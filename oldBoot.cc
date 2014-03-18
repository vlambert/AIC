#include "BootStrapping.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooTrace.h"
#define NUM_CPU 1

#define __DO_TRACE 0

using namespace std;

BootStrapping::BootStrapping(const TString& inputFileName)
{
  if(inputFileName != ""){
    inputFile = new TFile(inputFileName);
    w = ((RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace"));
  }
}


void BootStrapping::makeFits(RooWorkspace *wout) {
  RooAbsData *dc = w->data("Data_Combined");
  RooRealVar *mass = w->var("mass");
  
  int NData = dc->sumEntries();
  RooDataSet *Sample = new RooDataSet("Sample","",RooArgSet(*mass));
  int NSample = 0;
  while (NSample < NData) {
    int value = rand() % NData;
    double DataEntry = dc->get(value)->getRealValue("mass");
    mass->setVal(DataEntry);
    Sample->add(RooArgSet(*mass));
    NSample = Sample->sumEntries();
  }
  
  std::vector<RooExtendPdf*> ModelComp;
  for (int type=0; type<nModels; type++) {
    RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(type,mass);
  
    int k = MakeAICFits::Num_Params(type);
    RooRealVar *Nbkg = new RooRealVar("Nbkg","N Background Events", NData,0,1e9);
    RooExtendPdf *BkgModel = new RooExtendPdf("BKGFIT_bkgModel", "Background Model", *ModelShape, *Nbkg);
    
    RooFitResult *bkg_databkg = BkgModel->fitTo(*Sample, RooFit::Save(kTRUE), RooFit::Optimize(0));
    Double_t bkg_databkg_Nll = bkg_databkg->minNll();
    
    // set all fit parameters as constant
    RooArgSet* vars = BkgModel->getVariables();
    RooFIter iter = vars->fwdIterator();
    RooAbsArg* a;
    while( (a = iter.next()) ){
      if(string(a->GetName()).compare("mass")==0) continue;
      static_cast<RooRealVar*>(a)->setConstant(kTRUE);
    }
    // Add Model to Composite Array
    
    ModelComp.push_back(BkgModel);

    // Calculate AIC for each model
    LogLikelihood[type] += 2.*bkg_databkg_Nll;
    AIC_bkg_array[type] += 2.*(k + k*(k + 1.)/(NData - (k + 1.)) + bkg_databkg_Nll);
       
    // Clean up objects
    delete bkg_databkg;
    delete ModelShape;
    bkg_databkg_Nll = 0.;
  }
  //delete databkg;
  std::cout<<"Finished with fits for toy sample"<<std::endl;
  // Find minimum AIC value
  for (int type = 0; type<nModels;type++) {
    if (AIC_bkg_array[type] < minAIC) {
      minAIC = AIC_bkg_array[type];
    }
  }

  // Calculate Delta_i values and choose "best" model
  double sumExp = 0;
  for (int type = 0; type<nModels;type++) {
    DeltaI[type] = AIC_bkg_array[type] - minAIC;
    sumExp+= exp(-DeltaI[type]/2);
  }
  // Calculate AIC Weights
  for (int type = 0;type<nModels; type++) {
    AICweights[type] = exp(-DeltaI[type]/2)/sumExp;
    assert(ModelComp[type] != 0);
  }
  // Build Composite Model
  std::cout<<"Creating Composite Model"<<std::endl;
  //RooAbsPdf* singExp = ModelComp[0];
  //RooAbsPdf* doubExp = ModelComp[1];
  //RooAbsPdf* tripExp = ModelComp[2];
  //RooAbsPdf* modExp  = ModelComp[3];
  //RooAbsPdf* Poly    = ModelComp[4];
  //RooAbsPdf* singPow = ModelComp[5];
  //RooAbsPdf* doubPow = ModelComp[6];
  RooArgList* ModelSet = new RooArgList(*ModelComp[0],*ModelComp[1],*ModelComp[2],*ModelComp[3],*ModelComp[4],*ModelComp[5],*ModelComp[6]);
  //RooArgList* ModelSet = new RooArgList(*singExp,*doubExp,*tripExp,*modExp,*Poly,*singPow,*doubPow);

  for (int i=0; i<7;i++) {
    std::cout<<AICweights[i]<<std::endl;
  }
  
  RooRealVar* w1 = new RooRealVar("SingExpW","SingExpW",AICweights[0]);
  RooRealVar* w2 = new RooRealVar("DoubExpW","DoubExpW",AICweights[1]);
  RooRealVar* w3 = new RooRealVar("TripExpW","TripExpW",AICweights[2]);
  RooRealVar* w4 = new RooRealVar("ModExpW", "ModExpW" ,AICweights[3]);
  RooRealVar* w5 = new RooRealVar("PolyW", "PolyW",AICweights[4]);
  RooRealVar* w6 = new RooRealVar("SingPowW","SingPowW",AICweights[5]);
  RooRealVar* w7 = new RooRealVar("DoubPowW","DoubPowW",AICweights[6]);
  //RooArgList* Weights = new RooArgList(*w1,*w2,*w3,*w4,*w5,*w6,*w7);

  std::cout<<"All params set and made"<<std::endl;
  //RooAbsPdf *JointModel = new RooAddPdf("Composite", "Background Model",*ModelSet,*Weights);
  RooAbsPdf* JointModel = new RooAddPdf("C","B", RooArgList(*ModelComp[0],*ModelComp[1],*ModelComp[2],*ModelComp[3]),RooArgList(*w1,*w2,*w3));
  std::cout<<"Composite Model made"<<std::endl;
  CompositeModels.push_back(JointModel);
  wout->import(*JointModel);
  std::cout<<"Model stored in workspace"<<std::endl;
  delete JointModel;
}


void BootStrapping::run() {
  outfile.open("BootStrappingLogs.txt");
  for (int i = 0; i<NToys; i++) {
    RooWorkspace *wout = new RooWorkspace("AIC");
    makeFits(wout);
    std::cout<<"going to print"<<std::endl;
    print(i);
    wout->writeToFile(Form("AICModels_%d.root",i));
    
  }
  outfile.close();
  //plot();
}

void BootStrapping::plot() {
  RooRealVar *mass = w->var("mass");
  RooPlot* frame = mass->frame();
  TLegend *l = new TLegend(0.55,0.55,0.9,0.9, NULL, "NDC");
  for (auto model : CompositeModels) {
    model->plotOn(frame);
  }
  TCanvas *c = new TCanvas("", "", 800, 600);
  frame->Draw();
  l->Draw();
  c->Update();
  c->Print("CompositeModels.pdf");
}

void BootStrapping::print(int toy) {
  outfile<<"Toy Entry: " << toy << " :  Printing AIC Values" << std::endl;
  for (int type = 0; type<nModels; type++) {
    outfile<<"Log Likelihood for Model " << type << " ==== " << LogLikelihood[type] <<std::endl;
    outfile<<"AIC Value for Model: " << type << " ==== " << AIC_bkg_array[type]  << std::endl;
  }
  outfile<<"Minimum AIC: " << minAIC << std::endl;
  for (int type = 0; type<nModels; type++) {
    outfile<<"Delta AIC values : " << type << " ====" << DeltaI[type] <<std::endl;
  }
  for (int type = 0; type<nModels; type++) {
    outfile<<"Relative Likelihood AIC " << type << " ==== " <<exp(-DeltaI[type]/2)<<std::endl;
  }
  for (int type = 0; type<nModels; type++) {
    outfile<<"AIC Weights  : " << type << " =====" << AICweights[type] <<std::endl;
  }
}
  
  
  
  

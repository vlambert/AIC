// ========================================================== //
//           Bootstrapping analysis for stability of          //
//              H->gg background composite modeling           //     
//               Valere Lambert, Caltech 2014                 //                                  
// ========================================================== //
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


RooAbsPdf* BootStrapping::makeFits(int type , RooDataSet *Sample, RooRealVar *Nbkg) {
  // Returns a fitted model for a given data sample and model type
  // after calculating the AIC value                              

  // Build Model Shape
  RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(type,mass);
  
  // Set number of free parameters for the fit
  int k = MakeAICFits::Num_Params(type);

  // Extend Background PDF
  RooExtendPdf *BkgModel = new RooExtendPdf("BKGFIT_bkgModel", "Background Model", *ModelShape, *Nbkg);
    
  // Fit Model to Sample data and calculate min Log Likelihood
  int covariance = 0;
  RooFitResult *bkg_databkg;
  BkgModel->fitTo(*Sample,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  while (covariance < 2) {
    bkg_databkg = BkgModel->fitTo(*Sample, RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    covariance = bkg_databkg->covQual();
  }
  Double_t bkg_databkg_Nll = bkg_databkg->minNll();
    
  // Calculate AIC for each model
  LogLikelihood[type] = 2.*bkg_databkg_Nll;
  AIC_bkg_array[type] = 2.*(k + k*(k + 1.)/(NData - (k + 1.)) + bkg_databkg_Nll);
  
  // Clean up objects
  delete bkg_databkg;
  bkg_databkg_Nll = 0;
  
  return BkgModel;
}

void BootStrapping::run() {
  dc = w->data("Data_Combined");
  mass = w->var("mass");
  NData = dc->sumEntries();
  
  // Create weight objects
  RooRealVar* w1 = new RooRealVar("SingExpW","SingExpW",0);
  RooRealVar* w2 = new RooRealVar("DoubExpW","DoubExpW",0);
  RooRealVar* w3 = new RooRealVar("TripExpW","TripExpW",0);
  RooRealVar* w4 = new RooRealVar("ModExpW", "ModExpW" ,0);
  RooRealVar* w5 = new RooRealVar("PolyW", "PolyW",0);
  RooRealVar* w6 = new RooRealVar("SingPowW","SingPowW",0);
  RooRealVar* w7 = new RooRealVar("DoubPowW","DoubPowW",0);
  
  Double_t minAIC = 10000000000.0;  
  // Open text file for logs
  outfile.open("BootStrappingLogs.txt");
  
  // Set up for model overlay for each toy
  //TCanvas *cv = 0;
  //TLegend *l = 0;
  //RooPlot* f1 = 0;
  
  // Set up histogram for Signal Bias values
  TH1D * SigBias = new TH1D ("SignalBias", "",50,0.087,0.091);
  
  
  for (int i = 0; i<NToys; i++) {
    std::cout<< "New Toys being created " << "Toy : " << i << std::endl;
    w1->setVal(0);
    w2->setVal(0);
    w3->setVal(0);
    w4->setVal(0);
    w5->setVal(0);
    w6->setVal(0);
    w7->setVal(0);
    
    // Create new workspace to output model for toy sample
    //RooWorkspace *wout = new RooWorkspace("AIC");
    //cv = new TCanvas("","",800,600);
    //l = new TLegend(0.55,0.55,0.90,0.85);
    //l->SetTextSize(0.03);
    //f1 = mass->frame();
    
    // Generate data sample from original data set
    RooDataSet *Sample = new RooDataSet("Sample","",RooArgSet(*mass));
    int NSample = Sample->sumEntries();
    assert(NSample == 0);
    
    while (NSample < NData) {
      int value = rand() % NData;
      double DataEntry = dc->get(value)->getRealValue("mass");
      mass->setVal(DataEntry);
      Sample->add(RooArgSet(*mass));
      NSample = Sample->sumEntries();
    }
    //Sample->plotOn(f1);
    
    // Fit sample and produce pdfs for various models
    RooRealVar *Nbkg_0 = new RooRealVar("Nbkg_0","N Background Events",NData,0,1e9);
    RooAbsPdf* singExp = BootStrapping::makeFits(0, Sample, Nbkg_0);
    singExp->SetNameTitle(Form("singExp_%d",i),"Single Exponential");
    //singExp->plotOn(f1,RooFit::Name(Form("singExp_%d",i)),RooFit::LineColor(kRed));
    //l->AddEntry(f1->findObject(Form("singExp_%d",i)),"Single Exponential","l");
    
    RooRealVar *Nbkg_1 = new RooRealVar("Nbkg_1","N Background Events", NData, 0,1e9);
    RooAbsPdf* doubExp = BootStrapping::makeFits(1, Sample, Nbkg_1);
    doubExp->SetNameTitle(Form("doubExp_%d",i),"Double Exponential");
    //doubExp->plotOn(f1,RooFit::Name(Form("doubExp_%d",i)),RooFit::LineColor(kBlue));
    //l->AddEntry(f1->findObject(Form("doubExp_%d",i)),"Double Exponential","l");
    
    RooRealVar *Nbkg_2 = new RooRealVar("Nbkg_2","N Background Events", NData, 0,1e9);
    RooAbsPdf* tripExp = BootStrapping::makeFits(2, Sample, Nbkg_2);
    tripExp->SetNameTitle(Form("tripExp_%d",i),"Triple Exponential");
    //tripExp->plotOn(f1,RooFit::Name(Form("tripExp_%d",i)),RooFit::LineColor(38));
    //l->AddEntry(f1->findObject(Form("tripExp_%d",i)),"Triple Exponential","l");
    
    RooRealVar *Nbkg_3 = new RooRealVar("Nbkg_3","N Background Events", NData, 0,1e9);
    RooAbsPdf* modExp  = BootStrapping::makeFits(3, Sample, Nbkg_3);
    modExp->SetNameTitle(Form("modExp_%d",i),"Modified Exponential");
    //modExp->plotOn(f1,RooFit::Name(Form("modExp_%d",i)),RooFit::LineColor(kMagenta));
    //l->AddEntry(f1->findObject(Form("modExp_%d",i)),"Modified Exponential","l");
    
    RooRealVar *Nbkg_4 = new RooRealVar("Nbkg_4","N Background Events", NData, 0,1e9);
    RooAbsPdf* Poly    = BootStrapping::makeFits(4, Sample, Nbkg_4);
    Poly->SetNameTitle(Form("Poly_%d",i),"5th Order Polynomial");
    //Poly->plotOn(f1,RooFit::Name(Form("Poly_%d",i)),RooFit::LineColor(7));
    //l->AddEntry(f1->findObject(Form("Poly_%d",i)),"5th Order Polynomial","l");
    
    RooRealVar *Nbkg_5 = new RooRealVar("Nbkg_5","N Background Events", NData, 0,1e9);
    RooAbsPdf* singPow = BootStrapping::makeFits(5, Sample, Nbkg_5);
    singPow->SetNameTitle(Form("singPow_%d",i),"Single Power Law");
    //singPow->plotOn(f1,RooFit::Name(Form("doubPow_%d",i)),RooFit::LineColor(kGreen));
    //l->AddEntry(f1->findObject(Form("singPow_%d",i)),"Single Power Law","l");
    
    RooRealVar *Nbkg_6 = new RooRealVar("Nbkg_6","N Background Events", NData, 0,1e9);
    RooAbsPdf* doubPow = BootStrapping::makeFits(6, Sample, Nbkg_6);
    doubPow->SetNameTitle(Form("doubPow_%d",i),"Double Power Law");
    //doubPow->plotOn(f1,RooFit::Name(Form("doubPow_%d",i)),RooFit::LineColor(41));
    //l->AddEntry(f1->findObject(Form("doubPow_%d",i)),"Double Power Law","l");
    
    
    RooArgList* ModelSet = new RooArgList(*singExp,*doubExp,*tripExp,*modExp,*Poly,*singPow,*doubPow);
    std::cout<< "Making AIC"<<std::endl;
    minAIC = 1000000000.0;
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
    }   

    // Set AIC weight values
    w1->setVal(AICweights[0]);
    w2->setVal(AICweights[1]);
    w3->setVal(AICweights[2]);
    w4->setVal(AICweights[3]);
    w5->setVal(AICweights[4]);
    w6->setVal(AICweights[5]);
    w7->setVal(AICweights[6]);
    RooArgList* Weights = new RooArgList(*w1,*w2,*w3,*w4,*w5,*w6,*w7);

    //Build Composite Model
    RooAbsPdf *JointModel = new RooAddPdf("Composite", "Background Model", *ModelSet, *Weights);

    RooArgSet* vjoint = JointModel->getVariables();
    RooFIter iter = vjoint->fwdIterator();
    RooAbsArg* a;
    while ( ( a=iter.next()) ){
      if(string(a->GetName()).compare("mass")==0) continue;
      static_cast<RooRealVar*>(a)->setConstant(kTRUE);
    }
    // fit composite model to data
    RooRealVar *NBkg = new RooRealVar("NBkg","N Background Events",NData,0,1e9);
    RooExtendPdf *JointModelExt = new RooExtendPdf("","",*JointModel,*NBkg);
    JointModelExt->SetNameTitle(Form("Composite_%d",i),"Background Model");
    JointModelExt->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    RooFitResult *composite_databkg = JointModelExt->fitTo(*dc,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

    // Determine number of background events in signal region
    mass->setRange("SigRegion",123.0,127.0);
    mass->setRange("Full",mass->getMin(),mass->getMax());
    RooAbsReal *NSig = JointModelExt->createIntegral(*mass,RooFit::Range("SigRegion"),RooFit::NormSet(*mass));
    RooAbsReal *NAll = JointModelExt->createIntegral(*mass,RooFit::Range("Full"),RooFit::NormSet(*mass));
    NSigRange[i] = NSig->getVal() / NAll->getVal();
    NSigRangeE[i] = NSigRange[i] * sqrt(1/NAll->getVal() + pow(NSig->getPropagatedError(*composite_databkg)/NAll->getVal(),2));
    SigBias->Fill(NSigRange[i]);
    //SigBias->SetBinError(double(i),NSigRangeE[i]);

    CompositeM.push_back(JointModelExt);

    //Plot joint model and save it to workspace
    //JointModelExt->plotOn(f1,RooFit::Name("CompModel"),RooFit::LineColor(28));
    //l->AddEntry(f1->findObject("CompModel"),"Composite Model","l");

    LL[i] = 2*(composite_databkg->minNll());
    //wout->import(*JointModelExt);
    
    // Print logs and save workspace to a ROOT file
    print(i);
    
    //f1->Draw();
    //l->Draw();
    //cv->Update();   
    //cv->SaveAs(Form("Workspaces/overlayPlots_%d.pdf",i),"RECREATE");

    
    //outputFile = new TFile(Form("Workspaces/AICModels_%d.root",i),"RECREATE");
    //outputFile->cd();
    //wout->Write(wout->GetName(),TObject::kWriteDelete);
    //outputFile->Close();

    RooArgSet* vjointend = JointModelExt->getVariables();
    RooFIter iterend = vjointend->fwdIterator();
    RooAbsArg* aend;
    while ( ( aend=iterend.next()) ){
      if(string(aend->GetName()).compare("mass")==0) continue;
      static_cast<RooRealVar*>(aend)->setConstant(kFALSE);
    }
    
    delete composite_databkg;
    delete ModelSet;
    delete Weights;
    delete NAll;
    delete NSig;
    delete Sample;
    //delete wout;
  }
  
  //RooPlot *frame = mass->frame();
  TCanvas *c = 0;
  //c = new TCanvas("","",800,600);
  //RooDataHist *dh = new RooDataHist("dh","dh",RooArgSet(*mass),*dc);
  //dh->plotOn(frame,RooFit::Name("Data"));
  //std::vector<Int_t> color = {2,3,4,6,5,28,8,9,12,7};
  //for (int j = 0; j<NToys; j++) {
  //  CompositeM[j]->plotOn(frame,RooFit::Name(Form("CompModel_%d",j)),RooFit::LineColor(color[j]));
  //}
  
  //frame->Draw();
  for (int k = 0; k<NToys;k++) {
    //Chi2[k] = frame->chiSquare(Form("CompModel_%d",k),"Data");
    outfile<< "Chi2 for Toy Comp Model "<<k<< ":   " << Chi2[k] << std::endl; 
  }
  for (int k = 0; k<NToys; k++) {
    outfile<< "Signal Bias for Model " <<k<< ":  "<< NSigRange[k] << " +/- " << NSigRangeE[k] << std::endl;
  }

  //c->SaveAs("CompositeModels.pdf");


  c = new TCanvas("","",800,600);
  SigBias->Draw();
  TFile *histofile = new TFile("BootStrap_SigBiasHistos.root","RECREATE"); 
  SigBias->Write();
  histofile->Close();
  TF1 *func1 = new TF1("func1","gaus",0.087,0.091);
  TFitResultPtr r = SigBias->Fit(func1,"RS");
  gStyle->SetOptFit(0111);
  c->SaveAs("BootStrap_SigBias.pdf");
  
  outfile.close();
  std::cout<< "file closed"<<std::endl;
  inputFile->Close();
}

void BootStrapping::print(int toy) {
  outfile<<"Toy Entry: " << toy << " :  Printing AIC Values" << std::endl;
  /*for (int type = 0; type<nModels; type++) {
    outfile<<"Log Likelihood for Model " << type << " ==== " << LogLikelihood[type] <<std::endl;
    outfile<<"AIC Value for Model: " << type << " ==== " << AIC_bkg_array[type]  << std::endl;
  }
  for (int type = 0; type<nModels; type++) {
    outfile<<"Delta AIC values : " << type << " ====" << DeltaI[type] <<std::endl;
  }
  for (int type = 0; type<nModels; type++) {
    outfile<<"Relative Likelihood AIC " << type << " ==== " <<exp(-DeltaI[type]/2)<<std::endl;
  }
  for (int type = 0; type<nModels; type++) {
    outfile<<"AIC Weights  : " << type << " =====" << AICweights[type] <<std::endl;
  }
  outfile<<"Composite MinNll      : "  << LL[toy] <<std::endl;
  */
}
 
  
  
  

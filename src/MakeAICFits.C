// ========================================================== //
//    Multimodel Averaging for H->gg background modeling,     //
//    including construction of composite background model    //
//               Valere Lambert, Caltech 2014                 //
// ========================================================== //
#include "MakeAICFits.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooTrace.h"
#define NUM_CPU 1
#define __DO_TRACE 0

using namespace std;

MakeAICFits::MakeAICFits(const TString& inputFileName)
{
  if(inputFileName != ""){
    inputFile = new TFile(inputFileName);
    w = ((RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace"));
    // Get category labels from input workspace
    getLabels("evtcat",&catLabels,w);
  }
}

void MakeAICFits::getLabels(const char *varName, std::vector<TString> *lblVec, RooWorkspace* w){
  // Grab category labels from the data workspace
  RooCategory* labels = ((RooCategory*)w->obj(varName));
  lblVec->clear();
  if(labels==0) return;
  for(int i=0;i<labels->numBins("");i++){
    labels->setIndex(i);
    lblVec->push_back(labels->getLabel());
    std::cout<<lblVec->back() <<std::endl;
  }
}

int MakeAICFits::Num_Params(int type) {
  // Number of estimable parameters
  int Parameters;     
  switch(type){
  case kSingleExp:
    { 
      //single exponential
      Parameters = 2;
      break;
    }
  case kDoubleExp:
    {
      // double exponential
      Parameters = 4;
      break;
    }
  case kTripleExp:
    {
      // triple exponential
      Parameters = 6;
      break;
    }
  case kModifiedExp:
    { 
      // modified exponential 
      Parameters = 3; 
      break;
    } 
  case kPoly:
    { 
      // polynomial 
      Parameters = 6;
      break; 
    }
  case kPow:
    { 
      //Power 
      Parameters = 2;
      break; 
    }                                              
  case kDoublePow:
    { 
      // Double Power
      Parameters = 4;
      break;
    }
  default: 
    std::cout << "Invalid Background Model" << std::endl; 
    assert(false);
    break;  
  } 
  return Parameters; 
}        


RooAbsPdf* MakeAICFits::getBackgroundPdf(int type, RooRealVar* mass) {                  
  //background model
  RooAbsPdf* BkgShape;
  switch(type){
  case kSingleExp:
    { 
      //single exponential
      RooRealVar* alpha1 = new RooRealVar("alpha1","alpha1",-0.1,-1.,0.);                            
      BkgShape = new RooExponential("SingleExp","Background Model",*mass,*alpha1);                    
      break;
    }
  case kDoubleExp:                                                                                                                                           
    {                                                                                                                                                        
      //double exponential                                                                                                                                   
      RooRealVar* alpha1_2 = new RooRealVar("alpha1_2","alpha1_2",-0.1,-1.,0.);                            
      RooRealVar* alpha2_2 = new RooRealVar("alpha2_2","alpha2_22",-0.15,-1.,0.);                              
      RooRealVar* f1_bkg  = new RooRealVar("f1_bkg","f1_bkg",0.1,0,1);                                        
      RooExponential* exp1 = new RooExponential("exp1","exp1",*mass,*alpha1_2);                          
      RooExponential* exp2 = new RooExponential("exp2","exp2",*mass,*alpha2_2);

      BkgShape = new RooAddPdf("DoubleExp","Background Model",RooArgList(*exp1,*exp2),*f1_bkg);
      break;  
    }                                                                  
  case kTripleExp: 
    {
      //triple exponential  
      RooRealVar* alpha1_3 = new RooRealVar("alpha1_3","alpha1_3",-0.1,-1.,0.);                              
      RooRealVar* alpha2_3 = new RooRealVar("alpha2_3","alpha2_3",-0.15,-1.,0.);                              
      RooRealVar* alpha3_3 = new RooRealVar("alpha3_3","alpha3_3",-0.2,-1.,0.);                              
      RooRealVar* f1_2bkg  = new RooRealVar("f1_2bkg","f1_2bkg",0.1,0,1);                                     
      RooRealVar* f2_2bkg  = new RooRealVar("f2_2bkg","f2_2bkg",0.1,0,1);                                     
      RooExponential* exp1_3 = new RooExponential("exp1_3","exp1_3",*mass,*alpha1_3);                          
      RooExponential* exp2_3 = new RooExponential("exp2_3","exp2_3",*mass,*alpha2_3);                          
      RooExponential* exp3_3 = new RooExponential("exp3_3","exp3_3",*mass,*alpha3_3);                          
      
      BkgShape = new RooAddPdf("TripleExp","Background Model",                                          
                               RooArgList(*exp1_3,*exp2_3,*exp3_3),RooArgList(*f1_2bkg,*f2_2bkg));    
      break; 
    }
  case kModifiedExp:
    { 
      // pdf = e^(alpha1*m^alpha2)
      RooRealVar *alpha1_m = new RooRealVar("alpha1_m","",-1.,-10.,0.);                                    
      RooRealVar *alpha2_m = new RooRealVar("alpha2_m","",0.5,0.,10.);                                     
      BkgShape = new RooGenericPdf("ModifiedExp","Background Model","exp(@0*@1^@2)",RooArgList(*alpha1_m,*mass,*alpha2_m));   
      break;
    }                   
  case kPoly:
    { 
      //5th order polynomial
      RooRealVar *pC = new RooRealVar("pC","pC",1);                                                    
      RooRealVar *p0 = new RooRealVar("p0","p0",0,-10,10);                                             
      RooRealVar *p1 = new RooRealVar("p1","p1",0,-10,10);                                             
      RooRealVar *p2 = new RooRealVar("p2","p2",0,-10,10);                                             
      RooRealVar *p3 = new RooRealVar("p3","p3",0,-10,10);                                             
      RooRealVar *p4 = new RooRealVar("p4","p4",0,-10,10);                                             
      //enforce all coefficients positive 
      RooFormulaVar *pCmod = new RooFormulaVar("pCmod","pCmod","@0*@0",*pC);                                
      RooFormulaVar *p0mod = new RooFormulaVar("p0mod","p0mod","@0*@0",*p0);                                
      RooFormulaVar *p1mod = new RooFormulaVar("p1mod","p1mod","@0*@0",*p1);                                
      RooFormulaVar *p2mod = new RooFormulaVar("p2mod","p2mod","@0*@0",*p2);                                
      RooFormulaVar *p3mod = new RooFormulaVar("p3mod","p3mod","@0*@0",*p3);                                
      RooFormulaVar *p4mod = new RooFormulaVar("p4mod","p4mod","@0*@0",*p4);                                
      
      RooArgList *args; 
      args = new RooArgList(*pCmod,*p1mod,*p1mod,*p2mod,*p3mod,*p4mod);
      BkgShape = new RooBernstein("Polynomial","Background Model",*mass,*args);                          
      break;
    }
  case kPow:
    { 
      // pdf = m^alpha
      RooRealVar *alpha_p = new RooRealVar("alpha_p","alpha_p",-3.,-10.,0.);                                      
      BkgShape = new RooGenericPdf("Power","Background Model","@0^@1",RooArgList(*mass,*alpha_p));                    
      break; 
    }          
  case kDoublePow:
    { //pdf = f*m^alpha_1 + (1-f)*m^alpha_2  
      RooRealVar *alpha1_p2 = new RooRealVar("alpha1_p2","alpha1_p2",-4.0,-10.,0.);                                    
      RooRealVar *alpha2_p2 = new RooRealVar("alpha2_p2","alpha2_p2",-4.5,-10.,0.);                                   
      RooRealVar *f_pbkg  = new RooRealVar("f_pbkg","f_pbkg",0.1,0.05,0.95);                                       
      RooGenericPdf *pow1 = new RooGenericPdf("pow1","Power 1","@0^@1",RooArgList(*mass,*alpha1_p2));            
      RooGenericPdf *pow2 = new RooGenericPdf("pow2","Power 2","@0^@1",RooArgList(*mass,*alpha2_p2));
      
      BkgShape = new RooAddPdf("DoublePower","Background Model",RooArgList(*pow1,*pow2),*f_pbkg);                          
      break; 
    }
  default:
    std::cout << "INVALID BACKGROUND MODEL" << std::endl;
    assert(false); 
    break;
  }
  return BkgShape;
}         

RooAbsPdf* MakeAICFits::makeBackgroundFits(int type, RooRealVar* Nbkg) {
  // Build Fit model and calculate AIC
  RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(type,mass);
  int k = MakeAICFits::Num_Params(type);
  RooExtendPdf *BkgModel = new RooExtendPdf("BKGFIT_bkgModel", "Background Model", *ModelShape, *Nbkg);

  // Enforce that the fit converges
  int covariance = 0;
  RooFitResult *bkg_databkg;

  // covQual: -1 unknown, 0 not calculated , 1 approximation only, 2 full matrix but forced positive-definite, 3 full accurate covariance maxtrix
  while (covariance < 2) {
    BkgModel->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    bkg_databkg = BkgModel->fitTo(*dc, RooFit::Save(kTRUE), RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    covariance = bkg_databkg->covQual();
  }
  avgcov[type] = covariance;
  Double_t bkg_databkg_Nll = bkg_databkg->minNll();

  // Calculate AIC for each model
  LogLikelihood[type] = 2.*bkg_databkg_Nll;
  AIC_bkg_array[type] = 2.*(k + k*(k + 1.)/(SampleSize - (k + 1.)) + bkg_databkg_Nll);

  // Clean up objects
  delete bkg_databkg;
  bkg_databkg_Nll = 0.;

  return BkgModel;
}


void MakeAICFits::run() {
  //Set up categories
  for (auto catIt=catLabels.begin(); catIt !=catLabels.end();catIt++){
    runCategory(*catIt);
  }
  const TString &inclusive = TString("cat4");
  runCategory(inclusive);
}

void MakeAICFits::runCategory(const TString& catTag){
  if (catTag == "cat4"){
    dc = w->data("Data_Combined");
  }
  else{
    dc = w->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+catTag);
  }
  
  mass = w->var("mass");
  SampleSize = dc->sumEntries();
  
  // Create weight objects
  RooRealVar* w1 = new RooRealVar("SingExpW","SingExpW",0);
  RooRealVar* w2 = new RooRealVar("DoubExpW","DoubExpW",0);
  RooRealVar* w3 = new RooRealVar("TripExpW","TripExpW",0);
  RooRealVar* w4 = new RooRealVar("ModExpW", "ModExpW" ,0);
  RooRealVar* w5 = new RooRealVar("PolyW", "PolyW",0);
  RooRealVar* w6 = new RooRealVar("SingPowW","SingPowW",0);
  RooRealVar* w7 = new RooRealVar("DoubPowW","DoubPowW",0);
  
  // Open text file for logs
  outfile.open(TString("AICLogs_")+catTag+TString("_.txt"));
  
  // Set up for motivational plot with three model families and data
  RooPlot* frame = mass->frame();
  dc->plotOn(frame);
  TLegend *leg = new TLegend(0.55,0.55,0.9,0.9,NULL,"NDC");
  leg->AddEntry(dc,"Hgg Background","lep");
  
  // Set up for composite model plot
  RooPlot* compframe = mass->frame();
  TLegend *compleg = new TLegend(0.55,0.65,0.9,0.8,NULL,"NDC");
  
  // create workspace to output models
  RooWorkspace *wout = new RooWorkspace("AIC");
  
  // Fit sample and produce pdfs for various models
  RooRealVar *Nbkg_0 = new RooRealVar("Nbkg_0","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* singExp = MakeAICFits::makeBackgroundFits(0, Nbkg_0);
  singExp->SetNameTitle("singExp","Single Exponential");
  singExp->plotOn(frame,RooFit::Name("singExp"),RooFit::LineColor(kRed));
  leg->AddEntry(frame->findObject("singExp"),"Single Exponential","l");

  RooRealVar *Nbkg_1 = new RooRealVar("Nbkg_1","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* doubExp = MakeAICFits::makeBackgroundFits(1, Nbkg_1);
  doubExp->SetNameTitle("doubExp","Double Exponential");
  
  RooRealVar *Nbkg_2 = new RooRealVar("Nbkg_2","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* tripExp = MakeAICFits::makeBackgroundFits(2, Nbkg_2);
  tripExp->SetNameTitle("tripExp","Triple Exponential");
  
  RooRealVar *Nbkg_3 = new RooRealVar("Nbkg_3","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* modExp  = MakeAICFits::makeBackgroundFits(3, Nbkg_3);
  modExp->SetNameTitle("modExp","Modified Exponential");
  
  RooRealVar *Nbkg_4 = new RooRealVar("Nbkg_4","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* Poly    = MakeAICFits::makeBackgroundFits(4, Nbkg_4);
  Poly->SetNameTitle("Poly","5th Order Polynomial");
  Poly->plotOn(frame,RooFit::Name("Poly"),RooFit::LineColor(kBlue));
  leg->AddEntry(frame->findObject("Poly"),"5th Order Polynomial","l");
  
  RooRealVar *Nbkg_5 = new RooRealVar("Nbkg_5","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* singPow = MakeAICFits::makeBackgroundFits(5, Nbkg_5);
  singPow->SetNameTitle("singPow","Single Power Law");
  singPow->plotOn(frame,RooFit::Name("singPow"),RooFit::LineColor(kGreen));
  leg->AddEntry(frame->findObject("singPow"),"Single Power Law","l");
  
  RooRealVar *Nbkg_6 = new RooRealVar("Nbkg_6","N Background Events", SampleSize,0,1e9);
  RooAbsPdf* doubPow = MakeAICFits::makeBackgroundFits(6, Nbkg_6);
  doubPow->SetNameTitle("doubPow","Double Power Law");
  
  // Plot motivational plot with model families
  TCanvas *c = new TCanvas("","",800,600);
  frame->Draw();
  leg->Draw();
  c->Update();
  c->SaveAs(TString("HggData_")+catTag+TString(".pdf"));

  RooArgList* ModelSet = new RooArgList(*singExp,*doubExp,*tripExp,*modExp,*Poly,*singPow,*doubPow);
  // Find Minimum AIC value
  Double_t minAIC = 10000000000;
  for (int type = 0; type<nModels;type++) {
    if (AIC_bkg_array[type] < minAIC) {
      minAIC = AIC_bkg_array[type];
    }
  }
  // Calculate Delta_i values and choose "best" model
  Double_t sumExp = 0;
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
  RooRealVar *NBkg = new RooRealVar("NBkg","N Background Events",SampleSize,0,1e9);
  RooExtendPdf *JointModelExt = new RooExtendPdf("","",*JointModel,*NBkg);

  int CompCov = 0;
  RooFitResult *composite_databkg;
  JointModelExt->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  composite_databkg = JointModelExt->fitTo(*dc,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  CompCov = composite_databkg->covQual();
  
  //Plot joint model and save it to workspace
  RooDataHist *dh = new RooDataHist("dh","dh",RooArgSet(*mass),*dc);
  dh->plotOn(compframe,RooFit::Name("Data"));
  compleg->AddEntry(dh,"Hgg Background","lep");  
  JointModelExt->plotOn(compframe, RooFit::Name("CompModel"));
  compleg->AddEntry(compframe->findObject("CompModel"),"Composite Model","l");
  
  // Save composite model plot
  TCanvas *cv = new TCanvas("","",800,600);
  compframe->Draw();
  compleg->Draw();

  //Calculate Chi Square and Log Likelihood for Composite Model
  Double_t Chi2 = compframe->chiSquare("CompModel","Data");
  Double_t CompLL = 2*(composite_databkg->minNll());

  TPaveText *p1 = new TPaveText(0.55,0.6,0.9,0.65,"NDC");
  p1->AddText(Form("-2log(L) : %.3f",CompLL));
  p1->Draw();
  TPaveText *p2 = new TPaveText(0.55,0.55,0.9,0.6,"NDC");
  p2->AddText(Form("#chi^{2} : %.3f", Chi2));
  p2->Draw();
  cv->Update();
  cv->SaveAs(TString("HggAIC_composite_cat")+catTag+TString(".pdf"));
  
  wout->import(*JointModelExt);
  print();

  // Save workspace to file
  TFile *outputFile = new TFile(TString("Workspaces/HggAIC_workspace_")+catTag+TString(".root"),"RECREATE");
  outputFile->cd();
  wout->Write(wout->GetName(),TObject::kWriteDelete);
  outputFile->Close();

  // Print Composite results to category log file
  outfile <<"Composite Log Likelihood:  " << std::setprecision(10) << CompLL <<std::endl;
  outfile <<"Composite Chi^2:  "<< std::setprecision(5)<< Chi2 <<std::endl;
  outfile <<"Compositve covariance: "<< CompCov <<std::endl;
  outfile.close();

  // Release all composite model variables
  RooArgSet* vjointend = JointModelExt->getVariables();
  RooFIter iterend = vjointend->fwdIterator();
  RooAbsArg* aend;
  while ( ( aend=iterend.next()) ){
    if(string(aend->GetName()).compare("mass")==0) continue;
    delete static_cast<RooRealVar*>(aend);
  }

  // Clean up objects
  delete composite_databkg;
  delete frame;
  delete compframe;
  delete wout;
  delete vjoint;
  delete JointModel;
  delete JointModelExt;
  delete singExp;
  delete doubExp;
  delete tripExp;
  delete modExp;
  delete Poly;
  delete singPow;
  delete doubPow;
}


void MakeAICFits::print() {
  // Print results to log file
  outfile<<"Printing AIC Values" << std::endl;
  for (int type = 0; type<nModels; type++) {
    outfile<<"average covariance quality" << type <<" ===" << avgcov[type] <<std::endl;
    outfile<<"Log Likelihood for Model " << type << " ==== " << std::setprecision(10) << LogLikelihood[type] <<std::endl;
    outfile<<"AIC Value for Model: " << type << " ==== " << std::setprecision(10) << AIC_bkg_array[type] <<std::endl;
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
}
  
  
  
  

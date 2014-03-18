#include "MakeAICFits.h"
//#include "subtract.cc"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooTrace.h"
#define NUM_CPU 1

#define __DO_TRACE 0

using namespace std;


int MakeAICFits::Num_Params(int type) {                                                                                                                     
  int Parameters;                                                                                                                                           
  
  switch(type){                                                                                                                                              
  case kSingleExp:                                                                                                                                           
    {                                                                                                                                                        
      //single exponential                                                                                                                                   
      Parameters = 2; //1                                                                                                                                       
      break;                                                                                                                                                 
    }                                                                                                                                                        
  case kDoubleExp:                                                                                                                                           
    {                                                                                                                                                        
      // double exponential                                                                                                                                  
      Parameters = 4;    //3                                                                                                                                    
      break;                                                                                                                                                
    }                                                                                                                                                        
  case kTripleExp:                                                                                                                                           
    {                                                                                                                                                        
      // triple exponential
      Parameters = 6; //5
      break;                                                                                                                                                 
    }                                                                                                                                                        
  case kModifiedExp:                                                                                                                                         
    {                                                                                                                                                        
      // modified exponential                                                                                                                                
      Parameters = 3; //2                                                                                                                                        
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
      Parameters = 2; //1                                                                                                                                       
      break;                                                                                                                                                 
    }                                              
  case kDoublePow:                                                                                                                                           
    {                                                                                                                                                        
      // Double Power                                                                                                                                        
      Parameters = 4; //3                                                                                                                                        
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
      //std::cout<<"begin case single exponential"<<std::endl;
      RooRealVar* alpha1 = new RooRealVar("alpha1","alpha1",-0.1,-1.,0.);                            
      //std::cout<<"RooRealVar exists"<<std::endl;
      BkgShape = new RooExponential("BkgExp","Background Model",*mass,*alpha1);                    
      //std::cout<<"ending case single exponential"<<std::endl;
      break;                                                                                                                                                 
    }                                                                                                                                                        
  case kDoubleExp:                                                                                                                                           
    {                                                                                                                                                        
      //double exponential                                                                                                                                   
      RooRealVar* alpha1 = new RooRealVar("alpha1","alpha1",-0.1,-1.,0.);                            
      RooRealVar* alpha2 = new RooRealVar("alpha2","alpha2",-0.1,-1.,0.);                              
      RooRealVar* f_bkg  = new RooRealVar("f_bkg","f_bkg",0.1,0,1);                                        
      RooExponential* exp1 = new RooExponential("exp1","exp1",*mass,*alpha1);                          
      RooExponential* exp2 = new RooExponential("exp2","exp2",*mass,*alpha2);                          
                                                                                                                                                               
      BkgShape = new RooAddPdf("DoubleExp","Background Model",RooArgList(*exp1,*exp2),*f_bkg);                                                                                                
      break;                                                                                                                                                   
    }                                                                  
  case kTripleExp:                                                                                                                                             
    {                                                                                                                                                          
      //triple exponential                                                                                                                                     
      RooRealVar* alpha1 = new RooRealVar("alpha1","alpha1",-0.1,-1.,0.);                              
      RooRealVar* alpha2 = new RooRealVar("alpha2","alpha2",-0.1,-1.,0.);                              
      RooRealVar* alpha3 = new RooRealVar("alpha3","alpha3",-0.1,-1.,0.);                              
      RooRealVar* f1_bkg  = new RooRealVar("f1_bkg","f1_bkg",0.1,0,1);                                     
      RooRealVar* f2_bkg  = new RooRealVar("f2_bkg","f2_bkg",0.1,0,1);                                     
      RooExponential* exp1 = new RooExponential("exp1","exp1",*mass,*alpha1);                          
      RooExponential* exp2 = new RooExponential("exp2","exp2",*mass,*alpha2);                          
      RooExponential* exp3 = new RooExponential("exp3","exp3",*mass,*alpha3);                          
                                                                                                                                                               
      BkgShape = new RooAddPdf("TripleExp","Background Model",                                          
                               RooArgList(*exp1,*exp2,*exp3),RooArgList(*f1_bkg,*f2_bkg));                                                                     
      break;                                                                                                                                                   
    }                                                                                                                                                          
  case kModifiedExp:                                                                                                                                           
    { // pdf = e^(alpha1*m^alpha2)                                                                                                                             
      RooRealVar *alpha1 = new RooRealVar("alpha1","",-1.,-10.,0.);                                    
      RooRealVar *alpha2 = new RooRealVar("alpha2","",0.5,0.,10.);                                     
      BkgShape = new RooGenericPdf("ModifiedExp","Background Model","exp(@0*@1^@2)",RooArgList(*alpha1,*mass,*alpha2));   
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
    { // pdf = m^alpha                                                                                                                                         
      RooRealVar *alpha = new RooRealVar("alpha","alpha",-3.,-10.,0.);                                      
      BkgShape = new RooGenericPdf("Power","Background Model","@0^@1",RooArgList(*mass,*alpha));                    
      break;                                                                                                                                                   
    }          

  case kDoublePow:                                                                                                                                             
    { //pdf = f*m^alpha_1 + (1-f)*m^alpha_2                                                                                                                    
      RooRealVar *alpha1 = new RooRealVar("alpha1","alpha1",-4.,-10.,0.);                                    
      RooRealVar *alpha2 = new RooRealVar("alpha2","alpha2",-3.9,-10.,0.);                                   
      RooRealVar *f_bkg  = new RooRealVar("f_bkg","f_bkg",0.1,0.05,0.95);                                       
      RooGenericPdf *pow1 = new RooGenericPdf("pow1","Power 1","@0^@1",RooArgList(*mass,*alpha1));            
      RooGenericPdf *pow2 = new RooGenericPdf("pow2","Power 2","@0^@1",RooArgList(*mass,*alpha2));            
                                                                                                                                                               
      BkgShape = new RooAddPdf("DoublePower","Background Model",RooArgList(*pow1,*pow2),*f_bkg);                          
      break;                                                                                                                                                   
    }                                                                                                                                                          
                                                                                                                                                               
  default:                                                                                                                                                     
    std::cout << "INVALID BACKGROUND MODEL" << std::endl;                                                                                                      
    assert(false);                                                                                                                                             
    break;                                                                                                                                                     
  }                                                                                                                                                            
  return BkgShape;                                                                                                                                             
}         

MakeAICFits::MakeAICFits() {
  //RooRandom::randomGenerator()->SetSeed(314159);
  // Create Toy Data Set
  // Single Exponential
  RooRealVar *mass   = new RooRealVar("mass","mass", 50., 35., 65.); 
  //RooRealVar *alpha1 = new RooRealVar("Exp_alpha","Exp_alpha",-0.1,-1.,0.);                     
  //alpha1->setVal(-0.1);
  RooRealVar *alpha1 = new RooRealVar("Exp_alpha","Exp_alpha", -0.1);
  RooExponential bkg("exp","Background Exponential",*mass,*alpha1);             
  // KPower 
  //RooRealVar *alpha2 = new RooRealVar("Power_alpha","Power_alpha",-3.,-10.,0.);                             
  //alpha2->setVal(-3.);
  RooRealVar *alpha2 = new RooRealVar("Power_alpha", "Power_alpha", -3);
  RooGenericPdf pow("pow","Background Power","@0^@1",RooArgList(*mass,*alpha2));           
  //RooRealVar ratio("ratio","Background Ratio", 0.5, 0., 1.);
  //ratio.setVal(0.5);
  RooRealVar ratio("ratio", "Background Ratio", 0.5);


  //Create Background pdf
  //  RooAddPdf bkg("bkg","bkg",RooArgList(pow,expon),ratio);

  std::cout<<"==========  Data Model Made  ==========="<<std::endl;

  const Int_t nToys = 1;;
  //ratio.setVal(0);
  RooDataSet* data = bkg.generate(RooArgSet(*mass), 1000000);
  //ratio.setVal(0.20);


  std::cout<<"==========  Data Set Made    ==========="<<std::endl;
  
  // Make plain projection of data and pdf on mass
  //bkg.fitTo(*data);
  
  RooFitResult* bkg_data = bkg.fitTo(*data, RooFit::Save(kTRUE), RooFit::Optimize(0));
  Double_t bkg_data_Nll = bkg_data->minNll();
  std::cout<<" ======== Data fitted ==========="<<std::endl;
  RooArgSet* bkgParams = bkg.getParameters(*data);
  const RooArgList& fitbkgParams = bkg_data->floatParsFinal();
  std::cout<< "=======================  parameters done  ========================"<<std::endl;   

  Double_t LogLikelihood[8] = {0,0,0,0,0,0,0,0};
  Double_t AIC_bkg_array[8] = {0,0,0,0,0,0,0,0};
  Double_t AICc_bkg_array[8] = {0,0,0,0,0,0,0,0};
  //Double_t BIC_bkg_array[7] = {0,0,0,0,0,0,0};
  
  Double_t LogLikelihood_data = bkg_data_Nll;
  Int_t avgcov[8] = {0,0,0,0,0,0,0};
  std::cout<<"======================   Starting Toys  ==============="<<std::endl;
  for (Int_t i=0; i<nToys; i++) {
    if (i%10==0) {
      std::cout<<">> Processing Toy Trial " << i << "/" << nToys << std::endl;
    }

    RooPlot* frame = mass->frame();
    leg = new TLegend(0.55,0.55,0.9,0.9);
    mass->setConstant(kFALSE);
    alpha1->setConstant(kFALSE);
    alpha2->setConstant(kFALSE);
    ratio.setConstant(kFALSE);
    const RooArgList ranbkgParams = bkg_data->randomizePars();
    *bkgParams = ranbkgParams;
    
    Int_t SampleSize = 100000;
    RooDataSet* toybkg = bkg.generate(RooArgSet(*mass),SampleSize);
    toybkg->plotOn(frame);
    leg->AddEntry(toybkg,"Toy Background", "lep");
    *bkgParams = fitbkgParams;
    mass->setConstant(kTRUE);
    alpha1->setConstant(kTRUE);
    alpha2->setConstant(kTRUE);
    ratio.setConstant(kTRUE);
    for (int type=0; type<8; type++) {
      std::cout<<type<<endl;
    }
    for (int type=0; type<7; type++) {
      if (type<7) {
	//std::cout<<"Model Shape:    "<<type<<std::endl;
	RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(type,mass);
	//std::cout<<"Model Shape made"<<std::endl;
	int k = MakeAICFits::Num_Params(type);
	//std::cout<<"Params counted"<<std::endl;
      }
      if (type==7) {
	RooAbsPdf* Model1 = MakeAICFits::getBackgroundPdf(0,mass);
	RooAbsPdf* Model2 = MakeAICFits::getBackgroundPdf(1,mass);
	RooAbsPdf* Model3 = MakeAICFits::getBackgroundPdf(2,mass);
	RooAbsPdf* Model4 = MakeAICFits::getBackgroundPdf(3,mass);
	RooAbsPdf* Model5 = MakeAICFits::getBackgroundPdf(4,mass);
	int k = MakeAICFits::Num_Params(3);
	k+= MakeAICFits::Num_Params(1);
	k+= MakeAICFits::Num_Params(0);
	//k+= MakeAICFits::Num_Params(3);
	//k+= MakeAICFits::Num_Params(4);
	RooRealVar* modratio1 = new RooRealVar("modrat1", "modrat1", 0.62, 0.6, 0.7);
	RooRealVar* modratio2 = new RooRealVar("modrat2", "modrat2", 0.29, 0.25, 0.35);
	RooRealVar* modratio3 = new RooRealVar("modrat3", "modrat3", 0.01);
	//RooRealVar* modratio4 = new RooRealVar("modrat4", "modrat4", 0.25);
	RooAbsPdf* ModelShape = new RooAddPdf("Composite", "Background Model", RooArgList(*Model1, *Model4, *Model2), RooArgList(*modratio1, *modratio2));
	//RooAbsPdf* ModelShape = new RooAddPdf("Composite", "Background Model", RooArgList(*Model1, *Model4), *modratio1);
      }
      assert(ModelShape!=0);
      RooRealVar *Nbkg = new RooRealVar("Nbkg","N Background Events", SampleSize,0,1e9);
      RooExtendPdf *BkgModel = new RooExtendPdf("BKGFIT_bkgModel", "Background Model", *ModelShape, *Nbkg);
      TH1F* Model = new TH1F("Model", "Model", 100,0,100);
      //BkgModel->fitTo(*toybkg, RooFit::Strategy(0), RooFit::NumCPU(NUM_CPU), RooFit::Minos(kFALSE), RooFit::Extended(kTRUE));
      //RooFitResult *bkg_toybkg = BkgModel->fitTo(*toybkg,RooFit::Save(kTRUE), RooFit::Strategy(2), RooFit::NumCPU(NUM_CPU), RooFit::Minos(kFALSE), RooFit::Extended(kTRUE));
      RooFitResult *bkg_toybkg = BkgModel->fitTo(*toybkg, RooFit::Save(kTRUE), RooFit::Optimize(0));
      if (type == 0) {
	BkgModel->plotOn(frame, RooFit::LineColor(kBlue));
	Model->RooFit::SetLineColor(kBlue);
	leg->AddEntry(Model, "Exponential Model", "l");
      }
      if (type == 4) { 
	BkgModel->plotOn(frame, RooFit::LineColor(kRed));
	Model->RooFit::SetLineColor(kRed);
	leg->AddEntry(Model, "Polynomial Model", "l");
      }
      if (type == 5) { 
	BkgModel->plotOn(frame, RooFit::LineColor(kGreen));
	Model->RooFit::SetLineColor("kGreen");
	leg->AddEntry(Model, "Power Model", "l");
      }
      if (type == 7) {
	BkgModel->plotOn(frame, RooFit::LineColor(kMagenta));
	Model->RooFit::SetLineColor("kMagenta");
	leg->AddEntry(Model, "Composite Model", "l");
      }
      Double_t bkg_toybkg_Nll = bkg_toybkg->minNll();
      Int_t covariance = bkg_toybkg->covQual();
      avgcov[type] += covariance;
      //assert (covariance == 3);
      // Calculate AIC for each model
      LogLikelihood[type] += -bkg_toybkg_Nll;
      AICc_bkg_array[type] += 2.*(k + k*(k + 1.)/(SampleSize - (k + 1.)) + bkg_toybkg_Nll);
      //BIC_bkg_array[type]  += 2.*(k*log(SampleSize)/2. + bkg_toybkg_Nll);
      AIC_bkg_array[type] += 2.*(k + bkg_toybkg_Nll);
      // Clean up objects
      delete bkg_toybkg;
      bkg_toybkg_Nll = 0.;
    }
    delete toybkg;
    TCanvas *c = new TCanvas("", "", 800, 600);
    frame->Draw();
    leg->Draw();
    c->Update();
    c->Print("exp_plot_combined.pdf");
  }
  
  std::cout<<"Printing AIC Values" << std::endl;
  std::cout<<"Log Likelihood Data :    " << LogLikelihood_data <<std::endl;
  int display = 7;
  for (int type = 0; type<display; type++) {
    avgcov[type] = avgcov[type]/nToys;
    LogLikelihood[type] = LogLikelihood[type]/nToys;
    AIC_bkg_array[type] = AIC_bkg_array[type]/nToys;
    AICc_bkg_array[type] = AICc_bkg_array[type]/nToys;
    //BIC_bkg_array[type] = BIC_bkg_array[type]/nToys;
    std::cout<<"average covariance quality" << type <<" ===" << avgcov[type] <<std::endl;
    std::cout<<"Log Likelihood for Model " << type << " ==== " << LogLikelihood[type] <<std::endl;
    std::cout<<"AICc Value for Model: " << type << " ==== " << AICc_bkg_array[type] <<std::endl;
    std::cout<<"AIC Value for Model: " << type << " ==== " << AIC_bkg_array[type]  << std::endl;
    //std::cout<<"BIC Value for Model: " << type << " ==== " << BIC_bkg_array[type] <<std::endl;
  }
  double minAIC = 10000000000.;
  //double minBIC = 10000000000.;
  for (int type = 0; type<display; type++) {
    if (AICc_bkg_array[type] < minAIC) {
      minAIC = AICc_bkg_array[type];
    }
    //if (BIC_bkg_array[type] < minBIC) {
    //  minBIC = BIC_bkg_array[type];
    //}
  }
  std::cout<<"Minimum AIC: " << minAIC << std::endl;
  double DeltaIA[8];
  //double DeltaIB[7];
  double sumExpA=0;
  //double sumExpB=0;
  int bestmodelA;
  //int bestmodelB;
  for (int type = 0; type<display; type++) {
    DeltaIA[type] = AICc_bkg_array[type] - minAIC;
    //DeltaIB[type] = BIC_bkg_array[type] - minBIC;
    if (DeltaIA[type] == 0) {
      bestmodelA = type;
    }
    //if (DeltaIB[type] == 0) {
    //  bestmodelB = type;
    //}
    std::cout<<"Delta AIC values : " << type << " ====" << DeltaIA[type] <<std::endl;
    //std::cout<<"Delta BIC values : " << type << " ====" << DeltaIB[type] <<std::endl;
    sumExpA+= exp(-DeltaIA[type]/2);
    //sumExpB+= exp(-DeltaIB[type]/2);
  }
  double AICweights[8];
  //double BICweights[7];
  for (int type = 0; type<display; type++) {
    AICweights[type] = exp(-DeltaIA[type]/2)/sumExpA;
    //BICweights[type] = exp(-DeltaIB[type]/2)/sumExpB;
    std::cout<<"Relative Likelihood AIC " << type << " ==== " <<exp(-DeltaIA[type]/2)<<std::endl;
    std::cout<<"AIC Weights  : " << type << " =====" << AICweights[type] <<std::endl;
    //std::cout<<"Relative Likelihood BIC " << type << " ==== " <<exp(-DeltaIB[type]/2)<<std::endl;
    //std::cout<<"BIC Weights  : " << type << " =====" << BICweights[type] <<std::endl;
  }
  for (int type2 = 0; type2<display; type2++) {
    std::cout<< "AIC Ratio for:  " << "Model " << bestmodelA << " / " << "Model " << type2 << "  =  " << AICweights[bestmodelA]/AICweights[type2] <<std::endl;
    //std::cout<< "BIC Ratio for:  " << "Model " << bestmodelB << " / " << "Model " << type2 << "  =  " << BICweights[bestmodelB]/BICweights[type2] <<std::endl;
  }
}
  
  
  
  

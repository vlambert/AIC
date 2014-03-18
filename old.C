#include "MakeBiasStudy.h"
#include "MakeAICFits.h"
#include "assert.h"

using namespace std;

string MakeBiasStudy::Category(int type) {
  string Cat;
  switch(type){
  case kSingleExp:
    {
      Cat = "Single Exponential";
      break;
    }
  case kDoubleExp:
    {
      Cat = "Double Exponential";
      break;
    }
  case kTripleExp:
    {
      Cat = "Triple Exponential";
      break;
    }
  case kModifiedExp:
    {
      Cat = "Modified Exponential";
      break;
    }
  case kPoly:
    {
      Cat = "5th-Order Polynomial";
      break;
    }
  case kPow:
    {
      Cat = "Power Law";
      break;
    }
  case kDoublePow:
    {
      Cat = "Double Power Law";
      break;
    }
  default:
    std::cout << "Invalid Background Type" << std::endl;
    assert(false);
    break;
  }
  return Cat;
}
  

void MakeBiasStudy::biasStudy() {
  int Nmodels = 8;
  //RooRealVar* mass = ws->var("mass");
  TFile *f = TFile::Open("/mnt/hadoop/store/user/amott/Hgg2013/Hgg/workspaces/hgg_22Jan2013_R9_CIC.root");
  RooWorkspace *w = (RooWorkspace*) f->Get("cms_hgg_spin_workspace");
  RooDataSet *dc = w->data("Data_Combined");
  RooRealVar *mass = w->var("mass");
  RooRealVar* nBkgTruth("TruthNBkg","", 0,1e9);

  double Bias[Nmodels][Nmodels];
  double BiasE[Nmodels][Nmodels];

  for(int truthType = 0; truthType < Nmodels; truthType++){

    RooAbsPdf *truthPdf  = MakeAICFits::getBackgroundPdf(truthType,mass);
    RooExtendPdf truthExtendedPdf("truthExtendedPdf","",*truthPdf,nBkgTruth);
    truthExtendedPdf.fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    truthExtendedPdf.fitTo(*dc,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    
    RooArgSet *bkgpars = truthExtendedPdf.getVariables();
    RooFITer iter = bkg_pars->fwdIterator();
    RooAbsArg* a;
    while( (a = iter.next()) ) {
      if(string(a->GetName()).compare("mass")==0) continue;
      static_cast<RooRealVar*>(a)->setConstant(kTRUE);
    }
    double BiasWindow = 2.00;
    mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
    double TruthFrac = truthExtendedPdf.createIntegral(mass,RooFit::Range("biasRegion"),RooFit::NormSet(mass))->getVal();
    double NTruth = TruthFrac * nBkgTruth.getVal();
    double NTruthE = TruthFrac * nBkgTruth.getError();
    
    RooDataSet* truthbkg = truthPdf.generate(RooArgSet(*mass),NBkgTruth);

    for(int modelType = 0; modelType < Nmodels; modelType++){
      RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(modelType,mass);
      RooRealVar nBkgFit("FitNBkg", "", 0, 1e9);
      RooExtendPdf ModelExtendedPdf("ModelExtendedPdf", "",*ModelShape, nBkgFit);
      ModelExtendedPdf.fitTo(truthbkg, RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
      ModelExtendedPdf.fitTo(truthbkg, RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

      double FitFrac = ModelExtendedPdf.createIntegral(mass,RooFit::Range("biasRegion"),RooFit::NormSet(mass))->getVal();
      double NFit = FitFrac * nBkgFit.getVal();
      double NFitE = FitFrac * nBkgFit.getError();

      Bias[truthType][modelType] = fabs(NFit - NTruth);
      BiasE[truthType][modelType] = fabs(NFitE - NTruthE);
    }


  }

  for(int i = 0; i < Nmodels; i++) { 
    std::cout <<  "===== Truth Model : " << MakeBiasStudy::Category(i) << " ===== " << std::endl; 
    for (int j = 0; j < Nmodels; j++) { 
      std::cout << "Fit Model: " << MakeBiasStudy::Category(j) << "  , Bias = " << Bias[i][j] << " +/- " << BiasE[i][j] << std::endl; 
    }   
  }
  

}



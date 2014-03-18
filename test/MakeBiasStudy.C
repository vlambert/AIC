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
  

MakeBiasStudy::MakeBiasStudy() {
  int Nmodels = 8;
  //RooRealVar* mass = ws->var("mass");
  RooRealVar *mass = new RooRealVar("mass","mass", 100,180);
  RooRealVar *nBkgTruth = new RooRealVar("TruthNBkg","", 0,1e9);
  //  RooAbsData* realData = ws->data("Data_Combined")->reduce( Form("evtcat==evtcat::%s",cat.Data()) );

  double Bias[Nmodels][Nmodels];
  double BiasE[Nmodels][Nmodels];
  MakeAICFits MakeAIC_Fits;
  for(int truthType = 0; truthType < Nmodels; truthType++){

    RooAbsPdf *truthPdf = MakeAIC_Fits.getBackgroundPdf(truthType,mass);
    RooExtendPdf *truthExtendedPdf = new RooExtendPdf("truthExtendedPdf","",*truthPdf,*nBkgTruth);
    //truthExtendedPdf.fitTo(*realData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    //truthExtendedPdf.fitTo(*realData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    
    double BiasWindow = 2.00;
    mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
    double TruthFrac = truthExtendedPdf->createIntegral(mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass))->getVal();
    double NTruth = TruthFrac * nBkgTruth->getVal();
    double NTruthE = TruthFrac * nBkgTruth->getError();
    
    RooDataSet* truthbkg = truthPdf->generate(RooArgSet(*mass),nBkgTruth);

    for(int modelType = 0; modelType < Nmodels; modelType++){
      RooAbsPdf* ModelShape = MakeAIC_Fits.getBackgroundPdf(modelType,mass);
      RooRealVar *nBkgFit = new RooRealVar("FitNBkg", "", 0, 1e9);
      RooExtendPdf ModelExtendedPdf = new RooExtendPdf("ModelExtendedPdf", "",*ModelShape, *nBkgFit);
      ModelExtendedPdf.fitTo(truthbkg, RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
      ModelExtendedPdf.fitTo(truthbkg, RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

      double FitFrac = ModelExtendedPdf.createIntegral(mass,RooFit::Range("biasRegion"),RooFit::NormSet(mass))->getVal();
      double NFit = FitFrac * nBkgFit->getVal();
      double NFitE = FitFrac * nBkgFit->getError();

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



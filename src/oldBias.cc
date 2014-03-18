#include "MakeBiasStudy.h"

#include "assert.h"

MakeBiasStudy::MakeBiasStudy(const TString& inputFileName)
{
  //inputFile = TFile::Open("/mnt/hadoop/store/user/amott/Hgg2013/Hgg/workspaces/hgg_22Jan2013_R9_CIC.root");
  if(inputFileName !=""){
    inputFile = new TFile(inputFileName);
    w = ((RooWorkspace*) inputFile->Get("cms_hgg_spin_workspace")); 
    const char* icon = "Defaults";
    catLabels.push_back(icon);
  }
}
MakeBiasStudy::~MakeBiasStudy() {
}
  

void MakeBiasStudy::biasStudy(MakeAICFits::BkgFitType truthType, const TString& cat) {
  RooRealVar *mass = w->var("mass");
  RooAbsData *dc = w->data("Data_Combined");
  RooRealVar *nBkgTruth = new RooRealVar("TruthNBkg","", 0,1e9);
  
  RooAbsPdf *truthPdf = MakeAICFits::getBackgroundPdf(truthType,mass);
  RooExtendPdf *truthExtendedPdf = new RooExtendPdf("truthExtendedPdf","",*truthPdf,*nBkgTruth);
  truthExtendedPdf->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  truthExtendedPdf->fitTo(*dc,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    
  double BiasWindow = 2.00;
  mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
  mass->setRange("all",w->var("mass")->getMin(), w->var("mass")->getMax());
  double TruthRange = truthExtendedPdf->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass))->getVal();
  double TruthAll = truthExtendedPdf->createIntegral(*mass,RooFit::Range("all"),RooFit::NormSet(*mass))->getVal();
  double NTruth = TruthRange /TruthAll * nBkgTruth->getVal();
  double NTruthE = TruthRange /TruthAll * nBkgTruth->getError();
  
  RooDataSet* truthbkg = truthPdf->generate(RooArgSet(*mass),int(nBkgTruth->getVal()),RooFit::Extended(true));

  std::vector<std::vector<double>> Bias(FitTypes.size());
  std::vector<std::vector<double>> BiasErr(FitTypes.size());
  for (int i=0;i<NToys;i++) {
    auto fitIt = FitTypes.begin();
    auto BiasIt = Bias.begin();
    auto BiasEIt = BiasErr.begin();
    for(; fitIt!=FitTypes.end();fitIt++,BiasEIt++,BiasIt++) {
      std::tuple<double,double> results = getNFit(*truthbkg,*fitIt);
      BiasIt->push_back( ( fabs(std::get<0>(results) - NTruth)/std::get<1>(results)) );
      BiasEIt->push_back( fabs(std::get<0>(results)-NTruth) );
    }
  }
  auto fitIt = FitTypes.begin();
  auto BiasIt = Bias.begin();
  auto BiasEIt = BiasErr.begin();
  
  std::cout<< " ===== Truth Model : " << fitNameMap[truthType] << " ===== " << std::endl;
  for(; fitIt!=FitTypes.end();fitIt++,BiasEIt++,BiasIt++) {
    BiasBkg[cat][truthType].push_back( BiasIt->at( NToys/2 ));
    BiasBkgErr[cat][truthType].push_back( BiasEIt->at( NToys/2 ));

    std::cout << "  Fit: " << fitNameMap[*fitIt] << std::endl;
    for (auto bias: *BiasIt) {
      std::cout << bias << " ";
    }
    for (auto biasE: *BiasEIt) {
      std::cout << biasE << " ";
    }
   
    std::cout << std::endl;
  }
  delete truthPdf;; 
}
   
  
std::tuple<double,double> MakeBiasStudy::getNFit(RooAbsData& truthbkg, MakeAICFits::BkgFitType modelType) {
  RooRealVar* mass = w->var("mass");
  RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(modelType,mass);
  RooRealVar *nBkgFit = new RooRealVar("FitNBkg", "", 0, 1e9);
  RooExtendPdf *ModelExtendedPdf = new RooExtendPdf("ModelExtendedPdf", "",*ModelShape, *nBkgFit);
  ModelExtendedPdf->fitTo(truthbkg, RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  ModelExtendedPdf->fitTo(truthbkg, RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  
  double BiasWindow = 2.00;
  mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
  mass->setRange("all",w->var("mass")->getMin(),w->var("mass")->getMax());
  double FitRange = ModelExtendedPdf->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass))->getVal();
  double FitAll = ModelExtendedPdf->createIntegral(*mass,RooFit::Range("all"),RooFit::NormSet(*mass))->getVal();
  double NFit = FitRange / FitAll * nBkgFit->getVal();
  double NFitE = FitRange / FitAll * nBkgFit->getError();

  delete ModelShape;
  return std::tuple<double,double>(NFit,NFitE);
}

void MakeBiasStudy::run(){
  for (auto cat: catLabels) {
    for (auto truthIt : FitTypes) {
      biasStudy(truthIt,cat);
    }
  }
  print();
}

void MakeBiasStudy::print(){
  outfile.open("BiasOutput1.txt");
  std::cout << "\n \n ************************* BKG ERROR ***************************\n" << std::endl;
  printFormatted(BiasBkg);
  printFormatted(BiasBkgErr);
  outfile.close();
}

void MakeBiasStudy::printFormatted(const biasList& list) {
  const char *form = "%40s";
  for(auto catIt: catLabels) {
    for(auto truthIt: FitTypes) {
      outfile << fitNameMap[truthIt];
      for (auto BiasIt : list.at(catIt).at(truthIt) ) {
	outfile << Form("%40.2f \n",BiasIt);
      }
      outfile << std::endl;
    }
  }
}

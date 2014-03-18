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
    modelFile = TFile::Open("Workspaces/HggAIC_workspace.root");
    wmodels = ((RooWorkspace*) modelFile->Get("AIC")); 
  }
}
MakeBiasStudy::~MakeBiasStudy() {
}
  

void MakeBiasStudy::biasStudy(const TString& truthType, const TString& cat) {
  RooRealVar *mass = w->var("mass");
  RooAbsData *dc = w->data("Data_Combined");
  SampleSize = dc->sumEntries();
  RooRealVar *nBkgTruth = new RooRealVar("TruthNBkg","",SampleSize,0,1e9);

  RooAbsPdf *truthPdf = wmodels->pdf(truthType);
  RooExtendPdf *truthExtendedPdf = new RooExtendPdf("truthExtendedPdf","",*truthPdf,*nBkgTruth);

  truthExtendedPdf->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  truthExtendedPdf->fitTo(*dc,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

  double BiasWindow = 2.00;
  mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
  mass->setRange("FullRegion",mass->getMin(),mass->getMax());
  double TruthRange = truthExtendedPdf->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass))->getVal();
  double TruthAll = truthExtendedPdf->createIntegral(*mass,RooFit::Range("FullRegion"),RooFit::NormSet(*mass))->getVal();
  double NTruth = TruthRange/TruthAll * nBkgTruth->getVal();
  double NTruthE = TruthRange/TruthAll * nBkgTruth->getError();
  
  RooDataSet* truthbkg = truthExtendedPdf->generate(RooArgSet(*mass),int(nBkgTruth->getVal()),RooFit::Extended(true));

  std::vector<std::vector<double>> Bias(FitTypes.size());
  std::vector<std::vector<double>> BiasErr(FitTypes.size());
  std::vector<std::vector<double>> Slope(FitTypes.size());
  std::vector<std::vector<double>> SlopeNorm(FitTypes.size());

  for (int i=0;i<NToys;i++) {
    auto fitIt = FitTypes.begin();
    auto BiasIt = Bias.begin();
    auto BiasEIt = BiasErr.begin();
    auto SlopeIt = Slope.begin();
    auto SlopeNormIt = SlopeNorm.begin();
    for(; fitIt!=FitTypes.end();fitIt++,BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++) {
      std::tuple<double,double,double,double> results = getNFit(*truthbkg,*fitIt);
      BiasEIt->push_back( ( std::get<0>(results) - NTruth)/std::get<1>(results) );
      BiasIt->push_back( std::get<0>(results)-NTruth );
      SlopeIt->push_back( std::get<2>(results) );
      SlopeNormIt->push_back( std::get<3>(results) );
    }
  }
  auto fitIt = FitTypes.begin();
  auto BiasIt = Bias.begin();
  auto BiasEIt = BiasErr.begin();
  auto SlopeIt = Slope.begin();
  auto SlopeNormIt = SlopeNorm.begin();
  
  std::cout<< " ===== Truth Model : " << fitNameMap[truthType] << " ===== " << std::endl;
  for(; fitIt!=FitTypes.end();fitIt++,BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++) {
    std::sort( BiasIt->begin(), BiasIt->end() );
    std::sort( BiasEIt->begin(), BiasEIt->end() );
    std::sort( SlopeIt->begin(), SlopeIt->end() );
    std::sort( SlopeNormIt->begin(), SlopeNormIt->end() );

    BiasBkg[cat][truthType].push_back( BiasIt->at( NToys/2 ));
    BiasBkgErr[cat][truthType].push_back( fabs(BiasEIt->at( NToys/2 )));
    SlopeBkg[cat][truthType].push_back( SlopeIt->at( NToys/2 ));
    SlopeNormBkg[cat][truthType].push_back( SlopeNormIt->at( NToys/2 ));

    std::cout << "  Fit: " << fitNameMap[*fitIt] << std::endl;
    for (auto bias: *BiasIt) {
      std::cout << bias << " ";
    }
    for (auto biasE: *BiasEIt) {
      std::cout << biasE << " ";
    }
   
    std::cout << std::endl;
  }
}
   
  
std::tuple<double,double,double,double> MakeBiasStudy::getNFit(RooAbsData& truthbkg, const TString& modelType) {
  RooRealVar* mass = w->var("mass");

  RooAbsPdf *ModelShape = wmodels->pdf(modelType);
  RooRealVar *nBkgFit = new RooRealVar("FitNBkg", "",SampleSize, 0, 1e9);
  RooExtendPdf *ModelExtendedPdf = new RooExtendPdf("ModelExtendedPdf", "",*ModelShape, *nBkgFit);

  RooArgSet* vs = ModelExtendedPdf->getVariables();
  RooFIter it = vs->fwdIterator();
  RooAbsArg* a;
  while( (a = it.next()) ) {
    if(std::string(a->GetName()).compare("mass")==0) continue;
    static_cast<RooRealVar*>(a)->setConstant(kFALSE);
  }

  ModelExtendedPdf->fitTo(truthbkg, RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  ModelExtendedPdf->fitTo(truthbkg, RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

  double BiasWindow = 2.00;
  mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
  mass->setRange("FullRegion",mass->getMin(),mass->getMax() );
  double FitRange = ModelExtendedPdf->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass))->getVal();
  double FitAll = ModelExtendedPdf->createIntegral(*mass,RooFit::Range("FullRegion"),RooFit::NormSet(*mass))->getVal();
  double NFit = FitRange/FitAll * nBkgFit->getVal();
  double NFitE = FitRange/FitAll * nBkgFit->getError();

  mass->setRange("lowerend",mh-BiasWindow,mh);
  mass->setRange("upperend",mh,mh+BiasWindow);
  double LowInt = ModelExtendedPdf->createIntegral(*mass,RooFit::Range("lowerend"),RooFit::NormSet(*mass))->getVal();
  double HighInt = ModelExtendedPdf->createIntegral(*mass,RooFit::Range("upperend"),RooFit::NormSet(*mass))->getVal();
  double NLow = LowInt * nBkgFit->getVal();
  double NHigh = HighInt *nBkgFit->getVal();
  double Slope = (NHigh - NLow)/2;
  double SlopeNorm = (NHigh-NLow)/(2*FitRange*nBkgFit->getVal());

  return std::tuple<double,double,double,double>(NFit,NFitE,Slope,SlopeNorm);
}

void MakeBiasStudy::run(){
  for (auto cat: catLabels) {
    for (auto truthIt : FitTypes) {
      std::cout<<truthIt<<std::endl;
      biasStudy(truthIt,cat);
    }
  }
  print();
  std::cout<< "========== Finished Running Bias Study =========="<<std::endl;
}

void MakeBiasStudy::print(){
  outfile.open("BiasOutput.txt");
  outfile << "\n \n ************************* BKG ERROR *************************** \n" << std::endl;
  outfile << "  ======================= Bias =======================  " << std::endl;
  printFormatted(BiasBkg);
  outfile << "  ======================= Bias Error =======================  " <<std::endl;
  printFormatted(BiasBkgErr);
  outfile << "  ======================= Slope =======================  " <<std::endl;
  printFormatted(SlopeBkg);
  outfile <<"  =======================  Normalized Slope =======================  " << std::endl;
  printFormatted(SlopeNormBkg);
  outfile.close();
}

void MakeBiasStudy::printFormatted(const biasList& list) {
  /*
    Print out the bias results for each truth and fit model pair
                                                                    */

  for(auto catIt: catLabels) {
    for(auto truthIt: FitTypes) {
      outfile <<"Truth Type:  " <<  fitNameMap[truthIt] <<std::endl;
      for (auto BiasIt : list.at(catIt).at(truthIt) ) {
	outfile << Form("%40.3f \n",BiasIt) <<std::endl;
      }
      outfile << std::endl;
    }
  }
}

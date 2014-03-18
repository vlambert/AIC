#include "ShapeStudy.h"

#include "assert.h"

ShapeStudy::ShapeStudy(const TString& inputFileName)
{
  //inputFile = TFile::Open("/mnt/hadoop/store/user/amott/Hgg2013/Hgg/workspaces/hgg_22Jan2013_R9_CIC.root");
  if(inputFileName !=""){
    inputFile = new TFile(inputFileName);
    w = ((RooWorkspace*) inputFile->Get("cms_hgg_spin_workspace")); 
    const char* icon = "Defaults";
    catLabels.push_back(icon);
  }
}
ShapeStudy::~ShapeStudy() {
}
  

std::tuple<double,double> ShapeStudy::getBiasedMass(double k_value, RooDataSet* toydata) {
  RooRealVar* alpha = new RooRealVar("alpha","alpha",k_value);
  alpha->setVal(k_value);
  alpha->setConstant(kTRUE);
  RooExponential *Exp = new RooExponential("ExpFit","",*mass,*alpha);
  RooExtendPdf *ExpModel = new RooExtendPdf("ExpModel","",*Exp,*Nbkg);
  
  RooRealVar *mHFit = new RooRealVar("mHFit","mHFit",mh,120,130);
  RooRealVar *sHFit = new RooRealVar("sHFit","sHFit",1.5,0.2,5.0);
  RooGaussian *gaus = new RooGaussian("gaus","gaus",*mass, *mHFit,*sHFit);
  RooAddPdf *combModel = new RooAddPdf("combModel","",RooArgSet(*ExpModel,*gaus),RooArgSet(nBkgSB,nSigSB));
  
  combModel->fitTo(*toydata,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  combModel->fitTo(*toydata,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

  double mFit = mHFit->getVal();
  double mFitE = mHFit->getError();

  delete combModel;
  return std::tuple<double,double>(mFit,mFitE);
}

void ShapeStudy::run(){
  mass = w->var("mass");
  dc = w->data("Data_Combined");

  RooAbsPdf* Bkg = MakeAICFits::getBackgroundPdf(0,mass);
  SampleSize = dc->sumEntries();
  Nbkg = new RooRealVar("Nbkg","N Background Events",SampleSize,0,1e9);
  RooExtendPdf *BkgModel = new RooExtendPdf("BKGTrue_bkgModel", "Background Model",*Bkg,*Nbkg);

  RooRealVar *massH = new RooRealVar("massH","Higgs Mass",mh,123,127);
  RooRealVar *sigmatrue = new RooRealVar("sigma","",1.5,0.2,5.0);
  RooGaussian* Sig = new RooGaussian("Signal","Signal Peak",*mass,*massH,*sigmatrue);

  nBkgSB = new RooRealVar("NBkgSB","",0,1e9);
  nSigSB = new RooRealVar("NSigSB","",0,1e5);

  RooAddPdf *sumPdf = new RooAddPdf("sum","",RooArgSet(*BkgModel,*Sig),RooArgSet(nBkgSB,nSigSB));

  std::vector<std::vector<double,double>> HiggsM(NToys);
  std::vector<std::vector<double,double>> HiggsM_1s(NToys);
  std::vector<std::vector<double,double>> HiggsM_2s(NToys);
  std::vector<std::vector<double,double>> HiggsM_m1s(NToys);
  std::vector<std::vector<double,double>> HiggsM_m2s(NToys);

  for(int i = 0; i<NToys;i++){
    RooDataSet *toyData = sumPdf->generate(RooArgSet(*mass),int(Nbkg->getVal()),RooFit::Extended(true));

    RooRealVar* k = new RooRealVar("k","k",-0.1,-1,0);
    RooExponential *FitBkg = new RooExponential("Exp","",*mass,*k);
    RooExtendPdf *BkgFitModel = new RooExtendPdf("BKGFit_bkgModel","Background Model",*FitBkg,*Nbkg);

    RooRealVar *massHFit = new RooRealVar("massHFit","Fit Higgs Mass",mh,120,130);
    RooRealVar *SigmaFit = new RooRealVar("SigmaFit","Fit Sigma",1.5,0.2,5.0);
    RooGaussian *FitSig = new RooGaussian("FitSignal","Signal Peak",*mass,*massHFit,*SigmaFit);
    RooAddPdf *FitModel = new RooAddPdf("Fitsum","",RooArgSet(*BkgFitModel,*FitSig),RooArgSet(nBkgSB,nSigSB));
    FitModel->fitTo(*toyData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    FitModel->fitTo(*toyData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    
    double k_fit = k->getVal();
    double k_error = k->getError();
    double MassH = massHFit->getVal();
    double MassHE = massHFit->getError();

    std::tuple<double,double> MassH_k = getBiasedMass(k_fit + k_error,*toyData);
    std::tuple<double,double> MassH_2k = getBiasedMass(k_fit+2*k_error,*toyData);
    std::tuple<double,double> MassH_mk = getBiasedMass(k_fit-k_error,*toyData);
    std::tuple<double,double> MassH_m2k = getBiasedMass(k_fit-2*k_error,*toyData);
    HiggsM.push_back(std::tuple<double,double>(MassH,MassHE));
    HiggsM_1s.push_back(MassH_k);
    HiggsM_2s.push_back(MassH_2k);
    HiggsM_m1s.push_back(MassH_mk);
    HiggsM_m2s.push_back(MassH_m2k);
    delete FitModel;
    delete toyData;
  }
  
  print();
}

void ShapeStudy::print(){
  outfile.open("ShapeStudy.txt");
  outfile << "\n \n ************************* Shape Study *************************** \n" << std::endl;
  auto HIt = HiggsM.begin();
  auto HIt_1s = HiggsM_1s.begin();
  auto HIt_2s = HiggsM_2s.begin();
  auto HIt_m1s = HiggsM_m1s.begin();
  auto HIt_m2s = HiggsM_m2s.begin();
  outfile<<"Fit mH    :    k+sigma   :    k+2*sigma   :    k-sigma   :    k-2*sigma  "<<std::endl;
  for(; HIt!=HiggsM.end();HIt++,HIt_1s++,HIt_2s++,HIt_m1s++,HIt_m2s++) {
    outfile<< Form("%.2f + %.2f",std::get<0>(HIt),std::get<1>(HIt)) << Form("%.2f + %.2f",std::get<0>(HIt_1s),std::get<1>(HIt_1s)) << Form("%.2f + %.2f",std::get<0>(HIt_2s),std::get<1>(HIt_2s) << Form("%.2f + %.2f",std::get<0>(HIt_m1s),std::get<1>(HIt_m1s))) << Form("%.2f + %.2f",std::get<0>(HIt_m2s),std::get<1>(HIt_m2s)) << std::endl;
    outfile << << std::endl;
  }
  outfile.close();
}




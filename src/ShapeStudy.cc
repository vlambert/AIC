#include "ShapeStudy.h"

#include "assert.h"

ShapeStudy::ShapeStudy(const TString& inputFileName)
{
  if(inputFileName !=""){
    inputFile = new TFile(inputFileName);
    w = ((RooWorkspace*) inputFile->Get("cms_hgg_spin_workspace")); 
  }
}
ShapeStudy::~ShapeStudy() {
}
  
std::tuple<double,double> ShapeStudy::getBiasedMass(double k_value, RooAbsData& toydata) {
  
  // Create Exponential Background Model
  RooRealVar* alpha = new RooRealVar("alpha","alpha",k_value);
  alpha->setVal(k_value);
  alpha->setConstant(kTRUE);
  RooExponential *Exp = new RooExponential("ExpFit","",*mass,*alpha);
  RooExtendPdf *ExpModel = new RooExtendPdf("ExpModel","",*Exp,*Nbkg);
  
  // Create Signal Model
  RooRealVar *mHFit = new RooRealVar("mHFit","mHFit",mh,123,127);
  RooRealVar *sHFit = new RooRealVar("sHFit","sHFit",1.5,0.2,5.0);
  RooGaussian *gaus = new RooGaussian("gaus","gaus",*mass, *mHFit,*sHFit);

  // Combine Background and Signal Models
  RooAddPdf *combModel = new RooAddPdf("combModel","",RooArgSet(*ExpModel,*gaus),RooArgSet(*nBkgSB,*nSigSB));
  
  // Fit model to toy data
  combModel->fitTo(toydata,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  combModel->fitTo(toydata,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

  //Extract fit mass value and width
  double mFit = mHFit->getVal();
  double mFitE = mHFit->getError();
  
  alpha->setConstant(kFALSE);
  delete combModel;
  return std::tuple<double,double>(mFit,mFitE);
}

void ShapeStudy::run(){
  mass = w->var("mass");
  dc = w->data("Data_Combined");

  // Create Truth Exponential Model
  RooAbsPdf* Bkg = MakeAICFits::getBackgroundPdf(0,mass);
  SampleSize = dc->sumEntries();
  Nbkg = new RooRealVar("Nbkg","N Background Events",SampleSize,0,1e9);
  RooExtendPdf *BkgModel = new RooExtendPdf("BKGTrue_bkgModel", "Background Model",*Bkg,*Nbkg);

  // Create Truth Signal Model
  RooRealVar *massH = new RooRealVar("massH","Higgs Mass",mh,123,127);
  RooRealVar *sigmatrue = new RooRealVar("sigma","",1.5,0.2,5.0);
  RooGaussian* Sig = new RooGaussian("Signal","Signal Peak",*mass,*massH,*sigmatrue);

  nBkgSB = new RooRealVar("NBkgSB","",0,1e9);
  nSigSB = new RooRealVar("NSigSB","",0,1e5);

  // Combine Background and Signal Truth Models
  RooAddPdf *sumPdf = new RooAddPdf("sum","",RooArgSet(*BkgModel,*Sig),RooArgSet(*nBkgSB,*nSigSB));

  sumPdf->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  sumPdf->fitTo(*dc,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

  std::vector<std::tuple<double,double>> HiggsM;
  std::vector<std::tuple<double,double>> HiggsM_1s;
  std::vector<std::tuple<double,double>> HiggsM_2s;
  std::vector<std::tuple<double,double>> HiggsM_m1s;
  std::vector<std::tuple<double,double>> HiggsM_m2s;

  std::vector<std::tuple<double,double>> Diff_1s;
  std::vector<std::tuple<double,double>> Diff_2s;
  std::vector<std::tuple<double,double>> Diff_m1s;
  std::vector<std::tuple<double,double>> Diff_m2s;

  TH1D *p1s = new TH1D("1s","",50,122,128);
  TH1D *p2s = new TH1D("2s","",50,122,128);
  TH1D *m1s = new TH1D("m1s","",50,122,128);
  TH1D *m2s = new TH1D("m2s","",50,122,128);
  TH1D *norm = new TH1D("mH","",50,122,128);

  Double_t normE = 0.0;

  // Commence Toy Study
  for(int i = 0; i<NToys;i++){
    // Generate toy data from truth model
    RooDataSet *toyData = sumPdf->generate(RooArgSet(*mass),int(Nbkg->getVal()),RooFit::Extended(true));

    // Create fit Exponential model
    RooRealVar* k = new RooRealVar("k","k",-0.1,-1.0,0.0);
    RooExponential *FitBkg = new RooExponential("Exp","",*mass,*k);
    RooExtendPdf *BkgFitModel = new RooExtendPdf("BKGFit_bkgModel","Background Model",*FitBkg,*Nbkg);

    // Fit Exponential to toy data
    BkgFitModel->fitTo(*toyData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    BkgFitModel->fitTo(*toyData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    
    // Extract exponential power and set constant
    double k_fit = k->getVal();
    double k_error = k->getError();
    k->setConstant(kTRUE);
   
    // Combine Exponential with Signal Fit Model
    RooRealVar *massHFit = new RooRealVar("massHFit","Fit Higgs Mass",mh, 123,127);
    RooRealVar *SigmaFit = new RooRealVar("SigmaFit","Fit Sigma", 1.5,0.2, 5.0);
    RooGaussian *FitSig = new RooGaussian("FitSignal","Signal Peak",*mass, *massHFit,*SigmaFit);
    RooAddPdf *FitModel = new RooAddPdf("Fitsum","",RooArgSet(*BkgFitModel,*FitSig),RooArgSet(*nBkgSB,*nSigSB));

    // Fit toy data with combined model
    FitModel->fitTo(*toyData,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    FitModel->fitTo(*toyData,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

    // Extract Signal mass and width
    double MassH = massHFit->getVal();
    double MassHE = massHFit->getError();

    // Extract Signal mass and width for varied k values
    std::tuple<double,double> MassH_k   = getBiasedMass( k_fit +   k_error ,*toyData);
    std::tuple<double,double> MassH_2k  = getBiasedMass( k_fit + 2*k_error ,*toyData);
    std::tuple<double,double> MassH_mk  = getBiasedMass( k_fit -   k_error ,*toyData);
    std::tuple<double,double> MassH_m2k = getBiasedMass( k_fit - 2*k_error ,*toyData);
    HiggsM.push_back(std::tuple<double,double>(MassH,MassHE));
    normE += std::pow(MassHE,2);
    HiggsM_1s.push_back(MassH_k);
    HiggsM_2s.push_back(MassH_2k);
    HiggsM_m1s.push_back(MassH_mk);
    HiggsM_m2s.push_back(MassH_m2k);

    p1s->Fill(std::get<0>(MassH_k));
    p2s->Fill(std::get<0>(MassH_2k));
    m1s->Fill(std::get<0>(MassH_mk));
    m2s->Fill(std::get<0>(MassH_m2k));
    norm->Fill(MassH);

    std::tuple<double,double> D_1s  = std::tuple<double,double>(std::get<0>(MassH_k)-MassH, pow(pow(std::get<1>(MassH_k),2)+pow(MassHE,2),0.5));
    std::tuple<double,double> D_2s  = std::tuple<double,double>(std::get<0>(MassH_2k)-MassH, pow(pow(std::get<1>(MassH_2k),2)+pow(MassHE,2),0.5));
    std::tuple<double,double> D_m1s = std::tuple<double,double>(std::get<0>(MassH_mk)-MassH, pow(pow(std::get<1>(MassH_mk),2)+pow(MassHE,2),0.5));
    std::tuple<double,double> D_m2s = std::tuple<double,double>(std::get<0>(MassH_m2k)-MassH, pow(pow(std::get<1>(MassH_m2k),2)+pow(MassHE,2),0.5));  
    Diff_1s.push_back(D_1s);
    Diff_2s.push_back(D_2s);
    Diff_m1s.push_back(D_m1s);
    Diff_m2s.push_back(D_m2s);
    
    delete FitModel;
    delete toyData;
    std::cout<< "k_fit +/- k_error =  " << k_fit << "+/-" << k_error <<std::endl; 
  }
  
  Double_t avgD_1s  = 0.0;
  Double_t avgD_2s  = 0.0;
  Double_t avgD_m1s = 0.0;
  Double_t avgD_m2s = 0.0;
  
  Double_t D_1sE  = 0.0;
  Double_t D_2sE  = 0.0;
  Double_t D_m1sE = 0.0;
  Double_t D_m2sE = 0.0;
  
  std::cout<<"running averages"<<std::endl;
  auto DIt_1s  = Diff_1s.begin();
  auto DIt_2s  = Diff_2s.begin();
  auto DIt_m1s = Diff_m1s.begin();
  auto DIt_m2s = Diff_m2s.begin();
  // assert that vectors are all the same length
  std::cout<<"start loop"<<std::endl;
  for(; DIt_1s!=Diff_1s.end(); DIt_1s++, DIt_2s++, DIt_m1s++, DIt_m2s++) {
    
    std::cout<<"in loop"<<std::endl;
    
    avgD_1s+=std::get<0>(*DIt_1s);
    avgD_2s+=std::get<0>(*DIt_2s);
    avgD_m1s+=std::get<0>(*DIt_m1s);
    avgD_m2s+=std::get<0>(*DIt_m2s);
    
    std::cout<<" diff"<<std::endl;
    
    D_1sE += std::pow(std::get<1>(*DIt_1s),2);
    D_2sE += std::pow(std::get<1>(*DIt_2s),2);
    D_m1sE +=std::pow(std::get<1>(*DIt_m1s),2);
    D_m2sE +=std::pow(std::get<1>(*DIt_m2s),2);
  }
  
  std::cout<< "summations made" << std::endl;
  avgD_1s = avgD_1s/NToys;
  avgD_2s = avgD_2s/NToys;
  avgD_m1s = avgD_m1s/NToys;
  avgD_m2s = avgD_m2s/NToys;
  D_1sE = std::pow(D_1sE/NToys,0.5);
  D_2sE = std::pow(D_2sE/NToys,0.5);
  D_m1sE = std::pow(D_m1sE/NToys,0.5);
  D_m2sE = std::pow(D_m2sE/NToys,0.5);
  std::cout<<" up to norm"<<std::endl;
  normE = std::pow(normE/NToys,0.5);

  std::cout<<"making histogram"<<std::endl;
  TH1D *Diff = new TH1D("Mass Differences","",5,0,6);
  Diff->GetXaxis()->SetBinLabel(1.0,"k-2sigma");
  Diff->GetXaxis()->SetBinLabel(2.0,"k-2sigma");
  Diff->GetXaxis()->SetBinLabel(3.0,"k");
  Diff->GetXaxis()->SetBinLabel(4.0,"k+1sigma");
  Diff->GetXaxis()->SetBinLabel(5.0,"k+2sigma");
  Diff->Fill(1.0,avgD_m2s);
  Diff->SetBinError(1.0,D_m2sE);
  Diff->Fill(2.0,avgD_m1s);
  Diff->SetBinError(2.0,D_m1sE);
  Diff->Fill(3.0,0.0);
  Diff->SetBinError(3.0,normE);
  Diff->Fill(4.0,avgD_1s);
  Diff->SetBinError(4.0,D_1sE);
  Diff->Fill(5.0,avgD_2s);
  Diff->SetBinError(5.0,D_2sE);
  


  TCanvas *c = 0;
  std::cout<< "printing canvases"<<std::endl;
  c = new TCanvas("","",800,600);
  Diff->Draw();
  c->SaveAs("mH_differences.pdf");

  c = new TCanvas("","",800,600);
  norm->Draw();
  TF1* func1 = new TF1("func1","gaus",122,128);
  TFitResultPtr r = norm->Fit(func1,"RS");
  gStyle->SetOptFit(0111);
  c->SaveAs("MassHiggs.pdf");

  c = new TCanvas("","",800,600);
  p1s->Draw();
  TF1* func2 = new TF1("func2","gaus",122,128);
  TFitResultPtr r1 = p1s->Fit(func2,"RS");
  gStyle->SetOptFit(0111);
  c->SaveAs("MassHiggs_1k.pdf");

  c = new TCanvas("","",800,600);
  p2s->Draw();
  TF1* func3 = new TF1("func3","gaus",122,128);
  TFitResultPtr r2 = p2s->Fit(func3,"RS");
  gStyle->SetOptFit(0111);
  c->SaveAs("MassHiggs_2k.pdf");

  c = new TCanvas("","",800,600);
  m1s->Draw();
  TF1* func4  = new TF1("func4","gaus",122,128);
  TFitResultPtr r3 = m1s->Fit(func4,"RS");
  gStyle->SetOptFit(0111);
  c->SaveAs("MassHiggs_m1k.pdf");

  c = new TCanvas("","",800,600);
  m2s->Draw();
  TF1* func5= new TF1("func5","gaus",122,128);
  TFitResultPtr r4 = m2s->Fit(func5,"RS");
  gStyle->SetOptFit(0111);
  c->SaveAs("MassHiggs_m2k.pdf");
  std::cout<<"printing in general"<<std::endl;
  print(HiggsM,HiggsM_1s,HiggsM_2s,HiggsM_m1s,HiggsM_m2s,Diff_1s,Diff_2s,Diff_m1s,Diff_m2s);
  std::cout<<"done"<<std::endl;
  
}

void ShapeStudy::print(std::vector<std::tuple<double,double>> HiggsM,std::vector<std::tuple<double,double>> HiggsM_1s, std::vector<std::tuple<double,double>> HiggsM_2s, std::vector<std::tuple<double,double>> HiggsM_m1s, std::vector<std::tuple<double,double>> HiggsM_m2s, std::vector<std::tuple<double,double>> Diff_1s, std::vector<std::tuple<double,double>> Diff_2s, std::vector<std::tuple<double,double>> Diff_m1s, std::vector<std::tuple<double,double>> Diff_m2s ){
  outfile.open("ShapeStudy.txt");
  outfile << "\n \n *********************************************** Shape Study ***************************************************** \n" << std::endl;
  auto HIt = HiggsM.begin();
  auto HIt_1s = HiggsM_1s.begin();
  auto HIt_2s = HiggsM_2s.begin();
  auto HIt_m1s = HiggsM_m1s.begin();
  auto HIt_m2s = HiggsM_m2s.begin();

  auto DIt_1s = Diff_1s.begin();
  auto DIt_2s = Diff_2s.begin();
  auto DIt_m1s = Diff_m1s.begin();
  auto DIt_m2s = Diff_m2s.begin();

  outfile<<"     Fit mH             :           k+sigma          :           k+2*sigma           :           k-sigma          :          k-2*sigma    "<<std::endl;
  for(; HIt!=HiggsM.end();HIt++,HIt_1s++,HIt_2s++,HIt_m1s++,HIt_m2s++) {
    outfile<< Form(" (%.6f +/- %.6f)  ",std::get<0>(*HIt),std::get<1>(*HIt)) << Form("  (%.6f +/- %.6f)  ",std::get<0>(*HIt_1s),std::get<1>(*HIt_1s)) << Form("   (%.6f +/- %.6f)  ",std::get<0>(*HIt_2s),std::get<1>(*HIt_2s)) << Form("   (%.6f +/- %.6f)  ",std::get<0>(*HIt_m1s),std::get<1>(*HIt_m1s)) << Form("   (%.6f +/- %.6f)  ",std::get<0>(*HIt_m2s),std::get<1>(*HIt_m2s)) << std::endl;
    outfile <<"" << std::endl;
  }
  outfile<< "    Difference 1s       :              2s             :              -1s               :              -2s               "<<std::endl;
  for(; DIt_1s!=Diff_1s.end();DIt_1s++,DIt_2s++,DIt_m1s++,DIt_m2s++) {
    outfile<<Form(" (%.6f +/- %.6f) ",std::get<0>(*DIt_1s),std::get<1>(*DIt_1s)) << Form("  (%.6f +/- %.6f)  ",std::get<0>(*DIt_2s),std::get<1>(*DIt_2s)) << Form("%.6f +/ %.6f)  ",std::get<0>(*DIt_m1s),std::get<1>(*DIt_m1s)) << Form("   (%.6f +/ %.6f)  ",std::get<0>(*DIt_m2s),std::get<1>(*DIt_m2s)) <<std::endl;
  }
  outfile.close();
}




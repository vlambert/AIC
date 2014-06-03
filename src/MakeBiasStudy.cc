// ========================================================== //
//           Bias analysis H->gg background modeling          //
//        Including comparison of AIC Composite Result        //
//               Valere Lambert, Caltech 2014                 //
// ========================================================== //  
#include "MakeBiasStudy.h"
#include "assert.h"

MakeBiasStudy::MakeBiasStudy(const TString& inputFileName)
{
  if(inputFileName !=""){
    inputFile = new TFile(inputFileName);
    w = ((RooWorkspace*) inputFile->Get("cms_hgg_spin_workspace")); 
    MakeAICFits::getLabels("evtcat",&catLabels,w);
    catLabels.push_back("cat4");
  }
}
MakeBiasStudy::~MakeBiasStudy() {
}
  
void MakeBiasStudy::biasStudy(const TString& truthType, const TString& cat) {
  // Build truth model
  RooRealVar *mass = w->var("mass");
  RooAbsData *dc;
  if (cat == "cat4"){
    dc = w->data("Data_Combined");
  }
  else{
    dc = w->data("Data_Combined")->reduce(TString("evtcat==evtcat::")+cat);
  }
  SampleSize = dc->sumEntries();
  RooRealVar *nBkgTruth = new RooRealVar("TruthNBkg","",SampleSize,0,1e9);
  RooAbsPdf *truthPdf = wmodels->pdf(truthType);

  // Set truth parameters constant
  RooArgSet* vtruth = truthPdf->getVariables();
  RooFIter ittruth = vtruth->fwdIterator();
  RooAbsArg* atruth;
  while ( ( atruth=ittruth.next()) ){
    if(std::string(atruth->GetName()).compare("mass")==0) continue;
    static_cast<RooRealVar*>(atruth)->setConstant(kTRUE);
  }

  RooExtendPdf *truthExtendedPdf = new RooExtendPdf("truthExtendedPdf","",*truthPdf,*nBkgTruth);
  truthExtendedPdf->fitTo(*dc,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));

  // Require fit to converge
  int covariance = 0;
  RooFitResult* result;
  while (covariance < 2) {
    result =  truthExtendedPdf->fitTo(*dc,RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE),RooFit::Save(kTRUE));
    covariance = result->covQual();
  }

  // Determine number of truth background events within signal region
  double BiasWindow = 2.00;
  mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
  mass->setRange("FullRegion",mass->getMin(),mass->getMax());
  RooAbsReal *TruthRange = truthExtendedPdf->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass));
  RooAbsReal *TruthAll = truthExtendedPdf->createIntegral(*mass,RooFit::Range("FullRegion"),RooFit::NormSet(*mass));
  double NTruth = TruthRange->getVal() * nBkgTruth->getVal() / TruthAll->getVal();
  double NTruthE = NTruth * sqrt(1/TruthAll->getVal() + pow(TruthRange->getPropagatedError(*result)/TruthRange->getVal(),2) + 1/nBkgTruth->getVal());
  
  RooDataSet* truthbkg = truthExtendedPdf->generate(RooArgSet(*mass),int(nBkgTruth->getVal()),RooFit::Extended(true));

  std::vector<std::vector<double>> Bias(FitTypes.size());
  std::vector<std::vector<double>> BiasErr(FitTypes.size());
  std::vector<std::vector<double>> Slope(FitTypes.size());
  std::vector<std::vector<double>> SlopeNorm(FitTypes.size());

  for (int i=0;i<NToys;i++) {
    // Initialize model weights
    RooRealVar* w1 = new RooRealVar("SingExpW","SingExpW",0);
    RooRealVar* w2 = new RooRealVar("DoubExpW","DoubExpW",0);
    RooRealVar* w3 = new RooRealVar("TripExpW","TripExpW",0);
    RooRealVar* w4 = new RooRealVar("ModExpW", "ModExpW" ,0);
    RooRealVar* w5 = new RooRealVar("PolyW", "PolyW",0);
    RooRealVar* w6 = new RooRealVar("SingPowW","SingPowW",0);
    RooRealVar* w7 = new RooRealVar("DoubPowW","DoubPowW",0);
    
    auto BiasIt = Bias.begin();
    auto BiasEIt = BiasErr.begin();
    auto SlopeIt = Slope.begin();
    auto SlopeNormIt = SlopeNorm.begin();
    
    // Single Exponential Fit
    RooRealVar *Nbkg_0 = new RooRealVar("Nbkg_0","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r0 = MakeBiasStudy::makeBackgroundFits(*truthbkg, 0, Nbkg_0);
    RooAbsPdf* singExp = std::get<0>(r0);
    std::cout<< std::get<1>(r0) << "  "<< std::get<2>(r0) << "  "<<std::get<3>(r0) << "  "<< std::get<4>(r0)<<std::endl;
    singExp->SetNameTitle("singExp","Single Exponential");
    BiasEIt->push_back( ( std::get<1>(r0) - NTruth)/std::get<2>(r0) );
    BiasIt->push_back( std::get<1>(r0)-NTruth );
    if (std::get<3>(r0) != 0){
      SlopeIt->push_back( std::get<3>(r0) );
    }
    if (std::get<4>(r0) != 0){
    SlopeNormIt->push_back( std::get<4>(r0) );
    }
    BiasEIt++;
    BiasIt++;
    SlopeIt++;
    SlopeNormIt++;
    
    // Double Exponential Fit
    RooRealVar *Nbkg_1 = new RooRealVar("Nbkg_1","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r1 = MakeBiasStudy::makeBackgroundFits(*truthbkg,1, Nbkg_1);
    RooAbsPdf* doubExp = std::get<0>(r1);
    doubExp->SetNameTitle("doubExp","Double Exponential");
    BiasEIt->push_back( ( std::get<1>(r1) - NTruth)/std::get<2>(r1) );
    BiasIt->push_back( std::get<1>(r1)-NTruth );
    if (std::get<3>(r1) != 0){
      SlopeIt->push_back( std::get<3>(r1) );
    }
    if (std::get<4>(r1) != 0){
      SlopeNormIt->push_back( std::get<4>(r1) );
    }
    SlopeIt->push_back( std::get<3>(r1) );
    SlopeNormIt->push_back( std::get<4>(r1) );
    BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++;
    
    // Triple Exponential Fit
    RooRealVar *Nbkg_2 = new RooRealVar("Nbkg_2","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r2 = MakeBiasStudy::makeBackgroundFits(*truthbkg, 2, Nbkg_2);
    RooAbsPdf* tripExp = std::get<0>(r2);
    tripExp->SetNameTitle("tripExp","Triple Exponential");
    BiasEIt->push_back( ( std::get<1>(r2) - NTruth)/std::get<2>(r2) );
    BiasIt->push_back( std::get<1>(r2)-NTruth );
    if (std::get<3>(r2) != 0){
      SlopeIt->push_back( std::get<3>(r2) );
    }
    if (std::get<4>(r2) != 0){
      SlopeNormIt->push_back( std::get<4>(r2) );
    }
    BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++;
    
    // Modified Exponential Fit
    RooRealVar *Nbkg_3 = new RooRealVar("Nbkg_3","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r3  = MakeBiasStudy::makeBackgroundFits(*truthbkg, 3, Nbkg_3);
    RooAbsPdf* modExp  = std::get<0>(r3);
    modExp->SetNameTitle("modExp","Modified Exponential");
    BiasEIt->push_back( ( std::get<1>(r3) - NTruth)/std::get<2>(r3) );
    BiasIt->push_back( std::get<1>(r3)-NTruth );
    if (std::get<3>(r3) != 0){
      SlopeIt->push_back( std::get<3>(r3) );
    }
    if (std::get<4>(r3) != 0){
      SlopeNormIt->push_back( std::get<4>(r3) );
    }
    BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++;
    
    // 5th Order Polynomial Fit
    RooRealVar *Nbkg_4 = new RooRealVar("Nbkg_4","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r4  = MakeBiasStudy::makeBackgroundFits(*truthbkg, 4, Nbkg_4);
    RooAbsPdf* Poly    = std::get<0>(r4);
    Poly->SetNameTitle("Poly","5th Order Polynomial");
    BiasEIt->push_back( ( std::get<1>(r4) - NTruth)/std::get<2>(r4) );
    BiasIt->push_back( std::get<1>(r4)-NTruth );
    if (std::get<3>(r4) != 0){
      SlopeIt->push_back( std::get<3>(r4) );
    }
    if (std::get<4>(r4) != 0){
      SlopeNormIt->push_back( std::get<4>(r4) );
    }
    BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++;
    
    // Single Power Law Fit
    RooRealVar *Nbkg_5 = new RooRealVar("Nbkg_5","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r5 = MakeBiasStudy::makeBackgroundFits(*truthbkg, 5, Nbkg_5);
    RooAbsPdf* singPow = std::get<0>(r5);
    singPow->SetNameTitle("singPow","Single Power Law");
    BiasEIt->push_back( ( std::get<1>(r5) - NTruth)/std::get<2>(r5) );
    BiasIt->push_back( std::get<1>(r5)-NTruth );
    if (std::get<3>(r5) != 0){
      SlopeIt->push_back( std::get<3>(r5) );
    }
    if (std::get<4>(r5) != 0){
      SlopeNormIt->push_back( std::get<4>(r5) );
    }
    BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++;
    
    // Double Power Law Fit
    RooRealVar *Nbkg_6 = new RooRealVar("Nbkg_6","N Background Events", SampleSize,0,1e9);
    std::tuple<RooAbsPdf*,double,double,double,double> r6 = MakeBiasStudy::makeBackgroundFits(*truthbkg, 6, Nbkg_6);
    RooAbsPdf* doubPow = std::get<0>(r6);
    doubPow->SetNameTitle("doubPow","Double Power Law");  
    BiasEIt->push_back( ( std::get<1>(r6) - NTruth)/std::get<2>(r6) );
    BiasIt->push_back( std::get<1>(r6)-NTruth );
    if (std::get<3>(r6) != 0){
      SlopeIt->push_back( std::get<3>(r6) );
    }
    if (std::get<4>(r6) != 0){
      SlopeNormIt->push_back( std::get<4>(r6) );
    }
    BiasEIt++,BiasIt++,SlopeIt++,SlopeNormIt++;

    // Compare fit model AIC values
    RooArgList* ModelSet = new RooArgList(*singExp,*doubExp,*tripExp,*modExp,*Poly,*singPow,*doubPow);
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
    
    RooAbsPdf *JointModel = new RooAddPdf("Composite", "Background Model", *ModelSet, *Weights);
    // Fix all composite model parameters
    RooArgSet* vjoint = JointModel->getVariables();
    RooFIter iter = vjoint->fwdIterator();
    RooAbsArg* a;
    while ( ( a=iter.next()) ){
      if(std::string(a->GetName()).compare("mass")==0) continue;
      static_cast<RooRealVar*>(a)->setConstant(kTRUE);
    }
    
    // fit composite model to data
    RooRealVar *NBkg = new RooRealVar("NBkg","N Background Events",SampleSize,0,1e9);
    RooExtendPdf *JointModelExt = new RooExtendPdf("","",*JointModel,*NBkg);
    
    RooFitResult *composite_databkg;
    JointModelExt->fitTo(*truthbkg,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    composite_databkg = JointModelExt->fitTo(*truthbkg,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    
    // Calculate bias for composite model
    double BiasWindow = 2.00;
    mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
    mass->setRange("FullRegion",mass->getMin(),mass->getMax() );
    RooAbsReal *FitRange = JointModelExt->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass));
    RooAbsReal *FitAll = JointModelExt->createIntegral(*mass,RooFit::Range("FullRegion"),RooFit::NormSet(*mass));
    double NFit = FitRange->getVal() * NBkg->getVal()/FitAll->getVal();
    double NFitE = NFit * sqrt(1/FitAll->getVal() + pow(FitRange->getPropagatedError(*composite_databkg)/FitRange->getVal(),2) + 1/NBkg->getVal()); 

    // Calculate derivative for slope
    RooAbsReal* dmass = JointModelExt->derivative(*mass,1);
    RooPlot *dframe = mass->frame();
    dmass->plotOn(dframe,RooFit::Name("deriv"));
    TGraph *tg = (TGraph*)dframe->findObject("deriv");
    
    double SlopeCenter = tg->Eval(125)/FitAll->getVal();
    double SlopeLow = tg->Eval(123)/FitAll->getVal();
    double SlopeHigh = tg->Eval(127)/FitAll->getVal();
    double SlopeAvg = (SlopeHigh+SlopeLow)/2.0;
    
    BiasEIt->push_back( ( NFit - NTruth)/NFitE );
    BiasIt->push_back( NFit-NTruth );
    SlopeIt->push_back( SlopeCenter );
    SlopeNormIt->push_back( SlopeAvg );
    
    // clean up models
    delete composite_databkg;
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
    // don't fill if value is zero, find median of length of vector, if length is zero the fill with zero
    BiasBkg[cat][truthType].push_back( BiasIt->at( BiasIt->size()/2 ));
    BiasBkgErr[cat][truthType].push_back( fabs(BiasEIt->at( BiasEIt->size()/2 )));
    if (SlopeIt->size()/2 != 0){ 
      SlopeBkg[cat][truthType].push_back( SlopeIt->at( SlopeIt->size()/2 ));
    }
    else {
      SlopeBkg[cat][truthType].push_back(0);
    }
    if (SlopeNormIt->size()/2 !=0){
      SlopeNormBkg[cat][truthType].push_back( SlopeNormIt->at( SlopeNormIt->size()/2 ));
    }
    else{
      SlopeNormBkg[cat][truthType].push_back(0);
    }
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
   
  
std::tuple<RooAbsPdf*,double,double,double,double> MakeBiasStudy::makeBackgroundFits(RooAbsData& truthbkg, int type, RooRealVar *Nbkg) {
  // Build fit model
  RooRealVar* mass = w->var("mass"); 
  RooAbsPdf* ModelShape = MakeAICFits::getBackgroundPdf(type,mass);
  int k = MakeAICFits::Num_Params(type);
  RooExtendPdf *BkgModel = new RooExtendPdf("BKGFIT_bkgModel", "Background Model", *ModelShape, *Nbkg);

  // Require fit to converge
  int covariance = 0;
  RooFitResult *bkg_databkg;
  while (covariance < 2) {
    BkgModel->fitTo(truthbkg,RooFit::Strategy(0),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    bkg_databkg = BkgModel->fitTo(truthbkg, RooFit::Save(kTRUE), RooFit::Strategy(2),RooFit::NumCPU(NUM_CPU),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
    covariance = bkg_databkg->covQual();
  }
  Double_t bkg_databkg_Nll = bkg_databkg->minNll();

  // Calculate AIC for each model
  AIC_bkg_array[type] = 2.*(k + k*(k + 1.)/(SampleSize - (k + 1.)) + bkg_databkg_Nll);

  // Calculate bias within signal region
  double BiasWindow = 2.00;
  mass->setRange("biasRegion", mh-BiasWindow, mh+BiasWindow);
  mass->setRange("FullRegion",mass->getMin(),mass->getMax() );
  RooAbsReal *FitRange = BkgModel->createIntegral(*mass,RooFit::Range("biasRegion"),RooFit::NormSet(*mass));
  RooAbsReal *FitAll = BkgModel->createIntegral(*mass,RooFit::Range("FullRegion"),RooFit::NormSet(*mass));
  double NFit = FitRange->getVal() * Nbkg->getVal()/FitAll->getVal();
  double NFitE = sqrt(NFit);   
  //double NFitE = NFit *sqrt(1/FitAll->getVal() + pow(FitRange->getPropagatedError(*bkg_databkg)/FitRange->getVal(),2) + 1/Nbkg->getVal()); 

  // Determine derivative for slope
  RooAbsReal* dmass = BkgModel->derivative(*mass,1);
  RooPlot *dframe = mass->frame();
  dmass->plotOn(dframe,RooFit::Name("deriv"));
  TGraph *tg = (TGraph*)dframe->findObject("deriv");
  
  double SlopeCenter = (tg->Eval(125))/FitAll->getVal();
  double SlopeLow = (tg->Eval(123))/FitAll->getVal();
  double SlopeHigh = (tg->Eval(127))/FitAll->getVal();
  double SlopeAvg = (SlopeHigh+SlopeLow)/2.0;


  // clean up
  delete FitRange;
  delete FitAll;
  delete bkg_databkg;
  bkg_databkg_Nll = 0.;

  return std::tuple<RooAbsPdf*,double,double,double,double>(BkgModel, NFit, NFitE, SlopeCenter, SlopeAvg);
}

void MakeBiasStudy::run(){
  for (auto cat: catLabels) {
    std::cout<<cat<<std::endl;
  }
  for (auto cat: catLabels) {
    modelFile = TFile::Open(TString("Workspaces/HggAIC_workspace_")+cat+TString(".root"));
    wmodels = ((RooWorkspace*) modelFile->Get("AIC"));
    for (auto truthIt : FitTypes) {
      std::cout<<truthIt<<std::endl;
      biasStudy(truthIt,cat);
    }
    modelFile->Close();
  }

  print();
  inputFile->Close();
  std::cout<< "========== Finished Running Bias Study =========="<<std::endl;
}

void MakeBiasStudy::print(){
  outfile.open("BiasOutput.txt");
  outfile << "\n \n ************************* BKG ERROR *************************** \n" << std::endl;
  outfile << "  ======================= Bias =======================  " << std::endl;
  printFormatted(BiasBkg);
  outfile << "  ======================= Bias Error =======================  " <<std::endl;
  printFormatted(BiasBkgErr);
  outfile << "  ======================= Slope Center =======================  " <<std::endl;
  printFormatted(SlopeBkg);
  outfile <<"  =======================  Slope Average =======================  " << std::endl;
  printFormatted(SlopeNormBkg);
  outfile.close();
}

void MakeBiasStudy::printFormatted(const biasList& list) {
  // Print out the bias results for each truth and fit model pair
  for(auto truthIt: FitTypes) {
    outfile <<"Truth Type:  " <<  fitNameMap[truthIt] << std::endl;
    outfile << "Fit Type:          " <<"   SingExp   " << "   DoubExp   " <<  "   TripExp    "  << "   ModExp   " << "   Poly    " << "  SingPow  " << "  DoubPow   " << "   Composite  "  <<std::endl;
    for (auto catIt: catLabels){
      outfile<< catIt << "         ";
      for (auto BiasIt : list.at(catIt).at(truthIt) ) {
	outfile << Form(" & %10.7f ",BiasIt) << "   " ;
      }
      outfile << " \\\\ \\hline" <<std::endl;
    }
  }
}

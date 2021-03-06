#ifndef BootStrapping_h
#define BootStrapping_h

//! class to make the bootstrapping check

/*!
  

Author: Valere Lambert (Caltech)
Date: Nov 2013
*/

#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TMath.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooGlobalFunc.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include "RooStats/SPlot.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooBernstein.h"
#include "RooLinkedListIter.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooProdPdf.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"

#include "MakeAICFits.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <assert.h>

//#include <omp.h> //threading!


class BootStrapping{
 public:
  BootStrapping(const TString& inputFileName);  
  
  RooAbsPdf* makeFits(int type, RooDataSet *Sample, RooRealVar *Nbkg);
  void run();
  void print(int toy);
  enum BkgFitType{kSingleExp, kDoubleExp, kTripleExp, kModifiedExp, kPoly, kPow, kDoublePow};

 private:
  RooWorkspace *w;
  RooAbsData* dc;
  RooRealVar *mass;
  TFile* inputFile = 0;
  TFile* outputFile = 0;

  std::ofstream outfile;
  std::vector<MakeAICFits::BkgFitType> FitTypes = {MakeAICFits::kSingleExp,MakeAICFits::kDoubleExp,MakeAICFits::kTripleExp,MakeAICFits::kModifiedExp,MakeAICFits::kPoly,MakeAICFits::kPow,MakeAICFits::kDoublePow};
  
  int NData;
  int NToys = 200;
  int nModels = 7;
  Double_t LogLikelihood[100];
  Double_t AIC_bkg_array[100];
  Double_t DeltaI[100];
  Double_t AICweights[100];
  Double_t Chi2[100];
  Double_t LL[100];
  double NSigRange[100];
  double NSigRangeE[100];
  double NSig[100];
  std::vector<RooAbsPdf*> CompositeM;

  std::map<MakeAICFits::BkgFitType,TString> fitNameMap = {{MakeAICFits::kPow, "Single Power Law"}, {MakeAICFits::kDoublePow, "Double Power Law"}, {MakeAICFits::kSingleExp,"Single Exponential"},{MakeAICFits::kDoubleExp, "Double Exponential"},{MakeAICFits::kTripleExp,"Triple Exponential"},{MakeAICFits::kModifiedExp,"Modified Exponential"},{MakeAICFits::kPoly,"5th Order Polynomial"}};  
};

#endif

#ifndef MakeAICFits_h
#define MakeAICFits_h

//! class to make all the AIC Fits

/*!
  

Author: Valere Lambert (Caltech)
Date: Nov 2013
*/

#include <TPaveText.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TMath.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooDataHist.h>
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
//#include "RooIntegralMorph.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"

//#include "RooGaussianCorr.hh"

//#include <HggOutputReader2.h>
//#include <dataSetInfo.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <assert.h>

//#include <omp.h> //threading!


class MakeAICFits{
 public:
  MakeAICFits(const TString& inputFileName);
  
  static int Num_Params(int type);
  static RooAbsPdf* getBackgroundPdf(int type, RooRealVar* mass);
  static void getLabels(const char *varName, std::vector<TString> *lblVec,RooWorkspace *w);
  RooAbsPdf* makeBackgroundFits(int type, RooRealVar* Nbkg);
  void runCategory(const TString& catTag);
  void run();
  void print();

  enum BkgFitType{kSingleExp, kDoubleExp, kTripleExp, kModifiedExp, kPoly, kPow, kDoublePow};
 protected:
  std::vector<TString> catLabels;
  RooWorkspace *w;
  TFile* inputFile = 0;
  TFile* outputFile = 0;
  std::ofstream outfile;

  RooAbsData* dc;
  RooRealVar *mass;
  int nModels=7;
  int SampleSize;
  Double_t LogLikelihood[100];
  Double_t AIC_bkg_array[100];
  Double_t DeltaI[100];
  Double_t AICweights[100];
  Int_t avgcov[100];
};

#endif

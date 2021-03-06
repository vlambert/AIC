#ifndef MakeAICFits_h
#define MakeAICFits_h

//! class to make all the AIC Fits

/*!
  

Author: Valere Lambert (Caltech)
Date: Nov 2013
*/

#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
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
//#include "RooIntegralMorph.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"

//#include "RooGaussianCorr.hh"

//#include <HggOutputReader2.h>
//#include <dataSetInfo.h>

#include <iostream>
#include <map>
#include <vector>
#include <assert.h>

//#include <omp.h> //threading!


class MakeAICFits{
 public:
  MakeAICFits();
  
  //MakeAICFits(const TString& inputFileName, const TString& outputFileName);  
  int Num_Params(int type);
  static RooAbsPdf* getBackgroundPdf(int type, RooRealVar* mass);

  enum BkgFitType{kSingleExp, kDoubleExp, kTripleExp, kModifiedExp, kPoly, kPow, kDoublePow};
 protected:
  BkgFitType fitType;
};

#endif

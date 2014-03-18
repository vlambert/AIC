#ifndef MakeBiasStudy_h
#define MakeBiasStudy_h
/* class to do simple hypothesis testing over a mass range

Author: Valere Lambert (Caltech)
Date: Feb 2014
*/

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooLinkedListIter.h"
#include "RooAddPdf.h"

#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"

#include "MakeAICFits.h"
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

#define NUM_CPU 1

class MakeBiasStudy {
 public:
  MakeBiasStudy(const TString& inputFileName);
  virtual ~MakeBiasStudy();
  
  void biasStudy(const TString& truthType, const TString& cat);
  std::tuple<double,double,double,double> getNFit(RooAbsData& truthbkg, const TString& modelType);

  void run();
  void print();
  
  typedef std::map<TString,std::map<TString,std::vector<double> > > biasList;
 private:
  RooWorkspace *w;
  RooWorkspace *wmodels;

  TFile* inputFile = 0;
  TFile* modelFile = 0;
  std::ofstream outfile;
  std::vector<TString> catLabels;
  
  std::vector<TString> FitTypes = {"SingleExp","DoubleExp","TripleExp","ModifiedExp","Polynomial","Power","DoublePower","Composite"};
  //std::vector<TString> FitTypes = {"SingleExp","Power","Polynomial"};

  int SampleSize;
  biasList BiasBkg;
  biasList BiasBkgErr;
  biasList SlopeBkg;
  biasList SlopeNormBkg;
  float mh=125;
  int NToys = 100;

  std::map<TString,TString> fitNameMap = {{"SingleExp","Single Exponential"},{"DoubleExp","Double Exponential"},{"TripleExp","Triple Exponential"},{"ModifiedExp","Modified Exponential"},{"Polynomial","5th Order Polynomial"},{"Power","Single Power Law"},{"DoublePower","Double Power Law"},{"Composite","Composite Model"}}; 
  void printFormatted(const biasList& list);
};
#endif
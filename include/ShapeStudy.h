#ifndef ShapeStudy_h
#define ShapeStudy_h
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
#include "TF1.h"
#include "TH1.h"
#include <TFitResult.h>
#include <TStyle.h>
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

class ShapeStudy {
 public:
  ShapeStudy(const TString& inputFileName);
  virtual ~ShapeStudy();
  
  std::tuple<double,double> getBiasedMass(double k_value, RooAbsData& toydata);
  
  void run();
  // void print(std::vector<std::tuple<double,double>> HiggsM, std::vector<std::tuple<double,double>> HiggsM_1s, std::vector<std::tuple<double,double>> HiggsM_2s, std::vector<std::tuple<double,double>> HiggsM_m1s, std::vector<std::tuple<double,double>> HiggsM_m2s, std::vector<std::tuple<double,double>> Diff_1s, std::vector<std::tuple<double,double>> Diff_2s, std::vector<std::tuple<double,double>> Diff_m1s, std::vector<std::tuple<double,double>> Diff_m2s);

  void print(std::vector<std::tuple<double,double>> HiggsM,std::vector<std::tuple<double,double>> HiggsM_1s, std::vector<std::tuple<double,double>> HiggsM_2s, std::vector<std::tuple<double,double>> HiggsM_m1s, std::vector<std::tuple<double,double>> HiggsM_m2s, std::vector<std::tuple<double,double>> Diff_1s, std::vector<std::tuple<double,double>> Diff_2s, std::vector<std::tuple<double,double>> Diff_m1s, std::vector<std::tuple<double,double>> Diff_m2s );
  
 private:
  RooWorkspace *w;

  TFile* inputFile = 0;
  std::ofstream outfile;

  int SampleSize;
  RooRealVar *mass;
  RooAbsData *dc;
  RooRealVar *Nbkg;
  RooRealVar *nBkgSB;
  RooRealVar *nSigSB;
  float mh=125;
  int NToys = 100;


};
#endif

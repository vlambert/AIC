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

#define NUM_CPU 1

class MakeBiasStudy{
 public:
  MakeBiasStudy();

  string Category(int type);
  //void biasStudy();
  enum BkgFitType{kSingleExp, kDoubleExp, kTripleExp, kModifiedExp, kPoly, kPow, kDoublePow};

 private:
  float mh=125;


};

#endif

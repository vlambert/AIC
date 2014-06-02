#ifndef __CINT__
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#endif

#include <math.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <vector>
#include <map>



Int_t bkgColors[7] = { kRed, kBlue, kGreen+2, kMagenta, kCyan, kYellow, kBlack};
string modelLegendLabels[7] = { "Single Exp.",
				"Double Exp.",
				"Triple Exp.",
				"Modified Exp.",
				"Polynomial",
				"Single Power",
				"Double Power"};

string LatexLabels[7] = { "Single Exp.",
			  "Double Exp.",
			  "Triple Exp.",
			  "Modified Exp.",
			  "Polynomial",
			  "Single Power",
			  "Double Power"};
  
void BiasDiff () {
  vector<TH1F*> AIC;
  
  for (UInt_t i = 1 ; i <= 7; i++){
    AIC.push_back( new TH1F( Form("AIC_%d",i), " ", 9, 0,10 ));
    AIC[i-1]->GetXaxis()->SetBinLabel(1,"category 0");
    AIC[i-1]->GetXaxis()->SetBinLabel(2,"");
    AIC[i-1]->GetXaxis()->SetBinLabel(3,"category 1");
    AIC[i-1]->GetXaxis()->SetBinLabel(4,"");
    AIC[i-1]->GetXaxis()->SetBinLabel(5,"category 2");
    AIC[i-1]->GetXaxis()->SetBinLabel(6,"");
    AIC[i-1]->GetXaxis()->SetBinLabel(7,"category 3");
    AIC[i-1]->GetXaxis()->SetBinLabel(8,"");
    AIC[i-1]->GetXaxis()->SetBinLabel(9,"Inclusive");
    AIC[i-1]->SetFillColor(bkgColors[i-1]);
    AIC[i-1]->SetLineColor(bkgColors[i-1]);
  }
  
  AIC[0]->Fill(1.0,0.02);
  AIC[0]->Fill(3.0,0.3089);
  AIC[0]->Fill(5.0,0.62);
  AIC[0]->Fill(7.0,0.36);
  AIC[0]->Fill(9.0,0.03);
  
  AIC[1]->Fill(1.0,0.08);
  AIC[1]->Fill(3.0,0.21);
  AIC[1]->Fill(5.0,0.08);
  AIC[1]->Fill(7.0,0.16);
  AIC[1]->Fill(9.0,0.33);
  
  AIC[2]->Fill(1.0,0.01);
  AIC[2]->Fill(3.0,0.02);
  AIC[2]->Fill(5.0,0.01);
  AIC[2]->Fill(7.0,0.02);
  AIC[2]->Fill(9.0,0.04);
  
  AIC[3]->Fill(1.0,0.20);
  AIC[3]->Fill(3.0,0.38);
  AIC[3]->Fill(5.0,0.24);
  AIC[3]->Fill(7.0,0.451);
  AIC[3]->Fill(9.0,0.60);
  
  AIC[4]->Fill(1.0,0.03);
  AIC[4]->Fill(3.0,0.08);
  AIC[4]->Fill(5.0,0.017);
  AIC[4]->Fill(7.0,0.001);
  AIC[4]->Fill(9.0,0.0001);
  
  AIC[5]->Fill(1.0,0.58);
  AIC[5]->Fill(3.0,0.001);
  AIC[5]->Fill(5.0,0.03);
  AIC[5]->Fill(7.0,0.007);
  AIC[5]->Fill(9.0,0.000003);
  
  AIC[6]->Fill(1.0,0.08);
  AIC[6]->Fill(3.0,0.0001);
  AIC[6]->Fill(5.0,0.003);
  AIC[6]->Fill(7.0,0.001);
  AIC[6]->Fill(9.0,0.0000004);
  
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  
  // ===================================
  //           AIC
  // ===================================
  
  cv = new TCanvas("cv","cv",800,600);
  legend = new TLegend(0.64,0.64,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  
  THStack *stackAIC = new THStack();
  for (Int_t i = AIC.size()-1; i>=0; i--) {
    if (AIC[i]->Integral()>0) {
      stackAIC->Add(AIC[i]);
      legend->AddEntry(AIC[i],modelLegendLabels[i].c_str(),"F");
    }
  }
  stackAIC->Draw();
  stackAIC->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackAIC->GetHists()->At(0)))->GetXaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs("AICvalues.pdf");
  
  // =======================================================

  //double x[7] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
  double x[4] = {1.0,2.0,3.0,4.0};
  double comps[4],singEs[4],doubEs[4],tripEs[4],modEs[4],polys[4],pows[4],dpows[4];

 //Composite  
  //comps[0]  = ( (0.0042-0.0033)/0.0033+(0.0007-0.0036)/0.0036+(0.0180-0.012)/0.012+(0.0090-0.0023)/0.0023+ (0.0033-0.0031)/0.0031 ) /5;
  comps[0] = ( (0.0045-0.0043)/0.0043 + (0.0035-0.0018)/0.0018 + (0.0074-0.0084)/0.0084 + (0.0053-0.0065)/0.0065 + (0.0016-0.0029)/0.0029 ) /5;

  //comps[2] = ( (0.0036-0.0049)/0.0049 + (0.0075-0.0029)/0.0029 + (0.0125-0.0126)/0.0126 + (0.0041-0.0031)/0.0031 + (0.0053-0.0027)/0.0027 )/ 5;

  comps[1] = ( (0.0088-0.0122)/0.0122  + (0.0072-0.0064)/0.0064 + (0.0064-0.0057)/0.0057 ) / 3;

  comps[2] = ( (0.0142-0.0254)/0.0254 + (0.0118-0.0120)/0.0120 + (0.0024-0.0105)/0.0105 + (0.0057-0.0102)/0.0102 ) / 4;

  //comps[5] = ((0.0047 - 0.0043)/0.0043 + (0.0091 - 0.0083)/0.0083 + (0.0062 - 0.0071)/0.0071 + (0.0129-0.0120)/0.0120 + (0.0035 - 0.0014)/0.0014  ) / 5;

  comps[3] = ( (0.0040 - 0.0059)/0.0059 + (0.004 - 0.0009)/0.0009 + (0.0088 - 0.0046)/0.0046 + (0.0026 - 0.0033)/0.0033 + (0.0014 - 0.0022)/0.0022 ) / 5;


 //singE  

  singEs[0] = ( (0.0221-0.0043)/0.0043 + (0.0079-0.0018)/0.0018 + (0.0084-0.0084)/0.0084 + (0.0033-0.0065)/0.0065 + (0.0082-0.0029)/0.0029 ) /5;
 
  singEs[1] = ( (0.0050-0.0122)/0.0122  + (0.0084-0.0064)/0.0064 + (0.0114-0.0057)/0.0057 ) / 3;

  singEs[2] = ( (0.0196-0.0254)/0.0254 + (0.0024-0.0120)/0.0120 + (0.0008-0.0105)/0.0105 + (0.0003-0.0102)/0.0102 ) / 4;

  //singEs[3] = ( (0.0126 - 0.0059)/0.0059 + (0.0174 - 0.0009)/0.0009 + (0.0129 - 0.0046)/0.0046 + (0.0193 - 0.0033)/0.0033 + (0.0181 - 0.0022)/0.0022 ) / 5;
  singEs[3] = 3.0;
 //doubE  

  doubEs[0] = 0;
 
  doubEs[1] =   ( (0.0079-0.0122)/0.0122  + (0.0084-0.0064)/0.0064 + (0.0073-0.0057)/0.0057 ) / 3;

  doubEs[2] =  ( (0.0069-0.0254)/0.0254 + (0.0002-0.0120)/0.0120 + (0.0016-0.0105)/0.0105 + (0.0093-0.0102)/0.0102 ) / 4;

  doubEs[3] =  ( (0.0035 - 0.0059)/0.0059 + (0.004 - 0.0009)/0.0009 + (0.0097 - 0.0046)/0.0046 + (0.0041 - 0.0033)/0.0033 + (0.0052 - 0.0022)/0.0022 ) / 5;

  //tripE  

  tripEs[0] = ( (0.0045-0.0043)/0.0043 + (0.0035-0.0018)/0.0018 + (0.0074-0.0084)/0.0084 + (0.0053-0.0065)/0.0065 + (0.0016-0.0029)/0.0029 ) /5;
 
  tripEs[1] = ( (0.0079-0.0122)/0.0122  + (0.0084-0.0064)/0.0064 + (0.0073-0.0057)/0.0057 ) / 3;

  tripEs[2] = ( (0.0069-0.0254)/0.0254 + (0.0004-0.0120)/0.0120 + (0.0017-0.0105)/0.0105 + (0.0027-0.0102)/0.0102 ) / 4;

  tripEs[3] = ( (0.0039 - 0.0059)/0.0059 + (0.0043 - 0.0009)/0.0009 + (0.0092 - 0.0046)/0.0046 + (0.0015 - 0.0033)/0.0033 + (0.0039 - 0.0022)/0.0022 ) / 5;

//modE

  modEs[0] = ( (0.0076-0.0043)/0.0043 + (0.0022-0.0018)/0.0018 + (0.0022-0.0084)/0.0084 + (0.0071-0.0065)/0.0065 + (0.0008-0.0029)/0.0029 ) /5;

  modEs[1] = 0;
 
  modEs[2] = ( (0.0051-0.0254)/0.0254 + (0.001-0.0120)/0.0120 + (0.0017-0.0105)/0.0105 + (0.0041-0.0102)/0.0102 ) / 4;

  modEs[3] =  ( (0.0026 - 0.0059)/0.0059 + (0.0037 - 0.0009)/0.0009 + (0.0088 - 0.0046)/0.0046 + (0.0072 - 0.0033)/0.0033 + (0.0049 - 0.0022)/0.0022 ) / 5;


//poly  ****

  polys[0] = ( (0.0094-0.0043)/0.0043 + (0.0041-0.0018)/0.0018 + (0.0178-0.0084)/0.0084 + (0.0019-0.0065)/0.0065 + (0.0036-0.0029)/0.0029 ) /5;
 
  polys[1] = ( (0.0088-0.0122)/0.0122  + (0.0050-0.0064)/0.0064 + (0.0544-0.0057)/0.0057 ) / 3;

  polys[2] = 0;

  polys[3] = ( (0.0111 - 0.0059)/0.0059 + (0.0070 - 0.0009)/0.0009 + (0.0091 - 0.0046)/0.0046 + (0.0055 - 0.0033)/0.0033 + (0.0099 - 0.0022)/0.0022 ) / 5;

  //pow  

  pows[0] = ( (0.0042-0.0043)/0.0043 + (0.0093-0.0018)/0.0018 + (0.026-0.0084)/0.0084 + (0.0183-0.0065)/0.0065 + (0.0082-0.0029)/0.0029 ) /5;
 
  pows[1] = ( (0.0022-0.0122)/0.0122  + (0.0031-0.0064)/0.0064 + (0.0052-0.0057)/0.0057 ) / 3;

  pows[2] = ( (0.0018-0.0254)/0.0254 + (0.0153-0.0120)/0.0120 + (0.0194-0.0105)/0.0105 + (0.0149-0.0102)/0.0102 ) / 4;

  pows[3] = ( (0.0058 - 0.0059)/0.0059 + (0.0009 - 0.0009)/0.0009 + (0.0045 - 0.0046)/0.0046 + (0.0047 - 0.0033)/0.0033 + (0.0021 - 0.0022)/0.0022 ) / 5;

  //dpow   ****

  dpows[0] = ( (0.0005-0.0043)/0.0043 + (0.0093-0.0018)/0.0018 + (0.0261-0.0084)/0.0084 + (0.0183-0.0065)/0.0065 + (0.0082-0.0029)/0.0029 ) /5;
 
  dpows[1] = ( (0.0022-0.0122)/0.0122  + (0.0030-0.0064)/0.0064 + (0.0052-0.0057)/0.0057 ) / 3;

  dpows[2] = ( (0.0018-0.0254)/0.0254 + (0.0153-0.0120)/0.0120 + (0.0194-0.0105)/0.0105 + (0.0097-0.0102)/0.0102 ) / 4;

  dpows[3] = 0;

  TMultiGraph *mg = new TMultiGraph();

  TGraph *comp = new TGraph(4, x,comps);
  TGraph *singE = new TGraph(4,x,singEs);
  TGraph *doubE = new TGraph(4,x,doubEs);
  TGraph *tripE = new TGraph(4,x,tripEs);
  TGraph *modE = new TGraph(4,x,modEs);
  TGraph *poly = new TGraph(4,x,polys);
  TGraph *pow = new TGraph(4,x,pows);
  TGraph *dpow = new TGraph(4,x,dpows);
  legend = new TLegend(0.64,0.64,0.92,0.94);
  
  //comp->GetXaxis()->SetBinLabel(1,"Single Exp.");
  dpow->GetXaxis()->SetBinLabel(1,"Double Exp.");
  //comp->GetXaxis()->SetBinLabel(3,"Triple Exp.");
  dpow->GetXaxis()->SetBinLabel(2,"Modified Exp.");
  dpow->GetXaxis()->SetBinLabel(3,"Polynomial");
  //comp->GetXaxis()->SetBinLabel(6,"Single Power");
  dpow->GetXaxis()->SetBinLabel(4,"Double Power");
  singE->SetMarkerColor(bkgColors[0]);
  doubE->SetMarkerColor(bkgColors[1]);
  tripE->SetMarkerColor(bkgColors[2]);
  modE->SetMarkerColor(bkgColors[3]);
  poly->SetMarkerColor(bkgColors[4]);
  pow->SetMarkerColor(bkgColors[5]);
  dpow->SetMarkerColor(bkgColors[6]);
  comp->SetMarkerColor(28);

  singE->SetLineColor(bkgColors[0]);
  doubE->SetLineColor(bkgColors[1]);
  tripE->SetLineColor(bkgColors[2]);
  modE->SetLineColor(bkgColors[3]);
  poly->SetLineColor(bkgColors[4]);
  pow->SetLineColor(bkgColors[5]);
  dpow->SetLineColor(bkgColors[6]);
  comp->SetLineColor(28);

  legend->AddEntry(comp,"Composite","p");
  legend->AddEntry(singE,"Single Exponential","p");
  legend->AddEntry(doubE,"Double Exponential","p");
  legend->AddEntry(tripE,"Triple Exponential","p");
  legend->AddEntry(modE,"Modified Exponential","p");
  legend->AddEntry(poly,"Polynomial","p");
  legend->AddEntry(pow,"Power Law","p");
  legend->AddEntry(dpow,"Double Power Law","p");

  mg->Add(comp);
  mg->Add(singE);
  mg->Add(doubE);
  mg->Add(tripE);
  mg->Add(modE);
  mg->Add(poly);
  mg->Add(pow);
  mg->Add(dpow);





  //mg->GetYaxis()->SetRangeUser(-1,3);



  cv = new TCanvas("cv","cv",800,600);
  mg->Draw("apl");
  legend->Draw();
  cv->Update();
  cv->SaveAs("BiasAverages.pdf");
}

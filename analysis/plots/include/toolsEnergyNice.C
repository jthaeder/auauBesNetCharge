#include <fstream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TPad.h"
#include "TDirectory.h"
#include "TColor.h"
#include "TMath.h"
#include "TLatex.h"
#include "TGAxis.h"
#include "TLine.h"
#include "TGraphErrors.h"

TObjArray canA;

TCanvas *can;
TPad *pad;

TLegend *legExp;
TLegend *legTheo;

// -----------------------------------------------------------

      enum    namesType          {kNetQ, kNetK, kNetP}; 

const Int_t   nNames           = 3;
const Char_t *aNames[3]        = {"Net-Charge", "Net-Kaon", "Net-Proton"};

const Char_t *aNamesPt[3]      = { "0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0, |#eta| < 0.5",
				   "0.2 < #it{p}_{T} (GeV/#it{c}) < 1.6, |#it{y}| < 0.5",
				   "0.4 < #it{p}_{T} (GeV/#it{c}) < 2.0, |#it{y}| < 0.5" };  

// -----------------------------------------------------------

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma", "#kappa #sigma^{2}"};
const Char_t* aMomentsTitle2[] = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma/Skellam", "#kappa #sigma^{2}"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV"};
const Char_t* aMoments2[]      = {"C1","C2","C3","C4","VM","SDSk","KV"};
const Int_t   nMoments         = 7;

const Char_t* aDataSetsTitle[] = {"corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}"};
const Char_t* aDataSets[]      = {"twoeff_11"};
const int     nDataSets        = 1;

      Int_t   nCent            = 9;
const Char_t* cent[9]          = {"0005",  "0510",   "1020",   "2030",   "3040",   "4050",   "5060",   "6070",   "7080"};
const Char_t* cent1[9]         = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%"};

// -----------------------------------------------------------
// J. Cleymans et al., PRC 73, 034905 (2006).
//    mub(sqrt(s)) = 1.308 /￼(1+0.27× sqrt(s))
//    Tc(sqrt(s))  = 0.166−0.139 × mub^2 −0.053 mub^4
// !  Unit : √s, μB, T (GeV)
// -- energies

const Int_t   nEnergies        = 8;
const Char_t *energies[]       = {  "7",   "11",   "14",   "19",   "27",   "39",   "62",  "200"};
const Char_t *exactEnergies[]  = {"7.7", "11.5", "14.5", "19.6", "27.0", "39.0", "62.4", "200"};

const Double_t snn[]           = {7.7, 11.5, 14.5, 19.6,  27,  39, 62.4, 200};
const Double_t mub[]           = {422,  316,  266,  206, 156, 112,   73,  24};
const Double_t tc[]            = {140,  152,  156,  160, 163, 164,  165, 166};

// -----------------------------------------------------------

      Color_t aColors[]        = {kRed+2, kGreen+2, kOrange+9, kYellow, kAzure, kAzure, kCyan+2, kBlue+1, kBlack, kRed+2};

      Int_t   aMarkers[]       = {29, 25, 25, 25, 26, 28, 27, 30, 24, 24};

      Float_t aMinY[7]         = { -1,  -7, -3, -2200, -10, -9, -20};
      Float_t aMaxY[7]         = {  40, 250, 110, 800, 160, 11,  12};

      // Float_t aMinY[7]         = { 0, 0, -10, -2300, 5.5, -0.02, -12 };
      // Float_t aMaxY[7]         = { 40, 250, 90, 900, 9.8, 0.48, 6};

      Float_t aMinX            = 6;
      Float_t aMaxX            = 250;

      Float_t aMinRatioY[]     = {0.27, 0.27, 0.27, -1.8};
      Float_t aMaxRatioY[]     = {0.9, 0.9, 0.9,  2.6};

// -----------------------------------------------------------

TList lGraphStat;
TList lGraphSys; 
TList lGraphPoisson;
TList lGraphUrqmd; 
TList lGraph14; 

// -----------------------------------------------------------

// ______________________________________________________________________________________
void PrepareGraph(TGraphErrors* g) {
  // -- Prepare Graph

  if (!g)
    return;  
  
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");
  g->GetXaxis()->SetNdivisions(6, 5, 0);
  g->GetYaxis()->SetNdivisions(6, 5, 0);
  g->GetXaxis()->SetLabelSize(0.08);
  g->GetYaxis()->SetLabelSize(0.08);
  g->GetXaxis()->SetLabelOffset(0.008);
  g->GetXaxis()->SetNoExponent(kTRUE);
}

// ______________________________________________________________________________________
void PrepareGraphCumulants(TGraphErrors* g) {
  // -- Prepare Graph

  if (!g)
    return;

  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");
  g->GetXaxis()->SetNdivisions(6, 5, 0);
  g->GetYaxis()->SetNdivisions(6, 5, 0);
  g->GetXaxis()->SetLabelSize(0.062);
  g->GetYaxis()->SetLabelSize(0.062);
  g->GetXaxis()->SetLabelOffset(0.008);
  g->GetXaxis()->SetNoExponent(kTRUE);
}

// ______________________________________________________________________________________
void ConfigGraph(TGraphErrors* g, Int_t yLimitIdx, Int_t idxMarker, Int_t is = 0) {
  // -- Config Graph
  //    is = 1 :  URQMD
  //    is = 2 :  Sys
  //    is = 3 :  Poisson

  if (!g)
    return;
  
  g->SetMinimum(aMinY[yLimitIdx]);
  g->SetMaximum(aMaxY[yLimitIdx]);
  g->GetXaxis()->SetLimits(aMinX, aMaxX);  
  g->GetXaxis()->SetMoreLogLabels(kTRUE);

  g->SetMarkerStyle(aMarkers[idxMarker]);
  g->SetMarkerColor(aColors[idxMarker]);
  g->SetMarkerSize(1.4);
  if (aMarkers[idxMarker] == 29 || aMarkers[idxMarker] == 30)
    g->SetMarkerSize(1.8);

  g->SetLineColor(aColors[idxMarker]);
  if (is == 1) { // Poisson
    g->SetLineStyle(2);
    g->SetLineWidth(2);
  }
  if (is == 2) { // SYS
    g->SetFillColorAlpha(aColors[idxMarker], 0.35);
    g->SetFillColor(aColors[idxMarker]);
    g->SetFillStyle(3344);
  }
  if (is == 3) {  // UrQMD
    g->SetLineColor(kAzure);
    g->SetLineWidth(2);
    g->SetLineStyle(2);

    g->SetFillColorAlpha(kAzure, 0.35);
    g->SetFillColor(kAzure);
    g->SetFillStyle(3344);
  }
  if (is == 4) {
    g->SetMinimum(aMinRatioY[yLimitIdx]);
    g->SetMaximum(aMaxRatioY[yLimitIdx]);
  }

}

// ______________________________________________________________________________________
void ShiftGraphX(TGraphErrors* g, Double_t shift) {
  // -- shift datapoints

  for (Int_t idx = 0 ; idx < g->GetN(); ++idx) {
    Double_t x, y;
    g->GetPoint(idx, x, y);
    g->SetPoint(idx, x+shift, y);
  }

  PrepareGraph(g);
}

// ______________________________________________________________________________________
void DrawSet(TGraphErrors *gStat, TGraphErrors *gSys, TGraphErrors *gUrqmd, TGraphErrors *gPoisson, 
	     Int_t idxMoment, Int_t idxCent, TGraphErrors* g14 = NULL) {
  // -- Draw Set 

  PrepareGraph(gStat);
  PrepareGraph(gSys);
  PrepareGraph(gPoisson);
  PrepareGraph(gUrqmd);

  ConfigGraph(gStat,    idxMoment, idxCent);
  ConfigGraph(gSys,     idxMoment, idxCent, 2);
  ConfigGraph(gUrqmd,   idxMoment, idxCent, 3);
  ConfigGraph(gPoisson, idxMoment, idxCent, 1);
  
  if (g14) {
    PrepareGraph(g14);
    ConfigGraph(g14, idxMoment, idxCent);
    if (idxCent == 0)
      g14->SetMarkerStyle(29);
    else
      g14->SetMarkerStyle(20);
  }

  // -- draw box
  if (idxCent == 0) 
    gStat->Draw("AZP");

  // -- draw urqmd
  if (idxCent == 0) {
    if (idxMoment == 4)
      gUrqmd->Draw("L,same");
    gUrqmd->Draw("E3,same");
  }
  
  // -- draw poisson 
  if (idxMoment == 5 || idxMoment == 6) {
    if (idxCent == 0) {
      TLine *line1 = new TLine(aMinX, 1, aMaxX, 1);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(2);
      line1->SetLineWidth(2);
      line1->Draw();
    }
  }
  else if (idxMoment == 4)
    gPoisson->Draw("L,SAME");
  
  // -- draw datapoints
  //gSys->Draw("B2,SAME");
  gStat->Draw("ZP,SAME");
  gSys->Draw("[],SAME");

  if (g14)
    g14->Draw("ZP,SAME");

  if (idxMoment == 4) {
    legExp->AddEntry(gStat, Form("%s", cent1[idxCent]), "pl");
    legTheo->AddEntry(gPoisson, Form("%s Poisson", cent1[idxCent]), "l");
  }
}

// ______________________________________________________________________________________
TPad* SetupCanvas(const Char_t* canName, const Char_t *canTitle) {
  // -- setup canvas and pad
  
  canA.Add(new TCanvas(Form("can%s", canName), canTitle, 0, 0 , 420, 700));
  can = static_cast<TCanvas*>(canA.Last());
  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetBorderSize(0.0);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->cd();
  
  pad = new TPad("pad", "pad",0.05, 0.06, 0.99, 0.98);
  pad->SetBorderMode(0);
  pad->SetFillColor(0);
  pad->Draw();
  pad->cd();
  pad->Divide(1, 3, 0., 0., 0);

  can->cd();

  TLatex *texb_5 = new TLatex(0.45, 0.03, "#sqrt{#it{s}_{NN}} (GeV)");
  texb_5->SetTextSize(0.04);
  texb_5->Draw("same");
  
  TLatex *texb_6a = new TLatex(0.05,0.8, aMomentsTitle[4]);
  texb_6a->SetTextSize(0.04);
  texb_6a->SetTextAngle(90);
  texb_6a->Draw("same");
  
  TLatex *texb_6b = new TLatex(0.05,0.45, aMomentsTitle2[5]);
  texb_6b->SetTextSize(0.04);
  texb_6b->SetTextAngle(90);
  texb_6b->Draw("same");
  
  TLatex *texb_6c = new TLatex(0.05,0.2, aMomentsTitle[6]);
  texb_6c->SetTextSize(0.04);
  texb_6c->SetTextAngle(90);
  texb_6c->Draw("same");
  
  return pad;
}

// ______________________________________________________________________________________
void LabelCanvas(const Char_t* analysisTitle, const Char_t *analysisDetails) {
  // -- add Labels to canvas
 
  pad->Modified();
  pad->cd();

  TPad *pad2 = new TPad("pad2", "pad2",0.05, 0.06, 0.99, 0.98);
  pad2->SetBorderMode(0);
  pad2->SetFillColor(1182);
  pad2->Draw();
  pad2->cd();

  TLatex *texb_3b = new TLatex(0.08, 0.94, analysisTitle);
  texb_3b->SetTextSize(0.04);
  texb_3b->Draw("same");

  TLatex *texb_3 = new TLatex(0.08, 0.90, "Au+Au collisions at RHIC");
  texb_3->SetTextSize(0.035);
  texb_3->Draw("same");
  
  TLatex *texb_3a = new TLatex(0.08, 0.86, analysisDetails);
  texb_3a->SetTextSize(0.035);
  texb_3a->Draw("same");
  
  TLatex *texb_4 = new TLatex(0.08, 0.60, "STAR Preliminary");
  texb_4->SetTextSize(0.04);
  texb_4->Draw("same");


  legExp->Draw("lt");
  legTheo->Draw("lt");

  pad->Modified();
  can->Modified();
}

// ______________________________________________________________________________________
void CreateLegends(Int_t linesExp, Int_t linesTheo, Float_t leftX, Float_t topY) {
  // -- Create legends

  float widthX = 0.19;
  float widthY = 0.03;
  float spacer = 0.02;
  float textSize = 0.03;

  legExp = new TLegend(leftX, topY-(linesExp*widthY), leftX+widthX, topY);
  legExp->SetTextAlign(12);
  legExp->SetTextSize(textSize);
  legExp->SetTextFont(42);
  legExp->SetFillColor(0);
  legExp->SetLineColor(0);
  legExp->SetBorderSize(0);

  legTheo = new TLegend(leftX+widthX+spacer, topY-(linesTheo*widthY), leftX+(2*widthX+0.04)+spacer, topY);
  legTheo->SetTextAlign(12);
  legTheo->SetTextSize(textSize);
  legTheo->SetTextFont(42);
  legTheo->SetFillColor(0);
  legTheo->SetLineColor(0);
  legTheo->SetBorderSize(0);
}

// ______________________________________________________________________________________
void SaveCanvas(const Char_t* name, Bool_t isNice = kTRUE) {
  // -- Write out canvas
  
  gSystem->Exec(Form("mkdir -p results/nice/%s/png",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/pdf",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/eps",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/gif",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/root", name));
  if (isNice) {
    gSystem->Exec(Form("mkdir -p results/nice/pdf"));
    gSystem->Exec(Form("mkdir -p results/nice/png"));
  }

  // -----------------------------------------------------
  
  for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
    TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
    if (!c)
      continue;
    
    c->SaveAs(Form("results/nice/%s/png/%s.png",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/eps/%s.eps",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/gif/%s.gif",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/pdf/%s.pdf",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/root/%s.C",    name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/root/%s.root", name, c->GetName()));
    if (isNice) {
      c->SaveAs(Form("results/nice/pdf/%s.pdf",            c->GetName()));
      c->SaveAs(Form("results/nice/png/%s.png",            c->GetName()));
    }
  }
}

// ______________________________________________________________________________________
void SetupStyle() {
  // -- Setup Style

  gROOT->SetStyle("Plain");
  
  gStyle->SetHatchesSpacing(0.8);
  gStyle->SetHatchesLineWidth(1);
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFillStyle(1001);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);

  gStyle->SetLegendBorderSize(0);

  Int_t font = 42;
  
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  
  gStyle->SetTickLength(0.02,"xy");
  gStyle->SetEndErrorSize(3);
  
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.009,"xyz");

  TGaxis::SetMaxDigits(3);

  gStyle->SetTitleSize(0.06);  
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.12,"x"); 
  gStyle->SetTitleOffset(1.13,"yz");   
  gStyle->SetTitleSize(0.045,"xyz");  

  gStyle->SetMarkerSize(1.6);
  gStyle->SetPalette(1,0); 

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetLineWidth(2);

  gStyle->SetEndErrorSize(4);

  // -- Set plot styles for 2D colz
  const Int_t nRGBs = 5;
  const Int_t nCont = 255;
  
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);
}

// ______________________________________________________________________________________
void toolsEnergyNice(){;}
  

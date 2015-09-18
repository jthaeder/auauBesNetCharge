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
#include "TColor.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphErrors.h"

TObjArray canA;

// -----------------------------------------------------------

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma", "#kappa #sigma^{2}"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV"};
const Int_t   nMoments         = 7;

const Char_t* aDataSetsTitle[] = {"corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}"};
const Char_t* aDataSets[]      = {"twoeff_11"};
const int     nDataSets        = 1;

      Int_t   nCent            = 1;
const Char_t* cent[9]          = {"0005",  "0510",   "1020",   "2030",   "3040",   "4050",   "5060",   "6070",   "7080"};
const Char_t* cent1[9]         = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%"};

// -----------------------------------------------------------

      Color_t aColors[]        = {kBlack, kOrange+9, kOrange+9, kYellow, kAzure, kAzure, kCyan+2, kBlue+1, kRed+2, kRed+2};

      Int_t   aMarkers[]       = {24, 32, 25, 25, 26, 28, 27, 30, 32, 30};

      Float_t aMinY[7]         = { 0, 0, -10, -2300, 5.5, -0.02, -12 };
      Float_t aMaxY[7]         = { 40, 250, 90, 900, 9.8, 0.48, 6};

      Float_t aMinX            = -0.55;
      Float_t aMaxX            = 0.55;

// -----------------------------------------------------------


// ______________________________________________________________________________________
void PrepareGraph(TGraphErrors* g) {
  // -- Prepare Graph
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");
  g->GetXaxis()->SetNdivisions(6, 5, 0);
  g->GetYaxis()->SetNdivisions(6, 5, 0);
  g->GetXaxis()->SetLabelSize(0.07);
  g->GetYaxis()->SetLabelSize(0.07);
}

// ______________________________________________________________________________________
void ConfigGraph(TGraphErrors* g, Int_t yLimitIdx, Int_t idxMarker) {
  // -- Config Graph
  g->SetMinimum(aMinY[yLimitIdx]);
  g->SetMaximum(aMaxY[yLimitIdx]);
  g->GetXaxis()->SetLimits(aMinX, aMaxX);  
  g->SetMarkerSize(1.4);
  g->SetMarkerStyle(aMarkers[idxMarker]);
  g->SetMarkerColor(aColors[idxMarker]);
  g->SetLineColor(aColors[idxMarker]);
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
TPad* SetupCanvas(const Char_t* canName, const Char_t *canTitle, const Char_t *xTitle, Float_t xPosTitle) {
  // -- setup canvas and pad
  
  canA.Add(new TCanvas(canName, canTitle, 1200, 0 , 600, 1000));
  TCanvas *can = static_cast<TCanvas*>(canA.Last());
  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetBorderSize(0.0);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->cd();
  
  TPad* pad = new TPad("pad", "pad",0.05, 0.06, 0.99, 0.98);
  pad->SetBorderMode(0);
  pad->SetFillColor(0);
  pad->Draw();
  pad->cd();
  pad->Divide(1, 3, 0., 0., 0);

  can->cd();

  TLatex *texb_5 = new TLatex(xPosTitle, 0.03, xTitle);
  texb_5->SetTextSize(0.04);
  texb_5->Draw("same");
  
  TLatex *texb_6a = new TLatex(0.05,0.8, aMomentsTitle[4]);
  texb_6a->SetTextSize(0.04);
  texb_6a->SetTextAngle(90);
  texb_6a->Draw("same");
  
  TLatex *texb_6b = new TLatex(0.05,0.5, aMomentsTitle[5]);
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
void SaveCanvas(const Char_t* name) {
  // -- Write out canvas
  
  gSystem->Exec(Form("mkdir -p results/nice/%s/png",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/pdf",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/eps",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/gif",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/root", name));
  
  // -----------------------------------------------------
  
  for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
    TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
    if (!c)
      continue;
    
    c->SaveAs(Form("results/nice/%s/png/%s_14GeV.png",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/eps/%s_14GeV.eps",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/gif/%s_14GeV.gif",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/pdf/%s_14GeV.pdf",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/root/%s_14GeV.C",    name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/root/%s_14GeV.root", name, c->GetName()));
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

  gStyle->SetLineWidth(1);

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
void toolsEtaNice(){;}

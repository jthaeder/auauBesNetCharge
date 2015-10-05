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
#include "TGAxis.h"
#include "TLine.h"
#include "TArrow.h"
#include "TGraphErrors.h"
#include "data/BESIIerror.h"

TObjArray canA;

// -----------------------------------------------------------

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma", "#kappa #sigma^{2}"};

// ______________________________________________________________________________________
void PrepareGraph(TGraphErrors* g) {
  // -- Prepare Graph
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");
  g->GetXaxis()->SetNdivisions(6, 5, 0);
  g->GetYaxis()->SetNdivisions(6, 5, 0);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetLabelOffset(0.008);
  g->GetXaxis()->SetNoExponent(kTRUE);
}


// ______________________________________________________________________________________
void SaveCanvas(const Char_t* name) {
  // -- Write out canvas
  
  gSystem->Exec(Form("mkdir -p results/nice/%s/png",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/pdf",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/eps",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/gif",  name));
  gSystem->Exec(Form("mkdir -p results/nice/%s/root", name));
  gSystem->Exec(Form("mkdir -p results/nice/pdf"));
  gSystem->Exec(Form("mkdir -p results/nice/png"));
  
  // -----------------------------------------------------
  
  for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
    TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
    if (!c)
      continue;
    
    c->SaveAs(Form("results/nice/%s/png/%s.png",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/eps/%s.eps",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/gif/%s.gif",   name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/pdf/%s.pdf",   name, c->GetName()));
    c->SaveAs(Form("results/nice/pdf/%s.pdf",            c->GetName()));
    c->SaveAs(Form("results/nice/png/%s.png",            c->GetName()));
    c->SaveAs(Form("results/nice/%s/root/%s.C",    name, c->GetName()));
    c->SaveAs(Form("results/nice/%s/root/%s.root", name, c->GetName()));
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
void makeBESIIerror(){
    // -- make nice plot for BES II error prediction

  SetupStyle();
   
  TGraphErrors* data = new TGraphErrors(5, deltaYData, KVData, 0, eKVData); 
  TGraphErrors* esti = new TGraphErrors(8, deltaYEsti, KVEsti, 0, eKVEsti); 

  canA.Add(new TCanvas("canEstimatedBESIIerror", "estimated BES II error", 0, 0 , 600, 500));
  TCanvas* can = static_cast<TCanvas*>(canA.Last());
  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetBorderSize(0.0);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->cd();
  
  TPad* pad = new TPad("pad", "pad",0.05, 0.09, 0.99, 0.98);
  pad->SetBorderMode(0);
  pad->SetFillColor(1182);
  pad->Draw();
  pad->cd();

  can->cd();

  TLatex *texb_5 = new TLatex(0.50, 0.053, "#Delta#it{y}");
  texb_5->SetTextSize(0.07);
  texb_5->Draw("same");
    
  TLatex *texb_6c = new TLatex(0.07, 0.48, aMomentsTitle[6]);
  texb_6c->SetTextSize(0.07);
  texb_6c->SetTextAngle(90);
  texb_6c->Draw("same");
  
  pad->cd();

  PrepareGraph(data);
  PrepareGraph(esti);

  data->SetMaximum(13.);
  data->GetXaxis()->SetLimits(0, 1.8);

  data->SetMarkerStyle(29);
  data->SetMarkerColor(kRed+2);
  data->SetMarkerSize(1.8);
  data->SetLineColor(kRed+2);

  esti->SetLineColor(kAzure);
  esti->SetFillColorAlpha(kAzure, 0.35);
  esti->SetFillColor(kAzure);
  esti->SetFillStyle(3344);
  
  data->Draw("AP");

  TLine *line1 = new TLine(0, 1, 1.8, 1);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw();

  esti->Draw("E3,SAME");
  data->Draw("P,SAME");  

  // -- Add text
  // -----------------------------------------

  TLatex *texb_3b = new TLatex(0.08, 10.5, "Net-Proton, 0-5%");
  texb_3b->SetTextSize(0.05);
  texb_3b->Draw("same");

  TLatex *texb_3 = new TLatex(0.08, 9.5, "Au+Au collisions #sqrt{#it{s}_{NN}} = 7.7 GeV");
  texb_3->SetTextSize(0.045);
  texb_3->Draw("same");
  
  TLatex *texb_4 = new TLatex(0.08, 8.5, "STAR Preliminary");
  texb_4->SetTextSize(0.045);
  texb_4->Draw("same");


  TLegend* leg = new TLegend(0.12, 0.29, 0.42, 0.49);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);

  leg->AddEntry(data,  "0.4 < #it{p}_{T} (GeV/#it{c}) < 2.0", "pl");
  leg->AddEntry(esti,  "Estimated Error BES II", "f");  
  leg->AddEntry(line1, "Poisson Expectation", "l");

  leg->Draw("lt");

  TArrow *tpcA = new TArrow(0.01, 12.1, 1., 12.1 ,0.02,"<|>");
  tpcA->SetAngle(40);
  tpcA->SetLineWidth(2);
  tpcA->Draw();
  
  TArrow *itpcA = new TArrow(0.01, 11.7, 1.6, 11.7 ,0.02,"<|>");
  itpcA->SetAngle(40);
  itpcA->SetLineWidth(2);
  itpcA->Draw();

  TLatex *tpcT = new TLatex(0.45, 12.2, "TPC");
  tpcT->SetTextSize(0.035);
  tpcT->Draw("same");

  TLatex *itpcT = new TLatex(1.25, 11.8, "iTPC");
  itpcT->SetTextSize(0.035);
  itpcT->Draw("same");

  SaveCanvas("estimatedBESIIerror");
  
}

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


void SetGraphStyle(TGraphErrors* g, int color, int style) {
  g->SetMarkerColor(color);
  g->SetMarkerStyle(style);
}

void SetGraph(TGraphErrors* g, const char* title, const char* titleX, const char* titleY, int color, int style, float min, float max) {
  g->SetTitle(title);
  g->GetXaxis()->SetTitle(titleX);
  g->GetYaxis()->SetTitle(titleY);
  g->SetMaximum(max);
  g->SetMinimum(min);

  SetGraphStyle(g, color, style);
}

void Draw(TCanvas** can, TGraphErrors** graph, TGraphErrors** graph14, int idxPad) {
 
  graph[2]->GetXaxis()->SetNdivisions(9,5,0);

  for (Int_t jj = 0; jj < 2; ++jj) {
    can[jj]->cd(idxPad);
    gPad->SetLogx();
    if (idxPad == 1)
    gPad->SetLogy();
    
    graph[0]->Draw("AP");
    graph[1]->Draw("PSAME");
    if (jj == 1) {
      graph14[0]->Draw("PSAME");
      graph14[1]->Draw("PSAME");
    }
    
    can[jj+2]->cd(idxPad);
    if (idxPad == 1)
      gPad->SetLogy();
    
    graph[2]->Draw("AP");
    graph[3]->Draw("PSAME");
    if (jj == 1) {
      graph14[2]->Draw("PSAME");
      graph14[3]->Draw("PSAME");
    }
  }
}



void plotPublished(const Char_t* name = "2015-05-20_0.5", const Char_t* dataSet = "twoeff_11") {

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
  
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.009,"xyz");

  gStyle->SetTitleSize(0.06);  
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.11,"x");   // JMT 1.15	
  gStyle->SetTitleOffset(1.13,"yz");   // JMT 1.15	
  gStyle->SetTitleSize(0.044,"xyz");  


  gStyle->SetMarkerSize(1.4);  // JMT 1.1
  gStyle->SetPalette(1,0); 

  gStyle->SetOptDate(20);

  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("nemrks");
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

  // -----------------------------------------------------------------------
  // -- Fig.4
  // --  (Beam energy dependance of Variance/Mean, Skwness*Std. Deviation and Kurtosis*Variance)
  // -----------------------------------------------------------------------
  
  /*
    Variance/Mean (Sigma^2/M)                                               
    
    0-5% Centrality       
    70-80% Centrality 
  */ 
  
  // ENG (GeV)  sigma^2/M  Stat.Error  Sys.Error  NBD  Skellam
  Double_t publishedVM[2][7][6] = { 
    { {7.7,  3.72455, 0.0319487, 0.0159744, 6.72266, 5.4224},
      {11.5, 5.55127, 0.0366492, 0.0183246, 9.84451, 8.25644},
      {19.6, 8.95204, 0.0255873, 0.0127937, 15.7966, 13.3401},
      {27,   12.7145, 0.0299337, 0.0149669, 21.1496, 18.7789},
      {39,   17.4625, 0.0263059, 0.013153,  30.1741, 25.4531},
      {62.4, 27.3241, 0.0651797, 0.0325899, 49.3896, 39.0525},
      {200,  80.2318, 0.210213,  0.105106,  124.898, 108.848} },
    { {7.7,  7.47332, 0.102868, 0.051434, 14.6179, 9.20048},
      {11.5, 11.8146, 0.174244, 0.0871218, 22.8091, 15.4699},
      {19.6, 17.9191, 0.124656, 0.062328, 35.0745, 21.3709},
      {27,   26.9981, 0.217466, 0.108733, 54.9325, 37.9424},
      {39,   32.6472, 0.204098, 0.102049, 66.1672, 44.6294},
      {62.4, 47.2453, 0.62329,  0.311645, 94.9524, 64.815},
      {200,  99.0916, 1.56255,  0.781274, 201.971, 139.686} }};
  
  
  /*
    Skewness*Std. Deviation (S*sigma)                             
    
    0-5% Centrality       
    70-80% Centrality   
  */
  
  // ENG (GeV)  S*sigma  Stat.Error  Sys.Error  NBD  Skellam
  Double_t publishedSD[2][7][6] = { 
    { {7.7, 0.443332, 0.237605, 0.215502, 0.476933, 0.184454},
      {11.5, 0.503508, 0.208217, 0.174109, 0.290953, 0.121136},
      {19.6, 0.283159, 0.09437, 0.0777163, 0.168042, 0.0749629},
      {27, 0.107052, 0.0750091, 0.0562664, 0.111871, 0.053256},
      {39, 0.0133622, 0.0501263, 0.033736, 0.0768153, 0.0392905},
      {62.4, 0.145121, 0.067052, 0.0899881, 0.069852, 0.0255998},
      {200, -0.0188162, 0.0477387, 0.0409418, 0.0138754, 0.00919266} },
    { {7.7, 0.472407, 0.0347691, 0.0137378, 0.425598,  0.109716},
    {11.5, 0.321615, 0.0291255, 0.0610068, 0.299105,  0.0648058},
      {19.6, 0.204636, 0.0113238, 0.00967549, 0.189537,  0.0468747},
      {27, 0.124164, 0.00965298, 0.00534643, 0.136668,  0.0265044},
      {39, 0.103698, 0.00674335, 0.0388556, 0.103799,  0.0224682},
      {62.4, 0.0711058, 0.0118494, 0.0197959, 0.0674433,  0.0154511},
      {200, 0.0392179, 0.00761239, 0.021385, 0.0307338,  0.00717921} }};
  
  /*
    Kurtosis*Variance (KV)                             
    
    0-5%   Centrality  
    70-80% Centrality       
  */ 
  
  // ENG (GeV)  kppa*sigma^2  Stat.Error  Sys.Error  NBD  Skellam
  Double_t publishedKV[2][7][6] = { 
    { {7.7, -10.6236, 7.47779, 5.2201, 3.02669, 1},
      {11.5, 1.90761, 7.09204, 1.30255, 2.51506, 1},
      {19.6, 1.30334, 3.54197, 0.785453, 2.37491, 1},
      {27,  2.96161, 3.09955, 1.51754, 1.90713, 1},
      {39,  1.59113, 2.22224, 0.62209, 2.35608, 1},
      {62.4, 4.57877, 3.34581, 2.46885, 3.04649, 1},
      {200,  2.57939, 2.42202, 1.84051, 2.02919, 1} },
    { {7.7, 3.00154, 0.263386, 0.147408, 5.63941, 1},
      {11.5, 2.7137, 0.193646, 0.337091, 4.79646, 1},
      {19.6, 2.43159, 0.0922055, 0.0991678, 4.6705, 1},
      {27, 2.70718, 0.0839845, 0.124916, 5.00913, 1},
      {39, 2.64718, 0.0537805, 0.142572, 5.04093, 1},
      {62.4, 2.71035, 0.0964119, 0.316314, 5.06482, 1},
      {200, 2.33482, 0.0692845, 0.244065, 4.90635, 1} }};
  
  /*
    float snn[] = {  5, 7.7, 11.5, 14.5, 19.6,  27,  39, 62.4, 200,  2760};
    float mub[] = {553, 422,  316,  266,  206, 156, 112,   73,  24, 1.772};
    float  tc[] = {119, 140,  152,  156,  160, 163, 164,  165, 166,   166};
  */
  Double_t snn[] = { 7.7, 11.5, 19.6,  27,  39, 62.4, 200};
  Double_t mub[] = { 422,  316,  206, 156, 112,   73,  24};
  Double_t  tc[] = { 140,  152,  160, 163, 164,  165, 166};

  Double_t snn14[] = { 14.5};
  Double_t mub14[] = {  266};
  Double_t  tc14[] = {  156};
  
  
  
  // J. Cleymans et al., PRC 73, 034905 (2006).
  //    mub(sqrt(s)) = 1.308 /￼(1+0.27× sqrt(s))
  //    Tc(sqrt(s))  = 0.166−0.139 × mub^2 −0.053 mub^4
  // !  Unit : √s, μB, T (GeV)
  
  // -----------------------------------------------------------------------  

  TCanvas *can[4];
  for (int ii = 0; ii < 4; ++ii) {
    can[ii] = new TCanvas(Form("can_%d", ii), "", 0, 0 , 600, 1000);
    can[ii]->Divide(1,3, 0, 0.);
  }
  
  Int_t n = 7;
  Double_t y[2][7], ey[2][7];

  Int_t n14 = 1;
  Double_t y14[2][1], ey14[2][1];
  
  // -----------------------------------------------------------------------

  // -- Get new results
  TFile *fKVM = TFile::Open(Form("output/%s/%s/Moments_VM.root", name, dataSet));
  TGraphErrors * kVM = static_cast<TGraphErrors*>(fKVM->Get("VM"));

  TFile *fKSD = TFile::Open(Form("output/%s/%s/Moments_SD.root", name, dataSet));
  TGraphErrors * kSD = static_cast<TGraphErrors*>(fKSD->Get("SD"));

  TFile *fKKV = TFile::Open(Form("output/%s/%s/Moments_KV.root", name, dataSet));
  TGraphErrors * kKV = static_cast<TGraphErrors*>(fKKV->Get("KV"));

  int idxLow = 8;
  int idxHigh = 0;
  int aIdx14[2] = {0, 8};

  // -------------------------------------------------
  //  VM
  // -------------------------------------------------
  TGraphErrors *graphVM[4];
  TGraphErrors *graphVM14[4];
  
  for (Int_t jj = 0; jj < 2; ++jj) {
    for (Int_t ii = 0; ii < n; ++ii) {
      y[jj][ii]  = publishedVM[jj][ii][1];
      ey[jj][ii] = publishedVM[jj][ii][2];
    }

    graphVM[jj]   = new TGraphErrors(n, snn, y[jj], 0, ey[jj]);
    graphVM[jj+2] = new TGraphErrors(n, mub, y[jj], 0, ey[jj]);

    y14[jj][0]  = kVM->GetY()[aIdx14[jj]];
    ey14[jj][0] = kVM->GetEY()[aIdx14[jj]];
    
    graphVM14[jj]   = new TGraphErrors(n14, snn14, y14[jj], 0, ey14[jj]);
    graphVM14[jj+2] = new TGraphErrors(n14, mub14, y14[jj], 0, ey14[jj]);
  }

  SetGraph(graphVM[0], "#sigma^{2}/M - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#sqrt{s_{NN}} (GeV)", "#sigma^{2}/M", kRed+1, 24, 2, 250);
  SetGraphStyle(graphVM[1], kAzure, 25);
  SetGraph(graphVM[2], "#sigma^{2}/M - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#mu_{B} (MeV)",      "#sigma^{2}/M", kRed+1, 24, 2, 250);
  SetGraphStyle(graphVM[3], kAzure, 25);

  SetGraph(graphVM14[0], "#sigma^{2}/M - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#sqrt{s_{NN}} (GeV)", "#sigma^{2}/M", kRed+1, 20, 2, 250);
  SetGraphStyle(graphVM14[1], kAzure, 21);
  SetGraph(graphVM14[2], "#sigma^{2}/M - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#mu_{B} (MeV)",      "#sigma^{2}/M", kRed+1, 20, 2, 250);
  SetGraphStyle(graphVM14[3], kAzure, 21);

  Draw(can, graphVM, graphVM14, 1);

  // -------------------------------------------------
  //  SD
  // -------------------------------------------------
  TGraphErrors *graphSD[4];
  TGraphErrors *graphSD14[4];
  
  for (Int_t jj = 0; jj < 2; ++jj) {
    for (Int_t ii = 0; ii < n; ++ii) {
      y[jj][ii]  = publishedSD[jj][ii][1];
      ey[jj][ii] = publishedSD[jj][ii][2];
    }
    graphSD[jj]   = new TGraphErrors(n, snn, y[jj], 0, ey[jj]);
    graphSD[jj+2] = new TGraphErrors(n, mub, y[jj], 0, ey[jj]);

    y14[jj][0]  = kSD->GetY()[aIdx14[jj]];
    ey14[jj][0] = kSD->GetEY()[aIdx14[jj]];
    
    graphSD14[jj]   = new TGraphErrors(n14, snn14, y14[jj], 0, ey14[jj]);
    graphSD14[jj+2] = new TGraphErrors(n14, mub14, y14[jj], 0, ey14[jj]);
  }

  SetGraph(graphSD[0], "S #sigma - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#sqrt{s_{NN}} (GeV)", "S #sigma", kRed+2, 24, -0.1, 1);
  SetGraphStyle(graphSD[1], kAzure, 25);
  SetGraph(graphSD[2], "S #sigma - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#mu_{B} (MeV)",       "S #sigma", kRed+2, 24, -0.1, 1);
  SetGraphStyle(graphSD[3], kAzure, 25);

  SetGraph(graphSD14[0], "S #sigma - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#sqrt{s_{NN}} (GeV)", "S #sigma", kRed+2, 20, -0.1, 1);
  SetGraphStyle(graphSD14[1], kAzure, 21);
  SetGraph(graphSD14[2], "S #sigma - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#mu_{B} (MeV)",       "S #sigma", kRed+2, 20, -0.1, 1);
  SetGraphStyle(graphSD14[3], kAzure, 21);

  Draw(can, graphSD, graphSD14, 2);

  // -------------------------------------------------
  //  KV
  // -------------------------------------------------
  TGraphErrors *graphKV[4];
  TGraphErrors *graphKV14[4];
  
  for (Int_t jj = 0; jj < 2; ++jj) {
    for (Int_t ii = 0; ii < n; ++ii) {
      y[jj][ii]  = publishedKV[jj][ii][1];
      ey[jj][ii] = publishedKV[jj][ii][2];
    }
    graphKV[jj]   = new TGraphErrors(n, snn, y[jj], 0, ey[jj]);
    graphKV[jj+2] = new TGraphErrors(n, mub, y[jj], 0, ey[jj]);
    
    y14[jj][0]  = kKV->GetY()[aIdx14[jj]];
    ey14[jj][0] = kKV->GetEY()[aIdx14[jj]];
    
    graphKV14[jj]   = new TGraphErrors(n14, snn14, y14[jj], 0, ey14[jj]);
    graphKV14[jj+2] = new TGraphErrors(n14, mub14, y14[jj], 0, ey14[jj]);
  }

  SetGraph(graphKV[0], "#kappa #sigma^{2} - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#sqrt{s_{NN}} (GeV)", "#kappa #sigma^{2}", kRed+1, 24, -15, 10);
  SetGraphStyle(graphKV[1], kAzure, 25);
  SetGraph(graphKV[2], "#kappa #sigma^{2} - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#mu_{B} (MeV)",      "#kappa #sigma^{2}", kRed+1, 24, -15, 10);
  SetGraphStyle(graphKV[3], kAzure, 25);

  SetGraph(graphKV14[0], "#kappa #sigma^{2} - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#sqrt{s_{NN}} (GeV)", "#kappa #sigma^{2}", kRed+1, 20, -15, 10);
  SetGraphStyle(graphKV14[1], kAzure, 21);
  SetGraph(graphKV14[2], "#kappa #sigma^{2} - Net-Charge , 0.2 < #it{p}_{T} < 2.0 , |#eta| < 0.5", "#mu_{B} (MeV)",      "#kappa #sigma^{2}", kRed+1, 20, -15, 10);
  SetGraphStyle(graphKV14[3], kAzure, 21);

  Draw(can, graphKV, graphKV14, 3);
  
  TLegend * leg2 = new TLegend(0.13, 0.79, 0.75, 0.88);
  leg2->SetTextAlign(12);
  leg2->SetTextSize(0.035);
  leg2->SetTextFont(42);
  leg2->SetFillColor(1182);
  leg2->SetLineColor(0);
  leg2->SetBorderSize(0);

  leg2->AddEntry(graphKV[0], "0-5% - published", "p");
  leg2->AddEntry(graphKV[1], "70-80% - published", "p");

  can[0]->cd(2);
  leg2->Draw();
  can[2]->cd(2);
  leg2->Draw();


  TLegend * leg4 = new TLegend(0.13, 0.70, 0.75, 0.88);
  leg4->SetTextAlign(12);
  leg4->SetTextSize(0.035);
  leg4->SetTextFont(42);
  leg4->SetFillColor(1182);
  leg4->SetLineColor(0);
  leg4->SetBorderSize(0);

  leg4->AddEntry(graphKV[0],   "0-5% - published", "p");
  leg4->AddEntry(graphKV[1],   "70-80% - published", "p");
  leg4->AddEntry(graphKV14[0], "0-5% - corrected - 11.5 GeV #epsilon_{1}, #epsilon_{2}", "p");
  leg4->AddEntry(graphKV14[1], "70-80% - corrected - 11.5 GeV #epsilon_{1}, #epsilon_{2}", "p");

  can[1]->cd(2);
  leg4->Draw();
  can[3]->cd(2);
  leg4->Draw();

  can[0]->SaveAs(Form("results/%s/png/can_NetCharge_Ratio_snn.png",   name));
  can[0]->SaveAs(Form("results/%s/pdf/can_NetCharge_Ratio_snn.pdf",   name));
  can[0]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_snn.C",    name));
  can[0]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_snn.root", name));
  
  can[1]->SaveAs(Form("results/%s/png/can_NetCharge_Ratio_snn_14GeV_%s.png",   name, dataSet));
  can[1]->SaveAs(Form("results/%s/pdf/can_NetCharge_Ratio_snn_14GeV_%s.pdf",   name, dataSet));
  can[1]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_snn_14GeV_%s.C",    name, dataSet));
  can[1]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_snn_14GeV_%s.root", name, dataSet));

  can[2]->SaveAs(Form("results/%s/png/can_NetCharge_Ratio_mub.png",   name));
  can[2]->SaveAs(Form("results/%s/pdf/can_NetCharge_Ratio_mub.pdf",   name));
  can[2]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_mub.C",    name));
  can[2]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_mub.root", name));
  
  can[3]->SaveAs(Form("results/%s/png/can_NetCharge_Ratio_mub_14GeV_%s.png",   name, dataSet));
  can[3]->SaveAs(Form("results/%s/pdf/can_NetCharge_Ratio_mub_14GeV_%s.pdf",   name, dataSet));
  can[3]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_mub_14GeV_%s.C",    name, dataSet));
  can[3]->SaveAs(Form("results/%s/root/can_NetCharge_Ratio_mub_14GeV_%s.root", name, dataSet));
}




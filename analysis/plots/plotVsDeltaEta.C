
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

const Char_t* aEtaSets[]         = {"2015-06-21_delta_0.1_0",
				    "2015-05-20_0.1",
				    "2015-06-03_delta_0.3_8",
				    "2015-05-20_0.2", 
				    "2015-06-20_delta_0.5_6", 
				    "2015-05-20_0.3", "2015-06-20_delta_0.7_0", "2015-05-20_0.4", "2015-06-20_delta_0.9_0", "2015-05-20_0.5"};  

const int     nEtaSets           = 10; 

const float   aEtaSetBinCenter[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const float   aEtaSetBinWidth[]  = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

//const int     nEtaSuperSets           = 7; 
// const int     aEtaSuperSets[]         = {  9,   6,   1,   1,   1,   1,   1}; 
// const float   aEtaSuperSetsBinWidth[] = {0.3, 0.5, 1.0, 0.8, 0.6, 0.4, 0.2}; 

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma", "#kappa #sigma^{2}"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV"};
const int     nMoments         = 7;

// const Char_t* aDataSetsTitle[] = {"uncorrected",
// 				  "corrected - 11.5 GeV (#epsilon_{1} + #epsilon_{2}) / 2",
// 				  "corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}", 
// 				  "corrected - 19.6 GeV (#epsilon_{1} + #epsilon_{2}) / 2",
// 				  "corrected - 19.6 GeV #epsilon_{1} , #epsilon_{2}"};

// const Char_t* aDataSets[]      = {"effuncorr", "11", "twoeff_11", "19", "twoeff_19"};
// const int     nDataSets        = 5;


const Char_t* aDataSetsTitle[] = {"uncorrected",
				  "corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}"};
const Char_t* aDataSets[]      = {"effuncorr", "twoeff_11"};
const int     nDataSets        = 2;

const char *cent[9]  = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};
const char *cent1[9] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"};


TObjArray canA;

void plotVsDeltaEta(const Char_t* name = "deltaEtaScan") {

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
  gStyle->SetTitleOffset(1.12,"x");   // JMT 1.15	
  gStyle->SetTitleOffset(1.13,"yz");   // JMT 1.15	
  gStyle->SetTitleSize(0.045,"xyz");  


  gStyle->SetMarkerSize(1.6);  // JMT 1.1
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

  // -----------------------------------------------------


  TFile *inFiles[nEtaSets][nDataSets][nMoments];

  TGraphErrors *inGraphs[nEtaSets][nDataSets][nMoments];

  for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) 
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) { 
	inFiles[idxEta][idxDataSet][idxMoment] = TFile::Open(Form("output/%s/%s/Moments_%s.root", 
								  aEtaSets[idxEta], aDataSets[idxDataSet], aMoments[idxMoment]));
	inGraphs[idxEta][idxDataSet][idxMoment] = static_cast<TGraphErrors*>((inFiles[idxEta][idxDataSet][idxMoment]->Get(aMoments[idxMoment]))->Clone());
	if (inFiles[idxEta][idxDataSet][idxMoment])
	  (inFiles[idxEta][idxDataSet][idxMoment])->Close();
      }
  
  // -----------------------------------------------------
  
  TCanvas *can[nDataSets];
  TCanvas *canRat[nDataSets];

  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    canA.Add(new TCanvas(Form("canDeltaEta_Cumulants_%s", aDataSets[idxDataSet]), aDataSetsTitle[idxDataSet], 0, 0 , 1100, 600));
    can[idxDataSet] = static_cast<TCanvas*>(canA.Last());
    can[idxDataSet]->Divide(2,2, 0.002, 0.002);
    
    canA.Add(new TCanvas(Form("canDeltaEta_Ratio_%s", aDataSets[idxDataSet]), aDataSetsTitle[idxDataSet], 0, 0 , 1100, 600));
    canRat[idxDataSet] = static_cast<TCanvas*>(canA.Last());
    canRat[idxDataSet]->Divide(2,2, 0.002, 0.002);
  }

  // -----------------------------------------------------

  TLegend *leg[nDataSets];
  TLegend *legRat[nDataSets];

  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    leg[idxDataSet] = new TLegend(0.12, 0.12, 0.68, 0.50);
    leg[idxDataSet]->SetTextAlign(12);
    leg[idxDataSet]->SetTextSize(0.04);
    leg[idxDataSet]->SetTextFont(42);
    leg[idxDataSet]->SetFillColor(1182);
    leg[idxDataSet]->SetLineColor(0);
    leg[idxDataSet]->SetBorderSize(0);
    leg[idxDataSet]->SetHeader(aDataSetsTitle[idxDataSet]);
    
    legRat[idxDataSet] = new TLegend(0.12, 0.12, 0.70, 0.88);
    legRat[idxDataSet]->SetTextAlign(12);
    legRat[idxDataSet]->SetTextSize(0.04);
    legRat[idxDataSet]->SetTextFont(42);
    legRat[idxDataSet]->SetFillColor(1182);
    legRat[idxDataSet]->SetLineColor(0);
    legRat[idxDataSet]->SetBorderSize(0);
    legRat[idxDataSet]->SetHeader(aDataSetsTitle[idxDataSet]);
  }

  // -----------------------------------------------------

  TGraphErrors *etaGraphs[9][nDataSets][nMoments];

  for (int idxCent = 0; idxCent < 9; idxCent++) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	float y[nEtaSets];
	float ey[nEtaSets];
	for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) {	
	  y[idxEta]  = inGraphs[idxEta][idxDataSet][idxMoment]->GetY()[idxCent];
	  ey[idxEta] = inGraphs[idxEta][idxDataSet][idxMoment]->GetEY()[idxCent];
	} 
	  
	etaGraphs[idxCent][idxDataSet][idxMoment] = new TGraphErrors(nEtaSets, aEtaSetBinWidth, y, 0, ey);  

	TGraphErrors *g = etaGraphs[idxCent][idxDataSet][idxMoment];
	
	g->SetTitle(Form("%s - Net-Charge, 0.2 < #it{p}_{T} < 2.0, [%s]", aMomentsTitle[idxMoment], cent1[idxCent]));
	g->GetXaxis()->SetTitle("#Delta#eta");
	g->GetYaxis()->SetTitle(aMomentsTitle[idxMoment]);
	
      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
    } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  } // for (int idxCent = 0; idxCent < 9; idxCent++) {
  
  // -----------------------------------------------------

  for (int idxCent = 0; idxCent < 9; idxCent++) {

    if (idxCent != 0 && idxCent != 1 && idxCent != 8 && idxCent != 4 )
      continue;

    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	  
	if (idxMoment < 4) 
	  can[idxDataSet]->cd(idxMoment+1);
	else 
	  canRat[idxDataSet]->cd(idxMoment-3);
	
	TGraphErrors *g = etaGraphs[idxCent][idxDataSet][idxMoment];
	
	const Color_t aColorsRB[] = {kRed+1, kOrange+9, kOrange-3, kYellow, kSpring+4, kGreen+2, kCyan+2, kBlue+1, kViolet-1, kPink+4};
	
	const float aMinY[2][7] = { { -5, -5, -5, -200, 0, 0, -2 },
				    { -5, -5, -5, -2300, 2, 0, -12 } };
	
	const float aMaxY[2][7] = { { 50, 200, 40, 150, 30, 0.4, 4},
				    { 50, 250, 100, 900, 22, 0.5, 4} };
	
	g->SetMinimum(aMinY[idxDataSet][idxMoment]);
	g->SetMaximum(aMaxY[idxDataSet][idxMoment]);
	
	if (idxCent == 0) {
	  g->SetMarkerColor(aColorsRB[9]);
	  g->SetMarkerStyle(30);
	  g->SetMarkerSize(1.4);
	  g->GetXaxis()->SetLimits(0, 1.1);
	  g->Draw("AP");
	  TLine *line05 = new TLine( 0.5, aMinY[idxDataSet][idxMoment],  0.5, aMaxY[idxDataSet][idxMoment]);
	  line05->SetLineColor(kGray);
	  line05->SetLineStyle(3);
	  line05->Draw();
	  
	  TLine *line50 = new TLine(-0.5, aMinY[idxDataSet][idxMoment], -0.5, aMaxY[idxDataSet][idxMoment]);
	  line50->SetLineColor(kGray);
	  line50->SetLineStyle(3);
	  line50->Draw();
	  
	  TLine *line00 = new TLine(0., aMinY[idxDataSet][idxMoment],0, aMaxY[idxDataSet][idxMoment]);
	  line00->SetLineColor(kGray);
	  line00->SetLineStyle(3);
	  line00->Draw();
	} else if (idxCent == 1) {
	  g->SetMarkerColor(aColorsRB[8]);
	  g->SetMarkerStyle(32);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 2) {
	  g->SetMarkerColor(aColorsRB[7]);
	  g->SetMarkerStyle(26);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 3) {
	  g->SetMarkerColor(aColorsRB[6]);
	  g->SetMarkerStyle(25);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 4) {
	  g->SetMarkerColor(aColorsRB[5]);
	  g->SetMarkerStyle(28);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 5) {
	  g->SetMarkerColor(aColorsRB[4]);
	  g->SetMarkerStyle(29);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 6) {
	  g->SetMarkerColor(aColorsRB[3]);
	  g->SetMarkerStyle(30);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 7) {
	  g->SetMarkerColor(aColorsRB[2]);
	  g->SetMarkerStyle(32);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 8) {
	  g->SetMarkerColor(aColorsRB[1]);
	  g->SetMarkerStyle(33);
	  g->SetMarkerSize(1.4);
	} else if (idxCent == 9) {
	  g->SetMarkerColor(aColorsRB[0]);
	  g->SetMarkerStyle(20);
	  g->SetMarkerSize(1.4);
	}
	g->Draw("PSAME");
      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
      leg[idxDataSet]->AddEntry(etaGraphs[idxCent][idxDataSet][0],    
				Form("%s, 0.2< #it{p}_{T}< 2.0 (GeV/#it{c})", cent1[idxCent]), "pl");
      legRat[idxDataSet]->AddEntry(etaGraphs[idxCent][idxDataSet][0],    
				   Form("%s, 0.2< #it{p}_{T}< 2.0 (GeV/#it{c})", cent1[idxCent]), "pl");
    } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  } // for (int idxCent = 0; idxCent < 9; idxCent++) { 
  
  // -----------------------------------------------------

  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    //can[idxDataSet]->cd(4);
    //leg[idxDataSet]->Draw();
    
    canRat[idxDataSet]->cd(4);
    legRat[idxDataSet]->Draw();
  }
 
  // -----------------------------------------------------

  gSystem->Exec(Form("mkdir -p results/%s/png",  name));
  gSystem->Exec(Form("mkdir -p results/%s/pdf",  name));
  gSystem->Exec(Form("mkdir -p results/%s/root", name));

  // -----------------------------------------------------

  for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
    TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
    if (!c)
      continue;
    
    c->SaveAs(Form("results/%s/png/%s_14GeV.png",   name, c->GetName()));
    c->SaveAs(Form("results/%s/pdf/%s_14GeV.pdf",   name, c->GetName()));
    c->SaveAs(Form("results/%s/root/%s_14GeV.C",    name, c->GetName()));
    c->SaveAs(Form("results/%s/root/%s_14GeV.root", name, c->GetName()));
  }
}

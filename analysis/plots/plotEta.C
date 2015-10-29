
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

const Char_t* aEtaSets[]         = {"2015-06-03_delta_0.3_0", "2015-06-03_delta_0.3_1", "2015-06-03_delta_0.3_2", "2015-06-03_delta_0.3_3", 
				    "2015-06-03_delta_0.3_8",
				    "2015-06-03_delta_0.3_4", "2015-06-03_delta_0.3_5", "2015-06-03_delta_0.3_6", "2015-06-03_delta_0.3_7", 

				    "2015-06-03_delta_0.5_0", "2015-06-03_delta_0.5_1", "2015-06-03_delta_0.5_2", 
				    "2015-06-20_delta_0.5_6",
				    "2015-06-03_delta_0.5_3", "2015-06-03_delta_0.5_4", "2015-06-03_delta_0.5_5",
				    
				    "2015-06-21_delta_0.1_0",  "2015-05-20_0.1", "2015-06-03_delta_0.3_8", 
				    "2015-05-20_0.2", "2015-06-20_delta_0.5_6", "2015-05-20_0.3", 
				    "2015-06-20_delta_0.7_0", "2015-05-20_0.4", "2015-06-20_delta_0.9_0", "2015-05-20_0.5",

				    "2015-06-20_delta_0.5_6", "2015-06-21_delta_1.0_mult_0.5", "2015-05-20_0.5"};  





const int     nEtaSets           = 29; 

const float   aEtaSetBinCenter[] = {-0.35, -0.25, -0.15, -0.05, 0., 0.05,  0.15,  0.25,  0.35, 
 				    -0.25, -0.15, -0.05, 0, 0.05, 0.15, 0.25,
 				    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
				    0., 0., 0.};

const float   aEtaSetBinWidth[]  = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
 				    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
				    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
				    0.5, 1.0, 1.0};


const int     nEtaSuperSets           = 4; 
const int     aEtaSuperSets[]         = {  9,   7,  10,   3}; 
const float   aEtaSuperSetsBinWidth[] = {0.3, 0.5, 1.0, 1.0}; 
const Char_t* aEtaSuperSetsBinTitle[] = {"#Delta#eta = 0.3", "#Delta#eta = 0.5", "#Delta#eta = [0.1 - 1.0]", "#Delta#eta = 0.5/1.0, fixed mult"}; 

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

void plotEta(const Char_t* name = "etaVsNpart") {

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

  TGraphErrors *graphs[nEtaSets][nDataSets][nMoments];

  for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) 
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) { 
	inFiles[idxEta][idxDataSet][idxMoment] = TFile::Open(Form("output/%s/%s/Moments_%s.root", 
								  aEtaSets[idxEta], aDataSets[idxDataSet], aMoments[idxMoment]));
	graphs[idxEta][idxDataSet][idxMoment] = static_cast<TGraphErrors*>((inFiles[idxEta][idxDataSet][idxMoment]->Get(aMoments[idxMoment]))->Clone());
	TGraphErrors *g = graphs[idxEta][idxDataSet][idxMoment];
	g->SetTitle(Form("%s - Net-Charge, 0.2 < #it{p}_{T} < 2.0", aMomentsTitle[idxMoment]));
	g->GetXaxis()->SetTitle("<N_{part}>");
	g->GetYaxis()->SetTitle(aMomentsTitle[idxMoment]);

	if (inFiles[idxEta][idxDataSet][idxMoment])
	  (inFiles[idxEta][idxDataSet][idxMoment])->Close();
      }
  
  // -----------------------------------------------------
  
  TCanvas *can[nEtaSuperSets][nDataSets];
  TCanvas *canRat[nEtaSuperSets][nDataSets];
  
  for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      canA.Add(new TCanvas(Form("canEta_Cumulants_%s_%d", aDataSets[idxDataSet], idxEtaSuper), aDataSetsTitle[idxDataSet], 0, 0 , 1100, 600));
      can[idxEtaSuper][idxDataSet] = static_cast<TCanvas*>(canA.Last());
      can[idxEtaSuper][idxDataSet]->Divide(2,2, 0.002, 0.002);
      
      canA.Add(new TCanvas(Form("canEta_Ratio_%s_%d", aDataSets[idxDataSet], idxEtaSuper), aDataSetsTitle[idxDataSet], 0, 0 , 1100, 600));
      canRat[idxEtaSuper][idxDataSet] = static_cast<TCanvas*>(canA.Last());
      canRat[idxEtaSuper][idxDataSet]->Divide(2,2, 0.002, 0.002);
    }
  }

  // -----------------------------------------------------

  TLegend *leg[nEtaSuperSets][nDataSets];
  TLegend *legRat[nEtaSuperSets][nDataSets];

  for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      leg[idxEtaSuper][idxDataSet] = new TLegend(0.12, 0.12, 0.68, 0.50);
      leg[idxEtaSuper][idxDataSet]->SetTextAlign(12);
      leg[idxEtaSuper][idxDataSet]->SetTextSize(0.04);
      leg[idxEtaSuper][idxDataSet]->SetTextFont(42);
      leg[idxEtaSuper][idxDataSet]->SetFillColor(1182);
      leg[idxEtaSuper][idxDataSet]->SetLineColor(0);
      leg[idxEtaSuper][idxDataSet]->SetBorderSize(0);
      leg[idxEtaSuper][idxDataSet]->SetHeader(Form("%s %s", aDataSetsTitle[idxDataSet], aEtaSuperSetsBinTitle[idxEtaSuper]));
					      
      legRat[idxEtaSuper][idxDataSet] = new TLegend(0.12, 0.18, 0.70, 0.88);
      legRat[idxEtaSuper][idxDataSet]->SetTextAlign(12);
      legRat[idxEtaSuper][idxDataSet]->SetTextSize(0.04);
      legRat[idxEtaSuper][idxDataSet]->SetTextFont(42);
      legRat[idxEtaSuper][idxDataSet]->SetFillColor(1182);
      legRat[idxEtaSuper][idxDataSet]->SetLineColor(0);
      legRat[idxEtaSuper][idxDataSet]->SetBorderSize(0);
      legRat[idxEtaSuper][idxDataSet]->SetHeader(Form("%s %s", aDataSetsTitle[idxDataSet], aEtaSuperSetsBinTitle[idxEtaSuper]));
    }
  }

  // -----------------------------------------------------
  
  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
      
      int idxEta = -1;
      for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   

	if (idxMoment < 4) 
	  can[idxEtaSuper][idxDataSet]->cd(idxMoment+1);
	else 
	  canRat[idxEtaSuper][idxDataSet]->cd(idxMoment-3);
      
	for (int idxEtaSub = 0 ; idxEtaSub < aEtaSuperSets[idxEtaSuper]; ++idxEtaSub) {
	  ++idxEta;
	  
	  TGraphErrors *g = graphs[idxEta][idxDataSet][idxMoment];
	  
	  const Color_t aColorsRB[] = {kRed+1, kOrange+9, kOrange-3, kYellow, kSpring+4, kGreen+2, kCyan+2, kBlue+1, kViolet-1, kPink+4};
	  
	  if (idxMoment == 0) {
	    g->SetMaximum(40);
	  } else if (idxMoment == 1) {
	    g->SetMaximum(250);
	  } else if (idxMoment == 2) {
	    g->SetMaximum(100);
	  } else if (idxMoment == 3) {
	    g->SetMaximum(900);
	    g->SetMinimum(-2300);
	  } else if (idxMoment == 4) {
	    g->SetMaximum(25);
	    g->SetMinimum(4);
	  } else if (idxMoment == 5) {
	    g->SetMaximum(0.6);
	    g->SetMinimum(0.0);
	  } else if (idxMoment == 6) {
	    g->SetMaximum(5);
	    g->SetMinimum(-8);
	  }

	  g->SetMarkerColor(aColorsRB[idxEtaSub]);
	  g->SetMarkerSize(1.4);
	  
	  if (idxEtaSub == 0) {
	    g->SetMarkerStyle(30);
	  } else if (idxEtaSub == 1) {
	    g->SetMarkerStyle(32);
	  } else if (idxEtaSub == 2) {
	    g->SetMarkerStyle(26);
	  } else if (idxEtaSub == 3) {
	    g->SetMarkerStyle(25);
	  } else if (idxEtaSub == 4) {
	    g->SetMarkerStyle(27);
	  } else if (idxEtaSub == 5) {
	    g->SetMarkerStyle(28);
	  } else if (idxEtaSub == 6) {
	    g->SetMarkerStyle(29);
	  } else if (idxEtaSub == 7) {
	    g->SetMarkerStyle(30);
	  } else if (idxEtaSub == 8) {
	    g->SetMarkerStyle(32);
	  } else if (idxEtaSub == 9) {
	    g->SetMarkerStyle(20);
	  }

	  if (idxEtaSub == 0) 
	    g->Draw("AP");
	  else
	    g->Draw("PSAME");
	  

	  const char* multT = "";
	  if (idxEtaSuper == 3) {
	    if (idxEtaSub == 1){
	      multT = " (mult of #Delta#eta = 0.5)";
	      g->SetMarkerColor(kAzure);
	      g->SetMarkerStyle(24);
	    }
	  }

	  if (idxMoment == 0) 
	    leg[idxEtaSuper][idxDataSet]->AddEntry(graphs[idxEta][idxDataSet][0],    
						   Form("#Delta#eta = %.1f%s, 0.2<#it{p}_{T}<2.0 (GeV/#it{c})", aEtaSetBinWidth[idxEta], multT), "lp");
	  if (idxMoment == 4) 
	    legRat[idxEtaSuper][idxDataSet]->AddEntry(graphs[idxEta][idxDataSet][4], 
						      Form("#Delta#eta = %.1f%s, 0.2<#it{p}_{T}<2.0 (GeV/#it{c})", aEtaSetBinWidth[idxEta], multT), "pl");

	} // for (int idxEtaSub = 0 ; idxEtaSub < aEtaSuperSets[idxEtaSuper]; ++idxEtaSub) { 
	
      } //for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   
    } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
  } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  
  // -----------------------------------------------------
  
  for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      // can[idxEtaSuper][idxDataSet]->cd(4);
      // leg[idxEtaSuper][idxDataSet]->Draw();

      canRat[idxEtaSuper][idxDataSet]->cd(4);
      legRat[idxEtaSuper][idxDataSet]->Draw();
    }
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


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

      float   aEta[]           = {0.1, 0.2, 0.3, 0.4, 0.5};
const int     nEta             = 1;

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma", "#kappa #sigma^{2}"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV"};
const int     nMoments         = 7;

const Char_t* aDataSetsTitle[] = {"uncorrected",
				  "corrected - 11.5 GeV (#epsilon_{1} + #epsilon_{2}) / 2",
				  "corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}", 
				  "corrected - 19.6 GeV (#epsilon_{1} + #epsilon_{2}) / 2",
				  "corrected - 19.6 GeV #epsilon_{1} , #epsilon_{2}"};

const Char_t* aDataSets[]      = {"effuncorr", "11", "twoeff_11", "19", "twoeff_19"};
const int     nDataSets        = 5;

TObjArray canA;

void plotSet(const Char_t* name = "2015-05-20", float etaMax = 0.5) {

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
  
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");

  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.1,"x");   // JMT 1.15	
  gStyle->SetTitleOffset(1.1,"yz");   // JMT 1.15	
  gStyle->SetTitleSize(0.04,"xyz");  
  gStyle->SetTitleSize(0.04);  

  gStyle->SetMarkerSize(1.2);  // JMT 1.1
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

  aEta[0] = etaMax;

  // -----------------------------------------------------

  TFile *inFiles[nEta][nDataSets][nMoments];

  for (int idxEta = 0 ; idxEta < nEta; ++idxEta) 
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) 
	inFiles[idxEta][idxDataSet][idxMoment] = TFile::Open(Form("output/%s_%.1f/%s/Moments_%s.root", name, 
								  aEta[idxEta], aDataSets[idxDataSet], aMoments[idxMoment]));

  // -----------------------------------------------------
  
  TCanvas *can[2];
  TCanvas *canRat[2];

  canA.Add(new TCanvas(Form("can_Cumulants_eff_11"), "Efficiency 11.5 GeV", 0, 0 , 800, 600));
  can[0] = static_cast<TCanvas*>(canA.Last());
  can[0]->Divide(2,2, 0.002, 0.002);

  canA.Add(new TCanvas(Form("can_Cumulants_eff_11_19"), "Efficiency 11.5 and 19.6 GeV", 0, 0 , 800, 600));
  can[1] = static_cast<TCanvas*>(canA.Last());
  can[1]->Divide(2,2, 0.002, 0.002);

  canA.Add(new TCanvas(Form("can_Ratio_eff_11"), "Efficiency 11.5 GeV", 0, 0 , 1200, 400));
  canRat[0] = static_cast<TCanvas*>(canA.Last());
  canRat[0]->Divide(3,1, 0.002, 0.002);

  canA.Add(new TCanvas(Form("can_Ratio_eff_11_19"), "Efficiency 11.5 and 19.6 GeV", 0, 0 , 1200, 400));
  canRat[1] = static_cast<TCanvas*>(canA.Last());
  canRat[1]->Divide(3,1, 0.002, 0.002);

  // -----------------------------------------------------

  TLegend *leg[2];
  TLegend *legRat[2];

  for (int idx = 0 ; idx < 2; ++idx) {
    leg[idx] = new TLegend(0.12, 0.12, 0.68, 0.50);
    leg[idx]->SetTextAlign(12);
    leg[idx]->SetTextSize(0.04);
    leg[idx]->SetTextFont(42);
    leg[idx]->SetFillColor(1182);
    leg[idx]->SetLineColor(0);
    leg[idx]->SetBorderSize(0);
    leg[idx]->SetHeader(Form("Net-Charge, 0.2 #it{p}_{T}<2.0, |#eta|<%.1f", aEta[0]));
    
    legRat[idx] = new TLegend(0.12, 0.58, 0.70, 0.88);
    legRat[idx]->SetTextAlign(12);
    legRat[idx]->SetTextSize(0.04);
    legRat[idx]->SetTextFont(42);
    legRat[idx]->SetFillColor(1182);
    legRat[idx]->SetLineColor(0);
    legRat[idx]->SetBorderSize(0);
    legRat[idx]->SetHeader(Form("Net-Charge, 0.2 #it{p}_{T}<2.0, |#eta|<%.1f", aEta[0]));
  }

  // -----------------------------------------------------

  TGraphErrors *graphs[nEta][nDataSets][nMoments];

  for (int idxEta = 0 ; idxEta < nEta; ++idxEta) 
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	graphs[idxEta][idxDataSet][idxMoment] = static_cast<TGraphErrors*>(inFiles[idxEta][idxDataSet][idxMoment]->Get(aMoments[idxMoment]));

	TGraphErrors *g = graphs[idxEta][idxDataSet][idxMoment];
	g->SetTitle(Form("%s - Net-Charge, 0.2 < #it{p}_{T} < 2.0, |#eta|<%.1f", aMomentsTitle[idxMoment], aEta[idxEta]));
	g->GetXaxis()->SetTitle("<N_{part}>");
	g->GetYaxis()->SetTitle(aMomentsTitle[idxMoment]);
	
      }
  
  // -----------------------------------------------------
  
  for (int idxEta = 0 ; idxEta < nEta; ++idxEta) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {

	TGraphErrors *g = graphs[idxEta][idxDataSet][idxMoment];

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
	  g->SetMinimum(0);
	} else if (idxMoment == 5) {
	  g->SetMaximum(1);
	  g->SetMinimum(-0.1);
	} else if (idxMoment == 6) {
	  g->SetMaximum(10);
	  g->SetMinimum(-10);
	}

	if (idxDataSet == 0) {
	  g->SetMarkerColor(kSpring);
	  g->SetMarkerStyle(28);
	} else if (idxDataSet == 1) {
	  g->SetMarkerColor(kOrange);
	  g->SetMarkerStyle(27);
	} else if (idxDataSet == 2) {
	  g->SetMarkerColor(kBlack);
	  g->SetMarkerStyle(26);
	} else if (idxDataSet == 3) {
	  g->SetMarkerColor(kAzure);
	  g->SetMarkerStyle(25);
	} else if (idxDataSet == 4) {
	  g->SetMarkerColor(kRed+1);
	  g->SetMarkerStyle(24);
	}

	// -- plot 11 GeV only
	if (idxDataSet < 3) {
	  if (idxMoment < 4) 
	    can[0]->cd(idxMoment+1);
	  else 
	    canRat[0]->cd(idxMoment-3);
	  
	  if (idxDataSet == 0)
	    g->Draw("AP");
	  else
	    g->Draw("PSAME");
	}

	// -- plot 11 and 19 GeV
	if (idxMoment < 4) 
	  can[1]->cd(idxMoment+1);
	else 
	  canRat[1]->cd(idxMoment-3);

	if (idxDataSet == 0)
	  g->Draw("AP");
	else 
	  g->Draw("PSAME");

      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {

      for (int idx = 0 ; idx < 2; ++idx) {
	leg[idx]->AddEntry(graphs[idxEta][idxDataSet][0],    aDataSetsTitle[idxDataSet], "p");
	legRat[idx]->AddEntry(graphs[idxEta][idxDataSet][4], aDataSetsTitle[idxDataSet], "p");
      }
    } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  } // for (int idxEta = 0 ; idxEta < nEta; ++idxEta) {
      
  // -----------------------------------------------------

  for (int idx = 0 ; idx < 2; ++idx) {
    can[idx]->cd(4);
    leg[idx]->Draw();

    canRat[idx]->cd(2);
    legRat[idx]->Draw();
  }

  // -----------------------------------------------------

  gSystem->Exec(Form("mkdir -p results/%s_%.1f/png",  name, etaMax));
  gSystem->Exec(Form("mkdir -p results/%s_%.1f/pdf",  name, etaMax));
  gSystem->Exec(Form("mkdir -p results/%s_%.1f/root", name, etaMax));

  // -----------------------------------------------------

  for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
    TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
    if (!c)
      continue;

    c->SaveAs(Form("results/%s_%.1f/png/%s_14GeV.png",   name, etaMax, c->GetName()));
    c->SaveAs(Form("results/%s_%.1f/pdf/%s_14GeV.pdf",   name, etaMax, c->GetName()));
    c->SaveAs(Form("results/%s_%.1f/root/%s_14GeV.C",    name, etaMax, c->GetName()));
    c->SaveAs(Form("results/%s_%.1f/root/%s_14GeV.root", name, etaMax, c->GetName()));
  }

  // -----------------------------------------------------

  for (int idxEta = 0 ; idxEta < nEta; ++idxEta) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	
	cout << "\n" << aMoments[idxMoment] << " -- Eff: " << aDataSets[idxDataSet] << " \n" << endl;
	graphs[idxEta][idxDataSet][idxMoment]->Print();
	cout << "--------------------------------------------------" << endl;
      }
    }
  }
}

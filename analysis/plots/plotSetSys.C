
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


const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S#sigma", "#kappa#sigma^{2}", "(S#sigma - M/#sigma^{2})"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV", "VM"};
const int     nMoments         = 8;

// const Char_t* aDataSetsTitle[] = {"uncorrected",
// 				  "corrected - 11.5 GeV (#epsilon_{1} + #epsilon_{2}) / 2",
// 				  "corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}", 
// 				  "corrected - 19.6 GeV (#epsilon_{1} + #epsilon_{2}) / 2",
// 				  "corrected - 19.6 GeV #epsilon_{1} , #epsilon_{2}"};

// const Char_t* aDataSets[]      = {"effuncorr", "11", "twoeff_11", "19", "twoeff_19"};
// const int     nDataSets        = 5;

// const Char_t* aDataSetsTitle[] = {"uncorrected",
// 				  "corrected - #epsilon_{1} , #epsilon_{2}"};

//const int     nSysNames = 4;
//const Char_t* aSysNames[] = { "base", "sys_1_0", "sys_1_1", "sys_1_2", "sys_1_3" };
//const Char_t* aSysNames[] = { "base", "sys_0_0", "sys_0_1", "sys_0_2", "sys_0_3" };

const int     nSysSetNames = 5;
const int     aSysSets[nSysSetNames] = {1, 4, 4, 2, 2};  
const int     nSysNames    = 4;
const Char_t* aSysNames[nSysSetNames][nSysNames] = { {"base", "", "", ""},
						     {"sys_1_0", "sys_1_1", "sys_1_2", "sys_1_3"},
						     {"sys_0_0", "sys_0_1", "sys_0_2", "sys_0_3"},
						     {"base_plus5", "base_minus5", "", ""},
						     {"base_plus2minus2","base_minus2plus2", "", ""} };
const int     nSysNamesAll = 13;
const Char_t* aSysNamesAll[] = { "base", 
			      "sys_1_0", "sys_1_1", "sys_1_2", "sys_1_3",
			      "sys_0_0", "sys_0_1", "sys_0_2", "sys_0_3",
			      "base_plus5", "base_minus5", 
			      "base_plus2minus2","base_minus2plus2" };

const Char_t* aDataSets[]      = {"twoeff_11"};
const int     nDataSets        = 1;

const Char_t* aDataSetsTitle[] = {"corrected - #epsilon_{1} , #epsilon_{2}"};

//const Char_t* aSysTitle[] = {"default", "nHitPoints [16,24]"};
//const Char_t* aSysTitle[] = {"default", "eff (+5% && +5%)", "eff (-5% && -5%)", "eff (+5% && -5%)"};
const Char_t* aSysTitle[] = {"default", "eff (+5% && +5%)", "eff (-5% && -5%)", "eff (+2% && -2%)", "eff (-2% && +2%)"};

TObjArray canA;

void plotSetSys(const Char_t* name = "jobs_14.5") {

  float etaMax = 0.5;

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
  gStyle->SetTitleSize(0.045,"yz");  
  gStyle->SetTitleSize(0.045,"x");  

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

  // -----------------------------------------------------

  TFile *inFiles[nSysNamesAll][nDataSets][nMoments];
  TGraphErrors *graphs[nSysNamesAll][nDataSets][nMoments];

  for (int idxSys = 0; idxSys < nSysNamesAll; ++idxSys) 
    for (int idxDataSet = 0; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0; idxMoment < nMoments; ++idxMoment) {

	printf("output/%s_%s/%s/Moments_%s.root\n", name, aSysNamesAll[idxSys], aDataSets[idxDataSet], aMoments[idxMoment]);
	inFiles[idxSys][idxDataSet][idxMoment] = TFile::Open(Form("output/%s_%s/%s/Moments_%s.root", name, aSysNamesAll[idxSys],
								  aDataSets[idxDataSet], aMoments[idxMoment]));
	graphs[idxSys][idxDataSet][idxMoment] = static_cast<TGraphErrors*>(inFiles[idxSys][idxDataSet][idxMoment]->Get(aMoments[idxMoment])->Clone());
	
       	TGraphErrors *g = graphs[idxSys][idxDataSet][idxMoment];
	g->SetName(Form("gSys_%d", idxSys));
	g->SetTitle(Form("%s - Net-Charge, 0.2 < #it{p}_{T} < 2.0, |#eta|< 0.5", aMomentsTitle[idxMoment]));
	g->GetXaxis()->SetTitle("<N_{part}>");
	g->GetYaxis()->SetTitle(aMomentsTitle[idxMoment]);

	g->Print();
	inFiles[idxSys][idxDataSet][idxMoment]->Close();
    
	break;
      }
  //  return;
  // -----------------------------------------------------

  TCanvas *can[4];
  TCanvas *canRat[4];

  for (Int_t idx = 0; idx < 4; ++idx) {
    canA.Add(new TCanvas(Form("can_Sys_%d", idx), "", 0, 0 , 690, 600));
    can[idx] = static_cast<TCanvas*>(canA.Last());
    can[idx]->Divide(2,2, 0.002, 0.002);
    
    canA.Add(new TCanvas(Form("can_Ratio_Sys_%d", idx), "", 0, 0 , 1100, 600));
    canRat[idx] = static_cast<TCanvas*>(canA.Last());
    canRat[idx]->Divide(2,2, 0.002, 0.002);
  }

  // -----------------------------------------------------

  // TLegend *leg[2];
  // TLegend *legRat[2];

  // for (int idx = 0 ; idx < 2; ++idx) {
  //   leg[idx] = new TLegend(0.12, 0.12, 0.68, 0.50);
  //   leg[idx]->SetTextAlign(12);
  //   leg[idx]->SetTextSize(0.04);
  //   leg[idx]->SetTextFont(42);
  //   leg[idx]->SetFillColor(1182);
  //   leg[idx]->SetLineColor(0);
  //   leg[idx]->SetBorderSize(0);
  //   leg[idx]->SetHeader(Form("Net-Charge, 0.2 #it{p}_{T}<2.0, |#eta|<%.1f", aEta[0]));
    
  //   legRat[idx] = new TLegend(0.12, 0.60, 0.70, 0.88);
  //   legRat[idx]->SetTextAlign(12);
  //   legRat[idx]->SetTextSize(0.04);
  //   legRat[idx]->SetTextFont(42);
  //   legRat[idx]->SetFillColor(1182);
  //   legRat[idx]->SetLineColor(0);
  //   legRat[idx]->SetBorderSize(0);
  //   legRat[idx]->SetHeader(Form("Net-Charge, 0.2 #it{p}_{T}<2.0, |#eta|<%.1f", aEta[0]));
  // }

  // -----------------------------------------------------


  // -----------------------------------------------------
  Int_t idxSysAll = 0;
  for (int idxSysSet = 1; idxSysSet < nSysSetNames; ++idxSysSet) {
    for (int idxSys = 0; idxSys < aSysSets[idxSysSet]; ++idxSys) {
      ++idxSysAll;
      
      for (int idxDataSet = 0; idxDataSet < nDataSets; ++idxDataSet) {
	for (int idxMoment = 0; idxMoment < nMoments; ++idxMoment) {

	  TGraphErrors *gDefault = graphs[0][idxDataSet][idxMoment];	  
	  TGraphErrors *g        = graphs[idxSysAll][idxDataSet][idxMoment];
	  
	  
	  if (idxMoment == 0) {
	    gDefault->SetMaximum(40);
	  } else if (idxMoment == 1) {
	    gDefault->SetMaximum(250);
	  } else if (idxMoment == 2) {
	    gDefault->SetMaximum(150); // 100
	  } else if (idxMoment == 3) {
	    gDefault->SetMaximum(1800);  // 900
	    gDefault->SetMinimum(-2300);
	  } else if (idxMoment == 4) {
	    gDefault->SetMaximum(20);
	    gDefault->SetMinimum(4);
	  } else if (idxMoment == 5) {
	    gDefault->SetMaximum(0.8);
	    gDefault->SetMinimum(0.1);
	  } else if (idxMoment == 6) {
	    gDefault->SetMaximum(12); //5
	    gDefault->SetMinimum(-8);
	  } else if (idxMoment == 7) {
	    gDefault->SetMaximum(0.6);
	    gDefault->SetMinimum(-0.2);
	  }
	  
	  gDefault->SetMarkerColor(kRed+2);
	  gDefault->SetMarkerStyle(29);
	  gDefault->SetMarkerSize(1.6);
	  
	  if (idxSys == 0) {
	    g->SetMarkerColor(kAzure);
	    g->SetMarkerStyle(24);
	    g->SetMarkerSize(1.4);
	  } else if (idxSys == 1) {
	    g->SetMarkerColor(kBlack);
	    g->SetMarkerStyle(25);
	    g->SetMarkerSize(1.4);
	  } else if (idxSys == 2) {
	    g->SetMarkerColor(kOrange+4);
	    g->SetMarkerStyle(26);
	    g->SetMarkerSize(1.4);
	  } else if (idxSys == 3) {
	    g->SetMarkerColor(kGreen+2);
	    g->SetMarkerStyle(27);
	    g->SetMarkerSize(1.4);
	  } 
	  
	  if (idxMoment < 4) 
	    can[idxSysSet-1]->cd(idxMoment+1);
	  else 
	    canRat[idxSysSet-1]->cd(idxMoment-3);
	  
	  if (idxSys == 0)
	    gDefault->Draw("AP");
	  g->Draw("PSAME");
	  
	  
	} // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	
	// if (idxSys < nSysNames) {
        // 	for (int idx = 0 ; idx < 2; ++idx) {
        // 	  leg[idx]->AddEntry(graphs[idxSys][idxDataSet][0],    aSysTitle[idxSys], "p");
        // 	  legRat[idx]->AddEntry(graphs[idxSys][idxDataSet][4], aSysTitle[idxSys], "p");
        // 	}
        // }
      } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    } // for (int idxSys = 0 ; idxSys < nEta; ++idxSys) {
  }

}

#if 0
  // -----------------------------------------------------

  // for (int idx = 0 ; idx < 2; ++idx) {
  //   can[idx]->cd(4);
  //   leg[idx]->Draw();

  //   canRat[idx]->cd(2);
  //   legRat[idx]->Draw();
  // }

  // -----------------------------------------------------
  return;
  gSystem->Exec(Form("mkdir -p results/%s_%s_%.1f/png",  name, sys, etaMax));
  gSystem->Exec(Form("mkdir -p results/%s_%s_%.1f/pdf",  name, sys, etaMax));
  gSystem->Exec(Form("mkdir -p results/%s_%s_%.1f/root", name, sys, etaMax));

  // -----------------------------------------------------

  for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
    TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
    if (!c)
      continue;

    c->SaveAs(Form("results/%s_%s_%.1f/png/%s_14GeV.png",   name, sys, etaMax, c->GetName()));
    c->SaveAs(Form("results/%s_%s_%.1f/pdf/%s_14GeV.pdf",   name, sys, etaMax, c->GetName()));
    c->SaveAs(Form("results/%s_%s_%.1f/root/%s_14GeV.C",    name, sys, etaMax, c->GetName()));
    c->SaveAs(Form("results/%s_%s_%.1f/root/%s_14GeV.root", name, sys, etaMax, c->GetName()));
  }

  // -----------------------------------------------------

  for (int idxSys = 0 ; idxSys < nEta; ++idxSys) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	
	cout << "\n" << aMoments[idxMoment] << " -- Eff: " << aDataSets[idxDataSet] << " \n" << endl;
	graphs[idxSys][idxDataSet][idxMoment]->Print();
	cout << "--------------------------------------------------" << endl;
      }
    }
  }
}

#endif 

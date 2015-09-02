/* ***************************************************
 *  Fit efficiencies for identified particles   
 *  from embedding - pt dependent
 *
 * *************************************************** *
 *  Input: 
 *   - Output from makeEfficiency from embedding
 *     Seperate for all energies and particles
 *       <baseFolder>/efficiency/<particle><charge><energy>GeV.root
 *     eg.: <baseFolder>/efficiency/piminus11GeV.root
 *
 * *************************************************** *
 *  Output: 
 *    - Fits: 
 *        - results/fits/fit.root
 *    - Canvas: 
 *        - results/particles_pt
 *        - results/energy_cmp_pt
 *
 * *************************************************** *
 *  Macro adopted from Stephen Horvat
 *
 *  Latest changes: Jochen Thaeder <jmthader@lbl.gov>
 * 
 * ***************************************************
 *  Run :  root -l -b -q fitEfficiencyPt.C++
 * ***************************************************/

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

#include "setupStyle.C"

using namespace std;

const int   nCent   = 10;
const char* cent[]  = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};
const char* cent1[] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"};

const int   nNames   = 3; 
const char* names[]  = {"pion", "kaon", "proton"};
const char* names2[] = {"#pi^{+}", "K^{+}", "p"};
const char* names3[] = {"#pi^{-}", "K^{-}", "#bar{p}"};

const float namesPtMin[]  = {0.2, 0.2, 0.4};
const float namesPtMax[]  = {2.0, 1.6, 2.0};

const int namesPtBinMin[]  = {3, 3, 5};
const int namesPtBinMax[]  = {14, 13, 14};

// 1 -> 0 > 0.05 < 0.1
// 2 -> 0.1 > 0.15 < 0.2
// 3 -> 0.2 > 0.25 < 0.3
// 4 -> 0.3 > 0.35 < 0.4
// 5 -> 0.4 > 0.45 < 0.5
// 6 -> 0.5 > 0.55 < 0.6
// 7 -> 0.6 > 0.65 < 0.7
// 8 -> 0.7 > 0.75 < 0.8
// 9 -> 0.8 > 0.85 < 0.9
// 10 -> 0.9 > 0.95 < 1
// 11 -> 1 > 1.05 < 1.1
// 12 -> 1.1 > 1.25 < 1.4
// 13 -> 1.4 > 1.55 < 1.7
// 14 -> 1.7 > 1.85 < 2
// 15 -> 2 > 2.15 < 2.3

const int   nEnergies       = 7;
const float xPosLabel[]     = {0.15, 0.45, 0.72};
const char *energies[]      = {  "7",   "11",   "14",   "19",   "27",   "39",   "62",  "200"};
const char *exactEnergies[] = {"7.7", "11.5", "11.5", "19.6", "27.0", "39.0", "62.4", "62.4"};

double plateau(double *x,double *par) {
  return par[0]*x[0]+par[1];
}

// __________________________________________________________________________________
void makeEfficiencyStudy() {
  gROOT->Reset();
  gROOT->LoadMacro("./setupStyle.C");
  setupStyle();

  gSystem->Exec("mkdir -p ./results/effStudy/png  ./results/effStudy/pdf  ./results/effStudy/root  ./results/effStudy/root_macro");
  gSystem->Exec("mkdir -p ./results/effStudy/profile/png  ./results/effStudy/profile/pdf  ./results/effStudy/profile/root  ./results/effStudy/profile/root_macro");
  gSystem->Exec("mkdir -p ./results/effStudy/support/png  ./results/effStudy/support/pdf  ./results/effStudy/support/root  ./results/effStudy/support/root_macro");

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);

  TFile *fplus[3][nEnergies];
  TFile *fminus[3][nEnergies];

  TCanvas* can[nNames][nEnergies];
  TCanvas* can2D[2*nNames][nEnergies];
  TCanvas* canProfile[nNames][nEnergies]; 

  TH1D* histsPlus[nNames][nEnergies][nCent];
  TH1D* histsMinus[nNames][nEnergies][nCent];

  TH2D* hists2DPlus[nNames][nEnergies][nCent];
  TH2D* hists2DMinus[nNames][nEnergies][nCent];

  TH1D* histsProfilePlus[nNames][nEnergies][nCent];
  TH1D* histsProfileMinus[nNames][nEnergies][nCent];

  // ----------------------------------------------------------
  // -- read files
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    fplus[0][energyIdx]  = TFile::Open(Form("efficiencyStudy/piplus%sGeV.root",  energies[energyIdx]));
    fminus[0][energyIdx] = TFile::Open(Form("efficiencyStudy/piminus%sGeV.root", energies[energyIdx]));
    
    fplus[1][energyIdx]  = TFile::Open(Form("efficiencyStudy/kplus%sGeV.root",   energies[energyIdx]));
    fminus[1][energyIdx] = TFile::Open(Form("efficiencyStudy/kminus%sGeV.root",  energies[energyIdx]));
    
    fplus[2][energyIdx]  = TFile::Open(Form("efficiencyStudy/pplus%sGeV.root",   energies[energyIdx]));
    fminus[2][energyIdx] = TFile::Open(Form("efficiencyStudy/pminus%sGeV.root",  energies[energyIdx]));
  }

  // ----------------------------------------------------------
  // -- read histograms
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	histsPlus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fplus[idx][energyIdx]->Get(Form("hpt_width_%s",cent[centIdx-1])));
	histsPlus[idx][energyIdx][centIdx]->SetName(Form("hpt_width_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));
	histsPlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	
	histsMinus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fminus[idx][energyIdx]->Get(Form("hpt_width_%s",cent[centIdx-1])));
	histsMinus[idx][energyIdx][centIdx]->SetName(Form("hpt_width_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));	
	histsMinus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);

	hists2DPlus[idx][energyIdx][centIdx] = static_cast<TH2D*>(fplus[idx][energyIdx]->Get(Form("hpt_2D_%s",cent[centIdx-1])));
	hists2DPlus[idx][energyIdx][centIdx]->SetName(Form("hpt_2D_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));			

	hists2DMinus[idx][energyIdx][centIdx] = static_cast<TH2D*>(fminus[idx][energyIdx]->Get(Form("hpt_2D_%s",cent[centIdx-1])));
	hists2DMinus[idx][energyIdx][centIdx]->SetName(Form("hpt_2D_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));		
      }
    }
  }

  // ----------------------------------------------------------
  // -- Make Profiles
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	histsProfilePlus[idx][energyIdx][centIdx] = 
	  hists2DPlus[idx][energyIdx][centIdx]->ProjectionY(Form("%s_profile", hists2DPlus[idx][energyIdx][centIdx]->GetName()), namesPtBinMin[idx], namesPtBinMax[idx]);
	histsProfilePlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsProfilePlus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsProfilePlus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsProfilePlus[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsProfilePlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsProfilePlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	histsProfilePlus[idx][energyIdx][centIdx]->SetTitle("");
	histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetTitle("");
	histsProfilePlus[idx][energyIdx][centIdx]->GetYaxis()->SetTitle("");
	histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetLabelSize(0.07);
	histsProfilePlus[idx][energyIdx][centIdx]->GetYaxis()->SetLabelSize(0.05);
	histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetNdivisions(9,5,0);
	histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetRangeUser(-0.1, 1.1);
	histsProfilePlus[idx][energyIdx][centIdx]->GetYaxis()->SetRangeUser(0.4, 2000);



	histsProfileMinus[idx][energyIdx][centIdx] = 
	  hists2DMinus[idx][energyIdx][centIdx]->ProjectionY(Form("%s_profile", hists2DMinus[idx][energyIdx][centIdx]->GetName()), namesPtBinMin[idx], namesPtBinMax[idx]);
	histsProfileMinus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsProfileMinus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsProfileMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsProfileMinus[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsProfileMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsProfileMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
	histsProfileMinus[idx][energyIdx][centIdx]->SetTitle("");
	histsProfileMinus[idx][energyIdx][centIdx]->GetXaxis()->SetTitle("");
	histsProfileMinus[idx][energyIdx][centIdx]->GetYaxis()->SetTitle("");
	histsProfileMinus[idx][energyIdx][centIdx]->GetXaxis()->SetLabelSize(0.07);
	histsProfileMinus[idx][energyIdx][centIdx]->GetYaxis()->SetLabelSize(0.05);
	histsProfileMinus[idx][energyIdx][centIdx]->GetXaxis()->SetNdivisions(9,5,0);
      }
    }
  }

  // ----------------------------------------------------------
  // -- Fill canvas by name / energy -- eff width
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx <nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      can[idx][energyIdx] = new TCanvas(Form("can_%s_%s", names[idx], energies[energyIdx]), names[idx],0,0,1200,600);
      can[idx][energyIdx]->SetFillColor(0);
      can[idx][energyIdx]->SetBorderMode(0);
      can[idx][energyIdx]->SetBorderSize(0.0);
      can[idx][energyIdx]->SetFrameFillColor(0);
      can[idx][energyIdx]->SetFrameBorderMode(0);
      can[idx][energyIdx]->cd();

      TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
      pad->SetBorderMode(0);
      pad->SetFillColor(1182);
      pad->Draw();
      pad->cd();
      pad->Divide(5,2,0.,0.,0.);
      
      TLegend *leg1= new TLegend(0.15,0.25,0.8,0.40);
      leg1->SetName(Form("leg1_%s", names[idx]));
      leg1->SetTextAlign(12);
      leg1->SetTextSize(0.10);
      leg1->SetTextFont(42);
      leg1->SetFillColor(kWhite);
      leg1->SetLineColor(0);
      leg1->SetBorderSize(0);
      
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	pad->cd(centIdx);
	
	if (centIdx == 5) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 9) {
	  gPad->SetRightMargin(0.003);
	}
	
	TH2D *ff = new TH2D("","",20,0.009,2.09,20,0.00,0.19);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.05);
	ff->GetXaxis()->SetNdivisions(9,5,0);
	ff->Draw();

	TLine *line02 = new TLine(0.2,0,0.2,1);
	line02->SetLineColor(kGray+4);
	line02->SetLineStyle(3);
	line02->Draw();

	TLine *line20 = new TLine(2.,0,2.,1);
	line20->SetLineColor(kGray+4);
	line20->SetLineStyle(3);
	line20->Draw();

	histsPlus[idx][energyIdx][centIdx]->Draw("p,same");
	histsMinus[idx][energyIdx][centIdx]->Draw("p,same");
	
	TLatex *texb_Cent = new TLatex(0.8,0.1,cent1[centIdx-1]);
	texb_Cent->SetTextSize(0.09);
	texb_Cent->SetTextFont(42);
	texb_Cent->Draw("same");
	
	pad->Modified();
	can[idx][energyIdx]->cd();
	
	if (centIdx == 1) {
	  leg1->AddEntry(histsPlus[idx][energyIdx][centIdx],  names2[idx],"ple");
	  leg1->AddEntry((TObject*)0, "", "");
	  leg1->AddEntry(histsMinus[idx][energyIdx][centIdx], names3[idx],"ple");
	}
      }
      
      pad->cd(10);

      TH2D *ff = new TH2D("","", 20,0.009,2.09, 20,0.01,0.99);
      ff->GetXaxis()->SetLabelSize(0.07);
      ff->GetYaxis()->SetLabelSize(0.07);
      ff->SetNdivisions(95, "X");
    
      leg1->Draw("lt");
      
      TLatex *texb_1 = new TLatex(0.15,0.82,"Efficiency Width");
      texb_1->SetTextSize(0.1);
      texb_1->SetTextFont(42);
      texb_1->Draw("same");
      
      TLatex *texb_3 = new TLatex(0.15,0.7, Form("AuAu #sqrt{s_{NN}} = %s GeV", exactEnergies[energyIdx]));
      texb_3->SetTextSize(0.08);
      texb_3->SetTextFont(42);
      texb_3->Draw("same");
      
      TLatex *texb_4 = new TLatex(0.15,0.6,"|#eta| < 0.5");
      texb_4->SetTextSize(0.08);
      texb_4->SetTextFont(42);
      texb_4->Draw("same");
      
      TLatex *texb_4a = new TLatex(0.15,0.5,"sample size 20 events each");
      texb_4a->SetTextSize(0.07);
      texb_4a->SetTextFont(42);
      texb_4a->Draw("same");

      pad->Modified();
      can[idx][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.73,0.06,"#it{p}_{T} (GeV/#it{c})");
      texb_5->SetTextSize(0.03);
      texb_5->SetTextFont(42);
      texb_5->Draw("same");
      
      TLatex *texb_6 = new TLatex(0.03,0.4, Form("Efficiency Width- %s", names[idx]));
      texb_6->SetTextSize(0.03);
      texb_6->SetTextFont(42);
      texb_6->SetTextAngle(90);
      texb_6->Draw("same");

      pad->Modified();
      can[idx][energyIdx]->cd();
    } // end idx
  } // end energy idx

  // ----------------------------------------------------------
  // -- Fill canvas by name / energy -- 2D eff
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx <nEnergies; ++energyIdx) {
    for (int idx = 0; idx < 2*nNames; ++idx) {
      int idxName    = idx/2;
      int isNegative = idx%2;
      if (isNegative)
	can2D[idx][energyIdx] = new TCanvas(Form("can2D_%sNegative_%s", names[idxName], energies[energyIdx]), Form("%s negative", names[idxName]),0,0,1200,600);
      else
	can2D[idx][energyIdx] = new TCanvas(Form("can2D_%sPositive_%s", names[idxName], energies[energyIdx]), Form("%s positive", names[idxName]),0,0,1200,600);
      can2D[idx][energyIdx]->SetFillColor(0);
      can2D[idx][energyIdx]->SetBorderMode(0);
      can2D[idx][energyIdx]->SetBorderSize(0.0);
      can2D[idx][energyIdx]->SetFrameFillColor(0);
      can2D[idx][energyIdx]->SetFrameBorderMode(0);
      can2D[idx][energyIdx]->cd();

      TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
      pad->SetBorderMode(0);
      pad->SetFillColor(1182);
      pad->Draw();
      pad->cd();
      pad->Divide(5,2,0.,0.,0.);
           
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	pad->cd(centIdx);
	
	TH2D *hist = (isNegative) ? hists2DMinus[idxName][energyIdx][centIdx] :  hists2DPlus[idxName][energyIdx][centIdx];
	hist->SetTitle("");
	hist->GetXaxis()->SetTitle("");
	hist->GetYaxis()->SetTitle("");
	hist->GetXaxis()->SetRangeUser(0.009, 2.09);
	hist->GetYaxis()->SetRangeUser(0., 1.);
	hist->GetXaxis()->SetLabelSize(0.07);
	hist->GetYaxis()->SetLabelSize(0.05);
	hist->GetXaxis()->SetNdivisions(9,5,0);
	hist->Draw("COLZ");

	if (centIdx == 5) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 9) {
	  gPad->SetRightMargin(0.003);
	}
	
	TLine *line02 = new TLine(0.2,0,0.2,1);
	line02->SetLineColor(kGray+4);
	line02->SetLineStyle(3);
	line02->Draw();

	TLine *line20 = new TLine(2.,0,2.,1);
	line20->SetLineColor(kGray+4);
	line20->SetLineStyle(3);
	line20->Draw();

	TLatex *texb_Cent = new TLatex(0.8,0.1,cent1[centIdx-1]);
	texb_Cent->SetTextSize(0.09);
	texb_Cent->SetTextFont(42);
	texb_Cent->Draw("same");
	
	pad->Modified();
	can2D[idx][energyIdx]->cd();
      }
      
      pad->cd(10);

      TH2D *ff = new TH2D("","", 20,0.009,2.09, 20,0.01,0.99);
      ff->GetXaxis()->SetLabelSize(0.07);
      ff->GetYaxis()->SetLabelSize(0.07);
      ff->SetNdivisions(95, "X");
    
      TLatex *texb_1 = new TLatex(0.15,0.82,"Efficiency vs #it{p}_{T}");
      texb_1->SetTextSize(0.1);
      texb_1->SetTextFont(42);
      texb_1->Draw("same");
      
      TLatex *texb_3 = new TLatex(0.15,0.7, Form("AuAu #sqrt{s_{NN}} = %s GeV", exactEnergies[energyIdx]));
      texb_3->SetTextSize(0.08);
      texb_3->SetTextFont(42);
      texb_3->Draw("same");
      
      TLatex *texb_4 = new TLatex(0.15,0.6,"|#eta| < 0.5");
      texb_4->SetTextSize(0.08);
      texb_4->SetTextFont(42);
      texb_4->Draw("same");
      
      TLatex *texb_4a = new TLatex(0.15,0.5,"sample size 20 events each");
      texb_4a->SetTextSize(0.07);
      texb_4a->SetTextFont(42);
      texb_4a->Draw("same");

      pad->Modified();
      can2D[idx][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.73,0.06,"#it{p}_{T} (GeV/#it{c})");
      texb_5->SetTextSize(0.03);
      texb_5->SetTextFont(42);
      texb_5->Draw("same");
      
      TLatex *texb_6;
      if (isNegative)
	texb_6 = new TLatex(0.03,0.4, Form("Efficiency - %s", names2[idxName]));
      else
	texb_6 = new TLatex(0.03,0.4, Form("Efficiency - %s", names3[idxName]));
      texb_6->SetTextSize(0.03);
      texb_6->SetTextFont(42);
      texb_6->SetTextAngle(90);
      texb_6->Draw("same");

      pad->Modified();
      can2D[idx][energyIdx]->cd();

    } // end idx
  } // end energy idx

  // ----------------------------------------------------------
  // -- Fill canvas by name / energy -- eff width
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx <nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      canProfile[idx][energyIdx] = new TCanvas(Form("canProfile_%s_%s", names[idx], energies[energyIdx]), names[idx],0,0,1200,600);
      canProfile[idx][energyIdx]->SetFillColor(0);
      canProfile[idx][energyIdx]->SetBorderMode(0);
      canProfile[idx][energyIdx]->SetBorderSize(0.0);
      canProfile[idx][energyIdx]->SetFrameFillColor(0);
      canProfile[idx][energyIdx]->SetFrameBorderMode(0);
      canProfile[idx][energyIdx]->cd();

      TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
      pad->SetBorderMode(0);
      pad->SetFillColor(1182);
      pad->Draw();
      pad->cd();
      pad->Divide(5,2,0.,0.,0.);
      
      TLegend *leg1= new TLegend(0.15,0.15,0.8,0.30);
      leg1->SetName(Form("leg1_%s", names[idx]));
      leg1->SetTextAlign(12);
      leg1->SetTextSize(0.10);
      leg1->SetTextFont(42);
      leg1->SetFillColor(kWhite);
      leg1->SetLineColor(0);
      leg1->SetBorderSize(0);
      
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	pad->cd(centIdx);
	gPad->SetLogy();

	if (centIdx == 5) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 9) {
	  gPad->SetRightMargin(0.003);
	}
	
	histsProfilePlus[idx][energyIdx][centIdx]->Draw("e1,same");
	histsProfileMinus[idx][energyIdx][centIdx]->Draw("e1,same");
	
	TLatex *texb_Cent = new TLatex(0.05,500,cent1[centIdx-1]);
	texb_Cent->SetTextSize(0.09);
	texb_Cent->SetTextFont(42);
	texb_Cent->Draw("same");
	
	pad->Modified();
	canProfile[idx][energyIdx]->cd();
	
	if (centIdx == 1) {
	  leg1->AddEntry(histsProfilePlus[idx][energyIdx][centIdx],  names2[idx],"ple");
	  leg1->AddEntry((TObject*)0, "", "");
	  leg1->AddEntry(histsProfileMinus[idx][energyIdx][centIdx], names3[idx],"ple");
	}
      }
      
      pad->cd(10);
      
      TH2D *ff = new TH2D("","", 20,0.009,2.09, 20,0.01,0.99);
      ff->GetXaxis()->SetLabelSize(0.07);
      ff->GetYaxis()->SetLabelSize(0.07);
      ff->SetNdivisions(95, "X");
    
      leg1->Draw("lt");
      
      TLatex *texb_1 = new TLatex(0.15,0.82,"Efficiency Width");
      texb_1->SetTextSize(0.1);
      texb_1->SetTextFont(42);
      texb_1->Draw("same");
      
      TLatex *texb_3 = new TLatex(0.15,0.7, Form("AuAu #sqrt{s_{NN}} = %s GeV", exactEnergies[energyIdx]));
      texb_3->SetTextSize(0.08);
      texb_3->SetTextFont(42);
      texb_3->Draw("same");
      
      TLatex *texb_4 = new TLatex(0.15,0.6,"|#eta| < 0.5");
      texb_4->SetTextSize(0.08);
      texb_4->SetTextFont(42);
      texb_4->Draw("same");

      TLatex *texb_4a = new TLatex(0.15,0.5,Form("%.1f < #it{p}_{T} < %.1f (GeV/#it{c})", namesPtMin[idx], namesPtMax[idx]));
      texb_4a->SetTextSize(0.08);
      texb_4a->SetTextFont(42);
      texb_4a->Draw("same");
      
      TLatex *texb_4b = new TLatex(0.15,0.4,"sample size 20 events each");
      texb_4b->SetTextSize(0.07);
      texb_4b->SetTextFont(42);
      texb_4b->Draw("same");

      pad->Modified();
      canProfile[idx][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.73,0.06,"Efficiency");
      texb_5->SetTextSize(0.03);
      texb_5->SetTextFont(42);
      texb_5->Draw("same");
      
      TLatex *texb_6 = new TLatex(0.03,0.4, Form("Sample Counts - %s", names[idx]));
      texb_6->SetTextSize(0.03);
      texb_6->SetTextFont(42);
      texb_6->SetTextAngle(90);
      texb_6->Draw("same");

      pad->Modified();
      canProfile[idx][energyIdx]->cd();

    } // end idx
  } // end energy idx

  // ----------------------------------------------------------
  // -- Write out canvas
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      can[idx][energyIdx]->SaveAs(Form("./results/effStudy/root/effWidth_%s_%sGeV.root",    names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/effStudy/root_macro/effWidth_%s_%sGeV.C", names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/effStudy/png/effWidth_%s_%sGeV.png",      names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/effStudy/pdf/effWidth_%s_%sGeV.pdf",      names[idx], energies[energyIdx]));

      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/root/effProfile_%s_%sGeV.root",    names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/root_macro/effProfile_%s_%sGeV.C", names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/png/effProfile_%s_%sGeV.png",      names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/pdf/effProfile_%s_%sGeV.pdf",      names[idx], energies[energyIdx]));
    }

    for (int idx = 0; idx < 2*nNames; ++idx) {
      int idxName    = idx/2;
      int isNegative = idx%2;
      if (isNegative) {
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root/pt_eff_2D_%sNegative_%sGeV.root",    names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root_macro/pt_eff_2D_%sNegative_%sGeV.C", names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/png/pt_eff_2D_%sNegative_%sGeV.png",      names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/pdf/pt_eff_2D_%sNegative_%sGeV.pdf",      names[idxName], energies[energyIdx]));
      }
      else {
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root/pt_eff_2D_%sPositive_%sGeV.root",    names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root_macro/pt_eff_2D_%sPositive_%sGeV.C", names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/png/pt_eff_2D_%sPositive_%sGeV.png",      names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/pdf/pt_eff_2D_%sPositive_%sGeV.pdf",      names[idxName], energies[energyIdx]));
      }
    }
  }
}





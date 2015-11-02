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
const char* Names[]  = {"Pion", "Kaon", "Proton"};
const char* names2[] = {"#pi^{+}", "K^{+}", "p"};
const char* names3[] = {"#pi^{-}", "K^{-}", "#bar{p}"};

#if 1
const float namesPtMin[]  = {0.7, 1.0, 0.8};
const float namesPtMax[]  = {2.3, 2.3, 2.3};

const int namesPtBinMin[]  = { 8, 11,  9};
const int namesPtBinMax[]  = {15, 15, 15};
#else
const float namesPtMin[]  = {1.1, 1.1, 1.1};
const float namesPtMax[]  = {1.4, 1.4, 1.4};

const int namesPtBinMin[]  = {12, 12, 12};
const int namesPtBinMax[]  = {12, 12, 12};
#endif

// 1  -> 0   > 0.05 < 0.1
// 2  -> 0.1 > 0.15 < 0.2
// 3  -> 0.2 > 0.25 < 0.3
// 4  -> 0.3 > 0.35 < 0.4
// 5  -> 0.4 > 0.45 < 0.5
// 6  -> 0.5 > 0.55 < 0.6
// 7  -> 0.6 > 0.65 < 0.7
// 8  -> 0.7 > 0.75 < 0.8
// 9  -> 0.8 > 0.85 < 0.9
// 10 -> 0.9 > 0.95 < 1
// 11 -> 1   > 1.05 < 1.1
// 12 -> 1.1 > 1.25 < 1.4
// 13 -> 1.4 > 1.55 < 1.7
// 14 -> 1.7 > 1.85 < 2
// 15 -> 2   > 2.15 < 2.3

Double_t xbinsProton[4] = {0, 0.2, 0.8, 2.3}; 
Double_t xbinsPion[4]   = {0, 0.2, 0.7, 2.3}; 
Double_t xbinsKaon[4]   = {0, 0.2, 1.0, 2.3}; 


const int   nEnergies       = 8;
const char *energies[]      = {  "7",   "11",   "14",   "19", "27", "39",   "62", "200"};
const char *exactEnergies[] = {"7.7", "11.5", "14.5", "19.6", "27", "39", "62.4", "200"};

double plateau(double *x,double *par) {
  return par[0]*x[0]+par[1];
}

// __________________________________________________________________________________
void makeEfficiencyStudy() {
  gROOT->Reset();
  gROOT->LoadMacro("./setupStyle.C");
  setupStyle();

  gSystem->Exec("mkdir -p ./results/effStudy/png  ./results/effStudy/pdf ./results/effStudy/eps  ./results/effStudy/root  ./results/effStudy/root_macro");
  gSystem->Exec("mkdir -p ./results/effStudy/profile/png  ./results/effStudy/profile/pdf ./results/effStudy/profile/eps  ./results/effStudy/profile/root  ./results/effStudy/profile/root_macro");
  gSystem->Exec("mkdir -p ./results/effStudy/inset/png  ./results/effStudy/inset/pdf ./results/effStudy/inset/eps  ./results/effStudy/inset/root  ./results/effStudy/inset/root_macro");
  gSystem->Exec("mkdir -p ./results/effStudy/support/png  ./results/effStudy/support/pdf ./results/effStudy/support/eps  ./results/effStudy/support/root  ./results/effStudy/support/root_macro");

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);

  TFile *fplus[3][nEnergies];
  TFile *fminus[3][nEnergies];

  TCanvas* can[nNames][nEnergies];
  TCanvas* can2D[2*nNames][nEnergies];
  TCanvas* canProfile[nNames][nEnergies]; 
  TCanvas* canInset[nEnergies];

  TProfile* histsPlusEff[nNames][nEnergies][nCent];
  TProfile* histsMinusEff[nNames][nEnergies][nCent];

  TProfile* histsPlusEffRebin[nNames][nEnergies][nCent];
  TProfile* histsMinusEffRebin[nNames][nEnergies][nCent];


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
	histsPlus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fplus[idx][energyIdx]->Get(Form("hpt_widthSmeared_%s",cent[centIdx-1])));
	histsPlus[idx][energyIdx][centIdx]->SetName(Form("hpt_width_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));
	histsPlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	
	histsMinus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fminus[idx][energyIdx]->Get(Form("hpt_widthSmeared_%s",cent[centIdx-1])));
	histsMinus[idx][energyIdx][centIdx]->SetName(Form("hpt_width_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));	
	histsMinus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);


	histsPlusEff[idx][energyIdx][centIdx] = static_cast<TProfile*>(fplus[idx][energyIdx]->Get(Form("effProfileSmeared_%s",cent[centIdx-1])));
	histsPlusEff[idx][energyIdx][centIdx]->SetName(Form("effProfileSmeared_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));
	histsPlusEff[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsPlusEff[idx][energyIdx][centIdx]->SetMarkerStyle(24);
	histsPlusEff[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsPlusEff[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsPlusEff[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsPlusEff[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
       
	if (idx == 0)
	  histsPlusEffRebin[idx][energyIdx][centIdx] = static_cast<TProfile*>(histsPlusEff[idx][energyIdx][centIdx]->Rebin(3, Form("effProfileSmearedRebin_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]), xbinsPion));
	else if (idx == 1)
	  histsPlusEffRebin[idx][energyIdx][centIdx] = static_cast<TProfile*>(histsPlusEff[idx][energyIdx][centIdx]->Rebin(3, Form("effProfileSmearedRebin_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]), xbinsKaon));
	else if (idx == 2)
	  histsPlusEffRebin[idx][energyIdx][centIdx] = static_cast<TProfile*>(histsPlusEff[idx][energyIdx][centIdx]->Rebin(3, Form("effProfileSmearedRebin_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]), xbinsProton));

	histsMinusEff[idx][energyIdx][centIdx] = static_cast<TProfile*>(fminus[idx][energyIdx]->Get(Form("effProfileSmeared_%s",cent[centIdx-1])));
	histsMinusEff[idx][energyIdx][centIdx]->SetName(Form("effProfileSmeared_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));	
	histsMinusEff[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsMinusEff[idx][energyIdx][centIdx]->SetMarkerStyle(25);
	histsMinusEff[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsMinusEff[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsMinusEff[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsMinusEff[idx][energyIdx][centIdx]->SetLineColor(kAzure);
	histsMinusEff[idx][energyIdx][centIdx]->GetXaxis()->SetRangeUser(0.8, 2.3);
	histsMinusEff[idx][energyIdx][centIdx]->Sumw2();

	if (idx == 0)
	  histsMinusEffRebin[idx][energyIdx][centIdx] = static_cast<TProfile*>(histsMinusEff[idx][energyIdx][centIdx]->Rebin(3, Form("effProfileSmearedRebin_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]), xbinsPion));
	else if (idx == 1)
	  histsMinusEffRebin[idx][energyIdx][centIdx] = static_cast<TProfile*>(histsMinusEff[idx][energyIdx][centIdx]->Rebin(3, Form("effProfileSmearedRebin_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]), xbinsKaon));
	else if (idx == 2)
	  histsMinusEffRebin[idx][energyIdx][centIdx] = static_cast<TProfile*>(histsMinusEff[idx][energyIdx][centIdx]->Rebin(3, Form("effProfileSmearedRebin_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]), xbinsProton));

	hists2DPlus[idx][energyIdx][centIdx] = static_cast<TH2D*>(fplus[idx][energyIdx]->Get(Form("hpt_2D_%s",cent[centIdx-1])));
	hists2DPlus[idx][energyIdx][centIdx]->SetName(Form("hpt_2D_%s_plus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));			

	hists2DMinus[idx][energyIdx][centIdx] = static_cast<TH2D*>(fminus[idx][energyIdx]->Get(Form("hpt_2D_%s",cent[centIdx-1])));
	hists2DMinus[idx][energyIdx][centIdx]->SetName(Form("hpt_2D_%s_minus_%s_%s", names[idx], energies[energyIdx],cent[centIdx-1]));		
      }
    }
  }

  // ----------------------------------------------------------
  // -- print RMS of TProfile for high pT
  // ----------------------------------------------------------

  for (int idx = 0; idx < nNames; ++idx) {
    for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
      cout << names[idx] << " | " << exactEnergies[energyIdx] << " || ";
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	if (histsMinusEffRebin[idx][energyIdx][centIdx])
	  printf("%.3f ", histsMinusEffRebin[idx][energyIdx][centIdx]->GetBinError(3));

	  //	  cout << histsMinusEffRebin[idx][energyIdx][centIdx]->GetBinError(3) << " ";
	//	  cout << histsMinusEffRebin[idx][energyIdx][centIdx]->GetRMS(3) << " ";
      }
      cout << endl;
    }
  }
  return;
  // ----------------------------------------------------------
  // -- Make Profiles
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	histsProfilePlus[idx][energyIdx][centIdx] = 
	  hists2DPlus[idx][energyIdx][centIdx]->ProjectionY(Form("%s_profile", hists2DPlus[idx][energyIdx][centIdx]->GetName()), namesPtBinMin[idx], namesPtBinMax[idx]);
	histsProfilePlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsProfilePlus[idx][energyIdx][centIdx]->SetMarkerStyle(24);
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
	histsProfileMinus[idx][energyIdx][centIdx]->SetMarkerStyle(25);
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
      pad->SetFillColor(0);
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
	
	TH2D *ff = new TH2D("","",20,0.009,2.09,20,0.00,0.15);
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
      
      TLatex *texb_4a = new TLatex(0.15,0.5,"sample size 100 events each");
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
      pad->SetFillColor(0);
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
      
      TLatex *texb_4a = new TLatex(0.15,0.5,"sample size 100 events each");
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
      pad->SetFillColor(0);
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
	//	gPad->SetLogy();

	if (centIdx == 5) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 9) {
	  gPad->SetRightMargin(0.003);
	}
	
	histsProfilePlus[idx][energyIdx][centIdx]->Draw("e1,same");
	histsProfileMinus[idx][energyIdx][centIdx]->DrawCopy("e1,same");
	
	TLatex *texb_Cent = new TLatex(0.05, 500, cent1[centIdx-1]);
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
      
      TLatex *texb_4b = new TLatex(0.15,0.4,"sample size 100 events each");
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
  // -- Fill canvas by name / energy -- eff width
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx <nEnergies; ++energyIdx) {
    canInset[energyIdx] = new TCanvas(Form("canInset_%s", energies[energyIdx]), Form("canInset_%s", energies[energyIdx]), 0, 0, 1400, 700);
    canInset[energyIdx]->SetFillColor(0);
    canInset[energyIdx]->SetBorderMode(0);
    canInset[energyIdx]->SetBorderSize(0.0);
    canInset[energyIdx]->SetFrameFillColor(0);
    canInset[energyIdx]->SetFrameBorderMode(0);
    canInset[energyIdx]->cd();

    int centIdx = 1; 

    TPad* padContent = new TPad("padContent", "", 0.05, 0.1, 0.99, 0.99);
    padContent->SetBorderMode(0);
    padContent->SetFillColor(0);
    padContent->Draw();
    padContent->cd();
    padContent->Divide(3,1,0.,0.,0.);
    
    for (int idx = 0; idx < nNames; ++idx) {
      
      TLegend *leg1= new TLegend(0.75, 0.85, 1.2, 0.95);
      leg1->SetName(Form("leg1_%s", names[idx]));
      leg1->SetTextAlign(12);
      leg1->SetTextSize(0.06);
      leg1->SetTextFont(42);
      leg1->SetFillColor(1182);
      leg1->SetLineColor(0);
      leg1->SetBorderSize(0);
      
      leg1->AddEntry(histsPlusEff[idx][energyIdx][centIdx],  names2[idx],"ple");
      leg1->AddEntry((TObject*)0, "", "");
      leg1->AddEntry(histsMinusEff[idx][energyIdx][centIdx], names3[idx],"ple");
           
      padContent->cd(idx+1);
      
      double yMax = 1.18;
      TH2D *ff = new TH2D("","",20,0.009,2.09,20,0.00,yMax);
      ff->GetXaxis()->SetLabelSize(0.06);
      ff->GetYaxis()->SetLabelSize(0.06);
      ff->GetXaxis()->SetNdivisions(9,5,0);
      ff->Draw();
      
      yMax = 0.95;
      TLine *line02 = new TLine(namesPtMin[idx], 0, namesPtMin[idx], yMax);
      line02->SetLineColor(kGray+4);
      line02->SetLineStyle(3);
      line02->Draw();
      
      TLine *line20 = new TLine(namesPtMax[idx], 0, namesPtMax[idx], yMax);
      line20->SetLineColor(kGray+4);
      line20->SetLineStyle(3);
      line20->Draw();


      yMax = 1.18;

      TLine *line02a = new TLine(namesPtMin[idx], 1.12, namesPtMin[idx], yMax);
      line02a->SetLineColor(kGray+4);
      line02a->SetLineStyle(3);
      line02a->Draw();

      TLine *line20a = new TLine(namesPtMax[idx], 1.14, namesPtMax[idx], yMax);
      line20a->SetLineColor(kGray+4);
      line20a->SetLineStyle(3);
      line20a->Draw();



      histsPlusEff[idx][energyIdx][centIdx]->SetMarkerStyle(25);
      histsPlusEff[idx][energyIdx][centIdx]->SetMarkerSize(0.7);
      
      histsMinusEff[idx][energyIdx][centIdx]->SetMarkerStyle(24);
      histsMinusEff[idx][energyIdx][centIdx]->SetMarkerSize(0.7);

      histsPlusEff[idx][energyIdx][centIdx]->Draw("p,same");
      histsMinusEff[idx][energyIdx][centIdx]->Draw("p,same");
	
      leg1->Draw();

      TPad* padInset = new TPad("padInset", "padInset", 0.25, 0.13, 0.95, 0.50);
      padInset->SetBorderMode(0);
      padInset->SetFillColor(0);
      padInset->Draw();
      padInset->cd();

      TPad* padInsetContent = new TPad("padInsetContent", "", 0.05, 0.1, 0.99, 0.99);
      padInsetContent->SetRightMargin(0.002);
      padInsetContent->SetBorderMode(0);
      padInsetContent->SetFillColor(1182);
      padInsetContent->Draw();
      padInsetContent->cd();

      //      gPad->SetLogy();

      histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetRangeUser(0.15, 1.1);
      histsProfilePlus[idx][energyIdx][centIdx]->GetYaxis()->SetRangeUser(1, 300);

      histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetLabelSize(0.09);
      histsProfilePlus[idx][energyIdx][centIdx]->GetYaxis()->SetLabelSize(0.09);
      histsProfilePlus[idx][energyIdx][centIdx]->GetXaxis()->SetNdivisions(5,5,0);

      histsProfilePlus[idx][energyIdx][centIdx]->SetMarkerStyle(25);
      histsProfilePlus[idx][energyIdx][centIdx]->SetMarkerSize(0.7);

      histsProfileMinus[idx][energyIdx][centIdx]->SetMarkerStyle(24);
      histsProfileMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.7);

      histsProfilePlus[idx][energyIdx][centIdx]->Draw("e1");
      histsProfileMinus[idx][energyIdx][centIdx]->Draw("e1,same");
      
      TLatex *texbInset_3 = new TLatex(0.40, 80, Form("%0.1f < #it{p}_{T} (GeV/#it{c}) < %0.1f", namesPtMin[idx], namesPtMax[idx]));
      texbInset_3->SetTextSize(0.09);
      texbInset_3->SetTextFont(42);
      texbInset_3->Draw("same");

      padInset->cd();
      TLatex *texbInset_5 = new TLatex(0.47, 0.05,"Efficiency");
      texbInset_5->SetTextSize(0.09);
      texbInset_5->SetTextFont(42);
      texbInset_5->Draw("same");
      
      TLatex *texbInset_6 = new TLatex(0.05, 0.18, "Counts / 100 Events");
      texbInset_6->SetTextSize(0.09);
      texbInset_6->SetTextFont(42);
      texbInset_6->SetTextAngle(90);
      texbInset_6->Draw("same");
    } // end idx
    
    padContent->Modified();
    canInset[energyIdx]->cd();
    
    padContent->cd(1);
    TLatex *texb_3 = new TLatex(0.10, 1.08, Form("AuAu #sqrt{s_{NN}} = %s GeV", exactEnergies[energyIdx]));
    texb_3->SetTextSize(0.06);
    texb_3->SetTextFont(42);
    texb_3->Draw("same");
      
    for (int idx = 0; idx < nNames; ++idx) {
      padContent->cd(idx+1);
      TLatex *texb_4 = new TLatex(0.10,0.99,Form("%s, |#eta| < 0.5, %s", Names[idx], cent1[centIdx-1]));
      texb_4->SetTextSize(0.05);
      texb_4->SetTextFont(42);
      texb_4->Draw("same");
    }

    padContent->cd(2);
    TLatex *texb_2a = new TLatex(0.10,1.08, "vertical error bars represent RMS");
    texb_2a->SetTextSize(0.04);
    texb_2a->SetTextFont(42);
    texb_2a->Draw("same");
    
    padContent->cd(3);
    TLatex *texb_3a = new TLatex(0.10,1.08, "STAR Preliminary");
    texb_3a->SetTextSize(0.07);
    texb_3a->SetTextFont(42);
    texb_3a->Draw("same");
    
    padContent->Modified();
    canInset[energyIdx]->cd();
      
    TLatex *texb_5 = new TLatex(0.47, 0.07,"#it{p}_{T} (GeV/#it{c})");
    texb_5->SetTextSize(0.04);
    texb_5->SetTextFont(42);
    texb_5->Draw("same");
    
    TLatex *texb_6 = new TLatex(0.03,0.4, "Tracking Efficiency");
    texb_6->SetTextSize(0.04);
    texb_6->SetTextFont(42);
    texb_6->SetTextAngle(90);
    texb_6->Draw("same");
    
    padContent->Modified();
    canInset[energyIdx]->cd();
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
      can[idx][energyIdx]->SaveAs(Form("./results/effStudy/eps/effWidth_%s_%sGeV.eps",      names[idx], energies[energyIdx]));

      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/root/effProfile_%s_%sGeV.root",    names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/root_macro/effProfile_%s_%sGeV.C", names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/png/effProfile_%s_%sGeV.png",      names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/pdf/effProfile_%s_%sGeV.pdf",      names[idx], energies[energyIdx]));
      canProfile[idx][energyIdx]->SaveAs(Form("./results/effStudy/profile/eps/effProfile_%s_%sGeV.eps",      names[idx], energies[energyIdx]));
    }

    canInset[energyIdx]->SaveAs(Form("./results/effStudy/inset/root/effInset_%sGeV.root",    energies[energyIdx]));
    canInset[energyIdx]->SaveAs(Form("./results/effStudy/inset/root_macro/effInset_%sGeV.C", energies[energyIdx]));
    canInset[energyIdx]->SaveAs(Form("./results/effStudy/inset/png/effInset_%sGeV.png",      energies[energyIdx]));
    canInset[energyIdx]->SaveAs(Form("./results/effStudy/inset/pdf/effInset_%sGeV.pdf",      energies[energyIdx]));
    canInset[energyIdx]->SaveAs(Form("./results/effStudy/inset/eps/effInset_%sGeV.eps",      energies[energyIdx]));
    
    for (int idx = 0; idx < 2*nNames; ++idx) {
      int idxName    = idx/2;
      int isNegative = idx%2;
      if (isNegative) {
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root/pt_eff_2D_%sNegative_%sGeV.root",    names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root_macro/pt_eff_2D_%sNegative_%sGeV.C", names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/png/pt_eff_2D_%sNegative_%sGeV.png",      names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/pdf/pt_eff_2D_%sNegative_%sGeV.pdf",      names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/eps/pt_eff_2D_%sNegative_%sGeV.eps",      names[idxName], energies[energyIdx]));
      }
      else {
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root/pt_eff_2D_%sPositive_%sGeV.root",    names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/root_macro/pt_eff_2D_%sPositive_%sGeV.C", names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/png/pt_eff_2D_%sPositive_%sGeV.png",      names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/pdf/pt_eff_2D_%sPositive_%sGeV.pdf",      names[idxName], energies[energyIdx]));
	can2D[idx][energyIdx]->SaveAs(Form("./results/effStudy/support/eps/pt_eff_2D_%sPositive_%sGeV.eps",      names[idxName], energies[energyIdx]));
      }
    }
  }
}





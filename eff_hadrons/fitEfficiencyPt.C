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

const int   nEnergies       = 8;
const float xPosLabel[]     = { 0.08,   0.19,   0.30,   0.40,  0.51,    0.62,  0.73,   0.84};
const char *energies[]      = {  "7",   "11",   "14",   "19",   "27",   "39",   "62", "200"};
const char *exactEnergies[] = {"7.7", "11.5", "14.5", "19.6", "27.0", "39.0", "62.4", "200"};

const double fitRanges[2]        = {0.1, 3.};  // {0.2, 1.95};
const double fitRangesDelta[2]   = {0.1, 4.5}; 
const double fitRangesPlateau[2] = {1.,  3.};  // {1., 1.95};

double tpc(double *x,double *par) {
  return par[0]*TMath::Exp(-1*pow(par[1]/x[0],par[2]));
}

double plateau(double *x,double *par) {
  return par[0]*x[0]+par[1];
}

double deltaPt(double *x,double *par) {
  return par[0]+par[1]*exp(-pow(x[0]/par[2],par[3]))+par[4]*x[0];
}

// __________________________________________________________________________________
void fitEfficiencyPt() {
  gROOT->Reset();
  gROOT->LoadMacro("./setupStyle.C");
  setupStyle();

  gSystem->Exec("mkdir -p ./results/fits");
  gSystem->Exec("mkdir -p ./results/particles_pt/png  ./results/particles_pt/pdf ./results/particles_pt/eps  ./results/particles_pt/root  ./results/particles_pt/root_macro");
  gSystem->Exec("mkdir -p ./results/particles_deltaPt/png  ./results/particles_deltaPt/pdf ./results/particles_deltaPt/eps  ./results/particles_deltaPt/root  ./results/particles_deltaPt/root_macro");
  gSystem->Exec("mkdir -p ./results/energy_cmp_pt/png ./results/energy_cmp_pt/pdf ./results/energy_cmp_pt/eps ./results/energy_cmp_pt/root ./results/energy_cmp_pt/root_macro");

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);

  TFile *fplus[nNames][nEnergies];
  TFile *fminus[nNames][nEnergies];

  TCanvas* can[nNames][nEnergies];
  TCanvas* canDelta[nNames][nEnergies];

  TCanvas* canCent[nCent];
 
  TH1D* histsPlus[nNames][nEnergies][nCent];
  TH1D* histsMinus[nNames][nEnergies][nCent];

  TH1D* histsDeltaPlus[nNames][nEnergies][nCent];
  TH1D* histsDeltaMinus[nNames][nEnergies][nCent];

  TF1* funPlus[nNames][nEnergies][nCent];
  TF1* funMinus[nNames][nEnergies][nCent];

  TF1* funDeltaPlus[nNames][nEnergies][nCent];
  TF1* funDeltaMinus[nNames][nEnergies][nCent];

  TF1* funPlateauPlus[nNames][nEnergies][nCent];
  TF1* funPlateauMinus[nNames][nEnergies][nCent];

  // ----------------------------------------------------------
  // -- read files
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    fplus[0][energyIdx]  = TFile::Open(Form("efficiency/piplus%sGeV.root",  energies[energyIdx]));
    fminus[0][energyIdx] = TFile::Open(Form("efficiency/piminus%sGeV.root", energies[energyIdx]));
    
    fplus[1][energyIdx]  = TFile::Open(Form("efficiency/kplus%sGeV.root",   energies[energyIdx]));
    fminus[1][energyIdx] = TFile::Open(Form("efficiency/kminus%sGeV.root",  energies[energyIdx]));
    
    fplus[2][energyIdx]  = TFile::Open(Form("efficiency/pplus%sGeV.root",   energies[energyIdx]));
    fminus[2][energyIdx] = TFile::Open(Form("efficiency/pminus%sGeV.root",  energies[energyIdx]));
  }

  // ----------------------------------------------------------
  // -- read histograms
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	histsPlus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fplus[idx][energyIdx]->Get(Form("hpt_%s",cent[centIdx-1])));
	histsPlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	
	histsMinus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fminus[idx][energyIdx]->Get(Form("hpt_%s",cent[centIdx-1])));	
	histsMinus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);

	histsDeltaPlus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fplus[idx][energyIdx]->Get(Form("hdeltapt_%s",cent[centIdx-1])));		
	histsDeltaPlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsDeltaPlus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsDeltaPlus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsDeltaPlus[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsDeltaPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsDeltaPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);

	histsDeltaMinus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fminus[idx][energyIdx]->Get(Form("hdeltapt_%s",cent[centIdx-1])));
	histsDeltaMinus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsDeltaMinus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsDeltaMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsDeltaMinus[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsDeltaMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsDeltaMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
      }
    }
  }

  // ----------------------------------------------------------
  // -- Fit histograms
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {

	funPlateauPlus[idx][energyIdx][centIdx] = new TF1(Form("funEffPlateau_Plus_%s_%s_%s", names[idx], energies[energyIdx], cent[centIdx-1]), 
							  plateau, 0.1, 3, 3);
	funPlateauPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	funPlateauPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funPlateauPlus[idx][energyIdx][centIdx]->SetParameters(0.1, 0.7);

	histsPlus[idx][energyIdx][centIdx]->Fit(funPlateauPlus[idx][energyIdx][centIdx], "R", "", fitRangesPlateau[0], fitRangesPlateau[1]);

	// ----------------------------------------------------------
	
	funPlus[idx][energyIdx][centIdx] = new TF1(Form("funEff_Plus_%s_%s_%s", names[idx], energies[energyIdx], cent[centIdx-1]),  
						   tpc, 0.1, 3, 3);
	funPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	funPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funPlus[idx][energyIdx][centIdx]->SetParameters(0.7, 0.2, 2.0);

	histsPlus[idx][energyIdx][centIdx]->Fit(funPlus[idx][energyIdx][centIdx], "R", "", fitRanges[0], fitRanges[1]);

	// ----------------------------------------------------------
	
	funDeltaPlus[idx][energyIdx][centIdx] = new TF1(Form("funDeltaPt_Plus_%s_%s_%s", names[idx], energies[energyIdx], cent[centIdx-1]),  
	 						deltaPt, 0.1, 5., 5);
	funDeltaPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	funDeltaPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funDeltaPlus[idx][energyIdx][centIdx]->SetParameters(0.0016, 0.26, 0.086, 0.82, -0.0015);

	histsDeltaPlus[idx][energyIdx][centIdx]->Fit(funDeltaPlus[idx][energyIdx][centIdx], "ERsame", "", fitRangesDelta[0], fitRangesDelta[1]);
	funDeltaPlus[idx][energyIdx][centIdx]->SetRange(0.1,10);

	// ----------------------------------------------------------

	funPlateauMinus[idx][energyIdx][centIdx] = new TF1(Form("funEffPlateau_Minus_%s_%s_%s", names[idx], energies[energyIdx], cent[centIdx-1]), 
							   plateau, 0.1, 3, 3);
	funPlateauMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
	funPlateauMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funPlateauMinus[idx][energyIdx][centIdx]->SetParameters(0.1, 0.7);

	histsMinus[idx][energyIdx][centIdx]->Fit(funPlateauMinus[idx][energyIdx][centIdx], "R", "", fitRangesPlateau[0], fitRangesPlateau[1]);

	// ----------------------------------------------------------

	funMinus[idx][energyIdx][centIdx] = new TF1(Form("funEff_Minus_%s_%s_%s", names[idx], energies[energyIdx], cent[centIdx-1]), 
						    tpc, 0.1, 3, 3);
	funMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
	funMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funMinus[idx][energyIdx][centIdx]->SetParameters(0.7, 0.2, 2.0);
	
	histsMinus[idx][energyIdx][centIdx]->Fit(funMinus[idx][energyIdx][centIdx], "R", "", fitRanges[0], fitRanges[1]);
	
	// ----------------------------------------------------------
	
	funDeltaMinus[idx][energyIdx][centIdx] = new TF1(Form("funDeltaPt_Minus_%s_%s_%s", names[idx], energies[energyIdx], cent[centIdx-1]),  
							 deltaPt, 0.1, 5., 5);
	funDeltaMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
	funDeltaMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funDeltaMinus[idx][energyIdx][centIdx]->SetParameters(0.0016, 0.26, 0.086, 0.82, -0.0015);
	
	histsDeltaMinus[idx][energyIdx][centIdx]->Fit(funDeltaMinus[idx][energyIdx][centIdx], "ERsame", "", fitRangesDelta[0], fitRangesDelta[1]);
	funDeltaMinus[idx][energyIdx][centIdx]->SetRange(0.1,10);
      }
    }
  }

  // ----------------------------------------------------------
  // -- Fill canvas by name / energy -- eff
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
      
      TLegend *leg1= new TLegend(0.15,0.35,0.8,0.50);
      leg1->SetName(Form("leg1_%s", names[idx]));
      leg1->SetTextAlign(12);
      leg1->SetTextSize(0.10);
      leg1->SetTextFont(42);
      leg1->SetFillColor(kWhite);
      leg1->SetLineColor(0);
      leg1->SetBorderSize(0);
      
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	//for (int centIdx = nCent-1; centIdx > 0; centIdx--) {
	pad->cd(centIdx);
	
	if (centIdx == 5) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 9) {
	  gPad->SetRightMargin(0.003);
	}
	
	TH2D *ff = new TH2D("","",20,0.009,2.09,20,0.01,0.99);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.07);
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

	histsPlus[idx][energyIdx][centIdx]->Draw("psame");
	histsMinus[idx][energyIdx][centIdx]->Draw("psame");
	
	TLatex *texb_Cent = new TLatex(0.8,0.1,cent1[centIdx-1]);
	texb_Cent->SetTextSize(0.10);
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
      
      TLatex *texb_1 = new TLatex(0.15,0.82,"Tracking Efficiency");
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
      
      pad->Modified();
      can[idx][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.73,0.06,"#it{p}_{T} (GeV/#it{c})");
      texb_5->SetTextSize(0.03);
      texb_5->SetTextFont(42);
      texb_5->Draw("same");
      
      TLatex *texb_6 = new TLatex(0.03,0.4, Form("Efficiency - %s", names[idx]));
      texb_6->SetTextSize(0.03);
      texb_6->SetTextFont(42);
      texb_6->SetTextAngle(90);
      texb_6->Draw("same");

      pad->Modified();
      can[idx][energyIdx]->cd();

    } // end idx
  } // end energy idx

  // ----------------------------------------------------------
  // -- Fill canvas for Delta pt by name / energy 
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      canDelta[idx][energyIdx] = new TCanvas(Form("canDelta_%s_%s", names[idx], energies[energyIdx]), names[idx],0,0,1200,600);
      canDelta[idx][energyIdx]->SetFillColor(0);
      canDelta[idx][energyIdx]->SetBorderMode(0);
      canDelta[idx][energyIdx]->SetBorderSize(0.0);
      canDelta[idx][energyIdx]->SetFrameFillColor(0);
      canDelta[idx][energyIdx]->SetFrameBorderMode(0);
      canDelta[idx][energyIdx]->cd();

      TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
      pad->SetBorderMode(0);
      pad->SetFillColor(0);
      pad->Draw();
      pad->cd();
      pad->Divide(5,2,0.,0.,0.);
      
      TLegend *leg1= new TLegend(0.15,0.35,0.8,0.50);
      leg1->SetName(Form("legDelta1_%s", names[idx]));
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
	
	TH2D *ff = new TH2D("","",20, 0.009, 2.09, 21,-0.015, 0.05);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.06);
	ff->GetXaxis()->SetNdivisions(9,5,0);
	ff->Draw();

	TLine *line02 = new TLine(0.2, -0.015, 0.2, 0.05);
	line02->SetLineColor(kGray+4);
	line02->SetLineStyle(3);
	line02->Draw();

	TLine *line20 = new TLine(2., -0.015, 2., 0.05);
	line20->SetLineColor(kGray+4);
	line20->SetLineStyle(3);
	line20->Draw();

	histsDeltaPlus[idx][energyIdx][centIdx]->Draw("psame");
        histsDeltaMinus[idx][energyIdx][centIdx]->Draw("psame");
	
	TLatex *texb_Centx = new TLatex(0.95,0.035,cent1[centIdx-1]);
	texb_Centx->SetTextSize(0.10);
	texb_Centx->SetTextFont(42);
	texb_Centx->Draw("same");
	
	pad->Modified();
	canDelta[idx][energyIdx]->cd();
	
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
      
      TLatex *texb_1 = new TLatex(0.15,0.82,"#delta#it{p}_{T}  (#it{p}_{T,MC}-#it{p}_{T,rec})/#it{p}_{T,MC}");
      texb_1->SetTextSize(0.09);
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
      
      pad->Modified();
      canDelta[idx][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.73,0.06,"#it{p}_{T} (GeV/#it{c})");
      texb_5->SetTextSize(0.03);
      texb_5->SetTextFont(42);
      texb_5->Draw("same");
      
      TLatex *texb_6 = new TLatex(0.03,0.4, Form("#delta #it{p}_{T} (GeV/#it{c}) - %s", names[idx]));
      texb_6->SetTextSize(0.03);
      texb_6->SetTextFont(42);
      texb_6->SetTextAngle(90);
      texb_6->Draw("same");

      pad->Modified();
      canDelta[idx][energyIdx]->cd();
    } // end idx
  } // end energy idx

  // ----------------------------------------------------------
  // -- Fill canvas by cent
  // ----------------------------------------------------------
  TLegend *leg2[nNames];

  for (int centIdx = 1; centIdx < nCent; centIdx++) {
    canCent[centIdx] = new TCanvas(Form("canCent_%d", centIdx), Form("canCent_%d", centIdx),0,0,1400,800);
    canCent[centIdx]->SetFillColor(0);
    canCent[centIdx]->SetBorderMode(0);
    canCent[centIdx]->SetBorderSize(0.0);
    canCent[centIdx]->SetFrameFillColor(0);
    canCent[centIdx]->SetFrameBorderMode(0);
    canCent[centIdx]->cd();

    TPad* pad = new TPad("pad", "pad",0.05,0.08,0.94,0.94);
    pad->SetBorderMode(0);
    pad->SetFillColor(0);
    pad->Draw();
    pad->cd();
    pad->Divide(nEnergies, nNames, 0.,0.,0.);

    for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
      for (int idx = 0; idx < nNames; ++idx) {
	pad->cd(idx*nEnergies+energyIdx+1);
			
	TH2D *ff = new TH2D("","",20,0.009,2.09,20,0.01,0.99);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.07);
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

	histsPlus[idx][energyIdx][centIdx]->Draw("psame");
	histsMinus[idx][energyIdx][centIdx]->Draw("psame");

	if (centIdx == 1) {
	  if (idx != 2)
	    leg2[idx] = new TLegend(0.65,0.15,0.8,0.3);
	  else
	    leg2[idx] = new TLegend(0.65,0.25,0.8,0.4);
	  leg2[idx]->SetName(Form("leg1_%s", names[idx]));
	  leg2[idx]->SetTextAlign(12);
	  leg2[idx]->SetTextSize(0.08 );
	  leg2[idx]->SetTextFont(42);
	  leg2[idx]->SetFillColor(kWhite);
	  leg2[idx]->SetLineColor(0);
	  leg2[idx]->SetBorderSize(0);
      
	  leg2[idx]->AddEntry(histsPlus[idx][energyIdx][centIdx],  names2[idx],"ple");
	  leg2[idx]->AddEntry((TObject*)0, "", "");
	  leg2[idx]->AddEntry(histsMinus[idx][energyIdx][centIdx], names3[idx],"ple");
	}

	leg2[idx]->Draw("lt");
      } // end idx
    } // end energy idx

    pad->Modified();
    canCent[centIdx]->cd();

    TLatex *texb_1 = new TLatex(0.3, 0.97, Form("Tracking Efficiency #pi, K, p - [%s] - |#eta| < 0.5",cent1[centIdx-1]));
    texb_1->SetTextSize(0.03);
    texb_1->SetTextFont(42);
    texb_1->Draw("same");

    for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
      TLatex *texb_3 = new TLatex(xPosLabel[energyIdx], 0.925, Form("AuAu #sqrt{s_{NN}} = %s GeV", exactEnergies[energyIdx]));
      texb_3->SetTextSize(0.015);
      texb_3->SetTextFont(42);
      texb_3->Draw("same");
    }
    
    TLatex *texb_5 = new TLatex(0.45,0.06,"#it{p}_{T} (GeV/#it{c})");
    texb_5->SetTextSize(0.03);
    texb_5->SetTextFont(42);
    texb_5->Draw("same");
      
    TLatex *texb_6a = new TLatex(0.045,0.72, Form("Efficiency - %s", names[0]));
    texb_6a->SetTextSize(0.02);
    texb_6a->SetTextFont(42);
    texb_6a->SetTextAngle(90);
    texb_6a->Draw("same");

    TLatex *texb_6b = new TLatex(0.045,0.46, Form("Efficiency - %s", names[1]));
    texb_6b->SetTextSize(0.02);
    texb_6b->SetTextFont(42);
    texb_6b->SetTextAngle(90);
    texb_6b->Draw("same");

    TLatex *texb_6c = new TLatex(0.045,0.15, Form("Efficiency - %s", names[2]));
    texb_6c->SetTextSize(0.02);
    texb_6c->SetTextFont(42);
    texb_6c->SetTextAngle(90);
    texb_6c->Draw("same");
  } // cent idx

  // ----------------------------------------------------------
  // -- Write out fits
  // ----------------------------------------------------------
  TFile *fOutput = TFile::Open("./results/fits/fit.root", "RECREATE");
  fOutput->cd();
  
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (int idx = 0; idx < nNames; ++idx) 
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	funPlus[idx][energyIdx][centIdx]->Write();
	funMinus[idx][energyIdx][centIdx]->Write();
	funDeltaPlus[idx][energyIdx][centIdx]->Write();
	funDeltaMinus[idx][energyIdx][centIdx]->Write();
      }

  fOutput->Close();

  // ----------------------------------------------------------
  // -- Write out canvas
  // ----------------------------------------------------------
  for (int energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      can[idx][energyIdx]->SaveAs(Form("./results/particles_pt/root/pt_%s_%sGeV.root",    names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_pt/root_macro/pt_%s_%sGeV.C", names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_pt/png/pt_%s_%sGeV.png",      names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_pt/pdf/pt_%s_%sGeV.pdf",      names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_pt/eps/pt_%s_%sGeV.eps",      names[idx], energies[energyIdx]));
 
      canDelta[idx][energyIdx]->SaveAs(Form("./results/particles_deltaPt/root/pt_%s_%sGeV.root",    names[idx], energies[energyIdx]));
      canDelta[idx][energyIdx]->SaveAs(Form("./results/particles_deltaPt/root_macro/pt_%s_%sGeV.C", names[idx], energies[energyIdx]));
      canDelta[idx][energyIdx]->SaveAs(Form("./results/particles_deltaPt/png/pt_%s_%sGeV.png",      names[idx], energies[energyIdx]));
      canDelta[idx][energyIdx]->SaveAs(Form("./results/particles_deltaPt/pdf/pt_%s_%sGeV.pdf",      names[idx], energies[energyIdx]));
      canDelta[idx][energyIdx]->SaveAs(Form("./results/particles_deltaPt/eps/pt_%s_%sGeV.eps",      names[idx], energies[energyIdx]));
    }
  }

  for (int centIdx = 1; centIdx < nCent; centIdx++) {
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_pt/root/pt_cmp_%d.root",    centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_pt/root_macro/pt_cmp_%d.C", centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_pt/png/pt_cmp_%d.png",      centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_pt/pdf/pt_cmp_%d.pdf",      centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_pt/eps/pt_cmp_%d.eps",      centIdx));
  }
}





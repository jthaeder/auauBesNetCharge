/* ***************************************************
 *  Fit efficiencies for identified particles   
 *  from embedding - eta dependent
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
 *        - results/fits/fit_eta.root
 *    - Canvas: 
 *        - results/particles_eta
 *        - results/energy_cmp_eta
 *
 * *************************************************** *
 *  Macro adopted from Stephen Horvat
 *
 *  Latest changes: Jochen Thaeder <jmthader@lbl.gov>
 * 
 * ***************************************************
 *  Run :  root -l -b -q fitEfficiencyEta.C++
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

const int   nCent    = 10;
const char *cent[9]  = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};
const char *cent1[9] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"};

const int   nNames    = 3; 
const char* names[3]  = {"pion", "kaon", "proton"};
const char* names2[3] = {"#pi^{+}", "K^{+}", "p"};
const char* names3[3] = {"#pi^{-}", "K^{-}", "#bar{p}"};

const int   nEnergies        = 3; // 8
const char* energies[3]      = {"11",   "14",   "19"};
const char* exactEnergies[3] = {"11.5", "14.5", "19.6"};
const float xPosLabel[]      = {0.15, 0.45, 0.72};
//const char *energies[]       = {  "7",   "11",   "14",   "19",   "27",   "39",   "62",  "200"};
//const char *exactEnergies[] = {"7.7", "11.5", "11.5", "19.6", "27.0", "39.0", "62.4", "62.4"};

const double fitRanges[2]    = {-0.5, 0.5}; 

double plateau(double *x,double *par) {
  return  par[0]*x[0]+par[1];
}

void fitEfficiencyEta() {
  gROOT->Reset();
  gROOT->LoadMacro("./setupStyle.C");
  setupStyle();

  gSystem->Exec("mkdir -p ./results/fits");
  gSystem->Exec("mkdir -p ./results/particles_eta/png ./results/particles_eta/pdf ./results/particles_eta/root ./results/particles_eta/root_macro");
  gSystem->Exec("mkdir -p ./results/energy_cmp_eta/png ./results/energy_cmp_eta/pdf ./results/energy_cmp_eta/root ./results/energy_cmp_eta/root_macro");

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);

  TFile *fplus[nNames][nEnergies];
  TFile *fminus[nNames][nEnergies];

  TCanvas* can[nNames][nEnergies];
  TCanvas* canCent[nCent];
 
  TH1D* histsPlus[nNames][nEnergies][nCent];
  TH1D* histsMinus[nNames][nEnergies][nCent];

  TF1* funPlus[nNames][nEnergies][nCent];
  TF1* funMinus[nNames][nEnergies][nCent];

  // ----------------------------------------------------------
  // -- read files
  // ----------------------------------------------------------
  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
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
  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
	histsPlus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fplus[idx][energyIdx]->Get(Form("heta_%s",cent[centIdx-1])));
	histsPlus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsPlus[idx][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	histsPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	
	histsMinus[idx][energyIdx][centIdx] = static_cast<TH1D*>(fminus[idx][energyIdx]->Get(Form("heta_%s",cent[centIdx-1])));	
	histsMinus[idx][energyIdx][centIdx]->SetLineStyle(1);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerStyle(20);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerSize(0.5);
	histsMinus[idx][energyIdx][centIdx]->SetMarkerColor(kAzure);
	histsMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	histsMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
      }
    }
  }

  // ----------------------------------------------------------
  // -- Fit histograms
  // ----------------------------------------------------------
  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      for (int centIdx = 1; centIdx < nCent; centIdx++) {

	funPlus[idx][energyIdx][centIdx] = new TF1(Form("funEff_Plus_%d_%d_%d", idx, energyIdx,centIdx), plateau, 0.2, 2, 3);
	funPlus[idx][energyIdx][centIdx]->SetLineColor(kRed+3);
	funPlus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funPlus[idx][energyIdx][centIdx]->SetParameters(0.1, 0.7);

	histsPlus[idx][energyIdx][centIdx]->Fit(funPlus[idx][energyIdx][centIdx], "R", "", fitRanges[0], fitRanges[1]);

	// ----------------------------------------------------------

	funMinus[idx][energyIdx][centIdx] = new TF1(Form("funEff_Minus_%d_%d_%d", idx, energyIdx,centIdx), plateau, 0.2, 2, 3);
	funMinus[idx][energyIdx][centIdx]->SetLineColor(kAzure);
	funMinus[idx][energyIdx][centIdx]->SetLineWidth(1);
	funMinus[idx][energyIdx][centIdx]->SetParameters(0.1, 0.7);

	histsMinus[idx][energyIdx][centIdx]->Fit(funMinus[idx][energyIdx][centIdx], "R", "", fitRanges[0], fitRanges[1]);
      }
    }
  }

  // ----------------------------------------------------------
  // -- Create Canvas
  // ----------------------------------------------------------
  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      can[idx][energyIdx] = new TCanvas(Form("can_eta_%s_%s", names[idx], energies[energyIdx]), names[idx],0,0,1200,600);
      can[idx][energyIdx]->SetFillColor(0);
      can[idx][energyIdx]->SetBorderMode(0);
      can[idx][energyIdx]->SetBorderSize(0.0);
      can[idx][energyIdx]->SetFrameFillColor(0);
      can[idx][energyIdx]->SetFrameBorderMode(0);
      can[idx][energyIdx]->cd();
    }
  }

  // ----------------------------------------------------------
  // -- Fill canvas by name / energy 
  // ----------------------------------------------------------
  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      can[idx][energyIdx]->cd();

      TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
      pad->SetBorderMode(0);
      pad->SetFillColor(1182);
      pad->Draw();
      pad->cd();
      pad->Divide(5,2,0.,0.,0.);
      
      TLegend *leg1= new TLegend(0.35,0.35,0.55,0.50);
      leg1->SetName(Form("leg1_%s", names[idx]));
      leg1->SetTextAlign(12);
      leg1->SetTextSize(0.1);
      leg1->SetTextFont(42);
      leg1->SetFillColor(kWhite);
      leg1->SetLineColor(0);
      leg1->SetBorderSize(0);
      
      for (int centIdx = 1; centIdx < 10; centIdx++) {
	pad->cd(centIdx);
	
	if (centIdx == 5) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 9) {
	  gPad->SetRightMargin(0.003);
	}
	
	TH2D *ff = new TH2D("","",20,-1.1,1.1,20,0.01,0.99);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.07);
	ff->GetXaxis()->SetNdivisions(9,5,0);
	ff->Draw();

	TLine *line05 = new TLine(0.5,0,0.5,1);
	line05->SetLineColor(kGray);
	line05->SetLineStyle(3);
	line05->Draw();

	TLine *line50 = new TLine(-0.5,0,-0.5,1);
	line50->SetLineColor(kGray);
	line50->SetLineStyle(3);
	line50->Draw();

	histsPlus[idx][energyIdx][centIdx]->Draw("psame");
	histsMinus[idx][energyIdx][centIdx]->Draw("psame");
	
	TLatex *texb_Cent = new TLatex(-0.3,0.1,cent1[centIdx-1]);
	texb_Cent->SetTextSize(0.08);
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
      
      TLatex *texb_4 = new TLatex(0.15,0.6,"0.2 < #it{p}_{T} < 2.0 (GeV/#it{c})");
      texb_4->SetTextSize(0.08);
      texb_4->SetTextFont(42);
      texb_4->Draw("same");
      
      pad->Modified();
      can[idx][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.73,0.06,"#eta");
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
  // -- Get Canvas
  // ----------------------------------------------------------

  for (int centIdx = 1; centIdx < nCent; centIdx++) {
    canCent[centIdx] = new TCanvas(Form("canCent_%d", centIdx), Form("canCent_eta_%d", centIdx),0,0,1200,900);
    canCent[centIdx]->SetFillColor(0);
    canCent[centIdx]->SetBorderMode(0);
    canCent[centIdx]->SetBorderSize(0.0);
    canCent[centIdx]->SetFrameFillColor(0);
    canCent[centIdx]->SetFrameBorderMode(0);
    canCent[centIdx]->cd();
  }
  
  // ----------------------------------------------------------
  // -- Fill canvas by cent
  // ----------------------------------------------------------
  
  TLegend *leg2[3] = {NULL, NULL, NULL};

  for (int centIdx = 1; centIdx < nCent; centIdx++) {
    canCent[centIdx]->cd();

    TPad* pad = new TPad("pad", "pad",0.05,0.08,0.94,0.94);
    pad->SetBorderMode(0);
    pad->SetFillColor(1182);
    pad->Draw();
    pad->cd();
    pad->Divide(nEnergies, nNames, 0.,0.,0.);
    
    for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
      for (int idx = 0; idx < nNames; ++idx) {
	pad->cd(idx*3+energyIdx+1);
	
	TH2D *ff = new TH2D("","",20,-1.1,1.1,20,0.01,0.99);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.07);
	ff->GetXaxis()->SetNdivisions(9,5,0);
	ff->Draw();

	TLine *line05 = new TLine(0.5,0,0.5,1);
	line05->SetLineColor(kGray);
	line05->SetLineStyle(3);
	line05->Draw();

	TLine *line50 = new TLine(-0.5,0,-0.5,1);
	line50->SetLineColor(kGray);
	line50->SetLineStyle(3);
	line50->Draw();

	histsPlus[idx][energyIdx][centIdx]->Draw("psame");
	histsMinus[idx][energyIdx][centIdx]->Draw("psame");
	
	if (!leg2[idx]) {
	  leg2[idx] = new TLegend(0.35,0.15,0.55,0.3);
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

    TLatex *texb_1 = new TLatex(0.3, 0.97, Form("Tracking Efficiency #pi, K, p - [%s] - 0.2 < #it{p}_{T} < 2.0 (GeV/#it{c})",cent1[centIdx-1]));
    texb_1->SetTextSize(0.03);
    texb_1->SetTextFont(42);
    texb_1->Draw("same");
    
    for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
      TLatex *texb_3 = new TLatex(xPosLabel[energyIdx], 0.925, Form("AuAu #sqrt{s_{NN}} = %s GeV", exactEnergies[energyIdx]));
      texb_3->SetTextSize(0.02);
      texb_3->SetTextFont(42);
      texb_3->Draw("same");
    }

    TLatex *texb_5 = new TLatex(0.52,0.06,"#eta");
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
  TFile *fOutput = TFile::Open("./results/fits/fit_eta.root", "RECREATE");
  fOutput->cd();

  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) 
    for (int idx = 0; idx < nNames; ++idx) 
      for (int centIdx = 1; centIdx < nCent; centIdx++) {
   	funPlus[idx][energyIdx][centIdx]->Write();
  	funMinus[idx][energyIdx][centIdx]->Write();
      }
      fOutput->Close();
 
  // ----------------------------------------------------------
  // -- Write out canvas
  // ----------------------------------------------------------
  for (int energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
    for (int idx = 0; idx < nNames; ++idx) {
      can[idx][energyIdx]->SaveAs(Form("./results/particles_eta/root/eta_%s_%sGeV.root",    names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_eta/root_macro/eta_%s_%sGeV.C", names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_eta/png/eta_%s_%sGeV.png",      names[idx], energies[energyIdx]));
      can[idx][energyIdx]->SaveAs(Form("./results/particles_eta/pdf/eta_%s_%sGeV.pdf",      names[idx], energies[energyIdx]));
    }
  }

  for (int centIdx = 1; centIdx < nCent; centIdx++) {
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_eta/root/eta_cmp_%d.root",    centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_eta/root_macro/eta_cmp_%d.C", centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_eta/png/eta_cmp_%d.png",      centIdx));
    canCent[centIdx]->SaveAs(Form("./results/energy_cmp_eta/pdf/eta_cmp_%d.pdf",      centIdx));
  }
}





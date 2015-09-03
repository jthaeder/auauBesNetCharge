/* ***************************************************
 *  Make charged efficiencies form identified particles   
 *  from embedding, weighted by corrected spectra
 *
 * *************************************************** *
 *  Input: 
 *   - efficiency fits: fits/fit.root
 *   - spectra hists:   spectra/EvanSpectraMay14.root
 *
 * *************************************************** *
 *  Output:
 *    - Fits: 
 *        - fits/chargedEff.root
 *        - fits/chargedEffInput.root
 *    - Canvas: 
 *        - results/spectra
 *        - results/charged_eff
 *
 * *************************************************** *
 *  Macro adopted from Stephen Horvat
 *
 *  Latest changes: Jochen Thaeder <jmthader@lbl.gov>
 * 
 * ***************************************************
 *  Run :  root -l -b -q makeEfficiencyCharged.C++
 * ***************************************************/

#include <fstream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TPad.h"
#include "TColor.h"
#include "TMath.h"
#include "TLatex.h"
#include "TFitResultPtr.h"

#include "setupStyle.C"

using namespace std;

// -- input spectra naming ...
//    c0 - negative | c1 - positive
//    p0 - pion  | p1 - kaon | p2 - proton
//      cl0.ch1  60-80
//      cl2.ch3  40-60
//      cl4.ch5  20-40
//      cl6.ch6  10-20
//      cl7.ch7  5-10
//      cl8.ch8  0-5

// -- centralities
const int   nCent           = 9;
const char *cent[]          = {   "0005",    "0510",    "1020",    "2030",    "3040",    "4050",    "5060",    "6070",    "7080"};
const char *cent1[]         = {   "0-5%",   "5-10%",  "10-20%",  "20-30%",  "30-40%",  "40-50%",  "50-60%",  "60-70%",  "70-80%"};
const char *spectraCent[]   = {"cl8.ch8", "cl7.ch7", "cl6.ch6", "cl4.ch5", "cl4.ch5", "cl2.ch3", "cl2.ch3", "cl0.ch1", "cl0.ch1"};

// -- refmult2 [nEnergies][nCent]
const double refMult2[9][9] = { {407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},   
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},   
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12}, 
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},   
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12},
				{407, 336.4, 258.3, 178.3, 118.7, 75, 44.2, 24, 12} };

// -- fit particles
const int    nNames         = 3; 
const char  *names[]        = {   "pion",  "kaon",  "proton"};
const double masses[]       = {  0.13957, 0.49367, 0.938272 };
const char  *spectraNames[] = {     "p0",    "p1",      "p2"};
const char  *namesTitle[2][3] = { {"#pi^{-}", "K^{-}", "#bar{p}"},
				  {"#pi^{+}", "K^{+}",       "p"} };
// -- fit particle charges
const int   nCharges           = 2;
const char *chargedModes[]     = {"Minus", "Plus", "All"}; 
const char *spectraParticles[] = {   "c0",   "c1"};
const int   chargedColor[]     = { kAzure, kRed+2, kGray};

// -- energies
const int   nEnergies         = 7;
// const char *energies[]        = {"11",   "14",   "19"};
// const char *exactEnergies[]   = {"11.5", "14.5", "19.6"};
// const char *spectraEnergies[] = {"11.5", "11.5", "19.6"};
const char *energies[]        = {  "7",   "11",   "14",   "19",   "27",   "39",   "62",  "200"};
const char *exactEnergies[]   = {"7.7", "11.5", "14.5", "19.6", "27.0", "39.0", "62.4", "200"};
const char *spectraEnergies[] = {"7.7", "11.5", "11.5", "19.6", "27.0", "39.0", "62.4", "62.4"};

// -- pT range for average efficiencies
const double ptRange[] = {0.2, 2.0};

// -- input fit parameters
const double param[3][2][2]      = { { {1.1,   3.1}, {1.1,   0.2} },      // pi-  | pi+    check again 
				     { {1.1,   3.1}, {1.1,   3.1} },      // K-   | K+
				     { {1.,  100. }, {0.05, 10. } } };  // anti p   | p

// -- fit parameter limits for spectra fits
const double parLimits1[3][2][2] = { { {0.05, 10.}, {0.05, 10.} },    // pi-  | pi+
				     { {1.5,  10.}, {1.5,  10.} },    // K-   | K+
				     { {1.,  100.}, {0.05, 10.} } };  // anti p (.05, 10)  | p

const double parLimits2[3][2][2] = { { {0.09, 10.}, {0.09, 10.} },    // pi-  | pi+
				     { {0.05, 10.}, {0.05, 10.} },    // K-   | K+ 
				     { {1.05, 10.}, {0.05, 10.} } };  // anti p (.05, 0.5)  | p

// -- input efficiencies (input = fits)
TObjArray aEffFits;

// -- fitted efficiencies (input = fits)
TF1* funEff[nCharges][nNames][nEnergies][nCent];

// -- fitted Delta Pt
TF1* funEffDeltaPt[nCharges][nNames][nEnergies][nCent];

// -- spectras (histograms)
TH1D* hSpectra[nCharges][nNames][nEnergies][nCent];

// -- spectras (fits to histgrams)
TF1* funSpectra[nCharges][nNames][nEnergies][nCent];
TObjArray aSpectraFits;

// -- charged spectra
TF1* funSpectraCharged[3][nEnergies][nCent];
TF1* funSpectraChargedIn[3][nEnergies][nCent];

// -- charged weighted efficiencies
TF1* funEffWeightedCharged[3][nEnergies][nCent];

// -- charged efficiencies
TF1* funEffCharged[3][nEnergies][nCent];
TF1* funEffChargedIn[3][nEnergies][nCent];

// -- average charged efficiencies
Double_t avgEffCharged[3][nEnergies][nCent];

// -- average chaged efficiencies vs refMult2
TGraph* gAvgEffChargedRefMult2[3][nEnergies];
TF1*    funAvgEffChargedRefMult2[3][nEnergies];

// -- canvas
TCanvas *canSpectra[nCharges][nEnergies];
TCanvas *canSpectraCharged[nEnergies];
TCanvas *canEffCharged[nEnergies];

// -- index variables
int idxParticle, energyIdx, idx, centIdx;

// -- draw methods
void drawSpectra();
void drawSpectraCharged();
void drawEffCharged();

// -- function sumation / multiplication / division
double spectraChargedNeg(double *x, double *par);
double spectraChargedPos(double *x, double *par);
double spectraChargedAll(double *x, double *par);

double effWeightedChargedNeg(double *x, double *par);
double effWeightedChargedPos(double *x, double *par);
double effWeightedChargedAll(double *x, double *par);

double effChargedNeg(double *x, double *par);
double effChargedPos(double *x, double *par);
double effChargedAll(double *x, double *par);

// __________________________________________________________________
void makeEfficiencyCharged() {
  gROOT->Reset();
  gROOT->LoadMacro("./setupStyle.C");
  setupStyle();

  gSystem->Exec("mkdir -p ./results/chargedEff/png ./results/chargedEff/pdf ./results/chargedEff/root ./results/chargedEff/root_macro");
  gSystem->Exec("mkdir -p ./results/spectra/png    ./results/spectra/pdf    ./results/spectra/root    ./results/spectra/root_macro");
  gSystem->Exec("mkdir -p ./results/fits");

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);

  // ----------------------------------------------------------
  // -- read fits from file
  // ----------------------------------------------------------
  TFile *fitFile = TFile::Open("./results/fits/fit.root");

  for (idxParticle = 0; idxParticle < nCharges; ++idxParticle) 
    for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
      for (idx = 0; idx < nNames; ++idx) 
	for (centIdx = 0; centIdx < nCent; centIdx++) {
	  funEff[idxParticle][idx][energyIdx][centIdx] = static_cast<TF1*>(fitFile->Get(Form("funEff_%s_%s_%s_%s", 
											     chargedModes[idxParticle], names[idx], 
											     energies[energyIdx], cent[centIdx])));
	  aEffFits.Add(funEff[idxParticle][idx][energyIdx][centIdx]);
	  
	  funEffDeltaPt[idxParticle][idx][energyIdx][centIdx] = static_cast<TF1*>(fitFile->Get(Form("funEffDeltaPt_%s_%s_%s_%s", 
												    chargedModes[idxParticle], names[idx], 
												    energies[energyIdx], cent[centIdx])));
	  aEffFits.Add(funEffDeltaPt[idxParticle][idx][energyIdx][centIdx]);

	}

  // ----------------------------------------------------------
  // -- read spectra from file
  // ----------------------------------------------------------
  TFile *spectraFile = TFile::Open("spectra/EvanSpectraMay14.root");

  for (idxParticle = 0; idxParticle < nCharges; ++idxParticle) 
    for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
      for (idx = 0; idx < nNames; ++idx) 
	for (centIdx = 0; centIdx < nCent; centIdx++) {
	  hSpectra[idxParticle][idx][energyIdx][centIdx] = static_cast<TH1D*>(spectraFile->Get(Form("spectracorrected.%s.%s.%s.%s", 
												    spectraParticles[idxParticle],
												    spectraNames[idx], spectraCent[centIdx], 
												    spectraEnergies[energyIdx])));
	}
  
  // ----------------------------------------------------------
  // -- fit spectra
  // ----------------------------------------------------------
  for (idxParticle = 0; idxParticle < nCharges; ++idxParticle) 
    for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
      for (idx = 0; idx < nNames; ++idx) 
	for (centIdx = 0; centIdx < nCent; centIdx++) {
	  TH1D* hist = hSpectra[idxParticle][idx][energyIdx][centIdx];
	  
	  funSpectra[idxParticle][idx][energyIdx][centIdx] = new TF1(Form("funSpectra%s_%s_%s_%s", chargedModes[idxParticle],
									  names[idx], energies[energyIdx], cent[centIdx]),  
								     "[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))", 0.10, 10.0);
	  
	  aSpectraFits.Add(funSpectra[idxParticle][idx][energyIdx][centIdx]);

	  TF1* fun = funSpectra[idxParticle][idx][energyIdx][centIdx];
	  
	  fun->SetParameters(100.*hist->GetBinContent(10), 1.55836, 1.10682);

	  fun->SetParLimits(1, parLimits1[idx][idxParticle][0], parLimits1[idx][idxParticle][1]);
	  fun->SetParLimits(2, parLimits2[idx][idxParticle][0], parLimits2[idx][idxParticle][1]);
	  
	  int fitFlag = hist->Fit(fun, "QRsame", "", 0.10, 10.0);
	  TFitResultPtr fitResult = hist->Fit(fun, "QRSsame", "", 0.10, 10.0);
	  //	  cout << " ADD ID   #chi^{2} = " << (double) fitResult.Chi2() << endl;
	  cout << "fit flag           = " << fitFlag << endl;
	
	  // ----------------------------------------------------------

	  cout << "ID  = " << fun->GetName() << endl;
	  cout << "[0] = " << fun->GetParameter(0) << endl;
	  cout << "[1] = " << fun->GetParameter(1) << endl;
	  cout << "[2] = " << fun->GetParameter(2) << endl;

	} // for (centIdx = 0; centIdx < nCent; centIdx++) {

  // ----------------------------------------------------------
  // -- Set range for fitted spectra
  // ----------------------------------------------------------
    for (idxParticle = 0; idxParticle < nCharges; idxParticle++) 
    for (energyIdx = 0 ; energyIdx < nEnergies; energyIdx++) 
      for (idx = 0; idx < nNames; ++idx) 
	for (centIdx = 0; centIdx < nCent; centIdx++) 
	  funSpectra[idxParticle][idx][energyIdx][centIdx]->SetRange(0.1, 10.);

  // ----------------------------------------------------------
  // -- Get charged spectra functions
  // ----------------------------------------------------------
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx) {
      TString namePostFix(Form("%s_%s", energies[energyIdx], cent[centIdx]));
      funSpectraCharged[0][energyIdx][centIdx] = new TF1(Form("funSpectraCharged_%s_%s", chargedModes[0], namePostFix.Data()), spectraChargedNeg, 0.10, 10.0);
      funSpectraCharged[1][energyIdx][centIdx] = new TF1(Form("funSpectraCharged_%s_%s", chargedModes[1], namePostFix.Data()), spectraChargedPos, 0.10, 10.0);
      funSpectraCharged[2][energyIdx][centIdx] = new TF1(Form("funSpectraCharged_%s_%s", chargedModes[2], namePostFix.Data()), spectraChargedAll, 0.10, 10.0);
    }

  // ----------------------------------------------------------
  // -- Get charged weighted efficiency functions
  // ----------------------------------------------------------
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx) {
      TString namePostFix(Form("%s_%s", energies[energyIdx], cent[centIdx]));
      funEffWeightedCharged[0][energyIdx][centIdx] = new TF1(Form("funEffWeightedCharged_%s_%s", chargedModes[0], namePostFix.Data()), effWeightedChargedNeg, 0.10, 10.0);
      funEffWeightedCharged[1][energyIdx][centIdx] = new TF1(Form("funEffWeightedCharged_%s_%s", chargedModes[1], namePostFix.Data()), effWeightedChargedPos, 0.10, 10.0);
      funEffWeightedCharged[2][energyIdx][centIdx] = new TF1(Form("funEffWeightedCharged_%s_%s", chargedModes[2], namePostFix.Data()), effWeightedChargedAll, 0.10, 10.0);
    }

  // ----------------------------------------------------------
  // -- Get charged efficiency functions (weighted efficiency / spectra)
  // ----------------------------------------------------------
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx) {
      TString namePostFix(Form("%s_%s", energies[energyIdx], cent[centIdx]));
      funEffCharged[0][energyIdx][centIdx] = new TF1(Form("funEffCharged_%s_%s", chargedModes[0], namePostFix.Data()), effChargedNeg, 0.10, 10.0);
      funEffCharged[1][energyIdx][centIdx] = new TF1(Form("funEffCharged_%s_%s", chargedModes[1], namePostFix.Data()), effChargedPos, 0.10, 10.0);
      funEffCharged[2][energyIdx][centIdx] = new TF1(Form("funEffCharged_%s_%s", chargedModes[2], namePostFix.Data()), effChargedAll, 0.10, 10.0);
    }

  // ----------------------------------------------------------
  // -- Get charged average efficiency
  // ----------------------------------------------------------
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx)
      for (idx = 0; idx < 3; ++idx) {
	avgEffCharged[idx][energyIdx][centIdx] = funEffWeightedCharged[idx][energyIdx][centIdx]->Integral(ptRange[0], ptRange[1]) / 
	  funSpectraCharged[idx][energyIdx][centIdx]->Integral(ptRange[0], ptRange[1]); 
      }

  // ----------------------------------------------------------
  // -- Get average efficiency vs refmult
  // ----------------------------------------------------------
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (idx = 0; idx < 3; ++idx) {

      gAvgEffChargedRefMult2[idx][energyIdx] = new TGraph(nCent, refMult2[energyIdx], avgEffCharged[idx][energyIdx]);
      gAvgEffChargedRefMult2[idx][energyIdx]->SetNameTitle(Form("gAvgEffChargedRefMult2_%s_%s", chargedModes[idx], energies[energyIdx]),
							   Form("avg EffCharged vs refMult2 %s %s", chargedModes[idx], exactEnergies[energyIdx]));
      aEffFits.Add(gAvgEffChargedRefMult2[idx][energyIdx]);

      funAvgEffChargedRefMult2[idx][energyIdx] = new TF1(Form("funAvgEffChargedRefMult2_%s_%s", chargedModes[idx], energies[energyIdx]), "pol1", 0., 500.0);
      gAvgEffChargedRefMult2[idx][energyIdx]->Fit(funAvgEffChargedRefMult2[idx][energyIdx], "R", "", 0, 500);
    }

  // ----------------------------------------------------------
  // -- Write average efficiency
  // ----------------------------------------------------------
  ofstream outAvgEff(Form("./results/chargedEff/averageChargedEff.txt"));

  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (idx = 0; idx < 3; ++idx) {          
      outAvgEff << energies[energyIdx] << " \t" << chargedModes[idx];    
      for (centIdx = 0; centIdx < nCent; ++centIdx)
	outAvgEff << " \t" << avgEffCharged[idx][energyIdx][centIdx];
      outAvgEff << endl;
    }
  
  outAvgEff.close();
  
  // ----------------------------------------------------------
  // -- Write all functions and histograms
  // ----------------------------------------------------------
  TFile *outFileChargeEff = TFile::Open("./results/fits/chargedEff.root", "RECREATE");

  outFileChargeEff->cd();
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx) 
      for (idx = 0; idx < 3; idx++) 
	funEffCharged[idx][energyIdx][centIdx]->Write();

  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx) 
      for (idx = 0; idx < 3; idx++) 
	funEffWeightedCharged[idx][energyIdx][centIdx]->Write();

  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (centIdx = 0; centIdx < nCent; ++centIdx) 
      for (idx = 0; idx < 3; idx++) 
	funSpectraCharged[idx][energyIdx][centIdx]->Write();

  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) 
    for (idx = 0; idx < 3; ++idx)
      funAvgEffChargedRefMult2[idx][energyIdx]->Write();

  outFileChargeEff->Close();

  // ----------------------------------------------------------
  TFile *outFileInputs = TFile::Open("./results/fits/chargedEffInput.root", "RECREATE");
  outFileInputs->cd();
  aEffFits.Write();
  aSpectraFits.Write();
  outFileInputs->Close();

  // ----------------------------------------------------------
  // -- Draw functions - fitted spectra
  // ----------------------------------------------------------
  drawSpectra();

  // ----------------------------------------------------------
  // -- Draw functions - charged efficiencies
  // ----------------------------------------------------------
  drawSpectraCharged();
  drawEffCharged();

  // ----------------------------------------------------------
  // -- Write canvas by energy 
  // ----------------------------------------------------------
  for (energyIdx = 0; energyIdx < nEnergies; ++energyIdx) {
    canEffCharged[energyIdx]->SaveAs(Form("./results/chargedEff/root/chargedEff_%sGeV.root",    energies[energyIdx]));
    canEffCharged[energyIdx]->SaveAs(Form("./results/chargedEff/root_macro/chargedEff_%sGeV.C", energies[energyIdx]));
    canEffCharged[energyIdx]->SaveAs(Form("./results/chargedEff/png/chargedEff_%sGeV.png",      energies[energyIdx]));
    canEffCharged[energyIdx]->SaveAs(Form("./results/chargedEff/pdf/chargedEff_%sGeV.pdf",      energies[energyIdx]));
  }

  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (idxParticle = 0; idxParticle < nCharges; idxParticle++) {
      canSpectra[idxParticle][energyIdx]->SaveAs(Form("./results/spectra/root/spectra_%s_%sGeV.root",    chargedModes[idxParticle],energies[energyIdx]));
      canSpectra[idxParticle][energyIdx]->SaveAs(Form("./results/spectra/root_macro/spectra_%s_%sGeV.C", chargedModes[idxParticle],energies[energyIdx]));
      canSpectra[idxParticle][energyIdx]->SaveAs(Form("./results/spectra/png/spectra_%s_%sGeV.png",      chargedModes[idxParticle],energies[energyIdx]));
      canSpectra[idxParticle][energyIdx]->SaveAs(Form("./results/spectra/pdf/spectra_%s_%sGeV.pdf",      chargedModes[idxParticle],energies[energyIdx]));
    }
  }
}

// =========================================================================================================
// =========================================================================================================
// =========================================================================================================

// _________________________________________________________________
void drawEffCharged() {
  // -- Draw charged efficiency

  TFile *inFileChargeEff = TFile::Open("./results/fits/chargedEff.root");
 
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    canEffCharged[energyIdx] = new TCanvas(Form("canChargedEff_%s", energies[energyIdx]), Form("%s GeV - charged efficiency", exactEnergies[energyIdx]),0,0,1200,600);
    canEffCharged[energyIdx]->SetFillColor(0);
    canEffCharged[energyIdx]->SetBorderMode(0);
    canEffCharged[energyIdx]->SetBorderSize(0.0);
    canEffCharged[energyIdx]->SetFrameFillColor(0);
    canEffCharged[energyIdx]->SetFrameBorderMode(0);

    canEffCharged[energyIdx]->cd();

    TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
    pad->SetBorderMode(0);
    pad->SetFillColor(1182);
    pad->Draw();
    pad->cd();
    pad->Divide(5,2,0.,0.,0.);
    
    TLegend *leg1= new TLegend(0.15,0.25,0.8,0.50);
    leg1->SetName(Form("leg1_%s", names[idx]));
    leg1->SetTextAlign(12);
    leg1->SetTextSize(0.10);
    leg1->SetTextFont(42);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(0);
    leg1->SetBorderSize(0);
    
    for (centIdx = 0; centIdx < nCent; centIdx++) {
      pad->cd(centIdx+1);
      
      if (centIdx == 4) {
	gPad->SetRightMargin(0.002);
	gPad->SetBottomMargin(0.004);
      }
      if (centIdx == 8) {
	gPad->SetRightMargin(0.003);
      }
      
      TH2D *ff = new TH2D("","",20, 0.009, 2.09, 20, 0.01, 0.99);
      ff->GetXaxis()->SetLabelSize(0.07);
      ff->GetYaxis()->SetLabelSize(0.07);
      ff->GetXaxis()->SetNdivisions(9, 5, 0);
      ff->Draw();
      
      TLine *line02 = new TLine(ptRange[0], 0, ptRange[0], 1);
      line02->SetLineColor(kGray+4);
      line02->SetLineStyle(3);
      line02->Draw();
      
      TLine *line20 = new TLine(ptRange[1], 0, ptRange[1], 1);
      line20->SetLineColor(kGray+4);
      line20->SetLineStyle(3);
      line20->Draw();

      TString namePostFix(Form("%s_%s", energies[energyIdx], cent[centIdx]));
      for (idx = 0; idx < 3; ++idx) {
	funEffChargedIn[idx][energyIdx][centIdx] = static_cast<TF1*>(inFileChargeEff->Get(Form("funEffCharged_%s_%s", chargedModes[idx], namePostFix.Data())));
	funEffChargedIn[idx][energyIdx][centIdx]->SetLineColor(chargedColor[idx]);
	funEffChargedIn[idx][energyIdx][centIdx]->Draw("same");
      }
     
      TLine *lineEff[3];
      for (int idx = 0; idx < 3; ++idx) {
	lineEff[idx] = new TLine(ptRange[0], avgEffCharged[idx][energyIdx][centIdx], ptRange[1], avgEffCharged[idx][energyIdx][centIdx]);
	lineEff[idx]->SetLineColor(chargedColor[idx]);
	lineEff[idx]->SetLineStyle(7);
	lineEff[idx]->Draw();
      }

      TLatex *texb_Cent = new TLatex(0.8,0.1,cent1[centIdx]);
      texb_Cent->SetTextSize(0.10);
      texb_Cent->SetTextFont(42);
      texb_Cent->Draw("same");
      
      pad->Modified();
      canEffCharged[energyIdx]->cd();
      
      if (centIdx == 0) {
	leg1->AddEntry(funEffChargedIn[0][energyIdx][centIdx], "neg. Hadrons","l");
	leg1->AddEntry((TObject*)0, "", "");
	leg1->AddEntry(funEffChargedIn[1][energyIdx][centIdx], "pos. Hadrons","l");
	leg1->AddEntry((TObject*)0, "", "");
	leg1->AddEntry(funEffChargedIn[2][energyIdx][centIdx], "avg. Hadrons","l");
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
    canEffCharged[energyIdx]->cd();
    
    TLatex *texb_5 = new TLatex(0.45,0.06,"#it{p}_{T} (GeV/#it{c})");
    texb_5->SetTextSize(0.03);
    texb_5->SetTextFont(42);
    texb_5->Draw("same");
    
    TLatex *texb_6 = new TLatex(0.03,0.4, Form("Efficiency - charged Hadrons"));
    texb_6->SetTextSize(0.03);
    texb_6->SetTextFont(42);
    texb_6->SetTextAngle(90);
    texb_6->Draw("same");
    
    pad->Modified();
    canEffCharged[energyIdx]->cd();
  } // end energy idx
}

// _________________________________________________________________
void drawSpectra() {
  // -- Draw charged spectra
  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    for (idxParticle = 0; idxParticle < nCharges; idxParticle++) {
      canSpectra[idxParticle][energyIdx] = new TCanvas(Form("canSpectra_%s_%s", chargedModes[idxParticle], energies[energyIdx]), 
					Form("%s GeV - %s - spectra", exactEnergies[energyIdx], chargedModes[idxParticle]),0,0,1200,600);
      canSpectra[idxParticle][energyIdx]->SetFillColor(0);
      canSpectra[idxParticle][energyIdx]->SetBorderMode(0);
      canSpectra[idxParticle][energyIdx]->SetBorderSize(0.0);
      canSpectra[idxParticle][energyIdx]->SetFrameFillColor(0);
      canSpectra[idxParticle][energyIdx]->SetFrameBorderMode(0);
      
      canSpectra[idxParticle][energyIdx]->cd();

      TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
      pad->SetBorderMode(0);
      pad->SetFillColor(1182);
      pad->Draw();
      pad->cd();
      pad->Divide(5,2,0.,0.,0.);
      
      TLegend *leg1= new TLegend(0.15,0.25,0.8,0.50);
      leg1->SetName(Form("leg1_%s", names[idx]));
      leg1->SetTextAlign(12);
      leg1->SetTextSize(0.10);
      leg1->SetTextFont(42);
      leg1->SetFillColor(kWhite);
      leg1->SetLineColor(0);
      leg1->SetBorderSize(0);
      
      for (int centIdx = 0; centIdx < nCent; centIdx++) {
	pad->cd(centIdx+1);
	gPad->SetLogy();

	if (centIdx == 4) {
	  gPad->SetRightMargin(0.002);
	  gPad->SetBottomMargin(0.004);
	}
	if (centIdx == 8) {
	  gPad->SetRightMargin(0.003);
	}
	
	TH2D *ff = new TH2D("","",20,0.009,2.09, 20, 10e-5, 5e2);
	ff->GetXaxis()->SetLabelSize(0.07);
	ff->GetYaxis()->SetLabelSize(0.07);
	ff->GetXaxis()->SetNdivisions(9,5,0);
	ff->Draw();
	
	TLine *line02 = new TLine(ptRange[0], 10e-5, ptRange[0], 5e2);
	line02->SetLineColor(kGray+4);
	line02->SetLineStyle(3);
	line02->Draw();
      
	TLine *line20 = new TLine(ptRange[1], 10e-5, ptRange[1], 5e2);
	line20->SetLineColor(kGray+4);
	line20->SetLineStyle(3);
	line20->Draw();

	funSpectra[idxParticle][0][energyIdx][centIdx]->SetLineColor(kAzure);
	hSpectra[idxParticle][0][energyIdx][centIdx]->SetLineColor(kAzure);
	hSpectra[idxParticle][0][energyIdx][centIdx]->SetMarkerColor(kAzure);
	hSpectra[idxParticle][0][energyIdx][centIdx]->SetMarkerStyle(24);
	hSpectra[idxParticle][0][energyIdx][centIdx]->Draw("psame");
	funSpectra[idxParticle][0][energyIdx][centIdx]->Draw("psame");

	funSpectra[idxParticle][1][energyIdx][centIdx]->SetLineColor(kRed+2);
	hSpectra[idxParticle][1][energyIdx][centIdx]->SetLineColor(kRed+2);
	hSpectra[idxParticle][1][energyIdx][centIdx]->SetMarkerColor(kRed+2);
	hSpectra[idxParticle][1][energyIdx][centIdx]->SetMarkerStyle(25);
	hSpectra[idxParticle][1][energyIdx][centIdx]->Draw("psame");
	funSpectra[idxParticle][1][energyIdx][centIdx]->Draw("psame");

	funSpectra[idxParticle][2][energyIdx][centIdx]->SetLineColor(kGreen+2);
	hSpectra[idxParticle][2][energyIdx][centIdx]->SetLineColor(kGreen+2);
	hSpectra[idxParticle][2][energyIdx][centIdx]->SetMarkerColor(kGreen+2);
	hSpectra[idxParticle][2][energyIdx][centIdx]->SetMarkerStyle(26);
	hSpectra[idxParticle][2][energyIdx][centIdx]->Draw("psame");
	funSpectra[idxParticle][2][energyIdx][centIdx]->Draw("psame");
	
	TLatex *texb_Cent = new TLatex(0.3,0.001,cent1[centIdx]);
	texb_Cent->SetTextSize(0.10);
	texb_Cent->SetTextFont(42);
	texb_Cent->Draw("same");
	
	pad->Modified();
	canSpectra[idxParticle][energyIdx]->cd();
	
	if (centIdx == 0) {
	  leg1->AddEntry(hSpectra[idxParticle][0][energyIdx][centIdx], namesTitle[idxParticle][0],"ple");
	  leg1->AddEntry((TObject*)0, "", "");
	  leg1->AddEntry(hSpectra[idxParticle][1][energyIdx][centIdx], namesTitle[idxParticle][1],"ple");
	  leg1->AddEntry((TObject*)0, "", "");
	  leg1->AddEntry(hSpectra[idxParticle][2][energyIdx][centIdx], namesTitle[idxParticle][2],"ple");
	}
      }
      
      pad->cd(10);
      
      TH2D *ff = new TH2D("","", 20,0.009,2.09, 20,0.01,0.99);
      ff->GetXaxis()->SetLabelSize(0.07);
      ff->GetYaxis()->SetLabelSize(0.07);
      ff->SetNdivisions(95, "X");
      
      leg1->Draw("lt");
      
      TLatex *texb_1 = new TLatex(0.15,0.82,"#it{p}_{T} Spectra");
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
      canSpectra[idxParticle][energyIdx]->cd();
      
      TLatex *texb_5 = new TLatex(0.45,0.06,"#it{p}_{T} (GeV/#it{c})");
      texb_5->SetTextSize(0.03);
      texb_5->SetTextFont(42);
      texb_5->Draw("same");
      
      TLatex *texb_6 = new TLatex(0.03,0.4, Form("Spectra"));
      texb_6->SetTextSize(0.03);
      texb_6->SetTextFont(42);
      texb_6->SetTextAngle(90);
      texb_6->Draw("same");
      
      pad->Modified();
      canSpectra[idxParticle][energyIdx]->cd();
    } // end names idx
  } // end energy idx

}


// _________________________________________________________________
void drawSpectraCharged() {
  // -- Draw charged spectra

  TFile *inFileChargeEff = TFile::Open("./results/fits/chargedEff.root");

  for (energyIdx = 0 ; energyIdx < nEnergies; ++energyIdx) {
    canSpectraCharged[energyIdx] = new TCanvas(Form("canChargedSpectra_%s", energies[energyIdx]), Form("%s GeV - charged spectra", exactEnergies[energyIdx]),0,0,1200,600);
    canSpectraCharged[energyIdx]->SetFillColor(0);
    canSpectraCharged[energyIdx]->SetBorderMode(0);
    canSpectraCharged[energyIdx]->SetBorderSize(0.0);
    canSpectraCharged[energyIdx]->SetFrameFillColor(0);
    canSpectraCharged[energyIdx]->SetFrameBorderMode(0);

    canSpectraCharged[energyIdx]->cd();

    TPad* pad = new TPad("pad", "pad",0.05,0.1,0.99,0.99);
    pad->SetBorderMode(0);
    pad->SetFillColor(1182);
    pad->Draw();
    pad->cd();
    pad->Divide(5,2,0.,0.,0.);
    
    TLegend *leg1= new TLegend(0.15,0.25,0.8,0.50);
    leg1->SetName(Form("leg1_%s", names[idx]));
    leg1->SetTextAlign(12);
    leg1->SetTextSize(0.10);
    leg1->SetTextFont(42);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(0);
    leg1->SetBorderSize(0);
    
    for (centIdx = 0; centIdx < nCent; centIdx++) {
      pad->cd(centIdx+1);
      gPad->SetLogy();

      if (centIdx == 4) {
	gPad->SetRightMargin(0.002);
	gPad->SetBottomMargin(0.004);
      }
      if (centIdx == 8) {
	gPad->SetRightMargin(0.003);
      }
      
      TH2D *ff = new TH2D("","",20, 0.009, 2.09, 20, 10e-5, 5e2);
      ff->GetXaxis()->SetLabelSize(0.07);
      ff->GetYaxis()->SetLabelSize(0.07);
      ff->GetXaxis()->SetNdivisions(9, 5, 0);
      ff->Draw();
      
      TLine *line02 = new TLine(ptRange[0], 10e-5, ptRange[0], 5e2);
      line02->SetLineColor(kGray+4);
      line02->SetLineStyle(3);
      line02->Draw();
      
      TLine *line20 = new TLine(ptRange[1], 10e-5, ptRange[1], 5e2);
      line20->SetLineColor(kGray+4);
      line20->SetLineStyle(3);
      line20->Draw();

      TString namePostFix(Form("%s_%s", energies[energyIdx], cent[centIdx]));
      for (idx = 0; idx < 3; ++idx) {
	funSpectraChargedIn[idx][energyIdx][centIdx] = static_cast<TF1*>(inFileChargeEff->Get(Form("funSpectraCharged_%s_%s", chargedModes[idx], namePostFix.Data())));
	funSpectraChargedIn[idx][energyIdx][centIdx]->SetLineColor(chargedColor[idx]);
	funSpectraChargedIn[idx][energyIdx][centIdx]->Draw("same");
      }

      TLatex *texb_Cent = new TLatex(0.3,0.001,cent1[centIdx]);
      texb_Cent->SetTextSize(0.10);
      texb_Cent->SetTextFont(42);
      texb_Cent->Draw("same");
      
      pad->Modified();
      canSpectraCharged[energyIdx]->cd();
      
      if (centIdx == 0) {
	leg1->AddEntry(funSpectraChargedIn[0][energyIdx][centIdx], "neg. Hadrons","l");
	leg1->AddEntry((TObject*)0, "", "");
	leg1->AddEntry(funSpectraChargedIn[1][energyIdx][centIdx], "pos. Hadrons","l");
	leg1->AddEntry((TObject*)0, "", "");
	leg1->AddEntry(funSpectraChargedIn[2][energyIdx][centIdx], "avg. Hadrons","l");
      }
    }
    
    pad->cd(10);
     
    TH2D *ff = new TH2D("","", 20,0.009,2.09, 20,0.01,0.99);
    ff->GetXaxis()->SetLabelSize(0.07);
    ff->GetYaxis()->SetLabelSize(0.07);
    ff->SetNdivisions(95, "X");
    
    leg1->Draw("lt");
    
    TLatex *texb_1 = new TLatex(0.15,0.82,"#it{p}_{T} Spectra");
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
    canSpectraCharged[energyIdx]->cd();
    
    TLatex *texb_5 = new TLatex(0.45,0.06,"#it{p}_{T} (GeV/#it{c})");
    texb_5->SetTextSize(0.03);
    texb_5->SetTextFont(42);
    texb_5->Draw("same");
    
    TLatex *texb_6 = new TLatex(0.03,0.4, Form("Spectra"));
    texb_6->SetTextSize(0.03);
    texb_6->SetTextFont(42);
    texb_6->SetTextAngle(90);
    texb_6->Draw("same");
    
    pad->Modified();
    canSpectraCharged[energyIdx]->cd();
  } // end energy idx
}

// =========================================================================================================
// =========================================================================================================
// =========================================================================================================

// __________________________________________________________________
double spectraChargedNeg(double *x, double *par) {
  // -- get spectra for negative particles
  return ( funSpectra[0][0][energyIdx][centIdx]->EvalPar(x,par) + 
	   funSpectra[0][1][energyIdx][centIdx]->EvalPar(x,par) + 
	   funSpectra[0][2][energyIdx][centIdx]->EvalPar(x,par) );
}

// __________________________________________________________________
double spectraChargedPos(double *x, double *par) {
  // -- get spectra for positive particles
  return ( funSpectra[1][0][energyIdx][centIdx]->EvalPar(x,par) + 
	   funSpectra[1][1][energyIdx][centIdx]->EvalPar(x,par) + 
	   funSpectra[1][2][energyIdx][centIdx]->EvalPar(x,par) );
}

// __________________________________________________________________
double spectraChargedAll(double *x, double *par) {
  // -- get spectra for positive particles
  return ( funSpectraCharged[0][energyIdx][centIdx]->EvalPar(x,par) + 
	   funSpectraCharged[1][energyIdx][centIdx]->EvalPar(x,par) );
}

// __________________________________________________________________
double effWeightedChargedNeg(double *x, double *par) {
  // -- get weighted efficiency for negative particles
  return  ( (funSpectra[0][0][energyIdx][centIdx]->EvalPar(x,par)*funEff[0][0][energyIdx][centIdx]->EvalPar(x, par)) +
	    (funSpectra[0][1][energyIdx][centIdx]->EvalPar(x,par)*funEff[0][1][energyIdx][centIdx]->EvalPar(x, par)) +
	    (funSpectra[0][2][energyIdx][centIdx]->EvalPar(x,par)*funEff[0][2][energyIdx][centIdx]->EvalPar(x, par)) );
}

// __________________________________________________________________
double effWeightedChargedPos(double *x, double *par) {
  // -- get weighted efficiency for postive particles
  return  ( (funSpectra[1][0][energyIdx][centIdx]->EvalPar(x,par)*funEff[1][0][energyIdx][centIdx]->EvalPar(x, par)) +
	    (funSpectra[1][1][energyIdx][centIdx]->EvalPar(x,par)*funEff[1][1][energyIdx][centIdx]->EvalPar(x, par)) +
	    (funSpectra[1][2][energyIdx][centIdx]->EvalPar(x,par)*funEff[1][2][energyIdx][centIdx]->EvalPar(x, par)) );
}

// __________________________________________________________________
double effWeightedChargedAll(double *x, double *par) {
  // -- get weighted efficiency for postive+negative particles
  return ( funEffWeightedCharged[0][energyIdx][centIdx]->EvalPar(x,par) +
	   funEffWeightedCharged[1][energyIdx][centIdx]->EvalPar(x,par) );
}

// __________________________________________________________________
double effChargedNeg(double *x, double *par) {
  // -- get weighted efficiency for negative particles
  return ( funEffWeightedCharged[0][energyIdx][centIdx]->EvalPar(x,par) /
	   funSpectraCharged[0][energyIdx][centIdx]->EvalPar(x,par)  );
}

// __________________________________________________________________
double effChargedPos(double *x, double *par) {
  // -- get weighted efficiency for postive particles
 return ( funEffWeightedCharged[1][energyIdx][centIdx]->EvalPar(x,par) /
	  funSpectraCharged[1][energyIdx][centIdx]->EvalPar(x,par)  );
}

// __________________________________________________________________
double effChargedAll(double *x, double *par) {
  // -- get weighted efficiency for postive+negative particles
 return ( funEffWeightedCharged[2][energyIdx][centIdx]->EvalPar(x,par) /
	  funSpectraCharged[2][energyIdx][centIdx]->EvalPar(x,par)  );
}


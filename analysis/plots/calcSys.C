
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

      float   aEta[]           = {0.5};
const int     nEta             = 1;

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S#sigma", "#kappa#sigma^{2}", "(S#sigma - M/#sigma^{2})"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV", "VM"};
const int     nMoments         = 7;

const int     nSysNames = 13;
const Char_t* aSysNames[] = { "base", 
			      "sys_1_0", "sys_1_1", "sys_1_2", "sys_1_3",
			      "sys_0_0", "sys_0_1", "sys_0_2", "sys_0_3",
			      "base_plus5", "base_minus5", 
			      "base_plus2minus2","base_minus2plus2" };

const Char_t* aDataSets[]      = {"twoeff_11"};
const int     nDataSets        = 1;

const Char_t* aDataSetsTitle[] = {"corrected - #epsilon_{1} , #epsilon_{2}"};

const Char_t* aSysTitle[] = {"default", "eff (+5% && +5%)", "eff (-5% && -5%)", "eff (+2% && -2%)", "eff (-2% && +2%)"};

Bool_t useRatioSDSkellam = kTRUE;
//Bool_t useRatioSDSkellam = kFALSE;

TObjArray canA;
TObjArray sysErrorA;
TObjArray statErrorA;

// ______________________________________________________________________________________
void calcSys(const Char_t* name = "jobs_14.5") {
  // -- get systematic errors

  TFile *inFiles[nSysNames][nDataSets][nMoments];
  TGraphErrors *graphs[nSysNames][nDataSets][nMoments];

  for (int idxSys = 0; idxSys < nSysNames; ++idxSys) 
    for (int idxDataSet = 0; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0; idxMoment < nMoments; ++idxMoment) {
	inFiles[idxSys][idxDataSet][idxMoment] = TFile::Open(Form("output/%s_%s/%s/Moments_%s.root", name, aSysNames[idxSys],
								  aDataSets[idxDataSet], aMoments[idxMoment]));

	 graphs[idxSys][idxDataSet][idxMoment] = (useRatioSDSkellam && idxMoment == 5) ?
	   static_cast<TGraphErrors*>(inFiles[idxSys][idxDataSet][idxMoment]->Get(Form("%s_Poisson_ratio",aMoments[idxMoment]))->Clone()) :	   
	   static_cast<TGraphErrors*>(inFiles[idxSys][idxDataSet][idxMoment]->Get(aMoments[idxMoment])->Clone()) ;
	   
	 if (useRatioSDSkellam && idxMoment == 5) {
	   if (idxSys == 0)
	     graphs[idxSys][idxDataSet][idxMoment]->SetName(Form("%s_Poisson_ratio_stat", aMoments[idxMoment]));
	   else
	     graphs[idxSys][idxDataSet][idxMoment]->SetName(Form("%s_Poisson_ratio_sys_%d", aMoments[idxMoment], idxSys));
	 }
	 else {
	   if (idxSys == 0)
	     graphs[idxSys][idxDataSet][idxMoment]->SetName(Form("%s_stat", aMoments[idxMoment]));
	   else
	     graphs[idxSys][idxDataSet][idxMoment]->SetName(Form("%s_sys_%d", aMoments[idxMoment], idxSys));
	 }

	 
	 inFiles[idxSys][idxDataSet][idxMoment]->Close();
      }

  // -----------------------------------------------------

  TGraphErrors *graphsSysError[nDataSets][nMoments];

  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {

      statErrorA.Add(graphs[0][idxDataSet][idxMoment]);

      Double_t* defaultX  = graphs[0][idxDataSet][idxMoment]->GetX(); 
      Double_t* defaultY  = graphs[0][idxDataSet][idxMoment]->GetY(); 
      Double_t* defaultEY = graphs[0][idxDataSet][idxMoment]->GetEY(); 
            
      // "sys_1_0", "sys_1_1", "sys_1_2", "sys_1_3",
      // "sys_0_0", "sys_0_1", "sys_0_2", "sys_0_3",
      // "base_plus5", "base_minus5", 
      // "base_plus2minus2","base_minus2plus2" };
      
      Double_t values[4][4][9];
      for (Int_t hh = 0; hh < 4; ++hh) 
	for (Int_t ii = 0; ii < 4; ++ii) 
	  for (Int_t jj = 0; jj < 9; ++jj) 
	    values[hh][ii][jj] = 0;

      Double_t *a;

      for (Int_t ii = 1; ii < 5; ++ii) {
	a = graphs[ii][idxDataSet][idxMoment]->GetY(); 
	for (Int_t bin = 0; bin < 9; bin++) 
	  values[0][ii-1][bin] = a[bin];
      }
      for (Int_t ii = 5; ii < 9; ++ii) {
	a = graphs[ii][idxDataSet][idxMoment]->GetY(); 
	for (Int_t bin = 0; bin < 9; bin++) 
	  values[1][ii-5][bin] = a[bin];
      }
      for (Int_t ii = 9; ii < 11; ++ii) {
	a = graphs[ii][idxDataSet][idxMoment]->GetY(); 
	for (Int_t bin = 0; bin < 9; bin++) 
	  values[2][ii-9][bin] = a[bin]; 
      }
      for (Int_t ii = 11; ii < 13; ++ii) {
	a = graphs[ii][idxDataSet][idxMoment]->GetY(); 
	for (Int_t bin = 0; bin < 9; bin++) 
	  values[3][ii-11][bin] = a[bin];
      }

      Double_t sum[4][9];
      Double_t rms[4][9];
      for (Int_t ii = 0; ii < 4; ++ii) 
	for (Int_t jj = 0; jj < 9; ++jj) {
	    rms[ii][jj] = 0;
	    sum[ii][jj] = 0;
	}

      // --------------------------------------------
      Int_t idxSet = 0;
      Double_t n   = 4;
      for (Int_t bin = 0; bin < 9; bin++) {
	for (Int_t ii = 1; ii < 5; ++ii) {
	  Int_t idx = ii-1;
	  sum[idxSet][bin] += (values[idxSet][idx][bin]-defaultY[bin])*(values[idxSet][idx][bin]-defaultY[bin]);
	}
	rms[idxSet][bin] = TMath::Sqrt(sum[idxSet][bin]/n);
      }

      // --------------------------------------------
      ++idxSet;
      n = 4;
      for (Int_t bin = 0; bin < 9; bin++) {    
	for (Int_t ii = 5; ii < 9; ++ii) {
	  Int_t idx = ii-5;
    	  sum[idxSet][bin] += (values[idxSet][idx][bin]-defaultY[bin])*(values[idxSet][idx][bin]-defaultY[bin]);
	}
	rms[idxSet][bin] = TMath::Sqrt(sum[idxSet][bin]/n);
      }

      // --------------------------------------------
      ++idxSet;
      n = 2;
      for (Int_t bin = 0; bin < 9; bin++) {
	for (Int_t ii = 9; ii < 11; ++ii) {
	  Int_t idx = ii-9;
	  sum[idxSet][bin] += (values[idxSet][idx][bin]-defaultY[bin])*(values[idxSet][idx][bin]-defaultY[bin]);
	}
	rms[idxSet][bin] = TMath::Sqrt(sum[idxSet][bin]/n);
      }

      // --------------------------------------------

      ++idxSet;
      n = 2;
      for (Int_t bin = 0; bin < 9; bin++) {
	for (Int_t ii = 11; ii < 13; ++ii) {
	  Int_t idx = ii-11;
	  sum[idxSet][bin] += (values[idxSet][idx][bin]-defaultY[bin])*(values[idxSet][idx][bin]-defaultY[bin]);
	}
	rms[idxSet][bin] = TMath::Sqrt(sum[idxSet][bin]/n);
      }

      // --------------------------------------------
      
      Double_t sysSum[9];
      Double_t sysErr[9];
      for (Int_t jj = 0; jj < 9; ++jj) {
	sysSum[jj] = 0;
	sysErr[jj] = 0;
      }

      for (Int_t bin = 0; bin <9; bin++) {
	for (Int_t idx = 0 ; idx< 4; ++idx)
	  sysSum[bin] += rms[idx][bin]*rms[idx][bin];
	sysErr[bin] = TMath::Sqrt(sysSum[bin]);
      }

      // --------------------------------------------
      
      cout << "MOM: " << idxMoment << endl;
      for (Int_t idx = 0 ; idx < 4; ++idx) {
	cout << "    Set " << idx << " : " << rms[idx][0] << endl;
      }
      cout << "    SYS " << " : " << sysErr[0] << endl;
      
      cout << "--------------------------------------------" << endl;
      

      graphsSysError[idxDataSet][idxMoment] = new TGraphErrors(9, defaultX, defaultY, 0, sysErr);

      if (useRatioSDSkellam && idxMoment == 5)
	graphsSysError[idxDataSet][idxMoment]->SetName(Form("%s_Poisson_ratio_sys", aMoments[idxMoment]));
      else 
	graphsSysError[idxDataSet][idxMoment]->SetName(Form("%s_sys", aMoments[idxMoment]));

      sysErrorA.Add(graphsSysError[idxDataSet][idxMoment]);
    }
  }

  // -----------------------------------------------------

  gSystem->Exec("mkdir -p sysError");

  if (nDataSets > 1)
    cout << "N DataSets larger 1 ... check sys error name" << endl;

  TFile *sysErrFile = TFile::Open(Form("sysError/%s_sysError.root", name), "RECREATE");
  sysErrFile->cd();
  sysErrorA.Write();
  statErrorA.Write();
  sysErrFile->Close();
}

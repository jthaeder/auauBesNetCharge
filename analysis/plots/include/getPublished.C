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


// ---------------------------------------------------------------------------------
//                    [snn/mub][moments][cent]
TGraphErrors *graphStat[2][8][9];     // 8 = C1-C4 + VM + SDSk + KV + SD
TGraphErrors *graphSys[2][8][9];      // 8 = C1-C4 + VM + SDSk + KV + SD
TGraphErrors *graphPoisson[2][8][9];  // 7 = C1-C4 + VM + SDSk + KV 
TGraphErrors *graphUrqmd[2][8][9];    // 7 = C1-C4 + VM + SDSk + KV 
TGraphErrors *graph14[2][8][9];       // 8 = C1-C4 + VM + SDSk + KV + SD

// ---------------------------------------------------------------------------------

void getPublished() {

  // -----------------------------------------------------------------------
  // -- Fig.4
  // --  (Beam energy dependance of Variance/Mean, Skwness*Std. Deviation and Kurtosis*Variance)
  // -----------------------------------------------------------------------
  
  /* -----------------------------------------------------------------------
   * Variance/Mean (Sigma^2/M)                                               
   * 
   * 0-5% Centrality       
   * 70-80% Centrality 
   *
   *                                  ENG (GeV)  sigma^2/M  Stat.Error  Sys.Error  NBD  Skellam
   * ----------------------------------------------------------------------- */ 
  Double_t publishedVM[2][7][6] = { { {7.7,  3.72455, 0.0319487, 0.0159744, 6.72266, 5.4224},
				      {11.5, 5.55127, 0.0366492, 0.0183246, 9.84451, 8.25644},
				      {19.6, 8.95204, 0.0255873, 0.0127937, 15.7966, 13.3401},
				      {27,   12.7145, 0.0299337, 0.0149669, 21.1496, 18.7789},
				      {39,   17.4625, 0.0263059, 0.013153,  30.1741, 25.4531},
				      {62.4, 27.3241, 0.0651797, 0.0325899, 49.3896, 39.0525},
				      {200,  80.2318, 0.210213,  0.105106,  124.898, 108.848} },
				    
				    { {7.7,  7.47332, 0.102868,  0.051434,  14.6179, 9.20048},
				      {11.5, 11.8146, 0.174244,  0.0871218, 22.8091, 15.4699},
				      {19.6, 17.9191, 0.124656,  0.062328,  35.0745, 21.3709},
				      {27,   26.9981, 0.217466,  0.108733,  54.9325, 37.9424},
				      {39,   32.6472, 0.204098,  0.102049,  66.1672, 44.6294},
				      {62.4, 47.2453, 0.62329,   0.311645,  94.9524, 64.815},
				      {200,  99.0916, 1.56255,   0.781274,  201.971, 139.686} } };
  
  /* -----------------------------------------------------------------------  
   * Skewness*Std. Deviation (S*sigma)                                                                         
   * 
   * 0-5% Centrality       
   * 70-80% Centrality 
   *
   *                                  ENG (GeV)  S*sigma  Stat.Error  Sys.Error  NBD  Skellam
   * ----------------------------------------------------------------------- */ 
  Double_t publishedSD[2][7][6] = { { {7.7,  0.443332,  0.237605,   0.215502,   0.476933,  0.184454},
				      {11.5, 0.503508,  0.208217,   0.174109,   0.290953,  0.121136},
				      {19.6, 0.283159,  0.09437,    0.0777163,  0.168042,  0.0749629},
				      {27,   0.107052,  0.0750091,  0.0562664,  0.111871,  0.053256},
				      {39,   0.0133622, 0.0501263,  0.033736,   0.0768153, 0.0392905},
				      {62.4, 0.145121,  0.067052,   0.0899881,  0.069852,  0.0255998},
				      {200, -0.0188162, 0.0477387,  0.0409418,  0.0138754, 0.00919266} },

				    { {7.7,  0.472407,  0.0347691,  0.0137378,  0.425598,  0.109716},
				      {11.5, 0.321615,  0.0291255,  0.0610068,  0.299105,  0.0648058},
				      {19.6, 0.204636,  0.0113238,  0.00967549, 0.189537,  0.0468747},
				      {27,   0.124164,  0.00965298, 0.00534643, 0.136668,  0.0265044},
				      {39,   0.103698,  0.00674335, 0.0388556,  0.103799,  0.0224682},
				      {62.4, 0.0711058, 0.0118494,  0.0197959,  0.0674433, 0.0154511},
				      {200,  0.0392179, 0.00761239, 0.021385,   0.0307338, 0.00717921} } };

  /* -----------------------------------------------------------------------  
   * Kurtosis*Variance (KV)
   * 
   * 0-5% Centrality       
   * 70-80% Centrality 
   *
   *                                  ENG (GeV)  kappa*sigma^2  Stat.Error  Sys.Error  NBD  Skellam
   * ----------------------------------------------------------------------- */ 
  Double_t publishedKV[2][7][6] = { { {7.7, -10.6236, 7.47779,   5.2201,    3.02669, 1},
				      {11.5, 1.90761, 7.09204,   1.30255,   2.51506, 1},
				      {19.6, 1.30334, 3.54197,   0.785453,  2.37491, 1},
				      {27,   2.96161, 3.09955,   1.51754,   1.90713, 1},
				      {39,   1.59113, 2.22224,   0.62209,   2.35608, 1},
				      {62.4, 4.57877, 3.34581,   2.46885,   3.04649, 1},
				      {200,  2.57939, 2.42202,   1.84051,   2.02919, 1} },

				    { {7.7,  3.00154, 0.263386,  0.147408,  5.63941, 1},
				      {11.5, 2.7137,  0.193646,  0.337091,  4.79646, 1},
				      {19.6, 2.43159, 0.0922055, 0.0991678, 4.6705,  1},
				      {27,   2.70718, 0.0839845, 0.124916,  5.00913, 1},
				      {39,   2.64718, 0.0537805, 0.142572,  5.04093, 1},
				      {62.4, 2.71035, 0.0964119, 0.316314,  5.06482, 1},
				      {200,  2.33482, 0.0692845, 0.244065,  4.90635, 1} } };


  // ---------------------------------------------------------------------------------
  // J. Cleymans et al., PRC 73, 034905 (2006).
  //    mub(sqrt(s)) = 1.308 /￼(1+0.27× sqrt(s))
  //    Tc(sqrt(s))  = 0.166−0.139 × mub^2 −0.053 mub^4
  // !  Unit : √s, μB, T (GeV)
  Double_t snn[] = {7.7, 11.5, 14.5, 19.6,  27,  39, 62.4, 200};
  Double_t mub[] = {422,  316,  266,  206, 156, 112,   73,  24};
  Double_t tc[]  = {140,  152,  156,  160, 163, 164,  165, 166};

  Double_t snnRed[] = {7.7, 11.5, 19.6,  27,  39, 62.4, 200};
  Double_t mubRed[] = {422,  316,  206, 156, 112,   73,  24};

  Double_t snn14[] = {14.5};
  Double_t mub14[] = {266};

  // ---------------------------------------------------------------------------------
  
  const Int_t nEnergyMax = 8;
  
  // Double_t exSys[nEnergyMax] = {0.75, 3, 3, 3, 4, 4, 4, 10};

  // for (Int_t idxEnergy = 0; idxEnergy < nEnergyMax; ++idxEnergy)   
  //   exSys[idxEnergy] = 1.;

  // ---------------------------------------------------------------------------------
  
  for (Int_t idxMoment = 4; idxMoment < 8; ++idxMoment) {

    Double_t y[2][nEnergyMax], eyStat[2][nEnergyMax], eySys[2][nEnergyMax], yPoisson[2][nEnergyMax];

    for (Int_t ii = 0; ii < 2; ++ii) {
      Int_t idxCent = (ii == 1) ? 8 : 0;
      
      for (Int_t idxEnergy = 0; idxEnergy < nEnergyMax; ++idxEnergy) {
	Int_t energyArrayAccess = (idxEnergy >= 2) ? idxEnergy-1 : idxEnergy;
	
	if (idxEnergy != 2) {
	  //  -- VM
	  if (idxMoment == 4) {  
	    y[ii][idxEnergy]        = publishedVM[ii][energyArrayAccess][1];
	    eyStat[ii][idxEnergy]   = publishedVM[ii][energyArrayAccess][2];
	    eySys[ii][idxEnergy]    = publishedVM[ii][energyArrayAccess][3];
	    yPoisson[ii][idxEnergy] = publishedVM[ii][energyArrayAccess][5];
	  }
	  //  -- SDSk
	  else if (idxMoment == 5) {  
	    y[ii][idxEnergy]        = publishedSD[ii][energyArrayAccess][1]/publishedSD[ii][energyArrayAccess][5];
	    eyStat[ii][idxEnergy]   = publishedSD[ii][energyArrayAccess][2]/TMath::Abs(publishedSD[ii][energyArrayAccess][5]);
	    eySys[ii][idxEnergy]    = publishedSD[ii][energyArrayAccess][3]/TMath::Abs(publishedSD[ii][energyArrayAccess][5]);
	    yPoisson[ii][idxEnergy] = 1.;
	  }
	  //  -- KV
	  else if (idxMoment == 6) {  
	    y[ii][idxEnergy]        = publishedKV[ii][energyArrayAccess][1];
	    eyStat[ii][idxEnergy]   = publishedKV[ii][energyArrayAccess][2];
	    eySys[ii][idxEnergy]    = publishedKV[ii][energyArrayAccess][3];
	    yPoisson[ii][idxEnergy] = publishedKV[ii][energyArrayAccess][5];
	  }
	  //  -- SD
	  else if (idxMoment == 7) {  
	    y[ii][idxEnergy]        = publishedSD[ii][energyArrayAccess][1];
	    eyStat[ii][idxEnergy]   = publishedSD[ii][energyArrayAccess][2];
	    eySys[ii][idxEnergy]    = publishedSD[ii][energyArrayAccess][3];
	    yPoisson[ii][idxEnergy] = publishedSD[ii][energyArrayAccess][5];
	  }
	}
	else {
	  y[ii][idxEnergy]        = 1.;
	  eyStat[ii][idxEnergy]   = 1.;
	  eySys[ii][idxEnergy]    = 1.;
	  yPoisson[ii][idxEnergy] = 1.;
	}
      } // for (Int_t idxEnergy = 0; idxEnergy < nEnergyMax; ++idxEnergy) {
      
      graphStat[0][idxMoment][idxCent]    = new TGraphErrors(nEnergyMax, snn, y[ii], 0, eyStat[ii]);
      graphStat[1][idxMoment][idxCent]    = new TGraphErrors(nEnergyMax, mub, y[ii], 0, eyStat[ii]);    
      
      // graphSys[0][idxMoment][idxCent]     = new TGraphErrors(nEnergyMax, snn, y[ii], exSys, eySys[ii]);    
      // graphSys[1][idxMoment][idxCent]     = new TGraphErrors(nEnergyMax, mub, y[ii], exSys, eySys[ii]);
      
      graphSys[0][idxMoment][idxCent]     = new TGraphErrors(nEnergyMax, snn, y[ii], 0, eySys[ii]);    
      graphSys[1][idxMoment][idxCent]     = new TGraphErrors(nEnergyMax, mub, y[ii], 0, eySys[ii]);
      
      graphPoisson[0][idxMoment][idxCent] = new TGraphErrors(nEnergyMax, snn, yPoisson[ii], 0, 0);    
      graphPoisson[1][idxMoment][idxCent] = new TGraphErrors(nEnergyMax, mub, yPoisson[ii], 0, 0);
    } // for (Int_t ii = 0; ii < 2; ++i) {
  } // for (Int_t idxMoment = 4; idxMoment < 8; ++idxMoment) {

  // -----------------------------------------------------------------------

  for (Int_t idxMoment = 0; idxMoment < 8; ++idxMoment) {
    for (Int_t idxCent = 0; idxCent < 9; ++idxCent) {
      if (!graphStat[0][idxMoment][idxCent]) {
	graphStat[0][idxMoment][idxCent]    = new TGraphErrors(nEnergyMax, snn, 0, 0, 0);    
	graphStat[1][idxMoment][idxCent]    = new TGraphErrors(nEnergyMax, mub, 0, 0, 0);    
      }
      graphStat[0][idxMoment][idxCent]->SetName(Form("stat_snn_%d_%d", idxMoment, idxCent));
      graphStat[1][idxMoment][idxCent]->SetName(Form("stat_mub_%d_%d", idxMoment, idxCent));

      if (!graphSys[0][idxMoment][idxCent]) {
	graphSys[0][idxMoment][idxCent]     = new TGraphErrors(nEnergyMax, snn, 0, 0, 0);    
	graphSys[1][idxMoment][idxCent]     = new TGraphErrors(nEnergyMax, mub, 0, 0, 0);    
      }
      graphSys[0][idxMoment][idxCent]->SetName(Form("sys_snn_%d_%d", idxMoment, idxCent));
      graphSys[1][idxMoment][idxCent]->SetName(Form("sys_mub_%d_%d", idxMoment, idxCent));

      if (!graphPoisson[0][idxMoment][idxCent]) {
	graphPoisson[0][idxMoment][idxCent] = new TGraphErrors(nEnergyMax, snn, 0, 0, 0);    
	graphPoisson[1][idxMoment][idxCent] = new TGraphErrors(nEnergyMax, mub, 0, 0, 0);    
      }
      graphPoisson[0][idxMoment][idxCent]->SetName(Form("poisson_snn_%d_%d", idxMoment, idxCent));
      graphPoisson[1][idxMoment][idxCent]->SetName(Form("poisson_mub_%d_%d", idxMoment, idxCent));

      if (!graphUrqmd[0][idxMoment][idxCent]) {
	graphUrqmd[0][idxMoment][idxCent]   = new TGraphErrors(nEnergyMax-1, snnRed, 0, 0, 0);    
	graphUrqmd[1][idxMoment][idxCent]   = new TGraphErrors(nEnergyMax-1, mubRed, 0, 0, 0);    
      }
      graphUrqmd[0][idxMoment][idxCent]->SetName(Form("urqmd_snn_%d_%d", idxMoment, idxCent));
      graphUrqmd[1][idxMoment][idxCent]->SetName(Form("urqmd_mub_%d_%d", idxMoment, idxCent));

      if (!graph14[0][idxMoment][idxCent]) {
	graph14[0][idxMoment][idxCent]   = new TGraphErrors(1, snn14, 0, 0, 0);    
	graph14[1][idxMoment][idxCent]   = new TGraphErrors(1, mub14, 0, 0, 0);    
      }
      graph14[0][idxMoment][idxCent]->SetName(Form("14_snn_%d_%d", idxMoment, idxCent));
      graph14[1][idxMoment][idxCent]->SetName(Form("14_mub_%d_%d", idxMoment, idxCent));
    }
  }

  // -----------------------------------------------------------------------

}


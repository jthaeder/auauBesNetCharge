#include "TFile.h"
#include "TList.h"
#include "TGraphErrors.h"

void print_STAR_QM2015_Preliminary() {

  TFile *file = TFile::Open("STAR_QM2015_Preliminary.root");

  const Char_t *aNames[3]  = {"Net-Charge", "Net-Kaon", "Net-Proton"};
  const Char_t *aMom[3]    = {"VM", "SDSk", "KV"};
  const Char_t *aMom2[3]   = {"VM", "SD", "KV"};

  const Char_t *aTitle[3]  = {"#sigma^2 / M", "S #sigma / Skellam", "#kappa #sigma^2"};
  const Char_t *aTitle2[3] = {"#sigma^2 / M", "S #sigma", "#kappa #sigma^2"};

  const Char_t* aCent[9]   = {"0005",  "0510",   "1020",   "2030",   "3040",   "4050",   "5060",   "6070",   "7080"};

  const Double_t aEta[5]   = {0.1, 0.3, 0.5, 0.7, 1.0};

  printf("============================================================\n");
  printf("-- STAR QM15 Preliminary Data\n");
  printf("============================================================\n");
  printf("---------------------------------------------------\n");
  printf(" -- Kinematic region\n");
  printf("---------------------------------------------------\n");
  printf("   Net-Charge: 0.2 < p_T (GeV/c) < 2.0 , |#eta| < 0.5\n");
  printf("   Net-Kaon:   0.2 < p_T (GeV/c) < 1.6 , |y| < 0.5\n");
  printf("   Net-Proton: 0.4 < p_T (GeV/c) < 2.0 , |y| < 0.5\n");
  printf("---------------------------------------------------\n");

  // ======================================================================================
  // ======================================================================================
  for (int idxN = 0; idxN < 3; idxN++) {
    TList *list = static_cast<TList*>(file->Get(aNames[idxN]));
    printf("\n============================================================\n");
    for (int idxM = 0; idxM < 3; idxM++) {
      printf("\n---------------------------------------------------\n");
      printf("-- %s vs sqrt(s_NN) : %s\n", aNames[idxN], aTitle[idxM]);
      printf("---------------------------------------------------\n");
      
      for (int idxC = 0; idxC < 2; idxC++) {
	if (idxC == 1 && idxN != 2)
	  continue;
	TGraphErrors *gStat = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_stat", aNames[idxN], aMom[idxM], aCent[idxC])));
	TGraphErrors *gSys = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_sys",   aNames[idxN], aMom[idxM], aCent[idxC])));

	printf("\ncent ; sqrt(sNN);    value ;     stat ;      sys \n");
	for (Int_t idx = 0 ; idx < gStat->GetN(); ++idx) {
	  Double_t snn, value;
	  gStat->GetPoint(idx, snn, value);
	  Double_t eStat = gStat->GetErrorY(idx);	
	  Double_t eSys = gSys->GetErrorY(idx);	
	  printf("%s ; %8.1f ; %8f ; %8f ; %8f \n",  aCent[idxC], snn, value, eStat, eSys);
	} // for (Int_t idx = 0 ; idx < gStat->GetN(); ++idx) {
      } // for (int idxC = 0; idxC < 2; idxC++) {
    } // for (int idxM = 0; idxM < 3; idxM++) {
  } // for (int idxN = 0; idxN < 3; idxN++) {
  
  // ======================================================================================
  // ======================================================================================
  TList *list = static_cast<TList*>(file->Get("Net-Charge_VsDeltaEta"));

  printf("\n============================================================\n");
  for (int idxM = 0; idxM < 3; idxM++) {
    printf("\n---------------------------------------------------\n");
    printf("-- NetCharge (sqrt(s_NN) = 14.5 GeV) vs #DeltaEta : %s\n", aTitle2[idxM]);
    printf("   - statistical errors only\n");
    printf("---------------------------------------------------\n");
    
    for (int idxC = 0; idxC < 2; idxC++) {
      TGraphErrors *gStat = static_cast<TGraphErrors*>(list->FindObject(Form("Net-Charge_%s_DeltaEta_14.5GeV_%s_stat", aMom2[idxM], aCent[idxC])));
      
      printf("\ncent ; deltaEta;    value ;     stat\n");
      for (Int_t idx = 0 ; idx < gStat->GetN(); ++idx) {
	Double_t deltaEta, value;
	gStat->GetPoint(idx, deltaEta, value);
	Double_t eStat = gStat->GetErrorY(idx);	
	printf("%s ; %8f ; %8f ; %8f\n",  aCent[idxC], deltaEta, value, eStat);
      } // for (Int_t idx = 0 ; idx < gStat->GetN(); ++idx) {
    } // for (int idxC = 0; idxC < 2; idxC++) {
  } // for (int idxM = 0; idxM < 3; idxM++) {

  // ======================================================================================
  // ======================================================================================
  list = static_cast<TList*>(file->Get("Net-Charge_VsNpart"));

  printf("\n============================================================\n");
  for (int idxM = 0; idxM < 3; idxM++) {
    printf("\n---------------------------------------------------\n");
    printf("-- NetCharge (sqrt(s_NN) = 14.5 GeV) vs <N_part> : %s\n", aTitle2[idxM]);
    printf("   - statistical errors only\n");
    printf("---------------------------------------------------\n");
    
    for (int idxE = 0; idxE < 5; idxE++) {
      TGraphErrors *gStat = static_cast<TGraphErrors*>(list->FindObject(Form("Net-Charge_%s_Npart_14.5GeV_eta_%.1f_stat", aMom2[idxM], aEta[idxE])));
      printf("\ndeltaEta ; <N_part>;    value ;     stat\n");
      for (Int_t idx = 0 ; idx < gStat->GetN(); ++idx) {
	Double_t nPart, value;
	gStat->GetPoint(idx, nPart, value);
	Double_t eStat = gStat->GetErrorY(idx);	
	printf("%5.1f ; %8.1f ; %8f ; %8f\n",  aEta[idxE], nPart, value, eStat);
      } // for (Int_t idx = 0 ; idx < gStat->GetN(); ++idx) {
    } // for (int idxC = 0; idxC < 2; idxC++) {
  } // for (int idxM = 0; idxM < 3; idxM++) {
}

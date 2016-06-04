#include "include/toolsEnergyNice.C"
#include "include/getPublished.C"
#include "Net-Kaon/netKaon_sysErrors.h"

// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

  Float_t aMinYtmp[7] = { -5, -5, -5, -2300, -2, -0.3, -3.5};
  Float_t aMaxYtmp[7] = { 40, 250, 100, 900, 45, 2.3, 2.8};

  for (Int_t ii = 0; ii<7; ++ii) {
    aMinY[ii] = aMinYtmp[ii];
    aMaxY[ii] = aMaxYtmp[ii];
  }
}

// ______________________________________________________________________________________
void plotEnergyKaon(const Char_t* name = "ratioNetKaonVsEnergy") {

  Int_t idxNames = kNetK;

  gROOT->LoadMacro("include/toolsEnergyNice.C++");
  gROOT->LoadMacro("include/getPublished.C++");

  SetupStyle();

  SetGlobals();

  getPublished();

  // -----------------------------------------------------

  TFile *inFiles[nEnergies][nMoments];
  TFile *inFilesUrqmd[nEnergies];

  TGraphErrors *inGraphsStat[nEnergies][nMoments+1];
  TGraphErrors *inGraphsPoisson[nEnergies][nMoments];
  TGraphErrors *inGraphsUrqmd[nEnergies][nMoments];

  for (int idxEnergy = 0; idxEnergy < nEnergies; ++idxEnergy) { 

    if (idxEnergy != 2)
      inFilesUrqmd[idxEnergy]= TFile::Open(Form("URQMD/urqmd_kaon/AuAu%sGeV_Vz30_kpi1.root", exactEnergies[idxEnergy]));

    for (int idxMoment = 0; idxMoment < nMoments; ++idxMoment) { 
      // -- value and stat errors
      inFiles[idxEnergy][idxMoment] = TFile::Open(Form("%s/%sGeV/moments_ana/Moments_%s.root", 
						       aNames[idxNames], exactEnergies[idxEnergy], aMoments[idxMoment]));

      inGraphsStat[idxEnergy][idxMoment] = (idxMoment != 5) ? static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(aMoments[idxMoment]))->Clone()) :
      	static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(Form("%s_Poisson_ratio", aMoments[idxMoment])))->Clone());

      if (idxMoment == 5) // get SK also
	inGraphsStat[idxEnergy][nMoments] = static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(aMoments[idxMoment]))->Clone());

      // -- poisson
      inGraphsPoisson[idxEnergy][idxMoment] = static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(Form("%s_Poisson_base", aMoments[idxMoment])))->Clone());
      
      if (inFiles[idxEnergy][idxMoment])
       	inFiles[idxEnergy][idxMoment]->Close();
      
      // -- urqmd
      if (idxEnergy != 2) {
       	if (idxMoment == 4) 
       	  inGraphsUrqmd[idxEnergy][idxMoment] = static_cast<TGraphErrors*>((inFilesUrqmd[idxEnergy]->Get(Form("R21")))->Clone());
	else if (idxMoment == 5) 
       	  inGraphsUrqmd[idxEnergy][idxMoment] = static_cast<TGraphErrors*>((inFilesUrqmd[idxEnergy]->Get(Form("%s_Pos_ratio", aMoments[idxMoment])))->Clone());
       	else if (idxMoment == 6) 
       	  inGraphsUrqmd[idxEnergy][idxMoment] = static_cast<TGraphErrors*>((inFilesUrqmd[idxEnergy]->Get(aMoments[idxMoment]))->Clone());
      }
    }
  }

  // -----------------------------------------------------
  // -- Make graphs
  
  for (int idxEnergy = 0; idxEnergy < nEnergies; ++idxEnergy) { 
    for (int idxMoment = 4; idxMoment < nMoments+1; ++idxMoment) {   
      for (int idxCent = 0; idxCent < nCent; ++idxCent) {
	if (idxCent != 0 && idxCent != 8)
	  continue;

	Double_t xIn, yIn;
	inGraphsStat[idxEnergy][idxMoment]->GetPoint(idxCent, xIn, yIn);
	Double_t yErrorStatIn = inGraphsStat[idxEnergy][idxMoment]->GetErrorY(idxCent);

	Int_t accessCent = (idxCent == 0) ? 0 : 1;
	Int_t accessMoment = idxMoment - 4;
	Double_t yErrorSysIn  = sysErrors[accessMoment][accessCent][idxEnergy];

	for (Int_t idx = 0; idx < 2; ++idx) {
	  Double_t xBin = (idx == 0) ? snn[idxEnergy] : mub[idxEnergy];
	  graphStat[idx][idxMoment][idxCent]->SetPoint(idxEnergy, xBin, yIn);
	  graphStat[idx][idxMoment][idxCent]->SetPointError(idxEnergy, 0, yErrorStatIn);
	  graphSys[idx][idxMoment][idxCent]->SetPoint(idxEnergy, xBin, yIn);
	  graphSys[idx][idxMoment][idxCent]->SetPointError(idxEnergy, 0, yErrorSysIn);
	}
	
	if (idxMoment == nMoments)
	  continue;

	inGraphsPoisson[idxEnergy][idxMoment]->GetPoint(idxCent, xIn, yIn);    
	for (Int_t idx = 0; idx < 2; ++idx) {
	  Double_t xBin = (idx == 0) ? snn[idxEnergy] : mub[idxEnergy];
	  graphPoisson[idx][idxMoment][idxCent]->SetPoint(idxEnergy, xBin, yIn);
	}

	if (idxCent == 0 && idxEnergy != 2) {
	  Double_t yErrorUrqmdIn = inGraphsUrqmd[idxEnergy][idxMoment]->GetErrorY(idxCent);	

	  inGraphsUrqmd[idxEnergy][idxMoment]->GetPoint(idxCent, xIn, yIn);    
	  for (Int_t idx = 0; idx < 2; ++idx) {
	    Double_t xBin = (idx == 0) ? snn[idxEnergy] : mub[idxEnergy];
	    Int_t energyArrayAccess = (idxEnergy >= 2) ? idxEnergy-1 : idxEnergy;
	    graphUrqmd[idx][idxMoment][idxCent]->SetPoint(energyArrayAccess, xBin, yIn);
	    graphUrqmd[idx][idxMoment][idxCent]->SetPointError(energyArrayAccess, 0, yErrorUrqmdIn);
	  }
	}
      } // for (int idxCent = 0; idxCent < nCent; ++idxCent) {
    } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {   
  } // for (int idxEnergy = 0 ; idxEnergy < nEnergies; ++idxEnergy) { 

  // -----------------------------------------------------

  SetupCanvas(name, Form("%s Ratio energy dependence", aNames[idxNames]));
  CreateLegends(2, 3, 0.4, 0.09);

  // -----------------------------------------------------

  for (int idxMoment = 4; idxMoment < nMoments; ++idxMoment) {
    pad->cd(idxMoment-3);
    gPad->SetLogx();

    for (int idxCent = 0; idxCent < nCent; ++idxCent) {
      if (idxCent != 0 && idxCent != 8)
	continue;

      DrawSet(graphStat[0][idxMoment][idxCent],  graphSys[0][idxMoment][idxCent],
	      graphUrqmd[0][idxMoment][idxCent], graphPoisson[0][idxMoment][idxCent],
	      idxMoment, idxCent);
    } // for (int idxCent = 0; idxCent < nCent; ++idxCent) {

    graphStat[0][idxMoment][0]->Draw("ZP,SAME");
    graphSys[0][idxMoment][0]->Draw("[],SAME");
  } // for (int idxMoment = 4; idxMoment < nMoments; ++idxMoment) {

  legTheo->AddEntry(graphUrqmd[0][4][0], Form("%s UrQMD", cent1[0]), "f");
      
  // -----------------------------------------------------

  LabelCanvas(aNames[idxNames], aNamesPt[idxNames]);
  SaveCanvas(name);

  // -----------------------------------------------------

  TFile *fOut = TFile::Open("STAR_QM2015_Preliminary.root", "UPDATE");
  fOut->cd();
  
  TList* list = new TList;
  
  for (int idxMoment = 4; idxMoment < nMoments+1; ++idxMoment) {
    for (int idxCent = 0; idxCent < nCent; ++idxCent) {
      if (idxCent != 0) 
	continue;
      
      graphStat[0][idxMoment][idxCent]->SetName(Form("%s_%s_sNN_%s_stat", aNames[idxNames], aMoments2[idxMoment], cent[idxCent]));
      graphSys[0][idxMoment][idxCent]->SetName(Form("%s_%s_sNN_%s_sys",  aNames[idxNames], aMoments2[idxMoment], cent[idxCent]));
      
      list->Add(graphStat[0][idxMoment][idxCent]);
      list->Add(graphSys[0][idxMoment][idxCent]);
    }
  }
  list->Write(aNames[idxNames], TObject::kSingleKey);
  fOut->Close();

  // -----------------------------------------------------

  TFile *fOutAll = TFile::Open("STAR_Preliminary.root", "UPDATE");
  fOutAll->cd();
  
  TList* listAll = new TList;

  for (int idxMoment = 4; idxMoment < nMoments+1; ++idxMoment) {
    for (int idxCent = 0; idxCent < nCent; ++idxCent) {
      if (graphStat[0][idxMoment][idxCent] && graphStat[0][idxMoment][idxCent]->GetN() > 0) {
	graphStat[0][idxMoment][idxCent]->SetName(Form("%s_%s_sNN_%s_stat", aNames[idxNames], aMoments2[idxMoment], cent[idxCent]));
	listAll->Add(graphStat[0][idxMoment][idxCent]);
      }
      if (graphSys[0][idxMoment][idxCent] && graphSys[0][idxMoment][idxCent]->GetN() > 0) {
	graphSys[0][idxMoment][idxCent]->SetName(Form("%s_%s_sNN_%s_sys", aNames[idxNames], aMoments2[idxMoment], cent[idxCent]));
	listAll->Add(graphSys[0][idxMoment][idxCent]);
      }

      if (idxMoment == nMoments)
	continue;

      if (graphUrqmd[0][idxMoment][idxCent] && graphUrqmd[0][idxMoment][idxCent]->GetN() > 0) {
	graphUrqmd[0][idxMoment][idxCent]->SetName(Form("%s_%s_sNN_%s_urqmd", aNames[idxNames], aMoments2[idxMoment], cent[idxCent]));
	listAll->Add(graphUrqmd[0][idxMoment][idxCent]);
      }
      if (graphPoisson[0][idxMoment][idxCent] && graphPoisson[0][idxMoment][idxCent]->GetN() > 0) {
	graphPoisson[0][idxMoment][idxCent]->SetName(Form("%s_%s_sNN_%s_poisson", aNames[idxNames], aMoments2[idxMoment], cent[idxCent]));
	listAll->Add(graphPoisson[0][idxMoment][idxCent]);
      }
    }
  }
  listAll->Write(aNames[idxNames], TObject::kSingleKey);
  fOutAll->Close();
}

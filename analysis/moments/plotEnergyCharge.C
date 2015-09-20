
#include "include/toolsEnergyNice.C"
#include "include/getPublished.C"


// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

  // Float_t aMinYtmp[7] = { -5, -5, -5, -2300, 3.5, -0.02, -12 };
  // Float_t aMaxYtmp[7] = { 40, 250, 100, 900, 21,   0.62, 6};

  // for (Int_t ii = 0; ii<7; ++ii) {
  //   aMinY[ii] = aMinYtmp[ii];
  //   aMaxY[ii] = aMaxYtmp[ii];
  // }
}

// ______________________________________________________________________________________
void plotEnergyCharge(const Char_t* name = "ratioNetChargeVsEnergy") {

  gROOT->LoadMacro("include/toolsEnergyNice.C++");
  gROOT->LoadMacro("include/getPublished.C++");

  SetupStyle();

  SetGlobals();

  getPublished();

  // -----------------------------------------------------

  TFile *inFiles[nEnergies][nMoments];
  TFile *inFilesSys[nEnergies];

  TGraphErrors *inGraphsStat[nEnergies][nMoments];
  TGraphErrors *inGraphsSys[nEnergies][nMoments];
  TGraphErrors *inGraphsPoisson[nEnergies][nMoments];

  for (int idxEnergy  = 0 ; idxEnergy < nEnergies; ++idxEnergy) { 
    if (idxEnergy != 2) 
      continue;

    inFilesSys[idxEnergy] = TFile::Open(Form("sysError/jobs_%s_sysError.root", exactEnergies[idxEnergy]));
    inFilesSys[idxEnergy]->ls();

    for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) { 
      // -- value and stat errors
      inFiles[idxEnergy][idxMoment]  = TFile::Open(Form("output/jobs_14.5_base/%s/Moments_%s.root", aDataSets[0], aMoments[idxMoment]));

      inGraphsStat[idxEnergy][idxMoment] = (idxMoment != 5) ? static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(aMoments[idxMoment]))->Clone()) :
	static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(Form("%s_Poisson_ratio", aMoments[idxMoment])))->Clone());
      
      // -- poisson
      inGraphsPoisson[idxEnergy][idxMoment] = static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(Form("%s_Poisson_base", aMoments[idxMoment])))->Clone());
      
      if (inFiles[idxEnergy][idxMoment])
	inFiles[idxEnergy][idxMoment]->Close();
      
      // -- sysErrors
      // inGraphsSyS[idxEnergy][idxMoment] = (idxMoment != 5) ? static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(aMoments[idxMoment]))->Clone()) :
      // 	static_cast<TGraphErrors*>((inFiles[idxEnergy][idxMoment]->Get(Form("%s_Poisson_ratio", aMoments[idxMoment])))->Clone());
      

    }
  }
        
  // -----------------------------------------------------
  // -- Add 14.5 Data Points
  
  Int_t idxEnergy = 2;


  for (int idxMoment = 4 ; idxMoment < nMoments; ++idxMoment) {   
    for (int idxCent = 0; idxCent < nCent; ++idxCent) {
      if (idxCent != 0 && idxCent != 8)
	continue;
      
      Double_t xIn, yIn;
      inGraphsStat[idxEnergy][idxMoment]->GetPoint(idxCent, xIn, yIn);
      Double_t yErrorIn = inGraphsStat[idxEnergy][idxMoment]->GetErrorY(idxCent);

      for (Int_t idx = 0; idx < 2; ++idx) {
	Double_t xOut, yOut;
	graphStat[idx][idxMoment][idxCent]->GetPoint(idxEnergy, xOut, yOut);
	graphStat[idx][idxMoment][idxCent]->SetPoint(idxEnergy, xOut, yIn);
	graphStat[idx][idxMoment][idxCent]->SetPointError(idxEnergy, 0, yErrorIn);
      }

      inGraphsPoisson[idxEnergy][idxMoment]->GetPoint(idxEnergy, xIn, yIn);    
      for (Int_t idx = 0; idx < 2; ++idx) {
	Double_t xOut, yOut;
	graphPoisson[idx][idxMoment][idxCent]->GetPoint(idxEnergy, xOut, yOut);
	graphPoisson[idx][idxMoment][idxCent]->SetPoint(idxEnergy, xOut, yIn);
      }
    } // for (int idxCent = 0; idxCent < nCent; ++idxCent) {
  } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {   

    

  // -----------------------------------------------------


  TLegend *legRat = new TLegend(0.12, 0.12, 0.70, 0.57);
  legRat->SetTextAlign(12);
  legRat->SetTextSize(0.06);
  legRat->SetTextFont(42);
  legRat->SetFillColor(0);
  legRat->SetLineColor(0);
  legRat->SetBorderSize(0);

  // -----------------------------------------------------
  
  TPad *pad = SetupCanvas("canNetChargeRatioEnergy", "Net-Charge Ratio energy dependence");

  for (int idxMoment = 4; idxMoment < nMoments; ++idxMoment) {
    pad->cd(idxMoment-3);
    gPad->SetLogx();

    for (int idxCent = 0; idxCent < nCent; ++idxCent) {
      if (idxCent != 0 && idxCent != 8)
	continue;

      TGraphErrors *gStat    = graphStat[0][idxMoment][idxCent];
      TGraphErrors *gSys     = graphSys[0][idxMoment][idxCent];
      TGraphErrors *gPoisson = graphPoisson[0][idxMoment][idxCent];

      PrepareGraph(gStat);
      PrepareGraph(gSys);
      PrepareGraph(gPoisson);

      ConfigGraph(gStat, idxMoment, idxCent);
      ConfigGraph(gSys, idxMoment, idxCent);
      ConfigGraph(gPoisson, idxMoment, idxCent, 1);

      if (idxCent == 0) 
	gStat->Draw("AP");

      if (idxMoment == 5 || idxMoment == 6) {
	if (idxCent == 0) {
	  TLine *line1 = new TLine(aMinX, 1, aMaxX, 1);
	  line1->SetLineColor(kGray+1);
	  line1->SetLineStyle(2);
	  line1->SetLineWidth(2);
	  line1->Draw();
	}
      }
      else if (idxMoment == 4)
	gPoisson->Draw("L,SAME");

      gStat->Draw("P,SAME");
      gSys->Draw("[],SAME");

      // gSys->SetFillColor(kBlue);
      // gSys->Draw("2,SAME");

      if (idxMoment == 4) {
	// legRat->AddEntry(gStat, Form("%s", cent[idxCent]), "pl");
	// legRat->AddEntry(gPoisson, Form("%s Poisson", cent[idxCent]), "pl");
      }
    } // for (int idxCent = 0; idxCent < nCent; ++idxCent) {


  } // for (int idxMoment = 4; idxMoment < nMoments; ++idxMoment) {
      
      
  pad->cd(1);
  TLatex *texb_3 = new TLatex(130,19, "Au+Au collisions #sqrt{#it{s}_{NN}} = 14.5 GeV");
  texb_3->SetTextSize(0.07);
  texb_3->Draw("same");
  
  TLatex *texb_3a = new TLatex(130,17.5, "Net-Charge, 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0");
  texb_3a->SetTextSize(0.07);
  texb_3a->Draw("same");
  
  pad->cd(2);
  TLatex *texb_4 = new TLatex(15, 0.55, "STAR Preliminary");
  texb_4->SetTextSize(0.07);
  texb_4->Draw("same");
  
  pad->cd(3);
  //  legRat->Draw("lt");
  
  pad->Modified();
  
  // -----------------------------------------------------
  
  SaveCanvas(name);
}

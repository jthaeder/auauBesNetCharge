#include "include/toolsEnergyNice.C"
#include "include/getPublished.C"
#include "Net-Kaon/netKaon_sysErrors.h"

// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

  aMinX = 0;
  aMaxX = 380;

  aMarkers[2] = 30;
  aMarkers[1] = 24;
  aMarkers[0] = 25;

  aColors[1] = kAzure;
  aColors[0] = kBlack;
}

// ______________________________________________________________________________________
void plotCumulants(const Char_t* name = "cumulants") {

  gROOT->LoadMacro("include/toolsEnergyNice.C++");
  gROOT->LoadMacro("include/getPublished.C++");

  SetupStyle();

  SetGlobals();

  getPublished();

    // -----------------------------------------------------

  TFile *inFiles[3][nMoments];
  TFile *inFilesUnCorr[3][nMoments];
  TFile *inFilesSys[3];

  TGraphErrors *inGraphsStat[3][nMoments];
  TGraphErrors *inGraphsSys[3][nMoments];
  TGraphErrors *inGraphsUnCorr[3][nMoments];

  Int_t idxEnergy = 2;

  inFilesSys[0] = TFile::Open(Form("sysError/jobs_%s_sysError.root", exactEnergies[idxEnergy]));

  for (int idxMoment = 0 ; idxMoment < 4; ++idxMoment) { 

    // -- Net-Charge
    inFiles[0][idxMoment]       = TFile::Open(Form("output/jobs_%s_base/twoeff_11/Moments_%s.root", exactEnergies[idxEnergy], aMoments[idxMoment]));
    inFilesUnCorr[0][idxMoment] = TFile::Open(Form("output/jobs_%s_base/effuncorr/Moments_%s.root", exactEnergies[idxEnergy], aMoments[idxMoment]));
    
    inGraphsStat[0][idxMoment]   = static_cast<TGraphErrors*>((inFiles[0][idxMoment]->Get(aMoments[idxMoment]))->Clone());
    inGraphsUnCorr[0][idxMoment] = static_cast<TGraphErrors*>((inFilesUnCorr[0][idxMoment]->Get(aMoments[idxMoment]))->Clone());
    inGraphsSys[0][idxMoment]    = static_cast<TGraphErrors*>((inFilesSys[0]->Get(Form("%s_sys",aMoments[idxMoment])))->Clone());

    // -- Net-Kaon
    inFiles[1][idxMoment]       = TFile::Open(Form("Net-Kaon/%sGeV/moments_ana/Moments_%s.root", exactEnergies[idxEnergy], aMoments[idxMoment]));
    inFilesUnCorr[1][idxMoment] = TFile::Open(Form("Net-Kaon/%sGeV/moments_uncorr/Moments_%s.root", exactEnergies[idxEnergy], aMoments[idxMoment]));

    inGraphsStat[1][idxMoment]   = static_cast<TGraphErrors*>((inFiles[1][idxMoment]->Get(aMoments[idxMoment]))->Clone());
    inGraphsUnCorr[1][idxMoment] = static_cast<TGraphErrors*>((inFilesUnCorr[1][idxMoment]->Get(aMoments[idxMoment]))->Clone());
    inGraphsSys[1][idxMoment]    = static_cast<TGraphErrors*>((inFiles[1][idxMoment]->Get(aMoments[idxMoment]))->Clone());
    for (Int_t idxCent = 0; idxCent < 9; ++idxCent) 
      inGraphsSys[1][idxMoment]->SetPointError(idxCent, 0, sysError_kaon[idxMoment][idxCent]);

    // -- Net-Proton
    if (idxMoment == 0) {
      inFiles[2][0]       = TFile::Open(Form("Net-Proton/moments_%sGeV_preli.root", energies[idxEnergy]));
    }

    inGraphsStat[2][idxMoment]   = static_cast<TGraphErrors*>((inFiles[2][0]->Get(Form("%s_stat", aMoments[idxMoment])))->Clone());
    inGraphsSys[2][idxMoment]    = static_cast<TGraphErrors*>((inFiles[2][0]->Get(Form("%s_sys", aMoments[idxMoment])))->Clone());
    inGraphsUnCorr[2][idxMoment] = static_cast<TGraphErrors*>((inFiles[2][0]->Get(Form("%s_noEff_stat", aMoments[idxMoment])))->Clone());
  } // for (int idxMoment = 0 ; idxMoment < 4; ++idxMoment) { 

  // -----------------------------------------------------

  TGraphErrors *inGraphsRatio[3][nMoments];
  
  for (int idxMoment = 0 ; idxMoment < 4; ++idxMoment) { 
    for (Int_t idx = 0; idx < 3; ++idx) {
      TGraphErrors *g1 = inGraphsStat[idx][idxMoment];
      Double_t* aX   = g1->GetX(); 
      Double_t* g1Y  = g1->GetY(); 
      Double_t* g1EY = g1->GetEY(); 
      
      TGraphErrors *g2 = inGraphsUnCorr[idx][idxMoment];
      Double_t* g2Y  = g2->GetY(); 
      Double_t* g2EY = g2->GetEY(); 
      
      Double_t g1REY[9];
      Double_t g2REY[9];
      
      Double_t aY[9];
      Double_t aREY[9];
      Double_t aEY[9];

      for (int bin = 0; bin <9; bin++) {
	aY[bin]  = g2Y[bin] / g1Y[bin];
      
	g1REY[bin] = g1EY[bin]/g1Y[bin];
	g2REY[bin] = g2EY[bin]/g2Y[bin];
	
	aREY[bin] = sqrt( (g1REY[bin]*g1REY[bin])+(g2REY[bin]*g2REY[bin]));
	//	aEY[bin]  = sqrt( (g1REY[bin]/g1Y[bin])*(g1REY[bin]/g1Y[bin]) + g2EY[bin]*g2EY[bin]);
	aEY[bin]  = aREY[bin]*aY[bin];
      }

      inGraphsRatio[idx][idxMoment] = new TGraphErrors(9, aX, aY, 0, aEY);
      
    } // for (int idxMoment = 0 ; idxMoment < 4; ++idxMoment) { 
  }  // for (Int_t idx = 0; idx < 3; ++idx) {

  // -----------------------------------------------------

  canA.Add(new TCanvas(Form("can%s", name), "Cumulants", 0, 0 , 1200, 500));
  can = static_cast<TCanvas*>(canA.Last());
  can->SetFillColor(1182);
  can->SetBorderMode(0);
  can->SetBorderSize(0.0);
  can->SetFrameFillColor(1182);
  can->SetFrameBorderMode(0);
  can->cd();
  
  pad = new TPad("pad", "pad", 0.02, 0.0, 1., 1.);
  pad->SetTopMargin(0.);
  pad->SetBottomMargin(0.);
  pad->SetLeftMargin(0.);
  pad->SetRightMargin(0.);
  pad->SetBorderMode(0);
  pad->SetFillColor(1182);
  pad->Draw();
  pad->cd();


  TPad *padSet0[4];
  TPad *padSet1[4][2];
  TPad *padSet2[4];

  legExp = new TLegend(0.25, 0.45, 0.8, 0.65);
  legExp->SetTextAlign(12);
  legExp->SetTextSize(0.05);
  legExp->SetTextFont(42);
  legExp->SetFillColor(1182);
  legExp->SetLineColor(0);
  legExp->SetBorderSize(0);
  
  for (Int_t idxMoment = 0; idxMoment < 4; ++idxMoment) {
    pad->cd();
    padSet0[idxMoment] = new TPad(Form("pad_%d", idxMoment), "", idxMoment*0.25, 0.0, (idxMoment+1)*0.25, 1.);
    padSet0[idxMoment]->SetTopMargin(0.);
    padSet0[idxMoment]->SetBorderMode(0);
    padSet0[idxMoment]->SetRightMargin(0.);
    padSet0[idxMoment]->SetLeftMargin(0.);
    padSet0[idxMoment]->SetFillColor(1182);
    padSet0[idxMoment]->Draw();
    padSet0[idxMoment]->cd();

    padSet1[idxMoment][1] = new TPad(Form("pad_%d_1", idxMoment), "", 0.0, 0.4, 1.0, 1.0);
    padSet1[idxMoment][1]->SetTopMargin(0.1);
    padSet1[idxMoment][1]->SetBottomMargin(0.);
    padSet1[idxMoment][1]->SetRightMargin(0.01);
    padSet1[idxMoment][1]->SetBorderMode(0);
    padSet1[idxMoment][1]->SetFillColor(1182);
    padSet1[idxMoment][1]->Draw();

    padSet1[idxMoment][0] = new TPad(Form("pad_%d_0", idxMoment), "", 0.0, 0.06, 1.0, 0.4);
    padSet1[idxMoment][0]->SetTopMargin(0.);
    padSet1[idxMoment][0]->SetBorderMode(0);
    padSet1[idxMoment][0]->SetBottomMargin(0.1);
    padSet1[idxMoment][0]->SetRightMargin(0.01);
    padSet1[idxMoment][0]->SetFillColor(1182);
    padSet1[idxMoment][0]->Draw();

    padSet0[idxMoment]->cd();
    TLatex *texb_5 = new TLatex(0.48, 0.0175, "#LT#it{N}_{part}#GT");
    texb_5->SetTextSize(0.07);
    texb_5->Draw("same");
  }
  
  for (Int_t idxMoment = 0; idxMoment < 4; ++idxMoment) {
    for (Int_t idx = 0; idx < 3; ++idx) {
      padSet1[idxMoment][1]->cd();

      TGraphErrors* gStat      = inGraphsStat[idx][idxMoment];
      TGraphErrors* gSys       = inGraphsSys[idx][idxMoment];
      TGraphErrors* gRatioStat = inGraphsRatio[idx][idxMoment];
      
      PrepareGraphCumulants(gStat);
      ConfigGraph(gStat, idxMoment, idx);
      
      if (idx == 0)
	gStat->Draw("AZP");
      else
	gStat->Draw("ZP,SAME");

      PrepareGraph(gSys);
      ConfigGraph(gSys, idxMoment, idx, 2);
      gSys->Draw("[], SAME");
      
      padSet1[idxMoment][0]->cd();

      PrepareGraphCumulants(gRatioStat);
      ConfigGraph(gRatioStat, idxMoment, idx, 4);	
      gRatioStat->GetXaxis()->SetLabelSize(0.09);
      gRatioStat->GetYaxis()->SetLabelSize(0.09);
      gRatioStat->GetYaxis()->SetNdivisions(5, 5, 0);
      
      if (idx == 0)
	gRatioStat->Draw("AP");
      else
	gRatioStat->Draw("P,SAME");
    
    } // for (Int_t idx = 0; idx < 3; ++idx) {
  } // for (Int_t idxMoment = 0; idxMoment < 4; ++idxMoment) {

  // -----------------------------------------------------  

  legExp->AddEntry(inGraphsStat[0][0], "Net-Charge", "pl"); 
  legExp->AddEntry(inGraphsStat[1][0], "Net-Kaon", "pl"); 
  legExp->AddEntry(inGraphsStat[2][0], "Net-Proton", "pl"); 

  // -----------------------------------------------------

  for (Int_t idxMoment = 0; idxMoment < 4; ++idxMoment) {
    pad->cd();
    padSet2[idxMoment] = new TPad(Form("padCover_%d", idxMoment), "", idxMoment*0.25, 0.0, (idxMoment+1)*0.25, 1.);
    padSet2[idxMoment]->SetTopMargin(0.);
    padSet2[idxMoment]->SetBorderMode(0);
    padSet2[idxMoment]->SetRightMargin(0.);
    padSet2[idxMoment]->SetLeftMargin(0.);
    padSet2[idxMoment]->SetFillColor(1182);
    padSet2[idxMoment]->Draw();
    padSet2[idxMoment]->cd();
 
    TLatex *texb_6c = new TLatex(0.3,0.70, aMomentsTitle[idxMoment]);
    texb_6c->SetTextSize(0.1);
    texb_6c->Draw("same");
  }

  // -----------------------------------------------------

  can->cd();
  TLatex *texb_6b = new TLatex(0.017,0.13, "Corr./Uncorr.");
  texb_6b->SetTextSize(0.045);
  texb_6b->SetTextAngle(90);
  texb_6b->Draw("same");

  pad->cd();
  TLatex *texb_9 = new TLatex(0.03, 0.965, "Net-Charge: 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0, |#eta| < 0.5,  Net-Kaon: 0.2 < #it{p}_{T} (GeV/#it{c}) < 1.6, |#it{y}| < 0.5,  Net-Proton: 0.4 < #it{p}_{T} (GeV/#it{c}) < 2.0, |#it{y}| < 0.5");
  texb_9->SetTextSize(0.03);
  texb_9->Draw("same");

  padSet2[0]->cd();
  TLatex *texb_3 = new TLatex(0.15, 0.88, "Au+Au #sqrt{#it{s}_{NN}} =  14.5 GeV");
  texb_3->SetTextSize(0.065);
  texb_3->Draw("same");
  
  padSet2[1]->cd();
  TLatex *texb_4 = new TLatex(0.15, 0.88, "STAR Preliminary");
  texb_4->SetTextSize(0.065);
  texb_4->Draw("same");

  padSet2[3]->cd();
  legExp->Draw("lt");

  // -----------------------------------------------------

  SaveCanvas(name);
}

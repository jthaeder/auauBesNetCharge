
#include "include/toolsEtaNice.C"

const Char_t* aEtaSets[3][10]    = { {"2015-06-21_delta_0.1_0", "2015-05-20_0.1",
				      "2015-06-03_delta_0.3_8", "2015-05-20_0.2", 
				      "2015-06-20_delta_0.5_6", "2015-05-20_0.3",
				      "2015-06-20_delta_0.7_0", "2015-05-20_0.4", 
				      "2015-06-20_delta_0.9_0", "2015-05-20_0.5"},
				     {"jobs_14.5_chargeSep_1_-0.05_0.05", "jobs_14.5_chargeSep_1_-0.1_0.1",
				      "jobs_14.5_chargeSep_1_-0.15_0.15", "jobs_14.5_chargeSep_1_-0.2_0.2",
				      "jobs_14.5_chargeSep_1_-0.25_0.25", "jobs_14.5_chargeSep_1_-0.3_0.3",
				      "jobs_14.5_chargeSep_1_-0.35_0.35","jobs_14.5_chargeSep_1_-0.4_0.4",
				      "jobs_14.5_chargeSep_1_-0.45_0.45","jobs_14.5_chargeSep_1_-0.5_0.5"}, 
				     {"jobs_14.5_chargeSep_2_-0.05_0.05", "jobs_14.5_chargeSep_2_-0.1_0.1",
				      "jobs_14.5_chargeSep_2_-0.15_0.15", "jobs_14.5_chargeSep_2_-0.2_0.2",
				      "jobs_14.5_chargeSep_2_-0.25_0.25", "jobs_14.5_chargeSep_2_-0.3_0.3",
				      "jobs_14.5_chargeSep_2_-0.35_0.35","jobs_14.5_chargeSep_2_-0.4_0.4",
				      "jobs_14.5_chargeSep_2_-0.45_0.45","jobs_14.5_chargeSep_2_-0.5_0.5"} };

const char*   aSuperSetTitles[3] = { "", "pos: #eta>0, neg: #eta<0", "pos: #eta<0, neg: #eta>0"};  

const int     nEtaSets           = 10; 
const int     nEtaSuperSets      = 3; 

const float   aEtaSetBinCenter[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const float   aEtaSetBinWidth[]  = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

  nCent = 9;

  // Color_t aColorstmp[]  = {kRed+2, kAzure, kOrange+9, kYellow, kBlack, kAzure, kCyan+2, kBlue+1, kOrange+9, kOrange+9};
  // Int_t   aMarkerstmp[] = {30, 24, 25, 25, 25, 27, 27, 30, 27, 30};

  Color_t aColorstmp[]  = {kGray, kAzure, kRed+2,
			   kGray+2, kAzure+2, kRed+2, 
			   kCyan+2, kBlue+1, kOrange+9, 
			   kOrange+9, kOrange+9, kOrange+9};
  Int_t   aMarkerstmp[] = {20, 24, 24, 21, 25, 25, 
			   27, 30, 27, 30, 29, 29};

  for (Int_t ii = 0; ii < 12; ++ii) {
    aColors[ii] = aColorstmp[ii];
    aMarkers[ii] = aMarkerstmp[ii];
  }

  Float_t aMinYtmp[7] = { -5, -5, -5, -2300, 3.5, -0.02, -12 };
  Float_t aMaxYtmp[7] = { 40, 250, 100, 900, 52,   0.48, 6};

  for (Int_t ii = 0; ii<7; ++ii) {
    aMinY[ii] = aMinYtmp[ii];
    aMaxY[ii] = aMaxYtmp[ii];
  }

  aMinX = 0;
  aMaxX = 1.1;
}

// ______________________________________________________________________________________
void plotVsDeltaEta_chargeSep(const Char_t* name = "ratioVsDeltaEtaChargeSep") {

  gROOT->LoadMacro("include/toolsEtaNice.C++");

  SetupStyle();

  SetGlobals();

  // -----------------------------------------------------

  TFile *inFiles[nEtaSuperSets][nEtaSets][nDataSets][nMoments];

  TGraphErrors *graphs[nEtaSuperSets][nEtaSets][nDataSets][nMoments];

  for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) 
    for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) 
      for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) 
	for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) { 
	  inFiles[idxEtaSuper][idxEta][idxDataSet][idxMoment] = TFile::Open(Form("output/%s/%s/Moments_%s.root", 
								      aEtaSets[idxEtaSuper][idxEta], aDataSets[idxDataSet], aMoments[idxMoment]));
	  graphs[idxEtaSuper][idxEta][idxDataSet][idxMoment] = static_cast<TGraphErrors*>((inFiles[idxEtaSuper][idxEta][idxDataSet][idxMoment]->Get(aMoments[idxMoment]))->Clone());
	  if (inFiles[idxEtaSuper][idxEta][idxDataSet][idxMoment])
	    (inFiles[idxEtaSuper][idxEta][idxDataSet][idxMoment])->Close();
      }
  
  // -----------------------------------------------------
  
  TLegend *legRat[nDataSets][9];
  for (int idxCent = 0; idxCent < nCent; idxCent++) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      legRat[idxDataSet][idxCent] = new TLegend(0.12, 0.12, 0.70, 0.495);
      legRat[idxDataSet][idxCent]->SetTextAlign(12);
      legRat[idxDataSet][idxCent]->SetTextSize(0.06);
      legRat[idxDataSet][idxCent]->SetFillColor(1182);
      legRat[idxDataSet][idxCent]->SetLineColor(0);
      legRat[idxDataSet][idxCent]->SetBorderSize(0);
    }
  }

  TLegend *legCum[nDataSets][9];
  for (int idxCent = 0; idxCent < nCent; idxCent++) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      legCum[idxDataSet][idxCent] = new TLegend(0.12, 0.12, 0.70, 0.495);
      legCum[idxDataSet][idxCent]->SetTextAlign(12);
      legCum[idxDataSet][idxCent]->SetTextSize(0.06);
      legCum[idxDataSet][idxCent]->SetFillColor(1182);
      legCum[idxDataSet][idxCent]->SetLineColor(0);
      legCum[idxDataSet][idxCent]->SetBorderSize(0);
    }
  }


  // -----------------------------------------------------

  TGraphErrors *etaGraphs[nEtaSuperSets][9][nDataSets][nMoments];

  for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {
    for (int idxCent = 0; idxCent < nCent; idxCent++) {
      for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
	for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	  float y[nEtaSets];
	  float ey[nEtaSets];
	  for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) {	
	    y[idxEta]  = graphs[idxEtaSuper][idxEta][idxDataSet][idxMoment]->GetY()[idxCent];
	    ey[idxEta] = graphs[idxEtaSuper][idxEta][idxDataSet][idxMoment]->GetEY()[idxCent];
	  } 
	  
	  etaGraphs[idxEtaSuper][idxCent][idxDataSet][idxMoment] = new TGraphErrors(nEtaSets, aEtaSetBinWidth, y, 0, ey);
	  PrepareGraph(etaGraphs[idxEtaSuper][idxCent][idxDataSet][idxMoment]);
	} // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
      } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    } // for (int idxCent = 0; idxCent < 9; idxCent++) {
  } //  for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {

  // -----------------------------------------------------

  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {

    for (int idxCent = 0; idxCent < 9; idxCent++) {
      // if (idxCent != 0 && idxCent != 1 && idxCent != 8 && idxCent != 4 )
      // 	continue;
      
      TPad *padRat = SetupCanvas(Form("canDeltaEta_ChargeSep_Ratio_%s_%d", aDataSets[idxDataSet], idxCent), aDataSetsTitle[idxDataSet], "#Delta#eta", 0.48);
      TPad *padCum   = SetupCanvasCum(Form("canDeltaEta_ChargeSep_Cum_%s_%d",   aDataSets[idxDataSet], idxCent), aDataSetsTitle[idxDataSet], "#Delta#eta", 0.48);
      TPad *pad = NULL;

      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
	if (idxMoment < 4) {
	  pad = padCum;
	  pad->cd(idxMoment+1);
	}
	else {
	  pad = padRat;
	  pad->cd(idxMoment-3);
	}
	
	for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {      
	  TGraphErrors *g = etaGraphs[idxEtaSuper][idxCent][idxDataSet][idxMoment];

	  if (idxCent == 0) 
	    ShiftGraphX(g, -0.015);
	  else if (idxCent == 1) 
	    ShiftGraphX(g, -0.005);
	  else if (idxCent == 4) 
	    ShiftGraphX(g, 0.005);
	  else if (idxCent == 8) 
	    ShiftGraphX(g, 0.015);
	  
	  ConfigGraph(g, idxMoment, idxEtaSuper);

	  if (idxEtaSuper == 0) {
	    g->Draw("AP");
	    
	    if (idxMoment == 5) {
	      // TLine *line0 = new TLine(aMinX, 0, aMaxX, 0);
	      // line0->SetLineColor(kGray+1);
	      // line0->SetLineStyle(2);
	      // line0->SetLineWidth(2);
	      // line0->Draw();
	    }
	    else if (idxMoment == 6) {
	      TLine *line1 = new TLine(aMinX, 1, aMaxX, 1);
	      line1->SetLineColor(kGray+1);
	      line1->SetLineStyle(2);
	      line1->SetLineWidth(2);
	      line1->Draw();
	      legRat[idxDataSet][idxCent]->AddEntry(line1, "Poisson", "l");
	    }
	  }
	  
	  g->Draw("PSAME");
	  
	  if (idxMoment == 4) 
	    legRat[idxDataSet][idxCent]->AddEntry(etaGraphs[idxEtaSuper][idxCent][idxDataSet][4], Form("%s - %s", cent1[idxCent], aSuperSetTitles[idxEtaSuper]), "pl");
	  if (idxMoment == 0) 
	    legCum[idxDataSet][idxCent]->AddEntry(etaGraphs[idxEtaSuper][idxCent][idxDataSet][4], Form("%s - %s", cent1[idxCent], aSuperSetTitles[idxEtaSuper]), "pl");


	} // for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {      
      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
   
      padRat->cd(2);
      TLatex *texb_3 = new TLatex(0.05, 0.42, "Au+Au collisions #sqrt{#it{s}_{NN}} = 14.5 GeV");
      texb_3->SetTextSize(0.07);
      texb_3->Draw("same");
      
      TLatex *texb_3a = new TLatex(0.05, 0.38, "Net-Charge, 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0, 0-5%");
      texb_3a->SetTextSize(0.07);
      texb_3a->Draw("same");
      
      TLatex *texb_3b = new TLatex(0.05,0.34, "statistical errors only");
      texb_3b->SetTextSize(0.06);
      texb_3b->Draw("same");
      
      padRat->cd(1);
      TLatex *texb_4 = new TLatex(0.05, 47, "STAR Preliminary");
      texb_4->SetTextSize(0.07);
      texb_4->Draw("same");
      
      padRat->cd(3);
      legRat[idxDataSet][idxCent]->Draw();
    
      padRat->Modified();

      // --- --- --- --- --- --- --- --- --- --- --- --- 

      padCum->cd(3);
      TLatex *texc_3 = new TLatex(0.05, 82, "Au+Au collisions #sqrt{#it{s}_{NN}} = 14.5 GeV");
      texc_3->SetTextSize(0.07);
      texc_3->Draw("same");
      
      TLatex *texc_3a = new TLatex(0.05, 72, "Net-Charge, 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0, 0-5%");
      texc_3a->SetTextSize(0.07);
      texc_3a->Draw("same");
      
      TLatex *texc_3b = new TLatex(0.05, 62, "statistical errors only");
      texc_3b->SetTextSize(0.06);
      texc_3b->Draw("same");
      
      padCum->cd(1);
      TLatex *texc_4 = new TLatex(0.05, 35, "STAR Preliminary");
      texc_4->SetTextSize(0.07);
      texc_4->Draw("same");
      
      padCum->cd(4);
      legCum[idxDataSet][idxCent]->Draw();
    
      padCum->Modified();

    } // for (int idxCent = 0; idxCent < 9; idxCent++) { 
  } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  
  // -----------------------------------------------------

  SaveCanvas(name, kFALSE);

  // -----------------------------------------------------
  
  TFile *fOut = TFile::Open("STAR_QM2015_Preliminary.root", "UPDATE");
  fOut->cd();
  
  TList* list = new TList;

  for (int idxMoment = 4; idxMoment < nMoments; ++idxMoment) {
    for (int idxCent = 0; idxCent < nCent; ++idxCent) {
      if (idxCent > 1) 
	continue;
      for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   

	if (idxCent == 0) 
	  ShiftGraphX(etaGraphs[idxEtaSuper][idxCent][0][idxMoment], 0.015);
	else if (idxCent == 1) 
	  ShiftGraphX(etaGraphs[idxEtaSuper][idxCent][0][idxMoment], 0.005);
	else if (idxCent == 4) 
	  ShiftGraphX(etaGraphs[idxEtaSuper][idxCent][0][idxMoment], -0.005);
	else if (idxCent == 8) 
	  ShiftGraphX(etaGraphs[idxEtaSuper][idxCent][0][idxMoment], -0.015);
	
	etaGraphs[idxEtaSuper][idxCent][0][idxMoment]->SetName(Form("Net-Charge_%s_DeltaEta_ChargeSep_%d_14.5GeV_%s_stat", aMoments[idxMoment], idxEtaSuper, cent[idxCent]));
	list->Add(etaGraphs[idxEtaSuper][idxCent][0][idxMoment]);
      }
    }
  }
  
  list->Write("Net-Charge_VsDeltaEta_ChargeSep", TObject::kSingleKey);
  fOut->Close();
}

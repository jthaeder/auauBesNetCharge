
#include "include/toolsEtaNice.C"

const Char_t* aEtaSets[]         = {"2015-06-21_delta_0.1_0",
				    "2015-05-20_0.1",
				    "2015-06-03_delta_0.3_8",
				    "2015-05-20_0.2", 
				    "2015-06-20_delta_0.5_6", 
				    "2015-05-20_0.3", "2015-06-20_delta_0.7_0", "2015-05-20_0.4", "2015-06-20_delta_0.9_0", "2015-05-20_0.5"};  

const int     nEtaSets           = 10; 

const float   aEtaSetBinCenter[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const float   aEtaSetBinWidth[]  = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

  nCent = 9;

  Color_t aColorstmp[] = {kRed+2, kAzure, kOrange+9, kYellow, kBlack, kAzure, kCyan+2, kBlue+1, kOrange+9, kOrange+9};
  Int_t  aMarkerstmp[]   = {30, 24, 25, 25, 25, 27, 27, 30, 27, 30};

  for (Int_t ii = 0; ii < 10; ++ii) {
    aColors[ii] = aColorstmp[ii];
    aMarkers[ii] = aMarkerstmp[ii];
  }

  Float_t aMinYtmp[7] = { -5, -5, -5, -2300, 3.5, -0.02, -12 };
  Float_t aMaxYtmp[7] = { 40, 250, 100, 900, 21,   0.62, 6};

  for (Int_t ii = 0; ii<7; ++ii) {
    aMinY[ii] = aMinYtmp[ii];
    aMaxY[ii] = aMaxYtmp[ii];
  }

  aMinX = 0;
  aMaxX = 1.1;
}

// ______________________________________________________________________________________
void plotVsDeltaEta_nice(const Char_t* name = "ratioVsDeltaEta") {

  gROOT->LoadMacro("include/toolsEtaNice.C++");

  SetupStyle();

  SetGlobals();

  // -----------------------------------------------------

  TFile *inFiles[nEtaSets][nDataSets][nMoments];

  TGraphErrors *graphs[nEtaSets][nDataSets][nMoments];

  for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) 
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) 
      for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) { 
	inFiles[idxEta][idxDataSet][idxMoment] = TFile::Open(Form("output/%s/%s/Moments_%s.root", 
								  aEtaSets[idxEta], aDataSets[idxDataSet], aMoments[idxMoment]));
	graphs[idxEta][idxDataSet][idxMoment] = static_cast<TGraphErrors*>((inFiles[idxEta][idxDataSet][idxMoment]->Get(aMoments[idxMoment]))->Clone());
	if (inFiles[idxEta][idxDataSet][idxMoment])
	  (inFiles[idxEta][idxDataSet][idxMoment])->Close();
      }
  
  // -----------------------------------------------------
  
  TLegend *legRat[nDataSets];
  
  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    legRat[idxDataSet] = new TLegend(0.12, 0.12, 0.70, 0.495);
    legRat[idxDataSet]->SetTextAlign(12);
    legRat[idxDataSet]->SetTextSize(0.06);
    legRat[idxDataSet]->SetFillColor(1182);
    legRat[idxDataSet]->SetLineColor(0);
    legRat[idxDataSet]->SetBorderSize(0);
  }

  // -----------------------------------------------------

  TGraphErrors *etaGraphs[9][nDataSets][nMoments];

  for (int idxCent = 0; idxCent < nCent; idxCent++) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 4 ; idxMoment < nMoments; ++idxMoment) {
	float y[nEtaSets];
	float ey[nEtaSets];
	for (int idxEta = 0 ; idxEta < nEtaSets; ++idxEta) {	
	  y[idxEta]  = graphs[idxEta][idxDataSet][idxMoment]->GetY()[idxCent];
	  ey[idxEta] = graphs[idxEta][idxDataSet][idxMoment]->GetEY()[idxCent];
	} 
	
	etaGraphs[idxCent][idxDataSet][idxMoment] = new TGraphErrors(nEtaSets, aEtaSetBinWidth, y, 0, ey);
	PrepareGraph(etaGraphs[idxCent][idxDataSet][idxMoment]);
	
      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
    } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  } // for (int idxCent = 0; idxCent < 9; idxCent++) {

  // -----------------------------------------------------

  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    TPad *pad = SetupCanvas(Form("canDeltaEta_Ratio_%s", aDataSets[idxDataSet]), aDataSetsTitle[idxDataSet], "#Delta#eta", 0.45);
    
    for (int idxMoment = 4 ; idxMoment < nMoments; ++idxMoment) {
      pad->cd(idxMoment-3);
      
      for (int idxCent = 0; idxCent < 9; idxCent++) {
	if (idxCent != 0 && idxCent != 1 && idxCent != 8 && idxCent != 4 )
	  continue;
	
	TGraphErrors *g = etaGraphs[idxCent][idxDataSet][idxMoment];
	
	if (idxCent == 0) 
	  ShiftGraphX(g, -0.015);
	else if (idxCent == 1) 
	  ShiftGraphX(g, -0.005);
	else if (idxCent == 4) 
	  ShiftGraphX(g, 0.005);
	else if (idxCent == 8) 
	  ShiftGraphX(g, 0.015);
	
	ConfigGraph(g, idxMoment, idxCent);
	
	if (idxCent == 0) {
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
	    legRat[idxDataSet]->AddEntry(line1, "Poisson", "l");
	  }
	}
	
	g->Draw("PSAME");
	
	if (idxMoment == 4) 
	  legRat[idxDataSet]->AddEntry(etaGraphs[idxCent][idxDataSet][4], Form("%s", cent1[idxCent]), "pl");
	
      } // for (int idxCent = 0; idxCent < 9; idxCent++) { 
    } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
    
    pad->cd(2);
    TLatex *texb_3 = new TLatex(0.05, 0.55, "Au+Au collisions #sqrt{#it{s}_{NN}} = 14.5 GeV");
    texb_3->SetTextSize(0.07);
    texb_3->Draw("same");
    
    TLatex *texb_3a = new TLatex(0.05, 0.49, "Net-Charge, 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0, 0-5%");
    texb_3a->SetTextSize(0.07);
    texb_3a->Draw("same");
    
    TLatex *texb_3b = new TLatex(0.05,0.44, "statistical errors only");
    texb_3b->SetTextSize(0.06);
    texb_3b->Draw("same");
    
    pad->cd(1);
    TLatex *texb_4 = new TLatex(0.7, 19, "STAR Preliminary");
    texb_4->SetTextSize(0.07);
    texb_4->Draw("same");
    
    pad->cd(3);
    legRat[idxDataSet]->Draw();
    
    pad->Modified();
  } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  
  // -----------------------------------------------------

  SaveCanvas(name);
}

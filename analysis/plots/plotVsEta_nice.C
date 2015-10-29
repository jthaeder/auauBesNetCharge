
#include "include/toolsEtaNice.C"

const Char_t* aEtaSets[]         = {"2015-06-21_delta_0.1_0",
				    "2015-05-20_0.1",
				    "2015-06-03_delta_0.3_0", "2015-06-03_delta_0.3_1", "2015-06-03_delta_0.3_2", "2015-06-03_delta_0.3_3", 
				    "2015-06-03_delta_0.3_8",
				    "2015-06-03_delta_0.3_4", "2015-06-03_delta_0.3_5", "2015-06-03_delta_0.3_6", "2015-06-03_delta_0.3_7", 
				    "2015-05-20_0.2", 
				    "2015-06-03_delta_0.5_0", "2015-06-03_delta_0.5_1", "2015-06-03_delta_0.5_2", 
				    "2015-06-20_delta_0.5_6", 
				    "2015-06-03_delta_0.5_3", "2015-06-03_delta_0.5_4", "2015-06-03_delta_0.5_5", 
				    "2015-05-20_0.3", "2015-06-20_delta_0.7_0", "2015-05-20_0.4", "2015-06-20_delta_0.9_0", "2015-05-20_0.5"};  

const int     nEtaSets           = 24; 

const float   aEtaSetBinCenter[] = {0.,
				    0.,
				    -0.35, -0.25, -0.15, -0.05, 0., 0.05,  0.15,  0.25,  0.35, 
				    0,
				    -0.25, -0.15, -0.05, 0., 0.05, 0.15, 0.25,
				    0., 0., 0., 0., 0.};

const float   aEtaSetBinWidth[]  = {0.1,
				    0.2,
				    0., 0., 0., 0., 0.3, 0., 0., 0., 0.,
				    0.4,
				    0., 0., 0., 0.5, 0., 0., 0.,
				    0.6, 
				    0.7, 
				    0.8, 
				    0.9, 
				    1.0};

const int     nEtaSuperSets           = 10; 

const int     aEtaSuperSets[]         = { 1,   1,   9,   1,   7,   1,   1,   1,   1,   1}; 
const float   aEtaSuperSetsBinWidth[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; 

// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

}

// ______________________________________________________________________________________
void plotVsEta_nice(const Char_t* name = "ratioVsEta") {

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
    legRat[idxDataSet] = new TLegend(0.14, 0.12, 0.78, 0.27);
    legRat[idxDataSet]->SetTextAlign(12);
    legRat[idxDataSet]->SetTextSize(0.06);
    legRat[idxDataSet]->SetFillColor(1182);
    legRat[idxDataSet]->SetLineColor(0);
    legRat[idxDataSet]->SetBorderSize(0);
    legRat[idxDataSet]->SetNColumns(3);
  }

  // -----------------------------------------------------

  TGraphErrors *etaGraphs[9][nEtaSuperSets][nDataSets][nMoments];

  for (int idxCent = 0; idxCent < nCent; idxCent++) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      for (int idxMoment = 4 ; idxMoment < nMoments; ++idxMoment) {
	
	int idxEta = -1;
	for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {
	  float x[aEtaSuperSets[idxEtaSuper]];
	  float y[aEtaSuperSets[idxEtaSuper]];
	  float ex[aEtaSuperSets[idxEtaSuper]];
	  float ey[aEtaSuperSets[idxEtaSuper]];
	  
	  for (int idxEtaSub = 0 ; idxEtaSub < aEtaSuperSets[idxEtaSuper]; ++idxEtaSub) {
	    ++idxEta;
	    
	    x[idxEtaSub]  = aEtaSetBinCenter[idxEta];
	    ex[idxEtaSub] = aEtaSetBinWidth[idxEta]/2.;
	    y[idxEtaSub]  = graphs[idxEta][idxDataSet][idxMoment]->GetY()[idxCent];
	    ey[idxEtaSub] = graphs[idxEta][idxDataSet][idxMoment]->GetEY()[idxCent];
	  } // for (int idxEtaSub = 0 ; idxEtaSub < aEtaSuperSets[idxEtaSuper]; ++idxEtaSub) {
	  
	  etaGraphs[idxCent][idxEtaSuper][idxDataSet][idxMoment] = new TGraphErrors(aEtaSuperSets[idxEtaSuper], x, y, ex, ey);  
	  PrepareGraph(etaGraphs[idxCent][idxEtaSuper][idxDataSet][idxMoment]);
	} // for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {
	
      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
    } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  } // for (int idxCent = 0; idxCent < 9; idxCent++) {

  // -----------------------------------------------------

  for (int idxCent = 0; idxCent < 1; idxCent++) {
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      TPad *pad = SetupCanvas(Form("canEta_Ratio_%s_%s", aDataSets[idxDataSet], cent[idxCent]), aDataSetsTitle[idxDataSet], "#eta", 0.45);
    
      for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {
	
	if (idxEtaSuper == 1 || idxEtaSuper == 3 || idxEtaSuper == 5 || idxEtaSuper == 7 || idxEtaSuper == 8 )
	  continue;

	for (int idxMoment = 4 ; idxMoment < nMoments; ++idxMoment) {
	  pad->cd(idxMoment-3);
	  
	  TGraphErrors *g = etaGraphs[idxCent][idxEtaSuper][idxDataSet][idxMoment];

	  // if (idxCent == 0) 
	  //   ShiftGraphX(g, -0.015);
	  // else if (idxCent == 1) 
	  //   ShiftGraphX(g, -0.005);
	  // else if (idxCent == 4) 
	  //   ShiftGraphX(g, 0.005);
	  // else if (idxCent == 8) 
	  //   ShiftGraphX(g, 0.015);

	  ConfigGraph(g, idxMoment, idxEtaSuper);

	  if (idxEtaSuper == 0) {
	    g->Draw("AP");

	    TLine *line05 = new TLine( 0.5, aMinY[idxMoment],  0.5, aMaxY[idxMoment]);
	    line05->SetLineColor(kGray+1);
	    line05->SetLineStyle(3);
	    line05->Draw();
	    
	    TLine *line50 = new TLine(-0.5, aMinY[idxMoment], -0.5, aMaxY[idxMoment]);
	    line50->SetLineColor(kGray+1);
	    line50->SetLineStyle(3);
	    line50->Draw();
	    
	    TLine *line00 = new TLine(0., aMinY[idxMoment],0, aMaxY[idxMoment]);
	    line00->SetLineColor(kGray+1);
	    line00->SetLineStyle(3);
	    line00->Draw();

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
	    }
	  }
	  g->Draw("PSAME");

	} // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {

	if (idxCent == 0)  {
	  legRat[idxDataSet]->AddEntry(etaGraphs[0][idxEtaSuper][idxDataSet][4],    
	  			       Form("#Delta#eta = %.1f", aEtaSuperSetsBinWidth[idxEtaSuper]), "pl");
	}
      } // for (int idxEta = 0 ; idxEta < nEta; ++idxEta) {

      TLine *linex = new TLine(-0.55, 1, 0.55, 1);
      linex->SetLineColor(kGray+1);
      linex->SetLineStyle(2);
      linex->SetLineWidth(2);
      legRat[idxDataSet]->AddEntry(linex, "Poisson", "l");
      
      pad->cd(1);
      TLatex *texb_3 = new TLatex(-0.49, 9.3, "Au+Au collisions #sqrt{#it{s}_{NN}} = 14.5 GeV");
      texb_3->SetTextSize(0.07);
      texb_3->Draw("same");

      TLatex *texb_3a = new TLatex(-0.49,8.9, "Net-Charge, 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0, 0-5%");
      texb_3a->SetTextSize(0.07);
      texb_3a->Draw("same");

      pad->cd(2);
      TLatex *texb_3b = new TLatex(0.49,0.42, "statistical errors only");
      texb_3b->SetTextSize(0.06);
      texb_3b->SetTextAlign(31);
      texb_3b->Draw("same");

      TLatex *texb_3c = new TLatex(0.49,0.38, "horizontal error: #Delta#eta width");
      texb_3c->SetTextSize(0.06);
      texb_3c->SetTextAlign(31);
      texb_3c->Draw("same");

      TLatex *texb_4 = new TLatex(-0.49, 0.42, "STAR Preliminary");
      texb_4->SetTextSize(0.07);
      texb_4->Draw("same");

      pad->cd(3);
      legRat[idxDataSet]->Draw("lt");
      
      pad->Modified();
    } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
  } // for (int idxCent = 0; idxCent < 9; idxCent++) { 

  // -----------------------------------------------------

  SaveCanvas(name);
}

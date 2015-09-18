
#include "include/toolsEtaNice.C"

const Char_t* aEtaSets[]         = {"2015-06-03_delta_0.3_0", "2015-06-03_delta_0.3_1", "2015-06-03_delta_0.3_2", "2015-06-03_delta_0.3_3", 
				    "2015-06-03_delta_0.3_8",
				    "2015-06-03_delta_0.3_4", "2015-06-03_delta_0.3_5", "2015-06-03_delta_0.3_6", "2015-06-03_delta_0.3_7", 

				    "2015-06-03_delta_0.5_0", "2015-06-03_delta_0.5_1", "2015-06-03_delta_0.5_2", 
				    "2015-06-20_delta_0.5_6",
				    "2015-06-03_delta_0.5_3", "2015-06-03_delta_0.5_4", "2015-06-03_delta_0.5_5",
				    
				    "2015-06-21_delta_0.1_0",  "2015-05-20_0.1", "2015-06-03_delta_0.3_8", 
				    "2015-05-20_0.2", "2015-06-20_delta_0.5_6", "2015-05-20_0.3", 
				    "2015-06-20_delta_0.7_0", "2015-05-20_0.4", "2015-06-20_delta_0.9_0", "2015-05-20_0.5",

				    "2015-06-20_delta_0.5_6", "2015-06-21_delta_1.0_mult_0.5", "2015-05-20_0.5"};  

const int     nEtaSets           = 29; 

const float   aEtaSetBinCenter[] = {-0.35, -0.25, -0.15, -0.05, 0., 0.05,  0.15,  0.25,  0.35, 
 				    -0.25, -0.15, -0.05, 0, 0.05, 0.15, 0.25,
 				    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
				    0., 0., 0.};

const float   aEtaSetBinWidth[]  = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
 				    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
				    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
				    0.5, 1.0, 1.0};


const int     nEtaSuperSets           = 3; 
const int     aEtaSuperSets[]         = {  9,   7,  10,   3}; 
const float   aEtaSuperSetsBinWidth[] = {0.3, 0.5, 1.0, 1.0}; 
const Char_t* aEtaSuperSetsBinTitle[] = {"#Delta#eta = 0.3", "#Delta#eta = 0.5", "#Delta#eta = [0.1 - 1.0]", "#Delta#eta = 0.5/1.0, fixed mult"}; 

// ______________________________________________________________________________________
void SetGlobals() {
  // -- set globals

  Float_t aMinYtmp[7] = { -5, -5, -5, -2300, 3.5, -0.02, -12 };
  Float_t aMaxYtmp[7] = { 40, 250, 100, 900, 21,   0.62, 6};

  for (Int_t ii = 0; ii<7; ++ii) {
    aMinY[ii] = aMinYtmp[ii];
    aMaxY[ii] = aMaxYtmp[ii];
  }

  aMinX = 0;
  aMaxX = 375;
}

// ______________________________________________________________________________________
void plotEta_nice(const Char_t* name = "ratioEtaVsNpart") {

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

	TGraphErrors *g = graphs[idxEta][idxDataSet][idxMoment];
	PrepareGraph(graphs[idxEta][idxDataSet][idxMoment]);

	if (inFiles[idxEta][idxDataSet][idxMoment])
	  (inFiles[idxEta][idxDataSet][idxMoment])->Close();
      }
  
  // -----------------------------------------------------

  TLegend *legRat[nEtaSuperSets][nDataSets];

  for (int idxEtaSuper = 2 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   
    for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
      legRat[idxEtaSuper][idxDataSet] = new TLegend(0.12, 0.12, 0.70, 0.57);
      legRat[idxEtaSuper][idxDataSet]->SetTextAlign(12);
      legRat[idxEtaSuper][idxDataSet]->SetTextSize(0.06);
      legRat[idxEtaSuper][idxDataSet]->SetTextFont(42);
      legRat[idxEtaSuper][idxDataSet]->SetFillColor(0);
      legRat[idxEtaSuper][idxDataSet]->SetLineColor(0);
      legRat[idxEtaSuper][idxDataSet]->SetBorderSize(0);
    }
  }

  // -----------------------------------------------------
  
  for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {
    int idxEta = -1;
    for (int idxEtaSuper = 2 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   
      
      TPad *pad = SetupCanvas(Form("canEta_Ratio_%s_%d", aDataSets[idxDataSet], idxEtaSuper), aDataSetsTitle[idxDataSet], "<N_{part}>", 0.45);

      for (int idxMoment = 4 ; idxMoment < nMoments; ++idxMoment) {
	pad->cd(idxMoment-3);
	
	idxEta = aEtaSuperSets[0]+aEtaSuperSets[1]-1;
	for (int idxEtaSub = 0 ; idxEtaSub < aEtaSuperSets[idxEtaSuper]; ++idxEtaSub) {
	  ++idxEta;
	  
	  if (idxEtaSub == 1 || idxEtaSub == 3 || idxEtaSub == 5 || idxEtaSub == 7 || idxEtaSub == 8)
	    continue;

	  TGraphErrors *g = graphs[idxEta][idxDataSet][idxMoment];
	 
	  if (idxEtaSub == 0) 
	    ShiftGraphX(g, -10);
	  else if (idxEtaSub == 2) 
	    ShiftGraphX(g, -5);
	  else if (idxEtaSub == 4) 
	    ShiftGraphX(g, 0);
	  else if (idxEtaSub == 6) 
	    ShiftGraphX(g, 5 );
	  else if (idxEtaSub == 9) 
	    ShiftGraphX(g, 10 );

	  ConfigGraph(g, idxMoment, idxEtaSub);

	  if (idxEtaSub == 0) { 
	    g->Draw("AP");
	    
	    if (idxMoment == 5) {
	      TLine *line0 = new TLine(aMinX, 0, aMaxX, 0);
	      line0->SetLineColor(kGray+1);
	      line0->SetLineStyle(2);
	      line0->SetLineWidth(2);
	      line0->Draw();
	    }
	    else if (idxMoment == 6) {
	      TLine *line1 = new TLine(aMinX, 1, aMaxX, 1);
	      line1->SetLineColor(kGray+1);
	      line1->SetLineStyle(2);
	      line1->SetLineWidth(2);
	      line1->Draw();
	      legRat[idxEtaSuper][idxDataSet]->AddEntry(line1, "Poisson", "l");
	    }
	  }

	  g->Draw("PSAME");

	  if (idxMoment == 4) 
	    legRat[idxEtaSuper][idxDataSet]->AddEntry(graphs[idxEta][idxDataSet][4], Form("#Delta#eta = %.1f", aEtaSetBinWidth[idxEta]), "pl");
						      
	} // for (int idxEtaSub = 0 ; idxEtaSub < aEtaSuperSets[idxEtaSuper]; ++idxEtaSub) { 
      } // for (int idxMoment = 0 ; idxMoment < nMoments; ++idxMoment) {
      
      pad->cd(1);
      TLatex *texb_3 = new TLatex(130,19, "Au+Au collisions #sqrt{#it{s}_{NN}} = 14.5 GeV");
      texb_3->SetTextSize(0.07);
      texb_3->Draw("same");

      TLatex *texb_3a = new TLatex(130,17.5, "Net-Charge, 0.2 < #it{p}_{T} (GeV/#it{c}) < 2.0");
      texb_3a->SetTextSize(0.07);
      texb_3a->Draw("same");

      TLatex *texb_3b = new TLatex(130,16, "statistical errors only");
      texb_3b->SetTextSize(0.06);
      texb_3b->SetTextFont(42);
      texb_3b->Draw("same");

      pad->cd(2);
      TLatex *texb_4 = new TLatex(15, 0.55, "STAR Preliminary");
      texb_4->SetTextSize(0.07);
      texb_4->Draw("same");

      pad->cd(3);
      legRat[idxEtaSuper][idxDataSet]->Draw("lt");
      
      pad->Modified();

    } // for (int idxEtaSuper = 0 ; idxEtaSuper < nEtaSuperSets; ++idxEtaSuper) {   
  } // for (int idxDataSet = 0 ; idxDataSet < nDataSets; ++idxDataSet) {

  // -----------------------------------------------------

  SaveCanvas(name);
}

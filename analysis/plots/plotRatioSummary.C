#include "include/toolsEnergyNice.C"

const Int_t   nMom       = 3;
const Char_t *aMom[3]    = {"VM", "SDSk", "KV"};
const Char_t *aMom2[3]   = {"VM", "SD", "KV"};

const Char_t *aTitle[3]  = {"#sigma^{2}/M", "S#sigma/Skellam", "#kappa#sigma^{2}"};
const Char_t *aTitle2[3] = {"#sigma^{2}/M", "S#sigma", "#kappa#sigma^{2}"};

// ______________________________________________________________________________________
void SetGlobals(Int_t idx) {
  // -- set globals

  // Kaon
  Float_t aMinYtmp[3][3] = { { -10, -9, -19}, {-2, -0.3, -3.5}, { -0.5, 0.5, -0.5} };
  Float_t aMaxYtmp[3][3] = { {160, 11,  12},  {45, 2.3, 2.8},   {12, 1.18, 4.5} };

  for (Int_t ii = 0; ii<3; ++ii) {
    aMinY[ii] = aMinYtmp[idx][ii];
    aMaxY[ii] = aMaxYtmp[idx][ii];
  }

  // aMarkers[2] = 30;
  // aMarkers[1] = 24;
  // aMarkers[0] = 25;

  // aColors[1] = kAzure;
  // aColors[0] = kBlack;
}

// ______________________________________________________________________________________
void plotRatioSummary(const Char_t* name = "ratioSummary") {

  gROOT->LoadMacro("include/toolsEnergyNice.C++");

  SetupStyle();

  // -----------------------------------------------------
  
  TFile *file = TFile::Open("STAR_Preliminary.root");

  // -----------------------------------------------------

  TGraphErrors *graphStat[nNames][nMom][nCent];
  TGraphErrors *graphSys[nNames][nMom][nCent];
  TGraphErrors *graphPoisson[nNames][nMom][nCent];
  TGraphErrors *graphUrqmd[nNames][nMom][nCent];
  TGraphErrors *graph14[nNames][nMom][nCent];

  // -----------------------------------------------------

  for (int idxN = 0; idxN < nNames; idxN++) {
    TList *list = static_cast<TList*>(file->Get(aNames[idxN]));

    for (int idxMoment = 0; idxMoment < nMom; ++idxMoment) {
      for (int idxCent = 0; idxCent < nCent; ++idxCent) {
	graphStat[idxN][idxMoment][idxCent]      = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_stat",    aNames[idxN], aMom[idxMoment], cent[idxCent])));
      	graphSys[idxN][idxMoment][idxCent]       = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_sys",     aNames[idxN], aMom[idxMoment], cent[idxCent])));
	graphPoisson[idxN][idxMoment][idxCent]   = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_poisson", aNames[idxN], aMom[idxMoment], cent[idxCent])));
	graphUrqmd[idxN][idxMoment][idxCent]     = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_urqmd",   aNames[idxN], aMom[idxMoment], cent[idxCent])));
	graph14[idxN][idxMoment][idxCent]        = static_cast<TGraphErrors*>(list->FindObject(Form("%s_%s_sNN_%s_14",   aNames[idxN], aMom[idxMoment], cent[idxCent])));
      }
    }
  }

  // -----------------------------------------------------
  
  TLegend *leg[nMom];

  for (Int_t idxMoment = 0; idxMoment < nMom; ++idxMoment) {
    
    canA.Add(new TCanvas(Form("can%s_%s", name, aMom[idxMoment]), Form("CumulantRatio_%s", aTitle[idxMoment]), 0, 0 , 1000, 700));
    can = static_cast<TCanvas*>(canA.Last());
    can->SetFillColor(1182);
    can->SetBorderMode(0);
    can->SetBorderSize(0.0);
    can->SetFrameFillColor(1182);
    can->SetFrameBorderMode(0);
    can->cd();
  
    pad = new TPad("pad", "pad", 0.00, 0.0, 1., 1.);
    pad->SetTopMargin(0.);
    pad->SetBottomMargin(0.);
    pad->SetLeftMargin(0.);
    pad->SetRightMargin(0.);
    pad->SetBorderMode(0);
    pad->SetFillColor(1182);
    pad->Draw();
    pad->cd();
    
    TPad *padSet0[2];
    TPad *padSet1[3];
    TPad *padSet2[2][2];
    
    Double_t legTextSize = 0.06;
    Double_t legHeight = (idxMoment == 0) ? (legTextSize+0.02)*5 : (legTextSize+0.02)*4;
    leg[idxMoment] = new TLegend(0.15, 0.65-legHeight, 0.8, 0.65);
    leg[idxMoment]->SetTextAlign(12);
    leg[idxMoment]->SetTextSize(legTextSize);
    leg[idxMoment]->SetTextFont(42);
    leg[idxMoment]->SetFillColor(1182);
    leg[idxMoment]->SetLineColor(0);
    leg[idxMoment]->SetBorderSize(0);
    
    for (Int_t idx = 0; idx < 2; ++idx) {
      pad->cd();
      padSet0[idx] = new TPad(Form("pad_%d", idx), "", idx*0.5, 0.0, (idx+1)*0.5, 1.);
      padSet0[idx]->SetTopMargin(0.);
      padSet0[idx]->SetBorderMode(0);
      padSet0[idx]->SetRightMargin(0.);
      padSet0[idx]->SetLeftMargin(0.);
      padSet0[idx]->SetFillColor(1182);
      padSet0[idx]->Draw();
    }

    Double_t marginRight =  0.01;
    Double_t marginLeft  =  0.13;

    // left column
    Double_t heightBox = 0.5;
    Double_t margin    = 0.15; // percent

    Double_t marginHist = heightBox*margin;
 
    // right column
    Double_t marginNew    = marginHist / (2*heightBox);

    padSet0[0]->cd();
      
    // upper left
    padSet1[0] = new TPad(Form("pad_0_1"), "", 0.0, 1.0-heightBox, 1.0, 1.0);
    padSet1[0]->SetTopMargin(margin);
    padSet1[0]->SetBottomMargin(0);
    padSet1[0]->SetRightMargin(marginRight);
    padSet1[0]->SetLeftMargin(marginLeft);
    padSet1[0]->SetBorderMode(0);
    padSet1[0]->SetFillColor(1182);
    padSet1[0]->Draw();
    
    // lower left
    padSet1[1] = new TPad(Form("pad_0_0"), "", 0.0, 0., 1.0, heightBox);
    padSet1[1]->SetTopMargin(0.);
    padSet1[1]->SetBorderMode(0);
    padSet1[1]->SetBottomMargin(margin);
    padSet1[1]->SetRightMargin(marginRight);
    padSet1[1]->SetLeftMargin(marginLeft);
    padSet1[1]->SetFillColor(1182);
    padSet1[1]->Draw();
    
    padSet0[1]->cd();

    padSet1[2] = new TPad(Form("pad_1_1"), "", 0.0, 1.0-marginNew-heightBox, 1.0, 1.0-marginNew);
    padSet1[2]->SetTopMargin(0.01);
    padSet1[2]->SetBottomMargin(margin);
    padSet1[2]->SetRightMargin(marginRight);
    padSet1[2]->SetLeftMargin(marginLeft);
    padSet1[2]->SetBorderMode(0);
    padSet1[2]->SetFillColor(1182);
    padSet1[2]->Draw();

    for (int idxN = 0; idxN < nNames; idxN++) {
      padSet1[idxN]->cd();

      SetGlobals(idxN);      

      gPad->SetLogx();
      for (int idxCent = 0; idxCent < nCent; ++idxCent) {
       	if (idxCent != 0 && idxCent != 8)
	  continue;
	
	TGraphErrors* gStat    = graphStat[idxN][idxMoment][idxCent];
	TGraphErrors* gSys     = graphSys[idxN][idxMoment][idxCent];
	TGraphErrors* gUrqmd   = graphUrqmd[idxN][idxMoment][idxCent];
	TGraphErrors* gPoisson = graphPoisson[idxN][idxMoment][idxCent];

	PrepareGraphCumulants(gStat);
	PrepareGraphCumulants(gSys);
	PrepareGraphCumulants(gUrqmd);
	PrepareGraphCumulants(gPoisson);

        ConfigGraph(gStat,    idxMoment, idxCent);
        ConfigGraph(gSys,     idxMoment, idxCent, 2);
	ConfigGraph(gUrqmd,   idxMoment, idxCent, 3);
	ConfigGraph(gPoisson, idxMoment, idxCent, 1);

	gStat->GetXaxis()->SetNdivisions(2, 1, 0);

	// -- draw box
	if (idxCent == 0) 
	  gStat->Draw("AZP");

	// -- draw urqmd
	if (idxCent == 0) {
	  if (idxMoment == 0)
	    gUrqmd->Draw("L,same");
	  gUrqmd->Draw("E3,same");
	}

	// -- draw poisson 
	if (idxMoment > 0) {
	  if (idxCent == 0) {
	    TLine *line1 = new TLine(aMinX, 1, aMaxX, 1);
	    line1->SetLineColor(gPoisson->GetLineColor());
	    line1->SetLineStyle(2);
	    line1->SetLineWidth(2);
	    line1->Draw();
	  }
	}
	else if (idxMoment == 0)
	  gPoisson->Draw("L,SAME");
	
	// -- draw datapoints
	//gSys->Draw("B2,SAME");
	gStat->Draw("ZP,SAME");
	gSys->Draw("[],SAME");

      } // for (int idxCent = 0; idxCent < nCent; ++idxCent) {
    } // for (int idxN = 0; idxN < nNames; idxN++) {
    
    // -----------------------------------------------------  
    
    leg[idxMoment]->AddEntry(graphStat[0][idxMoment][0],      Form("%s STAR %s",    cent1[0], aTitle[idxMoment]), "pl");
    leg[idxMoment]->AddEntry(graphStat[0][idxMoment][8],      Form("%s STAR %s",    cent1[8], aTitle[idxMoment]), "pl");
    leg[idxMoment]->AddEntry(graphUrqmd[0][idxMoment][0],     Form("%s UrQMD %s",   cent1[0], aTitle[idxMoment]), "f");
    leg[idxMoment]->AddEntry(graphPoisson[0][idxMoment][0],   Form("%s Poisson %s", cent1[0], aTitle[idxMoment]), "l");
    if (idxMoment == 0)
      leg[idxMoment]->AddEntry(graphPoisson[0][idxMoment][8], Form("%s Poisson %s", cent1[8], aTitle[idxMoment]), "l");
    
    // -----------------------------------------------------

    for (Int_t idx1 = 0; idx1 < 2; ++idx1) {
      for (Int_t idx2 = 0; idx2 < 2; ++idx2) {
	pad->cd();
	padSet2[idx1][idx2] = new TPad(Form("padCover_%d_%d", idx1, idx2), "", idx1*0.5, idx2*0.5, (idx1+1)*0.5, (idx2+1)*0.5);
	padSet2[idx1][idx2]->SetTopMargin(0.);
	padSet2[idx1][idx2]->SetBorderMode(0);
	padSet2[idx1][idx2]->SetRightMargin(0.);
	padSet2[idx1][idx2]->SetLeftMargin(0.);
	padSet2[idx1][idx2]->SetFillColor(1182);
	padSet2[idx1][idx2]->Draw();
	padSet2[idx1][idx2]->cd();
      }
    }
    
    // [idxMoment][idxN]
    Double_t aTextPosX[3][3] = { {marginLeft+0.03, marginLeft+0.03, marginLeft+0.03},  // VM
				 {marginLeft+0.03, marginLeft+0.03, marginLeft+0.03},  // SK
				 {0.45, 0.45, 0.45} };  // KV

    Double_t aTextPosY[3][3] = { {0.76, 0.76, 0.76+margin},
				 {0.76, 0.76, 0.76+margin},
				 {0.20, 0.7 , 0.20+margin} };

    Double_t textSizeName = 0.07;

    TLatex *texb_6[3];
    TLatex *texb_7[3];
    TLatex *texb_8[3];

    for (Int_t ii = 0; ii < 3 ; ++ ii) {
      texb_6[ii] = new TLatex(aTextPosX[idxMoment][ii], aTextPosY[idxMoment][ii], aNames[ii]);
      texb_6[ii]->SetTextSize(textSizeName);

      texb_7[ii]= new TLatex(aTextPosX[idxMoment][ii], aTextPosY[idxMoment][ii]-textSizeName, aNamesPt[ii]);
      texb_7[ii]->SetTextSize(0.055);

      if (ii == 2)
	texb_8[ii] = new TLatex(0.055, ((1.0-margin)/2.)+margin, Form("%s %s", aNames[ii], aTitle[idxMoment]));
      else
	texb_8[ii] = new TLatex(0.055, ((1.0-margin)/2.), Form("%s %s", aNames[ii], aTitle[idxMoment]));
      texb_8[ii]->SetTextSize(0.08);
      texb_8[ii]->SetTextAngle(90);
      texb_8[ii]->SetTextAlign(21);
    }

    // -- proton
    Int_t idxN = 2;
    padSet2[0][0]->cd();
    
    texb_6[idxN]->Draw("same");
    texb_7[idxN]->Draw("same");
    texb_8[idxN]->Draw("same");

    TLatex *texb_5c = new TLatex(0.45, 0.03, "#sqrt{#it{s}_{NN}} (GeV)");
    texb_5c->SetTextSize(0.06);
    texb_5c->Draw("same");

    // -- charge
    idxN = 0;
    padSet2[0][1]->cd();

    texb_6[idxN]->Draw("same");
    texb_7[idxN]->Draw("same");
    texb_8[idxN]->Draw("same");

    // -- kaon
    idxN = 1;
    padSet2[1][1]->cd();

    texb_6[idxN]->Draw("same");
    texb_7[idxN]->Draw("same");
    texb_8[idxN]->Draw("same");

    // -- text
    padSet2[1][0]->cd();
    TLatex *texb_5d = new TLatex(0.45, 0.88, "#sqrt{#it{s}_{NN}} (GeV)");
    texb_5d->SetTextSize(0.06);
    texb_5d->Draw("same");

    TLatex *texb_3 = new TLatex(0.15, 0.70, Form("BES I - Net-Particle %s", aTitle2[idxMoment]));
    texb_3->SetTextSize(0.1);
    texb_3->Draw("same");
    
    TLatex *texb_4 = new TLatex(0.65, 0.05, "STAR Preliminary");
    texb_4->SetTextSize(0.065);
    texb_4->Draw("same");

    leg[idxMoment]->Draw("lt");

    // -----------------------------------------------------
    
  } // for (Int_t idxMoment = 0; idxMoment < nMom; ++idxMoment) {
  
  SaveCanvas(name);
}

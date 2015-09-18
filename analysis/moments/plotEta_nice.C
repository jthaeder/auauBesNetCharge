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

const Char_t* aMomentsTitle[]  = {"c_{1}","c_{2}","c_{3}","c_{4}","#sigma^{2}/M","S #sigma", "#kappa #sigma^{2}"};

const Char_t* aMoments[]       = {"C1","C2","C3","C4","VM","SD","KV"};
const int     nMoments         = 7;

const Char_t* aDataSetsTitle[] = {"corrected - 11.5 GeV #epsilon_{1} , #epsilon_{2}"};
const Char_t* aDataSets[]      = {"twoeff_11"};
const int     nDataSets        = 1;

const int   nCent    = 1;
const char *cent[9]  = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};
const char *cent1[9] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"};

const Color_t aColorsRB[] = {kBlack, kOrange+9, kOrange+9, kYellow, kAzure, kAzure, kCyan+2, kBlue+1, kRed+2, kRed+2};

const int  aMarkers[] = {24, 32, 25, 25, 26, 28, 27, 30, 32, 30};

const float aMinY[7] = { -5, -5, -5, -2300, 3.5, -0.02, -12 };
const float aMaxY[7] = { 40, 250, 100, 900, 21,   0.62, 6};

const float aMinX = 0;
const float aMaxX = 375;

TObjArray canA;
// ______________________________________________________________________________________
void PrepareGraph(TGraphErrors* g) {
  // -- Prepare Graph
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");
  g->GetXaxis()->SetNdivisions(6, 5, 0);
  g->GetYaxis()->SetNdivisions(6, 5, 0);
  g->GetXaxis()->SetLabelSize(0.07);
  g->GetYaxis()->SetLabelSize(0.07);
}

// ______________________________________________________________________________________
void ConfigGraph(TGraphErrors* g, Int_t yLimitIdx, Int_t idxMarker) {
  // -- Config Graph
  g->SetMinimum(aMinY[yLimitIdx]);
  g->SetMaximum(aMaxY[yLimitIdx]);
  g->GetXaxis()->SetLimits(aMinX, aMaxX);  
  g->SetMarkerSize(1.4);
  g->SetMarkerStyle(aMarkers[idxMarker]);
  g->SetMarkerColor(aColorsRB[idxMarker]);
  g->SetLineColor(aColorsRB[idxMarker]);
}

// ______________________________________________________________________________________
void ShiftGraphX(TGraphErrors* g, Double_t shift) {
  // -- shift datapoints

  for (Int_t idx = 0 ; idx < g->GetN(); ++idx) {
    Double_t x, y;
    g->GetPoint(idx, x, y);
    g->SetPoint(idx, x+shift, y);
  }

  PrepareGraph(g);
}

// ______________________________________________________________________________________
TPad* SetupCanvas(const Char_t* canName, const Char_t *canTitle, const Char_t *xTitle, Float_t xPosTitle) {
  // -- setup canvas and pad
  
  canA.Add(new TCanvas(canName, canTitle, 1200, 0 , 600, 1000));
  TCanvas *can = static_cast<TCanvas*>(canA.Last());
  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetBorderSize(0.0);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->cd();
  
  TPad* pad = new TPad("pad", "pad",0.05, 0.06, 0.99, 0.98);
  pad->SetBorderMode(0);
  pad->SetFillColor(0);
  pad->Draw();
  pad->cd();
  pad->Divide(1, 3, 0., 0., 0);

  can->cd();

  TLatex *texb_5 = new TLatex(xPosTitle, 0.03, xTitle);
  texb_5->SetTextSize(0.04);
  texb_5->Draw("same");
  
  TLatex *texb_6a = new TLatex(0.05,0.8, aMomentsTitle[4]);
  texb_6a->SetTextSize(0.04);
  texb_6a->SetTextAngle(90);
  texb_6a->Draw("same");
  
  TLatex *texb_6b = new TLatex(0.05,0.5, aMomentsTitle[5]);
  texb_6b->SetTextSize(0.04);
  texb_6b->SetTextAngle(90);
  texb_6b->Draw("same");
  
  TLatex *texb_6c = new TLatex(0.05,0.2, aMomentsTitle[6]);
  texb_6c->SetTextSize(0.04);
  texb_6c->SetTextAngle(90);
  texb_6c->Draw("same");
  
  return pad;
}

// ______________________________________________________________________________________
void SaveCanvas(const Char_t* name) {
  // -- Write out Canvas

   gSystem->Exec(Form("mkdir -p results/%s/png",  name));
   gSystem->Exec(Form("mkdir -p results/%s/pdf",  name));
   gSystem->Exec(Form("mkdir -p results/%s/eps",  name));
   gSystem->Exec(Form("mkdir -p results/%s/gif",  name));
   gSystem->Exec(Form("mkdir -p results/%s/root", name));

   // -----------------------------------------------------

   for (Int_t idx = 0; idx < canA.GetEntriesFast() ; ++idx) {
     TCanvas* c = static_cast<TCanvas*>(canA.At(idx));
     if (!c)
       continue;

     c->SaveAs(Form("results/%s/png/%s_14GeV.png",   name, c->GetName()));
     c->SaveAs(Form("results/%s/eps/%s_14GeV.eps",   name, c->GetName()));
     c->SaveAs(Form("results/%s/gif/%s_14GeV.gif",   name, c->GetName()));
     c->SaveAs(Form("results/%s/pdf/%s_14GeV.pdf",   name, c->GetName()));
     c->SaveAs(Form("results/%s/root/%s_14GeV.C",    name, c->GetName()));
     c->SaveAs(Form("results/%s/root/%s_14GeV.root", name, c->GetName()));
   }
 }
// ______________________________________________________________________________________
void SetupStyle() {
  // -- Setup Style

  gROOT->SetStyle("Plain");
  
  gStyle->SetHatchesSpacing(0.8);
  gStyle->SetHatchesLineWidth(1);
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFillStyle(1001);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);

  gStyle->SetLegendBorderSize(0);

  Int_t font = 42;
  
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  
  gStyle->SetTickLength(0.02,"xy");
  gStyle->SetEndErrorSize(3);
  
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.009,"xyz");

  gStyle->SetTitleSize(0.06);  
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.12,"x");   // JMT 1.15	
  gStyle->SetTitleOffset(1.13,"yz");   // JMT 1.15	
  gStyle->SetTitleSize(0.045,"xyz");  


  gStyle->SetMarkerSize(1.6);  // JMT 1.1
  gStyle->SetPalette(1,0); 

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetLineWidth(1);

  // -- Set plot styles for 2D colz

  const Int_t nRGBs = 5;
  const Int_t nCont = 255;
  
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);

  TColor *color = new TColor(1182, 1, 0, 0, " ", 0);
}

// ______________________________________________________________________________________
void plotEta_nice(const Char_t* name = "etaVsNpart_nice") {

  SetupStyle();

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

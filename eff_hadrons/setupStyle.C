//-*- Mode: C++ -*-

// Make Efficiciency and Contamination for NetParticle Studies
// - Setup Style Macro -
// Author: Jochen Thaeder <jochen@thaeder.de> 

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TColor.h"

// ____________________________________________________________________________
void setupStyle() {
  // Setup Style

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
  
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");

  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.3,"xyz");   // JMT 1.15	
  gStyle->SetTitleSize(0.04,"xyz");  
  gStyle->SetTitleSize(0.04);  

  gStyle->SetMarkerSize(1.2);  // JMT 1.1
  gStyle->SetPalette(1,0); 

  gStyle->SetOptDate(20);

  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("nemrks");
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

}

void DivideCanvas(TCanvas *can, Int_t nCols, Int_t nRows) {

  Int_t idxPad = 0;
  can->cd();

  Float_t widthRows = 1./nRows;
  Float_t widthCols = 1./nCols;

  for (Int_t idxRows = nRows; idxRows > 0; --idxRows) {
    Float_t idxY0 = widthRows*(idxRows-1);
    Float_t idxY1 = widthRows*idxRows;
    
    for (Int_t idxCols = nCols; idxCols > 0; --idxCols) {
      Float_t idxX0 = widthCols*(idxCols-1);
      Float_t idxX1 = widthCols*idxCols;
      
      ++idxPad;
      TPad *pad = new TPad(Form("%s_%d",can->GetName(), idxPad), "", idxX0, idxY0, idxX1, idxY1);
      pad->SetFillColor(idxPad);
      pad->Draw();
    }
  }
}

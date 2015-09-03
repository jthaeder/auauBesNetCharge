#include "TCanvas.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "math.h"
#include <iostream>

using namespace std;

// _______________________________________________________________
TH1D* getEffRatio(TH1D* hMC, TH1D* hRec, const char*name) {
  TH1D* hEff = static_cast<TH1D*>(hMC->Clone(name));

  for (int idx = 1; idx <= hEff->GetXaxis()->GetNbins(); idx++) {
    double binContent = hEff->GetBinContent(idx);
    double ratio = (binContent < 0.00001) ? 0. : hRec->GetBinContent(idx)/(double)binContent;
    double error = (binContent < 0.00001) ? 0. : sqrt((double)binContent*ratio*(1.0-ratio))/(double)binContent;
    hEff->SetBinContent(idx, ratio);
    hEff->SetBinError(idx, error);
  }
 
  return hEff;
}

// _______________________________________________________________
TH2D* getEffRatio(TH2D* hMC, TH2D* hRec, const char*name) {
  TH2D* hEff = static_cast<TH2D*>(hMC->Clone(name));

  for (int idxX = 1; idxX <= hEff->GetXaxis()->GetNbins(); idxX++) {
    for (int idxY = 1; idxY <= hEff->GetYaxis()->GetNbins(); idxY++) {
      double binContent = hEff->GetBinContent(idxX, idxY);
      double ratio = (binContent < 0.00001) ? 0. : hRec->GetBinContent(idxX, idxY)/(double)binContent;
      double ratio = (binContent < 0.00001) ? 0. : sqrt((double)binContent*ratio*(1.0-ratio))/(double)binContent);
      hEff->SetBinContent(idxX, idxY, ratio);
      hEff->SetBinError(idxX, idxY, error);
    }
  }
 
  return hEff;
}

// _______________________________________________________________
void make_efficiencies(const char* particle, int energy = 14) {
  static const Double_t pi = TMath::Pi();

  gStyle->SetOptStat(111111);
  gStyle->SetPalette(1);
  //gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameBorderMode(0);
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadRightMargin(0.05);
  //gStyle->SetPadBottomMargin(0.05);
  gStyle->SetPadLeftMargin(0.11);

  // --------------------------------------------------------------------------------
  // -- track cuts
  // --------------------------------------------------------------------------------
  
  float min_pt           = 0.2;
  float max_pt           = 2.0;

  float min_dedx         = 5;  // default - was 10 , try not to use  10
  float min_nCommon      = 10;
  float min_fitpts       = 20;//15 // 10
  float min_fitpts_nposs = .52;
  float max_eta          = 0.5;
  float max_dca          = 1.0;
  
  // --------------------------------------------------------------------------------
  // -- vertex cuts
  // --------------------------------------------------------------------------------

  // -- Vz cut
  float max_z = 30;
  if (     fabs(energy-200) < 0.01) max_z = 30;
  else if (fabs(energy-62)  < 0.01) max_z = 30;
  else if (fabs(energy-39)  < 0.01) max_z = 30;
  else if (fabs(energy-27)  < 0.01) max_z = 30;
  else if (fabs(energy-19)  < 0.01) max_z = 30;
  else if (fabs(energy-14)  < 0.01) max_z = 30;
  else if (fabs(energy-11)  < 0.01) max_z = 30;
  else if (fabs(energy-7)   < 0.01) max_z = 30;
  else {
    cout << "bad energy" << endl;
    return;
  }

  // -- Vr cut 
  float vyo, vxo, max_r;
  if(fabs(energy-14)<0.01){
    vyo = -0.89;
    vxo = 0;
    max_r = 1;
  }
  else {
    vyo = 0;
    vxo = 0;
    max_r = 2;
  }

  // --------------------------------------------------------------------------------
  TString NRG(Form("%d", energy));
  
  cout << "particle         = " << particle << endl;
  cout << "energy           = " << NRG << endl;
  cout << "vxo              = " << vxo << endl;
  cout << "vyo              = " << vyo << endl;
  cout << "max_r            = " << max_r << endl;
  cout << "max_z            = " << max_z << endl;
  cout << "min_dedx         = " << min_dedx << endl;
  cout << "min_fitpts       = " << min_fitpts << endl;
  cout << "min_fitpts_nposs = " << min_fitpts_nposs << endl;
  cout << "max_eta          = " << max_eta << endl;
  cout << "max_dca          = " << max_dca << endl;
  cout << "min_pt           = " << min_pt << endl;
  cout << "max_pt           = " << max_pt << endl;

  // --------------------------------------------------------------------------------
  // -- get file names
  // --------------------------------------------------------------------------------
  TString out_file, in_file;
  TString prefix("embeddingTrees");
  TString postfix("efficiency");

  gSystem->Exec("mkdir -p efficiency");

  int PID = 0;
  
  if(TString(particle) == TString("piplus")){
    in_file = prefix + TString("/SinglePiPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piplus") + NRG + TString("GeV.root");
    //    mass = .13975;
    PID = 8;
  }
  else if(TString(particle) == TString("piminus")){
    in_file = prefix + TString("/SinglePiMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piminus") + NRG + TString("GeV.root");
    //    mass = .13975;
    PID = 9;
  }
  else if(TString(particle) == TString("kaonplus")){
    in_file = prefix + TString("/SingleKPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kplus") + NRG + TString("GeV.root");
    //    mass = .493677;
    PID = 11;
  }
  else if(TString(particle) == TString("kaonminus")){
    in_file = prefix + TString("/SingleKMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kminus") + NRG + TString("GeV.root");
    //    mass = .493677;
    PID = 12;
  }
  else if(TString(particle) == TString("protonplus")){
    in_file = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pplus") + NRG + TString("GeV.root");
    //    mass = .93827;
    PID = 14;
  }
  else if(TString(particle) == TString("protonminus")){
    in_file = prefix + TString("/SinglePbarNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pminus") + NRG + TString("GeV.root");
    //    mass = .93827;
    PID = 15;
  }

  // --------------------------------------------------------------------------------
  // -- get input file
  // --------------------------------------------------------------------------------
  cout << "opening input file: " << in_file << endl;
  
  TFile *fIn = TFile::Open(in_file, "READ");
  TNtuple *MatchedPairs = static_cast<TNtuple*>(fIn->Get("MatchedPairs_NT"));
  TNtuple *McTrack      = static_cast<TNtuple*>(fIn->Get("McTrack_NT"));

  // --------------------------------------------------------------------------------
  // -- get matched ntuple
  float Dedx, RefMult, RefMultCorrected, CentralityWeight, Centrality16, VertexX, VertexY, VertexZ;
  float PtMc, PzMc, EtaMc, PhiMc, PtPr, PtGl, EtaPr, PhiPr, DcaGl, DcaZGl, DcaXYGl, Flag, FitPts, DedxPts;
  float AllPts, NPossible, ParentGeantId, GeantId,mErrP, NCommonHit;
  MatchedPairs->SetBranchAddress("Dedx", &Dedx);
  MatchedPairs->SetBranchAddress("RefMult", &RefMult);
  MatchedPairs->SetBranchAddress("RefMultCorrected", &RefMultCorrected);
  MatchedPairs->SetBranchAddress("CentralityWeight", &CentralityWeight);
  MatchedPairs->SetBranchAddress("Centrality16", &Centrality16);
  MatchedPairs->SetBranchAddress("VertexX", &VertexX);
  MatchedPairs->SetBranchAddress("VertexY", &VertexY);
  MatchedPairs->SetBranchAddress("VertexZ", &VertexZ);
  MatchedPairs->SetBranchAddress("PtMc", &PtMc);
  MatchedPairs->SetBranchAddress("PzMc", &PzMc);
  MatchedPairs->SetBranchAddress("EtaMc", &EtaMc);
  MatchedPairs->SetBranchAddress("PhiMc", &PhiMc);
  MatchedPairs->SetBranchAddress("PtPr", &PtPr);
  MatchedPairs->SetBranchAddress("mErrP", &mErrP);
  MatchedPairs->SetBranchAddress("PtGl", &PtGl);
  MatchedPairs->SetBranchAddress("EtaPr", &EtaPr);
  MatchedPairs->SetBranchAddress("PhiPr", &PhiPr);
  MatchedPairs->SetBranchAddress("DcaGl", &DcaGl);
  MatchedPairs->SetBranchAddress("DcaZGl", &DcaZGl);
  MatchedPairs->SetBranchAddress("DcaXYGl", &DcaXYGl);
  MatchedPairs->SetBranchAddress("Flag", &Flag);
  MatchedPairs->SetBranchAddress("FitPts", &FitPts);
  MatchedPairs->SetBranchAddress("DedxPts", &DedxPts);
  MatchedPairs->SetBranchAddress("AllPts", &AllPts);
  MatchedPairs->SetBranchAddress("NPossible", &NPossible);
  MatchedPairs->SetBranchAddress("ParentGeantId", &ParentGeantId);
  MatchedPairs->SetBranchAddress("GeantId", &GeantId);
  MatchedPairs->SetBranchAddress("NCommonHit", &NCommonHit);

  // --------------------------------------------------------------------------------
  // -- get MC ntuple
  float pRefMult, pRefMultCorrected, pCentralityWeight, pCentrality16, pVertexX, pVertexY, pVertexZ, pPtMc, pPzMc, pEtaMc, pPhiMc, pParentGeantId, pGeantId;
  McTrack->SetBranchAddress("RefMult", &pRefMult);
  McTrack->SetBranchAddress("RefMultCorrected", &pRefMultCorrected);
  McTrack->SetBranchAddress("CentralityWeight", &pCentralityWeight);
  McTrack->SetBranchAddress("Centrality16", &pCentrality16);
  McTrack->SetBranchAddress("VertexX", &pVertexX);
  McTrack->SetBranchAddress("VertexY", &pVertexY);
  McTrack->SetBranchAddress("VertexZ", &pVertexZ);
  McTrack->SetBranchAddress("PtMc", &pPtMc);
  McTrack->SetBranchAddress("PzMc", &pPzMc);
  McTrack->SetBranchAddress("EtaMc", &pEtaMc);
  McTrack->SetBranchAddress("PhiMc", &pPhiMc);
  McTrack->SetBranchAddress("ParentGeantId", &pParentGeantId);
  McTrack->SetBranchAddress("GeantId", &pGeantId);
  
  // --------------------------------------------------------------------------------
  // -- open output file
  // --------------------------------------------------------------------------------

  TFile *fOut = TFile::Open(out_file, "RECREATE");
  fOut->cd();
  
  // --------------------------------------------------------------------------------
  // -- Create histograms
  // --------------------------------------------------------------------------------

  const int   nCent       = 9;
  const char* cent[]      = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};

  const int    Nbins      = 20;
  const double pt_bin[21] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.4,1.7,2.0,2.3,2.7,3.5,4.0,4.5,5.0};

  // --------------------------------------------------------------------------------
  // -- Initialize pt hists
  TH1D* hpt_f[nCent];
  TH1D* hpt_mc[nCent];

  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    hpt_f[idxCent]  = new TH1D(Form("hpt_f_%s",  cent[idxCent]), "", Nbins, 0., 5.);
    hpt_mc[idxCent] = new TH1D(Form("hpt_mc_%s", cent[idxCent]), "", Nbins, 0., 5.);

    hpt_f[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
    hpt_mc[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
  }

  // --------------------------------------------------------------------------------
  // -- Initialize delta pt hists
  TH1D* hdeltapt[nCent];  
  TH1D* hNdeltapt[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    hdeltapt[idxCent]  = new TH1D(Form("hdeltapt_%s",  cent[idxCent]), "", 100, 0., 10.);
    hNdeltapt[idxCent] = new TH1D(Form("hNdeltapt_%s", cent[idxCent]), "", 100, 0., 10.);

    hdeltapt[idxCent]->Sumw2();
    hNdeltapt[idxCent]->Sumw2();
  }

  // --------------------------------------------------------------------------------
  // -- Initialize momemtum resolution hists
  TH2D* hmomenRes      = new TH2D("hmomenRes", "", Nbins, 0., 5.,20*Nbins,-3.,6.);
  hmomenRes->GetXaxis()->Set(Nbins, pt_bin);
  hmomenRes->Sumw2();
  TH2D* hmomenResDiff  = new TH2D("hmomenResDiff", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  hmomenResDiff->GetXaxis()->Set(Nbins, pt_bin);
  hmomenResDiff->Sumw2();

  TH2D* hmomenResDiffCent[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    hmomenResDiffCent[idxCent] = new TH2D(Form("hmomenResDiff_%s", cent[idxCent]), "", Nbins, 0., 5.,20*Nbins,-50.,50.);
    hmomenResDiffCent[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
    hmomenResDiffCent[idxCent]->Sumw2(); 
  }

  // --------------------------------------------------------------------------------
  // -- Initialize eta hists
  TH1D* heta_f[nCent];
  TH1D* heta_mc[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    heta_f[idxCent]  = new TH1D(Form("heta_f_%s",  cent[idxCent]), "", 31, -1.5, 1.5);
    heta_mc[idxCent] = new TH1D(Form("heta_mc_%s", cent[idxCent]), "", 31, -1.5, 1.5);
  }

  // --------------------------------------------------------------------------------
  // -- efficiency histograms

  TH1D* hcent_f     = new TH1D("hcent_f", "", 1000, 0.5, 1000.5);
  TH1D* hcent_mc    = new TH1D("hcent_mc", "", 1000, 0.5, 1000.5);
  TH1D* hphi_f      = new TH1D("hphi_f", "", 96,-pi,pi);
  TH1D* hphi_mc     = new TH1D("hphi_mc", "",96,-pi,pi);
  TH1D* hvz_f       = new TH1D("hvz_f", "", 100, -max_z, max_z);
  TH1D* hvz_mc      = new TH1D("hvz_mc", "", 100, -max_z, max_z);

  TH2D* hetacent_f  = new TH2D("hetacent_f", "", 8, -1., 1.,1000,0.5,1000.5);
  TH2D* hetacent_mc = new TH2D("hetacent_mc", "", 8, -1., 1.,1000,0.5,1000.5);
  TH2D* hetavz_f    = new TH2D("hetavz_f", "", 80, -1., 1.,110,-max_z,max_z);
  TH2D* hetavz_mc   = new TH2D("hetavz_mc", "", 80, -1., 1.,110,-max_z,max_z);
  TH2D* hcentvz_f   = new TH2D("hcentvz_f", "", 1000, 0.5, 1000.5,10,-max_z,max_z);
  TH2D* hcentvz_mc  = new TH2D("hcentvz_mc", "", 1000, 0.5, 1000.5,10,-max_z,max_z);

  TH2D* hpteta_f    = new TH2D("hpteta_f", "", Nbins, 0., 5.,8,-1.,1.);
  hpteta_f->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hpteta_mc   = new TH2D("hpteta_mc", "", Nbins, 0., 5.,8,-1.,1.);
  hpteta_mc->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptvz_f     = new TH2D("hptvz_f", "", Nbins, 0., 5.,10,-max_z,max_z);
  hptvz_f->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptvz_mc    = new TH2D("hptvz_mc", "", Nbins, 0., 5.,10,-max_z,max_z);
  hptvz_mc->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptphi_f    = new TH2D("hptphi_f", "", Nbins, 0., 5.,48,-pi,pi);
  hptphi_f  ->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptphi_mc   = new TH2D("hptphi_mc", "", Nbins, 0., 5.,48,-pi,pi);
  hptphi_f->GetXaxis()->Set(Nbins, pt_bin);

  TH2D* hptcent_f   = new TH2D("hptcent_f", "", Nbins, 0., 5.,1000,0.5,1000.5);
  hptcent_f->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptcent_mc  = new TH2D("hptcent_mc", "", Nbins, 0., 5.,1000,0.5,1000.5);
  hptcent_mc->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptcent2_f  = new TH2D("hptcent2_f", "", Nbins, 0., 5.,9,-0.5,8.5);
  hptcent2_f->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hptcent2_mc = new TH2D("hptcent2_mc", "", Nbins, 0., 5.,9,-0.5,8.5);
  hptcent2_mc->GetXaxis()->Set(Nbins, pt_bin);

  // --------------------------------------------------------------------------------
  // -- QA histograms

  TH1D* hDCA      = new TH1D("hDCA",    "", 100, -5., 5.);
  TH2D* hDCAphi   = new TH2D("hDCAphi", "", 100, -5., 5., 96, -pi, pi);
  TH2D* hphipt    = new TH2D("hphipt",  "", 200, -pi, pi, 200, 0.,10.);

  // --------------------------------------------------------------------------------
  // -- Loop over MC tracks
  // --------------------------------------------------------------------------------
  int nTracks = McTrack->GetEntries();

  int nTracksCut = 0;
  int nTracksCutEta = 0;
  
  for (int idxTrack = 0; idxTrack < nTracks; idxTrack++){
    McTrack->GetEntry(idxTrack);

    // --------------------------------------------------------------------------------
    // -- Get centrality 
    int idxCent = 0;

    if ( pCentrality16 == 15)     idxCent = 0;
    else if (pCentrality16 == 14) idxCent = 1;
    else if (pCentrality16 == 12 || pCentrality16 == 13) idxCent = 2;
    else if (pCentrality16 == 10 || pCentrality16 == 11) idxCent = 3;
    else if (pCentrality16 == 8  || pCentrality16 == 9)  idxCent = 4;
    else if (pCentrality16 == 6  || pCentrality16 == 7)  idxCent = 5;
    else if (pCentrality16 == 4  || pCentrality16 == 5)  idxCent = 6;
    else if (pCentrality16 == 2  || pCentrality16 == 3)  idxCent = 7;
    else if (pCentrality16 == 0  || pCentrality16 == 1)  idxCent = 8;
    else 
      continue;

    // --------------------------------------------------------------------------------
    // -- basic check cuts
    if (pParentGeantId != 0)
      continue;

    if (pGeantId != PID )
      continue;

    // -- within Vertex R
    if ( (pVertexX-vxo)*(pVertexX-vxo) + (pVertexY-vyo)*(pVertexY-vyo) > max_r*max_r )
      continue;

    // --------------------------------------------------------------------------------
    // -- check cuts

    // -- within Vertex Z
    int fVz = (pVertexZ < max_z && pVertexZ > -max_z) ? 0 : 1;	    

    // -- pt window
    int fPt = (pPtMc*pPtMc > min_pt*min_pt && pPtMc*pPtMc < max_pt*max_pt) ? 0 : 1;
   
    // -- with eta window
    int fEta = (fabs(pEtaMc) < max_eta) ? 0 : 1;

    // -- refMult usefull
    int fRef = (pRefMultCorrected > 0 && pRefMultCorrected < 1000) ? 0 : 1;

    // -- -- -- -- -- -- -- -- -- -- 

    if ((fVz+fEta+fPt) == 0) 
      hptcent2_mc->Fill(pPtMc,idxCent);

    // apply refmult cut
    if (fRef) 
      continue;

    if (fPt == 0)       
      hetavz_mc->Fill(pEtaMc,pVertexZ);

    if ((fEta+fPt) == 0)  {
      hvz_mc->Fill(pVertexZ);
      hptvz_mc->Fill(pPtMc,pVertexZ);
      hcentvz_mc->Fill(pRefMultCorrected,pVertexZ);
    }   

    // -- apply vz cut
    if (fVz) 
      continue;

    hpteta_mc->Fill(pPtMc, pEtaMc);

    // --------------------------------------------------------------------------------
    // --- Fill pt hists - MC
    if (fEta == 0) {
      ++nTracksCut;
      hpt_mc[idxCent]->Fill(pPtMc);

      hcent_mc->Fill(pRefMultCorrected);
      hptcent_mc->Fill(pPtMc, pRefMultCorrected);
      hptphi_mc->Fill(pPtMc, pPhiMc);
      hphi_mc->Fill(pPhiMc);
    }

    // --------------------------------------------------------------------------------
    // --- Fill eta hists - eta
    if (fPt == 0) {
      ++nTracksCutEta;
      heta_mc[idxCent]->Fill(pEtaMc);
      hetacent_mc->Fill(pEtaMc, pRefMultCorrected);
    }

    // --------------------------------------------------------------------------------
  } // for (int idxTrack = 0; idxTrack < nTracks; idxTrack++){

  // --------------------------------------------------------------------------------
  // -- Loop over matched rec tracks
  // --------------------------------------------------------------------------------
  int nMatched = MatchedPairs->GetEntries();
  int nMatchedCut = 0;
  int nMatchedCutEta = 0;

  for (int idxMatch = 0; idxMatch < nMatched; idxMatch++){
    MatchedPairs->GetEntry(idxMatch);

    // --------------------------------------------------------------------------------
    // -- Get centrality 
    int idxCent = 0;
    if      (Centrality16 == 15) idxCent = 0;
    else if (Centrality16 == 14) idxCent = 1;
    else if (Centrality16 == 12 || Centrality16 == 13) idxCent = 2;
    else if (Centrality16 == 10 || Centrality16 == 11) idxCent = 3;
    else if (Centrality16 == 8  || Centrality16 == 9)  idxCent = 4;
    else if (Centrality16 == 6  || Centrality16 == 7)  idxCent = 5;
    else if (Centrality16 == 4  || Centrality16 == 5)  idxCent = 6;
    else if (Centrality16 == 2  || Centrality16 == 3)  idxCent = 7;
    else if (Centrality16 == 0  || Centrality16 == 1)  idxCent = 8;
    else 
      continue;

    // --------------------------------------------------------------------------------
    // -- basic check cuts
    if (ParentGeantId != 0)
      continue;
   
    if (GeantId != PID )
      continue;

    if ((VertexX-vxo)*(VertexX-vxo) + (VertexY-vyo)*(VertexY-vyo) > max_r*max_r) 
      continue;
    
    if (Dedx < 0) 
      continue;

    if (DedxPts < min_dedx)
      continue;

    if (Flag < 0 || Flag > 699)
      continue;

    if (PtPr > 10./7.*PtGl || PtPr < 7./10.*PtGl) 
      continue;

    if (NCommonHit <= min_nCommon) 
      continue;

    if (FitPts <= min_fitpts)
      continue;
    
    if (FitPts/NPossible <= min_fitpts_nposs)
      continue;
    
    if (fabs(DcaGl) > max_dca )
      continue;

    // --------------------------------------------------------------------------------
    // -- check cuts

    int fVz   = (VertexZ < max_z && VertexZ > -max_z) ? 0 : 1;

    int fPt  = (PtMc*PtMc > min_pt*min_pt && PtMc*PtMc < max_pt*max_pt)  ? 0 : 1;
    int fEta = (fabs(EtaMc) < max_eta )  ? 0 : 1;

    // --------------------------------------------------------------------------------

    if (fPt == 0)
      hetavz_f->Fill(EtaMc,VertexZ);
    
    if ((fEta+fPt) == 0) {
      hvz_f->Fill(VertexZ);
      hptvz_f->Fill(PtMc,VertexZ);
      hcentvz_f->Fill(RefMultCorrected,VertexZ);
    }
    
    // -- apply vz cut
    if (fVz)
      continue;
    
    if ((fEta+fPt) == 0)
      hptcent2_f->Fill(PtMc,idxCent);
    
    hpteta_f->Fill(PtMc,EtaMc);

    // --------------------------------------------------------------------------------
    // --- Fill pt hists - matched
    if (fEta == 0) {
      ++nMatchedCut;
      hpt_f[idxCent]->Fill(PtMc);

      hcent_f->Fill(RefMultCorrected);
      hptcent_f->Fill(PtMc,RefMultCorrected);  
      hptphi_f->Fill(PtMc,PhiPr);
      hphi_f->Fill(PhiPr);

      hphipt->Fill(PhiPr, PtMc);
      hDCA->Fill(DcaXYGl);
      hDCAphi->Fill(DcaXYGl,PhiPr);

      hNdeltapt[idxCent]->Fill(PtMc);
      hdeltapt[idxCent]->Fill(PtMc, (PtMc-PtPr));

      if (FitPts == NCommonHit)
	hmomenRes->Fill(PtMc, 1./PtPr-1./PtMc);
      hmomenResDiff->Fill(PtMc, PtPr-PtMc);
      hmomenResDiffCent[idxCent]->Fill(PtMc,PtPr-PtMc);
    }
    
    // --------------------------------------------------------------------------------
    // --- Fill eta hists - matched
    if (fPt == 0) {
      ++nMatchedCutEta;
      heta_f[idxCent]->Fill(EtaMc);
      hetacent_f->Fill(EtaMc,RefMultCorrected);
    }
  } //   for (int idxMatch = 0; idxMatch < nMatched; idxMatch++){

  cout << "nMatched: " << nMatched << endl;
  cout << "nTracks:  " << nTracks << endl;

  cout << "nMatchedCut (pt) : " << nMatchedCut << endl;
  cout << "nMatchedCut (eta): " << nMatchedCutEta << endl;

  cout << "nTracksCut (pt)  : " << nTracksCut << endl;
  cout << "nTracksCut (eta) : " << nTracksCutEta << endl;
  
  // =============================================================================================
  // -- Make efficiencies 
  // =============================================================================================

  TH1D* hpt[nCent];
  for (int idxCent = 0; idxCent < nCent; idxCent++) 
    hpt[idxCent] = getEffRatio(hpt_mc[idxCent], hpt_f[idxCent], Form("hpt_%s", cent[idxCent]));
  
  TH1D* heta[nCent];
  for (int idxCent = 0; idxCent < nCent; idxCent++) 
    heta[idxCent] = getEffRatio(heta_mc[idxCent], heta_f[idxCent], Form("heta_%s", cent[idxCent]));
  
  TH1D* hcent    = getEffRatio(hcent_mc,    hcent_f,    "hcent");
  TH1D* hphi     = getEffRatio(hphi_mc,     hphi_f,     "hphi");
  TH1D* hvz      = getEffRatio(hvz_mc,      hvz_f,      "hvz");

  TH2D* hetacent = getEffRatio(hetacent_mc, hetacent_f, "hetacent");
  TH2D* hetavz   = getEffRatio(hetavz_mc,   hetavz_f,   "hetavz");
  TH2D* hcentvz  = getEffRatio(hcentvz_mc,  hcentvz_f,  "hcentvz");

  TH2D* hpteta   = getEffRatio(hpteta_mc,   hpteta_f,   "hpteta");
  TH2D* hptvz    = getEffRatio(hptvz_mc,    hptvz_f,    "hptvz");
  TH2D* hptphi   = getEffRatio(hptphi_mc,   hptphi_f,   "hptphi");

  TH2D* hptcent  = getEffRatio(hptcent_mc,  hptcent_f,  "hptcent");
  TH2D* hptcent2 = getEffRatio(hptcent2_mc, hptcent2_f, "hptcent2");

  // --------------------------------------------------------------------------------

  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) 
    hdeltapt[idxCent]->Divide(hNdeltapt[idxCent]);

  // =============================================================================================
  // -- Draw Canvas
  // =============================================================================================

  gSystem->Exec(Form("mkdir -p effQA/%d", energy));

  // --------------------------------------------------------------------------------

  TCanvas *cmomRes = new TCanvas("cmomRes", "momRes", 800, 600);
  hmomenRes->Draw("colz");
  cmomRes->SaveAs(Form("effQA/%d/hmomRes%d_%d.png",energy,energy,PID));

  // --------------------------------------------------------------------------------

  TCanvas *cdeltapt = new TCanvas("cdeltapt", "deltapt", 800, 600);
  cdeltapt->Divide(3,3,0,0);

  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    cdeltapt->cd(idxCent+1);
    hdeltapt[idxCent]->SetTitle(Form("Efficiency for refmult integrated p_{T} distribution for %s at %dGeV - cent %d",particle,energy,idxCent));
    hdeltapt[idxCent]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    hdeltapt[idxCent]->GetXaxis()->CenterTitle();
    hdeltapt[idxCent]->GetXaxis()->SetRangeUser(0.,5.);
    hdeltapt[idxCent]->GetYaxis()->SetTitle("<#Deltap_{T}>(GeV/c)");
    hdeltapt[idxCent]->GetYaxis()->SetTitleOffset(1.4);
    hdeltapt[idxCent]->GetYaxis()->CenterTitle();
    hdeltapt[idxCent]->SetMarkerStyle(20);
    hdeltapt[idxCent]->SetStats(kFALSE);
    hdeltapt[idxCent]->DrawCopy();
  }    
  cdeltapt->SaveAs(Form("effQA/%d/dhdeltapt%d_%d.png",energy,energy,PID));

  // --------------------------------------------------------------------------------

  TCanvas *cphi = new TCanvas("cphi", "phi", 800, 600);
  hphi->GetXaxis()->SetTitle("#phi");
  hphi->GetXaxis()->CenterTitle();
  hphi->GetYaxis()->SetTitle("Efficiency");
  hphi->GetYaxis()->CenterTitle();
  hphi->SetStats(kFALSE);
  hphi->DrawCopy();
  hphi->GetYaxis()->SetRangeUser(0.,1.2);
  cphi->SaveAs(Form("effQA/%d/hphi%d_%d.png",energy,energy,PID));
    
  TCanvas *cvz = new TCanvas("cvz", "vz", 800, 600);
  hvz->GetXaxis()->SetTitle("V_{Z}");
  hvz->GetXaxis()->CenterTitle();
  hvz->GetYaxis()->SetTitle("Efficiency");
  hvz->GetYaxis()->CenterTitle();
  hvz->SetStats(kFALSE);
  hvz->DrawCopy();
  cvz->SaveAs(Form("effQA/%d/hvz%d_%d.png",energy,energy,PID));
  
  TCanvas *ccent = new TCanvas("ccent", "cent", 800, 600);
  hcent->GetXaxis()->SetTitle("refMult");
  hcent->GetXaxis()->CenterTitle();
  hcent->GetYaxis()->SetTitle("Efficiency");
  hcent->GetYaxis()->CenterTitle();
  hcent->SetStats(kFALSE);
  hcent->DrawCopy();
  ccent->SaveAs(Form("effQA/%d/hcent%d_%d.png",energy,energy,PID));
  
  TCanvas *cptcent = new TCanvas("cptcent", "ptcent", 800, 600);
  hptcent->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptcent->GetXaxis()->CenterTitle();
  hptcent->GetXaxis()->SetRangeUser(0.2,4.5);
  hptcent->GetYaxis()->SetTitle("refMult");
  hptcent->GetYaxis()->CenterTitle();
  hptcent->SetStats(kFALSE);
  hptcent->DrawCopy("colz");
  cptcent->SaveAs(Form("effQA/%d/hptcent%d_%d.png",energy,energy,PID));

  TCanvas *cptcent2 = new TCanvas("cptcent2", "ptcent2", 800, 600);
  hptcent2->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptcent2->GetXaxis()->CenterTitle();
  hptcent2->GetXaxis()->SetRangeUser(0.2,4.5);
  hptcent2->GetYaxis()->SetTitle("centrality");
  hptcent2->GetYaxis()->CenterTitle();
  hptcent2->SetStats(kFALSE);
  hptcent2->DrawCopy("colz");
  cptcent2->SaveAs(Form("effQA/%d/hptcent2%d_%d.png",energy,energy,PID));
  
  TCanvas *cpteta = new TCanvas("cpteta", "pteta", 800, 600);
  hpteta->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hpteta->GetXaxis()->CenterTitle();
  hpteta->GetXaxis()->SetRangeUser(0.2,4.5);
  hpteta->GetYaxis()->SetTitle("#eta");
  hpteta->GetYaxis()->CenterTitle();
  hpteta->SetStats(kFALSE);
  hpteta->DrawCopy("colz");
  cpteta->SaveAs(Form("effQA/%d/hpteta%d_%d.png",energy,energy,PID));
  
  TCanvas *cptvz = new TCanvas("cptvz", "ptvz", 800, 600);
  hptvz->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptvz->GetXaxis()->CenterTitle();
  hptvz->GetXaxis()->SetRangeUser(0.2,4.5);
  hptvz->GetYaxis()->SetTitle("V_{Z}");
  hptvz->GetYaxis()->CenterTitle();
  hptvz->SetStats(kFALSE);
  hptvz->DrawCopy("colz");
  cptvz->SaveAs(Form("effQA/%d/hptvz%d_%d.png",energy,energy,PID));
  
  TCanvas *cetacent = new TCanvas("cetacent", "etacent", 800, 600);
  hetacent->GetXaxis()->SetTitle("#eta");
  hetacent->GetXaxis()->CenterTitle();
  hetacent->GetYaxis()->SetTitle("refMult");
  hetacent->GetYaxis()->CenterTitle();
  hetacent->SetStats(kFALSE);
  hetacent->DrawCopy("colz");
  cetacent->SaveAs(Form("effQA/%d/hetacent%d_%d.png",energy,energy,PID));
  
  TCanvas *cetavz = new TCanvas("cetavz", "etavz", 800, 600);
  hetavz->GetXaxis()->SetTitle("#eta");
  hetavz->GetXaxis()->CenterTitle();
  hetavz->GetYaxis()->SetTitle("V_{Z}");
  hetavz->GetYaxis()->CenterTitle();
  hetavz->SetStats(kFALSE);
  hetavz->DrawCopy("colz");
  cetavz->SaveAs(Form("effQA/%d/hetavz%d_%d.png",energy,energy,PID));
  
  TCanvas *ccentvz = new TCanvas("ccentvz", "centvz", 800, 600);
  hcentvz->GetXaxis()->SetTitle("refmult");
  hcentvz->GetXaxis()->CenterTitle();
  hcentvz->GetYaxis()->SetTitle("V_{Z}");
  hcentvz->GetYaxis()->CenterTitle();
  hcentvz->SetStats(kFALSE);
  hcentvz->DrawCopy("colz");
  ccentvz->SaveAs(Form("effQA/%d/hcentvz%d_%d.png",energy,energy,PID));
  
  TCanvas *cptphi = new TCanvas("cptphi", "ptphi", 800, 600);
  hptphi->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptphi->GetXaxis()->CenterTitle();
  hptphi->GetXaxis()->SetRangeUser(0.2,4.5);
  hptphi->GetYaxis()->SetTitle("#phi");
  hptphi->GetYaxis()->CenterTitle();
  hptphi->SetStats(kFALSE);
  hptphi->DrawCopy("colz");
  cptphi->SaveAs(Form("effQA/%d/hptphi%d_%d.png",energy,energy,PID));

  TCanvas *cDCA = new TCanvas("cDCA", "DCA", 800, 600);
  hDCA->GetXaxis()->SetTitle("DCA(cm)");
  hDCA->GetXaxis()->CenterTitle();
  hDCA->GetYaxis()->SetTitle("dN/dDCA(cm^{-1})");
  hDCA->GetYaxis()->CenterTitle();
  hDCA->SetStats(kFALSE);
  hDCA->DrawCopy("");
  cDCA->SaveAs(Form("effQA/%d/hDCA%d_%d.png",energy,energy,PID));

  TCanvas *cDCAphi = new TCanvas("cDCAphi", "DCAphi", 800, 600);
  hDCAphi->GetXaxis()->SetTitle("DCA(cm)");
  hDCAphi->GetXaxis()->CenterTitle();
  hDCAphi->GetYaxis()->SetTitle("#phi");
  hDCAphi->GetYaxis()->CenterTitle();
  hDCAphi->SetStats(kFALSE);
  hDCAphi->DrawCopy("colz");
  cDCAphi->SaveAs(Form("effQA/%d/hDCAphi%d_%d.png",energy,energy,PID));
  
  fOut->Write();
  fOut->Close();
}

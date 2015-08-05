#include "TCanvas.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "math.h"
#include <iostream>

using namespace std;

// _______________________________________________________________
TH1D* getEffRatio(TH1F* hMC, TH1D* hRec, const char*name) {
  TH1D* hEff = static_cast<TH1D*>(hMC->Clone(name));

  for (int idx = 1; idx <= hEff->GetXaxis()->GetNbins(); idx++) {
    double binContent = hEff->GetBinContent(idx);
    if (binContent < 0.00001)
      binContent = 1.;
    double ratio = hRec[idxCent]->GetBinContent(idx)/(double)binContent;
    hEff->SetBinContent(idx, ratio);
    hEff->SetBinError(idx, sqrt((double)binContent*ratio*(1.0-ratio))/(double)binContent);
  }
 
  return hEff;
}

// _______________________________________________________________
TH2D* getEffRatio(TH2F* hMC, TH2D* hRec, const char*name) {
  TH2D* hEff = static_cast<TH2D*>(hMC->Clone(name));

  for (int idxX = 1; idxY <= hEff->GetXaxis()->GetNbins(); idxX++) {
    for (int idxY = 1; idxY <= hEff->GetYaxis()->GetNbins(); idxY++) {
      double binContent = hEff->GetBinContent(idxX, idxY);
      if (binContent < <0.00001)
	binContent = 1.;
      double ratio = hRec[idxCent]->GetBinContent(idxX, idxY)/(double)binContent;
      hEff->SetBinContent(idxX, idxY, ratio);
      hEff->SetBinError(idxX, idxY, sqrt((double)binContent*ratio*(1.0-ratio))/(double)binContent);
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
  double  mass;
  float   PID = 0;
  
  if(TString(particle) == TString("piplus")){
    in_file = prefix + TString("/SinglePiPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piplus") + NRG + TString("GeV.root");
    mass = .13975;
    PID = 8;
  }
  else if(TString(particle) == TString("piminus")){
    in_file = prefix + TString("/SinglePiMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piminus") + NRG + TString("GeV.root");
    mass = .13975;
    PID = 9;
  }
  else if(TString(particle) == TString("kaonplus")){
    in_file = prefix + TString("/SingleKPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kplus") + NRG + TString("GeV.root");
    mass = .493677;
    PID = 11;
  }
  else if(TString(particle) == TString("kaonminus")){
    in_file = prefix + TString("/SingleKMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kminus") + NRG + TString("GeV.root");
    mass = .493677;
    PID = 12;
  }
  else if(TString(particle) == TString("protonplus")){
    in_file = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pplus") + NRG + TString("GeV.root");
    mass = .93827;
    PID = 14;
  }
  else if(TString(particle) == TString("protonminus")){
    in_file = prefix + TString("/SinglePbarNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pminus") + NRG + TString("GeV.root");
    mass = .93827;
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
  cout << "opening output file: " << out_file << endl;

  TFile *fOut = TFile::Open(out_file, "RECREATE");
  TFile* fmomenOut =TFile::Open("fmomenOut","UPDATE");   /// does what?
  fOut->cd();
  
  // --------------------------------------------------------------------------------
  // -- Create histograms
  // --------------------------------------------------------------------------------

  const int   nCent   = 10;
  const char* cent[]  = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};

  const int    Nbins      = 20;
  const double pt_bin[21] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,.0,1.1,1.4,1.7,2.0,2.3,2.7,3.5,4.0,4.5,5.0};

  // --------------------------------------------------------------------------------

  // -- Initialize pt hists
  TH1D* hpt_f[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    hpt_f[idxCent]  = new TH1D(Form("hpt_f_%s",  cent[idxCent]), "", Nbins, 0., 5.);
    hpt_mc[idxCent] = new TH1D(Form("hpt_mc_%s", cent[idxCent]), "", Nbins, 0., 5.);

    hpt_f[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
    hpt_mc[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
  }

  // --------------------------------------------------------------------------------

  // -- Initialize eta hists
  TH1D* heta_f[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    heta_f[idxCent]  = new TH1D(Form("heta_f_%s",  cent[idxCent]), "", 31, -1.5, 1.5);
    heta_mc[idxCent] = new TH1D(Form("heta_mc_%s", cent[idxCent]), "", 31, -1.5, 1.5);
  }

  // --------------------------------------------------------------------------------
  // -- efficiency histograms

  TH1D* hcent_f     = new TH1D("hcent_f", "", 1000, 0.5, 1000.5);
  TH1D* hcent_mc    = new TH1D("hcent_mc", "", 1000, 0.5, 1000.5);
  TH1D* heta_f      = new TH1D("heta_f", "", 44, -1.1, 1.1);
  TH1D* heta_mc     = new TH1D("heta_mc", "", 44, -1.1, 1.1);
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

  TH1D* hDCA      = new TH1D("hDCA", "", 100, -5., 5.);
  TH2D* hDCAphi   = new TH2D("hDCAphi", "", 100, -5., 5.,96, -pi, pi);

  // --------------------------------------------------------------------------------

  TH2D* hmomenRes      = new TH2D("hmomenRes", "", Nbins, 0., 5.,20*Nbins,-3.,6.);
  hmomenRes->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hmomenResDiff  = new TH2D("hmomenResDiff", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  hmomenResDiff->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hmomenResDiffC = new TH2D("hmomenResDiffC", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  hmomenResDiffC->GetXaxis()->Set(Nbins, pt_bin);
  TH2D* hmomenResDiffP = new TH2D("hmomenResDiffP", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  hmomenResDiffP->GetXaxis()->Set(Nbins, pt_bin);

  hmomenResDiff->Sumw2(); hmomenResDiffC->Sumw2(); hmomenResDiffP->Sumw2(); hmomenRes->Sumw2();

  // --------------------------------------------------------------------------------

  TH2D* hphipt     = new TH2D("hphipt","",200,-pi,pi,200,0.,10.);

  TH1D* hNdeltaptC = new TH1D("hNdeltaptC","",100,0.,10.);
  TH1D* hdeltaptC  = new TH1D("hdeltaptC","",100,0.,10.);
  TH1D* hNdeltaptP = new TH1D("hNdeltaptP","",100,0.,10.);
  TH1D* hdeltaptP  = new TH1D("hdeltaptP","",100,0.,10.);

  hNdeltaptC->Sumw2();hdeltaptC->Sumw2(); 
  hNdeltaptP->Sumw2();hdeltaptP->Sumw2();

  // --------------------------------------------------------------------------------

  Double_t xbins[38] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
			2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,4.5,5.0,6.0,7.0,10.0};

  TH1D* hMomenMCCent      = new TH1D("hMomenMCCent","",37,xbins);
  TH1D* hMomenMatchedCent = new TH1D("hMomenMatchedCent","",37,xbins);
  TH1D* hMomenMCPer       = new TH1D("hMomenMCPer","",37,xbins);
  TH1D* hMomenMatchedPer  = new TH1D("hMomenMatchedPer","",37,xbins);

  hMomenMatchedCent->Sumw2(); hMomenMCPer->Sumw2(); hMomenMatchedPer->Sumw2(); hMomenMCCent->Sumw2();
  
  cout << "histograms created" << endl;
    
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
   
    if (pGeantId == PID )
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

    if ((fVZ+fEta+fPt) == 0) 
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
    if (fVZ) 
      continue;

    hpteta_mc->Fill(pPtMc, pEtaMc);

    // --------------------------------------------------------------------------------
    // --- Fill pt hists - MC
    if (fEta == 0) {
      nTracksCut += 1;
      hpt_mc[idxCent]->Fill(pPtMc);

      hcent_mc->Fill(pRefMultCorrected);
      hptcent_mc->Fill(pPtMc, pRefMultCorrected);
      hptphi_mc->Fill(pPtMc, pPhiMc);
      hphi_mc->Fill(pPhiMc);
    }

    // --------------------------------------------------------------------------------
    // --- Fill eta hists - eta
    if (fPt == 0) {
      nTracksCutEta += 1;
      hpt_eta[idxCent]->Fill(pEtaMc);

      heta_mc->Fill(pEtaMc);
      hetacent_mc->Fill(pEtaMc, pRefMultCorrected);
    }

    // --------------------------------------------------------------------------------
  } // for (int idxTrack = 0; idxTrack < nTracks; idxTrack++){
  
  // --------------------------------------------------------------------------------
  // -- Loop over matched rec tracks
  // --------------------------------------------------------------------------------
  int nMatched = MatchedPairs->GetEntries();
  int nMatchedCut = 0;
  int nmatchedCutEta = 0;

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
   
    if (GeantId == PID )
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

    int VZ   = (VertexZ < max_z && VertexZ > -max_z) ? 0 : 1;

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
    if (fVZ)
      continue;
    
    if ((fEta+fPt) == 0)
      hptcent2_f->Fill(PtMc,idxCent);
    
    hpteta_f->Fill(PtMc,EtaMc);

    // --------------------------------------------------------------------------------
    // --- Fill pt hists - matched
    if (fEta == 0) {
      nMatchedCut += 1;
      hpt_f[idxCent]->Fill(PtMc);

      hcent_f->Fill(RefMultCorrected);
      hptcent_f->Fill(PtMc,RefMultCorrected);  
      hptphi_f->Fill(PtMc,PhiPr);
      hphi_f->Fill(PhiPr);

      hphipt->Fill(PhiPr, PtMc);
      hDCA->Fill(DcaXYGl);
      hDCAphi->Fill(DcaXYGl,PhiPr);

      if (FitPts == NCommonHit)
	hmomenRes->Fill(PtMc, 1./PtPr-1./PtMc);

      hmomenResDiff->Fill(PtMc, PtPr-PtMc);

      if (idxCent == 0) {
	hmomenResDiffC->Fill(PtMc,PtPr-PtMc);
	
	hNdeltaptC->Fill(PtMc);
	hdeltaptC->Fill(PtMc,(PtMc-PtPr));
      }

      if (idxCent == 7 || idxCent == 8) {
	hNdeltaptP->Fill(PtMc);
	hdeltaptP->Fill(PtMc,(PtMc-PtPr));
	hmomenResDiffP->Fill(PtMc,PtPr-PtMc);
      }
    }
    
    // --------------------------------------------------------------------------------
    // --- Fill eta hists - matched
    if (fPt == 0) {
      nmatchedcutEta += 1;
      heta_f[idxCent]->Fill(EtaMc);

      heta_f->Fill(EtaMc);
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
  // -- Mke efficiencies 
  // =============================================================================================

  TH1D* hpt[nCent];
  for (int idxCent = 0; idxCent < nCent; idxCent++) 
    hpt[idxCent] = getEffRatio(hpt_mc[idxCent], hpt_f[idxCent], Form("hpt_%s", cent[idxCent]));
  
  TH1D* heta[nCent];
  for (int idxCent = 0; idxCent < nCent; idxCent++) 
    heta[idxCent] = getEffRatio(heta_mc[idxCent], heta_f[idxCent], Form("heta_%s", cent[idxCent]));
  
  TH1D* hcent    = getEffRatio(hcent_mc,    hcent_f,    "hcent");
  TH1D* heta     = getEffRatio(heta_mc,     heta_f,     "heta");
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

  hdeltaptC->Divide(hNdeltaptC);
  hdeltaptP->Divide(hNdeltaptP);

  // =============================================================================================
  // -- Draw Canvas
  // =============================================================================================
  
  TCanvas *cmomentumResC = new TCanvas("cmomentumResC", "momentumResC", 800, 600);
  hMomenMatchedCent->Draw("C");

  TCanvas *cmomRes = new TCanvas("cmomRes", "momRes", 800, 600);
  hmomenRes->Draw("colz");

  TCanvas *cmomentumResP = new TCanvas("cmomentumResP", "momentumResP", 800, 600);
  hMomenMatchedPer->Draw("C");

  TCanvas *cdeltaptC = new TCanvas("cdeltaptC", "deltaptC", 800, 600);

  char buffer[100];
  sprintf(buffer,"Efficiency for refmult integrated p_{T} distribution for %s at %dGeV",particle,energy);
  //hdeltaptC->SetTitle(buffer);
  hdeltaptC->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hdeltaptC->GetXaxis()->CenterTitle();
  hdeltaptC->GetXaxis()->SetRangeUser(0.,5.);
  hdeltaptC->GetYaxis()->SetTitle("<#Deltap_{T}>(GeV/c)");
  hdeltaptC->GetYaxis()->SetTitleOffset(1.4);
  hdeltaptC->GetYaxis()->CenterTitle();
  hdeltaptC->SetMarkerStyle(20);
  hdeltaptC->SetStats(kFALSE);
  hdeltaptC->DrawCopy();
  sprintf(buffer, "effQA/hdeltaptC%d_%d.png",energy,PID);
  cdeltaptC->Print(buffer);

  TCanvas *cdeltaptP = new TCanvas("cdeltaptP", "deltaptP", 800, 600);
  sprintf(buffer,"Efficiency for refmult integrated p_{T} distribution for %s at %dGeV",particle,energy);
  //hdeltaptP->SetTitle(buffer);
  hdeltaptP->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hdeltaptP->GetXaxis()->CenterTitle();
  hdeltaptP->GetYaxis()->SetTitleOffset(1.4);
  hdeltaptP->GetXaxis()->SetRangeUser(0.,5.);
  hdeltaptP->GetYaxis()->SetTitle("<#Deltap_{T}>(GeV/c)");
  hdeltaptP->GetYaxis()->CenterTitle();
  hdeltaptP->SetMarkerStyle(20);
  hdeltaptP->SetStats(kFALSE);
  hdeltaptP->DrawCopy();
  sprintf(buffer, "effQA/hdeltaptP%d_%d.png",energy,PID);
  cdeltaptP->Print(buffer);

  TCanvas *cphi = new TCanvas("cphi", "phi", 800, 600);
  hphi->GetXaxis()->SetTitle("#phi");
  hphi->GetXaxis()->CenterTitle();
  hphi->GetYaxis()->SetTitle("Efficiency");
  hphi->GetYaxis()->CenterTitle();
  hphi->SetStats(kFALSE);
  hphi->DrawCopy();
  sprintf(buffer, "effQA/hphi%d_%d.png",energy,PID);
  hphi->GetYaxis()->SetRangeUser(0.,1.2);
  //cphi->Print(buffer);
  
  TCanvas *ceta = new TCanvas("ceta", "eta", 800, 600);
  heta->GetXaxis()->SetTitle("#eta");
  heta->GetXaxis()->CenterTitle();
  heta->GetYaxis()->SetTitle("Efficiency");
  heta->GetYaxis()->CenterTitle();
  heta->SetStats(kFALSE);
  heta->DrawCopy();
  sprintf(buffer, "effQA/heta%d_%d.png",energy,PID);
  ceta->Print(buffer);
  
  TCanvas *cvz = new TCanvas("cvz", "vz", 800, 600);
  hvz->GetXaxis()->SetTitle("V_{Z}");
  hvz->GetXaxis()->CenterTitle();
  hvz->GetYaxis()->SetTitle("Efficiency");
  hvz->GetYaxis()->CenterTitle();
  hvz->SetStats(kFALSE);
  hvz->DrawCopy();
  sprintf(buffer, "effQA/hvz%d_%d.png",energy,PID);
  cvz->Print(buffer);
  
  TCanvas *ccent = new TCanvas("ccent", "cent", 800, 600);
  hcent->GetXaxis()->SetTitle("refMult");
  hcent->GetXaxis()->CenterTitle();
  hcent->GetYaxis()->SetTitle("Efficiency");
  hcent->GetYaxis()->CenterTitle();
  hcent->SetStats(kFALSE);
  hcent->DrawCopy();
  sprintf(buffer, "effQA/hcent%d_%d.png",energy,PID);
  ccent->Print(buffer);
  
  TCanvas *cptcent = new TCanvas("cptcent", "ptcent", 800, 600);
  hptcent->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptcent->GetXaxis()->CenterTitle();
  hptcent->GetXaxis()->SetRangeUser(0.2,4.5);
  hptcent->GetYaxis()->SetTitle("refMult");
  hptcent->GetYaxis()->CenterTitle();
  hptcent->SetStats(kFALSE);
  hptcent->DrawCopy("colz");
  sprintf(buffer, "effQA/hptcent%d_%d.png",energy,PID);
  cptcent->Print(buffer);

  TCanvas *cptcent2 = new TCanvas("cptcent2", "ptcent2", 800, 600);
  hptcent2->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptcent2->GetXaxis()->CenterTitle();
  hptcent2->GetXaxis()->SetRangeUser(0.2,4.5);
  hptcent2->GetYaxis()->SetTitle("centrality");
  hptcent2->GetYaxis()->CenterTitle();
  hptcent2->SetStats(kFALSE);
  hptcent2->DrawCopy("colz");
  sprintf(buffer, "effQA/hptcent2%d_%d.png",energy,PID);
  cptcent2->Print(buffer);
  
  TCanvas *cpteta = new TCanvas("cpteta", "pteta", 800, 600);
  hpteta->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hpteta->GetXaxis()->CenterTitle();
  hpteta->GetXaxis()->SetRangeUser(0.2,4.5);
  hpteta->GetYaxis()->SetTitle("#eta");
  hpteta->GetYaxis()->CenterTitle();
  hpteta->SetStats(kFALSE);
  hpteta->DrawCopy("colz");
  sprintf(buffer, "effQA/hpteta%d_%d.png",energy,PID);
  cpteta->Print(buffer);
  
  TCanvas *cptvz = new TCanvas("cptvz", "ptvz", 800, 600);
  hptvz->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptvz->GetXaxis()->CenterTitle();
  hptvz->GetXaxis()->SetRangeUser(0.2,4.5);
  hptvz->GetYaxis()->SetTitle("V_{Z}");
  hptvz->GetYaxis()->CenterTitle();
  hptvz->SetStats(kFALSE);
  hptvz->DrawCopy("colz");
  sprintf(buffer, "effQA/hptvz%d_%d.png",energy,PID);
  cptvz->Print(buffer);
  
  TCanvas *cetacent = new TCanvas("cetacent", "etacent", 800, 600);
  hetacent->GetXaxis()->SetTitle("#eta");
  hetacent->GetXaxis()->CenterTitle();
  hetacent->GetYaxis()->SetTitle("refMult");
  hetacent->GetYaxis()->CenterTitle();
  hetacent->SetStats(kFALSE);
  hetacent->DrawCopy("colz");
  sprintf(buffer, "effQA/hetacent%d_%d.png",energy,PID);
  cetacent->Print(buffer);
  
  TCanvas *cetavz = new TCanvas("cetavz", "etavz", 800, 600);
  hetavz->GetXaxis()->SetTitle("#eta");
  hetavz->GetXaxis()->CenterTitle();
  hetavz->GetYaxis()->SetTitle("V_{Z}");
  hetavz->GetYaxis()->CenterTitle();
  hetavz->SetStats(kFALSE);
  hetavz->DrawCopy("colz");
  sprintf(buffer, "effQA/hetavz%d_%d.png",energy,PID);
  cetavz->Print(buffer);
  
  TCanvas *ccentvz = new TCanvas("ccentvz", "centvz", 800, 600);
  hcentvz->GetXaxis()->SetTitle("refmult");
  hcentvz->GetXaxis()->CenterTitle();
  hcentvz->GetYaxis()->SetTitle("V_{Z}");
  hcentvz->GetYaxis()->CenterTitle();
  hcentvz->SetStats(kFALSE);
  hcentvz->DrawCopy("colz");
  sprintf(buffer, "effQA/hcentvz%d_%d.png",energy,PID);
  ccentvz->Print(buffer);
  
  TCanvas *cptphi = new TCanvas("cptphi", "ptphi", 800, 600);
  hptphi->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  hptphi->GetXaxis()->CenterTitle();
  hptphi->GetXaxis()->SetRangeUser(0.2,4.5);
  hptphi->GetYaxis()->SetTitle("#phi");
  hptphi->GetYaxis()->CenterTitle();
  hptphi->SetStats(kFALSE);
  hptphi->DrawCopy("colz");
  sprintf(buffer, "effQA/hptphi%d_%d.png",energy,PID);
  cptphi->Print(buffer);

  TCanvas *cDCA = new TCanvas("cDCA", "DCA", 800, 600);
  hDCA->GetXaxis()->SetTitle("DCA(cm)");
  hDCA->GetXaxis()->CenterTitle();
  hDCA->GetYaxis()->SetTitle("dN/dDCA(cm^{-1})");
  hDCA->GetYaxis()->CenterTitle();
  hDCA->SetStats(kFALSE);
  hDCA->DrawCopy("");
  sprintf(buffer, "effQA/hDCA%d_%d.png",energy,PID);
  cDCA->Print(buffer);

  TCanvas *cDCAphi = new TCanvas("cDCAphi", "DCAphi", 800, 600);
  hDCAphi->GetXaxis()->SetTitle("DCA(cm)");
  hDCAphi->GetXaxis()->CenterTitle();
  hDCAphi->GetYaxis()->SetTitle("#phi");
  hDCAphi->GetYaxis()->CenterTitle();
  hDCAphi->SetStats(kFALSE);
  hDCAphi->DrawCopy("colz");
  sprintf(buffer, "effQA/hDCAphi%d_%d.png",energy,PID);
  cDCAphi->Print(buffer);
  
  fOut->Write();
  fOut->Close();
}

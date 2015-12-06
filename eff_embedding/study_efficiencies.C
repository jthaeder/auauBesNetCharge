#include "TCanvas.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "math.h"
#include "TRandom3.h"
#include <iostream>

using namespace std;

double binomial(double *x, double *par) {
  //  P(k; p,n) -> scale k = nx  -> P (x/n ; p,n) = 
  //  n = par[0]
  //  p = par[1]
  //
  //  return TMath::Binomial(par[0], x[0]*par[0])*pow(par[1], par[0]*x[0])*pow(1. - par[1], par[0]*(1. - x[0]));

  // not scaled
  return TMath::Binomial(par[0], x[0])*pow(par[1], x[0])*pow(1. - par[1], par[0] - x[0]);
}

void study_efficiencies(const char* particle = "piplus", int energy = 14) {

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

  int idxPart = 0;

  // --------------------------------------------------------------------------------
  // -- track cuts
  // --------------------------------------------------------------------------------
  
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

  // --------------------------------------------------------------------------------
  // -- get file names
  // --------------------------------------------------------------------------------
  TString out_file, in_file;
  TString prefix("embeddingTrees");
  TString postfix("efficiencyStudy");

  gSystem->Exec("mkdir -p efficiencyStudy");

  int PID = 0;
  
  if(TString(particle) == TString("piplus")){
    in_file = prefix + TString("/SinglePiPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piplus") + NRG + TString("GeV.root");
    PID = 8;
    idxPart = 0;
  }
  else if(TString(particle) == TString("piminus")){
    in_file = prefix + TString("/SinglePiMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piminus") + NRG + TString("GeV.root");
    PID = 9;
    idxPart = 0;
  }
  else if(TString(particle) == TString("kaonplus")){
    in_file = prefix + TString("/SingleKPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kplus") + NRG + TString("GeV.root");
    PID = 11;
    idxPart = 1;
  }
  else if(TString(particle) == TString("kaonminus")){
    in_file = prefix + TString("/SingleKMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kminus") + NRG + TString("GeV.root");
    PID = 12;
    idxPart = 1;
  }
  else if(TString(particle) == TString("protonplus")){
    in_file = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pplus") + NRG + TString("GeV.root");
    PID = 14;
    idxPart = 2;
  }
  else if(TString(particle) == TString("protonminus")){
    in_file = prefix + TString("/SinglePbarNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pminus") + NRG + TString("GeV.root");
    PID = 15;
    idxPart = 2;
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
  float AllPts, NPossible, ParentGeantId, GeantId,mErrP, NCommonHit, EventId, RunId;
  MatchedPairs->SetBranchAddress("EventId", &EventId);
  MatchedPairs->SetBranchAddress("RunId", &RunId);
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
  float pRefMult, pRefMultCorrected, pCentralityWeight, pCentrality16, pVertexX, pVertexY, pVertexZ, pPtMc, pPzMc, pEtaMc, pPhiMc, pParentGeantId, pGeantId, pEventId, pRunId;
  McTrack->SetBranchAddress("EventId", &pEventId);
  McTrack->SetBranchAddress("RunId", &pRunId);
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

  const int   nCent            = 9;
  const char* cent[]           = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};

  const int    Nbins           = 20;
  const double pt_bin[21]      = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.4,1.7,2.0,2.3,2.7,3.5,4.0,4.5,5.0};

  const int    NbinsRed        = 4;
  const double pt_binRed[3][5] = { {0.0, 0.5, 0.7, 2.3, 5.0},
				   {0.0, 0.5, 1.0, 2.3, 5.0},
				   {0.0, 0.5, 0.8, 2.3, 5.0} };

  const double smearBin        = 0.3;
    
  const int    max_events      = 100;

  const int    nTracksMcMax    = 25;

  // --------------------------------------------------------------------------------
  // -- Initialize pt hists
  TProfile *effProfile[nCent];
  TProfile *effProfileSmeared[nCent];

  TProfile *effProfileRed[nCent];
  TProfile *effProfileSmearedRed[nCent];

  TH1D* hpt_mc[nCent];
  TH1D* hpt_rec[nCent];

  TH1D* hnTracks_mc[nCent];
  TH1D* hnTracks_mc_sum = new TH1D("hnTracks_mc_sum", ";N_{Mc}; Events", nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);

  TH1D* hnTracks_rec[nCent];
  TH1D* hnTracks_rec_sum = new TH1D("hnTracks_rec_sum", ";N_{Rec}; Events", nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  
  TH1D* hnTracksRec[nCent][nTracksMcMax];
  TH1D* hnTracksRec_sum[nTracksMcMax];

  TF1* fnTracksRec[nCent][nTracksMcMax];
  TF1* fnTracksRec_sum[nTracksMcMax];

  for (Int_t idxTrack = 0; idxTrack < nTracksMcMax; idxTrack++) {
    hnTracksRec_sum[idxTrack] = new TH1D(Form("hnTracksRec_sum_%d", idxTrack), ";N_{Rec};Events",
					      nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);

      fnTracksRec_sum[idxTrack] = new TF1(Form("fnTracksRec_sum_%d", idxTrack), binomial, 0, nTracksMcMax, 2);
    
    for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
      hnTracksRec[idxCent][idxTrack] = new TH1D(Form("hnTracksRec_%s_%d", cent[idxCent], idxTrack), ";N_{Rec};Events",
						nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);
      fnTracksRec[idxCent][idxTrack] = new TF1(Form("fnTracksRec_%s_%d", cent[idxCent], idxTrack), binomial, 0, nTracksMcMax, 2);

    }
  }


  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  TH2D* hnTracks_mc_rec[nCent];
  TH2D* hnTracks_mc_rec_sum = new TH2D("hnTracks_mc_rec_sum", ";N_{Mc};N_{Rec}", 
				       nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5, nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);


  TH1D* hpt_width[nCent];
  TH1D* hpt_widthSmeared[nCent];
  TH1D* hpt_relativeWidth[nCent];

  TH2D* hpt_2D[nCent];

  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    effProfile[idxCent] = new TProfile(Form("effProfile_%s", cent[idxCent]), 
				       Form("effProfile_%s;#it{p}_{T} (GeV/#it{c});efficiency", cent[idxCent]), Nbins, 0., 5.);
    effProfile[idxCent]->GetXaxis()->Set(Nbins,pt_bin);

    effProfileSmeared[idxCent] = new TProfile(Form("effProfileSmeared_%s", cent[idxCent]), 
				       Form("effProfileSmeared_%s;#it{p}_{T} (GeV/#it{c});efficiency", cent[idxCent]), Nbins, 0., 5.);
    effProfileSmeared[idxCent]->GetXaxis()->Set(Nbins,pt_bin);

    effProfileRed[idxCent] = new TProfile(Form("effProfileRed_%s", cent[idxCent]), 
					  Form("effProfileRed_%s;#it{p}_{T} (GeV/#it{c});efficiency", cent[idxCent]), NbinsRed, 0., 5.);
    effProfileRed[idxCent]->GetXaxis()->Set(NbinsRed,pt_binRed[idxPart]);

    effProfileSmearedRed[idxCent] = new TProfile(Form("effProfileSmearedRed_%s", cent[idxCent]), 
						 Form("effProfileSmearedRed_%s;#it{p}_{T} (GeV/#it{c});efficiency", cent[idxCent]), NbinsRed, 0., 5.);
    effProfileSmearedRed[idxCent]->GetXaxis()->Set(NbinsRed,pt_binRed[idxPart]);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

    hpt_mc[idxCent] = new TH1D(Form("hpt_mc_%s",  cent[idxCent]), "", Nbins, 0., 5.);
    hpt_mc[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_rec[idxCent] = new TH1D(Form("hpt_rec_%s", cent[idxCent]), "", Nbins, 0., 5.);
    hpt_rec[idxCent]->GetXaxis()->Set(Nbins, pt_bin);


    hnTracks_mc[idxCent]  = new TH1D(Form("hnTracks_mc_%s",  cent[idxCent]), ";N_{Mc}; Events",  nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);
    hnTracks_rec[idxCent] = new TH1D(Form("hnTracks_rec_%s", cent[idxCent]), ";N_{Rec}; Events", nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);

    hnTracks_mc_rec[idxCent] = new TH2D(Form("hnTracks_mc_rec_%s",  cent[idxCent]), ";N_{Mc};N_{Rec}", 
					nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5, nTracksMcMax+1, -0.5, Double_t(nTracksMcMax)+0.5);

    

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

    hpt_width[idxCent] = new TH1D(Form("hpt_width_%s", cent[idxCent]), "Width;#it{p}_{T} (GeV/#it{c});width", Nbins, 0., 5.);
    hpt_width[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_widthSmeared[idxCent] = new TH1D(Form("hpt_widthSmeared_%s", cent[idxCent]), "Width;#it{p}_{T} (GeV/#it{c});width", Nbins, 0., 5.);
    hpt_widthSmeared[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_relativeWidth[idxCent] = new TH1D(Form("hpt_relativeWidth_%s", cent[idxCent]), "Relative Width;#it{p}_{T} (GeV/#it{c});width/mean", Nbins, 0., 5.);
    hpt_relativeWidth[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_2D[idxCent] = new TH2D(Form("hpt_2D_%s", cent[idxCent]), "distribution;#it{p}_{T} (GeV/#it{c});efficiency", Nbins, 0., 5., 40, 0., 1.);
    hpt_2D[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
  }

  // --------------------------------------------------------------------------------
  // -- Prepare syncing - get nEvents for MC and Rec
  // --------------------------------------------------------------------------------
  int nTracksMC  = McTrack->GetEntries();
  int nTracksRec = MatchedPairs->GetEntries();

  int lastEvent = -1;
  int nEventsMC  = 0;
  //  int nEventsRec = 0;

  for(int i = 0; i < nTracksMC; i++){
    McTrack->GetEntry(i);
    int current=Int_t(pEventId);
    
    if (current != lastEvent) {
      lastEvent = current;
      ++nEventsMC;
    }
  }
  
  // lastEvent = -1;
  // for(int i = 0; i < nTracksRec; i++){
  //   MatchedPairs->GetEntry(i);
  //   int current=Int_t(EventId);
    
  //   if (current != lastEvent) {
  //     lastEvent = current;
  //     ++nEventsRec;
  //   }
  // }

  printf("nEvents MC %d \n",  nEventsMC);
  //  printf("nEvents Rec %d \n", nEventsRec);

  // --------------------------------------------------------------------------------
  // -- Loop over events 
  //    - and associate tracks to them
  // --------------------------------------------------------------------------------
  int idxMC = 0;
  McTrack->GetEntry(0);
  int currentEventMC  = Int_t(pEventId); 
  int lastEventMC     = 0;

  int idxRec = 0;
  MatchedPairs->GetEntry(0);
  int currentEventRec = Int_t(EventId);
  int lastEventRec    = 0;

  int eventCounter[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++)
    eventCounter[idxCent] = 0;
  int nEventsCommon  = 0;

  // -----------------------------------------------

  int lastRun = -1;
  int nRuns   = 0;

  // -----------------------------------------------


  // -- event loop
  for (int idxEvent = 0; idxEvent < nEventsMC; ++idxEvent) {

    int nTracksEventMc    = 0;
    int nTracksEventRec    = 0;

    // -- MC loop
    // --------------------------------------------------------------------------------
    lastEventMC = currentEventMC;
    for (; idxMC < nTracksMC; idxMC++) {
      McTrack->GetEntry(idxMC);
      currentEventMC = Int_t(pEventId);

      if (currentEventMC != lastEventMC) 
	break;

      // --------------------------------------------------------------------------------
      // -- analyze current MC event
      //    - get Centrality
      //    - check cuts
      //    - fill histogram
      // --------------------------------------------------------------------------------
     
      int idxCent = -1;
      if (pCentrality16 == 15)      idxCent = 0;
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
      
      if (pParentGeantId != 0)
	continue;

      if (pGeantId != PID)
	continue;
      
      // -- within Vertex R
      if ((pVertexX-vxo)*(pVertexX-vxo) + (pVertexY-vyo)*(pVertexY-vyo) > max_r*max_r)
	continue;
      
      // -- within Vertex Z
      if (pVertexZ > max_z || pVertexZ < -max_z) 
	continue;
      
      // -- refMult usefull
      if (pRefMultCorrected < 0 || pRefMultCorrected > 1000)
	continue;
      
      // -- with eta window
      if (fabs(pEtaMc) > max_eta)
	continue;
      // --------------------------------------------------------------------------------
      ++nTracksEventMc;
      hpt_mc[idxCent]->Fill(pPtMc);
    } // for(; idxMC < nTracksMC; idxMC++) {
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

    // -- Rec loop
    // --------------------------------------------------------------------------------
    lastEventRec = currentEventRec;
    for (; idxRec < nTracksRec; idxRec++) {
      MatchedPairs->GetEntry(idxRec);
      currentEventRec = Int_t(EventId);

      if (currentEventRec != lastEventRec) 
	break;

      if (currentEventRec != lastEventMC)
	break;

      // --------------------------------------------------------------------------------
      // -- analyze current Rec event
      //    - get Centrality
      //    - check cuts
      //    - fill histogram
      // --------------------------------------------------------------------------------
   
      int idxCent = -1;
      if (Centrality16 == 15)      idxCent = 0;
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

      if (ParentGeantId != 0)
	continue;
      
      if (GeantId != PID)
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
      
      if (VertexZ > max_z || VertexZ < -max_z) 
	continue;
      
      if (fabs(EtaMc) > max_eta)
	continue;

      // --------------------------------------------------------------------------------
      ++nTracksEventRec;
      hpt_rec[idxCent]->Fill(PtMc);
    } // for(; idxRec < nTracksRec; idxRec++) {

    // -----------------------------------------------
    // -- get centrality
    int idxCent = -1;
    if (pCentrality16 == 15)      idxCent = 0;
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

    // -----------------------------------------------

    hnTracks_mc_sum->Fill(nTracksEventMc);
    hnTracks_mc[idxCent]->Fill(nTracksEventMc);

    hnTracks_rec_sum->Fill(nTracksEventRec);
    hnTracks_rec[idxCent]->Fill(nTracksEventRec);

    hnTracks_mc_rec_sum->Fill(nTracksEventMc, nTracksEventRec);
    hnTracks_mc_rec[idxCent]->Fill(nTracksEventMc, nTracksEventRec);


    hnTracksRec[idxCent][nTracksEventMc]->Fill(nTracksEventRec);
    hnTracksRec_sum[nTracksEventMc]->Fill(nTracksEventRec);
    
    // hnTracksPtRec[idxCent][nTracksEventMc]->Fill(nTracksEventRec);
    // hnTracksPtRec_sum[nTracksEventMc]->Fill(nTracksEventRec);
    


    // -----------------------------------------------

    ++nEventsCommon;
    ++(eventCounter[idxCent]);

    // -----------------------------------------------

    int currentRun = Int_t(pRunId);
    if (currentRun != lastRun) {
      lastRun = currentRun;
      ++nRuns;
    }

    //  if (nRuns > 20) {
      // -- pT dependent average over max_events
    if (eventCounter[idxCent] >= max_events) {
      for (int idx = 1; idx <= hpt_mc[idxCent]->GetXaxis()->GetNbins(); idx++) { 
	
	double binMC  = hpt_mc[idxCent]->GetBinContent(idx);
	double binRec = hpt_rec[idxCent]->GetBinContent(idx);

	if (binMC < 0.00001)
	  continue;

	double ratio = binRec / binMC;
	effProfile[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx), ratio);
	effProfileRed[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx), ratio);

	binMC  += ((2*smearBin*gRandom->Rndm()) - smearBin);
	binRec += ((2*smearBin*gRandom->Rndm()) - smearBin);
	ratio = binRec / binMC;
	effProfileSmeared[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx), ratio);
	effProfileSmearedRed[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx), ratio);
	
	hpt_2D[idxCent]->Fill(           hpt_mc[idxCent]->GetBinCenter(idx), ratio);
      }
      
      hpt_mc[idxCent]->Reset();
      hpt_rec[idxCent]->Reset();
      
      nRuns = 0;
      eventCounter[idxCent] = 0;
    } // if (eventCounter[idxCent] >= max_events) {
  } //  for (int idxEvent = 0 ; idxEvent < nEventsMC; ++idxEvent) {

  for (int idxCent = 0; idxCent < nCent; ++idxCent) 
    for (int idx = 1; idx <= hpt_mc[idxCent]->GetXaxis()->GetNbins(); idx++) {
      effProfile[idxCent]->SetBinError(idx, effProfile[idxCent]->GetRMS(idx));
      effProfileSmeared[idxCent]->SetBinError(idx, effProfileSmeared[idxCent]->GetRMS(idx));
    }

  for (int idxCent = 0; idxCent < nCent; ++idxCent) 
    for (int idx = 1; idx <= effProfileRed[idxCent]->GetXaxis()->GetNbins(); idx++) {
      effProfileRed[idxCent]->SetBinError(idx, effProfileRed[idxCent]->GetRMS(idx));
      effProfileSmearedRed[idxCent]->SetBinError(idx, effProfileSmearedRed[idxCent]->GetRMS(idx));
    }


  // --------------------------------------------------------------------------------

  printf("nEvents Common %d \n",  nEventsCommon);

  // --------------------------------------------------------------------------------

  for (int idxCent = 0; idxCent < nCent; ++idxCent) {
    for (int idx = 1; idx <= effProfile[idxCent]->GetXaxis()->GetNbins(); idx++) { 
      hpt_width[idxCent]->SetBinContent(idx, effProfile[idxCent]->GetBinError(idx));
      hpt_widthSmeared[idxCent]->SetBinContent(idx, effProfileSmeared[idxCent]->GetBinError(idx));

      double relativeWidth = (effProfile[idxCent]->GetBinContent(idx) != 0) ? 
	effProfile[idxCent]->GetBinError(idx)/effProfile[idxCent]->GetBinContent(idx) : 0.;
      hpt_relativeWidth[idxCent]->SetBinContent(idx, relativeWidth);
    }
  }

  // --------------------------------------------------------------------------------

  for (int idxCent = 0; idxCent < nCent; ++idxCent) {
    effProfile[idxCent]->GetYaxis()->SetRangeUser(0., 1.0);
    effProfile[idxCent]->SetMarkerStyle(21);
    effProfile[idxCent]->SetMarkerColor(kRed+2);
    effProfile[idxCent]->SetMarkerColor(kRed+2);

    effProfileSmeared[idxCent]->GetYaxis()->SetRangeUser(0., 1.0);
    effProfileSmeared[idxCent]->SetMarkerStyle(21);
    effProfileSmeared[idxCent]->SetMarkerColor(kRed+2);
    effProfileSmeared[idxCent]->SetMarkerColor(kRed+2);

    effProfileRed[idxCent]->GetYaxis()->SetRangeUser(0., 1.0);
    effProfileRed[idxCent]->SetMarkerStyle(21);
    effProfileRed[idxCent]->SetMarkerColor(kRed+2);
    effProfileRed[idxCent]->SetMarkerColor(kRed+2);

    effProfileSmearedRed[idxCent]->GetYaxis()->SetRangeUser(0., 1.0);
    effProfileSmearedRed[idxCent]->SetMarkerStyle(21);
    effProfileSmearedRed[idxCent]->SetMarkerColor(kRed+2);
    effProfileSmearedRed[idxCent]->SetMarkerColor(kRed+2);

    hpt_width[idxCent]->GetYaxis()->SetRangeUser(0., 0.2);
    hpt_width[idxCent]->SetMarkerStyle(20);
    hpt_width[idxCent]->SetMarkerColor(kRed+2);
    hpt_width[idxCent]->SetLineColor(kRed+2);

    hpt_widthSmeared[idxCent]->GetYaxis()->SetRangeUser(0., 0.2);
    hpt_widthSmeared[idxCent]->SetMarkerStyle(20);
    hpt_widthSmeared[idxCent]->SetMarkerColor(kRed+2);
    hpt_widthSmeared[idxCent]->SetLineColor(kRed+2);

    hpt_relativeWidth[idxCent]->GetYaxis()->SetRangeUser(0., 0.2);
    hpt_relativeWidth[idxCent]->SetMarkerStyle(21);
    hpt_relativeWidth[idxCent]->SetMarkerColor(kAzure);
    hpt_relativeWidth[idxCent]->SetLineColor(kAzure);
  }
  
  // --------------------------------------------------------------------------------

  for (Int_t idxTrack = 1; idxTrack < nTracksMcMax; idxTrack++) {
    Double_t nEntries = hnTracksRec_sum[idxTrack]->GetEntries();
    if (nEntries == 0){
      delete hnTracksRec_sum[idxTrack];
      hnTracksRec_sum[idxTrack] = NULL;
      continue;
    }

    Double_t mean = hnTracksRec_sum[idxTrack]->GetMean();
    Double_t eff  = mean / Double_t(idxTrack);

    hnTracksRec_sum[idxTrack]->Scale(1./hnTracksRec_sum[idxTrack]->Integral());

    fnTracksRec_sum[idxTrack]->SetParameter(0, nEntries);
    fnTracksRec_sum[idxTrack]->SetParameter(1, eff);
    // int flag = hnTracksRec_sum[idxTrack]->Fit(fnTracksRec_sum[idxTrack], "QRsame", "", 0.0, Double_t(idxTrack));

    fnTracksRec_sum[idxTrack]->SetLineColor(kRed+2);
    fnTracksRec_sum[idxTrack]->SetMarkerColor(kRed+2);
    fnTracksRec_sum[idxTrack]->SetMarkerStyle(20);
    fnTracksRec_sum[idxTrack]->Write();

    hnTracksRec_sum[idxTrack]->SetLineColor(kAzure);
    hnTracksRec_sum[idxTrack]->SetMarkerColor(kAzure);
    hnTracksRec_sum[idxTrack]->SetMarkerStyle(25);
  }

  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    delete hnTracksRec[idxCent][0];
    hnTracksRec[idxCent][0] = NULL;
  }

  for (Int_t idxTrack = 1; idxTrack < nTracksMcMax; idxTrack++) {
    for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
      Double_t nEntries = hnTracksRec[idxCent][idxTrack]->GetEntries();

      if (nEntries == 0){
	delete hnTracksRec[idxCent][idxTrack];
	hnTracksRec[idxCent][idxTrack] = NULL;
	continue;
      }

      Double_t mean = hnTracksRec[idxCent][idxTrack]->GetMean();
      Double_t eff  = mean / Double_t(idxTrack);
      hnTracksRec[idxCent][idxTrack]->Scale(1./hnTracksRec[idxCent][idxTrack]->Integral());

      fnTracksRec[idxCent][idxTrack]->SetParameter(0, nEntries);
      fnTracksRec[idxCent][idxTrack]->SetParameter(1, eff);
      // int flag = hnTracksRec[idxCent][idxTrack]->Fit(fnTracksRec[idxCent][idxTrack], "QRsame", "", 0.10, 1.0);
      
      fnTracksRec[idxCent][idxTrack]->SetLineColor(kRed+2);
      fnTracksRec[idxCent][idxTrack]->SetMarkerColor(kRed+2);
      fnTracksRec[idxCent][idxTrack]->SetMarkerStyle(20);
      fnTracksRec[idxCent][idxTrack]->Write();
      
      hnTracksRec[idxCent][idxTrack]->SetLineColor(kAzure);
      hnTracksRec[idxCent][idxTrack]->SetMarkerColor(kAzure);
      hnTracksRec[idxCent][idxTrack]->SetMarkerStyle(25);
    }
  }

  // --------------------------------------------------------------------------------

  fOut->Write();
  fOut->Close();
}

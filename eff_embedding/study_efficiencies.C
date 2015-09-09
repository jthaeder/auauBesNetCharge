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
#include "TRandom3.h"
#include <iostream>

using namespace std;

void study_efficiencies(const char* particle, int energy = 14) {

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
  
  float min_dedx         = 5;  // default - was 10 , try not to use  10
  float min_nCommon      = 10;
  float min_fitpts       = 20;//15 // 10
  float min_fitpts_nposs = .52;
  float max_eta          = 0.5;
  float max_dca          = 1.0;
  
  int   max_events       = 250;
  double smearBin        = 0.3;
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
  }
  else if(TString(particle) == TString("piminus")){
    in_file = prefix + TString("/SinglePiMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/piminus") + NRG + TString("GeV.root");
    PID = 9;
  }
  else if(TString(particle) == TString("kaonplus")){
    in_file = prefix + TString("/SingleKPlusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kplus") + NRG + TString("GeV.root");
    PID = 11;
  }
  else if(TString(particle) == TString("kaonminus")){
    in_file = prefix + TString("/SingleKMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kminus") + NRG + TString("GeV.root");
    PID = 12;
  }
  else if(TString(particle) == TString("protonplus")){
    in_file = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pplus") + NRG + TString("GeV.root");
    PID = 14;
  }
  else if(TString(particle) == TString("protonminus")){
    in_file = prefix + TString("/SinglePbarNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pminus") + NRG + TString("GeV.root");
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

  const int   nCent       = 9;
  const char* cent[]      = {"0005","0510","1020","2030","3040","4050","5060","6070","7080"};

  const int    Nbins      = 20;
  const double pt_bin[21] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.4,1.7,2.0,2.3,2.7,3.5,4.0,4.5,5.0};

  // --------------------------------------------------------------------------------
  // -- Initialize pt hists
  TProfile *effProfile[nCent];
  TProfile *effProfileSmeared[nCent];

  TH1D* hpt_mc[nCent];
  TH1D* hpt_rec[nCent];

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

    hpt_mc[idxCent] = new TH1D(Form("hpt_mc_%s",  cent[idxCent]), "", Nbins, 0., 5.);
    hpt_mc[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_rec[idxCent] = new TH1D(Form("hpt_rec_%s", cent[idxCent]), "", Nbins, 0., 5.);
    hpt_rec[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_width[idxCent] = new TH1D(Form("hpt_width_%s", cent[idxCent]), "Width;#it{p}_{T} (GeV/#it{c});width", Nbins, 0., 5.);
    hpt_width[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_widthSmeared[idxCent] = new TH1D(Form("hpt_widthSmeared_%s", cent[idxCent]), "Width;#it{p}_{T} (GeV/#it{c});width", Nbins, 0., 5.);
    hpt_widthSmeared[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_relativeWidth[idxCent] = new TH1D(Form("hpt_relativeWidth_%s", cent[idxCent]), "Relative Width;#it{p}_{T} (GeV/#it{c});width/mean", Nbins, 0., 5.);
    hpt_relativeWidth[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_2D[idxCent] = new TH2D(Form("hpt_2D_%s", cent[idxCent]), "distribution;#it{p}_{T} (GeV/#it{c});efficiency", Nbins, 0., 5., 101, 0., 1.);
    hpt_2D[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
  }

  TH2D* hpt_2DAll= new TH2D(Form("hpt_2DAll"), "distribution;#it{p}_{T} (GeV/#it{c});efficiency", Nbins, 0., 5., 101, 0., 1.);
  hpt_2DAll->GetXaxis()->Set(Nbins, pt_bin);

  // --------------------------------------------------------------------------------
  // -- Prepare syncing - get nEvents for MC and Rec
  // --------------------------------------------------------------------------------
  int nTracksMC  = McTrack->GetEntries();
  int nTracksRec = MatchedPairs->GetEntries();

  int lastEvent = -1;
  int nEventsMC  = 0;
  int nEventsRec = 0;

  for(int i = 0; i < nTracksMC; i++){
    McTrack->GetEntry(i);
    int current=Int_t(pEventId);
    
    if (current != lastEvent) {
      lastEvent = current;
      ++nEventsMC;
    }
  }
  
  lastEvent = -1;
  for(int i = 0; i < nTracksRec; i++){
    MatchedPairs->GetEntry(i);
    int current=Int_t(EventId);
    
    if (current != lastEvent) {
      lastEvent = current;
      ++nEventsRec;
    }
  }

  printf("nEvents MC %d \n",  nEventsMC);
  printf("nEvents Rec %d \n", nEventsRec);

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

	binMC  += (gRandom->Rndm() - smearBin);
	binRec += (gRandom->Rndm() - smearBin);
	ratio = binRec / binMC;
	effProfileSmeared[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx), ratio);
	
	hpt_2DAll->Fill(hpt_mc[idxCent]->GetBinCenter(idx),           ratio);
	hpt_2D[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx),     ratio);
      }
      
      hpt_mc[idxCent]->Reset();
      hpt_rec[idxCent]->Reset();
      
      nRuns = 0;
      eventCounter[idxCent] = 0;
    } // if (eventCounter[idxCent] >= max_events) {
  } //  for (int idxEvent = 0 ; idxEvent < nEventsMC; ++idxEvent) {

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
  
  fOut->Write();
  fOut->Close();
}

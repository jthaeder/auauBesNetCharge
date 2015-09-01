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

void study_efficiencies2(const char* particle, int energy = 14) {

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
  TString out_file, in_file1, in_file2;
  TString prefix("embeddingTrees");
  TString postfix("efficiencyStudy2");

  gSystem->Exec("mkdir -p efficiencyStudy2");

  int PID = 0;
  
  if(TString(particle) == TString("piplus") || TString(particle) == TString("piminus")){
    in_file1 = prefix + TString("/SinglePiPlusNT_Embed_") + NRG + TString("GeV.root");
    in_file2 = prefix + TString("/SinglePiMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pi") + NRG + TString("GeV.root");
    PID = 8;
  }
  else if(TString(particle) == TString("kaonplus") || TString(particle) == TString("kaonminus")){
    in_file1 = prefix + TString("/SingleKPlusNT_Embed_") + NRG + TString("GeV.root");
    in_file2 = prefix + TString("/SingleKMinusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/kplus") + NRG + TString("GeV.root");
    PID = 11;
  }
  else if(TString(particle) == TString("protonplus") || TString(particle) == TString("protonminus")){
    in_file1 = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    in_file2 = prefix + TString("/SinglePbarNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/pplus") + NRG + TString("GeV.root");
    PID = 14;
  }

  // --------------------------------------------------------------------------------
  // -- get input file
  // --------------------------------------------------------------------------------
  cout << "opening input files: " << in_file1 << "  " << in_file2<< endl;
  
  TFile *fIn1 = TFile::Open(in_file1, "READ");
  TFile *fIn2 = TFile::Open(in_file2, "READ");
  TNtuple *MatchedPairs1 = static_cast<TNtuple*>(fIn1->Get("MatchedPairs_NT"));
  TNtuple *MatchedPairs2 = static_cast<TNtuple*>(fIn2->Get("MatchedPairs_NT"));
  TNtuple *McTrack1      = static_cast<TNtuple*>(fIn1->Get("McTrack_NT"));
  TNtuple *McTrack2      = static_cast<TNtuple*>(fIn2->Get("McTrack_NT"));

  // --------------------------------------------------------------------------------
  // -- get matched ntuple 1
  float Dedx1, RefMult1, RefMultCorrected1, CentralityWeight1, Centrality161, VertexX1, VertexY1, VertexZ1;
  float PtMc1, PzMc1, EtaMc1, PhiMc1, PtPr1, PtGl1, EtaPr1, PhiPr1, DcaGl1, DcaZGl1, DcaXYGl1, Flag1, FitPts1, DedxPts1;
  float AllPts1, NPossible1, ParentGeantId1, GeantId1, mErrP1, NCommonHit1, EventId1;
  MatchedPairs1->SetBranchAddress("EventId", &EventId1);
  MatchedPairs1->SetBranchAddress("Dedx", &Dedx1);
  MatchedPairs1->SetBranchAddress("RefMult", &RefMult1);
  MatchedPairs1->SetBranchAddress("RefMultCorrected", &RefMultCorrected1);
  MatchedPairs1->SetBranchAddress("CentralityWeight", &CentralityWeight1);
  MatchedPairs1->SetBranchAddress("Centrality16", &Centrality161);
  MatchedPairs1->SetBranchAddress("VertexX", &VertexX1);
  MatchedPairs1->SetBranchAddress("VertexY", &VertexY1);
  MatchedPairs1->SetBranchAddress("VertexZ", &VertexZ1);
  MatchedPairs1->SetBranchAddress("PtMc", &PtMc1);
  MatchedPairs1->SetBranchAddress("PzMc", &PzMc1);
  MatchedPairs1->SetBranchAddress("EtaMc", &EtaMc1);
  MatchedPairs1->SetBranchAddress("PhiMc", &PhiMc1);
  MatchedPairs1->SetBranchAddress("PtPr", &PtPr1);
  MatchedPairs1->SetBranchAddress("mErrP", &mErrP1);
  MatchedPairs1->SetBranchAddress("PtGl", &PtGl1);
  MatchedPairs1->SetBranchAddress("EtaPr", &EtaPr1);
  MatchedPairs1->SetBranchAddress("PhiPr", &PhiPr1);
  MatchedPairs1->SetBranchAddress("DcaGl", &DcaGl1);
  MatchedPairs1->SetBranchAddress("DcaZGl", &DcaZGl1);
  MatchedPairs1->SetBranchAddress("DcaXYGl", &DcaXYGl1);
  MatchedPairs1->SetBranchAddress("Flag", &Flag1);
  MatchedPairs1->SetBranchAddress("FitPts", &FitPts1);
  MatchedPairs1->SetBranchAddress("DedxPts", &DedxPts1);
  MatchedPairs1->SetBranchAddress("AllPts", &AllPts1);
  MatchedPairs1->SetBranchAddress("NPossible", &NPossible1);
  MatchedPairs1->SetBranchAddress("ParentGeantId", &ParentGeantId1);
  MatchedPairs1->SetBranchAddress("GeantId", &GeantId1);
  MatchedPairs1->SetBranchAddress("NCommonHit", &NCommonHit1);

  // --------------------------------------------------------------------------------
  // -- get matched ntuple 2
  float Dedx2, RefMult2, RefMultCorrected2, CentralityWeight2, Centrality162, VertexX2, VertexY2, VertexZ2;
  float PtMc2, PzMc2, EtaMc2, PhiMc2, PtPr2, PtGl2, EtaPr2, PhiPr2, DcaGl2, DcaZGl2, DcaXYGl2, Flag2, FitPts2, DedxPts2;
  float AllPts2, NPossible2, ParentGeantId2, GeantId2, mErrP2, NCommonHit2, EventId2;
  MatchedPairs2->SetBranchAddress("EventId", &EventId2);
  MatchedPairs2->SetBranchAddress("Dedx", &Dedx2);
  MatchedPairs2->SetBranchAddress("RefMult", &RefMult2);
  MatchedPairs2->SetBranchAddress("RefMultCorrected", &RefMultCorrected2);
  MatchedPairs2->SetBranchAddress("CentralityWeight", &CentralityWeight2);
  MatchedPairs2->SetBranchAddress("Centrality16", &Centrality162);
  MatchedPairs2->SetBranchAddress("VertexX", &VertexX2);
  MatchedPairs2->SetBranchAddress("VertexY", &VertexY2);
  MatchedPairs2->SetBranchAddress("VertexZ", &VertexZ2);
  MatchedPairs2->SetBranchAddress("PtMc", &PtMc2);
  MatchedPairs2->SetBranchAddress("PzMc", &PzMc2);
  MatchedPairs2->SetBranchAddress("EtaMc", &EtaMc2);
  MatchedPairs2->SetBranchAddress("PhiMc", &PhiMc2);
  MatchedPairs2->SetBranchAddress("PtPr", &PtPr2);
  MatchedPairs2->SetBranchAddress("mErrP", &mErrP2);
  MatchedPairs2->SetBranchAddress("PtGl", &PtGl2);
  MatchedPairs2->SetBranchAddress("EtaPr", &EtaPr2);
  MatchedPairs2->SetBranchAddress("PhiPr", &PhiPr2);
  MatchedPairs2->SetBranchAddress("DcaGl", &DcaGl2);
  MatchedPairs2->SetBranchAddress("DcaZGl", &DcaZGl2);
  MatchedPairs2->SetBranchAddress("DcaXYGl", &DcaXYGl2);
  MatchedPairs2->SetBranchAddress("Flag", &Flag2);
  MatchedPairs2->SetBranchAddress("FitPts", &FitPts2);
  MatchedPairs2->SetBranchAddress("DedxPts", &DedxPts2);
  MatchedPairs2->SetBranchAddress("AllPts", &AllPts2);
  MatchedPairs2->SetBranchAddress("NPossible", &NPossible2);
  MatchedPairs2->SetBranchAddress("ParentGeantId", &ParentGeantId2);
  MatchedPairs2->SetBranchAddress("GeantId", &GeantId2);
  MatchedPairs2->SetBranchAddress("NCommonHit", &NCommonHit2);

  // --------------------------------------------------------------------------------
  // -- get MC ntuple 1 
  float pRefMult1, pRefMultCorrected1, pCentralityWeight1, pCentrality161, pVertexX1, pVertexY1, pVertexZ1, pPtMc1, pPzMc1, pEtaMc1, pPhiMc1, pParentGeantId1, pGeantId1, pEventId1;
  McTrack1->SetBranchAddress("EventId", &pEventId1);
  McTrack1->SetBranchAddress("RefMult", &pRefMult1);
  McTrack1->SetBranchAddress("RefMultCorrected", &pRefMultCorrected1);
  McTrack1->SetBranchAddress("CentralityWeight", &pCentralityWeight1);
  McTrack1->SetBranchAddress("Centrality16", &pCentrality161);
  McTrack1->SetBranchAddress("VertexX", &pVertexX1);
  McTrack1->SetBranchAddress("VertexY", &pVertexY1);
  McTrack1->SetBranchAddress("VertexZ", &pVertexZ1);
  McTrack1->SetBranchAddress("PtMc", &pPtMc1);
  McTrack1->SetBranchAddress("PzMc", &pPzMc1);
  McTrack1->SetBranchAddress("EtaMc", &pEtaMc1);
  McTrack1->SetBranchAddress("PhiMc", &pPhiMc1);
  McTrack1->SetBranchAddress("ParentGeantId", &pParentGeantId1);
  McTrack1->SetBranchAddress("GeantId", &pGeantId1);
  
  // --------------------------------------------------------------------------------
  // -- get MC ntuple 2 
  float pRefMult2, pRefMultCorrected2, pCentralityWeight2, pCentrality162, pVertexX2, pVertexY2, pVertexZ2, pPtMc2, pPzMc2, pEtaMc2, pPhiMc2, pParentGeantId2, pGeantId2, pEventId2;
  McTrack2->SetBranchAddress("EventId", &pEventId2);
  McTrack2->SetBranchAddress("RefMult", &pRefMult2);
  McTrack2->SetBranchAddress("RefMultCorrected", &pRefMultCorrected2);
  McTrack2->SetBranchAddress("CentralityWeight", &pCentralityWeight2);
  McTrack2->SetBranchAddress("Centrality16", &pCentrality162);
  McTrack2->SetBranchAddress("VertexX", &pVertexX2);
  McTrack2->SetBranchAddress("VertexY", &pVertexY2);
  McTrack2->SetBranchAddress("VertexZ", &pVertexZ2);
  McTrack2->SetBranchAddress("PtMc", &pPtMc2);
  McTrack2->SetBranchAddress("PzMc", &pPzMc2);
  McTrack2->SetBranchAddress("EtaMc", &pEtaMc2);
  McTrack2->SetBranchAddress("PhiMc", &pPhiMc2);
  McTrack2->SetBranchAddress("ParentGeantId", &pParentGeantId2);
  McTrack2->SetBranchAddress("GeantId", &pGeantId2);
  
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
  TProfile2D *effProfile2D[nCent];
  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    effProfile2D[idxCent] = new TProfile(Form("effProfile2D_%s", cent[idxCent]), 
					 Form("effProfile2D_%s", cent[idxCent]), Nbins, 0., 5., Nbins, 0., 5.);
    effProfile2D[idxCent]->GetXaxis()->Set(Nbins,pt_bin);
    effProfile2D[idxCent]->GetYaxis()->Set(Nbins,pt_bin);
  }

  // --------------------------------------------------------------------------------
  // -- Prepare syncing - get nEvents for MC and Rec - 1
  // --------------------------------------------------------------------------------
  int nTracksMC1  = McTrack1->GetEntries();
  int nTracksRec1 = MatchedPairs1->GetEntries();

  int lastEvent1 = -1;
  int nEventsMC1  = 0;
  int nEventsRec1 = 0;

  for(int i = 0; i < nTracksMC1; i++){
    McTrack1->GetEntry(i);
    int current=Int_t(pEventId1);
    
    if (current != lastEvent1) {
      lastEvent1 = current;
      ++nEventsMC1;
    }
  }
  
  lastEvent1 = -1;
  for(int i = 0; i < nTracksRec1; i++){
    MatchedPairs1->GetEntry(i);
    int current=Int_t(EventId1);
    
    if (current != lastEvent1) {
      lastEvent1 = current;
      ++nEventsRec1;
    }
  }

  printf("nEvents (1) MC %d \n",  nEventsMC1);
  printf("nEvents (1) Rec %d \n", nEventsRec1);


  // --------------------------------------------------------------------------------
  // -- Prepare syncing - get nEvents for MC and Rec - 2
  // --------------------------------------------------------------------------------
  int nTracksMC2  = McTrack2->GetEntries();
  int nTracksRec2 = MatchedPairs2->GetEntries();

  int lastEvent2 = -1;
  int nEventsMC2  = 0;
  int nEventsRec2 = 0;

  for(int i = 0; i < nTracksMC2; i++){
    McTrack2->GetEntry(i);
    int current=Int_t(pEventId2);
    
    if (current != lastEvent2) {
      lastEvent2 = current;
      ++nEventsMC2;
    }
  }
  
  lastEvent2 = -1;
  for(int i = 0; i < nTracksRec2; i++){
    MatchedPairs2->GetEntry(i);
    int current=Int_t(EventId2);
    
    if (current != lastEvent2) {
      lastEvent2 = current;
      ++nEventsRec2;
    }
  }

  printf("nEvents (2) MC %d \n",  nEventsMC2);
  printf("nEvents (2) Rec %d \n", nEventsRec2);

#if 
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

  // -- event loop
  for (int idxEvent = 0; idxEvent < nEventsMC; ++idxEvent) {
    int nTracksPerEventMC  = 0;
    int nTracksPerEventRec = 0;

    // -- MC loop
    // --------------------------------------------------------------------------------
    lastEventMC = currentEventMC;
    for (; idxMC < nTracksMC; idxMC++) {
      McTrack->GetEntry(idxMC);
      currentEventMC = Int_t(pEventId);

      if (currentEventMC != lastEventMC) 
	break;
#if 0
      // --------------------------------------------------------------------------------
      // -- analyze current MC event
      //    - check cuts
      //    - fill nTracksPerEventMC;
      // --------------------------------------------------------------------------------
      if (pParentGeantId != 0)
	continue;
      
      if (pGeantId != PID )
	continue;
      
      // -- within Vertex R
      if ( (pVertexX-vxo)*(pVertexX-vxo) + (pVertexY-vyo)*(pVertexY-vyo) > max_r*max_r )
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
#endif
      ++nTracksPerEventMC;
    } // for(; idxMC < nTracksMC; idxMC++) {
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

    bool skipEventMC = false;

    // -- Rec loop
    // --------------------------------------------------------------------------------
    lastEventRec = currentEventRec;
    for (; idxRec < nTracksRec; idxRec++) {
      MatchedPairs->GetEntry(idxRec);
      currentEventRec = Int_t(EventId);

      if (currentEventRec != lastEventRec) 
	break;

      if (currentEventRec != lastEventMC) {
	skipEventMC = true;
	break;
      }

      // --------------------------------------------------------------------------------
      // -- analyze current Rec event
      //    - check cuts
      //    - fill nTracksPerEventRec;
      // --------------------------------------------------------------------------------
#if 0
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
      
      if (VertexZ > max_z || VertexZ < -max_z) 
	continue;
      
      if (fabs(EtaMc) > max_eta )
	continue;
      // --------------------------------------------------------------------------------
#endif
      ++nTracksPerEventRec;
    } // for(; idxRec < nTracksRec; idxRec++) {

    float efficiencyPerEvent = 0.;
#if 0
    if (skipEventMC) {
      // -- No Rec event present
      continue;
    }    
    else 
      efficiencyPerEvent = (nTracksPerEventMC != 0) ? nTracksPerEventRec / Float_t(nTracksPerEventMC) : 0;
    
    // -----------------------------------------------
    // -- get centrality
    int idxCent = -1;
    if (pCentrality16 == 15)     idxCent = 0;
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

    effProfile[idxCent]->Fill(PtMc, efficiencyPerEvent);
#endif
  } //  for (int idxEvent = 0 ; idxEvent < nEventsMC; ++idxEvent) {


  // TCanvas *c1 = new TCanvas("c1", "", 10, 10, 1200, 800);
  // c1->Divide(3,3);
  // for (int idxCent = 0; idxCent < nCent; ++idxCent) {
  //   c1->cd(idxCent+1);
  //   effProfile[idxCent]->GetYaxis()->SetRangeUser(0., 0.8);
  //   effProfile[idxCent]->SetMarkerStyle(21);
  //   effProfile[idxCent]->SetMarkerColor(kRed+2);
  //   effProfile[idxCent]->SetLineColor(kRed+2);

  //   effProfile[idxCent]->Draw("e1");
  // }

  // --------------------------------------------------------------------------------
  return;
  fOut->Write();
  fOut->Close();
}

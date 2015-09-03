#include "TCanvas.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "math.h"
#include <iostream>
#include <map>
#include <utility>

using namespace std;

void study_efficiencies_both_particles(const char* particle, int energy = 14) {

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
  
  int   max_events       = 500;

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
  float AllPts1, NPossible1, ParentGeantId1, GeantId1, mErrP1, NCommonHit1, EventId1, RunId1;
  MatchedPairs1->SetBranchAddress("EventId", &EventId1);
  MatchedPairs1->SetBranchAddress("RunId", &RunId1);
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
  float AllPts2, NPossible2, ParentGeantId2, GeantId2, mErrP2, NCommonHit2, EventId2, RunId2;
  MatchedPairs2->SetBranchAddress("EventId", &EventId2);
  MatchedPairs2->SetBranchAddress("RunId", &RunId2);
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
  float pRefMult1, pRefMultCorrected1, pCentralityWeight1, pCentrality161, pVertexX1, pVertexY1, pVertexZ1, pPtMc1, pPzMc1, pEtaMc1, pPhiMc1, pParentGeantId1, pGeantId1, pEventId1, pRunId1;
  McTrack1->SetBranchAddress("EventId", &pEventId1);
  McTrack1->SetBranchAddress("RunId", &pRunId1);
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
  float pRefMult2, pRefMultCorrected2, pCentralityWeight2, pCentrality162, pVertexX2, pVertexY2, pVertexZ2, pPtMc2, pPzMc2, pEtaMc2, pPhiMc2, pParentGeantId2, pGeantId2, pEventId2, pRunId2;
  McTrack2->SetBranchAddress("EventId", &pEventId2);
  McTrack2->SetBranchAddress("RunId", &pRunId2);
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
  TProfile *effProfile[nCent];
  TProfile *effProfileAllPt[nCent];
  TH1D* hpt_mc[nCent];
  TH1D* hpt_rec[nCent];

  TH1D* hpt_width[nCent];
  TH1D* hpt_relativeWidth[nCent];

  for (Int_t idxCent = 0; idxCent < nCent; idxCent++) {
    effProfile2D[idxCent] = new TProfile2D(Form("effProfile2D_%s", cent[idxCent]), 
					 Form("effProfile2D_%s", cent[idxCent]), Nbins, 0., 5., Nbins, 0., 5.);
    effProfile2D[idxCent]->GetXaxis()->Set(Nbins,pt_bin);
    effProfile2D[idxCent]->GetYaxis()->Set(Nbins,pt_bin);

    effProfileAllPt[idxCent] = new TProfile(Form("effProfileAllPt_%s", cent[idxCent]), 
					    Form("effProfileAllPt_%s", cent[idxCent]), 1, 0., 1.);

    effProfile[idxCent] = new TProfile(Form("effProfile_%s", cent[idxCent]), 
				       Form("effProfile_%s", cent[idxCent]), Nbins, 0., 5.);
    effProfile[idxCent]->GetXaxis()->Set(Nbins,pt_bin);

    hpt_mc[idxCent] = new TH1D(Form("hpt_mc_%s",  cent[idxCent]), "", Nbins, 0., 5.);
    hpt_mc[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_rec[idxCent] = new TH1D(Form("hpt_rec_%s", cent[idxCent]), "", Nbins, 0., 5.);
    hpt_rec[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_width[idxCent] = new TH1D(Form("hpt_width_%s", cent[idxCent]), "Width", Nbins, 0., 5.);
    hpt_width[idxCent]->GetXaxis()->Set(Nbins, pt_bin);

    hpt_relativeWidth[idxCent] = new TH1D(Form("hpt_relativeWidth_%s", cent[idxCent]), "Relative Width", Nbins, 0., 5.);
    hpt_relativeWidth[idxCent]->GetXaxis()->Set(Nbins, pt_bin);
  }

  // --------------------------------------------------------------------------------
  // -- Prepare syncing - get nEvents for MC and Rec - 1
  // --------------------------------------------------------------------------------
  int nTracksMC1  = McTrack1->GetEntries();
  int nTracksRec1 = MatchedPairs1->GetEntries();

  int lastEvent1 = -1;
  int nEventsMC1  = 0;
  int nEventsRec1 = 0;

  typedef std::pair<int,int> keyEventRun;
  std::map<keyEventRun, int> eventMC1;
  std::map<keyEventRun, int> eventMC2;
  std::map<keyEventRun, int> eventRec1;
  std::map<keyEventRun, int> eventRec2;

  for(int i = 0; i < nTracksMC1; i++){
    McTrack1->GetEntry(i);
    int current=Int_t(pEventId1);

    if (current != lastEvent1) {
      lastEvent1 = current;
      ++nEventsMC1;
      eventMC1[std::make_pair(Int_t(pRunId1), Int_t(pEventId1))] = i;
      cout << " ======  " << i << " " << nEventsMC1 << " - "<< eventMC1.size() << " || " << Int_t(pRunId1) << " : " << current << endl;
    }
  }
  return; 
  lastEvent1 = -1;
  for(int i = 0; i < nTracksRec1; i++){
    MatchedPairs1->GetEntry(i);
    int current=Int_t(EventId1);
    
    if (current != lastEvent1) {
      lastEvent1 = current;
      ++nEventsRec1;
      eventRec1[std::make_pair(Int_t(RunId1), Int_t(EventId1))] = i;
    }
  }

  printf("nEvents (1) MC %d  - %d\n",  nEventsMC1, eventMC1.size());
  printf("nEvents (1) Rec %d - %d\n", nEventsRec1, eventRec1.size());


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
      eventMC2[std::make_pair(Int_t(pRunId2), Int_t(pEventId2))] = i;
    }
  }
  
  lastEvent2 = -1;
  for(int i = 0; i < nTracksRec2; i++){
    MatchedPairs2->GetEntry(i);
    int current=Int_t(EventId2);
    
    if (current != lastEvent2) {
      lastEvent2 = current;
      ++nEventsRec2;
      eventRec2[std::make_pair(Int_t(RunId2), Int_t(EventId2))] = i;
    }
  }

  printf("nEvents (2) MC %d  - %d\n",  nEventsMC2, eventMC2.size());
  printf("nEvents (2) Rec %d - %d\n", nEventsRec2, eventRec2.size());

  // --------------------------------------------------------------------------------
  int idxMC2 = -1;
  lastEvent1 = -1;     
  nEventsMC1 = 0;
  int commonEvents1 = 0;
  for(int i = 0; i < nTracksMC1; i++){
    McTrack1->GetEntry(i);
    int current1=Int_t(pEventId1);
    
    if (current1 != lastEvent1) {
      lastEvent1 = current1;
      ++nEventsMC1;
   
      std::map<keyEventRun, int>::iterator iter = eventMC2.find(std::make_pair(Int_t(pRunId1), Int_t(pEventId1)));
      if (iter != eventMC2.end()) {
	idxMC2 = eventMC2[std::make_pair(Int_t(pRunId1), Int_t(pEventId1))];
	//	cout << nEventsMC1 << " -- " << i << " -- " << idxMC2 << endl;
	++commonEvents1;
      }
    }
  }

  // --------------------------------------------------------------------------------
  int idxMC1 = -1;
  lastEvent2 = -1;     
  nEventsMC2 = 0;
  int commonEvents2 = 0;
  for(int i = 0; i < nTracksMC2; i++){
    McTrack2->GetEntry(i);
    int current2=Int_t(pEventId2);
    
    if (current2 != lastEvent2) {
      lastEvent2 = current2;
      ++nEventsMC2;

      std::map<keyEventRun, int>::iterator iter = eventMC1.find(std::make_pair(Int_t(pRunId2), Int_t(pEventId2)));
      if (iter != eventMC1.end()) {
	idxMC1 = eventMC1[std::make_pair(Int_t(pRunId2), Int_t(pEventId2))];
	//	cout << nEventsMC2 << " -- " << i << " -- " << idxMC1 << endl;
	++commonEvents2;
      }
    }
  }

  printf("nEvents (common1) %d - %d \n", commonEvents1, nEventsMC1);
  printf("nEvents (common2) %d - %d \n", commonEvents2, nEventsMC2);

#if 0
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
      //    - get Centrality
      //    - check cuts
      //    - fill nTracksPerEventMC;
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
#if 0
      // --------------------------------------------------------------------------------
      // -- analyze current Rec event
      //    - get Centrality
      //    - check cuts
      //    - fill nTracksPerEventRec;
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
#if 0

    // -----------------------------------------------

    float efficiencyPerEvent = (nTracksPerEventMC != 0) ? nTracksPerEventRec / Float_t(nTracksPerEventMC) : 0;
    
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

    // -- all pts
    effProfileAllPt[idxCent]->Fill(0.5, efficiencyPerEvent);

    // -- pT dependent average over max_events
    ++(eventCounter[idxCent]);
    if (eventCounter[idxCent] >= max_events) {
      for (int idx = 1; idx <= hpt_mc[idxCent]->GetXaxis()->GetNbins(); idx++) { 
	double binContent = hpt_mc[idxCent]->GetBinContent(idx);
	if (binContent < 0.00001)
	  binContent = 1.;
	double ratio = hpt_rec[idxCent]->GetBinContent(idx)/(double)binContent;

	effProfile[idxCent]->Fill(hpt_mc[idxCent]->GetBinCenter(idx), ratio);
      }
      
      hpt_mc[idxCent]->Reset();
      hpt_rec[idxCent]->Reset();
      
      eventCounter[idxCent] = 0;
    }
#endif
  } //  for (int idxEvent = 0 ; idxEvent < nEventsMC; ++idxEvent) {

#if 0

  // --------------------------------------------------------------------------------
 
  for (int idxCent = 0; idxCent < nCent; ++idxCent) {
    for (int idx = 1; idx <= effProfile[idxCent]->GetXaxis()->GetNbins(); idx++) { 
      hpt_width[idxCent]->SetBinContent(idx, effProfile[idxCent]->GetBinError(idx));
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

    effProfileAllPt[idxCent]->GetYaxis()->SetRangeUser(0., 1.0);
    effProfileAllPt[idxCent]->SetMarkerStyle(22);
    effProfileAllPt[idxCent]->SetMarkerColor(kAzure);
    effProfileAllPt[idxCent]->SetLineColor(kAzure);

    hpt_width[idxCent]->GetYaxis()->SetRangeUser(0., 0.2);
    hpt_width[idxCent]->SetMarkerStyle(20);
    hpt_width[idxCent]->SetMarkerColor(kRed+2);
    hpt_width[idxCent]->SetLineColor(kRed+2);

    hpt_relativeWidth[idxCent]->GetYaxis()->SetRangeUser(0., 0.2);
    hpt_relativeWidth[idxCent]->SetMarkerStyle(21);
    hpt_relativeWidth[idxCent]->SetMarkerColor(kAzure);
    hpt_relativeWidth[idxCent]->SetLineColor(kAzure);
  }
  
  // --------------------------------------------------------------------------------
#endif
  fOut->Write();
  fOut->Close();
#endif
}

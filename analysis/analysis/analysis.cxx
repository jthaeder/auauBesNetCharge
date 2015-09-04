#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions; 
#pragma link C++ class PlotFile;
#endif
#ifndef __CINT__
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <map>
#include <utility>
#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include "TTree.h" 
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h" 
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TRandom3.h"
#endif
#include <stdio.h>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "picoDST.h"
#include "Centrality.h"
#include "StRefMultCorr.h"
#include "cuts.h"

using namespace std;

#define TRACK_THN 0
#define EVENT_THN 0
#define USE_RANDOM_EFF 0

void   AddHistSetCent(const Char_t *name, const Char_t *title);
void   FillHistSetCent(const Char_t *name, Int_t idx, Int_t cent);

void   Reset();

void   InitializeEventStats();
void   InitializeMultiplicityStats();
void   InitializeEventHists();
void   InitializeTrackHists();

Int_t  GetCentrality(Int_t nRefMultTracksCorr);

Bool_t FillEventStats(Int_t *aEventCuts);
void   FillMultiplicityStats(Double_t *aMult, Int_t mode);
void   FillEventHists(Int_t runId, Double_t *aEvent, Int_t mode);
void   FillTrackHists(Double_t *aTrack, Int_t mode);

Int_t  RefMult2Correction(Double_t vz, Int_t refmult2);

Double_t NN(Double_t, Int_t);

// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------  

Int_t    binHnEvent[11] = {  10,   21,   21,    201,   21,   601,   601,   3001,    2,    201,   21};
Double_t minHnEvent[11] = {-0.5, -2.0, -3.0, -100.0,  0.0,   0.0,   0.0,    0.0, -0.5, -100.0,  0.0};
Double_t maxHnEvent[11] = { 9.5,  2.0,  1.0,  100.0,  2.0, 600.0, 600.0, 3000.0,  1.5,  100.0, 20.0};

Int_t    binHnUnCorr[12] = {  10,  34,   21,    3,  50,   51,   51,   51, 101,    2,    2, 81};
Double_t minHnUnCorr[12] = {-0.5, 0.1, -1.0, -1.5,   0,  0.0,  0.0,  0.0, 0.0, -0.5, -0.5, -4};
Double_t maxHnUnCorr[12] = { 9.5, 3.0,  1.0,  1.5,  10, 50.0, 50.0, 50.0, 1.0,  1.5,  1.5,  4};

Int_t nMultSets = 5;

const Char_t* multNames[]  = {
  "_base", 
  "_before_TOF_rejected",  "_before_TOF", 
  "_after_TOF_rejected", "_after_TOF"};

const Char_t* multTitles[] = {
  " (before)",
  " (before TOF cuts rejected)", " (before TOF cuts)", 
  " (after TOF cuts rejected)", " (after TOF cuts)"};

// ----------------------------------------------------------------------------  
// -- 
// ----------------------------------------------------------------------------  
const int   nEnergies       = 8;
const char *energies[]      = {  "7",   "11",   "14",   "19",   "27",   "39",   "62", "200"};
const char *exactEnergies[] = {"7.7", "11.5", "14.5", "19.6", "27.0", "39.0", "62.4", "200"};

const Char_t* name = "NetCharge";

enum particleCharge {kPOS, kNEG, kNET, kParticleCharge};

const Char_t* asParticleName[2]  =  {"neg",  "pos" };
const Char_t* asParticleTitle[2] =  {"neg.", "pos."};

const Int_t fHEventStatMax   = 10; 
const Char_t *aEventCutNames[]   = {"all", "bad run", "trigger", "#it{v}_{z} < #it{v}_{z}^{max}", 
				    "#it{v}_{z}- #it{v}_{z}^{vpd} < 3cm", "shifted vtx <1", 
				    "centrality", "nTOFMatch>2","nTOFMatch cut", "accepted"};

double randomEff[2][9];

// ----------------------------------------------------------------------------  
// -- Globals
// ----------------------------------------------------------------------------  

TList                *fOutList;               //! Output data container
// =======================================================================
Int_t                 fOrder = 8;             //  Max order of higher order distributions
// -----------------------------------------------------------------------
Int_t                 fNNp   = 3;             //  N sets of arrays of particle/anti-particle counts
                                              //   0 is all 
                                              //   1,2 for arbitrary subset
Int_t               **fNp;                    //  Array of particle/anti-particle counts
// =======================================================================
THnSparseD           *fHnTrackUnCorr;         //  THnSparseD : uncorrected probe particles
THnSparseD           *fHnEvent;               //  THnSparseD : event
// =======================================================================

Int_t                 energyIdx = -1;         // energyIdx

// ----------------------------------------------------------------------------  
// -- 
// ----------------------------------------------------------------------------  
Int_t main(int argc, char** argv) {

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // -----------------------------------------------------------------------
  // -- Read arguments
  // -----------------------------------------------------------------------
  if (argc > 1) {
    for (Int_t idx = 1; idx < argc; ++idx) {
      TString argument(argv[idx]);

      if (argument.BeginsWith("--energy=")) {
	TString parameter=argument.ReplaceAll("--energy=", "");
	parameter.Remove(TString::kLeading, ' '); 
	
	for (int enerIdx = 0; enerIdx < nEnergies; ++enerIdx) {
	  if (!parameter.CompareTo(exactEnergies[enerIdx]) ||
	      !parameter.CompareTo(Form("%sGeV",exactEnergies[enerIdx])) ||
	      !parameter.CompareTo(energies[enerIdx]) ||
	      !parameter.CompareTo(Form("%sGeV",energies[enerIdx]))) {
	    energyIdx = enerIdx;
	  }
	}
      } // if (argument.BeginsWith("--energy=")){
    } // for (Int_t idx = 1; idx < argc; ++idx) {
  } // if (argc > 1) {

  if (energyIdx == -1) {
    cout << "No energy selected" << endl;
    return 0;
  }

  // -----------------------------------------------------------------------
  // -- Set Globals
  // -----------------------------------------------------------------------

  // -- Add List
  fOutList = new TList;
  fOutList->SetName(Form("f%s", name));
  fOutList->SetOwner(kTRUE);
 
  // -------------------------------------------
  fNp = new Int_t*[fNNp];
  for (Int_t ii = 0 ; ii < fNNp; ++ii)
    fNp[ii] = new Int_t[2];

  // -----------------------------------------------------------------------
  // -- Get random efficiencies
  // -----------------------------------------------------------------------
#if USE_RANDOM_EFF
  TString basePath("/project/projectdirs/starprod/rnc/jthaeder/NetCharge/data");
  TFile* fileAll = TFile::Open(Form("%s/2015-05-20_0.5/Sum_NetCharge_AuAu14.5GeV_Vz50.root",         basePath.Data()));
  TFile* fileSub = TFile::Open(Form("%s/2015-06-20_delta_0.5_6/Sum_NetCharge_AuAu14.5GeV_Vz50.root", basePath.Data()));

  // -----------------------------------------------------------

  TList *listAll     = static_cast<TList*>(fileAll->Get("fNetCharge"));
  TList *distListAll = static_cast<TList*>(listAll->FindObject("fDist"));
  TH2F  *hDistposAll = static_cast<TH2F*>(distListAll->FindObject("hDistpos"));
  TH2F  *hDistnegAll = static_cast<TH2F*>(distListAll->FindObject("hDistneg"));

  // -----------------------------------------------------------

  TList *listSub     = static_cast<TList*>(fileSub->Get("fNetCharge"));
  TList *distListSub = static_cast<TList*>(listSub->FindObject("fDist"));
  TH2F  *hDistposSub = static_cast<TH2F*>(distListSub->FindObject("hDistpos"));
  TH2F  *hDistnegSub = static_cast<TH2F*>(distListSub->FindObject("hDistneg"));

  // -----------------------------------------------------------
  for (int idxCent = 0; idxCent <9 ; ++idxCent) {
    TH1D *prPosAll = hDistposAll->ProjectionY(Form("%s_%d", hDistposAll->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");
    TH1D *prPosSub = hDistposSub->ProjectionY(Form("%s_%d", hDistposSub->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");

    TH1D *prNegAll = hDistnegAll->ProjectionY(Form("%s_%d", hDistnegAll->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");
    TH1D *prNegSub = hDistnegSub->ProjectionY(Form("%s_%d", hDistnegSub->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");

    randomEff[0][idxCent] = prNegSub->GetMean()/ prNegAll->GetMean();
    randomEff[1][idxCent] = prPosSub->GetMean()/ prPosAll->GetMean();

    cout << "neg " <<  idxCent << "  " << randomEff[0][idxCent] << " -> " << prNegSub->GetMean() << " / " << prNegAll->GetMean() <<  endl;
    cout << "pos " <<  idxCent << "  " << randomEff[1][idxCent] << " -> " << prPosSub->GetMean() << " / " << prPosAll->GetMean() <<  endl;
  } // for (int idxCent = 0; idxCent <9 ; ++idxCent) {
#endif
  // ------------------------------------------------------------------
  // -- Get event container
  // ------------------------------------------------------------------

#if EVENT_THN
  // -- Event
  fOutList->Add(new THnSparseD("hnEvent", "cent:vx:vy:vz:shiftedVtx:nRefMult2:nRefMult2Corr:nTracks:isRejected", 9, 
			       binHnEvent, minHnEvent, maxHnEvent));  

  fHnEvent = static_cast<THnSparseD*>(fOutList->Last());
  fHnEvent->Sumw2(); 
  fHnEvent->GetAxis(0)->SetTitle("centrality");
  fHnEvent->GetAxis(1)->SetTitle("#it{v}_{x} (cm)");
  fHnEvent->GetAxis(2)->SetTitle("#it{v}_{y} (cm)");
  fHnEvent->GetAxis(3)->SetTitle("#it{v}_{z} (cm)");
  fHnEvent->GetAxis(4)->SetTitle("shifted vtx (cm)");
  fHnEvent->GetAxis(5)->SetTitle("refMult2");
  fHnEvent->GetAxis(6)->SetTitle("refMult2Corr");
  fHnEvent->GetAxis(7)->SetTitle("nTracks");
  fHnEvent->GetAxis(8)->SetTitle("isRejected");
  fHnEvent->GetAxis(9)->SetTitle("#it{v}_{z}^{vpd} (cm)");
  fHnEvent->GetAxis(10)->SetTitle("#Delta#it{v} (cm)");
#endif

  // ------------------------------------------------------------------
  // -- Get tracks container
  // ------------------------------------------------------------------

#if TRACK_THN
  // -- UnCorrected
  fOutList->Add(new THnSparseD("hnTrackUnCorr", "cent:pt:eta:sign:dca:nHitsDedx:nHitsFit:nFitPoss:ratio:isInRefMult:isTrackAccepted;nSigmaProton", 
			       12, binHnUnCorr, minHnUnCorr, maxHnUnCorr));  
  fHnTrackUnCorr = static_cast<THnSparseD*>(fOutList->Last());
  fHnTrackUnCorr->Sumw2(); 
  fHnTrackUnCorr->GetAxis(0)->SetTitle("centrality");
  fHnTrackUnCorr->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHnTrackUnCorr->GetAxis(2)->SetTitle("#eta");
  fHnTrackUnCorr->GetAxis(3)->SetTitle("sign");
  fHnTrackUnCorr->GetAxis(4)->SetTitle("dca");
  fHnTrackUnCorr->GetAxis(5)->SetTitle("nHitsDedx");
  fHnTrackUnCorr->GetAxis(6)->SetTitle("nHitsFit");
  fHnTrackUnCorr->GetAxis(7)->SetTitle("nFitPoss");
  fHnTrackUnCorr->GetAxis(8)->SetTitle("nHitsFit/nFitPoss");
  fHnTrackUnCorr->GetAxis(9)->SetTitle("isInRefMult");
  fHnTrackUnCorr->GetAxis(10)->SetTitle("isTrackAccepted");
  fHnTrackUnCorr->GetAxis(11)->SetTitle("nSigmaProton");
#endif

  // -----------------------------------------------------------------------
  // -- Get PicoDsts
  // -----------------------------------------------------------------------
  picoDST* pico = new picoDST();
  StRefMultCorr* refmult2Corr = new StRefMultCorr("refmult2");

  // -----------------------------------------------------------------------
  // -- Get refMult / refMult2 from muDst
  // -----------------------------------------------------------------------

  typedef std::pair<int,int> keyEventRun;
  std::map<keyEventRun, int> refMultMap;
  std::map<keyEventRun, int> refMult2Map;

  if (energyIdx != 2) {
    ifstream refMultFile;
    
    // -- open in working dir
    refMultFile.open("file.list.refMult");
    if (!refMultFile.is_open()) {
      cout << "Couldn't open refMult file" << endl;
      exit(EXIT_FAILURE);
    }
    
    if (refMultFile.is_open()) {
      int run, event, gRefMult, refMult, refMult2;
      
      refMultFile.ignore(10000,'\n');
      while(1) {
	refMultFile >> run >> event >> gRefMult >> refMult >> refMult2;
	
	// -- break at at of file
      if (refMultFile.eof())
	break;
      
      // -- break if error occured during reading
      if (!refMultFile.good()) 
	break;
      
      refMultMap[std::make_pair(run,event)] = refMult;
      refMult2Map[std::make_pair(run,event)] = refMult2;
      }
      
      refMultFile.close();
    }
  }

  // -----------------------------------------------------------------------
  // -- Define Histograms
  // -----------------------------------------------------------------------
  
  InitializeEventStats();
  InitializeMultiplicityStats();
  InitializeEventHists();
  InitializeTrackHists();

  //  TString sTitle(Form("|#eta|<%.1f #it{p}_{T} [%.1f,%.1f]", etaAbsRange[1], ptRange[0], ptRange[1]));
  TString sTitle(Form("%.1f<#eta<%.1f #it{p}_{T} [%.1f,%.1f]", etaAbsRange[0], etaAbsRange[1], ptRange[0], ptRange[1]));
  printf ("TITLE : .... %s", sTitle.Data());
  AddHistSetCent("Dist",       sTitle.Data());
  AddHistSetCent("Dist_lower", sTitle.Data());
  AddHistSetCent("Dist_upper", sTitle.Data());
  
  //*************************************************************************************************

  fOutList->Add(new TList);
  TList *fijkhList = static_cast<TList*>(fOutList->Last());
  fijkhList->SetName(Form("f%sFijkh", name));
  fijkhList->SetOwner(kTRUE);

  TProfile *fact[9][9][9][9][10];
  for(int l=0;l<=9;l++) 
    for(int i=0;i<=8;i++) 
      for(int j=0;j<=8;j++) 
	for(int k=0;k<=8;k++) 
	  for(int h=0;h<=8;h++) 
	    if((i+j+k+h)<=8) {
	      fijkhList->Add(new TProfile(Form("f%d%d%d%d_Cent%d",i,j,k,h,l),Form("f%d%d%d%d_Cent%d",i,j,k,h,l),1001,-0.5,1000.5));
	      fact[i][j][k][h][l] = static_cast<TProfile*>(fijkhList->Last());
	    }

  //*************************************************************************************************
  
  Long64_t nEvents = (static_cast<TTree*>(pico->fChain))->GetEntries();
  cout << nEvents << " Events to be done!" << endl;

  // ------------------------------------------------------------------
  // -- Event LOOP
  // ------------------------------------------------------------------
  for (Long64_t idxEvent = 0; idxEvent < nEvents; idxEvent++) {
    Long64_t ientry = pico->LoadTree(idxEvent);
    if (ientry < 0) 
      break;
    
    (static_cast<TTree*>(pico->fChain))->GetEntry(idxEvent);  

    if (idxEvent%5000 == 0)
      cout << "Event Loop: " << idxEvent << " events have been done!" << endl;

    // ------------------------------------------------------------------
    // -- Reset Event
    // ------------------------------------------------------------------

    Reset();
    
    Int_t runId   = pico->Event_mRunId[0];
    Int_t eventId = pico->Event_mEventId[0];

    // ------------------------------------------------------------------
    // -- Track loop - for refMult2 / TOF tracks - track multiplicity
    // ------------------------------------------------------------------
    Int_t nTracks         = pico->Tracks_;      
    Int_t nGlobalTracks   = pico->Event_mNumberOfGlobalTracks[0]; 
    Int_t nPrimaryTracks  = 0;
    
    Int_t nRefMultTracksPico  = pico->Event_mRefMultNeg[0] + pico->Event_mRefMultPos[0];
    Int_t nRefMultTracksMuDst = refMultMap[std::make_pair(runId, eventId)];
    Int_t nRefMultTracks  = (energyIdx == 2) ? nRefMultTracksPico : nRefMultTracksMuDst;

    Int_t nRefMult2TracksPico = pico->Event_mRefMult2NegEast[0] + pico->Event_mRefMult2PosEast[0] +
      pico->Event_mRefMult2NegWest[0] +  pico->Event_mRefMult2PosWest[0];
    Int_t nRefMult2TracksMuDst = refMult2Map[std::make_pair(runId, eventId)];
    Int_t nRefMult2Tracks = (energyIdx == 2) ? 0 : nRefMult2TracksMuDst;

    Int_t nTOFMatch = 0;
    
    if (energyIdx == 2) 
      nTOFMatch = pico->Event_mNBTOFMatch[0];
    else {
      for (int idx =0; idx < nTracks; ++idx) {
	if (pico->Tracks_mBTofMatchFlag[idx] > 0)
	  nTOFMatch++; 
      }
    }
    
    // ------------------------------------------------------------------

    for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
      Double_t pxcm     = pico->Tracks_mPMomentum_mX1[idxTrack];
      Double_t pycm     = pico->Tracks_mPMomentum_mX2[idxTrack];
      Double_t pzcm     = pico->Tracks_mPMomentum_mX3[idxTrack];
  
      Double_t pcm      = TMath::Sqrt(pxcm*pxcm + pycm*pycm + pzcm*pzcm);
      Double_t eta      = 0.5 * TMath::Log((pcm+pzcm)/(pcm-pzcm));
      
      Double_t DCA      = pico->Tracks_mGDca[idxTrack]/1000.;
      Int_t    nHitsFit = pico->Tracks_mNHitsFit[idxTrack];
      
      // -- count for refMult2 track -- only for 14.5 GeV
      // ------------------------------------------------------------------
      if (energyIdx == 2 && TMath::Abs(eta) > etaAbsRangeRefMult[0] && TMath::Abs(eta) <= etaAbsRangeRefMult[1]  && 
	  TMath::Abs(nHitsFit) > nHitsFitRefMult  && DCA < dcaMaxRefMult)
	++nRefMult2Tracks;

      // -- count for primary track
      // ------------------------------------------------------------------
      if (! (pxcm == 0 && pycm == 0 && pzcm == 0) )
	++nPrimaryTracks;
    }
    
    if (nRefMultTracksPico - nRefMultTracksMuDst != 0)
      cout << "  DELTA REFMULT " << nRefMultTracks  << " || Pico:" << nRefMultTracksPico  << " <> MuDst:" << nRefMultTracksMuDst  << endl;
    cout << "  REFMULT2      " << nRefMult2Tracks << " || Pico:" << nRefMult2TracksPico << " <> MuDst:" << nRefMult2TracksMuDst << endl;

    // ------------------------------------------------------------------
    // -- Event Cuts
    // ------------------------------------------------------------------

    Int_t *aEventCuts = new Int_t[fHEventStatMax];
    // set aEventCuts[ii] to 1 in case of reject
  
    for (Int_t ii=0;ii<fHEventStatMax; ++ii)
      aEventCuts[ii] = 0;
    
    Int_t iCut = 0;

    // -- 0 - Before Event Cuts
    aEventCuts[iCut] = 0;

    // -- bad runs / trigger 
    // ------------------------------------------------------------------

    // -- 1 - bad run
    ++iCut;
    if (refmult2Corr->isBadRun(runId)) 
      aEventCuts[iCut] = 1;

    // -- 2 - trigger - bit 5 || 6
    Int_t triggerId = pico->Event_mTriggerWord[0];
    ++iCut;

    if (energyIdx == 0) {
      if (!((triggerId&0x1)||((triggerId>>1)&0x1))) 
	aEventCuts[iCut] = 1;	
    }
    else if (energyIdx == 1) {
      if (!((triggerId&0x1)||((triggerId>>1)&0x1))) 
	aEventCuts[iCut] = 1;
    }
    else if (energyIdx == 2) {
      if (!(((triggerId>>5)&0x1)||((triggerId>>6)&0x1))) 
	aEventCuts[iCut] = 1;
    }
    else if (energyIdx == 3) {
      if(!((triggerId&0x1)||((triggerId>>1)&0x1)||((triggerId>>2)&0x1))) 
	aEventCuts[iCut] = 1;
    }
    else if (energyIdx == 4) {
      if (!((triggerId&0x1))) 
	aEventCuts[iCut] = 1;
    }
    else if (energyIdx == 5) {
      if (!((triggerId&0x1))) 
	aEventCuts[iCut] = 1;
    }
    else if (energyIdx == 6) {
      if(!((triggerId&0x1)||((triggerId>>1)&0x1)||((triggerId>>2)&0x1))) 
	aEventCuts[iCut] = 1;
    }
    else if (energyIdx == 7) {
      if(!(((triggerId)&0x1)||((triggerId>>1)&0x1)||((triggerId>>2)&0x1)||((triggerId>>3)&0x1))) 
	aEventCuts[iCut] = 1;
    }

    // -- Vertex cuts
    // ------------------------------------------------------------------

    Float_t vx = pico->Event_mPrimaryVertex_mX1[0];
    Float_t vy = pico->Event_mPrimaryVertex_mX2[0];
    Float_t vz = pico->Event_mPrimaryVertex_mX3[0];

    // -- 3 - vertexZ cut
    ++iCut;
    Double_t vzMax = (energyIdx == 0) ? 50 : 30; // cm
    if (TMath::Abs(vz) > vzMax) 
      aEventCuts[iCut] = 1;

    // -- 4 - vpd vertex cut
    ++iCut;  
    Double_t deltaVzMax = 3;
    Float_t  vpdVz      = Float_t(pico->Event_mVzVpd[0]);
    Double_t deltaVz    = TMath::Abs(vpdVz-vz);
    if (energyIdx > 4 && deltaVz > deltaVzMax)
      aEventCuts[iCut] = 1;
  
    // -- 5 - shifted vertex cut
    ++iCut;
    Double_t shiftedVtx = TMath::Sqrt(vx*vx + (vy+0.89)*(vy+0.89));
    if (energyIdx == 2 && shiftedVtx >= 1) 
      aEventCuts[iCut] = 1;
    
    // -- Centrality cuts 
    // ------------------------------------------------------------------

    Int_t nRefMult2TracksCorr, centrality;

    // -- 14.5 GeV
    if (energyIdx == 2) {
      // -- Correct refmult vz dependent 
      nRefMult2TracksCorr = RefMult2Correction(vz, nRefMult2Tracks);
      
      // -- Get centrality and fill stats
      centrality = GetCentrality(nRefMult2TracksCorr);
    }
    // -- other energies
    else {
      refmult2Corr->init(runId);
      refmult2Corr->initEvent(nRefMult2Tracks, vz);

      // -- Correct refmult vz dependent 
      nRefMult2TracksCorr = refmult2Corr->getRefMultCorr();

      centrality = 8 - refmult2Corr->getCentralityBin9();
    }

    // -- 6 - cut for centrality range
    ++iCut;
    if (centrality >= fNCentralityBins || centrality == -1)
      aEventCuts[iCut] = 1;

    // -- TOFMatch cuts -pileup 
    // ------------------------------------------------------------------

    // -- 7 - cut for nTOFMatch >2
    ++iCut;
    if (nTOFMatch <= 2)
      aEventCuts[iCut] = 1;
    
    // -- 8 - cut for nTOFMatch vs nRefMult
    //     y = 0.78 x - 10.2
    ++iCut;
    if (0.78*nRefMultTracks - 10.2 > nTOFMatch ) 
      aEventCuts[iCut] = 1;

    // -- Fill statistics
    // ------------------------------------------------------------------
    Bool_t isRejectedWithoutTOF = (aEventCuts[0] || aEventCuts[1] || aEventCuts[2] || aEventCuts[3] || 
				   aEventCuts[4] || aEventCuts[5] || aEventCuts[6]) ? kTRUE : kFALSE;

    Bool_t isRejected = FillEventStats(aEventCuts);

    // -- fill ThnSparse - events
    // ------------------------------------------------------------------
    Double_t aEvent[11] = {Double_t(centrality), vx, vy, vz, shiftedVtx, 
			   Double_t(nRefMult2Tracks), Double_t(nRefMult2TracksCorr), Double_t(nTracks), Double_t(isRejected),
			   vpdVz, deltaVz};
      
    Double_t aMult[7]  = {Double_t(centrality), Double_t(nRefMultTracks), 
			  Double_t(nRefMult2Tracks), Double_t(nRefMult2TracksCorr), 
			  Double_t(nGlobalTracks), Double_t(nPrimaryTracks),
			  Double_t(nTOFMatch) };

#if EVENT_THN
    fHnEvent->Fill(aEvent);
#endif
    
    // -- Fill eventHists
    FillEventHists(runId, aEvent, 0);
    
    // -- Fill multiplicty stats - before event cuts
    FillMultiplicityStats(aMult, 0);
      
    // -- Reject Event - but not TOF cuts
    // ------------------------------------------------------------------
    if (isRejectedWithoutTOF) {
      FillMultiplicityStats(aMult, 1);
      continue;
    }

    // -- Fill multiplicty stats - after event cuts (but not TOF cuts)
    FillMultiplicityStats(aMult, 2);

    // -- Reject Event
    // ------------------------------------------------------------------
    if (isRejected) {
      FillMultiplicityStats(aMult, 3);
      continue;
    }

    // if ( nRefMult2Tracks != nRefMult2TracksPico && TMath::Abs(nRefMult2Tracks - nRefMult2TracksPico) > 1)
    //   cout << " nRefMult2Tracks " << nRefMult2Tracks << " -> (pico) " << nRefMult2TracksPico << endl;
    // continue;

    // -- Fill multiplicty stats
    FillMultiplicityStats(aMult, 4);

    // -- Fill eventHists
    FillEventHists(runId, aEvent, 1);

    // ------------------------------------------------------------------
    // -- Track loop - track multiplicity
    // ------------------------------------------------------------------
    for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
      
      Double_t pxcm   = pico->Tracks_mPMomentum_mX1[idxTrack];
      Double_t pycm   = pico->Tracks_mPMomentum_mX2[idxTrack];
      Double_t pzcm   = pico->Tracks_mPMomentum_mX3[idxTrack];
      if (pxcm == 0 && pycm == 0 && pzcm == 0) 
	continue;
      
      Double_t pt     = TMath::Sqrt(pxcm*pxcm + pycm*pycm);
      Double_t pcm    = TMath::Sqrt(pt*pt + pzcm*pzcm);
      
      Double_t eta    = 0.5 * TMath::Log((pcm+pzcm)/(pcm-pzcm));
      Double_t DCA    = pico->Tracks_mGDca[idxTrack]/1000.;
      
      Int_t nHitsDedx = pico->Tracks_mNHitsDedx[idxTrack];
      Int_t nHitsFit  = pico->Tracks_mNHitsFit[idxTrack];
      Int_t nFitPoss  = pico->Tracks_mNHitsMax[idxTrack];
      
      Double_t ratio  =  (1+fabs(nHitsFit))/(1+nFitPoss);
      
      Double_t sign   = (nHitsFit > 0) ? +1 : -1;
      
      Float_t nSigmaProton = pico->Tracks_mNSigmaProton[idxTrack]/100.;
      
      // -- is in RefMult flag
      Bool_t isInRefMult = (TMath::Abs(eta) > etaAbsRangeRefMult[0] && TMath::Abs(eta) <= etaAbsRangeRefMult[1]  && 
			    TMath::Abs(nHitsFit) > nHitsFitRefMult  && DCA < dcaMaxRefMult) ? kTRUE : kFALSE;
      
      // -- is track accepted flag - kinematics
      Bool_t isTrackAcceptedKin = (eta > etaAbsRange[0] && eta < etaAbsRange[1] && 
				   pt > ptRange[0] && pt < ptRange[1]) ? kTRUE : kFALSE;
      
      // -- is track accepted flag - clusters/dca
      Bool_t isTrackAcceptedCut = (TMath::Abs(nHitsFit) > nHitsFitMin && DCA < dcaMax && 
				   nHitsDedx > nHitsDedxMin && ratio > ratioNHitsFitNFitPossMin) ? kTRUE : kFALSE;
      
      // -->> is track accepted  - clusters/dca && kinematics
      Bool_t isTrackAccepted = (isTrackAcceptedKin && isTrackAcceptedCut);

      
      // -- is track Spallation proton/anti-proton
      Bool_t isTrackSpallationProton = (isTrackAccepted && pt > ptRangeSpallation[0] && pt < ptRangeSpallation[1] &&
					TMath::Abs(nSigmaProton) < nSigmaProtonMaxSpallation)  ? kTRUE : kFALSE;
      
      // -- fill ThnSparse - tracks
      // ------------------------------------------------------------------
      Double_t aTrack[12] = {Double_t(centrality), pt, eta, sign, DCA,
			     Double_t(nHitsDedx), Double_t(nHitsFit), Double_t(nFitPoss), ratio, 
			     Double_t(isInRefMult), Double_t(isTrackAccepted), nSigmaProton};
      
#if TRACK_THN      
      fHnTrackUnCorr->Fill(aTrack);
#endif

      FillTrackHists(aTrack, 0);

      // -- reject track
      // ------------------------------------------------------------------
      if (!isTrackAccepted)
	continue;
      
      FillTrackHists(aTrack, 1);

      // -- reject spallation proton
      // ------------------------------------------------------------------
      if (isTrackSpallationProton)
	continue;
      
      FillTrackHists(aTrack, 2);

      // ------------------------------------------------------------------
      // -- Add up for event multiplicity
      // ------------------------------------------------------------------
      //  idxPart = 0 -> anti particle
      //  idxPart = 1 -> particle
      Int_t idxPart = (sign < 0) ? 0 : 1;

#if USE_RANDOM_EFF      
      // -- discard a random amount of tracks
      if (gRandom->Rndm() > randomEff[idxPart][centrality])
	continue;
#endif
 
     // -- in full pt Range
      fNp[0][idxPart] += 1;
      
      // -- divide in 2 parts
      if (pt < ptMidPoint)
	fNp[1][idxPart] += 1;
      else
	fNp[2][idxPart] += 1;
    } // for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
    
    // -- Fill histograms
    // ------------------------------------------------------------------
    FillHistSetCent("Dist",       0, centrality);
    FillHistSetCent("Dist_lower", 1, centrality);
    FillHistSetCent("Dist_upper", 2, centrality);

    // ------------------------------------------------------------------

    for(int i=0;i<=9;i++) 
      for(int j=0;j<=9;j++) 
	for(int k=0;k<=9;k++) 
	  for(int h=0;h<=9;h++) 
	    if((i+j+k+h)<=8) {
	      fact[i][j][k][h][centrality+1]->Fill(nRefMult2TracksCorr, NN(fNp[1][1],i) * NN(fNp[2][1],j) * NN(fNp[1][0],k) * NN(fNp[2][0],h));
	      fact[i][j][k][h][0]->Fill(           nRefMult2TracksCorr, NN(fNp[1][1],i) * NN(fNp[2][1],j) * NN(fNp[1][0],k) * NN(fNp[2][0],h));
	    }
    
  } // for (Long64_t idxEvent = 0; idxEvent < nEvents; idxEvent++) {

  cout << "Event Loop: " << nEvents << " events have been done!" << endl;

  // ------------------------------------------------------------------
  // -- Write Outlist
  // ------------------------------------------------------------------

  TH1::AddDirectory(oldStatus);

  TString collision(Form("Moments_hist_AuAu%sGeV_charge", exactEnergies[energyIdx]));
  TFile *outFile = new TFile(Form("%s.root", collision.Data()), "recreate");   
  outFile->cd();
  fOutList->Write(fOutList->GetName(), TObject::kSingleKey);
  outFile->Close();
  
  return 0;
 }

//________________________________________________________________________
Int_t RefMult2Correction(Double_t vz, Int_t refmult2) {
  // -- Correction of refmult2

  Double_t refMultZ = aRefMultCorrPar[0] + 
    aRefMultCorrPar[1]*vz + aRefMultCorrPar[2]*vz*vz + aRefMultCorrPar[3]*vz*vz*vz + 
    aRefMultCorrPar[4]*vz*vz*vz*vz + aRefMultCorrPar[5]*vz*vz*vz*vz*vz + aRefMultCorrPar[6]*vz*vz*vz*vz*vz*vz;

  Double_t Hovno    = (aRefMultCorrPar[0] + aRefMultCorrPar[7]) / refMultZ;

  Double_t refMultD = Double_t(refmult2) + gRandom->Rndm(); // random sampling over bin width -> avoid peak structures in corrected distribution
  
  return Int_t(refMultD * Hovno);
}

//________________________________________________________________________
double NN(double num, int order){
  return (order == 0) ? 1 : NN(num,order-1)*(num-order+1);
}

//________________________________________________________________________
void AddHistSetCent(const Char_t *name, const Char_t *title)  {
   // -- Add histogram sets for particle and anti-particle
   //    dependence : centrality

   TString sName(name);
   TString sTitle(title);

   // -- Add List
   fOutList->Add(new TList);
   TList *list = static_cast<TList*>(fOutList->Last());
   list->SetName(Form("f%s", name));
   list->SetOwner(kTRUE);

   // -- Create Titles
   TString sNetTitle(Form("N_{%s} - N_{%s}", asParticleTitle[1], asParticleTitle[0]));
   TString sSumTitle(Form("N_{%s} + N_{%s}", asParticleTitle[1], asParticleTitle[0]));
   
   // -----------------------------------------------------------------------------------------------
   
   // -- Add Particle / Anti-Particle Distributions
   for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
     list->Add(new TH2D(Form("h%s%s", name, asParticleName[idxPart]), 
			Form("N_{%s} : %s;Centrality;N_{%s}", asParticleTitle[idxPart], sTitle.Data(), asParticleTitle[idxPart]),
			fNCentralityBins, centBinRange[0], centBinRange[1], 601, -0.5, 600.49));
   } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
 
   // -- Add Particle vs Anti-Particle Distribution
   list->Add(new TH2D(Form("h%s%s%s", name, asParticleName[0],asParticleName[1]), 
		      Form("N_{%s} vs N_{%s} : %s;N_{%s};N_{%s}", asParticleTitle[0] ,asParticleTitle[1], sTitle.Data(), asParticleTitle[0], asParticleTitle[1]),
		      601, -0.5, 600.49,601, -0.5, 600.49));

   // -- Add NetParticle Distributions
   list->Add(new TH2D(Form("h%sNetCharge", name), 
		      Form("%s : %s;Centrality;%s", sNetTitle.Data(), sTitle.Data(), sNetTitle.Data()), 
		      fNCentralityBins, centBinRange[0], centBinRange[1], 601, -300.5, 300.49));

   // -- Add NetParticle vs SumParticle
   list->Add(new TH2D(Form("h%sNetChargeOverSum", name), 
		      Form("(%s)/(%s) : %s;Centrality;(%s)/(%s)", sNetTitle.Data(), sSumTitle.Data(), sTitle.Data(), sNetTitle.Data(), sSumTitle.Data()), 
		      fNCentralityBins, centBinRange[0], centBinRange[1], 41, -2.5, 2.49));   
   return;
}

//________________________________________________________________________
 void FillHistSetCent(const Char_t *name, Int_t idx, Int_t cent)  {
  // -- Fill histogram sets for particle and anti-particle
  //    dependence : centrality 
  
  // -- Get List
  TList *list = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));
  
  // -- Get Centrality Bin
  Float_t centralityBin = Float_t(cent);

  // -- Select MC or Data
  Int_t **np = fNp;

  // -----------------------------------------------------------------------------------------------

  Int_t sumNp   = np[idx][1]+np[idx][0];  // p + pbar
  Int_t deltaNp = np[idx][1]-np[idx][0];  // p - pbar

  // -- Fill Particle / Anti-Particle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%s%s", name, asParticleName[0]))))->Fill(centralityBin, np[idx][0]);
  (static_cast<TH2D*>(list->FindObject(Form("h%s%s", name, asParticleName[1]))))->Fill(centralityBin, np[idx][1]);

  (static_cast<TH2D*>(list->FindObject(Form("h%s%s%s", name, asParticleName[0], asParticleName[1]))))->Fill(np[idx][0],np[idx][1]);

  // -- Fill NetParticle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%sNetCharge",  name))))->Fill(centralityBin, deltaNp);

  // -- Fill NetParticle vs SumParticle
  Double_t deltaNpOverSumNp = (sumNp == 0.) ? 0. : deltaNp/Double_t(sumNp);
  (static_cast<TH2D*>(list->FindObject(Form("h%sNetChargeOverSum", name))))->Fill(centralityBin, deltaNpOverSumNp);

  // -----------------------------------------------------------------------------------------------

  return;
}

//________________________________________________________________________
void Reset() {
  // -- Reset eventwise

  // -- Reset N particles/anti-particles
  for (Int_t ii = 0; ii < fNNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fNp[ii][jj] = 0;
}

//________________________________________________________________________
void InitializeEventStats() {
  // -- Initialize event statistics histograms
  
  fOutList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", fHEventStatMax, -0.5, fHEventStatMax-0.5));
  TH1F *hEventStat0 = static_cast<TH1F*>(fOutList->Last());

  fOutList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", fHEventStatMax, -0.5, fHEventStatMax-0.5));
  TH1F *hEventStat1 = static_cast<TH1F*>(fOutList->Last());

  for (Int_t ii=0; ii < fHEventStatMax-1; ii++) {
    hEventStat0->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
    hEventStat1->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
  }

  hEventStat0->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", aCentralityMaxNames[9-1]));
  hEventStat1->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", aCentralityMaxNames[9-1]));
}

 
//________________________________________________________________________
Bool_t FillEventStats(Int_t *aEventCuts) {
  // -- Fill event / centrality statistics 
  
  Bool_t isRejected = kFALSE;
  
  // -- Fill event statistics
  for (Int_t idx = 0; idx < fHEventStatMax ; ++idx) {
    if (aEventCuts[idx])
      isRejected = kTRUE;
    else
      (static_cast<TH1F*>(fOutList->FindObject("hEventStat0")))->Fill(idx);
  }
  
  for (Int_t idx = 0; idx < fHEventStatMax; ++idx) {
    if (aEventCuts[idx])
      break;
    (static_cast<TH1F*>(fOutList->FindObject("hEventStat1")))->Fill(idx);
  }
  
  return isRejected;
}

//________________________________________________________________________
Int_t GetCentrality(Int_t nRefMultTracksCorr) {
  // -- Fill centrality statistics of accepted events
  //    first bin = 0
  
  Int_t centrality = -1;
  
  if (nRefMultTracksCorr >= aRefMult2[0])
    centrality = 0;   // 0-5

  else if (nRefMultTracksCorr < aRefMult2[fNCentralityBins-1])
    centrality = fNCentralityBins;   // 80-100    

  else {
    for (Int_t idx = 1; idx < fNCentralityBins; ++idx) {
      if (nRefMultTracksCorr >= aRefMult2[idx] && nRefMultTracksCorr < aRefMult2[idx-1])
	centrality = idx;  // 5-10  ... 70-80
    }
  }

  return centrality;
}

//________________________________________________________________________
void InitializeMultiplicityStats() {
  // -- Initialize trigger statistics histograms

  for (Int_t ii = 0 ; ii < nMultSets ; ++ii) {

    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_multHists%s", name, multNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH1F(Form("hCentralityStat%s", multNames[ii]),               Form("Centrality statistics%s;Centrality Bins;Events",multTitles[ii]),          fNCentralityBins,-0.5,fNCentralityBins-0.5));
    TH1F* hCentralityStat = static_cast<TH1F*>(list->Last());
    
    for (Int_t jj = 0; jj < fNCentralityBins; jj++)
      hCentralityStat->GetXaxis()->SetBinLabel(jj+1, aCentralityNames[jj]);
  
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    list->Add(new TH1F(Form("hRefMultStat%s", multNames[ii]),                    Form("RefMult  Statistics%s;RefMult;Events",multTitles[ii]),                                501, 0., 500.));
    list->Add(new TH1F(Form("hRefMult2Stat%s", multNames[ii]),                   Form("RefMult2 Statistics%s;RefMult2;Events",multTitles[ii]),                               501, 0., 500.));

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH2F(Form("hRefMult2_nGlobalTracks%s", multNames[ii]),         Form("RefMult2 vs nGlobalTracks%s;RefMult2;nGlobalTracks",multTitles[ii]),                  501, 0., 500., 2501, 0., 2500.));
    list->Add(new TH2F(Form("hRefMult2_nPrimaryTracks%s", multNames[ii]),        Form("RefMult2 vs nPrimaryTracks%s;RefMult2;nPrimaryTracks",multTitles[ii]),                501, 0., 500., 1001, 0., 1000.));
    list->Add(new TH2F(Form("hRefMult2_nTOFMatch%s", multNames[ii]),             Form("RefMult2 vs nTOFMatch%s;RefMult2;nTOFMatch",multTitles[ii]),                          501, 0., 500., 601, 0., 600.));

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH2F(Form("hRefMult_nGlobalTracks%s", multNames[ii]),          Form("RefMult vs nGlobalTracks%s;RefMult;nGlobalTracks",multTitles[ii]),                    501, 0., 500., 2501, 0., 2500.));
    list->Add(new TH2F(Form("hRefMult_nPrimaryTracks%s", multNames[ii]),         Form("RefMult vs nPrimaryTracks%s;RefMult;nPrimaryTracks",multTitles[ii]),                  501, 0., 500., 1001, 0., 1000.));
    list->Add(new TH2F(Form("hRefMult_nTOFMatch%s", multNames[ii]),              Form("RefMult vs nTOFMatch%s;RefMult;nTOFMatch",multTitles[ii]),                            501, 0., 500., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hRefMult2_RefMult2Corr%s", multNames[ii]),          Form("RefMult2 vs RefMult2Corr%s;RefMult2;RefMult2Corr",multTitles[ii]),                    501, 0., 500., 501, 0., 500.));
    list->Add(new TH2F(Form("hRefMult2_RefMult%s", multNames[ii]),               Form("RefMult2 vs RefMult%s;RefMult2;RefMult",multTitles[ii]),                              501, 0., 500., 501, 0., 500.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hnPrimaryTracks_nGlobalTracks%s", multNames[ii]),   Form("nPrimaryTracks vs nGlobalTracks%s;nPrimaryTracks;nGlobalTracks",multTitles[ii]),     1001, 0., 1000., 2501, 0., 2500.));
    list->Add(new TH2F(Form("hnPrimaryTracks_nTOFMatch%s", multNames[ii]),       Form("nPrimaryTracks vs nTOFMatch%s; nPrimaryTracks;nTOFMatch",multTitles[ii]),            1001, 0., 1000., 601, 0., 600.));    
    list->Add(new TH2F(Form("hnGlobalTracks_nTOFMatch%s", multNames[ii]),        Form("nGlobalTracks vs nTOFMatch%s;nGlobalTracks;nTOFMatch",multTitles[ii]),               2501, 0., 2500., 601, 0., 600.));
      
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    if ( ii == nMultSets-1 ) {
      list->Add(new TProfile("hRefMult2Mean",            "RefMult2 Mean vs Centrality;Centrality Bins;<RefMult2>",         fNCentralityBins,-0.5,fNCentralityBins-0.5));
      TH1F* hRefMult2Mean = static_cast<TH1F*>(list->Last());
      
      list->Add(new TProfile("hRefMult2CorrMean",        "RefMult2Corr Mean vs Centrality;Centrality Bins;<RefMult2Corr>", fNCentralityBins,-0.5,fNCentralityBins-0.5));
      TH1F* hRefMult2CorrMean = static_cast<TH1F*>(list->Last());
      
      // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      
      for (Int_t jj=0; jj < fNCentralityBins; jj++) {
	hRefMult2Mean->GetXaxis()->SetBinLabel(jj+1, aCentralityNames[jj]);
	hRefMult2CorrMean->GetXaxis()->SetBinLabel(jj+1, aCentralityNames[jj]);
      }
    }
  }    
}
 
//________________________________________________________________________
 void FillMultiplicityStats(Double_t *aMult, Int_t mode) {
  // -- Fill centrality statistics of accepted events

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_multHists%s", name, multNames[mode])));
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  (static_cast<TH1F*>(list->FindObject(Form("hCentralityStat%s",               multNames[mode]))))->Fill(aMult[0]);

  (static_cast<TH1F*>(list->FindObject(Form("hRefMultStat%s",                  multNames[mode]))))->Fill(aMult[1]);  
  (static_cast<TH1F*>(list->FindObject(Form("hRefMult2Stat%s",                 multNames[mode]))))->Fill(aMult[2]);

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  (static_cast<TH2F*>(list->FindObject(Form("hRefMult2_nGlobalTracks%s",       multNames[mode]))))->Fill(aMult[2], aMult[4]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult2_nPrimaryTracks%s",      multNames[mode]))))->Fill(aMult[2], aMult[5]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult2_nTOFMatch%s",           multNames[mode]))))->Fill(aMult[2], aMult[6]);

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  (static_cast<TH2F*>(list->FindObject(Form("hRefMult_nGlobalTracks%s",        multNames[mode]))))->Fill(aMult[1], aMult[4]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult_nPrimaryTracks%s",       multNames[mode]))))->Fill(aMult[1], aMult[5]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult_nTOFMatch%s",            multNames[mode]))))->Fill(aMult[1], aMult[6]);
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult2_RefMult2Corr%s",        multNames[mode]))))->Fill(aMult[2], aMult[3]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult2_RefMult%s",             multNames[mode]))))->Fill(aMult[2], aMult[1]);
  
  (static_cast<TH2F*>(list->FindObject(Form("hnPrimaryTracks_nGlobalTracks%s", multNames[mode]))))->Fill(aMult[5], aMult[4]);

  (static_cast<TH2F*>(list->FindObject(Form("hnPrimaryTracks_nTOFMatch%s",     multNames[mode]))))->Fill(aMult[5], aMult[6]);
  (static_cast<TH2F*>(list->FindObject(Form("hnGlobalTracks_nTOFMatch%s",      multNames[mode]))))->Fill(aMult[4], aMult[6]);
  
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  if (mode == nMultSets-1) {
    (static_cast<TProfile*>(list->FindObject("hRefMult2Mean")))->Fill(aMult[0], aMult[2]);
    (static_cast<TProfile*>(list->FindObject("hRefMult2CorrMean")))->Fill(aMult[0], aMult[3]);
  }
}

//________________________________________________________________________
void InitializeEventHists() {
  // -- Initialize trigger statistics histograms

  const Char_t *eventNames[]  = { "before", "after" };
  const Char_t *eventTitles[] = { " (before cuts)", " (after cuts)" };

  for (Int_t ii = 0 ; ii < 2 ; ++ii) {

    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_eventHists_%s", name, eventNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH1F(Form("vx_%s", eventNames[ii]),        Form("#it{v}_{x}%s;#it{v}_{x} (cm);Events", eventTitles[ii]),               binHnEvent[1], minHnEvent[1], maxHnEvent[1]));
    list->Add(new TH1F(Form("vy_%s", eventNames[ii]),        Form("#it{v}_{y}%s;#it{v}_{y} (cm);Events", eventTitles[ii]),               binHnEvent[2], minHnEvent[2], maxHnEvent[2]));
    list->Add(new TH1F(Form("vz_%s", eventNames[ii]),        Form("#it{v}_{z}%s;#it{v}_{z} (cm);Events", eventTitles[ii]),               binHnEvent[3], minHnEvent[3], maxHnEvent[3]));
    list->Add(new TH1F(Form("shiftedVr_%s", eventNames[ii]), Form("#it{v}_{r}^{shifted}%s;#it{v}_{r}^{shifted} (cm);Events", eventTitles[ii]), binHnEvent[4], minHnEvent[4], maxHnEvent[4]));
    list->Add(new TH2F(Form("vxvy_%s", eventNames[ii]),      Form("#it{v}_{x} vs #it{v}_{y}%s;#it{v}_{x} (cm); #it{v}_{y} (cm)", eventTitles[ii]), binHnEvent[1], minHnEvent[1], maxHnEvent[1], binHnEvent[2], minHnEvent[2], maxHnEvent[2]));

    list->Add(new TH1F(Form("vzVpd_%s", eventNames[ii]),     Form("#it{v}_{z}^{vpd}%s;#it{v}_{z}^{vpd} (cm);Events", eventTitles[ii]),         binHnEvent[9], minHnEvent[9], maxHnEvent[9]));

    list->Add(new TH1F(Form("deltaVz_%s", eventNames[ii]),   Form("#Delta#it{v}_{z}%s;#Delta#it{v}_{z} (cm);Events", eventTitles[ii]),         binHnEvent[10], minHnEvent[10], maxHnEvent[10]));
  }
}

//________________________________________________________________________
void FillEventHists(Int_t runId, Double_t *aEvent, Int_t mode) {
      
  const Char_t *eventNames[] = { "before", "after" };

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_eventHists_%s", name, eventNames[mode])));
  
  (static_cast<TH1F*>(list->FindObject(Form("vx_%s",eventNames[mode]))))->Fill(aEvent[1]);
  (static_cast<TH1F*>(list->FindObject(Form("vy_%s",eventNames[mode]))))->Fill(aEvent[2]);
  (static_cast<TH1F*>(list->FindObject(Form("vz_%s",eventNames[mode]))))->Fill(aEvent[3]);
  (static_cast<TH1F*>(list->FindObject(Form("shiftedVr_%s",eventNames[mode]))))->Fill(aEvent[3]);
  (static_cast<TH2F*>(list->FindObject(Form("vxvy_%s",eventNames[mode]))))->Fill(aEvent[1], aEvent[2]);
  (static_cast<TH1F*>(list->FindObject(Form("vzVpd_%s",eventNames[mode]))))->Fill(aEvent[9]);
  (static_cast<TH1F*>(list->FindObject(Form("deltaVz_%s",eventNames[mode]))))->Fill(aEvent[10]);
}

//________________________________________________________________________
void InitializeTrackHists() {
  // -- Initialize trigger statistics histograms

  const Char_t *trackNames[]  = { "before", "after", "final" };
  const Char_t *trackTitles[] = { " (before cuts)", " (after cuts)",  " (after final cuts)" };

  for (Int_t ii = 0 ; ii < 3 ; ++ii) {

    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_trackHists_%s", name, trackNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
 
    list->Add(new TH1F(Form("pt_%s", trackNames[ii]),                Form("#it{p}_{T}%s;#it{p}_{T} (GeV/#it{c});Tracks", trackTitles[ii]),  binHnUnCorr[1], minHnUnCorr[1], maxHnUnCorr[1]));
    list->Add(new TH1F(Form("eta_%s", trackNames[ii]),               Form("#eta%s;#eta;Tracks", trackTitles[ii]),  binHnUnCorr[2], minHnUnCorr[2], maxHnUnCorr[2]));
    list->Add(new TH1F(Form("dca_%s", trackNames[ii]),               Form("DCA%s;DCA (cm);Tracks", trackTitles[ii]),  binHnUnCorr[4], minHnUnCorr[4], maxHnUnCorr[4]));
    list->Add(new TH1F(Form("nHitsDedx_%s", trackNames[ii]),         Form("nHitsDedx%s;nHitsDedx;Tracks", trackTitles[ii]),  binHnUnCorr[5], minHnUnCorr[5], maxHnUnCorr[5]));
    list->Add(new TH1F(Form("nHitsFit_%s", trackNames[ii]),          Form("nHitsFit%s;nHitsFit;Tracks", trackTitles[ii]),  binHnUnCorr[6], minHnUnCorr[6], maxHnUnCorr[6]));
    list->Add(new TH1F(Form("nHitsFit_nFitPoss_%s", trackNames[ii]), Form("nHitsFit/nFitPoss%s;nHitsFit/nFitPoss;Tracks", trackTitles[ii]),  binHnUnCorr[8], minHnUnCorr[8], maxHnUnCorr[8]));
    list->Add(new TH1F(Form("nSigmaP_%s", trackNames[ii]),           Form("nSigmaProton%s;n#Sigma_{proton};Tracks", trackTitles[ii]),  binHnUnCorr[11], minHnUnCorr[11], maxHnUnCorr[11]));
    list->Add(new TH2F(Form("nSigmaP_pt_%s", trackNames[ii]),        Form("nSigmaProton vs #it{p}_{T}%s;n#sigma_{proton};#it{p}_{T} (GeV/#it{c})", trackTitles[ii]),  
		       binHnUnCorr[11], minHnUnCorr[11], maxHnUnCorr[11], binHnUnCorr[1], minHnUnCorr[1], maxHnUnCorr[1]));
  }
}

//________________________________________________________________________
void FillTrackHists(Double_t *aTrack, Int_t mode) {
      
  const Char_t *trackNames[] = { "before", "after", "final" };

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_trackHists_%s", name, trackNames[mode])));

  (static_cast<TH1F*>(list->FindObject(Form("pt_%s", trackNames[mode]))))->Fill(aTrack[1]);
  (static_cast<TH1F*>(list->FindObject(Form("eta_%s", trackNames[mode]))))->Fill(aTrack[2]);
  (static_cast<TH1F*>(list->FindObject(Form("dca_%s", trackNames[mode]))))->Fill(aTrack[4]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsDedx_%s", trackNames[mode]))))->Fill(aTrack[5]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsFit_%s", trackNames[mode]))))->Fill(aTrack[6]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsFit_nFitPoss_%s", trackNames[mode]))))->Fill(aTrack[8]);
  (static_cast<TH1F*>(list->FindObject(Form("nSigmaP_%s", trackNames[mode]))))->Fill(aTrack[11]);
  (static_cast<TH2F*>(list->FindObject(Form("nSigmaP_pt_%s", trackNames[mode]))))->Fill(aTrack[11], aTrack[1]);
}

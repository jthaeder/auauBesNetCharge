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
#define DEBUG 0

void   AddHistSetCent(const Char_t *name, const Char_t *title);
void   FillHistSetCent(const Char_t *name, Int_t idx, Int_t cent);

void   Reset();

void   InitializeEventStats();
void   InitializeMultiplicityStats();
void   InitializeEventHists();
void   InitializeTrackHists();

void   InitializeRunByRunEventHists();
void   InitializeRunByRunTrackHists();

Int_t  GetCentrality(Int_t nRefMultTracksCorr, Int_t anaIdx);

Bool_t FillEventStats(Int_t *aEventCuts);
void   FillMultiplicityStats(Double_t *aMult, Int_t mode);
void   FillEventHists(Double_t *aEvent, Int_t mode);
void   FillTrackHists(Double_t *aTrack, Int_t mode);

void   FillRunByRunEventHists(Double_t *aEvent, Int_t mode, Int_t runIdx);
void   FillRunByRunTrackHists(Double_t *aTrack, Int_t mode, Int_t runIdx);

Int_t  RefMultCorrection(Double_t vz, Int_t refmult2, Int_t anaIdx);

Double_t NN(Double_t, Int_t);

// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------  

const Double_t masses[3] = {0.1349764, 0.493677, 0.9382723128};

Int_t    binHnEvent[11] = {  10,   41,   41,    401,   21,   601,   601,   3001,    2,    201,   21};
Double_t minHnEvent[11] = {-0.5, -2.0, -3.0, -100.0,  0.0,   0.0,   0.0,    0.0, -0.5, -100.0,  0.0};
Double_t maxHnEvent[11] = { 9.5,  2.0,  1.0,  100.0,  2.0, 600.0, 600.0, 3000.0,  1.5,  100.0, 20.0};

Int_t    binHnUnCorr[12] = {  10,  34,   21,    3,  50,   51,   51,   51, 101,    2,    2, 81};
Double_t minHnUnCorr[12] = {-0.5, 0.1, -1.0, -1.5,   0,  0.0,  0.0,  0.0, 0.0, -0.5, -0.5, -4};
Double_t maxHnUnCorr[12] = { 9.5, 3.0,  1.0,  1.5,   5, 50.0, 50.0, 50.0, 1.0,  1.5,  1.5,  4};

Int_t nMultSets = 5;

const Char_t *qaNames[]  = { "before", "after"};
const Char_t *qaTitles[] = { " (before cuts)", " (after cuts)"};

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

const Char_t* name[3]  = {"NetCharge", "NetProton", "NetKaon"};
const Char_t* name2[3] = {"charge", "proton", "kaon"};

const Char_t* nameRefMult[3]      = {"refMult2", "refMult3", "refMult4"};
const Char_t* nameRefMultShort[3] = {"refmult2", "refmult3", "refmult4"};

enum particleCharge {kPOS, kNEG, kNET, kParticleCharge};

const Char_t* asParticleName[2]  =  {"neg",  "pos" };
const Char_t* asParticleTitle[2] =  {"neg.", "pos."};

const Int_t  fHEventStatMax   = 10; 
const Char_t *aEventCutNames[]   = {"all", "bad run", "trigger", "#it{v}_{z} < #it{v}_{z}^{max}", 
				    "#it{v}_{z}- #it{v}_{z}^{vpd} < 3cm", "shifted vtx <1", 
				    "centrality", "nTOFMatch>2","nTOFMatch cut", "accepted"};

double randomEff[2][9];

Int_t nRuns = 484;
Int_t runs[] = {15053001, 15053003, 15053004, 15053005, 15053007, 15053008, 15053009, 15053011, 15053012, 15053015, 15053016, 15053017, 15053019, 15053020, 15053021, 15053022, 15053023, 15053024, 15053025, 15053026, 15053047, 15053059, 15053060, 15053062, 15053064, 15053065, 15053067, 15054001, 15054002, 15054003, 15054004, 15054005, 15054006, 15054007, 15054008, 15054009, 15054010, 15054011, 15054012, 15054013, 15054014, 15054015, 15054016, 15054017, 15054018, 15054019, 15054020, 15054021, 15054023, 15054024, 15054025, 15054026, 15054028, 15054029, 15054030, 15054031, 15054037, 15054042, 15054043, 15054044, 15054046, 15054047, 15054048, 15054049, 15054050, 15054051, 15054052, 15055001, 15055002, 15055003, 15055004, 15055005, 15055006, 15055007, 15055008, 15055009, 15055011, 15055012, 15055013, 15055014, 15055015, 15055016, 15055017, 15055019, 15055020, 15055021, 15055135, 15055136, 15055138, 15055139, 15055140, 15055141, 15056001, 15056002, 15056003, 15056004, 15056005, 15056006, 15056007, 15056008, 15056009, 15056013, 15056014, 15056015, 15056016, 15056017, 15056018, 15056019, 15056020, 15056021, 15056022, 15056023, 15056024, 15056025, 15056026, 15056027, 15056028, 15056036, 15056037, 15056038, 15056039, 15058051, 15058052, 15058053, 15058054, 15058055, 15059001, 15059002, 15059003, 15059005, 15059006, 15059007, 15059009, 15059010, 15059011, 15059013, 15059014, 15059015, 15059016, 15059017, 15059018, 15059019, 15059020, 15059021, 15059024, 15059025, 15059026, 15059027, 15059028, 15059029, 15059033, 15059037, 15059038, 15059040, 15059041, 15059042, 15059055, 15059059, 15059060, 15059061, 15059063, 15059064, 15059068, 15059074, 15059076, 15059081, 15059082, 15059086, 15059087, 15059090, 15060001, 15060006, 15060007, 15060009, 15060011, 15060012, 15060014, 15060015, 15060016, 15060017, 15060018, 15060019, 15060020, 15060021, 15060022, 15060023, 15060024, 15060027, 15060028, 15060029, 15060031, 15060032, 15060033, 15060035, 15060036, 15060037, 15060044, 15060045, 15060046, 15060047, 15060048, 15060049, 15060050, 15060051, 15060052, 15060053, 15060067, 15060068, 15060069, 15060070, 15060071, 15061003, 15061004, 15061006, 15061007, 15061008, 15061010, 15061011, 15061012, 15061014, 15061015, 15061016, 15061018, 15061019, 15061021, 15061023, 15061024, 15061025, 15061026, 15061027, 15061028, 15061034, 15061036, 15061037, 15061038, 15061039, 15061041, 15061047, 15061048, 15061051, 15061053, 15061054, 15061055, 15061056, 15061059, 15061060, 15061061, 15061062, 15061063, 15061064, 15062003, 15062004, 15062005, 15062007, 15062008, 15062009, 15062010, 15062011, 15062012, 15062013, 15062015, 15062017, 15062018, 15062020, 15062022, 15062023, 15062024, 15062025, 15062026, 15062031, 15062032, 15062033, 15062034, 15062035, 15062036, 15062037, 15062038, 15062040, 15062041, 15062042, 15062066, 15062070, 15062071, 15062075, 15062076, 15062077, 15063001, 15063002, 15063003, 15063006, 15063008, 15063010, 15063011, 15063012, 15063013, 15063014, 15063016, 15063018, 15063020, 15063021, 15063029, 15063032, 15063036, 15063037, 15063040, 15063041, 15063042, 15063043, 15063045, 15063046, 15063048, 15063049, 15063053, 15063056, 15063059, 15063060, 15063061, 15063062, 15063063, 15063065, 15063066, 15063067, 15064001, 15064002, 15064003, 15064005, 15064006, 15064007, 15064008, 15064009, 15064010, 15064011, 15065011, 15065013, 15065017, 15065054, 15065056, 15065059, 15065060, 15065061, 15066005, 15066006, 15066010, 15066012, 15066014, 15066016, 15066019, 15066020, 15066021, 15066022, 15066023, 15066024, 15066025, 15066026, 15066064, 15066076, 15066077, 15066083, 15066085, 15066086, 15066088, 15067001, 15067002, 15067003, 15067005, 15067006, 15067008, 15067009, 15067011, 15067012, 15067013, 15067014, 15067015, 15067016, 15067018, 15067019, 15067020, 15067022, 15067023, 15067024, 15067025, 15067026, 15067032, 15067033, 15067034, 15067035, 15067036, 15067037, 15067038, 15067040, 15067041, 15067042, 15067043, 15067044, 15067046, 15067048, 15068001, 15068002, 15068003, 15068004, 15068005, 15068006, 15068007, 15068008, 15068009, 15068010, 15068012, 15068022, 15068023, 15068024, 15068025, 15068026, 15068027, 15068028, 15068031, 15068032, 15068033, 15068034, 15068035, 15068036, 15068037, 15068038, 15068039, 15068040, 15068041, 15068042, 15068043, 15068044, 15068045, 15068046, 15068047, 15068048, 15068049, 15069001, 15069002, 15069003, 15069004, 15069005, 15069006, 15069007, 15069008, 15069009, 15069010, 15069011, 15069012, 15069013, 15069014, 15069015, 15069016, 15069017, 15069018, 15069019, 15069020, 15069021, 15069022, 15069023, 15069030, 15069031, 15069034, 15069035, 15069038, 15069039, 15069040, 15069042, 15069043, 15069044, 15069045, 15069046, 15069047, 15069048, 15069049, 15069050, 15069051, 15070001, 15070002, 15070013, 15070014, 15070015, 15070016, 15070017, 15070018, 15070019, 15070020, 15070021} ;


// ----------------------------------------------------------------------------  
// -- Globals
// ----------------------------------------------------------------------------  

TList                *fOutList;                     //! Output data container
// =======================================================================
Int_t                 fOrder = 8;                   //  Max order of higher order distributions
// -----------------------------------------------------------------------
Int_t                 fNNp   = 3;                   //  N sets of arrays of particle/anti-particle counts
                                                    //   0 is all 
                                                    //   1,2 for arbitrary subset
Int_t               **fNp;                          //  Array of particle/anti-particle counts
// =======================================================================
THnSparseD           *fHnTrackUnCorr;               //  THnSparseD : uncorrected probe particles
THnSparseD           *fHnEvent;                     //  THnSparseD : event
// =======================================================================
Int_t                 energyIdx   = -1;             // energyIdx
// =======================================================================
Int_t                 analysisIdx = -1;             // 0 - net-charge
                                                    // 1 - net-proton
                                                    // 2 - net-kaon
// =======================================================================
Int_t                 useModeChargeSeparation = 0;  // 0 - off
                                                    // 1 - positive charge -> positive eta
                                                    // 2 - positive charge -> negative eta

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

      if (argument.BeginsWith("--analysis=")) {
	TString parameter = argument.ReplaceAll("--analysis=", "");
	parameter.Remove(TString::kLeading, ' '); 
	parameter.ToLower();
	if (!parameter.CompareTo("netcharge") || !parameter.CompareTo("net-charge") || !parameter.CompareTo("charge") )
	  analysisIdx = 0;
	else if (!parameter.CompareTo("netproton") || !parameter.CompareTo("net-proton") || !parameter.CompareTo("proton") )
	  analysisIdx = 1;
	else if (!parameter.CompareTo("netkaon") || !parameter.CompareTo("net-kaon") || !parameter.CompareTo("kaon") )
	  analysisIdx = 2;
      } // if (argument.BeginsWith("--analysis=")) {

      else if (argument.BeginsWith("--energy=")) {
	TString parameter = argument.ReplaceAll("--energy=", "");
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

      else if (argument.BeginsWith("--chargeSeparation=")) {
	TString parameter = argument.ReplaceAll("--chargeSeparation=", "");
	parameter.Remove(TString::kLeading, ' '); 
	useModeChargeSeparation = parameter.Atoi();
      } // else if (argument.BeginsWith("--chargeSeparation=")) {

      else if (argument.BeginsWith("--etaMin=")) {
	if (analysisIdx == -1) {
	  cout << "No analysis type selected, before eta Arguments" << endl;
	  return 0;
	}
	TString parameter = argument.ReplaceAll("--etaMin=", "");
	parameter.Remove(TString::kLeading, ' '); 
	etaAbsRange[analysisIdx][0] = parameter.Atof();
      } // else if (argument.BeginsWith("--etaMin=")) {

      else if (argument.BeginsWith("--etaMax=")) {
	if (analysisIdx == -1) {
	  cout << "No analysis type selected, before eta Arguments" << endl;
	  return 0;
	}
	TString parameter = argument.ReplaceAll("--etaMax=", "");
	parameter.Remove(TString::kLeading, ' '); 
	etaAbsRange[analysisIdx][1] = parameter.Atof();
      } // else if (argument.BeginsWith("--etaMax=")) {

    } // for (Int_t idx = 1; idx < argc; ++idx) {
  } // if (argc > 1) {

  if (energyIdx == -1) {
    cout << "No energy selected" << endl;
    return 0;
  }

  if (analysisIdx == -1) {
    cout << "No analysis type selected" << endl;
    return 0;
  }
  
  // -----------------------------------------------------------------------
  // -- Set Globals
  // -----------------------------------------------------------------------

  // -- Add List
  fOutList = new TList;
  fOutList->SetName(Form("f%s", name[analysisIdx]));
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
#if DEBUG
    cout << "neg " <<  idxCent << "  " << randomEff[0][idxCent] << " -> " << prNegSub->GetMean() << " / " << prNegAll->GetMean() <<  endl;
    cout << "pos " <<  idxCent << "  " << randomEff[1][idxCent] << " -> " << prPosSub->GetMean() << " / " << prPosAll->GetMean() <<  endl;
#endif
  } // for (int idxCent = 0; idxCent <9 ; ++idxCent) {
#endif
  // ------------------------------------------------------------------
  // -- Get event container
  // ------------------------------------------------------------------

#if EVENT_THN
  // -- Event
  fOutList->Add(new THnSparseD("hnEvent", "cent:vx:vy:vz:shiftedVtx:nRefMultX:nRefMultXCorr:nTracks:isRejected", 9, 
			       binHnEvent, minHnEvent, maxHnEvent));  

  fHnEvent = static_cast<THnSparseD*>(fOutList->Last());
  fHnEvent->Sumw2(); 
  fHnEvent->GetAxis(0)->SetTitle("centrality");
  fHnEvent->GetAxis(1)->SetTitle("#it{v}_{x} (cm)");
  fHnEvent->GetAxis(2)->SetTitle("#it{v}_{y} (cm)");
  fHnEvent->GetAxis(3)->SetTitle("#it{v}_{z} (cm)");
  fHnEvent->GetAxis(4)->SetTitle("shifted vtx (cm)");
  fHnEvent->GetAxis(5)->SetTitle("refMultX");
  fHnEvent->GetAxis(6)->SetTitle("refMultXCorr");
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
  StRefMultCorr* refmultCorr = new StRefMultCorr(nameRefMultShort[analysisIdx]);

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
	
	refMultMap[std::make_pair(run,event)]  = refMult;
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
  InitializeRunByRunEventHists();
  InitializeRunByRunTrackHists();

  TString sTitle("");
  sTitle += (!analysisIdx) ? Form("%.2f<#eta<%.2f", etaAbsRange[analysisIdx][0], etaAbsRange[analysisIdx][1]) : 
    Form("%.2f<y<%.2f", yAbsRange[analysisIdx][0], yAbsRange[analysisIdx][1]); 		   
  sTitle += Form(" #it{p}_{T} [%.1f,%.1f]", ptRange[analysisIdx][0], ptRange[analysisIdx][1]);
  
  cout << name[analysisIdx] << ": .... " << sTitle << " useChargeSeparation "<<  useModeChargeSeparation << endl;

  AddHistSetCent("Dist",       sTitle.Data());
  AddHistSetCent("Dist_lower", sTitle.Data());
  AddHistSetCent("Dist_upper", sTitle.Data());
  
  //*************************************************************************************************

  fOutList->Add(new TList);
  TList *fijkhList = static_cast<TList*>(fOutList->Last());
  fijkhList->SetName(Form("f%sFijkh", name[analysisIdx]));
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
    Int_t runIdx = 0;
    for (Int_t ii = 0; ii < nRuns; ++ii) {
      if (runId == runs[ii]) {
	runIdx = ii;
	break;
      }
    }
    
    Int_t eventId = pico->Event_mEventId[0];

    // ------------------------------------------------------------------
    // -- Track loop - for refMultX / TOF tracks - track multiplicity
    // ------------------------------------------------------------------
    Int_t nTracks         = pico->Tracks_;      
    Int_t nGlobalTracks   = pico->Event_mNumberOfGlobalTracks[0]; 
    Int_t nPrimaryTracks  = 0;
    
    Int_t nRefMultTracksPico  = pico->Event_mRefMultNeg[0] + pico->Event_mRefMultPos[0];
    Int_t nRefMultTracksMuDst = refMultMap[std::make_pair(runId, eventId)];

    Int_t nRefMultTracks      = (energyIdx == 2) ? nRefMultTracksPico : nRefMultTracksMuDst;
#if DEBUG
    // --- this will be used in after picoDST reproduction 
    Int_t nRefMult2TracksPico  = pico->Event_mRefMult2NegEast[0] + pico->Event_mRefMult2PosEast[0] +
      pico->Event_mRefMult2NegWest[0] +  pico->Event_mRefMult2PosWest[0];
#endif
    Int_t nRefMult2TracksMuDst = refMult2Map[std::make_pair(runId, eventId)];

    Int_t nRefMultXTracks[3], nRefMultXTracksCorr[3];
    nRefMultXTracks[0] = (energyIdx == 2) ? 0 : nRefMult2TracksMuDst;
    nRefMultXTracks[1] = (energyIdx == 2) ? 0 : nRefMult2TracksMuDst; // read out refMult3 from MUDSTS....
    nRefMultXTracks[2] = (energyIdx == 2) ? 0 : nRefMult2TracksMuDst; // read out refMult4 from MUDSTS....

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
      
      Double_t nSigmaProton = pico->Tracks_mNSigmaProton[idxTrack]/100.;
      Double_t nSigmaKaon   = pico->Tracks_mNSigmaKaon[idxTrack]/100.;
            
      Double_t beta         = pico->Tracks_mBTofBeta[idxTrack]/20000.;
      beta = (beta <= 0) ? -999 : beta;
            
      Double_t mSquare = (beta == -999. || beta <= 1.e-5) ? -999. : pow(pcm,2)*(pow(1/beta,2)-1);

      // -- count for refMult2 track -- only for 14.5 GeV
      // ------------------------------------------------------------------
      if (energyIdx == 2 
	  && TMath::Abs(eta) => etaAbsRangeRefMult[0][0] 
	  && TMath::Abs(eta) <  etaAbsRangeRefMult[0][1] 
	  && TMath::Abs(nHitsFit) > nHitsFitRefMult[0]  
	  && DCA < dcaMaxRefMult[0])
	++nRefMultXTracks[0];

      // -- count for refMult3 track -- only for 14.5 GeV
      // ------------------------------------------------------------------
      if (energyIdx == 2 
	  && TMath::Abs(eta) => etaAbsRangeRefMult[1][0] 
	  && TMath::Abs(eta) < etaAbsRangeRefMult[1][1] 
	  && TMath::Abs(nHitsFit) > nHitsFitRefMult[1]  
	  && DCA < dcaMaxRefMult[1]
	  && nSigmaProton < (-3) 
	  && mSquare < 0.4)
	++nRefMultXTracks[1];

      // -- count for refMult4 track -- only for 14.5 GeV
      // ------------------------------------------------------------------
      if (energyIdx == 2 
	  && TMath::Abs(eta) => etaAbsRangeRefMult[2][0] 
	  && TMath::Abs(eta) < etaAbsRangeRefMult[2][1] 
	  && TMath::Abs(nHitsFit) > nHitsFitRefMult[2]  
	  && DCA < dcaMaxRefMult[1]
	  && ( (mSquare == -999 && fabs(nSigmaKaon) > 3) || (mSquare != -999 && (mSquare > 0.6 || mSquare < 0.1))) )
	++nRefMultXTracks[2];

      // -- count for primary track
      // ------------------------------------------------------------------
      if (! (pxcm == 0 && pycm == 0 && pzcm == 0) )
	++nPrimaryTracks;
    }
    
#if DEBUG
    if (nRefMultTracksPico - nRefMultTracksMuDst != 0) 
      cout << "  DELTA REFMULT " << nRefMultTracks  << " || Pico:" << nRefMultTracksPico  << " <> MuDst:" << nRefMultTracksMuDst  << endl;
    cout << "  REFMULTX      " << nRefMultXTracks[analysisIdx] << " || Pico:" << nRefMult2TracksPico << " <> MuDst:" << nRefMult2TracksMuDst << endl;
#endif 
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
    if (refmultCorr->isBadRun(runId)) 
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
    Float_t  vpdVz      = Float_t(pico->Event_mVzVpd[0])/100;
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

    Int_t centrality;

    // -- 14.5 GeV
    if (energyIdx == 2) {
      // -- Correct refmult vz dependent 
      for (Int_t ii = 0; ii < 3; ++ii)
	nRefMultXTracksCorr[ii] = RefMultCorrection(vz, nRefMultXTracks[ii], ii);
      
      // -- Get centrality and fill stats - for p
      centrality = GetCentrality(nRefMultXTracksCorr[analysisIdx], analysisIdx);
    }
    // -- other energies
    else {
      refmultCorr->init(runId);
      refmultCorr->initEvent(nRefMultXTracks[analysisIdx], vz);

      // -- Correct refmult vz dependent 
      nRefMultXTracksCorr[analysisIdx] = refmultCorr->getRefMultCorr();

      centrality = 8 - refmultCorr->getCentralityBin9();
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
			   Double_t(nRefMultXTracks[analysisIdx]), Double_t(nRefMultXTracksCorr[analysisIdx]), Double_t(nTracks), Double_t(isRejected),
			   vpdVz, deltaVz};
      
    Double_t aMult[7]  = {Double_t(centrality), Double_t(nRefMultTracks), 
			  Double_t(nRefMultXTracks[analysisIdx]), Double_t(nRefMultXTracksCorr[analysisIdx]), 
			  Double_t(nGlobalTracks), Double_t(nPrimaryTracks),
			  Double_t(nTOFMatch) };

#if EVENT_THN
    fHnEvent->Fill(aEvent);
#endif

    // -- Fill eventHists
    FillEventHists(aEvent, 0);
    
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

#if DEBUG
    if ( nRefMultXTracks[analysisIdx] != nRefMult2TracksPico && TMath::Abs(nRefMultXTracks[analysisIdx] - nRefMult2TracksPico) > 1)
      cout << " nRefMultXTracks " << nRefMultXTracks[analysisIdx] << " -> (pico) " << nRefMult2TracksPico << endl;
#endif 

    // -- Fill multiplicty stats
    FillMultiplicityStats(aMult, 4);

    // -- Fill eventHists
    FillEventHists(aEvent, 1);

    // ------------------------------------------------------------------
    // -- Track loop - track multiplicity
    // ------------------------------------------------------------------
    for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
      
      Float_t pxcm    = pico->Tracks_mPMomentum_mX1[idxTrack];
      Float_t pycm    = pico->Tracks_mPMomentum_mX2[idxTrack];
      Float_t pzcm    = pico->Tracks_mPMomentum_mX3[idxTrack];
      if (pxcm == 0 && pycm == 0 && pzcm == 0) 
	continue;
      
      Float_t pt      = TMath::Sqrt(pxcm*pxcm + pycm*pycm);
      Float_t pcm     = TMath::Sqrt(pt*pt + pzcm*pzcm);
      Float_t ecm     = TMath::Sqrt(pcm*pcm + masses[analysisIdx]*masses[analysisIdx]);
      
      Float_t eta     = 0.5 * TMath::Log((pcm+pzcm)/(pcm-pzcm));
      Float_t DCA     = pico->Tracks_mGDca[idxTrack]/1000.;

      Float_t y       = 0.5 * TMath::Log((ecm+pzcm)/(ecm-pzcm));

      Int_t nHitsDedx = pico->Tracks_mNHitsDedx[idxTrack];
      Int_t nHitsFit  = pico->Tracks_mNHitsFit[idxTrack];
      Int_t nFitPoss  = pico->Tracks_mNHitsMax[idxTrack];
      
      Float_t ratio   =  (1+fabs(nHitsFit))/(1+nFitPoss);
      
      Float_t sign    = (nHitsFit > 0) ? +1 : -1;
      
      Float_t nSigma[3];
      nSigma[0]       = pico->Tracks_mNSigmaProton[idxTrack]/100.;
      nSigma[1]       = pico->Tracks_mNSigmaProton[idxTrack]/100.;
      nSigma[2]       = pico->Tracks_mNSigmaKaon[idxTrack]/100.;
           
      Float_t beta    = pico->Tracks_mBTofBeta[idxTrack]/20000.;
      beta = (beta <= 0) ? -999 : beta;
            
      Float_t mSquare = (beta == -999. || beta <= 1.e-5) ? -999. : pow(pcm,2)*(pow(1/beta,2)-1);

      // -- is in RefMult flag
      Bool_t isInRefMult = (TMath::Abs(eta) > etaAbsRangeRefMult[analysisIdx][0] 
			    && TMath::Abs(eta) <= etaAbsRangeRefMult[analysisIdx][1]  
			    && TMath::Abs(nHitsFit) > nHitsFitRefMult[analysisIdx]  
			    && DCA < dcaMaxRefMult[analysisIdx]) ? kTRUE : kFALSE;
      
      // -- is track accepted flag - kinematics
      Bool_t isTrackAcceptedKin = (pt > ptRange[analysisIdx][0] && pt < ptRange[analysisIdx][1]) ? kTRUE : kFALSE;

      if (isTrackAcceptedKin && analysisIdx == 0) 
	isTrackAcceptedKin = (eta > etaAbsRange[analysisIdx][0] && eta < etaAbsRange[analysisIdx][1]) ? kTRUE : kFALSE;
      else if (isTrackAcceptedKin && analysisIdx > 0) 
	isTrackAcceptedKin = (y > yAbsRange[analysisIdx][0] && y < yAbsRange[analysisIdx][1]) ? kTRUE : kFALSE;

      
      // -- is track accepted flag - clusters/dca
      Bool_t isTrackAcceptedCut = (TMath::Abs(nHitsFit) > nHitsFitMin[analysisIdx] 
				   && DCA < dcaMax[analysisIdx] 
				   && nHitsDedx > nHitsDedxMin[analysisIdx] 
				   && ratio > ratioNHitsFitNFitPossMin[analysisIdx]) ? kTRUE : kFALSE;
      
      // -- is track accepted flag - PID
      Bool_t isTrackAcceptedPid = kTRUE;
      
      if (analysisIdx > 0) {
	// -- PID for net-proton . net-kaon
	Bool_t isTrackAcceptedPidTPC = (TMath::Abs(nSigma[analysisIdx]) < nSigmaMax[analysisIdx]) ? kTRUE : kFALSE;
	Bool_t isTrackAcceptedPidTOF = kTRUE;

	if (pt > ptMidPoint[analysisIdx])
	  isTrackAcceptedPidTOF = (mSquare > mSquareRange[analysisIdx][0] 
				   && mSquare < mSquareRange[analysisIdx][1]) ? kTRUE : kFALSE;
	
	isTrackAcceptedPid = (isTrackAcceptedPidTPC && isTrackAcceptedPidTOF);
      }
      else {
	// -- is track Spallation proton/anti-proton
	Bool_t isTrackSpallationProton = (pt > ptRangeSpallation[0] && pt < ptRangeSpallation[1] 
					  && TMath::Abs(nSigma[0]) < nSigmaProtonMaxSpallation)  ? kTRUE : kFALSE;
	
	isTrackAcceptedPid = isTrackSpallationProton;
      }
      
      // -->> is track accepted  - clusters/dca && kinematics && PID
      Bool_t isTrackAccepted = (isTrackAcceptedKin && isTrackAcceptedCut && isTrackAcceptedPid);

      // -- fill ThnSparse - tracks
      // ------------------------------------------------------------------
      Double_t aTrack[12] = {Double_t(centrality), pt, eta, sign, DCA,
			     Double_t(nHitsDedx), Double_t(nHitsFit), Double_t(nFitPoss), ratio, 
			     Double_t(isInRefMult), Double_t(isTrackAccepted), nSigma[analysisIdx]};
      
#if TRACK_THN      
      fHnTrackUnCorr->Fill(aTrack);
#endif

      FillTrackHists(aTrack, 0);
      FillRunByRunTrackHists(aTrack, 0, runIdx);

      // -- reject track
      // ------------------------------------------------------------------
      if (!isTrackAccepted)
	continue;
      
      FillTrackHists(aTrack, 1);
      FillRunByRunTrackHists(aTrack, 1, runIdx);

      // ------------------------------------------------------------------
      // -- Add up for event multiplicity
      // ------------------------------------------------------------------
      //  idxPart = 0 -> anti particle
      //  idxPart = 1 -> particle
      Int_t idxPart    = (sign < 0) ? 0 : 1;
      Int_t idxEtaSign = (eta  < 0) ? 0 : 1;

#if USE_RANDOM_EFF      
      // -- discard a random amount of tracks
      if (gRandom->Rndm() > randomEff[idxPart][centrality])
	continue;
#endif

      // -- Apply Charge separation 
      //    -> default: off  => 0  
      //    -> positive particles from postive eta / negative particles from negative eta => 1
      //    -> positive particles from negative eta / negative particles from postive eta => 2

      Int_t count = 0;
      if (useModeChargeSeparation == 0) 
	++count;
      else if ( (useModeChargeSeparation == 1) && (idxPart == idxEtaSign) ) 
	++count;
      else if ( (useModeChargeSeparation == 2) && (idxPart != idxEtaSign) ) 
	++count;

      // -- in full pt Range
      fNp[0][idxPart] += count;
      
      // -- divide in 2 parts
      if (pt < ptMidPoint[analysisIdx])
	fNp[1][idxPart] += count;
      else
	fNp[2][idxPart] += count;
    } // for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
    
    // -- Fill histograms
    // ------------------------------------------------------------------
    FillHistSetCent("Dist",       0, centrality);
    FillHistSetCent("Dist_lower", 1, centrality);
    FillHistSetCent("Dist_upper", 2, centrality);

    // ------------------------------------------------------------------

    Double_t aEventRun[] = { Double_t(nRefMultTracks),  Double_t(nRefMultXTracks[analysisIdx]),
			     Double_t(fNp[0][1] - fNp[0][0]), 
			     Double_t(fNp[0][0]), Double_t (fNp[0][1])};

    FillRunByRunEventHists(aEventRun, 1, runIdx); 
  
    // ------------------------------------------------------------------

    for(int i=0;i<=9;i++) 
      for(int j=0;j<=9;j++) 
	for(int k=0;k<=9;k++) 
	  for(int h=0;h<=9;h++) 
	    if((i+j+k+h)<=8) {
	      fact[i][j][k][h][centrality+1]->Fill(nRefMultXTracksCorr[analysisIdx], NN(fNp[1][1],i) * NN(fNp[2][1],j) * NN(fNp[1][0],k) * NN(fNp[2][0],h));
	      fact[i][j][k][h][0]->Fill(           nRefMultXTracksCorr[analysisIdx], NN(fNp[1][1],i) * NN(fNp[2][1],j) * NN(fNp[1][0],k) * NN(fNp[2][0],h));
	    }
    
  } // for (Long64_t idxEvent = 0; idxEvent < nEvents; idxEvent++) {

  cout << "Event Loop: " << nEvents << " events have been done!" << endl;

  // ------------------------------------------------------------------
  // -- Write Outlist
  // ------------------------------------------------------------------

  TH1::AddDirectory(oldStatus);

  TString collision(Form("Moments_hist_AuAu%sGeV_%s", exactEnergies[energyIdx], name2[analysisIdx]));
  TFile *outFile = new TFile(Form("%s.root", collision.Data()), "recreate");   
  outFile->cd();
  fOutList->Write(fOutList->GetName(), TObject::kSingleKey);
  outFile->Close();
  
  return 0;
}

//________________________________________________________________________
Int_t RefMultCorrection(Double_t vz, Int_t refmultIn, Int_t anaIdx) {
  // -- Correction of refmultX

  Double_t refMultZ = aRefMultCorrPar[anaIdx][0] + 
    aRefMultCorrPar[anaIdx][1]*vz + aRefMultCorrPar[anaIdx][2]*vz*vz + aRefMultCorrPar[anaIdx][3]*vz*vz*vz + 
    aRefMultCorrPar[anaIdx][4]*vz*vz*vz*vz + aRefMultCorrPar[anaIdx][5]*vz*vz*vz*vz*vz + aRefMultCorrPar[anaIdx][6]*vz*vz*vz*vz*vz*vz;

  Double_t Hovno    = (aRefMultCorrPar[anaIdx][0] + aRefMultCorrPar[anaIdx][7]) / refMultZ;

  Double_t refMultD = Double_t(refmultIn) + gRandom->Rndm(); // random sampling over bin width -> avoid peak structures in corrected distribution
  
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
Int_t GetCentrality(Int_t nRefMultTracksCorr, Int_t anaIdx) {
  // -- Fill centrality statistics of accepted events
  //    first bin = 0
  
  Int_t centrality = -1;
  
  if (nRefMultTracksCorr >= aRefMult[anaIdx][0])
    centrality = 0;   // 0-5

  else if (nRefMultTracksCorr < aRefMult[anaIdx][fNCentralityBins-1])
    centrality = fNCentralityBins;   // 80-100    

  else {
    for (Int_t idx = 1; idx < fNCentralityBins; ++idx) {
      if (nRefMultTracksCorr >= aRefMult[anaIdx][idx] && nRefMultTracksCorr < aRefMult[anaIdx][idx-1])
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
    list->SetName(Form("f%s_multHists%s", name[analysisIdx], multNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH1F(Form("hCentralityStat%s", multNames[ii]),               Form("Centrality statistics%s;Centrality Bins;Events",multTitles[ii]),      
		       fNCentralityBins,-0.5,fNCentralityBins-0.5));

    TH1F* hCentralityStat = static_cast<TH1F*>(list->Last());    
    for (Int_t jj = 0; jj < fNCentralityBins; jj++)
      hCentralityStat->GetXaxis()->SetBinLabel(jj+1, aCentralityNames[jj]);
  
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    list->Add(new TH1F(Form("hRefMultStat%s", multNames[ii]),                    Form("RefMult  Statistics%s;RefMult;Events",multTitles[ii]),
		       601, 0., 600.));
    list->Add(new TH1F(Form("hRefMult2Stat%s", multNames[ii]),                   Form("RefMult2 Statistics%s;RefMult2;Events",multTitles[ii]),  
		       601, 0., 600.));

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH2F(Form("hRefMult2_nGlobalTracks%s", multNames[ii]),         Form("RefMult2 vs nGlobalTracks%s;RefMult2;nGlobalTracks",multTitles[ii]),
		       601, 0., 600., 2501, 0., 2500.));
    list->Add(new TH2F(Form("hRefMult2_nPrimaryTracks%s", multNames[ii]),        Form("RefMult2 vs nPrimaryTracks%s;RefMult2;nPrimaryTracks",multTitles[ii]),
		       601, 0., 600., 1001, 0., 1000.));
    list->Add(new TH2F(Form("hRefMult2_nTOFMatch%s", multNames[ii]),             Form("RefMult2 vs nTOFMatch%s;RefMult2;nTOFMatch",multTitles[ii]),
		       601, 0., 600., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hRefMult_nGlobalTracks%s", multNames[ii]),          Form("RefMult vs nGlobalTracks%s;RefMult;nGlobalTracks",multTitles[ii]),
		       601, 0., 600., 2501, 0., 2500.));
    list->Add(new TH2F(Form("hRefMult_nPrimaryTracks%s", multNames[ii]),         Form("RefMult vs nPrimaryTracks%s;RefMult;nPrimaryTracks",multTitles[ii]),
		       601, 0., 600., 1001, 0., 1000.));
    list->Add(new TH2F(Form("hRefMult_nTOFMatch%s", multNames[ii]),              Form("RefMult vs nTOFMatch%s;RefMult;nTOFMatch",multTitles[ii]),
		       601, 0., 600., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hRefMult2_RefMult2Corr%s", multNames[ii]),          Form("RefMult2 vs RefMult2Corr%s;RefMult2;RefMult2Corr",multTitles[ii]),
		       601, 0., 600., 601, 0., 600.));
    list->Add(new TH2F(Form("hRefMult2_RefMult%s", multNames[ii]),               Form("RefMult2 vs RefMult%s;RefMult2;RefMult",multTitles[ii]),
		       601, 0., 600., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hnPrimaryTracks_nGlobalTracks%s", multNames[ii]),   Form("nPrimaryTracks vs nGlobalTracks%s;nPrimaryTracks;nGlobalTracks",multTitles[ii]),
		       1001, 0., 1000., 2501, 0., 2500.));
    list->Add(new TH2F(Form("hnPrimaryTracks_nTOFMatch%s", multNames[ii]),       Form("nPrimaryTracks vs nTOFMatch%s; nPrimaryTracks;nTOFMatch",multTitles[ii]),
		       1001, 0., 1000., 601, 0., 600.));    
    list->Add(new TH2F(Form("hnGlobalTracks_nTOFMatch%s", multNames[ii]),        Form("nGlobalTracks vs nTOFMatch%s;nGlobalTracks;nTOFMatch",multTitles[ii]),
		       2501, 0., 2500., 601, 0., 600.));
      
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

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_multHists%s", name[analysisIdx], multNames[mode])));
  
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
  // -- Initialize event QA histograms

  for (Int_t ii = 0 ; ii < 2 ; ++ii) {

    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_eventHists_%s", name[analysisIdx], qaNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH1F(Form("vx_%s", qaNames[ii]),        Form("#it{v}_{x}%s;#it{v}_{x} (cm);Events", qaTitles[ii]), 
		       binHnEvent[1], minHnEvent[1], maxHnEvent[1]));

    list->Add(new TH1F(Form("vy_%s", qaNames[ii]),        Form("#it{v}_{y}%s;#it{v}_{y} (cm);Events", qaTitles[ii]),
		       binHnEvent[2], minHnEvent[2], maxHnEvent[2]));

    list->Add(new TH1F(Form("vz_%s", qaNames[ii]),        Form("#it{v}_{z}%s;#it{v}_{z} (cm);Events", qaTitles[ii]), 
		       binHnEvent[3], minHnEvent[3], maxHnEvent[3]));

    list->Add(new TH1F(Form("shiftedVr_%s", qaNames[ii]), Form("#it{v}_{r}^{shifted}%s;#it{v}_{r}^{shifted} (cm);Events", qaTitles[ii]), 
		       binHnEvent[4], minHnEvent[4], maxHnEvent[4]));

    list->Add(new TH2F(Form("vxvy_%s", qaNames[ii]),      Form("#it{v}_{x} vs #it{v}_{y}%s;#it{v}_{x} (cm); #it{v}_{y} (cm)", qaTitles[ii]), 
		       binHnEvent[1], minHnEvent[1], maxHnEvent[1], binHnEvent[2], minHnEvent[2], maxHnEvent[2]));
    
    list->Add(new TH1F(Form("vzVpd_%s", qaNames[ii]),     Form("#it{v}_{z}^{vpd}%s;#it{v}_{z}^{vpd} (cm);Events", qaTitles[ii]),
		       binHnEvent[9], minHnEvent[9], maxHnEvent[9]));
    
    list->Add(new TH1F(Form("deltaVz_%s", qaNames[ii]),   Form("#Delta#it{v}_{z}%s;#Delta#it{v}_{z} (cm);Events", qaTitles[ii]),
		       binHnEvent[10], minHnEvent[10], maxHnEvent[10]));
  }
}

//________________________________________________________________________
void InitializeRunByRunEventHists() {
  // -- Initialize run-by-run event distributions
  
  for (Int_t ii = 0 ; ii < 2 ; ++ii) {
    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_runByRunEventHists_%s", name[analysisIdx], qaNames[ii]));
    list->SetOwner(kTRUE);
    
    list->Add(new TProfile(Form("pRefMult_%s", qaNames[ii]), Form("<RefMult>%s;<RefMult>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    //    list->Add(new TProfile(Form("pRefMult2_%s", qaNames[ii]), Form("<RefMult2>%s;<RefMult2>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    list->Add(new TProfile(Form("pRefMult3_%s", qaNames[ii]), Form("<RefMult3>%s;<RefMult3>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    //    list->Add(new TProfile(Form("pRefMult4_%s", qaNames[ii]), Form("<RefMult4>%s;<RefMult4>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));

    list->Add(new TProfile(Form("pNetProton_%s", qaNames[ii]), Form("<NetProton>%s;<NetProton>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    list->Add(new TProfile(Form("pneg_%s", qaNames[ii]), Form("<neg>%s;<neg>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    list->Add(new TProfile(Form("ppos_%s", qaNames[ii]), Form("<pos>%s;<pos>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
  }
}

//________________________________________________________________________
void InitializeRunByRunTrackHists() {
  // -- Initialize run-by-run event distributions
  
  for (Int_t ii = 0 ; ii < 2 ; ++ii) {
    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_runByRunTrackHists_%s", name[analysisIdx], qaNames[ii]));
    list->SetOwner(kTRUE);
    
    list->Add(new TProfile(Form("pPt_neg_%s", qaNames[ii]), Form("<pT_neg>%s;<pT_neg>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    list->Add(new TProfile(Form("pPt_pos_%s", qaNames[ii]), Form("<pT_pos>%s;<pT_pos>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));

    list->Add(new TProfile(Form("pDca_%s", qaNames[ii]), Form("<dca>%s;<dca>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));
    list->Add(new TProfile(Form("pNSigmaProton_%s", qaNames[ii]), Form("<nSigmaProton>%s;<nSigmaProton>;Run", qaTitles[ii]), nRuns, -0.5, nRuns-0.5));

    list->Add(new TH2D(Form("hPt_neg_%s", qaNames[ii]), Form("pT_neg%s;Run;pT", qaTitles[ii]), 
		       nRuns, -0.5, nRuns-0.5,
		       binHnUnCorr[1], minHnUnCorr[1], maxHnUnCorr[1]));
    list->Add(new TH2D(Form("hPt_pos_%s", qaNames[ii]), Form("pT_pos%s;Run;pT", qaTitles[ii]), 
		       nRuns, -0.5, nRuns-0.5,
		       binHnUnCorr[1], minHnUnCorr[1], maxHnUnCorr[1]));

  }
}



//________________________________________________________________________
void FillEventHists(Double_t *aEvent, Int_t mode) {
  // -- Fill event QA histograms      

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_eventHists_%s", name[analysisIdx], qaNames[mode])));
  (static_cast<TH1F*>(list->FindObject(Form("vx_%s",qaNames[mode]))))->Fill(aEvent[1]);
  (static_cast<TH1F*>(list->FindObject(Form("vy_%s",qaNames[mode]))))->Fill(aEvent[2]);
  (static_cast<TH1F*>(list->FindObject(Form("vz_%s",qaNames[mode]))))->Fill(aEvent[3]);
  (static_cast<TH1F*>(list->FindObject(Form("shiftedVr_%s",qaNames[mode]))))->Fill(aEvent[3]);
  (static_cast<TH2F*>(list->FindObject(Form("vxvy_%s",qaNames[mode]))))->Fill(aEvent[1], aEvent[2]);
  (static_cast<TH1F*>(list->FindObject(Form("vzVpd_%s",qaNames[mode]))))->Fill(aEvent[9]);
  (static_cast<TH1F*>(list->FindObject(Form("deltaVz_%s",qaNames[mode]))))->Fill(aEvent[10]);
}

//________________________________________________________________________
void InitializeTrackHists() {
  // -- Initialize track QA histograms

  for (Int_t ii = 0 ; ii < 2 ; ++ii) {
    fOutList->Add(new TList);
    TList *list = static_cast<TList*>(fOutList->Last());
    list->SetName(Form("f%s_trackHists_%s", name[analysisIdx], qaNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    for (Int_t idxSign = 0 ; idxSign < 2; ++idxSign) { 
      list->Add(new TH1F(Form("pt_%s_%d", qaNames[ii], idxSign),                Form("#it{p}_{T}%s;#it{p}_{T} (GeV/#it{c});Tracks", qaTitles[ii]),  
			 binHnUnCorr[1], minHnUnCorr[1], maxHnUnCorr[1]));
      list->Add(new TH1F(Form("eta_%s_%d", qaNames[ii], idxSign),               Form("#eta%s;#eta;Tracks", qaTitles[ii]),  
			 binHnUnCorr[2], minHnUnCorr[2], maxHnUnCorr[2]));
      list->Add(new TH1F(Form("dca_%s_%d", qaNames[ii], idxSign),               Form("DCA%s;DCA (cm);Tracks", qaTitles[ii]),  
			 binHnUnCorr[4], minHnUnCorr[4], maxHnUnCorr[4]));
      list->Add(new TH1F(Form("nHitsDedx_%s_%d", qaNames[ii], idxSign),         Form("nHitsDedx%s;nHitsDedx;Tracks", qaTitles[ii]), 
			 binHnUnCorr[5], minHnUnCorr[5], maxHnUnCorr[5]));
      list->Add(new TH1F(Form("nHitsFit_%s_%d", qaNames[ii], idxSign),          Form("nHitsFit%s;nHitsFit;Tracks", qaTitles[ii]),  
			 binHnUnCorr[6], minHnUnCorr[6], maxHnUnCorr[6]));
      list->Add(new TH1F(Form("nHitsFit_nFitPoss_%s_%d", qaNames[ii], idxSign), Form("nHitsFit/nFitPoss%s;nHitsFit/nFitPoss;Tracks", qaTitles[ii]),  
			 binHnUnCorr[8], minHnUnCorr[8], maxHnUnCorr[8]));
      list->Add(new TH1F(Form("nSigmaP_%s_%d", qaNames[ii], idxSign),           Form("nSigmaProton%s;n#Sigma_{proton};Tracks", qaTitles[ii]),  
			 binHnUnCorr[11], minHnUnCorr[11], maxHnUnCorr[11]));
      list->Add(new TH2F(Form("nSigmaP_pt_%s_%d", qaNames[ii], idxSign),        Form("nSigmaProton vs #it{p}_{T}%s;n#sigma_{proton};#it{p}_{T} (GeV/#it{c})", qaTitles[ii]),
			 binHnUnCorr[11], minHnUnCorr[11], maxHnUnCorr[11], binHnUnCorr[1], minHnUnCorr[1], maxHnUnCorr[1]));
    }
  }
}

//________________________________________________________________________
void FillTrackHists(Double_t *aTrack, Int_t mode) {
  // -- Fill track QA histograms

  Int_t idxSign = (aTrack[3] < 0) ? 0 : 1;

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_trackHists_%s", name[analysisIdx], qaNames[mode])));
  (static_cast<TH1F*>(list->FindObject(Form("pt_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[1]);
  (static_cast<TH1F*>(list->FindObject(Form("eta_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[2]);
  (static_cast<TH1F*>(list->FindObject(Form("dca_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[4]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsDedx_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[5]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsFit_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[6]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsFit_nFitPoss_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[8]);
  (static_cast<TH1F*>(list->FindObject(Form("nSigmaP_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[11]);
  (static_cast<TH2F*>(list->FindObject(Form("nSigmaP_pt_%s_%d", qaNames[mode], idxSign))))->Fill(aTrack[11], aTrack[1]);
}

//________________________________________________________________________
void FillRunByRunEventHists(Double_t *aEvent, Int_t mode, Int_t runIdx) {
  // -- Fill run-by-run event distributions
  
  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_runByRunEventHists_%s", name[analysisIdx], qaNames[mode])));
  (static_cast<TProfile*>(list->FindObject(Form("pRefMult_%s",   qaNames[mode]))))->Fill(runIdx, aEvent[0]);
  (static_cast<TProfile*>(list->FindObject(Form("pRefMult3_%s",  qaNames[mode]))))->Fill(runIdx, aEvent[1]);
  (static_cast<TProfile*>(list->FindObject(Form("pNetProton_%s", qaNames[mode]))))->Fill(runIdx, aEvent[2]);
  (static_cast<TProfile*>(list->FindObject(Form("pneg_%s",       qaNames[mode]))))->Fill(runIdx, aEvent[3]);
  (static_cast<TProfile*>(list->FindObject(Form("ppos_%s",       qaNames[mode]))))->Fill(runIdx, aEvent[4]);
}

//________________________________________________________________________
void FillRunByRunTrackHists(Double_t *aTrack, Int_t mode, Int_t runIdx) {
  // -- Fill run-by-run track distributions

  Int_t idxSign = (aTrack[3] < 0) ? 0 : 1;  

  TList* list = static_cast<TList*>(fOutList->FindObject(Form("f%s_runByRunTrackHists_%s", name[analysisIdx], qaNames[mode])));
  
  if (idxSign == 0) {
    (static_cast<TProfile*>(list->FindObject(Form("pPt_neg_%s", qaNames[mode]))))->Fill(runIdx, aTrack[1]);
    (static_cast<TH2D*>(list->FindObject(Form("hPt_neg_%s", qaNames[mode]))))->Fill(runIdx, aTrack[1]);
  }
  else {
    (static_cast<TProfile*>(list->FindObject(Form("pPt_pos_%s", qaNames[mode]))))->Fill(runIdx, aTrack[1]);
    (static_cast<TH2D*>(list->FindObject(Form("hPt_pos_%s", qaNames[mode]))))->Fill(runIdx, aTrack[1]);
  }

  (static_cast<TProfile*>(list->FindObject(Form("pDca_%s", qaNames[mode]))))->Fill(runIdx, aTrack[4]);  
}

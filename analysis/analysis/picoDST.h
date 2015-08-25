//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 24 17:18:05 2015 by ROOT version 5.34/09
// from TChain PicoDst/
//////////////////////////////////////////////////////////

#ifndef picoDST_h
#define picoDST_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>


// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
//#include ".sl64_gcc447/obj/StRoot/StarClassLibrary/StThreeVector.hh"
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxEvent = 1;
const Int_t kMaxTracks = 20000;
const Int_t kMaxTrigger = 1;
const Int_t kMaxBTOWHit = 1;
const Int_t kMaxBTofHit = 20000;
const Int_t kMaxV0Ks = 1;
const Int_t kMaxV0L = 1;
const Int_t kMaxV0Lbar = 1;

class picoDST {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event_;
   Int_t           Event_mRunId[kMaxEvent];   //[Event_]
   Int_t           Event_mEventId[kMaxEvent];   //[Event_]
   UShort_t        Event_mFillId[kMaxEvent];   //[Event_]
   Float_t         Event_mBField[kMaxEvent];   //[Event_]
   Float_t         Event_mPrimaryVertex_mX1[kMaxEvent];   //[Event_]
   Float_t         Event_mPrimaryVertex_mX2[kMaxEvent];   //[Event_]
   Float_t         Event_mPrimaryVertex_mX3[kMaxEvent];   //[Event_]
   Float_t         Event_mSecondVertex_mX1[kMaxEvent];   //[Event_]
   Float_t         Event_mSecondVertex_mX2[kMaxEvent];   //[Event_]
   Float_t         Event_mSecondVertex_mX3[kMaxEvent];   //[Event_]
   UShort_t        Event_mTriggerWord[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultFtpcEast[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultFtpcWest[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultNeg[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultPos[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMult2NegEast[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMult2PosEast[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMult2NegWest[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMult2PosWest[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultHalfNegEast[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultHalfPosEast[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultHalfNegWest[kMaxEvent];   //[Event_]
   UShort_t        Event_mRefMultHalfPosWest[kMaxEvent];   //[Event_]
   UChar_t         Event_mNVpdHitsEast[kMaxEvent];   //[Event_]
   UChar_t         Event_mNVpdHitsWest[kMaxEvent];   //[Event_]
   UShort_t        Event_mNT0[kMaxEvent];   //[Event_]
   Short_t         Event_mVzVpd[kMaxEvent];   //[Event_]
   UShort_t        Event_mZDCx[kMaxEvent];   //[Event_]
   UShort_t        Event_mBBCx[kMaxEvent];   //[Event_]
   Float_t         Event_mBackgroundRate[kMaxEvent];   //[Event_]
   Float_t         Event_mBbcBlueBackgroundRate[kMaxEvent];   //[Event_]
   Float_t         Event_mBbcYellowBackgroundRate[kMaxEvent];   //[Event_]
   Float_t         Event_mBbcEastRate[kMaxEvent];   //[Event_]
   Float_t         Event_mBbcWestRate[kMaxEvent];   //[Event_]
   Float_t         Event_mZdcEastRate[kMaxEvent];   //[Event_]
   Float_t         Event_mZdcWestRate[kMaxEvent];   //[Event_]
   UShort_t        Event_mVpd[kMaxEvent][64];   //[Event_]
   UShort_t        Event_mZdcSumAdcEast[kMaxEvent];   //[Event_]
   UShort_t        Event_mZdcSumAdcWest[kMaxEvent];   //[Event_]
   UShort_t        Event_mZdcSmdEastHorizontal[kMaxEvent][8];   //[Event_]
   UShort_t        Event_mZdcSmdEastVertical[kMaxEvent][8];   //[Event_]
   UShort_t        Event_mZdcSmdWestHorizontal[kMaxEvent][8];   //[Event_]
   UShort_t        Event_mZdcSmdWestVertical[kMaxEvent][8];   //[Event_]
   Float_t         Event_mSpaceCharge[kMaxEvent];   //[Event_]
   UShort_t        Event_mbTofTrayMultiplicity[kMaxEvent];   //[Event_]
   UShort_t        Event_mNumberOfGlobalTracks[kMaxEvent];   //[Event_]
   Float_t         Event_mRanking[kMaxEvent];   //[Event_]
   UShort_t        Event_mNBEMCMatch[kMaxEvent];   //[Event_]
   UShort_t        Event_mNBTOFMatch[kMaxEvent];   //[Event_]
   UShort_t        Event_mBbcAdcEast[kMaxEvent][24];   //[Event_]
   UShort_t        Event_mBbcAdcWest[kMaxEvent][24];   //[Event_]
   UChar_t         Event_mHT_Th[kMaxEvent][3];   //[Event_]
   UChar_t         Event_mJP_Th[kMaxEvent][3];   //[Event_]
   Int_t           Tracks_;
   UShort_t        Tracks_mId[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mChi2[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mChi2Prob[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mGMomentum_mX1[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mGMomentum_mX2[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mGMomentum_mX3[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mPMomentum_mX1[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mPMomentum_mX2[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mPMomentum_mX3[kMaxTracks];   //[Tracks_]
   UChar_t         Tracks_mFlowFlag[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mOriginX[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mOriginY[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mOriginZ[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mGDca[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mDedx[kMaxTracks];   //[Tracks_]
   Char_t          Tracks_mNHitsFit[kMaxTracks];   //[Tracks_]
   UChar_t         Tracks_mNHitsMax[kMaxTracks];   //[Tracks_]
   UChar_t         Tracks_mNHitsDedx[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mNSigmaPion[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mNSigmaKaon[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mNSigmaProton[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mNSigmaElectron[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofCellId[kMaxTracks];   //[Tracks_]
   UChar_t         Tracks_mBTofMatchFlag[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mBTof[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mBTofBeta[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofYLocal[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofZLocal[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofHitPosX[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofHitPosY[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofHitPosZ[kMaxTracks];   //[Tracks_]
   Int_t           Trigger_;
   UInt_t          Trigger_fUniqueID[kMaxTrigger];   //[Trigger_]
   UInt_t          Trigger_fBits[kMaxTrigger];   //[Trigger_]
   UChar_t         Trigger_mFlag[kMaxTrigger];   //[Trigger_]
   UShort_t        Trigger_mId[kMaxTrigger];   //[Trigger_]
   UShort_t        Trigger_mAdc[kMaxTrigger];   //[Trigger_]
   Int_t           BTOWHit_;
   UInt_t          BTOWHit_fUniqueID[kMaxBTOWHit];   //[BTOWHit_]
   UInt_t          BTOWHit_fBits[kMaxBTOWHit];   //[BTOWHit_]
   UShort_t        BTOWHit_mId[kMaxBTOWHit];   //[BTOWHit_]
   UShort_t        BTOWHit_mAdc[kMaxBTOWHit];   //[BTOWHit_]
   Short_t         BTOWHit_mE[kMaxBTOWHit];   //[BTOWHit_]
   Int_t           BTofHit_;
   UInt_t          BTofHit_fUniqueID[kMaxBTofHit];   //[BTofHit_]
   UInt_t          BTofHit_fBits[kMaxBTofHit];   //[BTofHit_]
   Short_t         BTofHit_mId[kMaxBTofHit];   //[BTofHit_]
   Int_t           V0Ks_;
   Short_t         V0Ks_mIndex2Track[kMaxV0Ks][2];   //[V0Ks_]
   //StThreeVector<float> V0Ks_mMomentum[2][kMaxV0Ks];
   Float_t         V0Ks_mV0Pos_mX1[kMaxV0Ks];   //[V0Ks_]
   Float_t         V0Ks_mV0Pos_mX2[kMaxV0Ks];   //[V0Ks_]
   Float_t         V0Ks_mV0Pos_mX3[kMaxV0Ks];   //[V0Ks_]
   UShort_t        V0Ks_mDcaDaughters[kMaxV0Ks];   //[V0Ks_]
   UShort_t        V0Ks_mDca2Vtx[kMaxV0Ks];   //[V0Ks_]
   Float_t         V0Ks_mM[kMaxV0Ks];   //[V0Ks_]
   Int_t           V0L_;
   Short_t         V0L_mIndex2Track[kMaxV0L][2];   //[V0L_]
   //StThreeVector<float> V0L_mMomentum[2][kMaxV0L];
   Float_t         V0L_mV0Pos_mX1[kMaxV0L];   //[V0L_]
   Float_t         V0L_mV0Pos_mX2[kMaxV0L];   //[V0L_]
   Float_t         V0L_mV0Pos_mX3[kMaxV0L];   //[V0L_]
   UShort_t        V0L_mDcaDaughters[kMaxV0L];   //[V0L_]
   UShort_t        V0L_mDca2Vtx[kMaxV0L];   //[V0L_]
   Float_t         V0L_mM[kMaxV0L];   //[V0L_]
   Int_t           V0Lbar_;
   Short_t         V0Lbar_mIndex2Track[kMaxV0Lbar][2];   //[V0Lbar_]
 //  StThreeVector<float> V0Lbar_mMomentum[2][kMaxV0Lbar];
   Float_t         V0Lbar_mV0Pos_mX1[kMaxV0Lbar];   //[V0Lbar_]
   Float_t         V0Lbar_mV0Pos_mX2[kMaxV0Lbar];   //[V0Lbar_]
   Float_t         V0Lbar_mV0Pos_mX3[kMaxV0Lbar];   //[V0Lbar_]
   UShort_t        V0Lbar_mDcaDaughters[kMaxV0Lbar];   //[V0Lbar_]
   UShort_t        V0Lbar_mDca2Vtx[kMaxV0Lbar];   //[V0Lbar_]
   Float_t         V0Lbar_mM[kMaxV0Lbar];   //[V0Lbar_]

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_mRunId;   //!
   TBranch        *b_Event_mEventId;   //!
   TBranch        *b_Event_mFillId;   //!
   TBranch        *b_Event_mBField;   //!
   TBranch        *b_Event_mPrimaryVertex_mX1;   //!
   TBranch        *b_Event_mPrimaryVertex_mX2;   //!
   TBranch        *b_Event_mPrimaryVertex_mX3;   //!
   TBranch        *b_Event_mSecondVertex_mX1;   //!
   TBranch        *b_Event_mSecondVertex_mX2;   //!
   TBranch        *b_Event_mSecondVertex_mX3;   //!
   TBranch        *b_Event_mTriggerWord;   //!
   TBranch        *b_Event_mRefMultFtpcEast;   //!
   TBranch        *b_Event_mRefMultFtpcWest;   //!
   TBranch        *b_Event_mRefMultNeg;   //!
   TBranch        *b_Event_mRefMultPos;   //!
   TBranch        *b_Event_mRefMult2NegEast;   //!
   TBranch        *b_Event_mRefMult2PosEast;   //!
   TBranch        *b_Event_mRefMult2NegWest;   //!
   TBranch        *b_Event_mRefMult2PosWest;   //!
   TBranch        *b_Event_mRefMultHalfNegEast;   //!
   TBranch        *b_Event_mRefMultHalfPosEast;   //!
   TBranch        *b_Event_mRefMultHalfNegWest;   //!
   TBranch        *b_Event_mRefMultHalfPosWest;   //!
   TBranch        *b_Event_mNVpdHitsEast;   //!
   TBranch        *b_Event_mNVpdHitsWest;   //!
   TBranch        *b_Event_mNT0;   //!
   TBranch        *b_Event_mVzVpd;   //!
   TBranch        *b_Event_mZDCx;   //!
   TBranch        *b_Event_mBBCx;   //!
   TBranch        *b_Event_mBackgroundRate;   //!
   TBranch        *b_Event_mBbcBlueBackgroundRate;   //!
   TBranch        *b_Event_mBbcYellowBackgroundRate;   //!
   TBranch        *b_Event_mBbcEastRate;   //!
   TBranch        *b_Event_mBbcWestRate;   //!
   TBranch        *b_Event_mZdcEastRate;   //!
   TBranch        *b_Event_mZdcWestRate;   //!
   TBranch        *b_Event_mVpd;   //!
   TBranch        *b_Event_mZdcSumAdcEast;   //!
   TBranch        *b_Event_mZdcSumAdcWest;   //!
   TBranch        *b_Event_mZdcSmdEastHorizontal;   //!
   TBranch        *b_Event_mZdcSmdEastVertical;   //!
   TBranch        *b_Event_mZdcSmdWestHorizontal;   //!
   TBranch        *b_Event_mZdcSmdWestVertical;   //!
   TBranch        *b_Event_mSpaceCharge;   //!
   TBranch        *b_Event_mbTofTrayMultiplicity;   //!
   TBranch        *b_Event_mNumberOfGlobalTracks;   //!
   TBranch        *b_Event_mRanking;   //!
   TBranch        *b_Event_mNBEMCMatch;   //!
   TBranch        *b_Event_mNBTOFMatch;   //!
   TBranch        *b_Event_mBbcAdcEast;   //!
   TBranch        *b_Event_mBbcAdcWest;   //!
   TBranch        *b_Event_mHT_Th;   //!
   TBranch        *b_Event_mJP_Th;   //!
   TBranch        *b_Tracks_;   //!
   TBranch        *b_Tracks_mId;   //!
   TBranch        *b_Tracks_mChi2;   //!
   TBranch        *b_Tracks_mChi2Prob;   //!
   TBranch        *b_Tracks_mGMomentum_mX1;   //!
   TBranch        *b_Tracks_mGMomentum_mX2;   //!
   TBranch        *b_Tracks_mGMomentum_mX3;   //!
   TBranch        *b_Tracks_mPMomentum_mX1;   //!
   TBranch        *b_Tracks_mPMomentum_mX2;   //!
   TBranch        *b_Tracks_mPMomentum_mX3;   //!
   TBranch        *b_Tracks_mFlowFlag;   //!
   TBranch        *b_Tracks_mOriginX;   //!
   TBranch        *b_Tracks_mOriginY;   //!
   TBranch        *b_Tracks_mOriginZ;   //!
   TBranch        *b_Tracks_mGDca;   //!
   TBranch        *b_Tracks_mDedx;   //!
   TBranch        *b_Tracks_mNHitsFit;   //!
   TBranch        *b_Tracks_mNHitsMax;   //!
   TBranch        *b_Tracks_mNHitsDedx;   //!
   TBranch        *b_Tracks_mNSigmaPion;   //!
   TBranch        *b_Tracks_mNSigmaKaon;   //!
   TBranch        *b_Tracks_mNSigmaProton;   //!
   TBranch        *b_Tracks_mNSigmaElectron;   //!
   TBranch        *b_Tracks_mBTofCellId;   //!
   TBranch        *b_Tracks_mBTofMatchFlag;   //!
   TBranch        *b_Tracks_mBTof;   //!
   TBranch        *b_Tracks_mBTofBeta;   //!
   TBranch        *b_Tracks_mBTofYLocal;   //!
   TBranch        *b_Tracks_mBTofZLocal;   //!
   TBranch        *b_Tracks_mBTofHitPosX;   //!
   TBranch        *b_Tracks_mBTofHitPosY;   //!
   TBranch        *b_Tracks_mBTofHitPosZ;   //!
   TBranch        *b_Trigger_;   //!
   TBranch        *b_Trigger_fUniqueID;   //!
   TBranch        *b_Trigger_fBits;   //!
   TBranch        *b_Trigger_mFlag;   //!
   TBranch        *b_Trigger_mId;   //!
   TBranch        *b_Trigger_mAdc;   //!
   TBranch        *b_BTOWHit_;   //!
   TBranch        *b_BTOWHit_fUniqueID;   //!
   TBranch        *b_BTOWHit_fBits;   //!
   TBranch        *b_BTOWHit_mId;   //!
   TBranch        *b_BTOWHit_mAdc;   //!
   TBranch        *b_BTOWHit_mE;   //!
   TBranch        *b_BTofHit_;   //!
   TBranch        *b_BTofHit_fUniqueID;   //!
   TBranch        *b_BTofHit_fBits;   //!
   TBranch        *b_BTofHit_mId;   //!
   TBranch        *b_V0Ks_;   //!
   TBranch        *b_V0Ks_mIndex2Track;   //!
   TBranch        *b_V0Ks_mMomentum;   //!
   TBranch        *b_V0Ks_mV0Pos_mX1;   //!
   TBranch        *b_V0Ks_mV0Pos_mX2;   //!
   TBranch        *b_V0Ks_mV0Pos_mX3;   //!
   TBranch        *b_V0Ks_mDcaDaughters;   //!
   TBranch        *b_V0Ks_mDca2Vtx;   //!
   TBranch        *b_V0Ks_mM;   //!
   TBranch        *b_V0L_;   //!
   TBranch        *b_V0L_mIndex2Track;   //!
   TBranch        *b_V0L_mMomentum;   //!
   TBranch        *b_V0L_mV0Pos_mX1;   //!
   TBranch        *b_V0L_mV0Pos_mX2;   //!
   TBranch        *b_V0L_mV0Pos_mX3;   //!
   TBranch        *b_V0L_mDcaDaughters;   //!
   TBranch        *b_V0L_mDca2Vtx;   //!
   TBranch        *b_V0L_mM;   //!
   TBranch        *b_V0Lbar_;   //!
   TBranch        *b_V0Lbar_mIndex2Track;   //!
   TBranch        *b_V0Lbar_mMomentum;   //!
   TBranch        *b_V0Lbar_mV0Pos_mX1;   //!
   TBranch        *b_V0Lbar_mV0Pos_mX2;   //!
   TBranch        *b_V0Lbar_mV0Pos_mX3;   //!
   TBranch        *b_V0Lbar_mDcaDaughters;   //!
   TBranch        *b_V0Lbar_mDca2Vtx;   //!
   TBranch        *b_V0Lbar_mM;   //!

   picoDST(TTree *tree=0);
   virtual ~picoDST();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t    cent_deter(int , double *, int);
};

#endif

#ifdef picoDST_cxx
picoDST::picoDST(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("PicoDst",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("PicoDst","");
      //chain->Add("/project/projectdirs/star/starprod/picodsts/Run14/AuAu/15GeV/all/P14ii/047/15047104/st_physics_adc_15047104_raw_2500004.picoDst.root/PicoDst");
      
      char filename[400];

   ifstream filelist("file.list");


        while(filelist>>filename)
     {

        chain->Add(filename);
        std::cout<<"add file: "<<filename<<std::endl;
 //       if(chain->GetEntries()>100000) break;
                     }
      

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

picoDST::~picoDST()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t picoDST::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t picoDST::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void picoDST::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.mRunId", Event_mRunId, &b_Event_mRunId);
   fChain->SetBranchAddress("Event.mEventId", Event_mEventId, &b_Event_mEventId);
   fChain->SetBranchAddress("Event.mFillId", Event_mFillId, &b_Event_mFillId);
   fChain->SetBranchAddress("Event.mBField", Event_mBField, &b_Event_mBField);
   fChain->SetBranchAddress("Event.mPrimaryVertex.mX1", Event_mPrimaryVertex_mX1, &b_Event_mPrimaryVertex_mX1);
   fChain->SetBranchAddress("Event.mPrimaryVertex.mX2", Event_mPrimaryVertex_mX2, &b_Event_mPrimaryVertex_mX2);
   fChain->SetBranchAddress("Event.mPrimaryVertex.mX3", Event_mPrimaryVertex_mX3, &b_Event_mPrimaryVertex_mX3);
   fChain->SetBranchAddress("Event.mSecondVertex.mX1", Event_mSecondVertex_mX1, &b_Event_mSecondVertex_mX1);
   fChain->SetBranchAddress("Event.mSecondVertex.mX2", Event_mSecondVertex_mX2, &b_Event_mSecondVertex_mX2);
   fChain->SetBranchAddress("Event.mSecondVertex.mX3", Event_mSecondVertex_mX3, &b_Event_mSecondVertex_mX3);
   fChain->SetBranchAddress("Event.mTriggerWord", Event_mTriggerWord, &b_Event_mTriggerWord);
   fChain->SetBranchAddress("Event.mRefMultFtpcEast", Event_mRefMultFtpcEast, &b_Event_mRefMultFtpcEast);
   fChain->SetBranchAddress("Event.mRefMultFtpcWest", Event_mRefMultFtpcWest, &b_Event_mRefMultFtpcWest);
   fChain->SetBranchAddress("Event.mRefMultNeg", Event_mRefMultNeg, &b_Event_mRefMultNeg);
   fChain->SetBranchAddress("Event.mRefMultPos", Event_mRefMultPos, &b_Event_mRefMultPos);
   fChain->SetBranchAddress("Event.mRefMult2NegEast", Event_mRefMult2NegEast, &b_Event_mRefMult2NegEast);
   fChain->SetBranchAddress("Event.mRefMult2PosEast", Event_mRefMult2PosEast, &b_Event_mRefMult2PosEast);
   fChain->SetBranchAddress("Event.mRefMult2NegWest", Event_mRefMult2NegWest, &b_Event_mRefMult2NegWest);
   fChain->SetBranchAddress("Event.mRefMult2PosWest", Event_mRefMult2PosWest, &b_Event_mRefMult2PosWest);
   fChain->SetBranchAddress("Event.mRefMultHalfNegEast", Event_mRefMultHalfNegEast, &b_Event_mRefMultHalfNegEast);
   fChain->SetBranchAddress("Event.mRefMultHalfPosEast", Event_mRefMultHalfPosEast, &b_Event_mRefMultHalfPosEast);
   fChain->SetBranchAddress("Event.mRefMultHalfNegWest", Event_mRefMultHalfNegWest, &b_Event_mRefMultHalfNegWest);
   fChain->SetBranchAddress("Event.mRefMultHalfPosWest", Event_mRefMultHalfPosWest, &b_Event_mRefMultHalfPosWest);
   fChain->SetBranchAddress("Event.mNVpdHitsEast", Event_mNVpdHitsEast, &b_Event_mNVpdHitsEast);
   fChain->SetBranchAddress("Event.mNVpdHitsWest", Event_mNVpdHitsWest, &b_Event_mNVpdHitsWest);
   fChain->SetBranchAddress("Event.mNT0", Event_mNT0, &b_Event_mNT0);
   fChain->SetBranchAddress("Event.mVzVpd", Event_mVzVpd, &b_Event_mVzVpd);
   fChain->SetBranchAddress("Event.mZDCx", Event_mZDCx, &b_Event_mZDCx);
   fChain->SetBranchAddress("Event.mBBCx", Event_mBBCx, &b_Event_mBBCx);
   fChain->SetBranchAddress("Event.mBackgroundRate", Event_mBackgroundRate, &b_Event_mBackgroundRate);
   fChain->SetBranchAddress("Event.mBbcBlueBackgroundRate", Event_mBbcBlueBackgroundRate, &b_Event_mBbcBlueBackgroundRate);
   fChain->SetBranchAddress("Event.mBbcYellowBackgroundRate", Event_mBbcYellowBackgroundRate, &b_Event_mBbcYellowBackgroundRate);
   fChain->SetBranchAddress("Event.mBbcEastRate", Event_mBbcEastRate, &b_Event_mBbcEastRate);
   fChain->SetBranchAddress("Event.mBbcWestRate", Event_mBbcWestRate, &b_Event_mBbcWestRate);
   fChain->SetBranchAddress("Event.mZdcEastRate", Event_mZdcEastRate, &b_Event_mZdcEastRate);
   fChain->SetBranchAddress("Event.mZdcWestRate", Event_mZdcWestRate, &b_Event_mZdcWestRate);
   fChain->SetBranchAddress("Event.mVpd[64]", Event_mVpd, &b_Event_mVpd);
   fChain->SetBranchAddress("Event.mZdcSumAdcEast", Event_mZdcSumAdcEast, &b_Event_mZdcSumAdcEast);
   fChain->SetBranchAddress("Event.mZdcSumAdcWest", Event_mZdcSumAdcWest, &b_Event_mZdcSumAdcWest);
   fChain->SetBranchAddress("Event.mZdcSmdEastHorizontal[8]", Event_mZdcSmdEastHorizontal, &b_Event_mZdcSmdEastHorizontal);
   fChain->SetBranchAddress("Event.mZdcSmdEastVertical[8]", Event_mZdcSmdEastVertical, &b_Event_mZdcSmdEastVertical);
   fChain->SetBranchAddress("Event.mZdcSmdWestHorizontal[8]", Event_mZdcSmdWestHorizontal, &b_Event_mZdcSmdWestHorizontal);
   fChain->SetBranchAddress("Event.mZdcSmdWestVertical[8]", Event_mZdcSmdWestVertical, &b_Event_mZdcSmdWestVertical);
   fChain->SetBranchAddress("Event.mSpaceCharge", Event_mSpaceCharge, &b_Event_mSpaceCharge);
   fChain->SetBranchAddress("Event.mbTofTrayMultiplicity", Event_mbTofTrayMultiplicity, &b_Event_mbTofTrayMultiplicity);
   fChain->SetBranchAddress("Event.mNumberOfGlobalTracks", Event_mNumberOfGlobalTracks, &b_Event_mNumberOfGlobalTracks);
   fChain->SetBranchAddress("Event.mRanking", Event_mRanking, &b_Event_mRanking);
   fChain->SetBranchAddress("Event.mNBEMCMatch", Event_mNBEMCMatch, &b_Event_mNBEMCMatch);
   fChain->SetBranchAddress("Event.mNBTOFMatch", Event_mNBTOFMatch, &b_Event_mNBTOFMatch);
   fChain->SetBranchAddress("Event.mBbcAdcEast[24]", Event_mBbcAdcEast, &b_Event_mBbcAdcEast);
   fChain->SetBranchAddress("Event.mBbcAdcWest[24]", Event_mBbcAdcWest, &b_Event_mBbcAdcWest);
   fChain->SetBranchAddress("Event.mHT_Th[3]", Event_mHT_Th, &b_Event_mHT_Th);
   fChain->SetBranchAddress("Event.mJP_Th[3]", Event_mJP_Th, &b_Event_mJP_Th);
   fChain->SetBranchAddress("Tracks", &Tracks_, &b_Tracks_);
   fChain->SetBranchAddress("Tracks.mId", Tracks_mId, &b_Tracks_mId);
   fChain->SetBranchAddress("Tracks.mChi2", Tracks_mChi2, &b_Tracks_mChi2);
   fChain->SetBranchAddress("Tracks.mChi2Prob", Tracks_mChi2Prob, &b_Tracks_mChi2Prob);
   fChain->SetBranchAddress("Tracks.mGMomentum.mX1", Tracks_mGMomentum_mX1, &b_Tracks_mGMomentum_mX1);
   fChain->SetBranchAddress("Tracks.mGMomentum.mX2", Tracks_mGMomentum_mX2, &b_Tracks_mGMomentum_mX2);
   fChain->SetBranchAddress("Tracks.mGMomentum.mX3", Tracks_mGMomentum_mX3, &b_Tracks_mGMomentum_mX3);
   fChain->SetBranchAddress("Tracks.mPMomentum.mX1", Tracks_mPMomentum_mX1, &b_Tracks_mPMomentum_mX1);
   fChain->SetBranchAddress("Tracks.mPMomentum.mX2", Tracks_mPMomentum_mX2, &b_Tracks_mPMomentum_mX2);
   fChain->SetBranchAddress("Tracks.mPMomentum.mX3", Tracks_mPMomentum_mX3, &b_Tracks_mPMomentum_mX3);
   fChain->SetBranchAddress("Tracks.mFlowFlag", Tracks_mFlowFlag, &b_Tracks_mFlowFlag);
   fChain->SetBranchAddress("Tracks.mOriginX", Tracks_mOriginX, &b_Tracks_mOriginX);
   fChain->SetBranchAddress("Tracks.mOriginY", Tracks_mOriginY, &b_Tracks_mOriginY);
   fChain->SetBranchAddress("Tracks.mOriginZ", Tracks_mOriginZ, &b_Tracks_mOriginZ);
   fChain->SetBranchAddress("Tracks.mGDca", Tracks_mGDca, &b_Tracks_mGDca);
   fChain->SetBranchAddress("Tracks.mDedx", Tracks_mDedx, &b_Tracks_mDedx);
   fChain->SetBranchAddress("Tracks.mNHitsFit", Tracks_mNHitsFit, &b_Tracks_mNHitsFit);
   fChain->SetBranchAddress("Tracks.mNHitsMax", Tracks_mNHitsMax, &b_Tracks_mNHitsMax);
   fChain->SetBranchAddress("Tracks.mNHitsDedx", Tracks_mNHitsDedx, &b_Tracks_mNHitsDedx);
   fChain->SetBranchAddress("Tracks.mNSigmaPion", Tracks_mNSigmaPion, &b_Tracks_mNSigmaPion);
   fChain->SetBranchAddress("Tracks.mNSigmaKaon", Tracks_mNSigmaKaon, &b_Tracks_mNSigmaKaon);
   fChain->SetBranchAddress("Tracks.mNSigmaProton", Tracks_mNSigmaProton, &b_Tracks_mNSigmaProton);
   fChain->SetBranchAddress("Tracks.mNSigmaElectron", Tracks_mNSigmaElectron, &b_Tracks_mNSigmaElectron);
   fChain->SetBranchAddress("Tracks.mBTofCellId", Tracks_mBTofCellId, &b_Tracks_mBTofCellId);
   fChain->SetBranchAddress("Tracks.mBTofMatchFlag", Tracks_mBTofMatchFlag, &b_Tracks_mBTofMatchFlag);
   fChain->SetBranchAddress("Tracks.mBTof", Tracks_mBTof, &b_Tracks_mBTof);
   fChain->SetBranchAddress("Tracks.mBTofBeta", Tracks_mBTofBeta, &b_Tracks_mBTofBeta);
   fChain->SetBranchAddress("Tracks.mBTofYLocal", Tracks_mBTofYLocal, &b_Tracks_mBTofYLocal);
   fChain->SetBranchAddress("Tracks.mBTofZLocal", Tracks_mBTofZLocal, &b_Tracks_mBTofZLocal);
   fChain->SetBranchAddress("Tracks.mBTofHitPosX", Tracks_mBTofHitPosX, &b_Tracks_mBTofHitPosX);
   fChain->SetBranchAddress("Tracks.mBTofHitPosY", Tracks_mBTofHitPosY, &b_Tracks_mBTofHitPosY);
   fChain->SetBranchAddress("Tracks.mBTofHitPosZ", Tracks_mBTofHitPosZ, &b_Tracks_mBTofHitPosZ);
   fChain->SetBranchAddress("Trigger", &Trigger_, &b_Trigger_);
   fChain->SetBranchAddress("Trigger.fUniqueID", &Trigger_fUniqueID, &b_Trigger_fUniqueID);
   fChain->SetBranchAddress("Trigger.fBits", &Trigger_fBits, &b_Trigger_fBits);
   fChain->SetBranchAddress("Trigger.mFlag", &Trigger_mFlag, &b_Trigger_mFlag);
   fChain->SetBranchAddress("Trigger.mId", &Trigger_mId, &b_Trigger_mId);
   fChain->SetBranchAddress("Trigger.mAdc", &Trigger_mAdc, &b_Trigger_mAdc);
   fChain->SetBranchAddress("BTOWHit", &BTOWHit_, &b_BTOWHit_);
   fChain->SetBranchAddress("BTOWHit.fUniqueID", &BTOWHit_fUniqueID, &b_BTOWHit_fUniqueID);
   fChain->SetBranchAddress("BTOWHit.fBits", &BTOWHit_fBits, &b_BTOWHit_fBits);
   fChain->SetBranchAddress("BTOWHit.mId", &BTOWHit_mId, &b_BTOWHit_mId);
   fChain->SetBranchAddress("BTOWHit.mAdc", &BTOWHit_mAdc, &b_BTOWHit_mAdc);
   fChain->SetBranchAddress("BTOWHit.mE", &BTOWHit_mE, &b_BTOWHit_mE);
   fChain->SetBranchAddress("BTofHit", &BTofHit_, &b_BTofHit_);
   fChain->SetBranchAddress("BTofHit.fUniqueID", BTofHit_fUniqueID, &b_BTofHit_fUniqueID);
   fChain->SetBranchAddress("BTofHit.fBits", BTofHit_fBits, &b_BTofHit_fBits);
   fChain->SetBranchAddress("BTofHit.mId", BTofHit_mId, &b_BTofHit_mId);
   fChain->SetBranchAddress("V0Ks", &V0Ks_, &b_V0Ks_);
   fChain->SetBranchAddress("V0Ks.mIndex2Track[2]", &V0Ks_mIndex2Track, &b_V0Ks_mIndex2Track);
   //fChain->SetBranchAddress("V0Ks.mMomentum[2]", &V0Ks_mMomentum, &b_V0Ks_mMomentum);
   fChain->SetBranchAddress("V0Ks.mV0Pos.mX1", &V0Ks_mV0Pos_mX1, &b_V0Ks_mV0Pos_mX1);
   fChain->SetBranchAddress("V0Ks.mV0Pos.mX2", &V0Ks_mV0Pos_mX2, &b_V0Ks_mV0Pos_mX2);
   fChain->SetBranchAddress("V0Ks.mV0Pos.mX3", &V0Ks_mV0Pos_mX3, &b_V0Ks_mV0Pos_mX3);
   fChain->SetBranchAddress("V0Ks.mDcaDaughters", &V0Ks_mDcaDaughters, &b_V0Ks_mDcaDaughters);
   fChain->SetBranchAddress("V0Ks.mDca2Vtx", &V0Ks_mDca2Vtx, &b_V0Ks_mDca2Vtx);
   fChain->SetBranchAddress("V0Ks.mM", &V0Ks_mM, &b_V0Ks_mM);
   fChain->SetBranchAddress("V0L", &V0L_, &b_V0L_);
  // fChain->SetBranchAddress("V0L.mIndex2Track[2]", &V0L_mIndex2Track, &b_V0L_mIndex2Track);
   //fChain->SetBranchAddress("V0L.mMomentum[2]", &V0L_mMomentum, &b_V0L_mMomentum);
   fChain->SetBranchAddress("V0L.mV0Pos.mX1", &V0L_mV0Pos_mX1, &b_V0L_mV0Pos_mX1);
   fChain->SetBranchAddress("V0L.mV0Pos.mX2", &V0L_mV0Pos_mX2, &b_V0L_mV0Pos_mX2);
   fChain->SetBranchAddress("V0L.mV0Pos.mX3", &V0L_mV0Pos_mX3, &b_V0L_mV0Pos_mX3);
   fChain->SetBranchAddress("V0L.mDcaDaughters", &V0L_mDcaDaughters, &b_V0L_mDcaDaughters);
   fChain->SetBranchAddress("V0L.mDca2Vtx", &V0L_mDca2Vtx, &b_V0L_mDca2Vtx);
   fChain->SetBranchAddress("V0L.mM", &V0L_mM, &b_V0L_mM);
   fChain->SetBranchAddress("V0Lbar", &V0Lbar_, &b_V0Lbar_);
  // fChain->SetBranchAddress("V0Lbar.mIndex2Track[2]", &V0Lbar_mIndex2Track, &b_V0Lbar_mIndex2Track);
   //fChain->SetBranchAddress("V0Lbar.mMomentum[2]", &V0Lbar_mMomentum, &b_V0Lbar_mMomentum);
   fChain->SetBranchAddress("V0Lbar.mV0Pos.mX1", &V0Lbar_mV0Pos_mX1, &b_V0Lbar_mV0Pos_mX1);
   fChain->SetBranchAddress("V0Lbar.mV0Pos.mX2", &V0Lbar_mV0Pos_mX2, &b_V0Lbar_mV0Pos_mX2);
   fChain->SetBranchAddress("V0Lbar.mV0Pos.mX3", &V0Lbar_mV0Pos_mX3, &b_V0Lbar_mV0Pos_mX3);
   fChain->SetBranchAddress("V0Lbar.mDcaDaughters", &V0Lbar_mDcaDaughters, &b_V0Lbar_mDcaDaughters);
   fChain->SetBranchAddress("V0Lbar.mDca2Vtx", &V0Lbar_mDca2Vtx, &b_V0Lbar_mDca2Vtx);
   fChain->SetBranchAddress("V0Lbar.mM", &V0Lbar_mM, &b_V0Lbar_mM);
   Notify();
}

Bool_t picoDST::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void picoDST::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t picoDST::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}



Int_t  picoDST::cent_deter(int refmul,double *cent,int NumberOfCentrality)
{

        if(refmul>=cent[0]){ return 1;}
        if(refmul<cent[NumberOfCentrality-1]){return (NumberOfCentrality+1);}
        for(int m=0;m<=NumberOfCentrality-2;m++)
        {

                if(cent[m+1]<=refmul&&refmul<cent[m]){return (m+2);}

        }

 return 0;

}



#endif // #ifdef picoDST_cxx

#include "StMiniMcTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Riostream.h"
#include "TNtuple.h"
#include "../StRefMultCorr/StRefMultCorr.h"
#include <iostream>

using namespace std;

ClassImp(StMiniMcTree)

StMiniMcTree::StMiniMcTree(TString listfile, TTree *tree) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TChain * chain = new TChain("StMiniMcTree","");
    ifstream in; in.open(listfile.Data());
    Int_t nlines = 0; TString str;
    while (1) {
      in >> str;
      if (!in.good()) break;
      printf("%s\n",str.Data());
      chain->Add(str.Data());
      nlines++;
    }
    printf("%d files to analyze\n",nlines);
    tree = chain;
  }
    refmultCorrUtil = new StRefMultCorr();

  Init(tree);
}

StMiniMcTree::~StMiniMcTree() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StMiniMcTree::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t StMiniMcTree::LoadTree(Long64_t entry) {
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

void StMiniMcTree::Init(TTree *tree) {

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mEventId", &mEventId, &b_StMiniMcEvent_mEventId);
   fChain->SetBranchAddress("mRunId", &mRunId, &b_StMiniMcEvent_mRunId);
   fChain->SetBranchAddress("mOriginMult", &mOriginMult, &b_StMiniMcEvent_mOriginMult);
   fChain->SetBranchAddress("mCentralMult", &mCentralMult, &b_StMiniMcEvent_mCentralMult);
   fChain->SetBranchAddress("mCentrality", &mCentrality, &b_StMiniMcEvent_mCentrality);
   fChain->SetBranchAddress("mNUncorrectedNegativePrimaries", &mNUncorrectedNegativePrimaries, &b_StMiniMcEvent_mNUncorrectedNegativePrimaries);
   fChain->SetBranchAddress("mNUncorrectedPrimaries", &mNUncorrectedPrimaries, &b_StMiniMcEvent_mNUncorrectedPrimaries);
   fChain->SetBranchAddress("mNFtpcWUncorrectedPrimaries", &mNFtpcWUncorrectedPrimaries, &b_StMiniMcEvent_mNFtpcWUncorrectedPrimaries);
   fChain->SetBranchAddress("mNFtpcEUncorrectedPrimaries", &mNFtpcEUncorrectedPrimaries, &b_StMiniMcEvent_mNFtpcEUncorrectedPrimaries);
   fChain->SetBranchAddress("mMcMult", &mMcMult, &b_StMiniMcEvent_mMcMult);
   fChain->SetBranchAddress("mNMcNch", &mNMcNch, &b_StMiniMcEvent_mNMcNch);
   fChain->SetBranchAddress("mNMcFtpcWNch", &mNMcFtpcWNch, &b_StMiniMcEvent_mNMcFtpcWNch);
   fChain->SetBranchAddress("mNMcFtpcENch", &mNMcFtpcENch, &b_StMiniMcEvent_mNMcFtpcENch);
   fChain->SetBranchAddress("mNMcHminus", &mNMcHminus, &b_StMiniMcEvent_mNMcHminus);
   fChain->SetBranchAddress("mNMcGlobal", &mNMcGlobal, &b_StMiniMcEvent_mNMcGlobal);
   fChain->SetBranchAddress("mNMcGoodGlobal20", &mNMcGoodGlobal20, &b_StMiniMcEvent_mNMcGoodGlobal20);
   fChain->SetBranchAddress("mNRcGlobal", &mNRcGlobal, &b_StMiniMcEvent_mNRcGlobal);
   fChain->SetBranchAddress("mNRcGoodGlobal20", &mNRcGoodGlobal20, &b_StMiniMcEvent_mNRcGoodGlobal20);
   fChain->SetBranchAddress("mVertexX", &mVertexX, &b_StMiniMcEvent_mVertexX);
   fChain->SetBranchAddress("mVertexY", &mVertexY, &b_StMiniMcEvent_mVertexY);
   fChain->SetBranchAddress("mVertexZ", &mVertexZ, &b_StMiniMcEvent_mVertexZ);
   fChain->SetBranchAddress("mVertexCovMatrix[6]", mVertexCovMatrix, &b_StMiniMcEvent_mVertexCovMatrix);
   fChain->SetBranchAddress("mMcVertexX", &mMcVertexX, &b_StMiniMcEvent_mMcVertexX);
   fChain->SetBranchAddress("mMcVertexY", &mMcVertexY, &b_StMiniMcEvent_mMcVertexY);
   fChain->SetBranchAddress("mMcVertexZ", &mMcVertexZ, &b_StMiniMcEvent_mMcVertexZ);
   fChain->SetBranchAddress("mMagField", &mMagField, &b_StMiniMcEvent_mMagField);
   fChain->SetBranchAddress("mCenterOfMassEnergy", &mCenterOfMassEnergy, &b_StMiniMcEvent_mCenterOfMassEnergy);
   fChain->SetBranchAddress("mBackgroundRate", &mBackgroundRate, &b_StMiniMcEvent_mBackgroundRate);
   fChain->SetBranchAddress("mBeamMassNumberEast", &mBeamMassNumberEast, &b_StMiniMcEvent_mBeamMassNumberEast);
   fChain->SetBranchAddress("mBeamMassNumberWest", &mBeamMassNumberWest, &b_StMiniMcEvent_mBeamMassNumberWest);
   fChain->SetBranchAddress("mCtb", &mCtb, &b_StMiniMcEvent_mCtb);
   fChain->SetBranchAddress("mZdcE", &mZdcE, &b_StMiniMcEvent_mZdcE);
   fChain->SetBranchAddress("mZdcW", &mZdcW, &b_StMiniMcEvent_mZdcW);
   fChain->SetBranchAddress("mNMcTrack", &mNMcTrack, &b_StMiniMcEvent_mNMcTrack);
   fChain->SetBranchAddress("mNMatchedPair", &mNMatchedPair, &b_StMiniMcEvent_mNMatchedPair);
   fChain->SetBranchAddress("mNMergedPair", &mNMergedPair, &b_StMiniMcEvent_mNMergedPair);
   fChain->SetBranchAddress("mNSplitPair", &mNSplitPair, &b_StMiniMcEvent_mNSplitPair);
   fChain->SetBranchAddress("mNGhostPair", &mNGhostPair, &b_StMiniMcEvent_mNGhostPair);
   fChain->SetBranchAddress("mNContamPair", &mNContamPair, &b_StMiniMcEvent_mNContamPair);
   fChain->SetBranchAddress("mNMatGlobPair", &mNMatGlobPair, &b_StMiniMcEvent_mNMatGlobPair);
   fChain->SetBranchAddress("mMcTracks_", &mMcTracks_, &b_mMcTracks_);
   fChain->SetBranchAddress("mMcTracks.mIsValid", mMcTracks_mIsValid, &b_mMcTracks_mIsValid);
   fChain->SetBranchAddress("mMcTracks.mPtMc", mMcTracks_mPtMc, &b_mMcTracks_mPtMc);
   fChain->SetBranchAddress("mMcTracks.mPzMc", mMcTracks_mPzMc, &b_mMcTracks_mPzMc);
   fChain->SetBranchAddress("mMcTracks.mEtaMc", mMcTracks_mEtaMc, &b_mMcTracks_mEtaMc);
   fChain->SetBranchAddress("mMcTracks.mPhiMc", mMcTracks_mPhiMc, &b_mMcTracks_mPhiMc);
   fChain->SetBranchAddress("mMcTracks.mNHitMc", mMcTracks_mNHitMc, &b_mMcTracks_mNHitMc);
   fChain->SetBranchAddress("mMcTracks.mNSvtHitMc", mMcTracks_mNSvtHitMc, &b_mMcTracks_mNSvtHitMc);
   fChain->SetBranchAddress("mMcTracks.mNSsdHitMc", mMcTracks_mNSsdHitMc, &b_mMcTracks_mNSsdHitMc);
   fChain->SetBranchAddress("mMcTracks.mNFtpcHitMc", mMcTracks_mNFtpcHitMc, &b_mMcTracks_mNFtpcHitMc);
   fChain->SetBranchAddress("mMcTracks.mNBemcHitMc", mMcTracks_mNBemcHitMc, &b_mMcTracks_mNBemcHitMc);
   fChain->SetBranchAddress("mMcTracks.mNBprsHitMc", mMcTracks_mNBprsHitMc, &b_mMcTracks_mNBprsHitMc);
   fChain->SetBranchAddress("mMcTracks.mNBsmdeHitMc", mMcTracks_mNBsmdeHitMc, &b_mMcTracks_mNBsmdeHitMc);
   fChain->SetBranchAddress("mMcTracks.mNBsmdpHitMc", mMcTracks_mNBsmdpHitMc, &b_mMcTracks_mNBsmdpHitMc);
   fChain->SetBranchAddress("mMcTracks.mNEemcHitMc", mMcTracks_mNEemcHitMc, &b_mMcTracks_mNEemcHitMc);
   fChain->SetBranchAddress("mMcTracks.mNEprsHitMc", mMcTracks_mNEprsHitMc, &b_mMcTracks_mNEprsHitMc);
   fChain->SetBranchAddress("mMcTracks.mNEsmduHitMc", mMcTracks_mNEsmduHitMc, &b_mMcTracks_mNEsmduHitMc);
   fChain->SetBranchAddress("mMcTracks.mNEsmdvHitMc", mMcTracks_mNEsmdvHitMc, &b_mMcTracks_mNEsmdvHitMc);
   fChain->SetBranchAddress("mMcTracks.mGeantId", mMcTracks_mGeantId, &b_mMcTracks_mGeantId);
   fChain->SetBranchAddress("mMcTracks.mChargeMc", mMcTracks_mChargeMc, &b_mMcTracks_mChargeMc);
   fChain->SetBranchAddress("mMcTracks.mStopR", mMcTracks_mStopR, &b_mMcTracks_mStopR);
   fChain->SetBranchAddress("mMcTracks.mKey", mMcTracks_mKey, &b_mMcTracks_mKey);
   fChain->SetBranchAddress("mMcTracks.mParentKey", mMcTracks_mParentKey, &b_mMcTracks_mParentKey);
   fChain->SetBranchAddress("mMcTracks.mParentGeantId", mMcTracks_mParentGeantId, &b_mMcTracks_mParentGeantId);
   fChain->SetBranchAddress("mMcTracks.mEmcEnergyMcHit[3]", mMcTracks_mEmcEnergyMcHit, &b_mMcTracks_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mMcTracks.mEmcEnergyMcSum", mMcTracks_mEmcEnergyMcSum, &b_mMcTracks_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mMcTracks.mEmcSoftIdHiTowerMc[3]", mMcTracks_mEmcSoftIdHiTowerMc, &b_mMcTracks_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mMcTracks.mNAssocGl", mMcTracks_mNAssocGl, &b_mMcTracks_mNAssocGl);
   fChain->SetBranchAddress("mMcTracks.mNAssocPr", mMcTracks_mNAssocPr, &b_mMcTracks_mNAssocPr);
   fChain->SetBranchAddress("mMcTracks.mIsPrimary", mMcTracks_mIsPrimary, &b_mMcTracks_mIsPrimary);
   fChain->SetBranchAddress("mMcTracks.fUniqueID", mMcTracks_fUniqueID, &b_mMcTracks_fUniqueID);
   fChain->SetBranchAddress("mMcTracks.fBits", mMcTracks_fBits, &b_mMcTracks_fBits);
   fChain->SetBranchAddress("mMatchedPairs_", &mMatchedPairs_, &b_mMatchedPairs_);
   fChain->SetBranchAddress("mMatchedPairs.mNCommonHit", mMatchedPairs_mNCommonHit, &b_mMatchedPairs_mNCommonHit);
   fChain->SetBranchAddress("mMatchedPairs.mIsBestContam", mMatchedPairs_mIsBestContam, &b_mMatchedPairs_mIsBestContam);
   fChain->SetBranchAddress("mMatchedPairs.mDominatrack", mMatchedPairs_mDominatrack, &b_mMatchedPairs_mDominatrack);
   fChain->SetBranchAddress("mMatchedPairs.mDominCommonHit", mMatchedPairs_mDominCommonHit, &b_mMatchedPairs_mDominCommonHit);
   fChain->SetBranchAddress("mMatchedPairs.mAvgQuality", mMatchedPairs_mAvgQuality, &b_mMatchedPairs_mAvgQuality);
   fChain->SetBranchAddress("mMatchedPairs.mIsValid", mMatchedPairs_mIsValid, &b_mMatchedPairs_mIsValid);
   fChain->SetBranchAddress("mMatchedPairs.mPtMc", mMatchedPairs_mPtMc, &b_mMatchedPairs_mPtMc);
   fChain->SetBranchAddress("mMatchedPairs.mPzMc", mMatchedPairs_mPzMc, &b_mMatchedPairs_mPzMc);
   fChain->SetBranchAddress("mMatchedPairs.mEtaMc", mMatchedPairs_mEtaMc, &b_mMatchedPairs_mEtaMc);
   fChain->SetBranchAddress("mMatchedPairs.mPhiMc", mMatchedPairs_mPhiMc, &b_mMatchedPairs_mPhiMc);
   fChain->SetBranchAddress("mMatchedPairs.mNHitMc", mMatchedPairs_mNHitMc, &b_mMatchedPairs_mNHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNSvtHitMc", mMatchedPairs_mNSvtHitMc, &b_mMatchedPairs_mNSvtHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNSsdHitMc", mMatchedPairs_mNSsdHitMc, &b_mMatchedPairs_mNSsdHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNFtpcHitMc", mMatchedPairs_mNFtpcHitMc, &b_mMatchedPairs_mNFtpcHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNBemcHitMc", mMatchedPairs_mNBemcHitMc, &b_mMatchedPairs_mNBemcHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNBprsHitMc", mMatchedPairs_mNBprsHitMc, &b_mMatchedPairs_mNBprsHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNBsmdeHitMc", mMatchedPairs_mNBsmdeHitMc, &b_mMatchedPairs_mNBsmdeHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNBsmdpHitMc", mMatchedPairs_mNBsmdpHitMc, &b_mMatchedPairs_mNBsmdpHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNEemcHitMc", mMatchedPairs_mNEemcHitMc, &b_mMatchedPairs_mNEemcHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNEprsHitMc", mMatchedPairs_mNEprsHitMc, &b_mMatchedPairs_mNEprsHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNEsmduHitMc", mMatchedPairs_mNEsmduHitMc, &b_mMatchedPairs_mNEsmduHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mNEsmdvHitMc", mMatchedPairs_mNEsmdvHitMc, &b_mMatchedPairs_mNEsmdvHitMc);
   fChain->SetBranchAddress("mMatchedPairs.mGeantId", mMatchedPairs_mGeantId, &b_mMatchedPairs_mGeantId);
   fChain->SetBranchAddress("mMatchedPairs.mChargeMc", mMatchedPairs_mChargeMc, &b_mMatchedPairs_mChargeMc);
   fChain->SetBranchAddress("mMatchedPairs.mStopR", mMatchedPairs_mStopR, &b_mMatchedPairs_mStopR);
   fChain->SetBranchAddress("mMatchedPairs.mKey", mMatchedPairs_mKey, &b_mMatchedPairs_mKey);
   fChain->SetBranchAddress("mMatchedPairs.mParentKey", mMatchedPairs_mParentKey, &b_mMatchedPairs_mParentKey);
   fChain->SetBranchAddress("mMatchedPairs.mParentGeantId", mMatchedPairs_mParentGeantId, &b_mMatchedPairs_mParentGeantId);
   fChain->SetBranchAddress("mMatchedPairs.mEmcEnergyMcHit[3]", mMatchedPairs_mEmcEnergyMcHit, &b_mMatchedPairs_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mMatchedPairs.mEmcEnergyMcSum", mMatchedPairs_mEmcEnergyMcSum, &b_mMatchedPairs_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mMatchedPairs.mEmcSoftIdHiTowerMc[3]", mMatchedPairs_mEmcSoftIdHiTowerMc, &b_mMatchedPairs_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mMatchedPairs.mNAssocGl", mMatchedPairs_mNAssocGl, &b_mMatchedPairs_mNAssocGl);
   fChain->SetBranchAddress("mMatchedPairs.mNAssocPr", mMatchedPairs_mNAssocPr, &b_mMatchedPairs_mNAssocPr);
   fChain->SetBranchAddress("mMatchedPairs.mIsPrimary", mMatchedPairs_mIsPrimary, &b_mMatchedPairs_mIsPrimary);
   fChain->SetBranchAddress("mMatchedPairs.fUniqueID", mMatchedPairs_fUniqueID, &b_mMatchedPairs_fUniqueID);
   fChain->SetBranchAddress("mMatchedPairs.fBits", mMatchedPairs_fBits, &b_mMatchedPairs_fBits);
   fChain->SetBranchAddress("mMatchedPairs.mIsValidGl", mMatchedPairs_mIsValidGl, &b_mMatchedPairs_mIsValidGl);
   fChain->SetBranchAddress("mMatchedPairs.mPtPr", mMatchedPairs_mPtPr, &b_mMatchedPairs_mPtPr);
   fChain->SetBranchAddress("mMatchedPairs.mPzPr", mMatchedPairs_mPzPr, &b_mMatchedPairs_mPzPr);
   fChain->SetBranchAddress("mMatchedPairs.mEtaPr", mMatchedPairs_mEtaPr, &b_mMatchedPairs_mEtaPr);
   fChain->SetBranchAddress("mMatchedPairs.mPhiPr", mMatchedPairs_mPhiPr, &b_mMatchedPairs_mPhiPr);
   fChain->SetBranchAddress("mMatchedPairs.mDcaPr", mMatchedPairs_mDcaPr, &b_mMatchedPairs_mDcaPr);
   fChain->SetBranchAddress("mMatchedPairs.mDcaXYPr", mMatchedPairs_mDcaXYPr, &b_mMatchedPairs_mDcaXYPr);
   fChain->SetBranchAddress("mMatchedPairs.mDcaZPr", mMatchedPairs_mDcaZPr, &b_mMatchedPairs_mDcaZPr);
   fChain->SetBranchAddress("mMatchedPairs.mDcaXYPrMcV", mMatchedPairs_mDcaXYPrMcV, &b_mMatchedPairs_mDcaXYPrMcV);
   fChain->SetBranchAddress("mMatchedPairs.mDcaZPrMcV", mMatchedPairs_mDcaZPrMcV, &b_mMatchedPairs_mDcaZPrMcV);
   fChain->SetBranchAddress("mMatchedPairs.mCurvPr", mMatchedPairs_mCurvPr, &b_mMatchedPairs_mCurvPr);
   fChain->SetBranchAddress("mMatchedPairs.mTanLPr", mMatchedPairs_mTanLPr, &b_mMatchedPairs_mTanLPr);
   fChain->SetBranchAddress("mMatchedPairs.mErrP[5]", mMatchedPairs_mErrP, &b_mMatchedPairs_mErrP);
   fChain->SetBranchAddress("mMatchedPairs.mChi2Pr", mMatchedPairs_mChi2Pr, &b_mMatchedPairs_mChi2Pr);
   fChain->SetBranchAddress("mMatchedPairs.mFlag", mMatchedPairs_mFlag, &b_mMatchedPairs_mFlag);
   fChain->SetBranchAddress("mMatchedPairs.mDedx", mMatchedPairs_mDedx, &b_mMatchedPairs_mDedx);
   fChain->SetBranchAddress("mMatchedPairs.mPtGl", mMatchedPairs_mPtGl, &b_mMatchedPairs_mPtGl);
   fChain->SetBranchAddress("mMatchedPairs.mPzGl", mMatchedPairs_mPzGl, &b_mMatchedPairs_mPzGl);
   fChain->SetBranchAddress("mMatchedPairs.mEtaGl", mMatchedPairs_mEtaGl, &b_mMatchedPairs_mEtaGl);
   fChain->SetBranchAddress("mMatchedPairs.mPhiGl", mMatchedPairs_mPhiGl, &b_mMatchedPairs_mPhiGl);
   fChain->SetBranchAddress("mMatchedPairs.mDcaGl", mMatchedPairs_mDcaGl, &b_mMatchedPairs_mDcaGl);
   fChain->SetBranchAddress("mMatchedPairs.mDcaXYGl", mMatchedPairs_mDcaXYGl, &b_mMatchedPairs_mDcaXYGl);
   fChain->SetBranchAddress("mMatchedPairs.mDcaZGl", mMatchedPairs_mDcaZGl, &b_mMatchedPairs_mDcaZGl);
   fChain->SetBranchAddress("mMatchedPairs.mDcaXYGlMcV", mMatchedPairs_mDcaXYGlMcV, &b_mMatchedPairs_mDcaXYGlMcV);
   fChain->SetBranchAddress("mMatchedPairs.mDcaZGlMcV", mMatchedPairs_mDcaZGlMcV, &b_mMatchedPairs_mDcaZGlMcV);
   fChain->SetBranchAddress("mMatchedPairs.mCurvGl", mMatchedPairs_mCurvGl, &b_mMatchedPairs_mCurvGl);
   fChain->SetBranchAddress("mMatchedPairs.mTanLGl", mMatchedPairs_mTanLGl, &b_mMatchedPairs_mTanLGl);
   fChain->SetBranchAddress("mMatchedPairs.mErrG[5]", mMatchedPairs_mErrG, &b_mMatchedPairs_mErrG);
   fChain->SetBranchAddress("mMatchedPairs.mSeedQA", mMatchedPairs_mSeedQA, &b_mMatchedPairs_mSeedQA);
   fChain->SetBranchAddress("mMatchedPairs.mPidPion", mMatchedPairs_mPidPion, &b_mMatchedPairs_mPidPion);
   fChain->SetBranchAddress("mMatchedPairs.mPidProton", mMatchedPairs_mPidProton, &b_mMatchedPairs_mPidProton);
   fChain->SetBranchAddress("mMatchedPairs.mPidKaon", mMatchedPairs_mPidKaon, &b_mMatchedPairs_mPidKaon);
   fChain->SetBranchAddress("mMatchedPairs.mPidElectron", mMatchedPairs_mPidElectron, &b_mMatchedPairs_mPidElectron);
   fChain->SetBranchAddress("mMatchedPairs.mFirstZ", mMatchedPairs_mFirstZ, &b_mMatchedPairs_mFirstZ);
   fChain->SetBranchAddress("mMatchedPairs.mLastZ", mMatchedPairs_mLastZ, &b_mMatchedPairs_mLastZ);
   fChain->SetBranchAddress("mMatchedPairs.mFirstPadrow", mMatchedPairs_mFirstPadrow, &b_mMatchedPairs_mFirstPadrow);
   fChain->SetBranchAddress("mMatchedPairs.mLastPadrow", mMatchedPairs_mLastPadrow, &b_mMatchedPairs_mLastPadrow);
   fChain->SetBranchAddress("mMatchedPairs.mFirstFitPadrow", mMatchedPairs_mFirstFitPadrow, &b_mMatchedPairs_mFirstFitPadrow);
   fChain->SetBranchAddress("mMatchedPairs.mLastFitPadrow", mMatchedPairs_mLastFitPadrow, &b_mMatchedPairs_mLastFitPadrow);
   fChain->SetBranchAddress("mMatchedPairs.mFirstSector", mMatchedPairs_mFirstSector, &b_mMatchedPairs_mFirstSector);
   fChain->SetBranchAddress("mMatchedPairs.mLastSector", mMatchedPairs_mLastSector, &b_mMatchedPairs_mLastSector);
   fChain->SetBranchAddress("mMatchedPairs.mRecoKey", mMatchedPairs_mRecoKey, &b_mMatchedPairs_mRecoKey);
   fChain->SetBranchAddress("mMatchedPairs.mFitPts", mMatchedPairs_mFitPts, &b_mMatchedPairs_mFitPts);
   fChain->SetBranchAddress("mMatchedPairs.mFitSvt", mMatchedPairs_mFitSvt, &b_mMatchedPairs_mFitSvt);
   fChain->SetBranchAddress("mMatchedPairs.mFitSsd", mMatchedPairs_mFitSsd, &b_mMatchedPairs_mFitSsd);
   fChain->SetBranchAddress("mMatchedPairs.mFitFtpc", mMatchedPairs_mFitFtpc, &b_mMatchedPairs_mFitFtpc);
   fChain->SetBranchAddress("mMatchedPairs.mDedxPts", mMatchedPairs_mDedxPts, &b_mMatchedPairs_mDedxPts);
   fChain->SetBranchAddress("mMatchedPairs.mAllPts", mMatchedPairs_mAllPts, &b_mMatchedPairs_mAllPts);
   fChain->SetBranchAddress("mMatchedPairs.mCharge", mMatchedPairs_mCharge, &b_mMatchedPairs_mCharge);
   fChain->SetBranchAddress("mMatchedPairs.mNAssocMc", mMatchedPairs_mNAssocMc, &b_mMatchedPairs_mNAssocMc);
   fChain->SetBranchAddress("mMatchedPairs.mNPossible", mMatchedPairs_mNPossible, &b_mMatchedPairs_mNPossible);
   fChain->SetBranchAddress("mMatchedPairs.mIsValidPr", mMatchedPairs_mIsValidPr, &b_mMatchedPairs_mIsValidPr);
   fChain->SetBranchAddress("mMatchedPairs.mEmcTowerAdc[3]", mMatchedPairs_mEmcTowerAdc, &b_mMatchedPairs_mEmcTowerAdc);
   fChain->SetBranchAddress("mMatchedPairs.mEmcEnergyRcHit[3]", mMatchedPairs_mEmcEnergyRcHit, &b_mMatchedPairs_mEmcEnergyRcHit);
   fChain->SetBranchAddress("mMatchedPairs.mEmcSoftIdHiTowerRc[3]", mMatchedPairs_mEmcSoftIdHiTowerRc, &b_mMatchedPairs_mEmcSoftIdHiTowerRc);
   fChain->SetBranchAddress("mMergedPairs_", &mMergedPairs_, &b_mMergedPairs_);
#if 0
   fChain->SetBranchAddress("mMergedPairs.mNCommonHit", mMergedPairs_mNCommonHit, &b_mMergedPairs_mNCommonHit);
   fChain->SetBranchAddress("mMergedPairs.mIsBestContam", mMergedPairs_mIsBestContam, &b_mMergedPairs_mIsBestContam);
   fChain->SetBranchAddress("mMergedPairs.mDominatrack", mMergedPairs_mDominatrack, &b_mMergedPairs_mDominatrack);
   fChain->SetBranchAddress("mMergedPairs.mDominCommonHit", mMergedPairs_mDominCommonHit, &b_mMergedPairs_mDominCommonHit);
   fChain->SetBranchAddress("mMergedPairs.mAvgQuality", mMergedPairs_mAvgQuality, &b_mMergedPairs_mAvgQuality);
   fChain->SetBranchAddress("mMergedPairs.mIsValid", mMergedPairs_mIsValid, &b_mMergedPairs_mIsValid);
   fChain->SetBranchAddress("mMergedPairs.mPtMc", mMergedPairs_mPtMc, &b_mMergedPairs_mPtMc);
   fChain->SetBranchAddress("mMergedPairs.mPzMc", mMergedPairs_mPzMc, &b_mMergedPairs_mPzMc);
   fChain->SetBranchAddress("mMergedPairs.mEtaMc", mMergedPairs_mEtaMc, &b_mMergedPairs_mEtaMc);
   fChain->SetBranchAddress("mMergedPairs.mPhiMc", mMergedPairs_mPhiMc, &b_mMergedPairs_mPhiMc);
   fChain->SetBranchAddress("mMergedPairs.mNHitMc", mMergedPairs_mNHitMc, &b_mMergedPairs_mNHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNSvtHitMc", mMergedPairs_mNSvtHitMc, &b_mMergedPairs_mNSvtHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNSsdHitMc", mMergedPairs_mNSsdHitMc, &b_mMergedPairs_mNSsdHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNFtpcHitMc", mMergedPairs_mNFtpcHitMc, &b_mMergedPairs_mNFtpcHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNBemcHitMc", mMergedPairs_mNBemcHitMc, &b_mMergedPairs_mNBemcHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNBprsHitMc", mMergedPairs_mNBprsHitMc, &b_mMergedPairs_mNBprsHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNBsmdeHitMc", mMergedPairs_mNBsmdeHitMc, &b_mMergedPairs_mNBsmdeHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNBsmdpHitMc", mMergedPairs_mNBsmdpHitMc, &b_mMergedPairs_mNBsmdpHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNEemcHitMc", mMergedPairs_mNEemcHitMc, &b_mMergedPairs_mNEemcHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNEprsHitMc", mMergedPairs_mNEprsHitMc, &b_mMergedPairs_mNEprsHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNEsmduHitMc", mMergedPairs_mNEsmduHitMc, &b_mMergedPairs_mNEsmduHitMc);
   fChain->SetBranchAddress("mMergedPairs.mNEsmdvHitMc", mMergedPairs_mNEsmdvHitMc, &b_mMergedPairs_mNEsmdvHitMc);
   fChain->SetBranchAddress("mMergedPairs.mGeantId", mMergedPairs_mGeantId, &b_mMergedPairs_mGeantId);
   fChain->SetBranchAddress("mMergedPairs.mChargeMc", mMergedPairs_mChargeMc, &b_mMergedPairs_mChargeMc);
   fChain->SetBranchAddress("mMergedPairs.mStopR", mMergedPairs_mStopR, &b_mMergedPairs_mStopR);
   fChain->SetBranchAddress("mMergedPairs.mKey", mMergedPairs_mKey, &b_mMergedPairs_mKey);
   fChain->SetBranchAddress("mMergedPairs.mParentKey", mMergedPairs_mParentKey, &b_mMergedPairs_mParentKey);
   fChain->SetBranchAddress("mMergedPairs.mParentGeantId", mMergedPairs_mParentGeantId, &b_mMergedPairs_mParentGeantId);
   fChain->SetBranchAddress("mMergedPairs.mEmcEnergyMcHit[3]", mMergedPairs_mEmcEnergyMcHit, &b_mMergedPairs_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mMergedPairs.mEmcEnergyMcSum", mMergedPairs_mEmcEnergyMcSum, &b_mMergedPairs_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mMergedPairs.mEmcSoftIdHiTowerMc[3]", mMergedPairs_mEmcSoftIdHiTowerMc, &b_mMergedPairs_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mMergedPairs.mNAssocGl", mMergedPairs_mNAssocGl, &b_mMergedPairs_mNAssocGl);
   fChain->SetBranchAddress("mMergedPairs.mNAssocPr", mMergedPairs_mNAssocPr, &b_mMergedPairs_mNAssocPr);
   fChain->SetBranchAddress("mMergedPairs.mIsPrimary", mMergedPairs_mIsPrimary, &b_mMergedPairs_mIsPrimary);
   fChain->SetBranchAddress("mMergedPairs.fUniqueID", mMergedPairs_fUniqueID, &b_mMergedPairs_fUniqueID);
   fChain->SetBranchAddress("mMergedPairs.fBits", mMergedPairs_fBits, &b_mMergedPairs_fBits);
   fChain->SetBranchAddress("mMergedPairs.mIsValidGl", mMergedPairs_mIsValidGl, &b_mMergedPairs_mIsValidGl);
   fChain->SetBranchAddress("mMergedPairs.mPtPr", mMergedPairs_mPtPr, &b_mMergedPairs_mPtPr);
   fChain->SetBranchAddress("mMergedPairs.mPzPr", mMergedPairs_mPzPr, &b_mMergedPairs_mPzPr);
   fChain->SetBranchAddress("mMergedPairs.mEtaPr", mMergedPairs_mEtaPr, &b_mMergedPairs_mEtaPr);
   fChain->SetBranchAddress("mMergedPairs.mPhiPr", mMergedPairs_mPhiPr, &b_mMergedPairs_mPhiPr);
   fChain->SetBranchAddress("mMergedPairs.mDcaPr", mMergedPairs_mDcaPr, &b_mMergedPairs_mDcaPr);
   fChain->SetBranchAddress("mMergedPairs.mDcaXYPr", mMergedPairs_mDcaXYPr, &b_mMergedPairs_mDcaXYPr);
   fChain->SetBranchAddress("mMergedPairs.mDcaZPr", mMergedPairs_mDcaZPr, &b_mMergedPairs_mDcaZPr);
   fChain->SetBranchAddress("mMergedPairs.mDcaXYPrMcV", mMergedPairs_mDcaXYPrMcV, &b_mMergedPairs_mDcaXYPrMcV);
   fChain->SetBranchAddress("mMergedPairs.mDcaZPrMcV", mMergedPairs_mDcaZPrMcV, &b_mMergedPairs_mDcaZPrMcV);
   fChain->SetBranchAddress("mMergedPairs.mCurvPr", mMergedPairs_mCurvPr, &b_mMergedPairs_mCurvPr);
   fChain->SetBranchAddress("mMergedPairs.mTanLPr", mMergedPairs_mTanLPr, &b_mMergedPairs_mTanLPr);
   fChain->SetBranchAddress("mMergedPairs.mErrP[5]", mMergedPairs_mErrP, &b_mMergedPairs_mErrP);
   fChain->SetBranchAddress("mMergedPairs.mChi2Pr", mMergedPairs_mChi2Pr, &b_mMergedPairs_mChi2Pr);
   fChain->SetBranchAddress("mMergedPairs.mFlag", mMergedPairs_mFlag, &b_mMergedPairs_mFlag);
   fChain->SetBranchAddress("mMergedPairs.mDedx", mMergedPairs_mDedx, &b_mMergedPairs_mDedx);
   fChain->SetBranchAddress("mMergedPairs.mPtGl", mMergedPairs_mPtGl, &b_mMergedPairs_mPtGl);
   fChain->SetBranchAddress("mMergedPairs.mPzGl", mMergedPairs_mPzGl, &b_mMergedPairs_mPzGl);
   fChain->SetBranchAddress("mMergedPairs.mEtaGl", mMergedPairs_mEtaGl, &b_mMergedPairs_mEtaGl);
   fChain->SetBranchAddress("mMergedPairs.mPhiGl", mMergedPairs_mPhiGl, &b_mMergedPairs_mPhiGl);
   fChain->SetBranchAddress("mMergedPairs.mDcaGl", mMergedPairs_mDcaGl, &b_mMergedPairs_mDcaGl);
   fChain->SetBranchAddress("mMergedPairs.mDcaXYGl", mMergedPairs_mDcaXYGl, &b_mMergedPairs_mDcaXYGl);
   fChain->SetBranchAddress("mMergedPairs.mDcaZGl", mMergedPairs_mDcaZGl, &b_mMergedPairs_mDcaZGl);
   fChain->SetBranchAddress("mMergedPairs.mDcaXYGlMcV", mMergedPairs_mDcaXYGlMcV, &b_mMergedPairs_mDcaXYGlMcV);
   fChain->SetBranchAddress("mMergedPairs.mDcaZGlMcV", mMergedPairs_mDcaZGlMcV, &b_mMergedPairs_mDcaZGlMcV);
   fChain->SetBranchAddress("mMergedPairs.mCurvGl", mMergedPairs_mCurvGl, &b_mMergedPairs_mCurvGl);
   fChain->SetBranchAddress("mMergedPairs.mTanLGl", mMergedPairs_mTanLGl, &b_mMergedPairs_mTanLGl);
   fChain->SetBranchAddress("mMergedPairs.mErrG[5]", mMergedPairs_mErrG, &b_mMergedPairs_mErrG);
   fChain->SetBranchAddress("mMergedPairs.mSeedQA", mMergedPairs_mSeedQA, &b_mMergedPairs_mSeedQA);
   fChain->SetBranchAddress("mMergedPairs.mPidPion", mMergedPairs_mPidPion, &b_mMergedPairs_mPidPion);
   fChain->SetBranchAddress("mMergedPairs.mPidProton", mMergedPairs_mPidProton, &b_mMergedPairs_mPidProton);
   fChain->SetBranchAddress("mMergedPairs.mPidKaon", mMergedPairs_mPidKaon, &b_mMergedPairs_mPidKaon);
   fChain->SetBranchAddress("mMergedPairs.mPidElectron", mMergedPairs_mPidElectron, &b_mMergedPairs_mPidElectron);
   fChain->SetBranchAddress("mMergedPairs.mFirstZ", mMergedPairs_mFirstZ, &b_mMergedPairs_mFirstZ);
   fChain->SetBranchAddress("mMergedPairs.mLastZ", mMergedPairs_mLastZ, &b_mMergedPairs_mLastZ);
   fChain->SetBranchAddress("mMergedPairs.mFirstPadrow", mMergedPairs_mFirstPadrow, &b_mMergedPairs_mFirstPadrow);
   fChain->SetBranchAddress("mMergedPairs.mLastPadrow", mMergedPairs_mLastPadrow, &b_mMergedPairs_mLastPadrow);
   fChain->SetBranchAddress("mMergedPairs.mFirstFitPadrow", mMergedPairs_mFirstFitPadrow, &b_mMergedPairs_mFirstFitPadrow);
   fChain->SetBranchAddress("mMergedPairs.mLastFitPadrow", mMergedPairs_mLastFitPadrow, &b_mMergedPairs_mLastFitPadrow);
   fChain->SetBranchAddress("mMergedPairs.mFirstSector", mMergedPairs_mFirstSector, &b_mMergedPairs_mFirstSector);
   fChain->SetBranchAddress("mMergedPairs.mLastSector", mMergedPairs_mLastSector, &b_mMergedPairs_mLastSector);
   fChain->SetBranchAddress("mMergedPairs.mRecoKey", mMergedPairs_mRecoKey, &b_mMergedPairs_mRecoKey);
   fChain->SetBranchAddress("mMergedPairs.mFitPts", mMergedPairs_mFitPts, &b_mMergedPairs_mFitPts);
   fChain->SetBranchAddress("mMergedPairs.mFitSvt", mMergedPairs_mFitSvt, &b_mMergedPairs_mFitSvt);
   fChain->SetBranchAddress("mMergedPairs.mFitSsd", mMergedPairs_mFitSsd, &b_mMergedPairs_mFitSsd);
   fChain->SetBranchAddress("mMergedPairs.mFitFtpc", mMergedPairs_mFitFtpc, &b_mMergedPairs_mFitFtpc);
   fChain->SetBranchAddress("mMergedPairs.mDedxPts", mMergedPairs_mDedxPts, &b_mMergedPairs_mDedxPts);
   fChain->SetBranchAddress("mMergedPairs.mAllPts", mMergedPairs_mAllPts, &b_mMergedPairs_mAllPts);
   fChain->SetBranchAddress("mMergedPairs.mCharge", mMergedPairs_mCharge, &b_mMergedPairs_mCharge);
   fChain->SetBranchAddress("mMergedPairs.mNAssocMc", mMergedPairs_mNAssocMc, &b_mMergedPairs_mNAssocMc);
   fChain->SetBranchAddress("mMergedPairs.mNPossible", mMergedPairs_mNPossible, &b_mMergedPairs_mNPossible);
   fChain->SetBranchAddress("mMergedPairs.mIsValidPr", mMergedPairs_mIsValidPr, &b_mMergedPairs_mIsValidPr);
   fChain->SetBranchAddress("mMergedPairs.mEmcTowerAdc[3]", mMergedPairs_mEmcTowerAdc, &b_mMergedPairs_mEmcTowerAdc);
   fChain->SetBranchAddress("mMergedPairs.mEmcEnergyRcHit[3]", mMergedPairs_mEmcEnergyRcHit, &b_mMergedPairs_mEmcEnergyRcHit);
   fChain->SetBranchAddress("mMergedPairs.mEmcSoftIdHiTowerRc[3]", mMergedPairs_mEmcSoftIdHiTowerRc, &b_mMergedPairs_mEmcSoftIdHiTowerRc);
   fChain->SetBranchAddress("mSplitPairs_", &mSplitPairs_, &b_mSplitPairs_);
   fChain->SetBranchAddress("mSplitPairs.mNCommonHit", mSplitPairs_mNCommonHit, &b_mSplitPairs_mNCommonHit);
   fChain->SetBranchAddress("mSplitPairs.mIsBestContam", mSplitPairs_mIsBestContam, &b_mSplitPairs_mIsBestContam);
   fChain->SetBranchAddress("mSplitPairs.mDominatrack", mSplitPairs_mDominatrack, &b_mSplitPairs_mDominatrack);
   fChain->SetBranchAddress("mSplitPairs.mDominCommonHit", mSplitPairs_mDominCommonHit, &b_mSplitPairs_mDominCommonHit);
   fChain->SetBranchAddress("mSplitPairs.mAvgQuality", mSplitPairs_mAvgQuality, &b_mSplitPairs_mAvgQuality);
   fChain->SetBranchAddress("mSplitPairs.mIsValid", mSplitPairs_mIsValid, &b_mSplitPairs_mIsValid);
   fChain->SetBranchAddress("mSplitPairs.mPtMc", mSplitPairs_mPtMc, &b_mSplitPairs_mPtMc);
   fChain->SetBranchAddress("mSplitPairs.mPzMc", mSplitPairs_mPzMc, &b_mSplitPairs_mPzMc);
   fChain->SetBranchAddress("mSplitPairs.mEtaMc", mSplitPairs_mEtaMc, &b_mSplitPairs_mEtaMc);
   fChain->SetBranchAddress("mSplitPairs.mPhiMc", mSplitPairs_mPhiMc, &b_mSplitPairs_mPhiMc);
   fChain->SetBranchAddress("mSplitPairs.mNHitMc", mSplitPairs_mNHitMc, &b_mSplitPairs_mNHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNSvtHitMc", mSplitPairs_mNSvtHitMc, &b_mSplitPairs_mNSvtHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNSsdHitMc", mSplitPairs_mNSsdHitMc, &b_mSplitPairs_mNSsdHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNFtpcHitMc", mSplitPairs_mNFtpcHitMc, &b_mSplitPairs_mNFtpcHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNBemcHitMc", mSplitPairs_mNBemcHitMc, &b_mSplitPairs_mNBemcHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNBprsHitMc", mSplitPairs_mNBprsHitMc, &b_mSplitPairs_mNBprsHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNBsmdeHitMc", mSplitPairs_mNBsmdeHitMc, &b_mSplitPairs_mNBsmdeHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNBsmdpHitMc", mSplitPairs_mNBsmdpHitMc, &b_mSplitPairs_mNBsmdpHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNEemcHitMc", mSplitPairs_mNEemcHitMc, &b_mSplitPairs_mNEemcHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNEprsHitMc", mSplitPairs_mNEprsHitMc, &b_mSplitPairs_mNEprsHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNEsmduHitMc", mSplitPairs_mNEsmduHitMc, &b_mSplitPairs_mNEsmduHitMc);
   fChain->SetBranchAddress("mSplitPairs.mNEsmdvHitMc", mSplitPairs_mNEsmdvHitMc, &b_mSplitPairs_mNEsmdvHitMc);
   fChain->SetBranchAddress("mSplitPairs.mGeantId", mSplitPairs_mGeantId, &b_mSplitPairs_mGeantId);
   fChain->SetBranchAddress("mSplitPairs.mChargeMc", mSplitPairs_mChargeMc, &b_mSplitPairs_mChargeMc);
   fChain->SetBranchAddress("mSplitPairs.mStopR", mSplitPairs_mStopR, &b_mSplitPairs_mStopR);
   fChain->SetBranchAddress("mSplitPairs.mKey", mSplitPairs_mKey, &b_mSplitPairs_mKey);
   fChain->SetBranchAddress("mSplitPairs.mParentKey", mSplitPairs_mParentKey, &b_mSplitPairs_mParentKey);
   fChain->SetBranchAddress("mSplitPairs.mParentGeantId", mSplitPairs_mParentGeantId, &b_mSplitPairs_mParentGeantId);
   fChain->SetBranchAddress("mSplitPairs.mEmcEnergyMcHit[3]", mSplitPairs_mEmcEnergyMcHit, &b_mSplitPairs_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mSplitPairs.mEmcEnergyMcSum", mSplitPairs_mEmcEnergyMcSum, &b_mSplitPairs_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mSplitPairs.mEmcSoftIdHiTowerMc[3]", mSplitPairs_mEmcSoftIdHiTowerMc, &b_mSplitPairs_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mSplitPairs.mNAssocGl", mSplitPairs_mNAssocGl, &b_mSplitPairs_mNAssocGl);
   fChain->SetBranchAddress("mSplitPairs.mNAssocPr", mSplitPairs_mNAssocPr, &b_mSplitPairs_mNAssocPr);
   fChain->SetBranchAddress("mSplitPairs.mIsPrimary", mSplitPairs_mIsPrimary, &b_mSplitPairs_mIsPrimary);
   fChain->SetBranchAddress("mSplitPairs.fUniqueID", mSplitPairs_fUniqueID, &b_mSplitPairs_fUniqueID);
   fChain->SetBranchAddress("mSplitPairs.fBits", mSplitPairs_fBits, &b_mSplitPairs_fBits);
   fChain->SetBranchAddress("mSplitPairs.mIsValidGl", mSplitPairs_mIsValidGl, &b_mSplitPairs_mIsValidGl);
   fChain->SetBranchAddress("mSplitPairs.mPtPr", mSplitPairs_mPtPr, &b_mSplitPairs_mPtPr);
   fChain->SetBranchAddress("mSplitPairs.mPzPr", mSplitPairs_mPzPr, &b_mSplitPairs_mPzPr);
   fChain->SetBranchAddress("mSplitPairs.mEtaPr", mSplitPairs_mEtaPr, &b_mSplitPairs_mEtaPr);
   fChain->SetBranchAddress("mSplitPairs.mPhiPr", mSplitPairs_mPhiPr, &b_mSplitPairs_mPhiPr);
   fChain->SetBranchAddress("mSplitPairs.mDcaPr", mSplitPairs_mDcaPr, &b_mSplitPairs_mDcaPr);
   fChain->SetBranchAddress("mSplitPairs.mDcaXYPr", mSplitPairs_mDcaXYPr, &b_mSplitPairs_mDcaXYPr);
   fChain->SetBranchAddress("mSplitPairs.mDcaZPr", mSplitPairs_mDcaZPr, &b_mSplitPairs_mDcaZPr);
   fChain->SetBranchAddress("mSplitPairs.mDcaXYPrMcV", mSplitPairs_mDcaXYPrMcV, &b_mSplitPairs_mDcaXYPrMcV);
   fChain->SetBranchAddress("mSplitPairs.mDcaZPrMcV", mSplitPairs_mDcaZPrMcV, &b_mSplitPairs_mDcaZPrMcV);
   fChain->SetBranchAddress("mSplitPairs.mCurvPr", mSplitPairs_mCurvPr, &b_mSplitPairs_mCurvPr);
   fChain->SetBranchAddress("mSplitPairs.mTanLPr", mSplitPairs_mTanLPr, &b_mSplitPairs_mTanLPr);
   fChain->SetBranchAddress("mSplitPairs.mErrP[5]", mSplitPairs_mErrP, &b_mSplitPairs_mErrP);
   fChain->SetBranchAddress("mSplitPairs.mChi2Pr", mSplitPairs_mChi2Pr, &b_mSplitPairs_mChi2Pr);
   fChain->SetBranchAddress("mSplitPairs.mFlag", mSplitPairs_mFlag, &b_mSplitPairs_mFlag);
   fChain->SetBranchAddress("mSplitPairs.mDedx", mSplitPairs_mDedx, &b_mSplitPairs_mDedx);
   fChain->SetBranchAddress("mSplitPairs.mPtGl", mSplitPairs_mPtGl, &b_mSplitPairs_mPtGl);
   fChain->SetBranchAddress("mSplitPairs.mPzGl", mSplitPairs_mPzGl, &b_mSplitPairs_mPzGl);
   fChain->SetBranchAddress("mSplitPairs.mEtaGl", mSplitPairs_mEtaGl, &b_mSplitPairs_mEtaGl);
   fChain->SetBranchAddress("mSplitPairs.mPhiGl", mSplitPairs_mPhiGl, &b_mSplitPairs_mPhiGl);
   fChain->SetBranchAddress("mSplitPairs.mDcaGl", mSplitPairs_mDcaGl, &b_mSplitPairs_mDcaGl);
   fChain->SetBranchAddress("mSplitPairs.mDcaXYGl", mSplitPairs_mDcaXYGl, &b_mSplitPairs_mDcaXYGl);
   fChain->SetBranchAddress("mSplitPairs.mDcaZGl", mSplitPairs_mDcaZGl, &b_mSplitPairs_mDcaZGl);
   fChain->SetBranchAddress("mSplitPairs.mDcaXYGlMcV", mSplitPairs_mDcaXYGlMcV, &b_mSplitPairs_mDcaXYGlMcV);
   fChain->SetBranchAddress("mSplitPairs.mDcaZGlMcV", mSplitPairs_mDcaZGlMcV, &b_mSplitPairs_mDcaZGlMcV);
   fChain->SetBranchAddress("mSplitPairs.mCurvGl", mSplitPairs_mCurvGl, &b_mSplitPairs_mCurvGl);
   fChain->SetBranchAddress("mSplitPairs.mTanLGl", mSplitPairs_mTanLGl, &b_mSplitPairs_mTanLGl);
   fChain->SetBranchAddress("mSplitPairs.mErrG[5]", mSplitPairs_mErrG, &b_mSplitPairs_mErrG);
   fChain->SetBranchAddress("mSplitPairs.mSeedQA", mSplitPairs_mSeedQA, &b_mSplitPairs_mSeedQA);
   fChain->SetBranchAddress("mSplitPairs.mPidPion", mSplitPairs_mPidPion, &b_mSplitPairs_mPidPion);
   fChain->SetBranchAddress("mSplitPairs.mPidProton", mSplitPairs_mPidProton, &b_mSplitPairs_mPidProton);
   fChain->SetBranchAddress("mSplitPairs.mPidKaon", mSplitPairs_mPidKaon, &b_mSplitPairs_mPidKaon);
   fChain->SetBranchAddress("mSplitPairs.mPidElectron", mSplitPairs_mPidElectron, &b_mSplitPairs_mPidElectron);
   fChain->SetBranchAddress("mSplitPairs.mFirstZ", mSplitPairs_mFirstZ, &b_mSplitPairs_mFirstZ);
   fChain->SetBranchAddress("mSplitPairs.mLastZ", mSplitPairs_mLastZ, &b_mSplitPairs_mLastZ);
   fChain->SetBranchAddress("mSplitPairs.mFirstPadrow", mSplitPairs_mFirstPadrow, &b_mSplitPairs_mFirstPadrow);
   fChain->SetBranchAddress("mSplitPairs.mLastPadrow", mSplitPairs_mLastPadrow, &b_mSplitPairs_mLastPadrow);
   fChain->SetBranchAddress("mSplitPairs.mFirstFitPadrow", mSplitPairs_mFirstFitPadrow, &b_mSplitPairs_mFirstFitPadrow);
   fChain->SetBranchAddress("mSplitPairs.mLastFitPadrow", mSplitPairs_mLastFitPadrow, &b_mSplitPairs_mLastFitPadrow);
   fChain->SetBranchAddress("mSplitPairs.mFirstSector", mSplitPairs_mFirstSector, &b_mSplitPairs_mFirstSector);
   fChain->SetBranchAddress("mSplitPairs.mLastSector", mSplitPairs_mLastSector, &b_mSplitPairs_mLastSector);
   fChain->SetBranchAddress("mSplitPairs.mRecoKey", mSplitPairs_mRecoKey, &b_mSplitPairs_mRecoKey);
   fChain->SetBranchAddress("mSplitPairs.mFitPts", mSplitPairs_mFitPts, &b_mSplitPairs_mFitPts);
   fChain->SetBranchAddress("mSplitPairs.mFitSvt", mSplitPairs_mFitSvt, &b_mSplitPairs_mFitSvt);
   fChain->SetBranchAddress("mSplitPairs.mFitSsd", mSplitPairs_mFitSsd, &b_mSplitPairs_mFitSsd);
   fChain->SetBranchAddress("mSplitPairs.mFitFtpc", mSplitPairs_mFitFtpc, &b_mSplitPairs_mFitFtpc);
   fChain->SetBranchAddress("mSplitPairs.mDedxPts", mSplitPairs_mDedxPts, &b_mSplitPairs_mDedxPts);
   fChain->SetBranchAddress("mSplitPairs.mAllPts", mSplitPairs_mAllPts, &b_mSplitPairs_mAllPts);
   fChain->SetBranchAddress("mSplitPairs.mCharge", mSplitPairs_mCharge, &b_mSplitPairs_mCharge);
   fChain->SetBranchAddress("mSplitPairs.mNAssocMc", mSplitPairs_mNAssocMc, &b_mSplitPairs_mNAssocMc);
   fChain->SetBranchAddress("mSplitPairs.mNPossible", mSplitPairs_mNPossible, &b_mSplitPairs_mNPossible);
   fChain->SetBranchAddress("mSplitPairs.mIsValidPr", mSplitPairs_mIsValidPr, &b_mSplitPairs_mIsValidPr);
   fChain->SetBranchAddress("mSplitPairs.mEmcTowerAdc[3]", mSplitPairs_mEmcTowerAdc, &b_mSplitPairs_mEmcTowerAdc);
   fChain->SetBranchAddress("mSplitPairs.mEmcEnergyRcHit[3]", mSplitPairs_mEmcEnergyRcHit, &b_mSplitPairs_mEmcEnergyRcHit);
   fChain->SetBranchAddress("mSplitPairs.mEmcSoftIdHiTowerRc[3]", mSplitPairs_mEmcSoftIdHiTowerRc, &b_mSplitPairs_mEmcSoftIdHiTowerRc);
   fChain->SetBranchAddress("mGhostPairs_", &mGhostPairs_, &b_mGhostPairs_);
   fChain->SetBranchAddress("mGhostPairs.mNCommonHit", mGhostPairs_mNCommonHit, &b_mGhostPairs_mNCommonHit);
   fChain->SetBranchAddress("mGhostPairs.mIsBestContam", mGhostPairs_mIsBestContam, &b_mGhostPairs_mIsBestContam);
   fChain->SetBranchAddress("mGhostPairs.mDominatrack", mGhostPairs_mDominatrack, &b_mGhostPairs_mDominatrack);
   fChain->SetBranchAddress("mGhostPairs.mDominCommonHit", mGhostPairs_mDominCommonHit, &b_mGhostPairs_mDominCommonHit);
   fChain->SetBranchAddress("mGhostPairs.mAvgQuality", mGhostPairs_mAvgQuality, &b_mGhostPairs_mAvgQuality);
   fChain->SetBranchAddress("mGhostPairs.mIsValid", mGhostPairs_mIsValid, &b_mGhostPairs_mIsValid);
   fChain->SetBranchAddress("mGhostPairs.mPtMc", mGhostPairs_mPtMc, &b_mGhostPairs_mPtMc);
   fChain->SetBranchAddress("mGhostPairs.mPzMc", mGhostPairs_mPzMc, &b_mGhostPairs_mPzMc);
   fChain->SetBranchAddress("mGhostPairs.mEtaMc", mGhostPairs_mEtaMc, &b_mGhostPairs_mEtaMc);
   fChain->SetBranchAddress("mGhostPairs.mPhiMc", mGhostPairs_mPhiMc, &b_mGhostPairs_mPhiMc);
   fChain->SetBranchAddress("mGhostPairs.mNHitMc", mGhostPairs_mNHitMc, &b_mGhostPairs_mNHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNSvtHitMc", mGhostPairs_mNSvtHitMc, &b_mGhostPairs_mNSvtHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNSsdHitMc", mGhostPairs_mNSsdHitMc, &b_mGhostPairs_mNSsdHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNFtpcHitMc", mGhostPairs_mNFtpcHitMc, &b_mGhostPairs_mNFtpcHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNBemcHitMc", mGhostPairs_mNBemcHitMc, &b_mGhostPairs_mNBemcHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNBprsHitMc", mGhostPairs_mNBprsHitMc, &b_mGhostPairs_mNBprsHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNBsmdeHitMc", mGhostPairs_mNBsmdeHitMc, &b_mGhostPairs_mNBsmdeHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNBsmdpHitMc", mGhostPairs_mNBsmdpHitMc, &b_mGhostPairs_mNBsmdpHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNEemcHitMc", mGhostPairs_mNEemcHitMc, &b_mGhostPairs_mNEemcHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNEprsHitMc", mGhostPairs_mNEprsHitMc, &b_mGhostPairs_mNEprsHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNEsmduHitMc", mGhostPairs_mNEsmduHitMc, &b_mGhostPairs_mNEsmduHitMc);
   fChain->SetBranchAddress("mGhostPairs.mNEsmdvHitMc", mGhostPairs_mNEsmdvHitMc, &b_mGhostPairs_mNEsmdvHitMc);
   fChain->SetBranchAddress("mGhostPairs.mGeantId", mGhostPairs_mGeantId, &b_mGhostPairs_mGeantId);
   fChain->SetBranchAddress("mGhostPairs.mChargeMc", mGhostPairs_mChargeMc, &b_mGhostPairs_mChargeMc);
   fChain->SetBranchAddress("mGhostPairs.mStopR", mGhostPairs_mStopR, &b_mGhostPairs_mStopR);
   fChain->SetBranchAddress("mGhostPairs.mKey", mGhostPairs_mKey, &b_mGhostPairs_mKey);
   fChain->SetBranchAddress("mGhostPairs.mParentKey", mGhostPairs_mParentKey, &b_mGhostPairs_mParentKey);
   fChain->SetBranchAddress("mGhostPairs.mParentGeantId", mGhostPairs_mParentGeantId, &b_mGhostPairs_mParentGeantId);
   fChain->SetBranchAddress("mGhostPairs.mEmcEnergyMcHit[3]", mGhostPairs_mEmcEnergyMcHit, &b_mGhostPairs_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mGhostPairs.mEmcEnergyMcSum", mGhostPairs_mEmcEnergyMcSum, &b_mGhostPairs_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mGhostPairs.mEmcSoftIdHiTowerMc[3]", mGhostPairs_mEmcSoftIdHiTowerMc, &b_mGhostPairs_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mGhostPairs.mNAssocGl", mGhostPairs_mNAssocGl, &b_mGhostPairs_mNAssocGl);
   fChain->SetBranchAddress("mGhostPairs.mNAssocPr", mGhostPairs_mNAssocPr, &b_mGhostPairs_mNAssocPr);
   fChain->SetBranchAddress("mGhostPairs.mIsPrimary", mGhostPairs_mIsPrimary, &b_mGhostPairs_mIsPrimary);
   fChain->SetBranchAddress("mGhostPairs.fUniqueID", mGhostPairs_fUniqueID, &b_mGhostPairs_fUniqueID);
   fChain->SetBranchAddress("mGhostPairs.fBits", mGhostPairs_fBits, &b_mGhostPairs_fBits);
   fChain->SetBranchAddress("mGhostPairs.mIsValidGl", mGhostPairs_mIsValidGl, &b_mGhostPairs_mIsValidGl);
   fChain->SetBranchAddress("mGhostPairs.mPtPr", mGhostPairs_mPtPr, &b_mGhostPairs_mPtPr);
   fChain->SetBranchAddress("mGhostPairs.mPzPr", mGhostPairs_mPzPr, &b_mGhostPairs_mPzPr);
   fChain->SetBranchAddress("mGhostPairs.mEtaPr", mGhostPairs_mEtaPr, &b_mGhostPairs_mEtaPr);
   fChain->SetBranchAddress("mGhostPairs.mPhiPr", mGhostPairs_mPhiPr, &b_mGhostPairs_mPhiPr);
   fChain->SetBranchAddress("mGhostPairs.mDcaPr", mGhostPairs_mDcaPr, &b_mGhostPairs_mDcaPr);
   fChain->SetBranchAddress("mGhostPairs.mDcaXYPr", mGhostPairs_mDcaXYPr, &b_mGhostPairs_mDcaXYPr);
   fChain->SetBranchAddress("mGhostPairs.mDcaZPr", mGhostPairs_mDcaZPr, &b_mGhostPairs_mDcaZPr);
   fChain->SetBranchAddress("mGhostPairs.mDcaXYPrMcV", mGhostPairs_mDcaXYPrMcV, &b_mGhostPairs_mDcaXYPrMcV);
   fChain->SetBranchAddress("mGhostPairs.mDcaZPrMcV", mGhostPairs_mDcaZPrMcV, &b_mGhostPairs_mDcaZPrMcV);
   fChain->SetBranchAddress("mGhostPairs.mCurvPr", mGhostPairs_mCurvPr, &b_mGhostPairs_mCurvPr);
   fChain->SetBranchAddress("mGhostPairs.mTanLPr", mGhostPairs_mTanLPr, &b_mGhostPairs_mTanLPr);
   fChain->SetBranchAddress("mGhostPairs.mErrP[5]", mGhostPairs_mErrP, &b_mGhostPairs_mErrP);
   fChain->SetBranchAddress("mGhostPairs.mChi2Pr", mGhostPairs_mChi2Pr, &b_mGhostPairs_mChi2Pr);
   fChain->SetBranchAddress("mGhostPairs.mFlag", mGhostPairs_mFlag, &b_mGhostPairs_mFlag);
   fChain->SetBranchAddress("mGhostPairs.mDedx", mGhostPairs_mDedx, &b_mGhostPairs_mDedx);
   fChain->SetBranchAddress("mGhostPairs.mPtGl", mGhostPairs_mPtGl, &b_mGhostPairs_mPtGl);
   fChain->SetBranchAddress("mGhostPairs.mPzGl", mGhostPairs_mPzGl, &b_mGhostPairs_mPzGl);
   fChain->SetBranchAddress("mGhostPairs.mEtaGl", mGhostPairs_mEtaGl, &b_mGhostPairs_mEtaGl);
   fChain->SetBranchAddress("mGhostPairs.mPhiGl", mGhostPairs_mPhiGl, &b_mGhostPairs_mPhiGl);
   fChain->SetBranchAddress("mGhostPairs.mDcaGl", mGhostPairs_mDcaGl, &b_mGhostPairs_mDcaGl);
   fChain->SetBranchAddress("mGhostPairs.mDcaXYGl", mGhostPairs_mDcaXYGl, &b_mGhostPairs_mDcaXYGl);
   fChain->SetBranchAddress("mGhostPairs.mDcaZGl", mGhostPairs_mDcaZGl, &b_mGhostPairs_mDcaZGl);
   fChain->SetBranchAddress("mGhostPairs.mDcaXYGlMcV", mGhostPairs_mDcaXYGlMcV, &b_mGhostPairs_mDcaXYGlMcV);
   fChain->SetBranchAddress("mGhostPairs.mDcaZGlMcV", mGhostPairs_mDcaZGlMcV, &b_mGhostPairs_mDcaZGlMcV);
   fChain->SetBranchAddress("mGhostPairs.mCurvGl", mGhostPairs_mCurvGl, &b_mGhostPairs_mCurvGl);
   fChain->SetBranchAddress("mGhostPairs.mTanLGl", mGhostPairs_mTanLGl, &b_mGhostPairs_mTanLGl);
   fChain->SetBranchAddress("mGhostPairs.mErrG[5]", mGhostPairs_mErrG, &b_mGhostPairs_mErrG);
   fChain->SetBranchAddress("mGhostPairs.mSeedQA", mGhostPairs_mSeedQA, &b_mGhostPairs_mSeedQA);
   fChain->SetBranchAddress("mGhostPairs.mPidPion", mGhostPairs_mPidPion, &b_mGhostPairs_mPidPion);
   fChain->SetBranchAddress("mGhostPairs.mPidProton", mGhostPairs_mPidProton, &b_mGhostPairs_mPidProton);
   fChain->SetBranchAddress("mGhostPairs.mPidKaon", mGhostPairs_mPidKaon, &b_mGhostPairs_mPidKaon);
   fChain->SetBranchAddress("mGhostPairs.mPidElectron", mGhostPairs_mPidElectron, &b_mGhostPairs_mPidElectron);
   fChain->SetBranchAddress("mGhostPairs.mFirstZ", mGhostPairs_mFirstZ, &b_mGhostPairs_mFirstZ);
   fChain->SetBranchAddress("mGhostPairs.mLastZ", mGhostPairs_mLastZ, &b_mGhostPairs_mLastZ);
   fChain->SetBranchAddress("mGhostPairs.mFirstPadrow", mGhostPairs_mFirstPadrow, &b_mGhostPairs_mFirstPadrow);
   fChain->SetBranchAddress("mGhostPairs.mLastPadrow", mGhostPairs_mLastPadrow, &b_mGhostPairs_mLastPadrow);
   fChain->SetBranchAddress("mGhostPairs.mFirstFitPadrow", mGhostPairs_mFirstFitPadrow, &b_mGhostPairs_mFirstFitPadrow);
   fChain->SetBranchAddress("mGhostPairs.mLastFitPadrow", mGhostPairs_mLastFitPadrow, &b_mGhostPairs_mLastFitPadrow);
   fChain->SetBranchAddress("mGhostPairs.mFirstSector", mGhostPairs_mFirstSector, &b_mGhostPairs_mFirstSector);
   fChain->SetBranchAddress("mGhostPairs.mLastSector", mGhostPairs_mLastSector, &b_mGhostPairs_mLastSector);
   fChain->SetBranchAddress("mGhostPairs.mRecoKey", mGhostPairs_mRecoKey, &b_mGhostPairs_mRecoKey);
   fChain->SetBranchAddress("mGhostPairs.mFitPts", mGhostPairs_mFitPts, &b_mGhostPairs_mFitPts);
   fChain->SetBranchAddress("mGhostPairs.mFitSvt", mGhostPairs_mFitSvt, &b_mGhostPairs_mFitSvt);
   fChain->SetBranchAddress("mGhostPairs.mFitSsd", mGhostPairs_mFitSsd, &b_mGhostPairs_mFitSsd);
   fChain->SetBranchAddress("mGhostPairs.mFitFtpc", mGhostPairs_mFitFtpc, &b_mGhostPairs_mFitFtpc);
   fChain->SetBranchAddress("mGhostPairs.mDedxPts", mGhostPairs_mDedxPts, &b_mGhostPairs_mDedxPts);
   fChain->SetBranchAddress("mGhostPairs.mAllPts", mGhostPairs_mAllPts, &b_mGhostPairs_mAllPts);
   fChain->SetBranchAddress("mGhostPairs.mCharge", mGhostPairs_mCharge, &b_mGhostPairs_mCharge);
   fChain->SetBranchAddress("mGhostPairs.mNAssocMc", mGhostPairs_mNAssocMc, &b_mGhostPairs_mNAssocMc);
   fChain->SetBranchAddress("mGhostPairs.mNPossible", mGhostPairs_mNPossible, &b_mGhostPairs_mNPossible);
   fChain->SetBranchAddress("mGhostPairs.mIsValidPr", mGhostPairs_mIsValidPr, &b_mGhostPairs_mIsValidPr);
   fChain->SetBranchAddress("mGhostPairs.mEmcTowerAdc[3]", mGhostPairs_mEmcTowerAdc, &b_mGhostPairs_mEmcTowerAdc);
   fChain->SetBranchAddress("mGhostPairs.mEmcEnergyRcHit[3]", mGhostPairs_mEmcEnergyRcHit, &b_mGhostPairs_mEmcEnergyRcHit);
   fChain->SetBranchAddress("mGhostPairs.mEmcSoftIdHiTowerRc[3]", mGhostPairs_mEmcSoftIdHiTowerRc, &b_mGhostPairs_mEmcSoftIdHiTowerRc);
   fChain->SetBranchAddress("mContamPairs_", &mContamPairs_, &b_mContamPairs_);
   fChain->SetBranchAddress("mContamPairs.mParentGeantId", mContamPairs_mParentGeantId, &b_mContamPairs_mParentGeantId);
   fChain->SetBranchAddress("mContamPairs.mGeantProcess", mContamPairs_mGeantProcess, &b_mContamPairs_mGeantProcess);
   fChain->SetBranchAddress("mContamPairs.mPtMcParent", mContamPairs_mPtMcParent, &b_mContamPairs_mPtMcParent);
   fChain->SetBranchAddress("mContamPairs.mEtaMcParent", mContamPairs_mEtaMcParent, &b_mContamPairs_mEtaMcParent);
   fChain->SetBranchAddress("mContamPairs.mParentParentGeantId", mContamPairs_mParentParentGeantId, &b_mContamPairs_mParentParentGeantId);
   fChain->SetBranchAddress("mContamPairs.mPtMcParentParent", mContamPairs_mPtMcParentParent, &b_mContamPairs_mPtMcParentParent);
   fChain->SetBranchAddress("mContamPairs.mStartX", mContamPairs_mStartX, &b_mContamPairs_mStartX);
   fChain->SetBranchAddress("mContamPairs.mStartY", mContamPairs_mStartY, &b_mContamPairs_mStartY);
   fChain->SetBranchAddress("mContamPairs.mStartZ", mContamPairs_mStartZ, &b_mContamPairs_mStartZ);
   fChain->SetBranchAddress("mContamPairs.mNCommonHit", mContamPairs_mNCommonHit, &b_mContamPairs_mNCommonHit);
   fChain->SetBranchAddress("mContamPairs.mIsBestContam", mContamPairs_mIsBestContam, &b_mContamPairs_mIsBestContam);
   fChain->SetBranchAddress("mContamPairs.mDominatrack", mContamPairs_mDominatrack, &b_mContamPairs_mDominatrack);
   fChain->SetBranchAddress("mContamPairs.mDominCommonHit", mContamPairs_mDominCommonHit, &b_mContamPairs_mDominCommonHit);
   fChain->SetBranchAddress("mContamPairs.mAvgQuality", mContamPairs_mAvgQuality, &b_mContamPairs_mAvgQuality);
   fChain->SetBranchAddress("mContamPairs.mIsValid", mContamPairs_mIsValid, &b_mContamPairs_mIsValid);
   fChain->SetBranchAddress("mContamPairs.mPtMc", mContamPairs_mPtMc, &b_mContamPairs_mPtMc);
   fChain->SetBranchAddress("mContamPairs.mPzMc", mContamPairs_mPzMc, &b_mContamPairs_mPzMc);
   fChain->SetBranchAddress("mContamPairs.mEtaMc", mContamPairs_mEtaMc, &b_mContamPairs_mEtaMc);
   fChain->SetBranchAddress("mContamPairs.mPhiMc", mContamPairs_mPhiMc, &b_mContamPairs_mPhiMc);
   fChain->SetBranchAddress("mContamPairs.mNHitMc", mContamPairs_mNHitMc, &b_mContamPairs_mNHitMc);
   fChain->SetBranchAddress("mContamPairs.mNSvtHitMc", mContamPairs_mNSvtHitMc, &b_mContamPairs_mNSvtHitMc);
   fChain->SetBranchAddress("mContamPairs.mNSsdHitMc", mContamPairs_mNSsdHitMc, &b_mContamPairs_mNSsdHitMc);
   fChain->SetBranchAddress("mContamPairs.mNFtpcHitMc", mContamPairs_mNFtpcHitMc, &b_mContamPairs_mNFtpcHitMc);
   fChain->SetBranchAddress("mContamPairs.mNBemcHitMc", mContamPairs_mNBemcHitMc, &b_mContamPairs_mNBemcHitMc);
   fChain->SetBranchAddress("mContamPairs.mNBprsHitMc", mContamPairs_mNBprsHitMc, &b_mContamPairs_mNBprsHitMc);
   fChain->SetBranchAddress("mContamPairs.mNBsmdeHitMc", mContamPairs_mNBsmdeHitMc, &b_mContamPairs_mNBsmdeHitMc);
   fChain->SetBranchAddress("mContamPairs.mNBsmdpHitMc", mContamPairs_mNBsmdpHitMc, &b_mContamPairs_mNBsmdpHitMc);
   fChain->SetBranchAddress("mContamPairs.mNEemcHitMc", mContamPairs_mNEemcHitMc, &b_mContamPairs_mNEemcHitMc);
   fChain->SetBranchAddress("mContamPairs.mNEprsHitMc", mContamPairs_mNEprsHitMc, &b_mContamPairs_mNEprsHitMc);
   fChain->SetBranchAddress("mContamPairs.mNEsmduHitMc", mContamPairs_mNEsmduHitMc, &b_mContamPairs_mNEsmduHitMc);
   fChain->SetBranchAddress("mContamPairs.mNEsmdvHitMc", mContamPairs_mNEsmdvHitMc, &b_mContamPairs_mNEsmdvHitMc);
   fChain->SetBranchAddress("mContamPairs.mGeantId", mContamPairs_mGeantId, &b_mContamPairs_mGeantId);
   fChain->SetBranchAddress("mContamPairs.mChargeMc", mContamPairs_mChargeMc, &b_mContamPairs_mChargeMc);
   fChain->SetBranchAddress("mContamPairs.mStopR", mContamPairs_mStopR, &b_mContamPairs_mStopR);
   fChain->SetBranchAddress("mContamPairs.mKey", mContamPairs_mKey, &b_mContamPairs_mKey);
   fChain->SetBranchAddress("mContamPairs.mParentKey", mContamPairs_mParentKey, &b_mContamPairs_mParentKey);
//    fChain->SetBranchAddress("mContamPairs.mParentGeantId", mContamPairs_mParentGeantId, &b_mContamPairs_mParentGeantId);
   fChain->SetBranchAddress("mContamPairs.mEmcEnergyMcHit[3]", mContamPairs_mEmcEnergyMcHit, &b_mContamPairs_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mContamPairs.mEmcEnergyMcSum", mContamPairs_mEmcEnergyMcSum, &b_mContamPairs_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mContamPairs.mEmcSoftIdHiTowerMc[3]", mContamPairs_mEmcSoftIdHiTowerMc, &b_mContamPairs_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mContamPairs.mNAssocGl", mContamPairs_mNAssocGl, &b_mContamPairs_mNAssocGl);
   fChain->SetBranchAddress("mContamPairs.mNAssocPr", mContamPairs_mNAssocPr, &b_mContamPairs_mNAssocPr);
   fChain->SetBranchAddress("mContamPairs.mIsPrimary", mContamPairs_mIsPrimary, &b_mContamPairs_mIsPrimary);
   fChain->SetBranchAddress("mContamPairs.fUniqueID", mContamPairs_fUniqueID, &b_mContamPairs_fUniqueID);
   fChain->SetBranchAddress("mContamPairs.fBits", mContamPairs_fBits, &b_mContamPairs_fBits);
   fChain->SetBranchAddress("mContamPairs.mIsValidGl", mContamPairs_mIsValidGl, &b_mContamPairs_mIsValidGl);
   fChain->SetBranchAddress("mContamPairs.mPtPr", mContamPairs_mPtPr, &b_mContamPairs_mPtPr);
   fChain->SetBranchAddress("mContamPairs.mPzPr", mContamPairs_mPzPr, &b_mContamPairs_mPzPr);
   fChain->SetBranchAddress("mContamPairs.mEtaPr", mContamPairs_mEtaPr, &b_mContamPairs_mEtaPr);
   fChain->SetBranchAddress("mContamPairs.mPhiPr", mContamPairs_mPhiPr, &b_mContamPairs_mPhiPr);
   fChain->SetBranchAddress("mContamPairs.mDcaPr", mContamPairs_mDcaPr, &b_mContamPairs_mDcaPr);
   fChain->SetBranchAddress("mContamPairs.mDcaXYPr", mContamPairs_mDcaXYPr, &b_mContamPairs_mDcaXYPr);
   fChain->SetBranchAddress("mContamPairs.mDcaZPr", mContamPairs_mDcaZPr, &b_mContamPairs_mDcaZPr);
   fChain->SetBranchAddress("mContamPairs.mDcaXYPrMcV", mContamPairs_mDcaXYPrMcV, &b_mContamPairs_mDcaXYPrMcV);
   fChain->SetBranchAddress("mContamPairs.mDcaZPrMcV", mContamPairs_mDcaZPrMcV, &b_mContamPairs_mDcaZPrMcV);
   fChain->SetBranchAddress("mContamPairs.mCurvPr", mContamPairs_mCurvPr, &b_mContamPairs_mCurvPr);
   fChain->SetBranchAddress("mContamPairs.mTanLPr", mContamPairs_mTanLPr, &b_mContamPairs_mTanLPr);
   fChain->SetBranchAddress("mContamPairs.mErrP[5]", mContamPairs_mErrP, &b_mContamPairs_mErrP);
   fChain->SetBranchAddress("mContamPairs.mChi2Pr", mContamPairs_mChi2Pr, &b_mContamPairs_mChi2Pr);
   fChain->SetBranchAddress("mContamPairs.mFlag", mContamPairs_mFlag, &b_mContamPairs_mFlag);
   fChain->SetBranchAddress("mContamPairs.mDedx", mContamPairs_mDedx, &b_mContamPairs_mDedx);
   fChain->SetBranchAddress("mContamPairs.mPtGl", mContamPairs_mPtGl, &b_mContamPairs_mPtGl);
   fChain->SetBranchAddress("mContamPairs.mPzGl", mContamPairs_mPzGl, &b_mContamPairs_mPzGl);
   fChain->SetBranchAddress("mContamPairs.mEtaGl", mContamPairs_mEtaGl, &b_mContamPairs_mEtaGl);
   fChain->SetBranchAddress("mContamPairs.mPhiGl", mContamPairs_mPhiGl, &b_mContamPairs_mPhiGl);
   fChain->SetBranchAddress("mContamPairs.mDcaGl", mContamPairs_mDcaGl, &b_mContamPairs_mDcaGl);
   fChain->SetBranchAddress("mContamPairs.mDcaXYGl", mContamPairs_mDcaXYGl, &b_mContamPairs_mDcaXYGl);
   fChain->SetBranchAddress("mContamPairs.mDcaZGl", mContamPairs_mDcaZGl, &b_mContamPairs_mDcaZGl);
   fChain->SetBranchAddress("mContamPairs.mDcaXYGlMcV", mContamPairs_mDcaXYGlMcV, &b_mContamPairs_mDcaXYGlMcV);
   fChain->SetBranchAddress("mContamPairs.mDcaZGlMcV", mContamPairs_mDcaZGlMcV, &b_mContamPairs_mDcaZGlMcV);
   fChain->SetBranchAddress("mContamPairs.mCurvGl", mContamPairs_mCurvGl, &b_mContamPairs_mCurvGl);
   fChain->SetBranchAddress("mContamPairs.mTanLGl", mContamPairs_mTanLGl, &b_mContamPairs_mTanLGl);
   fChain->SetBranchAddress("mContamPairs.mErrG[5]", mContamPairs_mErrG, &b_mContamPairs_mErrG);
   fChain->SetBranchAddress("mContamPairs.mSeedQA", mContamPairs_mSeedQA, &b_mContamPairs_mSeedQA);
   fChain->SetBranchAddress("mContamPairs.mPidPion", mContamPairs_mPidPion, &b_mContamPairs_mPidPion);
   fChain->SetBranchAddress("mContamPairs.mPidProton", mContamPairs_mPidProton, &b_mContamPairs_mPidProton);
   fChain->SetBranchAddress("mContamPairs.mPidKaon", mContamPairs_mPidKaon, &b_mContamPairs_mPidKaon);
   fChain->SetBranchAddress("mContamPairs.mPidElectron", mContamPairs_mPidElectron, &b_mContamPairs_mPidElectron);
   fChain->SetBranchAddress("mContamPairs.mFirstZ", mContamPairs_mFirstZ, &b_mContamPairs_mFirstZ);
   fChain->SetBranchAddress("mContamPairs.mLastZ", mContamPairs_mLastZ, &b_mContamPairs_mLastZ);
   fChain->SetBranchAddress("mContamPairs.mFirstPadrow", mContamPairs_mFirstPadrow, &b_mContamPairs_mFirstPadrow);
   fChain->SetBranchAddress("mContamPairs.mLastPadrow", mContamPairs_mLastPadrow, &b_mContamPairs_mLastPadrow);
   fChain->SetBranchAddress("mContamPairs.mFirstFitPadrow", mContamPairs_mFirstFitPadrow, &b_mContamPairs_mFirstFitPadrow);
   fChain->SetBranchAddress("mContamPairs.mLastFitPadrow", mContamPairs_mLastFitPadrow, &b_mContamPairs_mLastFitPadrow);
   fChain->SetBranchAddress("mContamPairs.mFirstSector", mContamPairs_mFirstSector, &b_mContamPairs_mFirstSector);
   fChain->SetBranchAddress("mContamPairs.mLastSector", mContamPairs_mLastSector, &b_mContamPairs_mLastSector);
   fChain->SetBranchAddress("mContamPairs.mRecoKey", mContamPairs_mRecoKey, &b_mContamPairs_mRecoKey);
   fChain->SetBranchAddress("mContamPairs.mFitPts", mContamPairs_mFitPts, &b_mContamPairs_mFitPts);
   fChain->SetBranchAddress("mContamPairs.mFitSvt", mContamPairs_mFitSvt, &b_mContamPairs_mFitSvt);
   fChain->SetBranchAddress("mContamPairs.mFitSsd", mContamPairs_mFitSsd, &b_mContamPairs_mFitSsd);
   fChain->SetBranchAddress("mContamPairs.mFitFtpc", mContamPairs_mFitFtpc, &b_mContamPairs_mFitFtpc);
   fChain->SetBranchAddress("mContamPairs.mDedxPts", mContamPairs_mDedxPts, &b_mContamPairs_mDedxPts);
   fChain->SetBranchAddress("mContamPairs.mAllPts", mContamPairs_mAllPts, &b_mContamPairs_mAllPts);
   fChain->SetBranchAddress("mContamPairs.mCharge", mContamPairs_mCharge, &b_mContamPairs_mCharge);
   fChain->SetBranchAddress("mContamPairs.mNAssocMc", mContamPairs_mNAssocMc, &b_mContamPairs_mNAssocMc);
   fChain->SetBranchAddress("mContamPairs.mNPossible", mContamPairs_mNPossible, &b_mContamPairs_mNPossible);
   fChain->SetBranchAddress("mContamPairs.mIsValidPr", mContamPairs_mIsValidPr, &b_mContamPairs_mIsValidPr);
   fChain->SetBranchAddress("mContamPairs.mEmcTowerAdc[3]", mContamPairs_mEmcTowerAdc, &b_mContamPairs_mEmcTowerAdc);
   fChain->SetBranchAddress("mContamPairs.mEmcEnergyRcHit[3]", mContamPairs_mEmcEnergyRcHit, &b_mContamPairs_mEmcEnergyRcHit);
   fChain->SetBranchAddress("mContamPairs.mEmcSoftIdHiTowerRc[3]", mContamPairs_mEmcSoftIdHiTowerRc, &b_mContamPairs_mEmcSoftIdHiTowerRc);
   fChain->SetBranchAddress("mMatGlobPairs_", &mMatGlobPairs_, &b_mMatGlobPairs_);
   fChain->SetBranchAddress("mMatGlobPairs.mNCommonHit", mMatGlobPairs_mNCommonHit, &b_mMatGlobPairs_mNCommonHit);
   fChain->SetBranchAddress("mMatGlobPairs.mIsBestContam", mMatGlobPairs_mIsBestContam, &b_mMatGlobPairs_mIsBestContam);
   fChain->SetBranchAddress("mMatGlobPairs.mDominatrack", mMatGlobPairs_mDominatrack, &b_mMatGlobPairs_mDominatrack);
   fChain->SetBranchAddress("mMatGlobPairs.mDominCommonHit", mMatGlobPairs_mDominCommonHit, &b_mMatGlobPairs_mDominCommonHit);
   fChain->SetBranchAddress("mMatGlobPairs.mAvgQuality", mMatGlobPairs_mAvgQuality, &b_mMatGlobPairs_mAvgQuality);
   fChain->SetBranchAddress("mMatGlobPairs.mIsValid", mMatGlobPairs_mIsValid, &b_mMatGlobPairs_mIsValid);
   fChain->SetBranchAddress("mMatGlobPairs.mPtMc", mMatGlobPairs_mPtMc, &b_mMatGlobPairs_mPtMc);
   fChain->SetBranchAddress("mMatGlobPairs.mPzMc", mMatGlobPairs_mPzMc, &b_mMatGlobPairs_mPzMc);
   fChain->SetBranchAddress("mMatGlobPairs.mEtaMc", mMatGlobPairs_mEtaMc, &b_mMatGlobPairs_mEtaMc);
   fChain->SetBranchAddress("mMatGlobPairs.mPhiMc", mMatGlobPairs_mPhiMc, &b_mMatGlobPairs_mPhiMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNHitMc", mMatGlobPairs_mNHitMc, &b_mMatGlobPairs_mNHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNSvtHitMc", mMatGlobPairs_mNSvtHitMc, &b_mMatGlobPairs_mNSvtHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNSsdHitMc", mMatGlobPairs_mNSsdHitMc, &b_mMatGlobPairs_mNSsdHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNFtpcHitMc", mMatGlobPairs_mNFtpcHitMc, &b_mMatGlobPairs_mNFtpcHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNBemcHitMc", mMatGlobPairs_mNBemcHitMc, &b_mMatGlobPairs_mNBemcHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNBprsHitMc", mMatGlobPairs_mNBprsHitMc, &b_mMatGlobPairs_mNBprsHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNBsmdeHitMc", mMatGlobPairs_mNBsmdeHitMc, &b_mMatGlobPairs_mNBsmdeHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNBsmdpHitMc", mMatGlobPairs_mNBsmdpHitMc, &b_mMatGlobPairs_mNBsmdpHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNEemcHitMc", mMatGlobPairs_mNEemcHitMc, &b_mMatGlobPairs_mNEemcHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNEprsHitMc", mMatGlobPairs_mNEprsHitMc, &b_mMatGlobPairs_mNEprsHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNEsmduHitMc", mMatGlobPairs_mNEsmduHitMc, &b_mMatGlobPairs_mNEsmduHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNEsmdvHitMc", mMatGlobPairs_mNEsmdvHitMc, &b_mMatGlobPairs_mNEsmdvHitMc);
   fChain->SetBranchAddress("mMatGlobPairs.mGeantId", mMatGlobPairs_mGeantId, &b_mMatGlobPairs_mGeantId);
   fChain->SetBranchAddress("mMatGlobPairs.mChargeMc", mMatGlobPairs_mChargeMc, &b_mMatGlobPairs_mChargeMc);
   fChain->SetBranchAddress("mMatGlobPairs.mStopR", mMatGlobPairs_mStopR, &b_mMatGlobPairs_mStopR);
   fChain->SetBranchAddress("mMatGlobPairs.mKey", mMatGlobPairs_mKey, &b_mMatGlobPairs_mKey);
   fChain->SetBranchAddress("mMatGlobPairs.mParentKey", mMatGlobPairs_mParentKey, &b_mMatGlobPairs_mParentKey);
   fChain->SetBranchAddress("mMatGlobPairs.mParentGeantId", mMatGlobPairs_mParentGeantId, &b_mMatGlobPairs_mParentGeantId);
   fChain->SetBranchAddress("mMatGlobPairs.mEmcEnergyMcHit[3]", mMatGlobPairs_mEmcEnergyMcHit, &b_mMatGlobPairs_mEmcEnergyMcHit);
   fChain->SetBranchAddress("mMatGlobPairs.mEmcEnergyMcSum", mMatGlobPairs_mEmcEnergyMcSum, &b_mMatGlobPairs_mEmcEnergyMcSum);
   fChain->SetBranchAddress("mMatGlobPairs.mEmcSoftIdHiTowerMc[3]", mMatGlobPairs_mEmcSoftIdHiTowerMc, &b_mMatGlobPairs_mEmcSoftIdHiTowerMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNAssocGl", mMatGlobPairs_mNAssocGl, &b_mMatGlobPairs_mNAssocGl);
   fChain->SetBranchAddress("mMatGlobPairs.mNAssocPr", mMatGlobPairs_mNAssocPr, &b_mMatGlobPairs_mNAssocPr);
   fChain->SetBranchAddress("mMatGlobPairs.mIsPrimary", mMatGlobPairs_mIsPrimary, &b_mMatGlobPairs_mIsPrimary);
   fChain->SetBranchAddress("mMatGlobPairs.fUniqueID", mMatGlobPairs_fUniqueID, &b_mMatGlobPairs_fUniqueID);
   fChain->SetBranchAddress("mMatGlobPairs.fBits", mMatGlobPairs_fBits, &b_mMatGlobPairs_fBits);
   fChain->SetBranchAddress("mMatGlobPairs.mIsValidGl", mMatGlobPairs_mIsValidGl, &b_mMatGlobPairs_mIsValidGl);
   fChain->SetBranchAddress("mMatGlobPairs.mPtPr", mMatGlobPairs_mPtPr, &b_mMatGlobPairs_mPtPr);
   fChain->SetBranchAddress("mMatGlobPairs.mPzPr", mMatGlobPairs_mPzPr, &b_mMatGlobPairs_mPzPr);
   fChain->SetBranchAddress("mMatGlobPairs.mEtaPr", mMatGlobPairs_mEtaPr, &b_mMatGlobPairs_mEtaPr);
   fChain->SetBranchAddress("mMatGlobPairs.mPhiPr", mMatGlobPairs_mPhiPr, &b_mMatGlobPairs_mPhiPr);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaPr", mMatGlobPairs_mDcaPr, &b_mMatGlobPairs_mDcaPr);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaXYPr", mMatGlobPairs_mDcaXYPr, &b_mMatGlobPairs_mDcaXYPr);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaZPr", mMatGlobPairs_mDcaZPr, &b_mMatGlobPairs_mDcaZPr);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaXYPrMcV", mMatGlobPairs_mDcaXYPrMcV, &b_mMatGlobPairs_mDcaXYPrMcV);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaZPrMcV", mMatGlobPairs_mDcaZPrMcV, &b_mMatGlobPairs_mDcaZPrMcV);
   fChain->SetBranchAddress("mMatGlobPairs.mCurvPr", mMatGlobPairs_mCurvPr, &b_mMatGlobPairs_mCurvPr);
   fChain->SetBranchAddress("mMatGlobPairs.mTanLPr", mMatGlobPairs_mTanLPr, &b_mMatGlobPairs_mTanLPr);
   fChain->SetBranchAddress("mMatGlobPairs.mErrP[5]", mMatGlobPairs_mErrP, &b_mMatGlobPairs_mErrP);
   fChain->SetBranchAddress("mMatGlobPairs.mChi2Pr", mMatGlobPairs_mChi2Pr, &b_mMatGlobPairs_mChi2Pr);
   fChain->SetBranchAddress("mMatGlobPairs.mFlag", mMatGlobPairs_mFlag, &b_mMatGlobPairs_mFlag);
   fChain->SetBranchAddress("mMatGlobPairs.mDedx", mMatGlobPairs_mDedx, &b_mMatGlobPairs_mDedx);
   fChain->SetBranchAddress("mMatGlobPairs.mPtGl", mMatGlobPairs_mPtGl, &b_mMatGlobPairs_mPtGl);
   fChain->SetBranchAddress("mMatGlobPairs.mPzGl", mMatGlobPairs_mPzGl, &b_mMatGlobPairs_mPzGl);
   fChain->SetBranchAddress("mMatGlobPairs.mEtaGl", mMatGlobPairs_mEtaGl, &b_mMatGlobPairs_mEtaGl);
   fChain->SetBranchAddress("mMatGlobPairs.mPhiGl", mMatGlobPairs_mPhiGl, &b_mMatGlobPairs_mPhiGl);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaGl", mMatGlobPairs_mDcaGl, &b_mMatGlobPairs_mDcaGl);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaXYGl", mMatGlobPairs_mDcaXYGl, &b_mMatGlobPairs_mDcaXYGl);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaZGl", mMatGlobPairs_mDcaZGl, &b_mMatGlobPairs_mDcaZGl);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaXYGlMcV", mMatGlobPairs_mDcaXYGlMcV, &b_mMatGlobPairs_mDcaXYGlMcV);
   fChain->SetBranchAddress("mMatGlobPairs.mDcaZGlMcV", mMatGlobPairs_mDcaZGlMcV, &b_mMatGlobPairs_mDcaZGlMcV);
   fChain->SetBranchAddress("mMatGlobPairs.mCurvGl", mMatGlobPairs_mCurvGl, &b_mMatGlobPairs_mCurvGl);
   fChain->SetBranchAddress("mMatGlobPairs.mTanLGl", mMatGlobPairs_mTanLGl, &b_mMatGlobPairs_mTanLGl);
   fChain->SetBranchAddress("mMatGlobPairs.mErrG[5]", mMatGlobPairs_mErrG, &b_mMatGlobPairs_mErrG);
   fChain->SetBranchAddress("mMatGlobPairs.mSeedQA", mMatGlobPairs_mSeedQA, &b_mMatGlobPairs_mSeedQA);
   fChain->SetBranchAddress("mMatGlobPairs.mPidPion", mMatGlobPairs_mPidPion, &b_mMatGlobPairs_mPidPion);
   fChain->SetBranchAddress("mMatGlobPairs.mPidProton", mMatGlobPairs_mPidProton, &b_mMatGlobPairs_mPidProton);
   fChain->SetBranchAddress("mMatGlobPairs.mPidKaon", mMatGlobPairs_mPidKaon, &b_mMatGlobPairs_mPidKaon);
   fChain->SetBranchAddress("mMatGlobPairs.mPidElectron", mMatGlobPairs_mPidElectron, &b_mMatGlobPairs_mPidElectron);
   fChain->SetBranchAddress("mMatGlobPairs.mFirstZ", mMatGlobPairs_mFirstZ, &b_mMatGlobPairs_mFirstZ);
   fChain->SetBranchAddress("mMatGlobPairs.mLastZ", mMatGlobPairs_mLastZ, &b_mMatGlobPairs_mLastZ);
   fChain->SetBranchAddress("mMatGlobPairs.mFirstPadrow", mMatGlobPairs_mFirstPadrow, &b_mMatGlobPairs_mFirstPadrow);
   fChain->SetBranchAddress("mMatGlobPairs.mLastPadrow", mMatGlobPairs_mLastPadrow, &b_mMatGlobPairs_mLastPadrow);
   fChain->SetBranchAddress("mMatGlobPairs.mFirstFitPadrow", mMatGlobPairs_mFirstFitPadrow, &b_mMatGlobPairs_mFirstFitPadrow);
   fChain->SetBranchAddress("mMatGlobPairs.mLastFitPadrow", mMatGlobPairs_mLastFitPadrow, &b_mMatGlobPairs_mLastFitPadrow);
   fChain->SetBranchAddress("mMatGlobPairs.mFirstSector", mMatGlobPairs_mFirstSector, &b_mMatGlobPairs_mFirstSector);
   fChain->SetBranchAddress("mMatGlobPairs.mLastSector", mMatGlobPairs_mLastSector, &b_mMatGlobPairs_mLastSector);
   fChain->SetBranchAddress("mMatGlobPairs.mRecoKey", mMatGlobPairs_mRecoKey, &b_mMatGlobPairs_mRecoKey);
   fChain->SetBranchAddress("mMatGlobPairs.mFitPts", mMatGlobPairs_mFitPts, &b_mMatGlobPairs_mFitPts);
   fChain->SetBranchAddress("mMatGlobPairs.mFitSvt", mMatGlobPairs_mFitSvt, &b_mMatGlobPairs_mFitSvt);
   fChain->SetBranchAddress("mMatGlobPairs.mFitSsd", mMatGlobPairs_mFitSsd, &b_mMatGlobPairs_mFitSsd);
   fChain->SetBranchAddress("mMatGlobPairs.mFitFtpc", mMatGlobPairs_mFitFtpc, &b_mMatGlobPairs_mFitFtpc);
   fChain->SetBranchAddress("mMatGlobPairs.mDedxPts", mMatGlobPairs_mDedxPts, &b_mMatGlobPairs_mDedxPts);
   fChain->SetBranchAddress("mMatGlobPairs.mAllPts", mMatGlobPairs_mAllPts, &b_mMatGlobPairs_mAllPts);
   fChain->SetBranchAddress("mMatGlobPairs.mCharge", mMatGlobPairs_mCharge, &b_mMatGlobPairs_mCharge);
   fChain->SetBranchAddress("mMatGlobPairs.mNAssocMc", mMatGlobPairs_mNAssocMc, &b_mMatGlobPairs_mNAssocMc);
   fChain->SetBranchAddress("mMatGlobPairs.mNPossible", mMatGlobPairs_mNPossible, &b_mMatGlobPairs_mNPossible);
   fChain->SetBranchAddress("mMatGlobPairs.mIsValidPr", mMatGlobPairs_mIsValidPr, &b_mMatGlobPairs_mIsValidPr);
   fChain->SetBranchAddress("mMatGlobPairs.mEmcTowerAdc[3]", mMatGlobPairs_mEmcTowerAdc, &b_mMatGlobPairs_mEmcTowerAdc);
   fChain->SetBranchAddress("mMatGlobPairs.mEmcEnergyRcHit[3]", mMatGlobPairs_mEmcEnergyRcHit, &b_mMatGlobPairs_mEmcEnergyRcHit);
   fChain->SetBranchAddress("mMatGlobPairs.mEmcSoftIdHiTowerRc[3]", mMatGlobPairs_mEmcSoftIdHiTowerRc, &b_mMatGlobPairs_mEmcSoftIdHiTowerRc);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_StMiniMcEvent_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_StMiniMcEvent_fBits);
#endif
   Notify();
}

Bool_t StMiniMcTree::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void StMiniMcTree::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t StMiniMcTree::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void StMiniMcTree::Loop(Long64_t nentries, TString outdir, TString particle, Int_t energy) {
  if (fChain == 0) return;

  //RunIdEventsDb* runevdb = new RunIdEventsDb();
  //runevdb->ReadDb("./StRoot/RunIdEventsDb/",energy,"");
  //Int_t NrRunIds = runevdb->GetTotalNrRunIds(); 

  TFile* out = new TFile(Form("%s/Single%sNT_Embed_%dGeV.root",outdir.Data(),particle.Data(),energy),"recreate");
  // histograms
  TH1F* hEvents = new TH1F("hEvents",";;# events",5,0,5);
  TH1F* hRefMult = new TH1F("hRefMult",";RefMult;counts",1000,0,1000);
  TH1F* hVzMc = new TH1F("hVzMc",";VzMc (cm);",200,-50,50);
  TH2F* hVyMcVxMc = new TH2F("hVyMcVxMc",";VxMc (cm);VyMc (cm)",100,-1,1,100,-1,1);
  TH2F* hNrMc_RefMult = new TH2F("hNrMc_RefMult",";RefMult;# of MC particles",200,0,1000,30,0,30);
  TH1F* hViRc_Mc[3];
  for ( Int_t i = 0; i < 3; ++i ) hViRc_Mc[i] = new TH1F(Form("hViRc_Mc[%d]",i),";ViRc - ViMc (cm)",100,-4,4);
  // mctrack ntuple
  static TNtuple *McTrack_NT;
  static Float_t McTrack_NTDataArray[27];
  McTrack_NT = new TNtuple("McTrack_NT",Form("McTrack_NT Ntuple single %s",particle.Data()),"EventId:RunId:RefMult:RefMultCorrected:CentralMult:CentralityWeight:Centrality16:Centrality:NMcNch:NMcHminus:VertexX:VertexY:VertexZ:McVertexX:McVertexY:McVertexZ:NMcTrack:NMatchedPair:PtMc:PzMc:EtaMc:PhiMc:NHitMc:GeantId:Key:ParentKey:ParentGeantId");
  McTrack_NT ->SetAutoSave( 5000000 );
  // matchedpairs ntuple
  static TNtuple *MatchedPairs_NT;
  static Float_t MatchedPairs_NTDataArray[45];
  MatchedPairs_NT = new TNtuple("MatchedPairs_NT",Form("McTrack_NT Ntuple single %s",particle.Data()),"EventId:RunId:RefMult:RefMultCorrected:CentralMult:CentralityWeight:Centrality16:Centrality:NMcNch:NMcHminus:VertexX:VertexY:VertexZ:McVertexX:McVertexY:McVertexZ:NMcTrack:NMatchedPair:PtMc:PzMc:EtaMc:PhiMc:NHitMc:GeantId:Key:ParentKey:ParentGeantId:NCommonHit:PtPr:PtGl:EtaPr:PhiPr:DcaGl:DcaXYGl:DcaZGl:DcaXYGlMcV:DcaZGlMcV:Flag:Dedx:RecoKey:FitPts:DedxPts:AllPts:NPossible:mErrP");
  MatchedPairs_NT ->SetAutoSave( 5000000 );

  if ( nentries < 0 ) nentries = fChain->GetEntries();
  cout << "nentries = " << nentries << endl;
  Int_t ntrks = 0;
  Int_t ntrks2 = 0;
  setbuf(stdout, NULL);

  Float_t EventId, RunId, RefMult, RefMultCorrected, CentralMult, CentralityWeight, Centrality16, Centrality, NMcNch, NMcHminus, VertexX, VertexY, VertexZ;
  Float_t McVertexX, McVertexY, McVertexZ, NMcTrack, NMatchedPair; 
  Float_t PtMc, PzMc, EtaMc, PhiMc, NHitMc, GeantId, Key, ParentKey, ParentGeantId;
  Float_t NCommonHit, PtPr, PtGl, EtaPr, PhiPr, DcaGl, DcaXYGl, DcaZGl, DcaXYGlMcV, DcaZGlMcV, Flag, Dedx, ErrP; 
  Float_t RecoKey, FitPts, DedxPts, AllPts, NPossible;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry); 

    hEvents->Fill(1);

    VertexX      = (Float_t)mVertexX; 
    VertexY      = (Float_t)mVertexY;
    VertexZ      = (Float_t)mVertexZ;

    if ( VertexZ > -100 && VertexZ < 100 && VertexX*VertexX+VertexY*VertexY < 4 ) {

      hEvents->Fill(2);

      RefMult   = (Float_t)mNUncorrectedPrimaries;

      refmultCorrUtil->init(mRunId);
      refmultCorrUtil->initEvent(RefMult, VertexZ);

      if ( ! refmultCorrUtil->isBadRun(mRunId) ) {

        hEvents->Fill(3);

        EventId      = (Float_t)mEventId;
        RunId        = (Float_t)mRunId;
        CentralMult  = (Float_t)mCentralMult;
        RefMultCorrected = (Float_t) refmultCorrUtil->getRefMultCorr();
        CentralityWeight   = (Float_t)refmultCorrUtil->getWeight();
        Centrality16   = (Float_t)refmultCorrUtil->getCentralityBin16();
        Centrality   = (Float_t)mCentrality;
        NMcNch       = (Float_t)mNMcNch;
        NMcHminus    = (Float_t)mNMcHminus;
        McVertexX    = (Float_t)mMcVertexX;
        McVertexY    = (Float_t)mMcVertexY;
        McVertexZ    = (Float_t)mMcVertexZ;
        NMcTrack     = (Float_t)mNMcTrack;
        NMatchedPair = (Float_t)mNMatchedPair;

        hRefMult->Fill(RefMult);
        hVyMcVxMc->Fill(McVertexX,McVertexY);
        hVzMc->Fill(McVertexZ);
        hViRc_Mc[0]->Fill(mVertexX-McVertexX);
        hViRc_Mc[1]->Fill(mVertexY-McVertexY);
        hViRc_Mc[2]->Fill(mVertexZ-McVertexZ);
        hNrMc_RefMult->Fill(RefMult,NMcTrack);

        // loop and fill McTracks variables
        for ( Int_t entry = 0; entry < mNMcTrack; ++entry ) {
          ntrks++;
          if ( !(ntrks%1000) ) printf(".");
          if ( !(ntrks%10000) ) printf(" %dk\n",(Int_t)ntrks/10000);

          PtMc           = (Float_t)mMcTracks_mPtMc[entry];
          PzMc           = (Float_t)mMcTracks_mPzMc[entry];
          EtaMc          = (Float_t)mMcTracks_mEtaMc[entry];
          PhiMc          = (Float_t)mMcTracks_mPhiMc[entry];
          NHitMc         = (Float_t)mMcTracks_mNHitMc[entry];
          GeantId        = (Float_t)mMcTracks_mGeantId[entry];
          Key            = (Float_t)mMcTracks_mKey[entry];
          ParentKey      = (Float_t)mMcTracks_mParentKey[entry];
          ParentGeantId  = (Float_t)mMcTracks_mParentGeantId[entry];
      
//*********************by xuj*****************************************
 //           if(fabs(EtaMc)<1&&fabs(EtaMc)>0.5) 


          McTrack_NTDataArray[0]  = EventId       ;
          McTrack_NTDataArray[1]  = RunId         ;
          McTrack_NTDataArray[2]  = RefMult    ;
          McTrack_NTDataArray[3]  = RefMultCorrected    ;
          McTrack_NTDataArray[4]  = CentralMult   ;
          McTrack_NTDataArray[5]  = CentralityWeight    ;
          McTrack_NTDataArray[6]  = Centrality16    ;
          McTrack_NTDataArray[7]  = Centrality    ;
          McTrack_NTDataArray[8]  = NMcNch        ;
          McTrack_NTDataArray[9]  = NMcHminus     ;
          McTrack_NTDataArray[10]  = VertexX       ;
          McTrack_NTDataArray[11]  = VertexY       ;
          McTrack_NTDataArray[12]  = VertexZ       ;
          McTrack_NTDataArray[13] = McVertexX     ;
          McTrack_NTDataArray[14] = McVertexY     ;
          McTrack_NTDataArray[15] = McVertexZ     ;
          McTrack_NTDataArray[16] = NMcTrack      ;
          McTrack_NTDataArray[17] = NMatchedPair  ;
          McTrack_NTDataArray[18] = PtMc           ;
          McTrack_NTDataArray[19] = PzMc           ;
          McTrack_NTDataArray[20] = EtaMc          ;
          McTrack_NTDataArray[21] = PhiMc          ;
          McTrack_NTDataArray[22] = NHitMc         ;
          McTrack_NTDataArray[23] = GeantId        ;
          McTrack_NTDataArray[24] = Key            ;
          McTrack_NTDataArray[25] = ParentKey      ;
          McTrack_NTDataArray[26] = ParentGeantId  ;
          McTrack_NT->Fill(McTrack_NTDataArray);
        }

        // loop and fill MatchedPairs variables
        for ( Int_t entry = 0; entry < mNMatchedPair; ++entry ) {
          ntrks2++;
          if ( !(ntrks2%1000) ) printf(".");
          if ( !(ntrks2%10000) ) printf(" %dk\n",(Int_t)ntrks2/10000);

          PtMc           = (Float_t)mMatchedPairs_mPtMc[entry];
          PzMc           = (Float_t)mMatchedPairs_mPzMc[entry];
          EtaMc          = (Float_t)mMatchedPairs_mEtaMc[entry];
          PhiMc          = (Float_t)mMatchedPairs_mPhiMc[entry];
          NHitMc         = (Float_t)mMatchedPairs_mNHitMc[entry];
          GeantId        = (Float_t)mMatchedPairs_mGeantId[entry];
          Key            = (Float_t)mMatchedPairs_mKey[entry];
          ParentKey      = (Float_t)mMatchedPairs_mParentKey[entry];
          ParentGeantId  = (Float_t)mMatchedPairs_mParentGeantId[entry];
          NCommonHit  = (Float_t)mMatchedPairs_mNCommonHit[entry];
          PtPr        = (Float_t)mMatchedPairs_mPtPr[entry];
          PtGl        = (Float_t)mMatchedPairs_mPtGl[entry];
          EtaPr       = (Float_t)mMatchedPairs_mEtaPr[entry];
          PhiPr       = (Float_t)mMatchedPairs_mPhiPr[entry];
          DcaGl       = (Float_t)mMatchedPairs_mDcaGl[entry];
          DcaXYGl     = (Float_t)mMatchedPairs_mDcaXYGl[entry];
          DcaZGl      = (Float_t)mMatchedPairs_mDcaZGl[entry];
          DcaXYGlMcV  = (Float_t)mMatchedPairs_mDcaXYGlMcV[entry];
          DcaZGlMcV   = (Float_t)mMatchedPairs_mDcaZGlMcV[entry];
          Flag        = (Float_t)mMatchedPairs_mFlag[entry];
          Dedx        = (Float_t)mMatchedPairs_mDedx[entry];
          RecoKey     = (Float_t)mMatchedPairs_mRecoKey[entry];
          FitPts      = (Float_t)mMatchedPairs_mFitPts[entry];
          DedxPts     = (Float_t)mMatchedPairs_mDedxPts[entry];
          AllPts      = (Float_t)mMatchedPairs_mAllPts[entry];
          NPossible   = (Float_t)mMatchedPairs_mNPossible[entry];
	  ErrP        = (Float_t)mMatchedPairs_mErrP[entry][3];

          MatchedPairs_NTDataArray[0]  = EventId       ;
          MatchedPairs_NTDataArray[1]  = RunId         ;
          MatchedPairs_NTDataArray[2]  = RefMult    ;
          MatchedPairs_NTDataArray[3]  = RefMultCorrected    ;
          MatchedPairs_NTDataArray[4]  = CentralMult   ;
          MatchedPairs_NTDataArray[5]  = CentralityWeight    ;
          MatchedPairs_NTDataArray[6]  = Centrality16    ;
          MatchedPairs_NTDataArray[7]  = Centrality    ;
          MatchedPairs_NTDataArray[8]  = NMcNch        ;
          MatchedPairs_NTDataArray[9]  = NMcHminus     ;
          MatchedPairs_NTDataArray[10]  = VertexX       ;
          MatchedPairs_NTDataArray[11]  = VertexY       ;
          MatchedPairs_NTDataArray[12]  = VertexZ       ;
          MatchedPairs_NTDataArray[13] = McVertexX     ;
          MatchedPairs_NTDataArray[14] = McVertexY     ;
          MatchedPairs_NTDataArray[15] = McVertexZ     ;
          MatchedPairs_NTDataArray[16] = NMcTrack      ;
          MatchedPairs_NTDataArray[17] = NMatchedPair  ;
          MatchedPairs_NTDataArray[18] = PtMc           ;
          MatchedPairs_NTDataArray[19] = PzMc           ;
          MatchedPairs_NTDataArray[20] = EtaMc          ;
          MatchedPairs_NTDataArray[21] = PhiMc          ;
          MatchedPairs_NTDataArray[22] = NHitMc         ;
          MatchedPairs_NTDataArray[23] = GeantId        ;
          MatchedPairs_NTDataArray[24] = Key            ;
          MatchedPairs_NTDataArray[25] = ParentKey      ;
          MatchedPairs_NTDataArray[26] = ParentGeantId  ;
          MatchedPairs_NTDataArray[27] = NCommonHit;
          MatchedPairs_NTDataArray[28] = PtPr      ;
          MatchedPairs_NTDataArray[29] = PtGl      ;
          MatchedPairs_NTDataArray[30] = EtaPr     ;
          MatchedPairs_NTDataArray[31] = PhiPr     ;
          MatchedPairs_NTDataArray[32] = DcaGl     ;
          MatchedPairs_NTDataArray[33] = DcaXYGl   ;
          MatchedPairs_NTDataArray[34] = DcaZGl    ;
          MatchedPairs_NTDataArray[35] = DcaXYGlMcV;
          MatchedPairs_NTDataArray[36] = DcaZGlMcV ;
          MatchedPairs_NTDataArray[37] = Flag      ;
          MatchedPairs_NTDataArray[38] = Dedx      ;
          MatchedPairs_NTDataArray[39] = RecoKey   ;
          MatchedPairs_NTDataArray[40] = FitPts    ;
          MatchedPairs_NTDataArray[41] = DedxPts   ;
          MatchedPairs_NTDataArray[42] = AllPts    ;
          MatchedPairs_NTDataArray[43] = NPossible ;
          MatchedPairs_NTDataArray[44] = ErrP      ;
          MatchedPairs_NT->Fill(MatchedPairs_NTDataArray);
        }

      }
    }


  }
  cout << "   ntrks mc = " << ntrks << endl;
  cout << "   ntrks matched = " << ntrks2 << endl;


  out->cd();
  hEvents->Write();
  hRefMult->Write();
  hVzMc->Write();
  hVyMcVxMc->Write();
  for ( Int_t i = 0; i < 3; ++i ) hViRc_Mc[i]->Write();
  hNrMc_RefMult->Write();
  McTrack_NT->AutoSave("SaveSelf");
  MatchedPairs_NT->AutoSave("SaveSelf");
  out->WriteStreamerInfo();
  out->Close();

  cout << "output written." << endl;
}

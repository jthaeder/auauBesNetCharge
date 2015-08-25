

void getRandomEffFactors() {

  TString basePath("/global/homes/j/jthaeder/analysis/15GeV/NetCharge/data");


  TFile* fileAll = TFile::Open(Form("%s/2015-05-20_0.5/Sum_NetCharge_AuAu14.5GeV_Vz50.root", basePath.Data()));
  TFile* fileSub = TFile::Open(Form("%s/2015-06-03_delta_0.5_0/Sum_NetCharge_AuAu14.5GeV_Vz50.root", basePath.Data()));

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


  double randomEff[2][9];

  // -----------------------------------------------------------
  for (int idxCent = 0; idxCent <9 ; ++idxCent) {

    TH1D *prPosAll = hDistposAll->ProjectionY(Form("%s_%d", hDistposAll->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");
    TH1D *prPosSub = hDistposSub->ProjectionY(Form("%s_%d", hDistposSub->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");

    TH1D *prNegAll = hDistnegAll->ProjectionY(Form("%s_%d", hDistnegAll->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");
    TH1D *prNegSub = hDistnegSub->ProjectionY(Form("%s_%d", hDistnegSub->GetName(), idxCent+1), idxCent+1, idxCent+1, "e");

    randomEff[0][idxCent] = prNegSub->GetMean()/ prNegAll->GetMean();
    randomEff[1][idxCent] = prPosSub->GetMean()/ prPosAll->GetMean();
  } // for (int idxCent = 0; idxCent <9 ; ++idxCent) {

  // -----------------------------------------------------------

  // -----------------------------------------------------------

}

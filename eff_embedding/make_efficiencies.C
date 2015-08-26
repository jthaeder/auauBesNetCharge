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

  int Nbins = 20;
  double pt_bin[21]={0.0,0.1,0.2,0.3,0.4,
		     0.5,0.6,0.7,0.8,0.9,
		     1.0,1.1,1.4,1.7,2.0,
		     2.3,2.7,3.5,4.0,4.5,5.0};

  float max_z = 30;

  float min_pt = 0.2;
  float max_pt = 2.0;

  float min_dedx    = 5;  // default - was 10 , try no 10
  float min_nCommon = 10;
  float min_fitpts  = 20;//15 // 10
  float min_fitpts_nposs = .52;
  float max_eta     = 0.5;
  float max_dca     = 1.0;
  float PID = 0;
  char cbuff[10];
  char buffer[100];
  float vyo,vxo;

  float max_r;

  //if(fabs(energy-200)>0.01){
  //  cout << "change dca cut" << endl;
  //  return;
  //}

  if(fabs(energy-200)<0.01)     max_z = 30;
  else if(fabs(energy-62)<0.01) max_z = 30;
  else if(fabs(energy-39)<0.01) max_z = 30;
  else if(fabs(energy-27)<0.01) max_z = 30;
  else if(fabs(energy-19)<0.01) max_z = 30;
  else if(fabs(energy-14)<0.01) max_z = 30;
  else if(fabs(energy-11)<0.01) max_z = 30;
  else if(fabs(energy-7)<0.01)  max_z = 30;
  else {
    cout << "bad energy" << endl;
    return;
  }

#if 0
  if(fabs(energy-200)<0.01){
    vyo = 0.05141;
    vxo = 0.4738; 
  }
  else if(fabs(energy-62)<0.01){
    vyo = -0.04843;
    vxo = 0.3912;
  }
  else if(fabs(energy-39)<0.01){
    vyo = -0.07631;
    vxo = 0.2935;
  }
  else if(fabs(energy-27)<0.01){
    vyo = -0.07332;
    vxo = 0.3123;
  }
  else if(fabs(energy-19)<0.01){
    vyo = -0.03506;
    vxo = 0.375;
  }
  else if(fabs(energy-14)<0.01){
    vyo = -0.89;
    vxo = 0;
  }
  else if(fabs(energy-11)<0.01){
    vyo = -0.2277;
    vxo = 0.1044;
  }
  else if(fabs(energy-7)<0.01){
    vyo = -0.13;
    vxo = 0.2812;
  }
  else {
    cout << "bad energy" << endl;
    return;
  }
#else
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
#endif


  sprintf(cbuff, "%d",energy);
  TString NRG = TString(cbuff);
  
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

  TString out_file, in_file;
  TString prefix; 
  TString postfix; 
  prefix = TString("embeddingTrees");
  postfix = TString("efficiency");
  double mass;

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
  
  cout << "opening input file: " << in_file << endl;
  TFile *f = new TFile(in_file, "READ");
  TNtuple *MatchedPairs = (TNtuple*) f->Get("MatchedPairs_NT");
  TNtuple *McTrack = (TNtuple*) f->Get("McTrack_NT");
  
  float Dedx, RefMult, RefMultCorrected, CentralityWeight, Centrality16, VertexX, VertexY, VertexZ, PtMc, PzMc, EtaMc, PhiMc, PtPr, PtGl, EtaPr, PhiPr, DcaGl, DcaZGl, DcaXYGl, Flag, FitPts, DedxPts, AllPts, NPossible, ParentGeantId, GeantId,mErrP, NCommonHit;
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
  
  cout << "opening output file: " << out_file << endl;
  TFile *f_out = new TFile(out_file, "RECREATE");
  TFile* fmomenOut = new TFile("fmomenOut","UPDATE");
  f_out->cd();
  
  TH1D* hDCA = new TH1D("hDCA", "", 100, -5., 5.);
  TH2D* hDCAphi = new TH2D("hDCAphi", "", 100, -5., 5.,96, -pi, pi);
  //TH2D* hDCAphiCent = new TH2D("hDCAphiCent", "", 100, -5., 5.,100, -pi, pi);

  TH1D* hpt_f_0005 = new TH1D("hpt_f_0005", "", Nbins, 0., 5.);
  TH1D* hpt_f_0510 = new TH1D("hpt_f_0510", "", Nbins, 0., 5.);
  TH1D* hpt_f_1020 = new TH1D("hpt_f_1020", "", Nbins, 0., 5.);
  TH1D* hpt_f_2030 = new TH1D("hpt_f_2030", "", Nbins, 0., 5.);
  TH1D* hpt_f_3040 = new TH1D("hpt_f_3040", "", Nbins, 0., 5.);
  TH1D* hpt_f_4050 = new TH1D("hpt_f_4050", "", Nbins, 0., 5.);
  TH1D* hpt_f_5060 = new TH1D("hpt_f_5060", "", Nbins, 0., 5.);
  TH1D* hpt_f_6070 = new TH1D("hpt_f_6070", "", Nbins, 0., 5.);
  TH1D* hpt_f_7080 = new TH1D("hpt_f_7080", "", Nbins, 0., 5.);
  TH1D* hpt_mc_0005 = new TH1D("hpt_mc_0005", "", Nbins, 0., 5.);
  TH1D* hpt_mc_0510 = new TH1D("hpt_mc_0510", "", Nbins, 0., 5.);
  TH1D* hpt_mc_1020 = new TH1D("hpt_mc_1020", "", Nbins, 0., 5.);
  TH1D* hpt_mc_2030 = new TH1D("hpt_mc_2030", "", Nbins, 0., 5.);
  TH1D* hpt_mc_3040 = new TH1D("hpt_mc_3040", "", Nbins, 0., 5.);
  TH1D* hpt_mc_4050 = new TH1D("hpt_mc_4050", "", Nbins, 0., 5.);
  TH1D* hpt_mc_5060 = new TH1D("hpt_mc_5060", "", Nbins, 0., 5.);
  TH1D* hpt_mc_6070 = new TH1D("hpt_mc_6070", "", Nbins, 0., 5.);
  TH1D* hpt_mc_7080 = new TH1D("hpt_mc_7080", "", Nbins, 0., 5.);
  hpt_f_0005->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_0510->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_1020->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_2030->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_3040->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_4050->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_5060->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_6070->GetXaxis()->Set(Nbins,pt_bin);
  hpt_f_7080->GetXaxis()->Set(Nbins,pt_bin);
  
  hpt_mc_0005->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_0510->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_1020->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_2030->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_3040->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_4050->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_5060->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_6070->GetXaxis()->Set(Nbins,pt_bin);
  hpt_mc_7080->GetXaxis()->Set(Nbins,pt_bin);

  TH1D* heta_f_0005 = new TH1D("heta_f_0005", "", 31, -1.5, 1.5);
  TH1D* heta_f_0510 = new TH1D("heta_f_0510", "", 31, -1.5, 1.5);
  TH1D* heta_f_1020 = new TH1D("heta_f_1020", "", 31, -1.5, 1.5);
  TH1D* heta_f_2030 = new TH1D("heta_f_2030", "", 31, -1.5, 1.5);
  TH1D* heta_f_3040 = new TH1D("heta_f_3040", "", 31, -1.5, 1.5);
  TH1D* heta_f_4050 = new TH1D("heta_f_4050", "", 31, -1.5, 1.5);
  TH1D* heta_f_5060 = new TH1D("heta_f_5060", "", 31, -1.5, 1.5);
  TH1D* heta_f_6070 = new TH1D("heta_f_6070", "", 31, -1.5, 1.5);
  TH1D* heta_f_7080 = new TH1D("heta_f_7080", "", 31, -1.5, 1.5);
  TH1D* heta_mc_0005 = new TH1D("heta_mc_0005", "", 31, -1.5, 1.5);
  TH1D* heta_mc_0510 = new TH1D("heta_mc_0510", "", 31, -1.5, 1.5);
  TH1D* heta_mc_1020 = new TH1D("heta_mc_1020", "", 31, -1.5, 1.5);
  TH1D* heta_mc_2030 = new TH1D("heta_mc_2030", "", 31, -1.5, 1.5);
  TH1D* heta_mc_3040 = new TH1D("heta_mc_3040", "", 31, -1.5, 1.5);
  TH1D* heta_mc_4050 = new TH1D("heta_mc_4050", "", 31, -1.5, 1.5);
  TH1D* heta_mc_5060 = new TH1D("heta_mc_5060", "", 31, -1.5, 1.5);
  TH1D* heta_mc_6070 = new TH1D("heta_mc_6070", "", 31, -1.5, 1.5);
  TH1D* heta_mc_7080 = new TH1D("heta_mc_7080", "", 31, -1.5, 1.5);


  TH1D* heffcent_f[9];
  TH1D* heffcent_mc[9];
  //for(int i=0;i<9;i++){
  //  sprintf(buffer, "heffcent_f_%d",i);
  //  heffcent_f[i] = new TH1D(buffer, "", Nbins, 0., 5.);
  //  sprintf(buffer, "heffcent_mc_%d",i);
  //  heffcent_mc[i] = new TH1D(buffer, "", Nbins, 0., 5.);
  //}
  TH2D* hptcent_f = new TH2D("hptcent_f", "", Nbins, 0., 5.,1000,0.5,1000.5);
  TH2D* hptcent_mc = new TH2D("hptcent_mc", "", Nbins, 0., 5.,1000,0.5,1000.5);
  TH2D* hptcent2_f = new TH2D("hptcent2_f", "", Nbins, 0., 5.,9,-0.5,8.5);
  TH2D* hptcent2_mc = new TH2D("hptcent2_mc", "", Nbins, 0., 5.,9,-0.5,8.5);
  TH1D* hcent_f = new TH1D("hcent_f", "", 1000, 0.5, 1000.5);
  TH1D* hcent_mc = new TH1D("hcent_mc", "", 1000, 0.5, 1000.5);
  TH1D* heta_f = new TH1D("heta_f", "", 44, -1.1, 1.1);
  TH1D* heta_mc  = new TH1D("heta_mc", "", 44, -1.1, 1.1);
  TH1D* hvz_f = new TH1D("hvz_f", "", 100, -max_z, max_z);
  TH1D* hvz_mc = new TH1D("hvz_mc", "", 100, -max_z, max_z);
  TH2D* hpteta_f = new TH2D("hpteta_f", "", Nbins, 0., 5.,8,-1.,1.);
  TH2D* hpteta_mc = new TH2D("hpteta_mc", "", Nbins, 0., 5.,8,-1.,1.);
  TH2D* hptvz_f = new TH2D("hptvz_f", "", Nbins, 0., 5.,10,-max_z,max_z);
  TH2D* hptvz_mc = new TH2D("hptvz_mc", "", Nbins, 0., 5.,10,-max_z,max_z);
  TH2D* hetacent_f = new TH2D("hetacent_f", "", 8, -1., 1.,1000,0.5,1000.5);
  TH2D* hetacent_mc = new TH2D("hetacent_mc", "", 8, -1., 1.,1000,0.5,1000.5);
  TH2D* hetavz_f = new TH2D("hetavz_f", "", 80, -1., 1.,110,-max_z,max_z);
  TH2D* hetavz_mc = new TH2D("hetavz_mc", "", 80, -1., 1.,110,-max_z,max_z);
  TH2D* hcentvz_f = new TH2D("hcentvz_f", "", 1000, 0.5, 1000.5,10,-max_z,max_z);
  TH2D* hcentvz_mc = new TH2D("hcentvz_mc", "", 1000, 0.5, 1000.5,10,-max_z,max_z);
  TH2D* hptphi_f = new TH2D("hptphi_f", "", Nbins, 0., 5.,48,-pi,pi);
  TH2D* hptphi_mc = new TH2D("hptphi_mc", "", Nbins, 0., 5.,48,-pi,pi);
  TH2D* hmomenRes = new TH2D("hmomenRes", "", Nbins, 0., 5.,20*Nbins,-3.,6.);
  TH2D* hmomenResDiff = new TH2D("hmomenResDiff", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  TH2D* hmomenResDiffC = new TH2D("hmomenResDiffC", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  TH2D* hmomenResDiffP = new TH2D("hmomenResDiffP", "", Nbins, 0., 5.,20*Nbins,-50.,50.);
  TH1D* hphi_f = new TH1D("hphi_f", "", 96,-pi,pi);
  TH1D* hphi_mc = new TH1D("hphi_mc", "",96,-pi,pi);
  TH2D* hphipt = new TH2D("hphipt","",200,-pi,pi,200,0.,10.);
  TH1D* hNdeltaptC = new TH1D("hNdeltaptC","",100,0.,10.);
  TH1D* hdeltaptC = new TH1D("hdeltaptC","",100,0.,10.);
  TH1D* hNdeltaptP = new TH1D("hNdeltaptP","",100,0.,10.);
  TH1D* hdeltaptP = new TH1D("hdeltaptP","",100,0.,10.);
  Double_t xbins[38] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
			1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
			2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
			3.0,3.5,4.0,4.5,5.0,6.0,7.0,10.0};
  TH1D* hMomenMCCent = new TH1D("hMomenMCCent","",37,xbins);
  TH1D* hMomenMatchedCent = new TH1D("hMomenMatchedCent","",37,xbins);
  TH1D* hMomenMCPer = new TH1D("hMomenMCPer","",37,xbins);
  TH1D* hMomenMatchedPer = new TH1D("hMomenMatchedPer","",37,xbins);
  
  hNdeltaptC->Sumw2();hdeltaptC->Sumw2(); hMomenMCCent->Sumw2();
  hNdeltaptP->Sumw2();hdeltaptP->Sumw2();hMomenMatchedCent->Sumw2();
  hMomenMCPer->Sumw2();hMomenMatchedPer->Sumw2();hmomenRes->Sumw2();
  hmomenResDiff->Sumw2();hmomenResDiffC->Sumw2();hmomenResDiffP->Sumw2();
  
  
  cout << "histograms created" << endl;
  
  int ntracks = McTrack->GetEntries();
  int ntrackscut = 0, nmatchedcut = 0;
  int ntrackscutEta = 0, nmatchedcutEta = 0;
  int fP,fG,fR,fVR,fVZ,fE,fF,fD,fDP,fPt,fFP,fFPP,fDCA,fC,fM,fCO;
  int mycent = 0;

  for(int i = 0; i < ntracks; i++){
    fP=1;fG=1;fR=1;fVR=1;fVZ=1;fE=1;fPt=1;fC=0;
    McTrack->GetEntry(i);
    if(pCentrality16==15)mycent=0;
    else if(pCentrality16==14)mycent=1;
    else if(pCentrality16==12 ||pCentrality16==13)mycent=2;
    else if(pCentrality16==10 ||pCentrality16==11)mycent=3;
    else if(pCentrality16==8 ||pCentrality16==9)mycent=4;
    else if(pCentrality16==6 ||pCentrality16==7)mycent=5;
    else if(pCentrality16==4 ||pCentrality16==5)mycent=6;
    else if(pCentrality16==2 ||pCentrality16==3)mycent=7;
    else if(pCentrality16==0 ||pCentrality16==1)mycent=8;
    else fC=1;

    // float   pcm   =  sqrt(pPtMc*pPtMc+pPzMc*pPzMc);
    // float  Ecm=sqrt(pcm*pcm+0.493677*0.493677);
    // float   Ycm = 0.5*log((Ecm+pPzMc)/(Ecm-pPzMc));
    // if(Ycm!=Ycm) Ycm=-999;
    
    if( pParentGeantId == 0 )fP=0;
   
    if( pGeantId == PID )fG=0;
    
    if( pPtMc*pPtMc > min_pt*min_pt && pPtMc*pPtMc < max_pt*max_pt) fPt=0;
    if( (pVertexX-vxo)*(pVertexX-vxo) + (pVertexY-vyo)*(pVertexY-vyo) < max_r*max_r ) fVR=0;

    if( pRefMultCorrected > 0 && pRefMultCorrected < 1000)fR=0;
    if( pVertexZ < max_z && pVertexZ > -max_z )fVZ=0;	    
  
    if( fabs(pEtaMc) < max_eta )fE=0;

    // -- -- -- -- -- -- -- -- -- -- 

    if((fP+fG+fR+fVR+fVZ+fPt)==0) heta_mc->Fill(pEtaMc);
    if((fP+fG+fR+fVR+fVZ+fPt)==0) hpteta_mc->Fill(pPtMc,pEtaMc);
    if((fP+fG+fR+fVR+fVZ+fPt)==0) hetacent_mc->Fill(pEtaMc,pRefMultCorrected);
    if((fP+fG+fR+fVR+fE+fPt)==0)  hvz_mc->Fill(pVertexZ);
    if((fP+fG+fR+fVR+fE+fPt)==0)  hptvz_mc->Fill(pPtMc,pVertexZ);
    if((fP+fG+fR+fVR+fE+fPt)==0)  hcentvz_mc->Fill(pRefMultCorrected,pVertexZ);
    if((fP+fG+fR+fVR+fPt)==0)     hetavz_mc->Fill(pEtaMc,pVertexZ);
    if((fP+fG+fC+fVR+fVZ+fE+fPt)==0) hptcent2_mc->Fill(pPtMc,mycent);

    // --- Fill PT HISTOGRAMS
    if((fP+fG+fR+fVR+fVZ+fE)==0){
      ntrackscut += 1;
      //*****************************************************
      if(mycent==0) hpt_mc_0005->Fill(pPtMc);
      if(mycent==1) hpt_mc_0510->Fill(pPtMc);
      if(mycent==2) hpt_mc_1020->Fill(pPtMc);
      if(mycent==3) hpt_mc_2030->Fill(pPtMc);
      if(mycent==4) hpt_mc_3040->Fill(pPtMc);
      if(mycent==5) hpt_mc_4050->Fill(pPtMc);
      if(mycent==6) hpt_mc_5060->Fill(pPtMc);
      if(mycent==7) hpt_mc_6070->Fill(pPtMc);
      if(mycent==8) hpt_mc_7080->Fill(pPtMc);

      //*********************************************************
      hcent_mc->Fill(pRefMultCorrected);
      hptcent_mc->Fill(pPtMc,pRefMultCorrected);
      hptphi_mc->Fill(pPtMc,pPhiMc);
      hphi_mc->Fill(pPhiMc);
    }

    // --- Fill ETA HISTOGRAMS
    if((fP+fG+fR+fVR+fVZ+fPt)==0){
      ntrackscutEta += 1;
      //*****************************************************
      if(mycent==0) heta_mc_0005->Fill(pEtaMc);
      if(mycent==1) heta_mc_0510->Fill(pEtaMc);
      if(mycent==2) heta_mc_1020->Fill(pEtaMc);
      if(mycent==3) heta_mc_2030->Fill(pEtaMc);
      if(mycent==4) heta_mc_3040->Fill(pEtaMc);
      if(mycent==5) heta_mc_4050->Fill(pEtaMc);
      if(mycent==6) heta_mc_5060->Fill(pEtaMc);
      if(mycent==7) heta_mc_6070->Fill(pEtaMc);
      if(mycent==8) heta_mc_7080->Fill(pEtaMc);
    }
  }
  
  int nmatched = MatchedPairs->GetEntries();
  for(int i = 0; i < nmatched; i++){
    fP=1;fG=1;fVR=1;fVZ=1;fE=1;fF=1;fD=1;fDP=1;fPt=1;fFP=1;fFPP=1;fDCA=1;fC=0;fM=1;fCO=1;
    MatchedPairs->GetEntry(i);
    if(Centrality16==15)mycent=0;
    else if(Centrality16==14)mycent=1;
    else if(Centrality16==12 ||Centrality16==13)mycent=2;
    else if(Centrality16==10 ||Centrality16==11)mycent=3;
    else if(Centrality16==8 ||Centrality16==9)mycent=4;
    else if(Centrality16==6 ||Centrality16==7)mycent=5;
    else if(Centrality16==4 ||Centrality16==5)mycent=6;
    else if(Centrality16==2 ||Centrality16==3)mycent=7;
    else if(Centrality16==0 ||Centrality16==1)mycent=8;
    else fC=1;
	 
    // float   pcm   =  sqrt(PtMc*PtMc+PzMc*PzMc);
    // float   Ecm=sqrt(pcm*pcm+0.493677*0.493677);
    // float   Ycm = 0.5*log((Ecm+PzMc)/(Ecm-PzMc));
    // if(Ycm!=Ycm) Ycm=-999;

    if( ParentGeantId == 0 )fP=0;
    if( GeantId == PID )fG=0;
    if( Dedx >= 0 )fD=0;
    if( Flag >= 0 && Flag <=699 )fF=0;
    if( PtPr<10./7.*PtGl&&PtPr>7./10.*PtGl)fM=0;
    if( VertexZ < max_z && VertexZ > -max_z )fVZ=0;

    if( PtMc*PtMc > min_pt*min_pt && PtMc*PtMc < max_pt*max_pt) fPt=0;
    if( (VertexX-vxo)*(VertexX-vxo) + (VertexY-vyo)*(VertexY-vyo) < max_r*max_r )fVR=0;

    if( DedxPts > min_dedx )fDP=0;
    if( FitPts > min_fitpts )fFP=0;
    if( FitPts/NPossible > min_fitpts_nposs )fFPP=0;
    //if(PtMc<1.&& fabs(DcaGl) < max_dca )fDCA=0;
    if( fabs(DcaGl) < max_dca )fDCA=0;
    //if(PtMc>=1.&&PtMc<=2.&&fabs(DcaGl)<2.-0.5*PtMc)fDCA=0;
    //if(PtMc>2.&&fabs(DcaGl) < 1.5 )fDCA=0;

    if( fabs(EtaMc) < max_eta ) fE=0;

    if( NCommonHit>min_nCommon)fCO=0;

    // -- -- -- -- -- -- -- -- -- -- 

    if((fP+fG+fVR+fF+fD+fDP+fPt+fFP+fCO+fFPP+fDCA+fM)==0)
      hetavz_f->Fill(EtaMc,VertexZ);
    
    if((fP+fG+fVR+fE+fF+fD+fDP+fPt+fFP+fCO+fFPP+fDCA+fM)==0){
      hvz_f->Fill(VertexZ);
      hptvz_f->Fill(PtMc,VertexZ);
      hcentvz_f->Fill(RefMultCorrected,VertexZ);
    }
    
    if((fP+fG+fVR+fVZ+fF+fD+fDP+fPt+fFP+fFPP+fCO+fDCA+fM)==0){
      heta_f->Fill(EtaMc);
      hpteta_f->Fill(PtMc,EtaMc);
      hetacent_f->Fill(EtaMc,RefMultCorrected);
    }
    if((fP+fG+fVR+fVZ+fE+fF+fD+fDP+fPt+fFP+fFPP+fCO+fDCA+fC+fM)==0)hptcent2_f->Fill(PtMc,mycent);

    // --- Fill PT HISTOGRAMS
    if((fP+fG+fVR+fVZ+fE+fF+fD+fDP+fFP+fFPP+fDCA+fCO+fM)==0){
      nmatchedcut += 1;
      //********************************************************************
      if(mycent==0) hpt_f_0005->Fill(PtMc);
      if(mycent==1) hpt_f_0510->Fill(PtMc);
      if(mycent==2) hpt_f_1020->Fill(PtMc);
      if(mycent==3) hpt_f_2030->Fill(PtMc);
      if(mycent==4) hpt_f_3040->Fill(PtMc);
      if(mycent==5) hpt_f_4050->Fill(PtMc);
      if(mycent==6) hpt_f_5060->Fill(PtMc);
      if(mycent==7) hpt_f_6070->Fill(PtMc);
      if(mycent==8) hpt_f_7080->Fill(PtMc);
      //**********************************************************************
      if(FitPts==NCommonHit)hmomenRes->Fill(PtMc,1./PtPr-1./PtMc);
      hmomenResDiff->Fill(PtMc,PtPr-PtMc);
      //if(PtPr>0.&&PtPr<10.&&PtPr<10./7.*PtMc&&PtPr>7./10.*PtMc){
      if(mycent==0){
	//Double_t undptC = (double)hdeltaptC->GetBinError(hdeltaptC->FindBin(PtMc));
	hmomenResDiffC->Fill(PtMc,PtPr-PtMc);
	
	hNdeltaptC->Fill(PtMc);
	hdeltaptC->Fill(PtMc,(PtMc-PtPr));
      }
    }
    
    // --- Fill PT HISTOGRAMS
    if((fP+fG+fVR+fVZ+fF+fD+fDP+fPt+fFP+fFPP+fDCA+fCO+fM)==0){
      nmatchedcutEta += 1;
      //********************************************************************
      if(mycent==0) heta_f_0005->Fill(EtaMc);
      if(mycent==1) heta_f_0510->Fill(EtaMc);
      if(mycent==2) heta_f_1020->Fill(EtaMc);
      if(mycent==3) heta_f_2030->Fill(EtaMc);
      if(mycent==4) heta_f_3040->Fill(EtaMc);
      if(mycent==5) heta_f_4050->Fill(EtaMc);
      if(mycent==6) heta_f_5060->Fill(EtaMc);
      if(mycent==7) heta_f_6070->Fill(EtaMc);
      if(mycent==8) heta_f_7080->Fill(EtaMc);
      }
      if(mycent==7||mycent==8){
	hNdeltaptP->Fill(PtMc);
	hdeltaptP->Fill(PtMc,(PtMc-PtPr));
	hmomenResDiffP->Fill(PtMc,PtPr-PtMc);
      }
      
      hcent_f->Fill(RefMultCorrected);
      hptcent_f->Fill(PtMc,RefMultCorrected);  
      hptphi_f->Fill(PtMc,PhiPr);
      hphi_f->Fill(PhiPr);
      hphipt->Fill(PhiPr,PtMc);
      hDCA->Fill(DcaXYGl);
      hDCAphi->Fill(DcaXYGl,PhiPr);
  }
  cout << "nmatched: " << nmatched << endl;
  cout << "ntracks: " << ntracks << endl;
  cout << "nmatchedcut (pt) : " << nmatchedcut << endl;
  cout << "nmatchedcut (eta): " << nmatchedcutEta << endl;
  cout << "ntrackscut (pt)  : " << ntrackscut << endl;
  cout << "ntrackscut (eta) : " << ntrackscutEta << endl;
  
  // =============================================================================================
  // =============================================================================================
  // =============================================================================================
  // =============================================================================================

  TH1D* hpt_0005 = (TH1D*)hpt_mc_0005->Clone("hpt_0005");
  TH1D* hpt_0510 = (TH1D*)hpt_mc_0510->Clone("hpt_0510");
  TH1D* hpt_1020 = (TH1D*)hpt_mc_1020->Clone("hpt_1020");
  TH1D* hpt_2030 = (TH1D*)hpt_mc_2030->Clone("hpt_2030");
  TH1D* hpt_3040 = (TH1D*)hpt_mc_3040->Clone("hpt_3040");
  TH1D* hpt_4050 = (TH1D*)hpt_mc_4050->Clone("hpt_4050");
  TH1D* hpt_5060 = (TH1D*)hpt_mc_5060->Clone("hpt_5060");
  TH1D* hpt_6070 = (TH1D*)hpt_mc_6070->Clone("hpt_6070");
  TH1D* hpt_7080 = (TH1D*)hpt_mc_7080->Clone("hpt_7080");

  TH1D* heta_0005 = (TH1D*)heta_mc_0005->Clone("heta_0005");
  TH1D* heta_0510 = (TH1D*)heta_mc_0510->Clone("heta_0510");
  TH1D* heta_1020 = (TH1D*)heta_mc_1020->Clone("heta_1020");
  TH1D* heta_2030 = (TH1D*)heta_mc_2030->Clone("heta_2030");
  TH1D* heta_3040 = (TH1D*)heta_mc_3040->Clone("heta_3040");
  TH1D* heta_4050 = (TH1D*)heta_mc_4050->Clone("heta_4050");
  TH1D* heta_5060 = (TH1D*)heta_mc_5060->Clone("heta_5060");
  TH1D* heta_6070 = (TH1D*)heta_mc_6070->Clone("heta_6070");
  TH1D* heta_7080 = (TH1D*)heta_mc_7080->Clone("heta_7080");

  TH1D* hcent = (TH1D*)hcent_mc->Clone("hcent");
  TH1D* heta = (TH1D*)heta_mc->Clone("heta");
  TH1D* hvz = (TH1D*)hvz_mc->Clone("hvz");
  TH2D* hptcent = (TH2D*)hptcent_mc->Clone("hptcent");
  TH2D* hpteta = (TH2D*)hpteta_mc->Clone("hpteta");
  TH2D* hptvz = (TH2D*)hptvz_mc->Clone("hptvz");
  TH2D* hetacent = (TH2D*)hetacent_mc->Clone("hetacent");
  TH2D* hetavz = (TH2D*)hetavz_mc->Clone("hetavz");
  TH2D* hcentvz = (TH2D*)hcentvz_mc->Clone("hcentvz");
  TH2D* hptphi = (TH2D*)hptphi_mc->Clone("hptphi");
  TH1D* hphi = (TH1D*)hphi_mc->Clone("hphi");
  TH2D* hptcent2 = (TH2D*)hptcent2_mc->Clone("hptcent2");
  
  for(int i = 1; i <= Nbins; i++) {
    double n0 = hpt_0005->GetBinContent(i);
    double n1 = hpt_0510->GetBinContent(i);
    double n2 = hpt_1020->GetBinContent(i);
    double n3 = hpt_2030->GetBinContent(i);
    double n4 = hpt_3040->GetBinContent(i);
    double n5 = hpt_4050->GetBinContent(i);
    double n6 = hpt_5060->GetBinContent(i);
    double n7 = hpt_6070->GetBinContent(i);
    double n8 = hpt_7080->GetBinContent(i);
    if( n0 == 0 ) n0 = 1.;
    double p = hpt_f_0005->GetBinContent(i)/(double)n0;
    hpt_0005->SetBinContent(i, p);
    hpt_0005->SetBinError(i, sqrt((double)n0*p*(1.0-p))/(double)n0);
    if( n1 == 0 ) n1 = 1.;
    double p = hpt_f_0510->GetBinContent(i)/(double)n1;
    hpt_0510->SetBinContent(i, p);
    hpt_0510->SetBinError(i, sqrt((double)n1*p*(1.0-p))/(double)n1);
    if( n2 == 0 ) n2 = 1.;
    double p = hpt_f_1020->GetBinContent(i)/(double)n2;
    hpt_1020->SetBinContent(i, p);
    hpt_1020->SetBinError(i, sqrt((double)n2*p*(1.0-p))/(double)n2);
    if( n3 == 0 ) n3 = 1.;
    double p = hpt_f_2030->GetBinContent(i)/(double)n3;
    hpt_2030->SetBinContent(i, p);
    hpt_2030->SetBinError(i, sqrt((double)n3*p*(1.0-p))/(double)n3);
    if( n4 == 0 ) n4 = 1.;
    double p = hpt_f_3040->GetBinContent(i)/(double)n4;
    hpt_3040->SetBinContent(i, p);
    hpt_3040->SetBinError(i, sqrt((double)n4*p*(1.0-p))/(double)n4);
    if( n5 == 0 ) n5 = 1.;
    double p = hpt_f_4050->GetBinContent(i)/(double)n5;
    hpt_4050->SetBinContent(i, p);
    hpt_4050->SetBinError(i, sqrt((double)n5*p*(1.0-p))/(double)n5);
    if( n6 == 0 ) n6 = 1.;
    double p = hpt_f_5060->GetBinContent(i)/(double)n6;
    hpt_5060->SetBinContent(i, p);
    hpt_5060->SetBinError(i, sqrt((double)n6*p*(1.0-p))/(double)n6);
    if( n7 == 0 ) n7 = 1.;
    double p = hpt_f_6070->GetBinContent(i)/(double)n7;
    hpt_6070->SetBinContent(i, p);
    hpt_6070->SetBinError(i, sqrt((double)n7*p*(1.0-p))/(double)n7);
    if( n8 == 0 ) n8 = 1.;
    double p = hpt_f_7080->GetBinContent(i)/(double)n8;
    hpt_7080->SetBinContent(i, p);
    hpt_7080->SetBinError(i, sqrt((double)n8*p*(1.0-p))/(double)n8);    
  }
  // -- -- -- -- -- -- -- -- -- 
  for(int i = 1; i <= 31; i++) {
    double n0 = heta_0005->GetBinContent(i);
    double n1 = heta_0510->GetBinContent(i);
    double n2 = heta_1020->GetBinContent(i);
    double n3 = heta_2030->GetBinContent(i);
    double n4 = heta_3040->GetBinContent(i);
    double n5 = heta_4050->GetBinContent(i);
    double n6 = heta_5060->GetBinContent(i);
    double n7 = heta_6070->GetBinContent(i);
    double n8 = heta_7080->GetBinContent(i);
    if( n0 == 0 ) n0 = 1.;
    double p = heta_f_0005->GetBinContent(i)/(double)n0;
    heta_0005->SetBinContent(i, p);
    heta_0005->SetBinError(i, sqrt((double)n0*p*(1.0-p))/(double)n0);
    if( n1 == 0 ) n1 = 1.;
    double p = heta_f_0510->GetBinContent(i)/(double)n1;
    heta_0510->SetBinContent(i, p);
    heta_0510->SetBinError(i, sqrt((double)n1*p*(1.0-p))/(double)n1);
    if( n2 == 0 ) n2 = 1.;
    double p = heta_f_1020->GetBinContent(i)/(double)n2;
    heta_1020->SetBinContent(i, p);
    heta_1020->SetBinError(i, sqrt((double)n2*p*(1.0-p))/(double)n2);
    if( n3 == 0 ) n3 = 1.;
    double p = heta_f_2030->GetBinContent(i)/(double)n3;
    heta_2030->SetBinContent(i, p);
    heta_2030->SetBinError(i, sqrt((double)n3*p*(1.0-p))/(double)n3);
    if( n4 == 0 ) n4 = 1.;
    double p = heta_f_3040->GetBinContent(i)/(double)n4;
    heta_3040->SetBinContent(i, p);
    heta_3040->SetBinError(i, sqrt((double)n4*p*(1.0-p))/(double)n4);
    if( n5 == 0 ) n5 = 1.;
    double p = heta_f_4050->GetBinContent(i)/(double)n5;
    heta_4050->SetBinContent(i, p);
    heta_4050->SetBinError(i, sqrt((double)n5*p*(1.0-p))/(double)n5);
    if( n6 == 0 ) n6 = 1.;
    double p = heta_f_5060->GetBinContent(i)/(double)n6;
    heta_5060->SetBinContent(i, p);
    heta_5060->SetBinError(i, sqrt((double)n6*p*(1.0-p))/(double)n6);
    if( n7 == 0 ) n7 = 1.;
    double p = heta_f_6070->GetBinContent(i)/(double)n7;
    heta_6070->SetBinContent(i, p);
    heta_6070->SetBinError(i, sqrt((double)n7*p*(1.0-p))/(double)n7);
    if( n8 == 0 ) n8 = 1.;
    double p = heta_f_7080->GetBinContent(i)/(double)n8;
    heta_7080->SetBinContent(i, p);
    heta_7080->SetBinError(i, sqrt((double)n8*p*(1.0-p))/(double)n8);    
  }

  for(int i = 1; i <= Nbins; i++) {
    for(int j =1; j<=1000;j++)
      {
	double nn = hptcent->GetBinContent(i,j);
	//if(j<10)double nn2 = hptcent2->GetBinContent(i,j);
	if( nn <0.00001  ) nn = 1.;
	//if( nn2 <0.00001  ) nn2 = 1.;
	double pp = hptcent_f->GetBinContent(i,j)/(double)nn;
	//if(j<10)double pp2 = hptcent2_f->GetBinContent(i,j)/(double)nn2;
	hptcent->SetBinContent(i,j, pp);
	//if(j<10)hptcent2->SetBinContent(i,j, pp2);
	hptcent->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
	//if(j<10)hptcent2->SetBinError(i,j, sqrt((double)nn2*pp2*(1.0-pp2))/(double)nn2);
      }
    for(int j =1; j<=9;j++){
      //double nn = hptcent->GetBinContent(i,j);
      double nn2 = hptcent2->GetBinContent(i,j);
      //if( nn <0.00001  ) nn = 1.;
      if( nn2 <0.00001  ) nn2 = 1.;
      //double pp = hptcent_f->GetBinContent(i,j)/(double)nn;
      double pp2 = hptcent2_f->GetBinContent(i,j)/(double)nn2;
      //hptcent->SetBinContent(i,j, pp);
      hptcent2->SetBinContent(i,j, pp2);
      //hptcent->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
      hptcent2->SetBinError(i,j, sqrt((double)nn2*pp2*(1.0-pp2))/(double)nn2);
    }
    
    for(int j =1; j<=8;j++){
      double nn = hpteta->GetBinContent(i,j);
      if( nn == 0 ) nn = 1.;
      double pp = hpteta_f->GetBinContent(i,j)/(double)nn;
      hpteta->SetBinContent(i,j, pp);
      hpteta->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
    }
    for(int j =1; j<=10;j++){
      double nn = hptvz->GetBinContent(i,j);
      if( nn == 0 ) nn = 1.;
      double pp = hptvz_f->GetBinContent(i,j)/(double)nn;
      hptvz->SetBinContent(i,j, pp);
      hptvz->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
    }
    for(int j =1; j<=48;j++){
      double nn = hptphi->GetBinContent(i,j);
      if( nn < 0.001 ) nn = 1.;
      double pp = hptphi_f->GetBinContent(i,j)/(double)nn;
      hptphi->SetBinContent(i,j, pp);
      hptphi->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
    }
  }
  
  for(int j =1; j<=1000;j++){
    double nn = hcent->GetBinContent(j);
    if( nn == 0 ) nn = 1.;
    double pp = hcent_f->GetBinContent(j)/(double)nn;
    hcent->SetBinContent(j, pp);
    hcent->SetBinError(j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
  }
  
  for(int j =1; j<=44;j++){
    double nn = heta->GetBinContent(j);
    if( nn == 0 ) nn = 1.;
    double pp = heta_f->GetBinContent(j)/(double)nn;
    heta->SetBinContent(j, pp);
    heta->SetBinError(j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
  }
  
  for(int j =1; j<=96;j++){
    double nn = hphi->GetBinContent(j);
    if( nn < 0.001 ) nn = 1.;
    double pp = hphi_f->GetBinContent(j)/(double)nn;
    hphi->SetBinContent(j, pp);
    hphi->SetBinError(j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
  }
  
  for(int i =1; i<=8;i++){
    for(int j = 1; j <= 1000; j++){
      double nn = hetacent->GetBinContent(i,j);
      if( nn == 0 ) nn = 1.;
      double pp = hetacent_f->GetBinContent(i,j)/(double)nn;
      hetacent->SetBinContent(i,j, pp);
      hetacent->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
    }
  }
  
  for(int i =1; i<=80;i++){
    for(int j = 1; j <= 100; j++){
      double nn = hetavz->GetBinContent(i,j);
      if( nn == 0 ) nn = 1.;
      double pp = hetavz_f->GetBinContent(i,j)/(double)nn;
      hetavz->SetBinContent(i,j, pp);
      hetavz->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
    }
  }
  
  for(int i =1; i<=1000;i++){
    for(int j = 1; j <= 10; j++){
      double nn = hcentvz->GetBinContent(i,j);
      if( nn == 0 ) nn = 1.;
      double pp = hcentvz_f->GetBinContent(i,j)/(double)nn;
      hcentvz->SetBinContent(i,j, pp);
      hcentvz->SetBinError(i,j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
    }
  }
  
  for(int j =1; j<=100;j++){
    double nn = hvz->GetBinContent(j);
    if( nn == 0 ) nn = 1.;
    double pp = hvz_f->GetBinContent(j)/(double)nn;
    hvz->SetBinContent(j, pp);
    hvz->SetBinError(j, sqrt((double)nn*pp*(1.0-pp))/(double)nn);
  }
  
  TCanvas *cmomentumResC = new TCanvas("cmomentumResC", "momentumResC", 800, 600);
  hMomenMatchedCent->Draw("C");

  TCanvas *cmomRes = new TCanvas("cmomRes", "momRes", 800, 600);
  hmomenRes->Draw("colz");

  TCanvas *cmomentumResP = new TCanvas("cmomentumResP", "momentumResP", 800, 600);
  hMomenMatchedPer->Draw("C");

  hdeltaptC->Divide(hNdeltaptC);
  hdeltaptP->Divide(hNdeltaptP);

  TCanvas *cdeltaptC = new TCanvas("cdeltaptC", "deltaptC", 800, 600);
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
  //return;
  /*
    TCanvas *cpt = new TCanvas("cpt", "pt", 800, 600);
    sprintf(buffer,"Efficiency for refmult integrated p_{T} distribution for %s at %dGeV",particle,energy);
    hpt->SetTitle(buffer);
    hpt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    hpt->GetXaxis()->CenterTitle();
    hpt->GetXaxis()->SetRangeUser(0.2,4.5);
    hpt->GetYaxis()->SetTitle("Efficiency");
    hpt->GetYaxis()->CenterTitle();
    //hpt->SetMarkerStyle(20);
    hpt->SetStats(kFALSE);
    hpt->DrawCopy();
    sprintf(buffer, "effQA/hpt%d_%d.png",energy,PID);
    cpt->Print(buffer);
  */
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
  
  f_out->Write();
  f_out->Close();
  //f->Close();
}

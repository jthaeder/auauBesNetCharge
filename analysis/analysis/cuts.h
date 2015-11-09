
// ----------------------------------------------------------------------------  
// -- Ranges / Cuts
//    Net-Charge / Net-Proton / Net-Kaon
// ----------------------------------------------------------------------------  
Double_t centBinRange[2]             = {0, 9};

Double_t ptRange[3][2]               = { {0.2, 2.0}, {0.4, 2.0}, {0.2, 1.6}};
Double_t ptMidPoint[3]               = { 1.0, 0.8, 0.4};    //  pseudo (net-charge), to get 2 region pt for f(i,j,k,l)

Double_t etaAbsRange[3][2]           = { {-0.5, 0.5}, {-999.0, 999.0}, {-999., 999.} };
Double_t yAbsRange[3][2]             = { {-999., 999.}, {-999.0, 999.0}, {-999., 999.} };
Double_t dcaMax[3]                   = { 1.  ,  1.  ,  1.};

Double_t nHitsFitMin[3]              = {20.  , 20.  , 20};
Double_t nHitsDedxMin[3]             = { 5.  ,  5.  ,  5.};
Double_t ratioNHitsFitNFitPossMin[3] = { 0.52,  0.52,  0.52};

Double_t dcaMaxRefMult[3]            = { 3.,  3.,  3.};
Double_t nHitsFitRefMult[3]          = {10., 10., 15.};


Double_t etaAbsRangeRefMult[3][2]    = { {0.5, 1.0}, {0., 1.0}, {0., 1.0} }; 

// -- For net-charge only
Double_t ptRangeSpallation[2]        = {0.2, 0.4};
Double_t nSigmaProtonMaxSpallation   = 2;

// -- For net-proton and net-kaon usage
Double_t nSigmaMax[3]                = { -1., 2.0,  2.0};
Double_t mSquareRange[3][2]          = { {-1, -1}, { 0.6, 1.2}, { 0.6, 1.2} };


//***********systematic error cut condition***
/* double dca_cut=1; */
/* double nfit_cut=20; */
/* double Nsigma_cut=2; */
/* double nfitratio_cut=0.52; */
/* double pT_cut[20]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1}; */
/* const int NumberOfPts=10; */
/* //\*********************************** */
/* double Ycut_p=0.5; */
/* double Ycut_k=0.5; */
/* double Ycut_pi=0.5; */
/* //double eta_cut=0.5; */
/* //\***************************** */

/* double Ptcut_p[2]={0.4,0.8}; */
/* double Ptcut_k[2]={0.2,0.6}; */
/* double Ptcut_pi[2]={0.2,1}; */
/* double Ptcut_ch[2]={0.15,1}; */

/* //\********************************** */

/* //\*****eta dependence************* */
/* const int num=1; */
/* double eta_cut[1]={1}; */
/* double y_cut[26]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0}; */
/* const int num_y=26; */

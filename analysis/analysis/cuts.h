
// ----------------------------------------------------------------------------  
// -- Ranges / Cuts
// ----------------------------------------------------------------------------  
Double_t centBinRange[2]           = {0, 9};

Double_t ptRange[2]                = {0.2, 2.0};
Double_t ptMidPoint                = 1.0;    // -- pseudo to get 2 region pt for f(i,j,k,l)

//Double_t etaAbsRange[2]            = {-0.45, 0.45};
Double_t etaAbsRange[2]            = {-0.5, 0.5};
Double_t dcaMax                    = 1.;

Double_t nHitsFitMin               = 20.;
Double_t nHitsDedxMin              = 5.;
Double_t ratioNHitsFitNFitPossMin  = 0.52;

Double_t dcaMaxRefMult             = 3.;
Double_t nHitsFitRefMult           = 10.;
Double_t etaAbsRangeRefMult[2]     = {0.5, 1.0};

Double_t ptRangeSpallation[2]      = {0.2, 0.4};
Double_t nSigmaProtonMaxSpallation = 2;


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

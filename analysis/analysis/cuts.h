
// ----------------------------------------------------------------------------  
// -- Ranges / Cuts
//    Net-Charge / Net-Proton / Net-Kaon
// ----------------------------------------------------------------------------  
Double_t centBinRange[2]             = {0, 9};

Double_t ptRange[3][2]               = { {0.2, 2.0}, {0.4, 2.0}, {0.2, 1.6}};
Double_t ptMidPoint[3]               = { 1.0, 0.8, 0.4};    //  pseudo (net-charge), to get 2 region pt for f(i,j,k,l)

Double_t etaAbsRange[3][2]           = { {-0.5, 0.5},   {-999.0, 999.0}, {-999., 999.} };
Double_t yAbsRange[3][2]             = { {-999., 999.}, {-0.5, 0.5},     {0.5, 0.5} };
Double_t dcaMax[3]                   = { 1.  ,  1.  ,  1.};

Double_t nHitsFitMin[3]              = {20.  , 20.  , 20};
Double_t nHitsDedxMin[3]             = { 5.  ,  5.  ,  5.};  // 10 !!!
Double_t ratioNHitsFitNFitPossMin[3] = { 0.52,  0.52,  0.52};

Double_t dcaMaxRefMult[3]            = { 3.,  3.,  3.};
Double_t nHitsFitRefMult[3]          = {10., 10., 15.};


Double_t etaAbsRangeRefMult[3][2]    = { {0.5, 1.0}, {0., 1.0}, {0., 1.0} }; 

// -- For net-charge only
Double_t ptRangeSpallation[2]        = {0.2, 0.4};
Double_t nSigmaProtonMaxSpallation   = 2;

// -- For net-proton and net-kaon usage
Double_t nSigmaMax[3]                = { -1., 2.0, 2.0};
Double_t mSquareRange[3][2]          = { {-1, -1}, { 0.6, 1.2}, { 0.15, 0.4} };



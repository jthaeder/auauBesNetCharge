void minimc_macro(Long64_t nentries, TString list, TString particle, Int_t energy) {
  gSystem ->Load("StMiniMcTree.so");
  gSystem ->Load("StRefMultCorr.so");

  StMiniMcTree* minimc = new StMiniMcTree(list);
  TString outdir = "embeddingTrees";
  minimc->Loop(nentries, outdir, particle, energy);
}

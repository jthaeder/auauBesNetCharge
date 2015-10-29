void makeFact4(const Char_t* path, const Char_t* name) {


  TFile *fin = TFile::Open(Form("data/%s/%s", path, name));

  TList *l = static_cast<TList*>(fin->Get("fNetCharge"));
  TList *list = static_cast<TList*>(l->FindObject("fNetChargeFijkh"));

  gSystem->Exec(Form("mkdir -p output/%s", path));
  TFile *fout = TFile::Open(Form("output/%s/Fact4.root", path), "RECREATE");
  fout->cd();
  
  for (Int_t cent = 0 ; cent < 10 ; ++cent) {
    for (int ii = 0 ; ii < 9 ; ++ii) {
      for (int jj = 0 ; jj < 9 ; ++jj) {
	for (int kk = 0 ; kk < 9 ; ++kk) {
	  for (int ll = 0 ; ll < 9 ; ++ll) {
	    TProfile * p = static_cast<TProfile*>(list->FindObject(Form("f%d%d%d%d_Cent%d",ii,jj,kk,ll,cent)));
	    if (p)
	      p->Write();
	  }
	}
      }
    }
  }
  
  fout->Close();
}

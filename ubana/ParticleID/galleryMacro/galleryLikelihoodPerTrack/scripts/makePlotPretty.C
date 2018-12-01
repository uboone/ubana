void makePlotPretty(){

  TFile* f = new TFile("likelihoodMaps.root", "READ");

  if (!_file0->IsOpen()) {
    printf("<E> Cannot open input file") ;
    exit(1) ;
  }

  TList* list = _file0->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TObject* obj ;

  TGraph* grp = (TGraph*)_file0->Get("Graph;1");
  TGraph* grmu = (TGraph*)_file0->Get("Graph;2");

  for (int i = 0; i < list->GetSize(); i++){
  
    TKey* key = (TKey*)list->At(i) ;
    obj = key->ReadObj();
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
        && (!obj->InheritsFrom("TH2"))
        && (!obj->InheritsFrom("TH1")) 
       ) {
      printf("<W> Object %s is not 1D or 2D histogram : "
          "will not be converted\n",obj->GetName()) ;
      continue;
    }
    printf("Histo name:%s title:%s\n",obj->GetName(),obj->GetTitle());

    std::cout << "1" << std::endl;
    TH2D* h = (TH2D*)_file0->Get(obj->GetName());

    h->GetXaxis()->SetRangeUser(0,30);
    h->SetMarkerStyle(2);

    std::cout << "2" << std::endl;
    grp->SetLineColor(kRed);
    grp->SetLineWidth(2);

    grmu->SetLineColor(kAzure+1);
    grmu->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
    c1->cd();

    TH2D* h1 = (TH2D*)f->Get("h_combined");
    h1->Draw("colz");
    h->SetMarkerColor(kWhite);
    h->Draw("same");
    //grp->Draw("same");
    //grmu->Draw("same");

    TString st = Form("%s.png", obj->GetName());

    c1->SaveAs(st);


  }
}

void CompareFileContentsHists(std::string file1, std::string file2){

  // Get TFiles
  TFile *f1 = new TFile(file1.c_str(),"open");
  TFile *f2 = new TFile(file2.c_str(),"open");

  if (!f1->IsOpen()) {
    std::cout << "Cannot open input file" << file1 << std::endl;
    exit(1) ;
  }

  if (!f2->IsOpen()) {
    std::cout << "Cannot open input file" << file2 << std::endl;
    exit(1) ;
  }

  TDirectory *d1 = (TDirectory*)f1->Get("pidvalid");
  TDirectory *d2 = (TDirectory*)f2->Get("pidvalid");

  TList* list = d1->GetListOfKeys() ;
  if (!list) { std::cout <<  "No keys found in file" << std::endl ; exit(1) ; }
  TObject* obj ;

  for (int i = 0; i < list->GetSize(); i++){

    TKey* key = (TKey*)list->At(i) ;
    obj = key->ReadObj();
    if ( (!obj->InheritsFrom("TH2"))
        && (!obj->InheritsFrom("TH1"))
       ) {
      std::cout << "Object " << obj->GetName() << " is not 1D or 2D histogram : will not be compared" << std::endl;
      continue;
    }

    std::cout << "Comparing: " << obj->GetName() << std::endl;
    if (obj->InheritsFrom("TH2")){
      TH2F* h1 = (TH2F*)d1->Get(obj->GetName());
      TH2F* h2 = (TH2F*)d2->Get(obj->GetName());

      if (!h2){
        std::cout << "ERROR " << obj->GetName() << " does not exist in file " << file2.c_str() << std::endl;
        continue;
      }

      if (h1->GetSize() != h2->GetSize()) std::cout << "ERROR h1 has " << h1->GetSize() << "-2 bins, h2 has " << h2->GetSize() << "-2 bins" << std::endl;

      for (int i_bin=1; i_bin <= h1->GetSize(); i_bin++){
        if (h1->GetBinContent(i_bin) != h2->GetBinContent(i_bin)) std::cout << "ERROR bin " << i_bin << ": h1 bin content = " << h1->GetBinContent(i_bin) << ", h2 bin content = " << h2->GetBinContent(i_bin) << std::endl;
      }
    }
    else if (obj->InheritsFrom("TH1")){
      TH1F* h1 = (TH1F*)d1->Get(obj->GetName());
      TH1F* h2 = (TH1F*)d2->Get(obj->GetName());

      if (!h2){
        std::cout << "ERROR " << obj->GetName() << " does not exist in file " << file2.c_str() << std::endl;
        continue;
      }

      if (h1->GetSize() != h2->GetSize()) std::cout << "ERROR h1 has " << h1->GetSize() << "-2 bins, h2 has " << h2->GetSize() << "-2 bins" << std::endl;

      for (int i_bin=1; i_bin <= h1->GetSize(); i_bin++){
        if (h1->GetBinContent(i_bin) != h2->GetBinContent(i_bin)) std::cout << "ERROR bin " << i_bin << ": h1 bin content = " << h1->GetBinContent(i_bin) << ", h2 bin content = " << h2->GetBinContent(i_bin) << std::endl;
      }
    }

  }
}

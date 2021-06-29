void getSimPot(){

  double sr_pot;

  TTree* tree = (TTree*)_file0->Get("pidvalid/potTree");

  tree->SetBranchAddress("sr_pot", &sr_pot);

  double totPot = 0;
  for (int i = 0; i < tree->GetEntries(); i++){

    tree->GetEntry(i);

    totPot += sr_pot;
  }

  std::cout << "Total POT: " <<  totPot << std::endl;

}

void dumpRunSubRun(){

  int run;
  int sub_run;

  TTree* tree = (TTree*)_file0->Get("pidvalid/potTree");

  tree->SetBranchAddress("sr_run", &run);
  tree->SetBranchAddress("sr_sub_run", &sub_run);

  std::vector<std::pair<int, int>> rsr_pairs;

  for (int i = 0; i < tree->GetEntries(); i++){

    tree->GetEntry(i);

    std::pair<int, int> rsr_pair;
    rsr_pair.first = run;
    rsr_pair.second = sub_run;
    rsr_pairs.push_back(rsr_pair);

  }

  // open output txt file
  std::ofstream ofs("run_subrun.txt", std::ofstream::out);

  std::vector<std::pair<int, int>> outvec;
  for (int i = 0; i < rsr_pairs.size(); i++){

    // get run subrun info
    int runt1 = rsr_pairs.at(i).first;
    int subrunt1 = rsr_pairs.at(i).second;

    bool isWritten = false;
    for (int j = 0; j < outvec.size(); j++){

      int runt2 = outvec.at(j).first;
      int subrunt2 = outvec.at(j).second;
  
      if (runt1 == runt2 && subrunt1 == subrunt2){
        isWritten = true;
      }

    }
    if (isWritten == false) outvec.push_back(rsr_pairs.at(i));

  }

  for (int i = 0; i < outvec.size(); i++){

    ofs << outvec.at(i).first << " " << outvec.at(i).second << std::endl;

  }

  ofs.close();

}

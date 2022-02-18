
#include <iostream>
#include <fstream>

//Root Includes
#include "TFile.h" 
#include "TTree.h" 

int GetPOT(const char *_file1){

    std::cout << "File: " << _file1 << std::endl;

    // First we need to open the root file
    TFile * f = new TFile(_file1);
    if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }

    // Try the second method of getting the POT
    TTree * pot_tree2 = (TTree*)f->Get("GetPOT/pottree");

    if (pot_tree2 == NULL){
        std::cout << "help can't get the branch so exiting..." << std::endl;
        gSystem->Exit(1);
    }
    double pot_sum2 = 0;
    double pot2;
    int run;
    int subrun;
    pot_tree2->SetBranchAddress("pot", &pot2);
    pot_tree2->SetBranchAddress("run", &run);
    pot_tree2->SetBranchAddress("subrun", &subrun);

    for (int i = 0; i < pot_tree2->GetEntries(); i++){
        pot_tree2->GetEntry(i);
        pot_sum2 += pot2;
    }

    std::cout << "Total POT: " << pot_sum2 << std::endl;
    return 0;
}
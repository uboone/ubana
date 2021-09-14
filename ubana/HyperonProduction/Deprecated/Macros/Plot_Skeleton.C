#define OutputTree_cxx
#include "OutputTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//#include "Functions.h"
//#include "Errors.h"

#include "Alg/FV.h"
#include "Alg/Muon_ID.h"
#include "Alg/Track_Lengths.h"



void OutputTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L OutputTree.C
//      root> OutputTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;


int nSignal = 0;
int nSelected = 0;
int nSignalSelected = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


if(IsSignal) nSignal++;

//apply preselection from module
//if(!SelectedEvent) continue;

/////////////////////////////////////////////////////////

//apply preselection here if you like 

//FV cut
if(!inFV(*RecoPrimaryVertex)) continue;



//3 tracks zero shows
if(NPrimaryTrackDaughters < 3) continue;
if(NPrimaryShowerDaughters > 0) continue;



std::vector<double> lengths;
std::vector<double> PIDs;

for(int i_tr=0;i_tr<TracklikePrimaryDaughters_;i_tr++){

lengths.push_back(TracklikePrimaryDaughters_TrackLength[i_tr]);
PIDs.push_back(TracklikePrimaryDaughters_TrackPID[i_tr]);

}



//muon pid
int i_muon = Muon_ID(lengths,PIDs);
if(i_muon < 0) continue;


//track length cut
if(!Track_Length_Cut(lengths,i_muon)) continue;

/////////////////////////////////////////////////////////


nSelected++;
if(IsSignal) nSignalSelected++;


   }

std::cout << "Signal Events: " << nSignal << std::endl;
std::cout << "Selected Events: " << nSelected << std::endl;
std::cout << "Selected Signal: " << nSignalSelected << std::endl;

double Efficiency = (double)nSignalSelected/nSignal;
double Purity = (double)nSignalSelected/nSelected;
double ExP = Efficiency*Purity;

std::cout << "Efficiency = " << nSignalSelected << "/" << nSignal << " = " << Efficiency << std::endl;
std::cout << "Purity = " << nSignalSelected << "/" << nSelected << " = " << Purity << std::endl;
std::cout << "Efficiency x Purity = " << ExP << std::endl;
std::cout << "Signal/Sqrt(Background) = " << nSignalSelected/sqrt(nSelected-nSignalSelected) << std::endl;




}

#define pfpprofile_ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void pfpprofile_ana::Loop()
{
  if (fChain == 0) return;

  // plot options
  int nbin_longdist = 100;
  double lbin_longdist = 0.0;
  double hbin_longdist = 25.0; // [X0 = 14 cm]
  
  int nbin_trandist = 16;
  double lbin_trandist = 0.0;
  double hbin_trandist = 8.0; // [cm]

  int nbin_energy = 10;
  double lbin_energy = 0.0;
  double hbin_energy = 1000.0; // [MeV]

  int n_particles = 5;
  int n_energies = 10;
  int n_profiles = 2;

  const char * particle_options[5] = {"electron", "photon", "proton", "pion", "others"}; // 0: electron; 1: photon; 2: proton; 3: pion; 4: others
  const char * energy_options[10] = {"0-100MeV", "100-200MeV", "200-300MeV","300-400MeV", "400-500MeV","500-600MeV", "600-700MeV","700-800MeV", "800-900MeV","900-1000MeV"};
  const char * profile_options[7] = {"long", "trans", "trans_1", "trans_2", "trans_3", "trans_4","trans_5"};

  double cut_purity = 0.0;
  double cut_completeness = 0.9;
  
  // output
  TFile output("output_ana.root","recreate");

  TProfile *hlong = new TProfile("hlong","Longitudinal profile", nbin_longdist, lbin_longdist, hbin_longdist);
  TProfile *htran = new TProfile("htran","Transverse profile", nbin_trandist, lbin_trandist, hbin_trandist);

  TProfile *hprofile[5][10][2]; // [particle][energy][profile]
  for (int i=0; i<n_particles; i++) {
    for (int j=0; j<n_energies; j++) {
      for (int k=0; k<n_profiles; k++) {
        if (k == 0) {
          hprofile[i][j][k] = new TProfile(TString::Format("hprofile_%s_%s_%s", particle_options[i], energy_options[j], profile_options[k]), "longitudinal profile; t [X0=14cm]; Q", nbin_longdist, lbin_longdist, hbin_longdist);
        }
        else {
          hprofile[i][j][k] = new TProfile(TString::Format("hprofile_%s_%s_%s", particle_options[i], energy_options[j], profile_options[k]), "transverse profile; d [cm]; Q", nbin_longdist, lbin_longdist, hbin_longdist);
        }
      }
    }
  }
 
  // loop over entries
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries: " << nentries << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    for (int i = 0; i< npfps; ++i){
      for (int p=0; p<3; p++) {
        if (p != 2) continue;

        if (TotalCharge[i][p]) {
          int index_particle = -1;
          int index_energy = -1;

          if ( abs(pfppdg_truth[i]) == 11 ) index_particle = 0;
          else if ( abs(pfppdg_truth[i]) == 22 ) index_particle = 1;
          else if ( abs(pfppdg_truth[i]) == 2212 ) index_particle = 2;
          else if ( abs(pfppdg_truth[i]) == 211 ) index_particle = 3;
          else index_particle = 4;

          for (int ienergy=0; ienergy<nbin_energy; ienergy++) {
            if ( (pfpenergy_recon[i][p] >= ((hbin_energy-lbin_energy)/nbin_energy)*ienergy+lbin_energy) && (pfpenergy_recon[i][p] < ((hbin_energy-lbin_energy)/nbin_energy)*(ienergy+1)+lbin_energy) ) {
              index_energy = ienergy;  
            }
          }
          
          if (index_energy == -1 || index_particle == -1) {
            cout << "index_particle: " << index_particle << endl;
            cout << "index_energy: " << index_energy << endl;
          }
          for (int j=0; j<nbin_longdist; j++) {
            hlong->Fill(hlong->GetBinCenter(j+1), LongProf[i][p][j]);
            if (index_energy != -1) hprofile[index_particle][index_energy][0]->Fill(hlong->GetBinCenter(j+1), LongProf[i][p][j]);
          }
          for (int j=0; j<nbin_trandist; j++) {
            htran->Fill(htran->GetBinCenter(j+1), TranProf[i][p][j]);
            if (index_energy != -1) hprofile[index_particle][index_energy][1]->Fill(htran->GetBinCenter(j+1), TranProf[i][p][j]);
          }
        }
      } // loop over planes p
    } // loop over npfps i

  } // loop over nentries jentry

  // output
  output.Write();
  output.Close();
}

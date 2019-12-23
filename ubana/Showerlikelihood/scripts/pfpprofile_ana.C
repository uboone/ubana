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

  int nbin_resolution_energy = 40;
  double lbin_resolution_energy = -1.0;
  double hbin_resolution_energy = 1.0;

  int nbin_energy_mc = 50;
  double lbin_energy_mc = 0.0;
  double hbin_energy_mc = 5000.0; // [MeV]
  
  int nbin_energy_reco = 50;
  double lbin_energy_reco = 0.0;
  double hbin_energy_reco = 5000.0; // [MeV]

  int n_particles = 6;
  int n_energies = 10;
  int n_profiles = 7;

  const char * particle_options[6] = {"electron", "photon", "muon", "proton", "pion", "others"}; // 0: electron; 1: photon; 2:muon; 3: proton; 4: pion; 5: others
  const char * energy_options[10] = {"0_100MeV", "100_200MeV", "200_300MeV","300_400MeV", "400_500MeV","500_600MeV", "600_700MeV","700_800MeV", "800_900MeV","900_1000MeV"};
  const char * profile_options[7] = {"long", "trans", "trans_1", "trans_2", "trans_3", "trans_4","trans_5"};

  double cut_datahits = 0.1;
  double cut_purity = 0.9;
  double cut_completeness = 0.9;
  
  // output
  TFile output("output_ana.root","recreate");

  TProfile *hlong = new TProfile("hlong","Longitudinal profile", nbin_longdist, lbin_longdist, hbin_longdist);
  TProfile *htran = new TProfile("htran","Transverse profile", nbin_trandist, lbin_trandist, hbin_trandist);

  TProfile *hprofile[6][10][7]; // [particle][energy][profile]

  TH1D *h1_energy_resolution[6][3]; // [particle][plane]
  TH2D *h2_energy_resolution[6][3]; // [particle][plane]

  for (int i=0; i<n_particles; i++) {
    for (int j=0; j<n_energies; j++) {
      for (int k=0; k<n_profiles; k++) {
        if (k == 0) {
          hprofile[i][j][k] = new TProfile(TString::Format("hprofile_%s_%s_%s", particle_options[i], energy_options[j], profile_options[k]), "longitudinal profile; t [X0=14cm]; Q", nbin_longdist, lbin_longdist, hbin_longdist);
        }
        else {
          hprofile[i][j][k] = new TProfile(TString::Format("hprofile_%s_%s_%s", particle_options[i], energy_options[j], profile_options[k]), "transverse profile; d [cm]; Q", nbin_trandist, lbin_trandist, hbin_trandist);
        }
      }
    }

    for (int p=0; p<3; p++) {
      h1_energy_resolution[i][p] = new TH1D(TString::Format("h1_energy_resolution_%s_%d", particle_options[i], p), TString::Format("purity >= %.1f and completeness >= %.1f; Energy resolution (#frac{E_{Reco}-E_{MC}}{E_{MC}}); Arbitrary units", cut_purity, cut_completeness), nbin_resolution_energy, lbin_resolution_energy, hbin_resolution_energy); 
      h2_energy_resolution[i][p] = new TH2D(TString::Format("h2_energy_resolution_%s_%d", particle_options[i], p), TString::Format("purity >= %.1f and completeness >= %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", cut_purity, cut_completeness), nbin_energy_mc, lbin_energy_mc, hbin_energy_mc, nbin_energy_reco, lbin_energy_reco, hbin_energy_reco); 
    }
  }
 
  // energy correction
  int npoint_energy = 12;
  int nbin_energy_correction = 10;
  double step_energy = 200.0; // [MeV]
  TH1D *h1_energy_electron[3][12]; // for energy correction functions: [plane][energy]
  for (int p=0; p<3; p++) {
    for (int n=0; n<npoint_energy; n++) {
      h1_energy_electron[p][n] = new TH1D(TString::Format("h1_energy_electron_%d_%d",p, n), "Electron pfp energy", nbin_energy_correction*(n+1), 0, (n+1)*step_energy);
    }
  }
  
  double energy_corr_slope[3] = {0.654342, 0.624423, 0.824178};
  double energy_corr_intercept[3] ={1.23691, 17.7857, -3.23663};
  bool UseEnergyCorrection = true;

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
      if (pfppurity[i] < cut_purity || pfpcompleteness[i] < cut_completeness) continue;
      if (pfp_datahits[i][0]>= cut_datahits || pfp_datahits[i][1]>=cut_datahits || pfp_datahits[i][2]>=cut_datahits) continue;
      
      // resolution
      for (int p=0; p<3; p++) {
        if (TotalCharge[i][p]) {

          int index_particle = -1;
          if ( abs(pfppdg_truth[i]) == 11 ) index_particle = 0;
          else if ( abs(pfppdg_truth[i]) == 22 ) index_particle = 1;
          else if ( abs(pfppdg_truth[i]) == 13 ) index_particle = 2;
          else if ( abs(pfppdg_truth[i]) == 2212 ) index_particle = 3;
          else if ( abs(pfppdg_truth[i]) == 211 ) index_particle = 4;
          else index_particle = 5;

          if (! UseEnergyCorrection) {
            h1_energy_resolution[index_particle][p]->Fill((pfpenergy_recon[i][p] - pfpenergy_mc[i])/ pfpenergy_mc[i]);
            h2_energy_resolution[index_particle][p]->Fill(pfpenergy_mc[i], pfpenergy_recon[i][p]);
          }
          else {
            h1_energy_resolution[index_particle][p]->Fill(( (pfpenergy_recon[i][p]-energy_corr_intercept[p])/energy_corr_slope[p] - pfpenergy_mc[i])/ pfpenergy_mc[i]);
            h2_energy_resolution[index_particle][p]->Fill(pfpenergy_mc[i], (pfpenergy_recon[i][p]-energy_corr_intercept[p])/energy_corr_slope[p]);
          }

          // energy correction
          if (index_particle == 0 && pfp_primary_e[i] == 1) {
            for (int n=0; n<npoint_energy; n++) {
              if ( (pfpenergy_mc[i] >= n*step_energy) && (pfpenergy_mc[i] <= (n+1)*step_energy) ) {
                h1_energy_electron[p][n]->Fill(pfpenergy_recon[i][p]);
              }
            }
          }
        }
      }
      
      // profile
      for (int p=0; p<3; p++) {
        if (p != 2) continue;

        if (TotalCharge[i][p]) {
          int index_particle = -1;
          int index_energy = -1;

          if ( abs(pfppdg_truth[i]) == 11 ) index_particle = 0;
          else if ( abs(pfppdg_truth[i]) == 22 ) index_particle = 1;
          else if ( abs(pfppdg_truth[i]) == 13 ) index_particle = 2;
          else if ( abs(pfppdg_truth[i]) == 2212 ) index_particle = 3;
          else if ( abs(pfppdg_truth[i]) == 211 ) index_particle = 4;
          else index_particle = 5;

          for (int ienergy=0; ienergy<nbin_energy; ienergy++) {
            if (!UseEnergyCorrection) {
              if ( (pfpenergy_recon[i][p] >= ((hbin_energy-lbin_energy)/nbin_energy)*ienergy+lbin_energy) && (pfpenergy_recon[i][p] < ((hbin_energy-lbin_energy)/nbin_energy)*(ienergy+1)+lbin_energy) ) {
                index_energy = ienergy;  
              }
            }
            else {
              if ( ( (pfpenergy_recon[i][p]-energy_corr_intercept[p])/energy_corr_slope[p] >= ((hbin_energy-lbin_energy)/nbin_energy)*ienergy+lbin_energy) && ( (pfpenergy_recon[i][p]-energy_corr_intercept[p])/energy_corr_slope[p] < ((hbin_energy-lbin_energy)/nbin_energy)*(ienergy+1)+lbin_energy) ) {
                index_energy = ienergy;  
              }
            
            }
          }
         
          for (int j=0; j<nbin_longdist; j++) {
            hlong->Fill(hlong->GetBinCenter(j+1), LongProf[i][p][j]);
            if (index_energy != -1) hprofile[index_particle][index_energy][0]->Fill(hlong->GetBinCenter(j+1), LongProf[i][p][j]);
          }
          for (int j=0; j<nbin_trandist; j++) {
            htran->Fill(htran->GetBinCenter(j+1), TranProf[i][p][j]);
            if (index_energy != -1) {
              hprofile[index_particle][index_energy][1]->Fill(hprofile[index_particle][index_energy][1]->GetBinCenter(j+1), TranProf[i][p][j]);
              hprofile[index_particle][index_energy][2]->Fill(hprofile[index_particle][index_energy][2]->GetBinCenter(j+1), TranProf_1[i][p][j]);
              hprofile[index_particle][index_energy][3]->Fill(hprofile[index_particle][index_energy][3]->GetBinCenter(j+1), TranProf_2[i][p][j]);
              hprofile[index_particle][index_energy][4]->Fill(hprofile[index_particle][index_energy][4]->GetBinCenter(j+1), TranProf_3[i][p][j]);
              hprofile[index_particle][index_energy][5]->Fill(hprofile[index_particle][index_energy][5]->GetBinCenter(j+1), TranProf_4[i][p][j]);
              hprofile[index_particle][index_energy][6]->Fill(hprofile[index_particle][index_energy][6]->GetBinCenter(j+1), TranProf_5[i][p][j]);
            }
          }
        }
      } // loop over planes p
    } // loop over npfps i

  } // loop over nentries jentry


  // ------- energy correction ------
  std::vector<std::vector<double>> energy_recon;
  std::vector<std::vector<double>> energy_recon_error;
  std::vector<std::vector<double>> energy_true;
  std::vector<std::vector<double>> energy_true_error;
  energy_recon.clear();
  energy_recon_error.clear();
  energy_true.clear();
  energy_true_error.clear();

  TCanvas *c1_electron_energy[3];
  for (int p=0; p<3; p++) {
    std::vector<double> energy_recon_per_plane;
    std::vector<double> energy_recon_per_plane_error;
    std::vector<double> energy_true_per_plane;
    std::vector<double> energy_true_per_plane_error;
    energy_recon_per_plane.clear();
    energy_recon_per_plane_error.clear();
    energy_true_per_plane.clear();
    energy_true_per_plane_error.clear();

    c1_electron_energy[p] = new TCanvas(TString::Format("c1_electron_energy_%d",p), TString::Format("c1_electron_energy_%d",p), 1200, 900);
    c1_electron_energy[p]->Divide(4,3);

    for (int n=0; n<npoint_energy; n++) {
      c1_electron_energy[p]->cd(n+1);
      h1_energy_electron[p][n]->Draw();
      if (n > 0 && n < 7) {  
        double max_fit = h1_energy_electron[p][n]->GetBinCenter(h1_energy_electron[p][n]->GetMaximumBin()); 
        double left_fit = max_fit - 0.2*step_energy-n*0.05*step_energy;
        double right_fit = max_fit + 0.2*step_energy+n*0.05*step_energy;

        auto myfunc = new TF1("myfunc", "gaus", left_fit, right_fit);
        double par[3];
        h1_energy_electron[p][n]->Fit(myfunc, "", "", left_fit, right_fit);
        myfunc->GetParameters(&par[0]);
        cout << par[0] << " +/-  " << myfunc->GetParError(0) << endl;
        cout << par[1] << " +/-  " << myfunc->GetParError(1) << endl;
        cout << par[2] << " +/-  " << myfunc->GetParError(2) << endl;
        energy_recon_per_plane.push_back(par[1]);
        energy_recon_per_plane_error.push_back(myfunc->GetParError(1));
        energy_true_per_plane.push_back((n+0.5)*step_energy);
        energy_true_per_plane_error.push_back(0.5*step_energy);
      }
    }
    energy_recon.push_back(energy_recon_per_plane);
    energy_recon_error.push_back(energy_recon_per_plane_error);
    energy_true.push_back(energy_true_per_plane);
    energy_true_error.push_back(energy_true_per_plane_error);
    
    c1_electron_energy[p]->Write(TString::Format("hist_electron_energy_%d",p));
    c1_electron_energy[p]->SaveAs(TString::Format("hist_electron_energy_%d.pdf",p));
    delete c1_electron_energy[p];
  }

  // --------------------------------
  TCanvas *c1_grErr_energy_correction[3];
  TGraphErrors *grErr_energy_corr[3];
  for (int p = 0; p < 3; p++) {
    c1_grErr_energy_correction[p] = new TCanvas(TString::Format("c1_grErr_energy_correction_plane%d",p), "c1_grErr_energy_correction", 800,600);
    TPaveText *pt_energy_corr = new TPaveText(400.0,700.0,700.0,800.0);
    grErr_energy_corr[p] = new TGraphErrors(energy_recon[p].size(), &energy_true[p][0],&energy_recon[p][0], &energy_true_error[p][0], &energy_recon_error[p][0]);
    auto f1 = new TF1("f1", "[0]+[1]*x", 150.0, 1500.0);
    f1->SetParameters(0.0,1.0);
    double par_corr[2];
    grErr_energy_corr[p]->SetTitle(TString::Format("Shower energy correction (plane %d)", p));
    grErr_energy_corr[p]->GetXaxis()->SetTitle("True electron energy [MeV]");
    grErr_energy_corr[p]->GetYaxis()->SetTitle("Recon shower energy [MeV]");
    grErr_energy_corr[p]->Draw("AP");
    grErr_energy_corr[p]->Fit(f1,"","",150.0,1500.0);
    f1->GetParameters(&par_corr[0]);
    cout << par_corr[0] << endl;
    cout << par_corr[1] << endl;
    pt_energy_corr->AddText(TString::Format("Slope: %f",par_corr[1]));
    pt_energy_corr->AddText(TString::Format("Intercept: %f",par_corr[0]));
    pt_energy_corr->Draw();
    c1_grErr_energy_correction[p]->Write(TString::Format("grErr_energy_corr_plane_%d",p));
    c1_grErr_energy_correction[p]->SaveAs(TString::Format("grErr_energy_corr_plane_%d.pdf",p));

    delete c1_grErr_energy_correction[p];
  }
  

  // output
  output.Write();
  output.Close();
}

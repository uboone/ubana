#define ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void selection_ana::Loop()
{
  if (fChain == 0) return;

  // output
  TFile output("ouput_ana.root","recreate");
  ofstream myfile("events_for_scan.txt");

  // -------- for 3 planes with option whether plotting the induction planes ----------
  bool plotInductionPlanes = true;
  const int  plane[3] = {0, 1, 2};
  const char * shower_particle[6] = {"primary_e", "photon", "proton", "pion", "other_e", "all_others"}; // 0: primary electron shower; 1: photon shower; 2: proton shower; 3:pion shower; 4: other electron shower; 5: all the other shower
  const int color_options[3] = {2, 3, 4}; // red, green, blue
  const char * draw_options[3] = {"", "SAME", "SAME"};
  const char * eff_draw_options[4] = {"", "SAME", "SAME", "SAME"};
  const int eff_color_options[4] = {2, 3, 4, 6}; // red, green, blue
  const char * cut_options[4] = {"no cuts", "p>0.7; c>0.7", "p>0.8; c>0.8", "p>0.9; c>0.9"};

  double E_bins[21] = {0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 8000.0, 10000.0};

  //TH1D *h_Ee_num = new TH1D("h_Ee_num", "Electron Energy; E_{e} [MeV]; Shower Reconstruction Efficiency", 20, E_bins);
  //TH1D *h_Ee_num[4];
  TH1D *h_Ee_num[3][4];
  TH1D *h_Ee_den = new TH1D("h_Ee_den", "Electron Energy; E_{e} [MeV]; Shower Reconstruction Efficiency", 20, E_bins);
  h_Ee_den->Sumw2();
  for (int p=0; p<3; p++){
    for (int iplot = 0; iplot < 4; iplot++) {
      h_Ee_num[p][iplot] = new TH1D(TString::Format("h_Ee_num_plane%d_%d",p,iplot), "Efficiency; E_{e} [MeV]; #epsilon", 20, E_bins);
      h_Ee_num[p][iplot]->Sumw2();
    }
  }


  TEfficiency* h_Eff_Ee = new TEfficiency("h_Eff_Ee", "Electron Energy; E_{e} [MeV]; Shower Reconstruction Efficiency", 50, 0, 5000.0); 

  // correction for reco shower energy
  double step_energy = 200.0; // [MeV]
  int npoints_energy = 12;
  int nbins_energy = 10;
  TH1D *h_energy_shower_primary_e[12];
  TH1D *h_primary_e_shower_recoEnergy[3][12];
  for (int n = 0; n < npoints_energy; n++) {
    h_energy_shower_primary_e[n] = new TH1D(TString::Format("h_energy_shower_primary_e_%d",n), "Primary e shower energy", nbins_energy*(n+1), 0, (n+1)*step_energy);
  }
  for (int p = 0; p < 3; p++) {
    for (int n = 0; n < npoints_energy; n++) {
      h_primary_e_shower_recoEnergy[p][n] = new TH1D(TString::Format("h_primary_e_shower_recoEnergy_plane%d_%d",p,n), "Primary e shower energy", nbins_energy*(n+1), 0, (n+1)*step_energy);
    }
  }


  double energy_corr_slope[3] = {0.644342, 0.590877, 0.813455};
  double energy_corr_intercept[3] ={-24.662, 12.8708, -17.4905};

  // -------------------------------------------

  TH1D * h_shower_hasPrimary_e[3];

  TH1D *h_shower_numbers[3];
  TH1D *h_shower_purity[3][6];
  TH1D *h_shower_completeness[3][6];
  TH1D *h_shower_energy_resolution[3][6][4];
  TH2D *h2_shower_energy[3][6][4];

  TH1D *h_shower_vertex_resolution_x[3][6][4];
  TH1D *h_shower_vertex_resolution_y[3][6][4];
  TH1D *h_shower_vertex_resolution_z[3][6][4];

  TH1D *h_shower_direction_resolution[3][6][4];
  int nBins_shower_direction = 20;
  int lowBin_shower_direciton = -1.0;
  //int highBin_shower_direciton = 3.1415927;
  int highBin_shower_direciton = 1.0;

  int nPlots_shower_particles = 6;
  int nPlots_shower_purity_completeness_cuts = 4;

  double shower_cut_purity[4];
  double shower_cut_completeness[4];
  double shower_cutStep_purity = 0.1;
  double shower_cutStep_completeness = 0.1;
  double shower_cutLow_purity = 0.7;
  double shower_cutLow_completeness = 0.7;

  int nBins_energy_resolution = 50;
  double lowBin_energy_resolution = -1.5;
  double highBin_energy_resolution = 1.0;

  int nBins_shower_energyMC = 55;
  int nBins_shower_energyReco = 45;
  double binwidth_shower_energyMC = 100.0; // [MeV]
  double binwidth_shower_energyReco = 100.0; // [MeV]
  double lowBin_shower_energyMC = 0.0;
  double lowBin_shower_energyReco = 0.0;

  double highBin_shower_energyMC = nBins_shower_energyMC * binwidth_shower_energyMC + lowBin_shower_energyMC;
  double highBin_shower_energyReco = nBins_shower_energyReco * binwidth_shower_energyReco + lowBin_shower_energyReco;

  int nBins_vertex_resolution_x = 40;
  double lowBin_vertex_resolution_x = -10.0;
  double highBin_vertex_resolution_x =10.0;

  int nBins_vertex_resolution_y = 40;
  double lowBin_vertex_resolution_y = -10.0;
  double highBin_vertex_resolution_y = 10.0;

  int nBins_vertex_resolution_z = 40;
  double lowBin_vertex_resolution_z = -10.0;
  double highBin_vertex_resolution_z = 10.0;

  for (int i = 0; i < 3; i++) {
    h_shower_hasPrimary_e[i] = new TH1D(TString::Format("h_shower_hasPrimary_e_plane_%d",plane[i]), "h_shower_hasPrimary_e", 2, 0, 2);
    h_shower_hasPrimary_e[i]->GetXaxis()->SetTitle("shower from primary electron: 0 [N]; 1[Y]");

    h_shower_numbers[i] = new TH1D("h_shower_numbers", "Number of showers", 8, 0, 8);

    for (int j = 0; j < nPlots_shower_particles; j++) {
      h_shower_purity[i][j] = new TH1D(TString::Format("h_%s_shower_purity_plane_%d", shower_particle[j], plane[i]), TString::Format("h_%s_shower_purity", shower_particle[j]), 110, 0.0, 1.1);
      h_shower_completeness[i][j] = new TH1D(TString::Format("h_%s_shower_completeness_plane_%d", shower_particle[j], plane[i]), TString::Format("h_%s_shower_completeness", shower_particle[j]), 110, 0.0, 1.1);

      for (int k = 0; k < nPlots_shower_purity_completeness_cuts; k++) {
        if (k == 0) {
          shower_cut_purity[k] = 0;
          shower_cut_completeness[k] = 0;
          h_shower_energy_resolution[i][j][k] = new TH1D(TString::Format("h_%s_shower_energy_plane_%d_%d",shower_particle[j], plane[i], k), TString::Format("%s: no cuts on purity and completeness", shower_particle[j]), nBins_energy_resolution, lowBin_energy_resolution, highBin_energy_resolution);

          h2_shower_energy[i][j][k] = new TH2D(TString::Format("h_%s_shower_energy_plane_%d_%d", shower_particle[j], plane[i], k), TString::Format("%s: no cuts on purity and completeness; E_{MC} [MeV]; E_{Reco} [MeV]", shower_particle[j]), nBins_shower_energyMC, lowBin_shower_energyMC, highBin_shower_energyMC, nBins_shower_energyReco, lowBin_shower_energyReco, highBin_shower_energyReco);

          h_shower_vertex_resolution_x[i][j][k] = new TH1D(TString::Format("h_%s_shower_vertex_x_%d_%d",shower_particle[j], plane[i], k), TString::Format("%s: no cuts on purity and completeness", shower_particle[j]), nBins_vertex_resolution_x, lowBin_vertex_resolution_x, highBin_vertex_resolution_x);
          h_shower_vertex_resolution_y[i][j][k] = new TH1D(TString::Format("h_%s_shower_vertex_y_%d_%d",shower_particle[j], plane[i], k), TString::Format("%s: no cuts on purity and completeness", shower_particle[j]), nBins_vertex_resolution_y, lowBin_vertex_resolution_y, highBin_vertex_resolution_y);
          h_shower_vertex_resolution_z[i][j][k] = new TH1D(TString::Format("h_%s_shower_vertex_z_%d_%d",shower_particle[j], plane[i], k), TString::Format("%s: no cuts on purity and completeness", shower_particle[j]), nBins_vertex_resolution_z, lowBin_vertex_resolution_z, highBin_vertex_resolution_z);

          h_shower_direction_resolution[i][j][k] = new TH1D(TString::Format("h_%s_shower_direction_%d_%d",shower_particle[j], plane[i], k), TString::Format("%s: no cuts on purity and completeness", shower_particle[j]), nBins_shower_direction, lowBin_shower_direciton, highBin_shower_direciton);
        }
        else {
          shower_cut_purity[k] = shower_cutLow_purity + (k - 1) * shower_cutStep_purity;
          shower_cut_completeness[k] = shower_cutLow_completeness + (k - 1) * shower_cutStep_completeness;
          h_shower_energy_resolution[i][j][k] = new TH1D(TString::Format("h_%s_shower_energy_plane_%d_%d",shower_particle[j],plane[i], k), TString::Format("%s: purity > %.1f completeness > %.1f", shower_particle[j], shower_cut_purity[k], shower_cut_completeness[k]), nBins_energy_resolution, lowBin_energy_resolution, highBin_energy_resolution);

          h2_shower_energy[i][j][k] = new TH2D(TString::Format("h_%s_shower_energy_plane_%d_%d", shower_particle[j], plane[i], k), TString::Format("%s: purity > %.1f completeness > %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", shower_particle[j], shower_cut_purity[k], shower_cut_completeness[k]), nBins_shower_energyMC, lowBin_shower_energyMC, highBin_shower_energyMC, nBins_shower_energyReco, lowBin_shower_energyReco, highBin_shower_energyReco);

          h_shower_vertex_resolution_x[i][j][k] = new TH1D(TString::Format("h_%s_shower_vertex_x_%d_%d",shower_particle[j], plane[i], k),TString::Format("%s: purity > %.1f completeness > %.1f", shower_particle[j], shower_cut_purity[k], shower_cut_completeness[k]), nBins_vertex_resolution_x, lowBin_vertex_resolution_x, highBin_vertex_resolution_x);
          h_shower_vertex_resolution_y[i][j][k] = new TH1D(TString::Format("h_%s_shower_vertex_y_%d_%d",shower_particle[j], plane[i], k),TString::Format("%s: purity > %.1f completeness > %.1f", shower_particle[j], shower_cut_purity[k], shower_cut_completeness[k]), nBins_vertex_resolution_y, lowBin_vertex_resolution_y, highBin_vertex_resolution_y);
          h_shower_vertex_resolution_z[i][j][k] = new TH1D(TString::Format("h_%s_shower_vertex_z_%d_%d",shower_particle[j], plane[i], k),TString::Format("%s: purity > %.1f completeness > %.1f", shower_particle[j], shower_cut_purity[k], shower_cut_completeness[k]), nBins_vertex_resolution_z, lowBin_vertex_resolution_z, highBin_vertex_resolution_z);

          h_shower_direction_resolution[i][j][k] = new TH1D(TString::Format("h_%s_shower_direction_%d_%d",shower_particle[j], plane[i], k), TString::Format("%s: purity > %.1f completeness > %.1f", shower_particle[j], shower_cut_purity[k], shower_cut_completeness[k]), nBins_shower_direction, lowBin_shower_direciton, highBin_shower_direciton);
        }

        h_shower_energy_resolution[i][j][k]->GetXaxis()->SetTitle("shower energy resolution (#frac{E_{Reco}-E_{MC}}{E_{MC}})");
        h_shower_vertex_resolution_x[i][j][k]->GetXaxis()->SetTitle("shower vertex x resolution (x_{Reco}-x_{MC}) [cm]");
        h_shower_vertex_resolution_y[i][j][k]->GetXaxis()->SetTitle("shower vertex y resolution (y_{Reco}-y_{MC}) [cm]");
        h_shower_vertex_resolution_z[i][j][k]->GetXaxis()->SetTitle("shower vertex z resolution (z_{Reco}-z_{MC}) [cm]");

        h_shower_direction_resolution[i][j][k]->GetXaxis()->SetTitle("cos#theta");

      } // k
    } // j

  } // i

  TH1D *h_numOfShowers = new TH1D("h_numOfShowers", "h_numOfShowers", 50, 0, 50);
  //TH1D *h_num_showers = new TH1D("h_num_showers","h_num_showers", 8, 0, 8); // 1: primary electron shower; 2: photon shower; 3: proton shower; 4:pion shower; 5: other electron shower; 6: all the other shower

  TH1D *h_sh_purity = new TH1D("h_sh_purity", "h_sh_purity", 110, 0.0, 1.1);
  TH1D *h_sh_purity_electron = new TH1D("h_sh_purity_electron", "h_sh_purity_electron", 110, 0.0, 1.1);
  TH1D *h_sh_purity_photon = new TH1D("h_sh_purity_photon", "h_sh_purity_photon", 110, 0.0, 1.1);
  TH1D *h_sh_purity_proton = new TH1D("h_sh_purity_proton", "h_sh_purity_proton", 110, 0.0, 1.1);
  TH1D *h_sh_purity_pion = new TH1D("h_sh_purity_pion", "h_sh_purity_pion", 110, 0.0, 1.1);
  TH1D *h_sh_purity_electron_other = new TH1D("h_sh_purity_electron_other", "h_sh_purity_electron_other", 110, 0.0, 1.1);
  TH1D *h_sh_purity_all_other = new TH1D("h_sh_purity_all_other", "h_sh_purity_all_other", 110, 0.0, 1.1);

  TH1D *h_sh_completeness = new TH1D("h_sh_completeness", "h_sh_completeness", 110, 0.0, 1.1);
  TH1D *h_sh_completeness_electron = new TH1D("h_sh_completeness_electron", "h_sh_completeness_electron", 110, 0.0, 1.1);
  TH1D *h_sh_completeness_photon = new TH1D("h_sh_completeness_photon", "h_sh_completeness_photon", 110, 0.0, 1.1);
  TH1D *h_sh_completeness_proton = new TH1D("h_sh_completeness_proton", "h_sh_completeness_proton", 110, 0.0, 1.1);
  TH1D *h_sh_completeness_pion = new TH1D("h_sh_completeness_pion", "h_sh_completeness_pion", 110, 0.0, 1.1);
  TH1D *h_sh_completeness_electron_other = new TH1D("h_sh_completeness_electron_other", "h_sh_completeness_electron_other", 110, 0.0, 1.1);
  TH1D *h_sh_completeness_all_other = new TH1D("h_sh_completeness_all_other", "h_sh_completeness_all_other", 110, 0.0, 1.1);

  int nPlots = 4;
  double cutStep_purity = 0.1;
  double cutStep_completeness = 0.1;
  double cutLow_purity = 0.7;
  double cutLow_completenss = 0.7;

  // -------------- energy resolution -------------
  TH1D *h_sh_E_resoluton = new TH1D("h_sh_E_resoluton", "h_sh_E_resoluton",90, -35500, 4500); // only the diff: Reco - MC for all showers
  TH1D* h_sh_energy_resolution_electron[4];
  TH1D* h_sh_energy_resolution_photon[4];
  TH1D* h_sh_energy_resolution_proton[4];
  TH1D* h_sh_energy_resolution_pion[4];

  int nBins_energyResolution = 50;
  double lowBin_energyResolution = -1.5;
  double highBin_energyResolution = 1.0;

  TH2D* hist_sh_energy_electron[4];
  TH2D* hist_sh_energy_photon[4];
  TH2D* hist_sh_energy_proton[4];
  TH2D* hist_sh_energy_pion[4];
  TH2D* hist_sh_energy_electron_other[4];
  int nBins_energyMC = 55;
  int nBins_energyReco = 45;
  double binwidth_energyMC = 100.0;
  double binwidth_energyReco = 100.0;
  double highBin_energyMC = nBins_energyMC*binwidth_energyMC;
  double highBin_energyReco = nBins_energyReco*binwidth_energyReco;

  int cout_electron_with_purity_completeless[4];
  double cut_purity[4];
  double cut_completeness[4];
  for (int i = 0; i < nPlots; ++i) {
    cout_electron_with_purity_completeless[i] = 0;

    if (i == 0) {
      //h_Ee_num[i] = new TH1D(TString::Format("h_Ee_num_%d",i), "Efficiency; E_{e} [MeV]; #epsilon ", 20, E_bins);
      //h_Ee_num[i]->Sumw2();

      h_sh_energy_resolution_electron[i] = new TH1D(TString::Format("h_electron_%d",i),"electron: no cuts on purity and completeness", nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);
      h_sh_energy_resolution_photon[i] = new TH1D(TString::Format("h_photon_%d",i),"photon: no cuts on purity and completeness", nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);
      h_sh_energy_resolution_proton[i] = new TH1D(TString::Format("h_proton_%d",i),"proton: no cuts on purity and completeness", nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);
      h_sh_energy_resolution_pion[i] = new TH1D(TString::Format("h_pion_%d",i),"pion: no cuts on purity and completeness", nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);

      hist_sh_energy_electron[i] = new TH2D(TString::Format("hist_electron_%d",i),"electron: no cuts on purity and completeness; E_{MC} [MeV]; E_{Reco} [MeV]", nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_photon[i] = new TH2D(TString::Format("hist_photon_%d",i),"photon: no cuts on purity and completeness; E_{MC} [MeV]; E_{Reco} [MeV]", nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_proton[i] = new TH2D(TString::Format("hist_proton_%d",i),"proton: no cuts on purity and completeness; E_{MC} [MeV]; E_{Reco} [MeV]", nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_pion[i] = new TH2D(TString::Format("hist_pion_%d",i),"pion: no cuts on purity and completeness; E_{MC} [MeV]; E_{Reco} [MeV]", nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_electron_other[i] = new TH2D(TString::Format("hist_electron_other_%d",i),"electron (other): no cuts on purity and completeness; E_{MC} [MeV]; E_{Reco} [MeV]", nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      cut_purity[i] = 0;
      cut_completeness[i] = 0;
    }
    else {
      cut_purity[i] = cutLow_purity+(i-1)*cutStep_purity;
      cut_completeness[i] = cutLow_completenss+(i-1)*cutStep_completeness;

      //h_Ee_num[i] = new TH1D(TString::Format("h_Ee_num_%d",i), "Efficiency; E_{e} [MeV]; #epsilon", 20, E_bins);
      //h_Ee_num[i]->Sumw2();

      h_sh_energy_resolution_electron[i] = new TH1D(TString::Format("h_electron_%d",i),TString::Format("electron: purity > %.1f completeness > %.1f", cut_purity[i], cut_completeness[i]), nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);
      h_sh_energy_resolution_photon[i] = new TH1D(TString::Format("h_photon_%d",i),TString::Format("photon: purity > %.1f completeness > %.1f", cut_purity[i], cut_completeness[i]), nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);
      h_sh_energy_resolution_proton[i] = new TH1D(TString::Format("h_proton_%d",i),TString::Format("proton: purity > %.1f completeness > %.1f", cut_purity[i], cut_completeness[i]), nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);
      h_sh_energy_resolution_pion[i] = new TH1D(TString::Format("h_pion_%d",i),TString::Format("pion: purity > %.1f completeness > %.1f", cut_purity[i], cut_completeness[i]), nBins_energyResolution, lowBin_energyResolution, highBin_energyResolution);

      hist_sh_energy_electron[i] = new TH2D(TString::Format("hist_electron_%d",i),TString::Format("electron: purity > %.1f completeness > %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", cut_purity[i], cut_completeness[i]), nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_photon[i] = new TH2D(TString::Format("hist_photon_%d",i),TString::Format("photon: purity > %.1f completeness > %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", cut_purity[i], cut_completeness[i]), nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_proton[i] = new TH2D(TString::Format("hist_proton_%d",i),TString::Format("proton: purity > %.1f completeness > %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", cut_purity[i], cut_completeness[i]), nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_pion[i] = new TH2D(TString::Format("hist_pion_%d",i),TString::Format("pion: purity > %.1f completeness > %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", cut_purity[i], cut_completeness[i]), nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);
      hist_sh_energy_electron_other[i] = new TH2D(TString::Format("hist_electron_other_%d",i),TString::Format("electron (other): purity > %.1f completeness > %.1f; E_{MC} [MeV]; E_{Reco} [MeV]", cut_purity[i], cut_completeness[i]), nBins_energyMC, 0, highBin_energyMC, nBins_energyReco, 0, highBin_energyReco);

    }
    cout << cut_purity[i] << endl;
    cout << cut_completeness[i] << endl;
  }



  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries: " << nentries << endl;


  int sh_count = 0;
  int sh_count_primary_electron[3] = {0,0,0}; // 1: primary electron shower
  int sh_count_photon[3] = {0,0,0};  // 2: photon shower
  int sh_count_proton[3] = {0,0,0}; // 3: proton shower
  int sh_count_pion[3] = {0,0,0}; // 4: pion shower
  int sh_count_other_electron[3] = {0,0,0}; // 5: other electron shower
  int sh_count_all_other[3] = {0,0,0}; // 6: shower other than 1-5

  int test = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    h_numOfShowers->Fill(numOfShowers); // this contains all showers; 

    int shower_id = 0; // 0: undefined; 1: primary electron shower; 2: photon shower; 3: proton shower; 4:pion shower; 5: other electron shower; 6: all the other shower
    std::vector<int> showerID = {0,0,0}; // for three planes; 0: undefined; 1: primary electron shower; 2: photon shower; 3: proton shower; 4:pion shower; 5: other electron shower; 6: all the other shower

    //int sh_primary_e_count_per_event[4] = {0, 0, 0, 0}; 
    int shower_primary_e_count_per_event[3][4];
    for (int p = 0; p < 3; p++) {
      for (int iplot = 0; iplot < nPlots; iplot++) {
        shower_primary_e_count_per_event[p][iplot] = 0;
      }
    }

    // here, we look at showers in a slice if it exists
    for (int k = 0; k < sh_start_X->size(); k++) {
      
      // all showers
      sh_count += 1;
      //cout << showerEnergyReco->size() << endl;   
      
      for (int p = 0; p < 3; p++) {

        if ( (*showerParticlePDG)[k][p] == -99999 || (*sh_purity)[k][p] <= 0.0 || (*sh_completeness)[k][p] <= 0.0 || (*showerEnergyReco)[k][p] < 1.0e-20) {
          h_shower_numbers[p]->Fill(6); // label the showers that cannot be associated with a MC particle as all other showers (log scale)
          continue; // skip showers that cannot be associated with a MC particle.
        }
        //if ((*showerEnergyReco)[k][p] <=0) {
        //  cout << ".........www.........." << endl;
        //  cout << "(*showerEnergyReco)[k][p]: " << (*showerEnergyReco)[k][p] << endl;
        //}
        //else 
        if ( (((*showerEnergyReco)[k][p] - (*showerEnergyMC)[k][p]) / (*showerEnergyMC)[k][p]) < -0.95 ) {
          myfile << "event: " << event << " \n";
          myfile << "run: " << run << "\n";
          myfile << "subrun: " << subrun << "\n";
          myfile << "plane " << p <<  "\n";
          myfile << "(*slice_id)[k]: " << (*slice_id)[k] <<  "\n";
          myfile << "(*sh_start_X)[k]: " << (*sh_start_X)[k] << "\n";
          myfile << "(*sh_start_Y)[k]: " << (*sh_start_Y)[k] << "\n";
          myfile << "(*sh_start_Z)[k]: " << (*sh_start_Z)[k] << "\n";
          myfile << "(*showerParticlePDG)[k][0]: " << (*showerParticlePDG)[k][0] << "\n";
          myfile << "(*showerParticlePDG)[k][1]: " << (*showerParticlePDG)[k][1] << "\n";
          myfile << "(*showerParticlePDG)[k][2]: " << (*showerParticlePDG)[k][2] << "\n";
          myfile << "(*showerEnergyReco)[k][0]: " << (*showerEnergyReco)[k][0] << "\n";
          myfile << "(*showerEnergyReco)[k][1]: " << (*showerEnergyReco)[k][1] << "\n";
          myfile << "(*showerEnergyReco)[k][2]: " << (*showerEnergyReco)[k][2] << "\n";
          myfile << "(*showerEnergyMC)[k][0]: " << (*showerEnergyMC)[k][0] << "\n";
          myfile << "(*showerEnergyMC)[k][1]: " << (*showerEnergyMC)[k][1] << "\n";
          myfile << "(*showerEnergyMC)[k][2]: " << (*showerEnergyMC)[k][2] << "\n";
        }

        h_shower_hasPrimary_e[p]->Fill( (*sh_hasPrimary_e)[k][p] );

        // primary electron showers
        if ( (*sh_hasPrimary_e)[k][p] == 1 ) {
          shower_id = 1;
          sh_count_primary_electron[p] += 1;
        }
        else {
          // photon showers
          if ( abs((*showerParticlePDG)[k][p]) == 22 ) {
            shower_id = 2;
            sh_count_photon[p] += 1;
          }
          // proton showers
          else if ( abs((*showerParticlePDG)[k][p]) == 2212 ) {
            shower_id = 3;
            sh_count_proton[p] += 1;
          }
          // pion showers
          else if ( abs((*showerParticlePDG)[k][p]) == 211 ) {
            shower_id = 4;
            sh_count_pion[p] += 1;
          }
          // other electron showers
          else if ( abs((*showerParticlePDG)[k][p]) == 11 ) {
            shower_id = 5;
            sh_count_other_electron[p] += 1;
          }
          // all the other showers
          else {
            shower_id = 6;
            sh_count_all_other[p] += 1;
          } 
        }  // if ... else ... for different showers

        h_shower_numbers[p]->Fill(shower_id);
        h_shower_purity[p][shower_id-1]->Fill((*sh_purity)[k][p] );
        h_shower_completeness[p][shower_id-1]->Fill((*sh_completeness)[k][p] );

        for (int iplot = 0; iplot < nPlots; ++iplot) {
          if ( (*sh_purity)[k][p] >= cut_purity[iplot] && (*sh_completeness)[k][p] >= cut_completeness[iplot] ) {
            //h_shower_energy_resolution[p][shower_id-1][iplot]->Fill( ((*showerEnergyReco)[k][p] - (*showerEnergyMC)[k][p]) / (*showerEnergyMC)[k][p] ); 
            h_shower_energy_resolution[p][shower_id-1][iplot]->Fill( ( ((*showerEnergyReco)[k][p] - energy_corr_intercept[p]) / energy_corr_slope[p] - (*showerEnergyMC)[k][p]) / (*showerEnergyMC)[k][p] ); 
            //h2_shower_energy[p][shower_id-1][iplot]->Fill( (*showerEnergyMC)[k][p], (*showerEnergyReco)[k][p]) ;
            h2_shower_energy[p][shower_id-1][iplot]->Fill( (*showerEnergyMC)[k][p],((*showerEnergyReco)[k][p] - energy_corr_intercept[p]) / energy_corr_slope[p]) ;
            h_shower_vertex_resolution_x[p][shower_id-1][iplot]->Fill( (*sh_start_X)[k] - (*sh_start_X_MC)[k][p] );
            h_shower_vertex_resolution_y[p][shower_id-1][iplot]->Fill( (*sh_start_Y)[k] - (*sh_start_Y_MC)[k][p] );
            h_shower_vertex_resolution_z[p][shower_id-1][iplot]->Fill( (*sh_start_Z)[k] - (*sh_start_Z_MC)[k][p] );

            TVector3 shReconDirection((*sh_direction_X)[k], (*sh_direction_Y)[k], (*sh_direction_Z)[k]);
            TVector3 shMCDirection((*sh_direction_X_MC)[k][p], (*sh_direction_Y_MC)[k][p], (*sh_direction_Z_MC)[k][p]);
            h_shower_direction_resolution[p][shower_id-1][iplot]->Fill(cos(shMCDirection.Angle(shReconDirection)));

            if (shower_id == 1 &&  p == 2) {
              //sh_primary_e_count_per_event[iplot] += 1;

              // for correction of reco shower energy; no cut on purity and completeness
              // here, we use the primary electron shower on the collection plane
              if (iplot == 0) {
                for (int n1 = 0; n1 < npoints_energy; n1++) {
                  if ( (leptonEnergy_truth * 1000.0 > n1*step_energy) && (leptonEnergy_truth * 1000.0 <= (n1+1)*step_energy)) {
                    h_energy_shower_primary_e[n1]->Fill((*showerEnergyReco)[k][p]);
                  }
                }
              }
            } // end of if (shower_id == 1 &&  p == 2)
            
            // count primary e shower for each plane with/without cuts
            if (shower_id == 1) {
              shower_primary_e_count_per_event[p][iplot] += 1;
              if (iplot == 0) {
                for (int n1 = 0; n1 < npoints_energy; n1++) {
                  if ( (leptonEnergy_truth * 1000.0 > n1*step_energy) && (leptonEnergy_truth * 1000.0 <= (n1+1)*step_energy)) { 
                    h_primary_e_shower_recoEnergy[p][n1]->Fill((*showerEnergyReco)[k][p]);
                  }
                }
              }
            } // end of if (shower_id == 1) 

          }
        }


      } // end of loop p


    } // end of loop sh_start_X->size() k


    for (int p = 0; p < 3; p++) {
      for (int iplot = 0; iplot < nPlots; ++iplot) {
        if (shower_primary_e_count_per_event[p][iplot] > 0) h_Ee_num[p][iplot]->Fill(leptonEnergy_truth * 1000.0); // [MeV]
      }
    }

    if (ccnc_truth == 0 && abs(leptonPDG_truth) == 11) {
      h_Ee_den->Fill(leptonEnergy_truth * 1000.0); // [MeV] 
    }

  } // end of loop nentries jentry

  cout << "sh_count: " << sh_count << endl;
  for (int p = 0; p < 3; p++) {
    cout << "total shower on plane " << p << ": " << sh_count_primary_electron[p] + sh_count_photon[p] + sh_count_proton[p] + sh_count_pion[p] + sh_count_other_electron[p] + sh_count_all_other[p] << endl;
    cout << "sh_count_primary_electron: " << sh_count_primary_electron[p] << endl;
    cout << "sh_count_photon: " << sh_count_photon[p] << endl;
    cout << "sh_count_proton: " << sh_count_proton[p] << endl;
    cout << "sh_count_pion: " << sh_count_pion[p] << endl;
    cout << "sh_count_other_electron: " << sh_count_other_electron[p] << endl;
    cout << "sh_count_all_other: " << sh_count_all_other[p] << endl;

  }

  // --------------------
  TCanvas *c1_hasPrimary_e = new TCanvas("c1_hasPrimary_e","c1_hasPrimary_e",800,600);
  auto lg_hasPrimary_e = new TLegend(0.7,0.7,0.9,0.9);
  for (int p = 0; p < 3; p++) {
    h_shower_hasPrimary_e[p]->Draw(draw_options[p]);
    if (p == 0) h_shower_hasPrimary_e[p]->SetStats(0);
    h_shower_hasPrimary_e[p]->SetLineColor(color_options[p]);
    lg_hasPrimary_e->AddEntry(h_shower_hasPrimary_e[p], TString::Format("plane %d",plane[p]),"l");
  }
  lg_hasPrimary_e->Draw();
  c1_hasPrimary_e->Write("h_shower_hasPrimary_e");
  c1_hasPrimary_e->SaveAs("h_shower_hasPrimary_e.pdf");

  // --------------------
  double maxVal_shower_numbers = -999.0;
  for (int p = 0; p < 3; p++) {
    //cout << "........" << h_shower_numbers[p]->GetMaximum() << endl;
    if (h_shower_numbers[p]->GetMaximum() > maxVal_shower_numbers) {
      maxVal_shower_numbers = h_shower_numbers[p]->GetMaximum();
    }
  }

  TCanvas *c1_shower_numbers = new TCanvas("c1_shower_numbers","c1_shower_numbers",800,600);
  // --- normal y scale
  //auto lg_shower_numbers = new TLegend(0.7,0.7,0.9,0.9);
  //TPaveText *pt_shower_numbers = new TPaveText(3.5,0.5*maxVal_shower_numbers,5.5,0.9*maxVal_shower_numbers); // (x1,y1,x2,y2): coordinate, not relative to the canvas, use the x,y values
  // --- log y scale
  auto lg_shower_numbers = new TLegend(0.2,0.75,0.4,0.9);
  TPaveText *pt_shower_numbers = new TPaveText(3.5,0.05*maxVal_shower_numbers,5.5,0.9*maxVal_shower_numbers); // (x1,y1,x2,y2): coordinate, not relative to the canvas, use the x,y values
  for (int p = 0; p < 3; p++) {
    h_shower_numbers[p]->Draw(draw_options[p]);
    gPad->SetLogy(); // for log y scale
    if (p == 0) {
      h_shower_numbers[p]->SetStats(0);
      h_shower_numbers[p]->SetMaximum(1.05*maxVal_shower_numbers);
      pt_shower_numbers->AddText("1: primary electron showers");
      pt_shower_numbers->AddText("2: photon showers");
      pt_shower_numbers->AddText("3: proton showers");
      pt_shower_numbers->AddText("4: pion showers");
      pt_shower_numbers->AddText("5: other electron showers");
      pt_shower_numbers->AddText("6: all the other showers");
      pt_shower_numbers->Draw();
    }
    h_shower_numbers[p]->SetLineColor(color_options[p]);
    lg_shower_numbers->AddEntry(h_shower_numbers[p], TString::Format("plane %d",plane[p]),"l");
  }
  lg_shower_numbers->Draw();
  c1_shower_numbers->Write("h_shower_numbers");
  c1_shower_numbers->SaveAs("h_shower_numbers.pdf");

  // -------------------------------------------------------
  double maxVal_shower_purity[6];
  double maxVal_shower_completeness[6];
  for (int i = 0; i < 6; i++) {
    maxVal_shower_purity[i] = -999.0;
    maxVal_shower_completeness[i] = -999.0;
    for (int p = 0; p < 3; p++) {
      if (h_shower_purity[p][i]->GetMaximum() > maxVal_shower_purity[i]) {
        maxVal_shower_purity[i] = h_shower_purity[p][i]->GetMaximum();
      }
      if (h_shower_completeness[p][i]->GetMaximum() > maxVal_shower_completeness[i]) {
        maxVal_shower_completeness[i] = h_shower_completeness[p][i]->GetMaximum();
      }
    }
  }

  // -------------------------------
  TCanvas *c1_shower_purity = new TCanvas("c1_shower_purity", "c1_shower_purity", 1200, 800);
  c1_shower_purity->Divide(3,2);
  for (int i = 0; i < 6; i++) {
    auto lg_shower_purity = new TLegend(0.4,0.6,0.6,0.8);
    c1_shower_purity->cd(i+1);
    for (int p = 0; p < 3; p++) {
      h_shower_purity[p][i]->Draw(draw_options[p]);
      if (p == 0) {
        h_shower_purity[p][i]->SetStats(0);
        h_shower_purity[p][i]->SetMaximum(1.05*maxVal_shower_purity[i]);
      }
      h_shower_purity[p][i]->SetLineColor(color_options[p]);
      lg_shower_purity->AddEntry(h_shower_purity[p][i], TString::Format("plane %d", p), "l");
    }
    lg_shower_purity->Draw();
  }
  c1_shower_purity->Write("h_shower_purity");
  c1_shower_purity->SaveAs("h_shower_purity.pdf");

  // ----------------------------
  TCanvas *c1_shower_completeness = new TCanvas("c1_shower_completeness", "c1_shower_completeness", 1200, 800);
  c1_shower_completeness->Divide(3,2);
  for (int i = 0; i < 6; i++) {
    auto lg_shower_completeness = new TLegend(0.4,0.6,0.6,0.8);
    c1_shower_completeness->cd(i+1);
    for (int p = 0; p < 3; p++) {
      h_shower_completeness[p][i]->Draw(draw_options[p]);
      if (p == 0) {
        h_shower_completeness[p][i]->SetStats(0);
        h_shower_completeness[p][i]->SetMaximum(1.05*maxVal_shower_completeness[i]);
      }
      h_shower_completeness[p][i]->SetLineColor(color_options[p]);
      lg_shower_completeness->AddEntry(h_shower_completeness[p][i], TString::Format("plane %d", p), "l");
    }
    lg_shower_completeness->Draw();
  }
  c1_shower_completeness->Write("h_shower_completeness");
  c1_shower_completeness->SaveAs("h_shower_completeness.pdf");

  // ----------------------------------------------------
  double maxVal_shower_energy_resolution[6][4];
  double maxVal_shower_vertex_resolution_x[6][4];
  double maxVal_shower_vertex_resolution_y[6][4];
  double maxVal_shower_vertex_resolution_z[6][4];
  double maxVal_shower_direction_resolution[6][4];
  for (int j = 0; j < 6; ++j) {
    for (int iplot = 0; iplot < nPlots; ++iplot) {
      maxVal_shower_energy_resolution[j][iplot] = -999.0;
      maxVal_shower_vertex_resolution_x[j][iplot] = -999.0;
      maxVal_shower_vertex_resolution_y[j][iplot] = -999.0;
      maxVal_shower_vertex_resolution_z[j][iplot] = -999.0;
      maxVal_shower_direction_resolution[j][iplot] = -999.0;
      for (int p = 0; p < 3; ++p){
        if (h_shower_energy_resolution[p][j][iplot]->GetMaximum() > maxVal_shower_energy_resolution[j][iplot]) {
          maxVal_shower_energy_resolution[j][iplot] = h_shower_energy_resolution[p][j][iplot]->GetMaximum();
        }
        if (h_shower_vertex_resolution_x[p][j][iplot]->GetMaximum() > maxVal_shower_vertex_resolution_x[j][iplot]) {
          maxVal_shower_vertex_resolution_x[j][iplot] = h_shower_vertex_resolution_x[p][j][iplot]->GetMaximum();
        }
        if (h_shower_vertex_resolution_y[p][j][iplot]->GetMaximum() > maxVal_shower_vertex_resolution_y[j][iplot]) {
          maxVal_shower_vertex_resolution_y[j][iplot] = h_shower_vertex_resolution_y[p][j][iplot]->GetMaximum();
        }
        if (h_shower_vertex_resolution_z[p][j][iplot]->GetMaximum() > maxVal_shower_vertex_resolution_z[j][iplot]) {
          maxVal_shower_vertex_resolution_z[j][iplot] = h_shower_vertex_resolution_z[p][j][iplot]->GetMaximum();
        }
        if (h_shower_direction_resolution[p][j][iplot]->GetMaximum() > maxVal_shower_direction_resolution[j][iplot]) {
          maxVal_shower_direction_resolution[j][iplot] = h_shower_direction_resolution[p][j][iplot]->GetMaximum();
        }
      }
    }
  }
  // ----------------------------
  TCanvas *c1_shower_energy_resolution[6];
  for (int j = 0; j < 6; ++j) {
    c1_shower_energy_resolution[j] = new TCanvas(TString::Format("c1_shower_energy_resolution_%s",shower_particle[j]), TString::Format("c1_shower_energy_resolution_%s",shower_particle[j]), 1000, 800);
    c1_shower_energy_resolution[j]->Divide(2,2);
    for (int iplot = 0; iplot < nPlots; ++iplot){
      auto lg_shower_energy_resolution = new TLegend(0.7,0.7,0.9,0.9);
      c1_shower_energy_resolution[j]->cd(iplot+1);
      for (int p = 0; p < 3; p++) {
        h_shower_energy_resolution[p][j][iplot]->Draw(draw_options[p]);
        if (p==0) {
          h_shower_energy_resolution[p][j][iplot]->SetStats(0);
          h_shower_energy_resolution[p][j][iplot]->SetMaximum(1.05*maxVal_shower_energy_resolution[j][iplot]);
        }
        h_shower_energy_resolution[p][j][iplot]->SetLineColor(color_options[p]);
        lg_shower_energy_resolution->AddEntry(h_shower_energy_resolution[p][j][iplot], TString::Format("plane %d", p), "l");
      }
      lg_shower_energy_resolution->Draw();
    }
    c1_shower_energy_resolution[j]->Write(TString::Format("h_shower_energy_resolution_%s", shower_particle[j]));
    c1_shower_energy_resolution[j]->SaveAs(TString::Format("h_shower_energy_resolution_%s.pdf", shower_particle[j]));
  }

  // --------------------
  TCanvas *c1_shower_energy[18];
  for (int p = 0; p < 3; p++) {
    for (int j = 0; j < 6; ++j) {
      c1_shower_energy[j+6*p] = new TCanvas(TString::Format("c1_shower_energy_%s_%d",shower_particle[j], plane[p]), TString::Format("c1_shower_energy_%s_%d",shower_particle[j],plane[p]), 1000, 800);
      c1_shower_energy[j+6*p]->Divide(2,2);
      for (int iplot = 0; iplot < nPlots; ++iplot){
        c1_shower_energy[j+6*p]->cd(iplot+1);
        h2_shower_energy[p][j][iplot]->Draw("COLZ");
        gStyle->SetStatY(0.9);// Set y-position (fraction of pad size)
        gStyle->SetStatX(0.45);// Set x-position (fraction of pad size)
        gStyle->SetStatW(0.3);// Set width of stat-box (fraction of pad size)
        gStyle->SetStatH(0.2);// Set height of stat-box (fraction of pad size)
        //TLine myline(0,0,4000.0,4000.0);
        auto myline = new TLine(0.0,0.0,4500.0,4500.0);
        myline->SetLineColor(kRed);
        myline->SetLineWidth(3);
        myline->Draw();
      }
      c1_shower_energy[j+6*p]->Write(TString::Format("h2_shower_energy_%s_%d",shower_particle[j], plane[p]));
      c1_shower_energy[j+6*p]->SaveAs(TString::Format("h2_shower_energy_%s_%d.pdf",shower_particle[j], plane[p]));
    }
  }

  // -----------------------------------------------------
  // -------------------------
  TCanvas *c1_shower_vertex_resolution[6];
  for (int j = 0; j < 6; ++j) {
    c1_shower_vertex_resolution[j] = new TCanvas(TString::Format("c1_shower_vertex_resolution_%s",shower_particle[j]), TString::Format("c1_shower_vertex_resolution_%s",shower_particle[j]), 1200, 900);
    c1_shower_vertex_resolution[j]->Divide(4,3);
    for (int i = 0; i < 3; i++) {
      for (int iplot = 0; iplot < nPlots; ++iplot){
        auto lg_shower_vertex_resolution = new TLegend(0.7,0.7,0.9,0.9);
        c1_shower_vertex_resolution[j]->cd(iplot + nPlots*i+1);
        if (i == 0) {
          for (int p = 0; p < 3; p++) {
            h_shower_vertex_resolution_x[p][j][iplot]->Draw(draw_options[p]);
            if (p == 0) {
              h_shower_vertex_resolution_x[p][j][iplot]->SetStats(0);
              h_shower_vertex_resolution_x[p][j][iplot]->SetMaximum(1.05*maxVal_shower_vertex_resolution_x[j][iplot]);
            }
            h_shower_vertex_resolution_x[p][j][iplot]->SetLineColor(color_options[p]);
            lg_shower_vertex_resolution->AddEntry(h_shower_vertex_resolution_x[p][j][iplot], TString::Format("plane %d", p), "l");
          }
        }
        else if (i == 1) {
          for (int p = 0; p < 3; p++) {
            h_shower_vertex_resolution_y[p][j][iplot]->Draw(draw_options[p]);
            if (p == 0) {
              h_shower_vertex_resolution_y[p][j][iplot]->SetStats(0);
              h_shower_vertex_resolution_y[p][j][iplot]->SetMaximum(1.05*maxVal_shower_vertex_resolution_y[j][iplot]);
            }
            h_shower_vertex_resolution_y[p][j][iplot]->SetLineColor(color_options[p]);
            lg_shower_vertex_resolution->AddEntry(h_shower_vertex_resolution_y[p][j][iplot], TString::Format("plane %d", p), "l");
          }
        }
        else {
          for (int p = 0; p < 3; p++) {
            h_shower_vertex_resolution_z[p][j][iplot]->Draw(draw_options[p]);
            if (p == 0) {
              h_shower_vertex_resolution_z[p][j][iplot]->SetStats(0);
              h_shower_vertex_resolution_z[p][j][iplot]->SetMaximum(1.05*maxVal_shower_vertex_resolution_z[j][iplot]);
            }
            h_shower_vertex_resolution_z[p][j][iplot]->SetLineColor(color_options[p]);
            lg_shower_vertex_resolution->AddEntry(h_shower_vertex_resolution_z[p][j][iplot], TString::Format("plane %d", p), "l");
          }
        }
        lg_shower_vertex_resolution->Draw();
      }
    }
    c1_shower_vertex_resolution[j]->Write(TString::Format("h_shower_vertex_resolution_%s", shower_particle[j]));
    c1_shower_vertex_resolution[j]->SaveAs(TString::Format("h_shower_vertex_resolution_%s.pdf", shower_particle[j]));
  }

  // --------------------------------------------
  TCanvas *c1_shower_direction_resolution[6];
  for (int j=0; j<6; ++j) {
    c1_shower_direction_resolution[j] = new TCanvas(TString::Format("c1_shower_direction_resolution_%s",shower_particle[j]), TString::Format("c1_shower_direction_resolution_%s",shower_particle[j]), 1000, 800);
    c1_shower_direction_resolution[j]->Divide(2,2);
    for (int iplot = 0; iplot < nPlots; ++iplot) {
      auto lg_shower_direction_resolution = new TLegend(0.1,0.7,0.3,0.9);
      c1_shower_direction_resolution[j]->cd(iplot+1);
      for (int p = 0; p < 3; p++) {
        h_shower_direction_resolution[p][j][iplot]->Draw(draw_options[p]);
        if (p==0) {
          h_shower_direction_resolution[p][j][iplot]->SetStats(0);
          h_shower_direction_resolution[p][j][iplot]->SetMaximum(1.05*maxVal_shower_direction_resolution[j][iplot]);
        } 
        h_shower_direction_resolution[p][j][iplot]->SetLineColor(color_options[p]);
        lg_shower_direction_resolution->AddEntry(h_shower_direction_resolution[p][j][iplot], TString::Format("plane %d", p), "l");
      }
      lg_shower_direction_resolution->Draw();
    }
    c1_shower_direction_resolution[j]->Write(TString::Format("h_shower_direction_resolution_%s", shower_particle[j]));
    c1_shower_direction_resolution[j]->SaveAs(TString::Format("h_shower_direction_resolution_%s.pdf", shower_particle[j]));
  }

/*
  if(TEfficiency::CheckConsistency(*h_Ee_num,*h_Ee_den)){
    h_Eff_Ee = new TEfficiency("h_Eff_Ee", *h_Ee_num,*h_Ee_den);
    h_Eff_Ee->Write("h_Eff_Ee");
    h_Eff_Ee->SaveAs("h_Eff_Ee.pdf");
  }
*/
  // ---------------------
  TCanvas *c1_grEff_Ee[3];
  //TCanvas *c1_grEff_Ee = new TCanvas("c1_grEff_Ee","c1_grEff_Ee", 800,600);
  TGraphAsymmErrors *grEff_Ee[3][4];
  for (int p=0; p< 3; p++) {
    c1_grEff_Ee[p] = new TCanvas(TString::Format("c1_grEff_Ee_plane%d",p),"c1_grEff_Ee", 800,600);
    auto lg_grEff_Ee = new TLegend(0.7,0.2,0.9,0.4);
    for (int iplot = 0; iplot < nPlots; iplot++) {
      grEff_Ee[p][iplot] = new TGraphAsymmErrors(h_Ee_num[p][iplot], h_Ee_den);
      grEff_Ee[p][iplot]->Draw(eff_draw_options[iplot]);
      if (iplot == 0) {
        grEff_Ee[p][iplot]->SetTitle(TString::Format("Efficiency (plane %d)",p));
        grEff_Ee[p][iplot]->GetXaxis()->SetTitle("E_{e} [MeV]");
        grEff_Ee[p][iplot]->GetYaxis()->SetTitle("#epsilon");
        grEff_Ee[p][iplot]->GetYaxis()->SetRangeUser(0, 1.05);
      }
      grEff_Ee[p][iplot]->SetLineColor(eff_color_options[iplot]);
      lg_grEff_Ee->AddEntry(grEff_Ee[p][iplot], cut_options[iplot], "l");
    }
    lg_grEff_Ee->Draw();
    c1_grEff_Ee[p]->Write(TString::Format("grEff_Ee_plane%d",p));
    c1_grEff_Ee[p]->SaveAs(TString::Format("grEff_Ee_plane%d.pdf",p));
  }
  
  // ---------------------------
  std::vector<double> recon_energy;
  std::vector<double> recon_energy_error;
  std::vector<double> true_energy;
  std::vector<double> true_energy_error;
  recon_energy.clear();
  recon_energy_error.clear();
  true_energy.clear();
  true_energy_error.clear();
  TCanvas *c1_h_energy_shower_primary_e = new TCanvas("c1_h_energy_shower_primary_e","c1_h_energy_shower_primary_e", 1200,900);
  c1_h_energy_shower_primary_e->Divide(4,3);
  for (int n = 0; n < npoints_energy; n++) {
    c1_h_energy_shower_primary_e->cd(n+1);
    h_energy_shower_primary_e[n]->Draw();
    gStyle->SetStatY(0.9);// Set y-position (fraction of pad size)
    gStyle->SetStatX(0.85);// Set x-position (fraction of pad size)
    gStyle->SetStatW(0.25);// Set width of stat-box (fraction of pad size)
    gStyle->SetStatH(0.2);// Set height of stat-box (fraction of pad size)
    if (n > 1 && n < 9) {
      double left_shift = 0.2*n - 0.8;
      double right_shift = 0.5-0.1*n;
      auto myfunc = new TF1("myfunc","gaus",(n-1-left_shift)*step_energy,(n+right_shift)*step_energy);
      double par[3];
      h_energy_shower_primary_e[n]->Fit(myfunc, "", "", (n-1-left_shift)*step_energy,(n+right_shift)*step_energy);
      myfunc->GetParameters(&par[0]);
      cout << par[0] << " +/-  " << myfunc->GetParError(0) << endl;
      cout << par[1] << " +/-  " << myfunc->GetParError(1) << endl;
      cout << par[2] << " +/-  " << myfunc->GetParError(2) << endl;
      recon_energy.push_back(par[1]);
      recon_energy_error.push_back(myfunc->GetParError(1));
      true_energy.push_back((n+0.5)*step_energy);
      true_energy_error.push_back(0.5*step_energy);
    }
  }
  c1_h_energy_shower_primary_e->Write("h_energy_shower_primary_e");
  c1_h_energy_shower_primary_e->SaveAs("h_energy_shower_primary_e.pdf");

  // --------------------------------------------
  std::vector<std::vector<double>> energy_recon;
  std::vector<std::vector<double>> energy_recon_error;
  std::vector<std::vector<double>> energy_true;
  std::vector<std::vector<double>> energy_true_error;
  energy_recon.clear();
  energy_recon_error.clear();
  energy_true.clear();
  energy_true_error.clear();
  TCanvas *c1_h_primary_e_shower_recoEnergy[3];
  for (int p = 0; p < 3; p++) {
    std::vector<double> energy_recon_per_plane;
    std::vector<double> energy_recon_per_plane_error;
    std::vector<double> energy_true_per_plane;
    std::vector<double> energy_true_per_plane_error;
    energy_recon_per_plane.clear();
    energy_recon_per_plane_error.clear();
    energy_true_per_plane.clear();
    energy_true_per_plane_error.clear();

    c1_h_primary_e_shower_recoEnergy[p] = new TCanvas(TString::Format("c1_h_primary_e_shower_recoEnergy_plane%d",p),"c1_h_primary_e_shower_recoEnergy", 1200,900);
    c1_h_primary_e_shower_recoEnergy[p]->Divide(4,3);
    for (int n = 0; n < npoints_energy; n++) {
      c1_h_primary_e_shower_recoEnergy[p]->cd(n+1);
      h_primary_e_shower_recoEnergy[p][n]->Draw();
      if (n > 1 & n < 9) {
        double left_shift;
        double right_shift;
        if (p == 0) {
          left_shift = 0.4*n - 0.55;
          right_shift = 0.5-0.25*n;
        }
        else if (p == 1) {
          left_shift = 0.5*n - 0.8;
          right_shift = 0.5-0.25*n;
        }
        else {
          left_shift = 0.2*n - 0.8;
          right_shift = 0.5-0.1*n;
        }
        auto myfunc = new TF1("myfunc","gaus",(n-1-left_shift)*step_energy,(n+right_shift)*step_energy);
        double par[3];
        h_primary_e_shower_recoEnergy[p][n]->Fit(myfunc, "", "", (n-1-left_shift)*step_energy,(n+right_shift)*step_energy);
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
    c1_h_primary_e_shower_recoEnergy[p]->Write(TString::Format("h_primary_e_shower_recoEnergy_plane%d",p));
    c1_h_primary_e_shower_recoEnergy[p]->SaveAs(TString::Format("h_primary_e_shower_recoEnergy_plane%d.pdf",p));
  }

  // ---------------------------------------
  TCanvas *c1_grErr_shower_energy_correction[3];
  TGraphErrors *grErr_shower_energy_corr[3];
  for (int p = 0; p < 3; p++) {
    c1_grErr_shower_energy_correction[p] = new TCanvas(TString::Format("c1_grErr_shower_energy_correction_plane%d",p), "c1_grErr_shower_energy_correction", 800,600);
    TPaveText *pt_shower_energy_corr = new TPaveText(400.0,800.0,850.0,1000.0);
    grErr_shower_energy_corr[p] = new TGraphErrors(energy_recon[p].size(), &energy_true[p][0],&energy_recon[p][0], &energy_true_error[p][0], &energy_recon_error[p][0]);
    auto f1 = new TF1("f1", "[0]+[1]*x", 400.0, 1800.0);
    f1->SetParameters(0.0,1.0);
    double par_corr[2];
    grErr_shower_energy_corr[p]->SetTitle(TString::Format("Shower energy correction (plane %d)", p));
    grErr_shower_energy_corr[p]->GetXaxis()->SetTitle("True electron energy [MeV]");
    grErr_shower_energy_corr[p]->GetYaxis()->SetTitle("Recon shower energy [MeV]");
    grErr_shower_energy_corr[p]->Draw("AP");
    grErr_shower_energy_corr[p]->Fit(f1,"","",400.0,1800.0);
    f1->GetParameters(&par_corr[0]);
    cout << par_corr[0] << endl;
    cout << par_corr[1] << endl;
    pt_shower_energy_corr->AddText(TString::Format("Slope: %f #pm %f",par_corr[1], f1->GetParError(1)));
    pt_shower_energy_corr->AddText(TString::Format("Intercept: %f #pm %f",par_corr[0], f1->GetParError(0)));
    pt_shower_energy_corr->Draw();
    c1_grErr_shower_energy_correction[p]->Write(TString::Format("grErr_shower_energy_corr_plane%d",p));
    c1_grErr_shower_energy_correction[p]->SaveAs(TString::Format("grErr_shower_energy_corr_plane%d.pdf",p));
  }


  // --------------------
  int npoints_select = recon_energy.size();
  //TGraph *gr_shower_energy_correction = new TGraph(npoints_select, &true_energy[0], &recon_energy[0]);
  TGraphErrors *gr_shower_energy_correction = new TGraphErrors(npoints_select, &true_energy[0], &recon_energy[0], &true_energy_error[0], &recon_energy_error[0]);
  auto f1 = new TF1("f1", "[0]+[1]*x", 400.0, 1800.0);
  f1->SetParameters(0.0,1.0);
  double par_corr[2];

  TCanvas *c1 = new TCanvas("c1","c1", 800,600);
  TPaveText *pt_energy_corr = new TPaveText(400.0,1000.0,850.0,1200.0);
  gr_shower_energy_correction->GetXaxis()->SetTitle("True electron energy [MeV]");
  gr_shower_energy_correction->GetYaxis()->SetTitle("Recon shower energy [MeV]");
  gr_shower_energy_correction->SetTitle("Shower energy correction");
  gr_shower_energy_correction->SetMarkerStyle(21);
  gr_shower_energy_correction->SetMarkerColor(4);
  gr_shower_energy_correction->SetMarkerSize(0.5);
  gr_shower_energy_correction->Draw("AP");
  gr_shower_energy_correction->Fit(f1,"","",400.0,1800.0);
  f1->GetParameters(&par_corr[0]);
  cout << par_corr[0] << endl;
  cout << par_corr[1] << endl;
  pt_energy_corr->AddText(TString::Format("Slope: %f #pm %f ",par_corr[1], f1->GetParError(1)));
  pt_energy_corr->AddText(TString::Format("Intercept: %f #pm %f",par_corr[0], f1->GetParError(0)));
  pt_energy_corr->Draw();

  c1->SaveAs("gr_shower_energy_correction.pdf");


  output.Write();
  output.Close();

  myfile.close();
}

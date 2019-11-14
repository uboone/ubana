
using namespace std;

void plot_template() {
  
  // input
  TFile *file_a = TFile::Open("output_ShowerTemplateMaker.root");

  // plot options
  int nSelectParticles = 4;
  int nSelectParticlesAll = 5;
  int nSelectEnergies = 11;
  int nSelectParticleProfiles = 7;

  const char * shower_particle[5] = {"electron", "photon", "proton", "pion", "others"}; // 0: primary electron; 1: photon; 2: proton; 3: pion; 4: others
  const char * energy_option[11] = {"0-10GeV", "0-100MeV", "100-200MeV", "200-300MeV","300-400MeV", "400-500MeV","500-600MeV", "600-700MeV","700-800MeV", "800-900MeV","900-1000MeV"};
  const char * profile_option[7] = {"long", "trans", "trans_1", "trans_2", "trans_3", "trans_4","trans_5"};

  const int color_options[4] = {1, 2, 3, 4}; // red, green, blue
  const int style_options[4] = {1, 2, 3, 4}; // red, green, blue
  const char * draw_options[4] = {"", "SAME", "SAME", "SAME"};


  // profile histogram
  TProfile *showerProfile[5][11][7]; // [particle][energy][profile]
  for (int i = 0; i < nSelectParticlesAll; i++) {
    for (int j = 0; j < nSelectEnergies; j++) {
      for (int k = 0; k < nSelectParticleProfiles; k++) {
        showerProfile[i][j][k] = (TProfile*)file_a->Get(TString::Format("showerProfile_%s_%s_%s",shower_particle[i], energy_option[j], profile_option[k]) );
      }
    }
  }

  
  // --------------------
  TCanvas *c1_profile[11];
  for (int j = 0; j < nSelectEnergies; j++) {
    c1_profile[j] = new TCanvas(TString::Format("c1_profile_%d",j),TString::Format("c1_profile_%d",j),1200,800);
    c1_profile[j]->Divide(3,2);
    for (int k = 0; k < nSelectParticleProfiles; k++) {
      if (k == 1) continue;
      if (k == 0)c1_profile[j]->cd(k+1);
      if (k > 1) c1_profile[j]->cd(k);
     
      double maxVal_charge = -999.0;
      for (int i = 0; i < nSelectParticles; i++) {
        if (showerProfile[i][j][k]->GetMaximum() > maxVal_charge) {
          maxVal_charge = showerProfile[i][j][k]->GetMaximum();
        }
      }

      auto lg_profile = new TLegend(0.7,0.7,0.9,0.9);
      for (int i = 0; i < nSelectParticles; i++) {
        showerProfile[i][j][k]->Draw(draw_options[i]);
        if (i == 0) showerProfile[i][j][k]->SetMaximum(1.05*maxVal_charge);
        showerProfile[i][j][k]->SetLineColor(color_options[i]);
        showerProfile[i][j][k]->SetMarkerColor(color_options[i]);
        //showerProfile[i][j][k]->SetMarkerStyle(style_options[i]);
        lg_profile->AddEntry(showerProfile[i][j][k], shower_particle[i]);
      }
      lg_profile->Draw();
    }
    c1_profile[j]->SaveAs(TString::Format("profile_%s.pdf",energy_option[j]));
  }
  
}

#define ana_cxx
#include "ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void template_ana::Loop()
{
  // output file
  TFile output("output_ShowerTemplateMaker.root","recreate");

  if (fChain == 0) return;
  
  // plot options
  int nSelectParticles = 4;
  int nSelectParticlesAll = 5;
  int nSelectEnergies = 11;
  int nSelectParticleProfiles = 7;

  const char * shower_particle[5] = {"electron", "photon", "proton", "pion", "others"}; // 0: primary electron; 1: photon; 2: proton; 3: pion; 4: others
  const char * energy_option[11] = {"0-10GeV", "0-100MeV", "100-200MeV", "200-300MeV","300-400MeV", "400-500MeV","500-600MeV", "600-700MeV","700-800MeV", "800-900MeV","900-1000MeV"};
  const char * profile_option[7] = {"long", "trans", "trans_1", "trans_2", "trans_3", "trans_4","trans_5"};

  // cut options
  double cutShowerPurity = 0.7;
  double cutShowerCompleteness = cutShowerPurity;
 
  // dimention options
  int nBins_energy = 100;
  double lowBin_energy = 0.0;
  double highBin_energy = 10000.0; // [MeV]

  int nBins_longDist = 100;
  double lowBin_longDist = 0.0;
  double highBin_longDist = 25.0; // [t]; radiation length 1 t = 14 cm

  int nBins_tranDist = 32; // todo: rethink about the bin width
  double lowBin_tranDist = -8.0;
  double highBin_tranDist = 8.0;

  int nBins_charge = 100;
  double lowBin_charge = 0.0;
  double highBin_charge = 150000.0; // [ADC]; todo: switch to # electron

  // profile histogram
  TProfile *showerProfile[5][11][7]; // [particle][energy][profile]
  for (int i = 0; i < nSelectParticlesAll; i++) {
    for (int j = 0; j < nSelectEnergies; j++) {
      for (int k = 0; k < nSelectParticleProfiles; k++) {
        if (k == 0) {
          showerProfile[i][j][k] = new TProfile(TString::Format("showerProfile_%s_%s_%s",shower_particle[i], energy_option[j], profile_option[k]),TString::Format("showerProfile_%s_%s_%s;t(X0=14cm);Q",shower_particle[i], energy_option[j], profile_option[k]), nBins_longDist, lowBin_longDist, highBin_longDist);
        }
        else {
          showerProfile[i][j][k] = new TProfile(TString::Format("showerProfile_%s_%s_%s",shower_particle[i], energy_option[j], profile_option[k]),TString::Format("showerProfile_%s_%s_%s;dist(1cm);Q",shower_particle[i], energy_option[j], profile_option[k]), nBins_tranDist, lowBin_tranDist, highBin_tranDist);
        }
      }
    }
  }

  // loop over entries
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries: " << nentries << endl;

  int nShowerReco = 0;
  int nShowerHitInfo = 0;
  int nShowerCuts = 0;
  int nShowerNoCuts = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  
    nShowerReco += showerEnergyReco->size();
    nShowerHitInfo += showerHitLongDist->size();
    
    // only for checking output purpose
    if (jentry < 5) {
      cout << "showerEnergyReco->size(): " << showerEnergyReco->size() << endl;
      cout << "showerHitTranDist->size(): " << showerHitTranDist->size() << endl;
    }

    // showers per event 
    int n_temp_shower = 0;
    for (int ishower = 0; ishower < showerEnergyReco->size(); ishower++) {
      nShowerNoCuts += 1;
      
      // purity and completeness cut: calculated only based on MC; in fact, this is not accurate when dealing when overlay sample due to the case that a shower may have hits from data.
      // we also require the total hit charge on plane 2 is bigger than 0, which is equivalent to "E_reco_plane2 > 0". slope: 0.813455; intercept: -17.490544; 21.501551 = (0+17.490544)/0.813455
      if ((*sh_purity)[ishower] > cutShowerPurity && (*sh_completeness)[ishower] > cutShowerCompleteness && (*showerEnergyReco)[ishower][2] > 21.51) {
        nShowerCuts += 1;
        n_temp_shower += 1;

        // for each effective shower, we save the 1D histogram
        auto *ltemp = new TH1F("ltemp","ltemp",nBins_longDist, lowBin_longDist, highBin_longDist);
        auto *ttemp = new TH1F("ttemp","ttemp", nBins_tranDist, lowBin_tranDist, highBin_tranDist);
        auto *ttemp_1 = new TH1F("ttemp_1","ttemp_1", nBins_tranDist, lowBin_tranDist, highBin_tranDist);
        auto *ttemp_2 = new TH1F("ttemp_2","ttemp_2", nBins_tranDist, lowBin_tranDist, highBin_tranDist);
        auto *ttemp_3 = new TH1F("ttemp_3","ttemp_3", nBins_tranDist, lowBin_tranDist, highBin_tranDist);
        auto *ttemp_4 = new TH1F("ttemp_4","ttemp_4", nBins_tranDist, lowBin_tranDist, highBin_tranDist);
        auto *ttemp_5 = new TH1F("ttemp_5","ttemp_5", nBins_tranDist, lowBin_tranDist, highBin_tranDist);
          
        // shower hits
        for (int ihit = 0; ihit < (*showerHitLongDist)[n_temp_shower-1].size(); ihit++) {
          ltemp->Fill((*showerHitLongDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          ttemp->Fill( (*showerHitTranDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          if ( (*showerHitLongDist)[n_temp_shower-1][ihit] < 1 ) {
            ttemp_1->Fill( (*showerHitTranDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          }
          else if ( (*showerHitLongDist)[n_temp_shower-1][ihit] < 2 ) {
            ttemp_2->Fill( (*showerHitTranDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          }
          else if ( (*showerHitLongDist)[n_temp_shower-1][ihit] < 3 ) {
            ttemp_3->Fill( (*showerHitTranDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          }
          else if ( (*showerHitLongDist)[n_temp_shower-1][ihit] < 4 ) {
            ttemp_4->Fill( (*showerHitTranDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          }
          else if ( (*showerHitLongDist)[n_temp_shower-1][ihit] < 5 ) {
            ttemp_5->Fill( (*showerHitTranDist)[n_temp_shower-1][ihit], (*showerHitChargeQ)[n_temp_shower-1][ihit] );
          }
        } // end of for ihit

        // find the particle index
        int particle_index = -99;
        if ( abs((*showerParticlePDG)[ishower]) == 11 && (*sh_hasPrimary_e)[ishower] == 1) {
            particle_index = 0;  
        }
        else if ( abs((*showerParticlePDG)[ishower]) == 22 ) {
          particle_index = 1;
        }
        else if ( abs((*showerParticlePDG)[ishower]) == 2212 ) {
          particle_index = 2;
        }
        else if ( abs((*showerParticlePDG)[ishower]) == 211 ) {
          particle_index = 3;
        }
        else {
          particle_index = 4;
        }

        // find the energy index
        int energy_index = -99;
        for (int ienergy = 0; ienergy < nSelectEnergies-1; ienergy++) {
          if ( ( (*showerEnergyReco)[ishower][2] > ((highBin_energy-lowBin_energy) / nBins_energy) * ienergy + lowBin_energy ) && ( (*showerEnergyReco)[ishower][2] <= ((highBin_energy-lowBin_energy) / nBins_energy) * (ienergy+1) + lowBin_energy ) ) {
            energy_index = ienergy + 1;
          } 
        }
        if (particle_index == -99) {
        cout << "energy_index: " << energy_index << endl;
        cout << "particle_index: " << particle_index << endl;
          cout << "(*showerParticlePDG)[ishower]: " << (*showerParticlePDG)[ishower] << endl;
        }
        
        // fill long profile
        for (int ilong = 0; ilong < nBins_longDist; ilong++) {
          showerProfile[particle_index][0][0]->Fill( ltemp->GetBinCenter(ilong+1), ltemp->GetBinContent(ilong+1) );
          if ( energy_index > 0 && (energy_index < nSelectEnergies) ) {
            showerProfile[particle_index][energy_index][0]->Fill( ltemp->GetBinCenter(ilong+1), ltemp->GetBinContent(ilong+1) );
          }
        }
        // fill trans profile
        for (int itrans = 0; itrans < nBins_tranDist; itrans++) {
          showerProfile[particle_index][0][1]->Fill( ttemp->GetBinCenter(itrans+1), ttemp->GetBinContent(itrans+1) );
          showerProfile[particle_index][0][2]->Fill( ttemp_1->GetBinCenter(itrans+1), ttemp_1->GetBinContent(itrans+1) );
          showerProfile[particle_index][0][3]->Fill( ttemp_2->GetBinCenter(itrans+1), ttemp_2->GetBinContent(itrans+1) );
          showerProfile[particle_index][0][4]->Fill( ttemp_3->GetBinCenter(itrans+1), ttemp_3->GetBinContent(itrans+1) );
          showerProfile[particle_index][0][5]->Fill( ttemp_4->GetBinCenter(itrans+1), ttemp_4->GetBinContent(itrans+1) );
          showerProfile[particle_index][0][6]->Fill( ttemp_5->GetBinCenter(itrans+1), ttemp_5->GetBinContent(itrans+1) );
          
          if ( energy_index > 0 && (energy_index < nSelectEnergies) ) {
            showerProfile[particle_index][energy_index][1]->Fill( ttemp->GetBinCenter(itrans+1), ttemp->GetBinContent(itrans+1) );
            showerProfile[particle_index][energy_index][2]->Fill( ttemp_1->GetBinCenter(itrans+1), ttemp_1->GetBinContent(itrans+1) );
            showerProfile[particle_index][energy_index][3]->Fill( ttemp_2->GetBinCenter(itrans+1), ttemp_2->GetBinContent(itrans+1) );
            showerProfile[particle_index][energy_index][4]->Fill( ttemp_3->GetBinCenter(itrans+1), ttemp_3->GetBinContent(itrans+1) );
            showerProfile[particle_index][energy_index][5]->Fill( ttemp_4->GetBinCenter(itrans+1), ttemp_4->GetBinContent(itrans+1) );
            showerProfile[particle_index][energy_index][6]->Fill( ttemp_5->GetBinCenter(itrans+1), ttemp_5->GetBinContent(itrans+1) );
          }
        }

        delete ltemp;
        delete ttemp;
        delete ttemp_1;
        delete ttemp_2;
        delete ttemp_3;
        delete ttemp_4;
        delete ttemp_5;

      } // end of if p > cutShowerPurity && c > cutShowerCompleteness && E_reco_plane2> 0

    } // end of loop ishower


  } // end of loop for jentry

  cout << "nShowerReco: " << nShowerReco << endl;
  cout << "nShowerHitInfo: " << nShowerHitInfo << endl;
  cout << "nShowerNoCuts: " << nShowerNoCuts << endl;
  cout << "nShowerCuts: " << nShowerCuts << endl;

  // ---------------------------------
  TCanvas *c1_profile[11][7];
  for (int j = 0; j < nSelectEnergies; j++) {
    for (int k = 0; k < nSelectParticleProfiles; k++) {
      c1_profile[j][k] = new TCanvas(TString::Format("c1_profile_%s_%s",energy_option[j],profile_option[k]), TString::Format("c1_profile_%s_%s",energy_option[j],profile_option[k]), 1000, 800);
      c1_profile[j][k]->Divide(2,2);
      for (int i = 0; i < nSelectParticles; i++) {
        c1_profile[j][k]->cd(i+1);
        showerProfile[i][j][k]->Draw();
      }
      c1_profile[j][k]->Write(TString::Format("profile_%s_%s",energy_option[j], profile_option[k]));
      c1_profile[j][k]->SaveAs(TString::Format("profile_%s_%s.pdf",energy_option[j], profile_option[k]));
    }
  }

  // output
  output.Write();
  output.Close();
}

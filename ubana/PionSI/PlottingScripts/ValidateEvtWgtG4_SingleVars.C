// This script is intended to validate that the Geant4WeightCalc eventweight calculator is working correctly
// It is intended to be used with unisims so we can look at them individually and check that they are right
// This script will not work with multisims (write another script for that!)

void ValidateEvtWgtG4_SingleVars(std::string infilename, std::string outfilename){

   gInterpreter->GenerateDictionary("vector<map<string,double> >", "vector");
   gSystem->Load("AutoDict_vector_map_string_double____cxx.so");

   // Get trees from file
   TFile *fin = new TFile(infilename.c_str(),"open");
   TTree *tpart = (TTree*)fin->Get("eventweightreint/ByParticleValidTree_211");
   TTree *tmctr = (TTree*)fin->Get("eventweightreint/ByMCTruthValidTree_211");

   // Variables to read from trees
   int p_event_num;
   tpart->SetBranchAddress("event_num",&p_event_num);
   int p_run_num;
   tpart->SetBranchAddress("run_num",&p_run_num);
   int p_subrun_num;
   tpart->SetBranchAddress("subrun_num",&p_subrun_num);
   std::vector<std::map<std::string, double>> *UniverseVals=nullptr;
   tpart->SetBranchAddress("UniverseVals", &UniverseVals);
   int p_pdg_to_reweight;
   tpart->SetBranchAddress("pdg_to_reweight",&p_pdg_to_reweight);
   std::vector<double> *p_inel_weight=nullptr;
   tpart->SetBranchAddress("inelastic_weight",&p_inel_weight);
   std::vector<double> *p_elastic_weight=nullptr;
   tpart->SetBranchAddress("elastic_weight",&p_elastic_weight);
   double p_track_length;
   tpart->SetBranchAddress("track_length",&p_track_length);
   int p_PDG;
   tpart->SetBranchAddress("PDG",&p_PDG);
   std::string *p_final_proc=nullptr;
   tpart->SetBranchAddress("final_proc",&p_final_proc);
   double p_init_momentum;
   tpart->SetBranchAddress("init_momentum",&p_init_momentum);
   double p_final_momentum;
   tpart->SetBranchAddress("final_momentum",&p_final_momentum);
   std::vector<double> *p_energies_inel=nullptr;
   tpart->SetBranchAddress("energies_inel",&p_energies_inel);
   std::vector<int> *p_sliceInts_inel=nullptr;
   tpart->SetBranchAddress("sliceInts_inel",&p_sliceInts_inel);
   std::vector<double> *p_energies_el=nullptr;
   tpart->SetBranchAddress("energies_el",&p_energies_el);
   std::vector<int> *p_sliceInts_el=nullptr;
   tpart->SetBranchAddress("sliceInts_el",&p_sliceInts_el);
   int p_nElasticScatters;
   tpart->SetBranchAddress("nElasticScatters",&p_nElasticScatters);

   int e_event_num;
   tmctr->SetBranchAddress("event_num",&e_event_num);
   int e_run_num;
   tmctr->SetBranchAddress("run_num",&e_run_num);
   int e_subrun_num;
   tmctr->SetBranchAddress("subrun_num",&e_subrun_num);
   int e_pdg_to_reweight;
   tmctr->SetBranchAddress("pdg_to_reweight",&e_pdg_to_reweight);
   std::vector<std::vector<double>> *e_inel_weight=nullptr;
   tmctr->SetBranchAddress("inelastic_weight",&e_inel_weight);
   std::vector<std::vector<double>> *e_elastic_weight=nullptr;
   tmctr->SetBranchAddress("elastic_weight",&e_elastic_weight);

   // UniverseVals is filled once and should be identical for all entries in both trees -- get it from entry 0
   tpart->GetEntry(0);

   // This script expects that each parameter will have a +/- 1 sigma variation stored in the tree. Work out the indices for those variations by looking through UniverseVals
   std::map<std::string,std::pair<int,double>> VarsDown;
   std::map<std::string,std::pair<int,double>> VarsUp;

   for (int i_univ=0; i_univ<UniverseVals->size(); i_univ++){
      int n_var = 0;
      double val = -999;
      std::string par;
      std::map<std::string, double>:: iterator it_par;
      for (it_par = UniverseVals->at(i_univ).begin(); it_par != UniverseVals->at(i_univ).end(); it_par++){
         // std::cout << it_par->first << ": " << it_par->second << std::endl;
         if (it_par->second != 1){
            n_var++;
            val = it_par->second;
            par = it_par->first;
         }
      } // it_par
      // Double-check that this parameter set only had a single parameter varied. If so, save this to make plots
      if (n_var != 1){
         std::cout << "Error: Parameter set for Universe " << i_univ << " has " << n_var << " varied parameters. Not plotting this universe." << std::endl;
      }
      else{
         std::pair<int,double> tmp = std::make_pair(i_univ,val);
         if (val<1){
            VarsDown[par] = tmp;
         }
         else{
            VarsUp[par] = tmp;
         }
      }
   } // i_univ

   // Check parameters we're plotting
   std::cout << "Plotting -1sigma variations for " << VarsDown.size() << " parameters: " << std::endl;
   std::map<std::string,std::pair<int,double>>:: iterator p;
   for (p = VarsDown.begin(); p != VarsDown.end(); p++)
      std::cout << "    " << p->first << ": " << p->second.second << " (index " << p->second.first << ")" << ")\n";
   std::cout << "Plotting +1sigma variations for " << VarsUp.size() << " parameters: " << std::endl;
   for (p = VarsUp.begin(); p != VarsUp.end(); p++)
      std::cout << "    " << p->first << ": " << p->second.second << " (index " << p->second.first << ")" << ")\n";

   if (VarsDown.size() != VarsUp.size()) std::cout << "WARNING: VarsDown.size() = " << VarsDown.size() << ", VarsUp.size() = " << VarsUp.size() << ". This script uses VarsDown to determine which parameters to plot -- that may cause you problems!" << std::endl;

   // Get number of parameter variations to plot and instantiate histograms
   // There are always 2 histograms: one incident and one interacting (for calculating the effective cross section)
   int n_vars = VarsDown.size();
   TH1D *h_part_nom_inel[2];
   TH1D *h_part_nom_el[2];
   TH1D *h_part_up[2][n_vars];
   TH1D *h_part_down[2][n_vars];

   h_part_nom_inel[0] = new TH1D("h_part_inc_nom",";Momentum (MeV/c);Inelastic cross section (mb)",20,0,2000);
   h_part_nom_inel[1] = new TH1D("h_part_int_nom",";Momentum (MeV/c);Inelastic cross section (mb)",20,0,2000);
   h_part_nom_el[0] = new TH1D("h_part_inc_nom",";Momentum (MeV/c);Elastic cross section (mb)",20,0,2000);
   h_part_nom_el[1] = new TH1D("h_part_int_nom",";Momentum (MeV/c);Elastic cross section (mb)",20,0,2000);

   int i_var = 0;
   std::map<int, std::string> i_var_to_str;
   for (p = VarsDown.begin(); p != VarsDown.end(); p++){
      if (p->first.find(std::string("Elast")) != std::string::npos){
         h_part_up[0][i_var] = new TH1D(Form("h_part_inc_up_%i",i_var),Form("%s;Momentum (MeV/c);Elastic cross section (mb)",p->first.c_str()),20,0,2000);
         h_part_down[0][i_var] = new TH1D(Form("h_part_inc_down_%i",i_var),Form("%s;Momentum (MeV/c);Elastic cross section (mb)",p->first.c_str()),20,0,2000);
         h_part_up[1][i_var] = new TH1D(Form("h_partint__up_%i",i_var),Form("%s;Momentum (MeV/c);Elastic cross section (mb)",p->first.c_str()),20,0,2000);
         h_part_down[1][i_var] = new TH1D(Form("h_part_int_down_%i",i_var),Form("%s;Momentum (MeV/c);Elastic cross section (mb)",p->first.c_str()),20,0,2000);
      }
      else{
         h_part_up[0][i_var] = new TH1D(Form("h_part_inc_up_%i",i_var),Form("%s;Momentum (MeV/c);Inelastic cross section (mb)",p->first.c_str()),20,0,2000);
         h_part_down[0][i_var] = new TH1D(Form("h_part_inc_down_%i",i_var),Form("%s;Momentum (MeV/c);Inelastic cross section (mb)",p->first.c_str()),20,0,2000);
         h_part_up[1][i_var] = new TH1D(Form("h_partint__up_%i",i_var),Form("%s;Momentum (MeV/c);Inelastic cross section (mb)",p->first.c_str()),20,0,2000);
         h_part_down[1][i_var] = new TH1D(Form("h_part_int_down_%i",i_var),Form("%s;Momentum (MeV/c);Inelastic cross section (mb)",p->first.c_str()),20,0,2000);
      }


      i_var_to_str[i_var] = p->first;
      i_var++;
   }


   TH1D *h_evt_wgt[2][UniverseVals->size()];
   TH1D *h_evt_bypart_wgt[2][UniverseVals->size()];

   for (int i_univ=0; i_univ<UniverseVals->size(); i_univ++){
      double val = -999;
      std::string par;
      std::map<std::string, double>:: iterator it_par;
      for (it_par = UniverseVals->at(i_univ).begin(); it_par != UniverseVals->at(i_univ).end(); it_par++){
         if (it_par->second != 1){
            val = it_par->second;
            par = it_par->first;
         }
      } // it_par

      h_evt_wgt[0][i_univ] = new TH1D(Form("h_evt_wgt_inel_%i",i_univ),Form("InelWgt_%s=%f;Weight;No. Entries",par.c_str(), val),50,-2,3);
      h_evt_bypart_wgt[0][i_univ] = new TH1D(Form("h_evt_bypart_wgt_inel_%i",i_univ),Form("InelWgt_%s=%f;Weight;No. Entries",par.c_str(), val),50,-2,3);

      h_evt_wgt[1][i_univ] = new TH1D(Form("h_evt_wgt_elastic_%i",i_univ),Form("ElasticWgt_%s=%f;Weight;No. Entries",par.c_str(), val),50,-2,3);
      h_evt_bypart_wgt[1][i_univ] = new TH1D(Form("h_evt_bypart_wgt_elastic_%i",i_univ),Form("ElasticWgt_%s=%f;Weight;No. Entries",par.c_str(), val),50,-2,3);
   } // i_univ



   std::map<int,std::string>:: iterator p1;
   for (p1 = i_var_to_str.begin(); p1 != i_var_to_str.end(); p1++)
      std::cout << "    " << p1->first << ": " << p1->second << ")\n";

   // Ok, now we're ready to make some plots!

   // ------------------------------------------------------ //
   //              Validation by particle                    //
   // ------------------------------------------------------ //

   TCanvas *c1 = new TCanvas();
   int run=-999, subrun=-999, event=-999;
   int n_entries = tpart->GetEntries();
   // if (n_entries>50000) n_entries = 50000;
   for (int i_p=0; i_p<n_entries; i_p++){
      if (i_p%1000==0) std::cout << i_p << "/" << n_entries << std::endl;
      tpart->GetEntry(i_p);

      if (p_run_num!=run || p_subrun_num!=subrun || p_event_num!=event){
         run = p_run_num;
         subrun = p_subrun_num;
         event = p_event_num;
         // std::cout << "Run " << run << " Subrun " << subrun << " Event " << event << std::endl;
      }
      // std::cout << "   " << "pdg = " << p_PDG << ", weights = ";
      // for (int i_univ=0; i_univ<p_inel_weight->size(); i_univ++){
      //    std::cout << p_inel_weight->at(i_univ) << "/" << p_elastic_weight->at(i_univ) << " ";
      // }
      // std::cout << std::endl;




      // Check if this is the particle type we're trying to reweight -- if not, make sure it has a weight of 1, but then don't include it in the validation plots
      // std::cout << "Particle type " << p_PDG << ", when should be reweighting particles of type " << p_pdg_to_reweight << std::endl;
      if (p_PDG != p_pdg_to_reweight){
         for (int i_wgt=0; i_wgt<p_inel_weight->size();i_wgt++){
            // if (p_inel_weight->at(i_wgt)!=1) std::cout << "ERROR: see inelastic weight of " << p_inel_weight->at(i_wgt) << " for particle type " << p_PDG << ", when should be reweighting particles of type " << p_pdg_to_reweight << std::endl;
         }
         for (int i_wgt=0; i_wgt<p_elastic_weight->size();i_wgt++){
            // if (p_elastic_weight->at(i_wgt)!=1) std::cout << "ERROR: see elastic weight of " << p_elastic_weight->at(i_wgt) << " for particle type " << p_PDG << ", when should be reweighting particles of type " << p_pdg_to_reweight << std::endl;
         }
         continue;
      } // if p_PDG != p_pdg_to_reweight

      // Now fill histograms! Loop through p_energies
      // Note that p_energies_el and p_energies_inel are identical. Depending on whether we're looking at an inelastic or elastic parameter, we need to decide whether to look at sliceInts_el or sliceInts_inel
      // std::cout << p_energies_inel->size() << "   " << p_sliceInts->size() << std::endl;

      // Get particle mass so we can plot kinetic energy (not total energy as stored in the tree)
      double mass = -999;
      if (TMath::Abs(p_PDG) == 211) mass = 139.570;
      else if (p_PDG == 2212) mass = 938.272;
      else{
         std::cout << "Unknown mass for PDG " << p_PDG << ", not filling kinetic energy plots" << std::endl;
      }
      for (int in=0; in<p_energies_inel->size(); in++){
         // std::cout << "    " << p_sliceInts_inel->at(in) << "  " << p_sliceInts_el->at(in) << "  " << p_energies_inel->at(in) << std::endl;

         double p = TMath::Sqrt(p_energies_inel->at(in)*p_energies_inel->at(in)-mass*mass);

         h_part_nom_inel[0]->Fill(p);
         h_part_nom_el[0]->Fill(p);

         if (p_sliceInts_inel->at(in)==1){
            h_part_nom_inel[1]->Fill(p);
         }
         if (p_sliceInts_el->at(in)==1){
            h_part_nom_el[1]->Fill(p);
         }

         for (int i_par=0; i_par<i_var_to_str.size(); i_par++){
            std::string par = i_var_to_str[i_par];
            int i_univ = VarsDown[par].first;
            double wgt = p_inel_weight->at(i_univ);
            wgt *= p_elastic_weight->at(i_univ);
            h_part_down[0][i_par]->Fill(p,wgt);

            if (par.find(std::string("Elast")) != std::string::npos && p_sliceInts_el->at(in)==1){
               h_part_down[1][i_par]->Fill(p,wgt);
            }
            else if (par.find(std::string("Elast")) == std::string::npos && p_sliceInts_inel->at(in)==1){
               h_part_down[1][i_par]->Fill(p,wgt);
            }


            i_univ = VarsUp[par].first;
            wgt = p_inel_weight->at(i_univ);
            wgt *= p_elastic_weight->at(i_univ);
            h_part_up[0][i_par]->Fill(p,wgt);

            if (par.find(std::string("Elast")) != std::string::npos && p_sliceInts_el->at(in)==1){
               h_part_up[1][i_par]->Fill(p,wgt);
            }
            else if (par.find(std::string("Elast")) == std::string::npos && p_sliceInts_inel->at(in)==1){
               h_part_up[1][i_par]->Fill(p,wgt);
            }
         }
      }

   } // end loop over particles/entries in tpart (i_p)

   // Divide histograms to get cross sections
   h_part_nom_inel[1]->Sumw2();
   h_part_nom_el[1]->Sumw2();
   h_part_nom_inel[1]->Divide(h_part_nom_inel[0]);
   h_part_nom_el[1]->Divide(h_part_nom_el[0]);
   for (int i_par=0; i_par<i_var_to_str.size(); i_par++){
      h_part_down[1][i_par]->Sumw2();
      h_part_up[1][i_par]->Sumw2();
      h_part_down[1][i_par]->Divide(h_part_down[0][i_par]);
      h_part_up[1][i_par]->Divide(h_part_up[0][i_par]);
   }

   // std::cout << "h_part_nom_inel[0] integral = " << h_part_nom_inel[0]->Integral() << std::endl;
   // std::cout << "h_part_nom_inel[1] integral = " << h_part_nom_inel[1]->Integral() << std::endl;
   // std::cout << "h_part_nom_el[0] integral = " << h_part_nom_el[0]->Integral() << std::endl;
   // std::cout << "h_part_nom_el[1] integral = " << h_part_nom_el[1]->Integral() << std::endl;
   // for (int i_par=0; i_par<i_var_to_str.size(); i_par++){
   //    std::cout << i_var_to_str.at(i_par) << " h_part_up[0] integral = " << h_part_up[0][i_par]->Integral() << std::endl;
   //    std::cout << i_var_to_str.at(i_par) << " h_part_up[1] integral = " << h_part_up[1][i_par]->Integral() << std::endl;
   //    std::cout << i_var_to_str.at(i_par) << " h_part_down[0] integral = " << h_part_down[0][i_par]->Integral() << std::endl;
   //    std::cout << i_var_to_str.at(i_par) << " h_part_down[1] integral = " << h_part_down[1][i_par]->Integral() << std::endl;
   // }


   // Now make plots
   h_part_nom_inel[1]->SetLineColor(kBlack);
   h_part_nom_inel[1]->SetLineStyle(1);
   h_part_nom_inel[1]->SetLineWidth(2);
   h_part_nom_el[1]->SetLineColor(kBlack);
   h_part_nom_el[1]->SetLineStyle(1);
   h_part_nom_el[1]->SetLineWidth(2);
   TH1D *h_part_nom_inel_ratio = (TH1D*)h_part_nom_inel[1]->Clone("h_part_nom_inel_ratio");
   TH1D *h_part_nom_el_ratio = (TH1D*)h_part_nom_el[1]->Clone("h_part_nom_el_ratio");
   h_part_nom_inel_ratio->Sumw2();
   h_part_nom_el_ratio->Sumw2();
   h_part_nom_inel_ratio->Divide(h_part_nom_inel[1]);
   h_part_nom_el_ratio->Divide(h_part_nom_el[1]);
   for (int i_par=0; i_par<i_var_to_str.size(); i_par++){
      h_part_up[1][i_par]->SetLineColor(kRed);
      h_part_up[1][i_par]->SetLineStyle(2);
      h_part_up[1][i_par]->SetLineWidth(2);
      h_part_down[1][i_par]->SetLineColor(kBlue);
      h_part_down[1][i_par]->SetLineStyle(2);
      h_part_down[1][i_par]->SetLineWidth(2);

      h_part_up[1][i_par]->Draw("hist");
      if(i_var_to_str.at(i_par).find(std::string("Elast")) != std::string::npos){
         h_part_nom_el[1]->Draw("hist e1 same");
      }
      else{
         h_part_nom_inel[1]->Draw("hist e1 same");
      }
      h_part_down[1][i_par]->Draw("hist same");
      c1->Print(Form("ValidPlot_%s.pdf",h_part_up[1][i_par]->GetTitle()));

      if(i_var_to_str.at(i_par).find(std::string("Elast")) != std::string::npos){
         h_part_up[1][i_par]->Divide(h_part_nom_el[1]);
         h_part_down[1][i_par]->Divide(h_part_nom_el[1]);
      }
      else{
         h_part_up[1][i_par]->Divide(h_part_nom_inel[1]);
         h_part_down[1][i_par]->Divide(h_part_nom_inel[1]);
      }

      h_part_up[1][i_par]->Draw("hist");
      h_part_up[1][i_par]->GetYaxis()->SetRangeUser(0.6,1.4);
      if(i_var_to_str.at(i_par).find(std::string("Elast")) != std::string::npos){
         h_part_nom_el_ratio->Draw("hist same");
      }
      else{
         h_part_nom_inel_ratio->Draw("hist same");
      }
      h_part_down[1][i_par]->Draw("hist same");
      c1->Print(Form("ValidPlot_Ratio_%s.pdf",h_part_up[1][i_par]->GetTitle()));
   }


   // ------------------------------------------------------ //
   //              Validation by event                       //
   // ------------------------------------------------------ //

   // Plan here is to check that the weights we get by event are the same as what we get when multiplying individual particle weights in the event

   // Fill by-event plots
   // n_entries = tmctr->GetEntries();
   // if (n_entries > 5000) n_entries = 5000;
   // for (int i_e=0; i_e<n_entries; i_e++){
   //    if (i_e%500==0) std::cout << i_e << "/" << n_entries << std::endl;
   //    tmctr->GetEntry(i_e);
   //
   //    if (e_inel_weight->size()>0){
   //       for (int i_truth=0; i_truth<e_inel_weight->size(); i_truth++){
   //          for (int i_univ=0; i_univ<e_inel_weight->at(i_truth).size(); i_univ++){
   //                h_evt_wgt[0][i_univ]->Fill(e_inel_weight->at(i_truth).at(i_univ));
   //             }
   //       }
   //    }
   //
   //    if (e_elastic_weight->size()>0){
   //       for (int i_truth=0; i_truth<e_elastic_weight->size(); i_truth++){
   //          for (int i_univ=0; i_univ<e_elastic_weight->at(i_truth).size(); i_univ++){
   //                h_evt_wgt[1][i_univ]->Fill(e_elastic_weight->at(i_truth).at(i_univ));
   //             }
   //       }
   //    }
   // } // end loop over events
   //
   // // Fill plots by event but calculated per track
   // run=-999;
   // subrun=-999;
   // event=-999;
   // std::vector<double> inel_weight(p_inel_weight->size(),1.0);
   // std::vector<double> elastic_weight(p_elastic_weight->size(),1.0);
   // for (int i_p=0; i_p<tpart->GetEntries(); i_p++){
   //    tpart->GetEntry(i_p);
   //
   //    if (p_run_num!=run || p_subrun_num!=subrun || p_event_num!=event){
   //       run = p_run_num;
   //       subrun = p_subrun_num;
   //       event = p_event_num;
   //
   //       // Fill histograms for the previous run/subrub/event and reset for the next
   //       for (int i_univ=0; i_univ<p_inel_weight->size(); i_univ++){
   //          h_evt_bypart_wgt[0][i_univ]->Fill(inel_weight.at(i_univ));
   //          inel_weight.at(i_univ) = 1.;
   //       }
   //       for (int i_univ=0; i_univ<p_elastic_weight->size(); i_univ++){
   //          h_evt_bypart_wgt[1][i_univ]->Fill(elastic_weight.at(i_univ));
   //          elastic_weight.at(i_univ) = 1.;
   //       }
   //    }
   //
   //    for (int i_univ=0; i_univ<p_inel_weight->size(); i_univ++){
   //       inel_weight.at(i_univ) *= p_inel_weight->at(i_univ);
   //    }
   //    for (int i_univ=0; i_univ<p_elastic_weight->size(); i_univ++){
   //       elastic_weight.at(i_univ) *= p_elastic_weight->at(i_univ);
   //    }
   // } // end loop over particles
   //
   // // Now draw comparisons
   // for (int i_univ=0; i_univ<UniverseVals->size(); i_univ++){
   //    h_evt_wgt[0][i_univ]->SetLineColor(kBlack);
   //    h_evt_wgt[0][i_univ]->SetLineWidth(2);
   //    h_evt_wgt[0][i_univ]->SetLineStyle(1);
   //    h_evt_bypart_wgt[0][i_univ]->SetLineColor(kBlue);
   //    h_evt_bypart_wgt[0][i_univ]->SetLineWidth(2);
   //    h_evt_bypart_wgt[0][i_univ]->SetLineStyle(2);
   //
   //    h_evt_wgt[1][i_univ]->SetLineColor(kBlack);
   //    h_evt_wgt[1][i_univ]->SetLineWidth(2);
   //    h_evt_wgt[1][i_univ]->SetLineStyle(1);
   //    h_evt_bypart_wgt[1][i_univ]->SetLineColor(kRed);
   //    h_evt_bypart_wgt[1][i_univ]->SetLineWidth(2);
   //    h_evt_bypart_wgt[1][i_univ]->SetLineStyle(2);
   //
   //    h_evt_wgt[0][i_univ]->Draw("hist");
   //    h_evt_bypart_wgt[0][i_univ]->Draw("hist same");
   //    c1->Print(Form("%s.pdf",h_evt_wgt[0][i_univ]->GetTitle()));
   //
   //    h_evt_wgt[1][i_univ]->Draw("hist");
   //    h_evt_bypart_wgt[1][i_univ]->Draw("hist same");
   //    c1->Print(Form("%s.pdf",h_evt_wgt[1][i_univ]->GetTitle()));
   // }
}

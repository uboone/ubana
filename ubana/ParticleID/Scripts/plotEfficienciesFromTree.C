#include "plotFromTreeHeader.h"

// What variables do we want these plots as a function of?
std::vector<std::vector<double>> GetPIDvarstoplot(treevars *vars){
  std::vector<std::vector<double>> varstoplot;
  for (size_t i=0; i<3; i++){
    varstoplot.push_back({
      vars->track_likelihood_p->at(i),
      vars->track_likelihood_mu->at(i),
      vars->track_likelihood_pi->at(i),
      vars->track_likelihood_k->at(i),
      vars->track_likelihood_mip->at(i),
      vars->track_likelihood_maxmumip->at(i),
      vars->track_chi2mu->at(i),
      vars->track_chi2p->at(i),
      vars->track_chi2pi->at(i),
      vars->track_chi2k->at(i),
      vars->track_PIDA_kde->at(i),
      vars->track_PIDA_median->at(i),
      vars->track_PIDA_mean->at(i),
      vars->track_likelihood_muoverp->at(i),
      vars->track_likelihood_mipoverp->at(i),
      vars->track_lnlikelihood_mipoverp->at(i),
      vars->track_likelihood_maxmumipoverp->at(i),
      vars->track_chi2_muminusp->at(i),
      vars->track_Lmu_0to1->at(i),
      vars->track_Lmip_0to1->at(i),
      vars->track_Lpi_0to1->at(i),
      vars->track_Lk_0to1->at(i),
      vars->track_Lp_0to1->at(i),
      vars->track_Lmumip_0to1->at(i),
      vars->track_Lmumippi_0to1->at(i),
      vars->track_Lmumip_0to1_nopionkaon->at(i),
      vars->track_Lmuovermip->at(i),
      vars->track_Lmumipoverpi->at(i),
      vars->track_depE_minus_rangeE_mu->at(i),
      vars->track_depE_minus_rangeE_p->at(i)
    });
  }

  // Add sum over all three planes
  std::vector<double> vars_sumplanes = varstoplot.at(0);
  for (size_t i=1; i<varstoplot.size(); i++){
    for (size_t ivar=0; ivar<vars_sumplanes.size(); ivar++){
      vars_sumplanes.at(ivar)+= varstoplot.at(i).at(ivar);
    }
  }
  // Normalise to number of planes (so we're using the average, not the sum)
  for (size_t ivar=0; ivar<vars_sumplanes.size(); ivar++){
    vars_sumplanes.at(ivar)/=varstoplot.size();
  }

  varstoplot.push_back(vars_sumplanes);

  return varstoplot;
};

// Binning (nbins, binlow, binhigh) in the same order as the vector above
std::vector<std::vector<double>> bins = {
                    {20,0,0.6}, // track_likelihood_p
                    {40,0,1.0}, // track_likelihood_mu
                    {40,0,1.0}, // track_likelihood_pi
                    {40,0,0.6}, // track_likelihood_k
                    {40,0,1.0}, // track_likelihood_mip
                    {40,0,1.0}, // track_likelihood_minmumip
                    {40,0,125}, // track_chi2mu
                    {50,0,400}, // track_chi2p
                    {40,0,125}, // track_chi2pi
                    {40,0,300}, // track_chi2k
                    {40,0,30}, // track_PIDA_kde
                    {40,0,30}, // track_PIDA_median
                    {40,0,30}, // track_PIDA_mean
                    {60,0,60}, // track_likelihood_muoverp
                    {60,0,60}, // track_likelihood_mipoverp
                    {60,-10,10}, // track_lnlikelihood_mipoverp
                    {60,0,60}, // track_likelihood_minmumipoverp
                    {50,-400,100}, // track_chi2_muminusp
                    {50,0,1}, // track_Lmu_0to1
                    {50,0,1}, // track_Lmip_0to1
                    {50,0,1}, // track_Lpi_0to1
                    {50,0,1}, // track_Lk_0to1
                    {50,0,1}, // track_Lp_0to1
                    {50,0,1}, // track_Lmumip_0to1
                    {50,0,1}, // track_Lmumippi_0to1
                    {50,0,1}, // track_Lmimip_0to1_nopionnokaon
                    {50,0,3}, // track_Lmuovermip
                    {50,0,3}, // track_Lmumipoverpi
                    {50,-100,150}, // track_depE_minus_rangeE_mu
                    {50,-100,100} // track_depE_minus_rangeE_p
                    };

// Histogram titles in the same order as the vector above
std::vector<std::string> histtitles = {
                    ";L_{p};",
                    ";L_{#mu};",
                    ";L_{#pi};",
                    ";L_{K};",
                    ";L_{MIP};",
                    ";L_{#mu/MIP};",
                    ";#chi^{2}_{#mu};",
                    ";#chi^{2}_{p};",
                    ";#chi^{2}_{#pi};",
                    ";#chi^{2}_{K};",
                    ";PIDa (by KDE);",
                    ";PIDa (by median);",
                    ";PIDa (by mean);",
                    ";(L_{#mu})/(L_{p});",
                    ";(L_{MIP})/(L_{p});",
                    ";ln(L_{MIP}/L_{p})",
                    ";(L_{#mu/MIP})/(L_{p});",
                    ";#chi^{2}_{#mu}-#chi^{2}_{p};",
                    ";L_{#mu}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{MIP}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{#pi}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{K}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{p}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP}+L_{#pi})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p});",
                    ";(L_{#mu}/L_{MIP});",
                    ";L_{#mu/MIP}/L_{#pi};",
                    ";Dep. E - E by range (muon assumption) [MeV];",
                    ";Dep. E - E by range (proton assumption) [MeV];"
                  };

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
                  "effpur_Lp",
                  "effpur_Lmu",
                  "effpur_Lpi",
                  "effpur_Lk",
                  "effpur_Lmip",
                  "effpur_Lmumip",
                  "effpur_chi2mu",
                  "effpur_chi2p",
                  "effpur_chi2pi",
                  "effpur_chi2k",
                  "effpur_pida_kde",
                  "effpur_pida_median",
                  "effpur_pida_mean",
                  "effpur_Lmuoverp",
                  "effpur_Lmipoverp",
                  "effpur_lnLmipoverp",
                  "effpur_Lmumipoverp",
                  "effpur_chi2muminusp",
                  "effpur_Lmu0to1",
                  "effpur_Lmip0to1",
                  "effpur_Lpi0to1",
                  "effpur_Lk0to1",
                  "effpur_Lp0to1",
                  "effpur_Lmumip0to1",
                  "effpur_Lmumippi0to1",
                  "effpur_Lmumip0to1nopionkaon",
                  "effpur_Lmuovermip",
                  "effpur_Lmumipoverpi",
                  "effpur_depErangeEmu",
                  "effpur_depErangeEp"
                };

// For efficiency/purity we need to know whether MIPs are supposed to be low or high. In the same order as the vector above
std::vector<bool> MIPlow = {
                    true, // track_likelihood_p
                    false, // track_likelihood_mu
                    false, // track_likelihood_pi
                    true, // track_likelihood_k
                    false, // track_likelihood_mip
                    false, // track_likelihood_minmumip
                    true, // track_chi2mu
                    false, // track_chi2p
                    true, // track_chi2pi
                    false, // track_chi2k
                    true, // track_PIDA_kde
                    true, // track_PIDA_median
                    true, // track_PIDA_mean
                    false, // track_likelihood_muminusp
                    false, // track_likelihood_mipminusp
                    false, // track_lnlikelihood_mipoverp
                    false, // track_likelihood_minmumipminusp
                    true, // track_chi2_muminusp
                    false, // track_Lmu_0to1
                    false, // track_Lmip_0to1
                    false, // track_Lpi_0to1
                    true, // track_Lk_0to1
                    true, // track_Lp_0to1
                    false, // track_Lmumip_0to1
                    false, // track_Lmumippi_0to1
                    false, // track_Lmumip0to1nopionkaon
                    false, // track_Lmuovermip
                    true, // track_Lmumipoverpi
                    true, // track_depE_minus_rangeE_mu
                    true // track_depE_minus_rangeE_p
                    };

// cut values for evaluating 2D efficiency/purity plots in terms of track length and theta
std::vector<double> cutValues = {
                    -999, // track_likelihood_p
                    -999, // track_likelihood_mu
                    -999, // track_likelihood_pi
                    -999, // track_likelihood_k
                    -999, // track_likelihood_mip
                    -999, // track_likelihood_minmumip
                    -999, // track_chi2mu
                    88.0, // track_chi2p
                    -999, // track_chi2pi
                    -999, // track_chi2k
                    -999, // track_PIDA_kde
                    -999, // track_PIDA_median
                    10,// 12.5, // track_PIDA_mean
                    -999, // track_likelihood_muminusp
                    2.0, // track_likelihood_mipminusp
                    1.0, // track_lnlikelihood_mipoverp
                    -999, // track_likelihood_minmumipminusp
                    -90.0, // track_chi2_muminusp
                    -999, // track_Lmu_0to1
                    -999, // track_Lmip_0to1
                    -999, // track_Lpi_0to1
                    -999, // track_Lk_0to1
                    -999, // track_Lp_0to1
                    -999, // track_Lmumip_0to1
                    -999, // track_Lmumippi_0to1
                    0.68, // track_Lmumip0to1nopionkaon
                    -999, // track_Lmuovermip
                    -999, // track_Lmumipoverpi
                    0.0, // track_depE_minus_rangeE_mu
                    -999 // track_depE_minus_rangeE_p
};

std::vector<string> twodEffPurVars = {
                    "track_length",
                    "track_theta"
};
std::vector<std::vector<double>> twodEffPurVars_bins = {
                    {50, 0, 700},
                    {50, 0, 3.15}
};

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void plotEfficienciesFromTree(std::string mcfile, double POTscaling=0., std::string onbeamdatafile="", std::string offbeamdatafile="", double offbeamscaling=0., bool onminusoffbeam=true, std::string outfile=NULL){

  // Make a file to store output
  TFile *fout = nullptr;
  if (outfile!=NULL){
    fout = new TFile(outfile.c_str(),"recreate");
  }

  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0.);

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get("pidvalid/pidTree");
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  TFile *f_onbeam=nullptr;
  TTree *t_onbeam=nullptr;
  treevars onbeam_vars;
  if (onbeamdatafile!=""){
    std::cout << "Making data-MC comparisons" << std::endl;
    f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
    t_onbeam = (TTree*)f_onbeam->Get("pidvalid/pidTree");
    settreevars(t_onbeam,&onbeam_vars);
  }

  TFile *f_offbeam=nullptr;
  TTree *t_offbeam=nullptr;
  treevars offbeam_vars;
  if (offbeamdatafile!=""){
    f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
    t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree");
    settreevars(t_offbeam,&offbeam_vars);
  }

  // Sanity check: the plot vectors should be the same size
  t_bnbcos->GetEntry(0);
  CalcPIDvars(&mc_vars, false);
  std::vector<std::vector<double>> PIDvarstoplot_dummy = GetPIDvarstoplot(&mc_vars);
  // if (PIDvarstoplot_dummy.size() != bins.size()) std::cout << "WARNING PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;
  std::cout << "PIDvarstoplot.size() = " << PIDvarstoplot_dummy.size() << std::endl;
  std::cout << "PIDvarstoplot.at(0).size() = " << PIDvarstoplot_dummy.at(0).size() << std::endl;
  std::cout << "bins.size() = " << bins.size() << std::endl;
  std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
  std::cout << "histnames.size() = " << histnames.size() << std::endl;
  std::cout << "MIPlow.size() = " << MIPlow.size() << std::endl;

  // ----------------- MC

  // set branch addresses for variables of interest
  double twodvar1 = -999;
  double twodvar2 = -999;

  // t_bnbcos->SetBranchAddress(twodEffPurVars.at(0).c_str(), &twodvar1);
  // t_bnbcos->SetBranchAddress(twodEffPurVars.at(1).c_str(), &twodvar2);

  // Make histograms to fill
  const size_t nplanes = PIDvarstoplot_dummy.size();
  const size_t nplots = PIDvarstoplot_dummy.at(0).size();
  hist1D *mc_hists[nplanes][nplots];
  hist1D *mc_passingfrac_hists[nplanes][nplots];
  hist2D *mc_2dhists[nplanes][nplots];
  for (int i_pl=0; i_pl<nplanes; i_pl++){
    for (int i_h=0; i_h<nplots; i_h++){
      mc_hists[i_pl][i_h] = new hist1D(std::string("h_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));

      mc_passingfrac_hists[i_pl][i_h] = new hist1D(std::string("h_pf_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),2,0,2);

      mc_2dhists[i_pl][i_h] = new hist2D(std::string("h2d_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h)+";"+twodEffPurVars.at(0)+";"+twodEffPurVars.at(1), twodEffPurVars_bins.at(0).at(0), twodEffPurVars_bins.at(0).at(1), twodEffPurVars_bins.at(0).at(2), twodEffPurVars_bins.at(1).at(0), twodEffPurVars_bins.at(1).at(1), twodEffPurVars_bins.at(1).at(2), bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }
  }

  // Loop through MC tree and fill plots
  for (int i = 0; i < t_bnbcos->GetEntries(); i++){
    t_bnbcos->GetEntry(i);
    CalcPIDvars(&mc_vars, false);
    std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&mc_vars);

    if (!(mc_vars.track_theta_x > 0 && mc_vars.track_theta_x < 45)) continue;

    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (size_t i_h = 0; i_h < nplots; i_h++){

        FillHist(mc_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),mc_vars.true_PDG);

        Fill2DHist(mc_2dhists[i_pl][i_h], twodvar1, twodvar2, PIDvarstoplot.at(i_pl).at(i_h), mc_vars.true_PDG);

        double passes = 0;
        if (PIDvarstoplot.at(i_pl).at(i_h)>cutValues.at(i_h)) passes = 1;
        FillHist(mc_passingfrac_hists[i_pl][i_h],passes,mc_vars.true_PDG);
      }
    }



  } // end loop over entries in tree

  // ----------------- On-beam data

  // Make histograms to fill
  hist1D *onb_hists[nplanes][nplots];
  if (t_onbeam){
    for (int i_pl=0; i_pl<nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){
        onb_hists[i_pl][i_h] = new hist1D(std::string("h_ondat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),2,0,2);
      }
    }

    // Loop through on-beam data tree and fill plots
    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);
      CalcPIDvars(&onbeam_vars, false);
      std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&onbeam_vars);

      if (!(onbeam_vars.track_theta_x > 0 && onbeam_vars.track_theta_x < 45)) continue;

      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        for (size_t i_h = 0; i_h < nplots; i_h++){
          double passes = 0;
          if (PIDvarstoplot.at(i_pl).at(i_h)>cutValues.at(i_h)) passes = 1;
          FillHist(onb_hists[i_pl][i_h],passes,0);
        }
      }
    } // end loop over entries in tree
  } // end if (t_onbeam)

  // ----------------- Off-beam data

  // Make histograms to fill
  hist1D *offb_hists[nplanes][nplots];
  if (t_offbeam){
    for (int i_pl=0; i_pl<nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){
        offb_hists[i_pl][i_h] = new hist1D(std::string("h_offdat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),2,0,2);
      }
    }

    // Loop through on-beam data tree and fill plots
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
      t_offbeam->GetEntry(i);
      CalcPIDvars(&offbeam_vars, false);
      std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&offbeam_vars);

      if (!(offbeam_vars.track_theta_x > 0 && offbeam_vars.track_theta_x < 45)) continue;

      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        for (size_t i_h = 0; i_h < nplots; i_h++){
          double passes = 0;
          if (PIDvarstoplot.at(i_pl).at(i_h)>cutValues.at(i_h)) passes = 1;
          FillHist(offb_hists[i_pl][i_h],passes,0);
        }
      }
    } // end loop over entries in tree
  } // end if (t_onbeam)


  // -------------------- Now make all the plots

  for (size_t i_pl=0; i_pl < nplanes; i_pl++){
    for (size_t i_h=0; i_h < nplots; i_h++){
      // Efficiency and purity: MC only
      // -- 1D
      TCanvas *c1 = new TCanvas();
      DrawMCEffPur(c1, mc_hists[i_pl][i_h],MIPlow.at(i_h),cutValues.at(i_h),fout);
      c1->Print(std::string(histnames[i_h]+std::string("_plane")+std::to_string(i_pl)+".png").c_str());
      delete c1;
      // -- 2D
      TCanvas *c2 = new TCanvas();
      TCanvas *c3 = new TCanvas();
      TCanvas *c4 = new TCanvas();
      DrawMCEffPur2D(c2, c3, c4, mc_2dhists[i_pl][i_h], MIPlow.at(i_h), cutValues.at(i_h), twodvar1, twodvar2, twodEffPurVars, fout);
      if (cutValues.at(i_h) != -999){
        c2->Print(std::string("eff_"+histnames[i_h]+std::string("_plane")+"_"+twodEffPurVars.at(0)+"_"+twodEffPurVars.at(1)+"_"+std::to_string(i_pl)+".png").c_str());
        c3->Print(std::string("pur_"+histnames[i_h]+std::string("_plane")+"_"+twodEffPurVars.at(0)+"_"+twodEffPurVars.at(1)+"_"+std::to_string(i_pl)+".png").c_str());
        c4->Print(std::string("effpur_"+histnames[i_h]+std::string("_plane")+"_"+twodEffPurVars.at(0)+"_"+twodEffPurVars.at(1)+"_"+std::to_string(i_pl)+".png").c_str());
      }
      delete c2;
      delete c3;
      delete c4;

      // Now plot passing fractions: can be MC and data
      TCanvas *c5 = new TCanvas();
      if (onminusoffbeam){
        DrawMC(mc_passingfrac_hists[i_pl][i_h],POTscaling,-999);
        if (f_onbeam && f_offbeam){
          OverlayOnMinusOffData(c5,onb_hists[i_pl][i_h],offb_hists[i_pl][i_h],offbeamscaling,POTscaling);
          TString e_str("h_err");
          TString o_str("h_ondat_"+histnames[i_h]+"_plane"+std::to_string(i_pl)+"_all");
          OverlayChi2(c5, e_str, o_str);
        }
      }

      else{
        if (f_onbeam && f_offbeam){
          DrawMCPlusOffbeam(mc_passingfrac_hists[i_pl][i_h],offb_hists[i_pl][i_h],POTscaling,offbeamscaling,-999);
          OverlayOnBeamData(c5,onb_hists[i_pl][i_h]);
          TString e_str("h_err");
          TString o_str("h_ondat_"+histnames[i_h]+"_plane"+std::to_string(i_pl)+"_all");
          OverlayChi2(c5, e_str, o_str);
        }
        else{
          DrawMC(mc_hists[i_pl][i_h],POTscaling,-999);
        }
      }
      if (cutValues.at(i_h)!=-999)
        c5->Print(std::string(std::string("PassingFrac_")+histnames[i_h]+std::string("_plane")+std::to_string(i_pl)+".png").c_str());
      delete c5;
    }
  }
  if (outfile!=NULL) fout->Close();
  delete f_bnbcos;
  if (offbeamdatafile!="") delete f_offbeam;
  if (onbeamdatafile!="") delete f_onbeam;
}

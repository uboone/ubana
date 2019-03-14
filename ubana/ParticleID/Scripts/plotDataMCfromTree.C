#include "/uboone/app/users/afurmans/CCInclusive_PID_Jan/srcs/ubobj/ubobj/UBXSec/UBXSecEvent.h"
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
      vars->track_PIDA_mean->at(i),
      vars->track_PIDA_median->at(i),
      vars->track_likelihood_muoverp->at(i),
      vars->track_likelihood_mipoverp->at(i),
      vars->track_lnlikelihood_mipoverp->at(i),
      vars->track_likelihood_maxmumipoverp->at(i),
      vars->track_chi2_muminusp->at(i),
      vars->track_Lmu_0to1->at(i),
      vars->track_Lmip_0to1->at(i),
      vars->track_Lpi_0to1->at(i),
      vars->track_Lp_0to1->at(i),
      vars->track_Lk_0to1->at(i),
      vars->track_Lmumip_0to1->at(i),
      vars->track_Lmumippi_0to1->at(i),
      vars->track_Lmumip_0to1_nopionkaon->at(i),
      vars->track_Lmumip_0to1_nopionkaon->at(i),
      vars->track_Lmuovermip->at(i),
      vars->track_Lmumipoverpi->at(i),
      vars->track_depE_minus_rangeE_mu->at(i),
      vars->track_depE_minus_rangeE_p->at(i),
      vars->track_Lmip_atstart->at(i),
      vars->track_lnLmip_atstart->at(i),
      vars->track_dEdx_mean_atstart->at(i),
      vars->track_shift_bestL->at(i)
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
                    {80,0,1.0}, // track_likelihood_mip
                    {40,0,1.0}, // track_likelihood_minmumip
                    {40,0,125}, // track_chi2mu
                    {50,0,400}, // track_chi2p
                    {40,0,125}, // track_chi2pi
                    {40,0,300}, // track_chi2k
                    {40,0,30}, // track_PIDA_kde
                    {40,0,30}, // track_PIDA_mean
                    {40,0,30}, // track_PIDA_median
                    {60,0,60}, // track_likelihood_muoverp
                    {60,0,60}, // track_likelihood_mipoverp
                    {60,-10,10}, // track_lnlikelihood_mipoverp
                    {60,0,60}, // track_likelihood_minmumipoverp
                    {50,-400,100}, // track_chi2_muminusp
                    {50,0,1}, // track_Lmu_0to1
                    {50,0,1}, // track_Lmip_0to1
                    {50,0,1}, // track_Lpi_0to1
                    {50,0,1}, // track_Lp_0to1
                    {50,0,1}, // track_Lk_0to1
                    {50,0,1}, // track_Lmumip_0to1
                    {50,0,1}, // track_Lmumippi_0to1
                    {50,0,1}, // track_Lmumip_0to1_nopionkaon
                    {50,0,1}, // track_Lmumip_0to1_nopionkaon_zoom
                    {50,0,3}, // track_Lmuovermip,
                    {50,0,3}, // track_Lmumipoverpi,
                    {50,-150,150}, // track_depE_minus_rangeE_mu
                    {50,-300,100}, // track_depE_minus_rangeE_p
                    {100,0,2}, // track_Lmip_atstart
                    {100,-10,10}, // track_lnLmip_atstart
                    {50,0,10}, // track_dEdx_mean_atstart
                    {20,-2,2} // track_shift_bestL
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
                    ";PIDa (by mean);",
                    ";PIDa (by median);",
                    ";(L_{#mu})/(L_{p});",
                    ";(L_{MIP})/(L_{p});",
                    ";ln(L_{MIP}/L_{p});",
                    ";(L_{#mu/MIP})/(L_{p});",
                    ";#chi^{2}_{#mu}-#chi^{2}_{p};",
                    ";L_{#mu}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{MIP}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{#pi}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{p}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{k}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP}+L_{#pi})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p});",
                    ";(L_{#mu}/L_{MIP});",
                    ";(L_{#mu/MIP}/L_{#pi});",
                    ";Dep. E - E. by range (muon assumption) [MeV];",
                    ";Dep. E - E. by range (proton assumption) [MeV];",
                    ";L_{MIP} at start of track;",
                    ";ln(L_{MIP}) at start of track;",
                    ";Truncated Mean dE/dx at start of track;",
                    ";Preferred shift for maximum likelihood [cm];"
                  };
//
// // What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
                  "Lp",
                  "Lmu",
                  "Lpi",
                  "Lk",
                  "Lmip",
                  "Lmumip",
                  "chi2mu",
                  "chi2p",
                  "chi2pi",
                  "chi2k",
                  "pida_kde",
                  "pida_mean",
                  "pida_median",
                  "Lmuoverp",
                  "Lmipoverp",
                  "lnLmipoverp",
                  "Lmumipoverp",
                  "chi2muminusp",
                  "Lmu0to1",
                  "Lmip0to1",
                  "Lpi0to1",
                  "Lp0to1",
                  "Lk0to1",
                  "Lmumip0to1",
                  "Lmumippi0to1",
                  "Lmumip0to1nopionkaon",
                  "Lmumip0to1nopionkaon_zoom",
                  "Lmuovermip",
                  "Lmumipoverpi",
                  "depErangeEmu",
                  "depErangeEp",
                  "Lmip_atstart",
                  "lnLmip_atstart",
                  "dEdx_truncmean_atstart",
                  "shift_bestL"
                };

// Set y-axis range to zoom in if we need to
// -999 means no zoom
std::vector<double> yrange = {
                  -999, // track_likelihood_p
                  -999, // track_likelihood_mu
                  -999, // track_likelihood_pi
                  -999, // track_likelihood_k
                  -999, // track_likelihood_mip
                  -999, // track_likelihood_minmumip
                  -999, // track_chi2mu
                  -999, // track_chi2p
                  -999, // track_chi2pi
                  -999, // track_chi2k
                  -999, // track_PIDA_kde
                  -999, // track_PIDA_mean
                  -999, // track_PIDA_median
                  -999, // track_likelihood_muoverp
                  -999, // track_likelihood_mipoverp
                  -999, // track_lnlikelihood_mipoverp
                  -999, // track_likelihood_minmumipoverp
                  -999, // track_chi2_muminusp
                  -999, // track_Lmu_0to1
                  -999, // track_Lmip_0to1
                  -999, // track_Lpi_0to1
                  -999, // track_Lp_0to1
                  -999, // track_Lk_0to1
                  -999, // track_Lmumip_0to1
                  -999, // track_Lmumippi_0to1
                  -999, // track_Lmumip_0to1_nopionkaon
                  1000, // track_Lmumip_0to1_nopionkaon_zoom
                  -999, // track_Lmuovermip,
                  -999, // track_Lmumipoverpi,
                  -999, // track_depE_minus_rangeE_mu
                  -999, // track_depE_minus_rangeE_p
                  -999, // track_Lmip_atstart
                  -999, // track_lnLmip_atstart
                  -999, // track_dEdx_mean_atstart
                  -999 // track_shift_bestL
                };

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void plotDataMCfromTree(std::string treename, std::string mcfile, double POTscaling=0., std::string onbeamdatafile="", std::string offbeamdatafile="", double offbeamscaling=0., bool onminusoffbeam=true, bool templatefit=false){
//void plotDataMCfromTree(){ // in case you want to hard-code things, you can uncomment this block
//  std::string treename="pidvalidcaliSCE/pidTree";
//  std::string mcfile="/uboone/data/users/sfehlber/Mar/Mar11/mc/mc.root"; // MC pure
////  std::string mcfile="/uboone/data/users/sfehlber/Mar/Mar11_test/overlay1/overlay1.root"; // Overlay
//  std::string onbeamdatafile="/uboone/data/users/sfehlber/Mar/Mar11/bnb/bnb.root"; // onbeam
//  std::string offbeamdatafile="/uboone/data/users/sfehlber/Mar/Mar11/extbnb/ext.root"; // onbeam
//  double offbeamscaling=0.647;
//  double POTscaling =1.118; //MC pure
//  //double POTscaling =6.202; Overlay
//  bool onminusoffbeam=false;
//  bool templatefit=false;


  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0.);
  gStyle->SetOptStat(0);

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get(treename.c_str());
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  TFile *f_onbeam=nullptr;
  TTree *t_onbeam=nullptr;
  treevars onbeam_vars;
  if (onbeamdatafile!=""){
    std::cout << "Making data-MC comparisons" << std::endl;
    f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
    t_onbeam = (TTree*)f_onbeam->Get(treename.c_str());
    settreevars(t_onbeam,&onbeam_vars);
  }

  TFile *f_offbeam=nullptr;
  TTree *t_offbeam=nullptr;
  treevars offbeam_vars;
  if (offbeamdatafile!=""){
    f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
    t_offbeam = (TTree*)f_offbeam->Get(treename.c_str());
    settreevars(t_offbeam,&offbeam_vars);
  }


  // Sanity check: the plot vectors should be the same size
  t_bnbcos->GetEntry(0);
  CalcPIDvars(&mc_vars, true);
  std::vector<std::vector<double>> PIDvarstoplot_dummy = GetPIDvarstoplot(&mc_vars);
  if (PIDvarstoplot_dummy.size() != 3 && PIDvarstoplot_dummy.size() != 4) std::cout << "WARNING PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.size() << ", should be 3 or 4." << std::endl;
  if (PIDvarstoplot_dummy.at(0).size() != bins.size()) std::cout << "WARNING PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;
  std::cout << "PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.at(0).size() << std::endl;
  std::cout << "bins.size() = " << bins.size() << std::endl;
  std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
  std::cout << "histnames.size() = " << histnames.size() << std::endl;
  std::cout << "yrange.size() = " << yrange.size() << std::endl;

  // ----------------- MC

  // Make histograms to fill
  const size_t nplanes = PIDvarstoplot_dummy.size();
  const size_t nplots = PIDvarstoplot_dummy.at(0).size();
  hist1D *mc_hists[nplanes][nplots];
  hist1D *mc_hists_MIPdEdx[nplanes];
  hist1D *mc_hists_dEdxSlices[nplanes][10];
  for (int i_pl=0; i_pl<nplanes; i_pl++){
    for (int i_h=0; i_h<nplots; i_h++){
      mc_hists[i_pl][i_h] = new hist1D(std::string("h_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }
    for (int slice(0); slice <= 10; slice++){
      mc_hists_dEdxSlices[i_pl][slice] = new hist1D("h_mc_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,7);
    }
    // Extra custom plots
    mc_hists_MIPdEdx[i_pl] = new hist1D(std::string("h_")+std::string("MIPregiondEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
  }


  // Loop through MC tree and fill plots
  for (int i = 0; i < t_bnbcos->GetEntries(); i++){
    t_bnbcos->GetEntry(i);
    CalcPIDvars(&mc_vars, true);
    std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&mc_vars);
    
    
    // if (mc_vars.track_theta_x > 75 && mc_vars.track_theta_x < 90){
      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        for (size_t i_h = 0; i_h < nplots; i_h++){
          FillHist(mc_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),mc_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<mc_vars.track_resrange_perhit_u->size(); i_hit++){
        if (mc_vars.track_resrange_perhit_u->at(i_hit)>100 && mc_vars.track_resrange_perhit_u->at(i_hit)<150){
          FillHist(mc_hists_MIPdEdx[0],mc_vars.track_dEdx_perhit_u->at(i_hit),mc_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<mc_vars.track_resrange_perhit_v->size(); i_hit++){
        if (mc_vars.track_resrange_perhit_v->at(i_hit)>100 && mc_vars.track_resrange_perhit_v->at(i_hit)<150){
          FillHist(mc_hists_MIPdEdx[1],mc_vars.track_dEdx_perhit_v->at(i_hit),mc_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<mc_vars.track_resrange_perhit_y->size(); i_hit++){
        if (mc_vars.track_resrange_perhit_y->at(i_hit)>100 && mc_vars.track_resrange_perhit_y->at(i_hit)<150){
          FillHist(mc_hists_MIPdEdx[2],mc_vars.track_dEdx_perhit_y->at(i_hit),mc_vars.true_PDG);
        }
      }
    // }
    double rr;
    double dedx;
    int slice;
    int nhits;
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      if (i_pl==0){
        nhits=mc_vars.track_resrange_perhit_u->size();}
      if (i_pl==1){
        nhits=mc_vars.track_resrange_perhit_v->size();}
      if (i_pl==2){
        nhits=mc_vars.track_resrange_perhit_y->size();}

      for (size_t i_hit=0; i_hit<nhits; i_hit++){
        if (i_pl==0){
          rr = mc_vars.track_resrange_perhit_u->at(i_hit);
          dedx = mc_vars.track_dEdx_perhit_u->at(i_hit);
          slice=rr/5;
        }
        if (i_pl==1){
          rr = mc_vars.track_resrange_perhit_v->at(i_hit);
          dedx = mc_vars.track_dEdx_perhit_v->at(i_hit);
          slice=rr/5;
        }
        if (i_pl==2){
          rr = mc_vars.track_resrange_perhit_y->at(i_hit);
          dedx = mc_vars.track_dEdx_perhit_y->at(i_hit);
          slice=rr/5;
        }
        if (slice>10){slice=10;}
        if (slice<0){continue;}
        FillHist(mc_hists_dEdxSlices[i_pl][slice],dedx,mc_vars.true_PDG);
      }
    }

  } // end loop over entries in tree
  std::cout << "Done MC, now on-beam" << std::endl;
  // ----------------- On-beam data
  hist1D *onb_hists[nplanes][nplots];
  hist1D *onb_hists_MIPdEdx[nplanes];
  hist1D *onb_hists_dEdxSlices[nplanes][10];
  if (t_onbeam){
    // Make histograms to fill
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){
        onb_hists[i_pl][i_h] = new hist1D(std::string("h_ondat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      }
      for (int slice(0); slice <=10; slice++){
        onb_hists_dEdxSlices[i_pl][slice] = new hist1D("h_onb_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,7);
      }
      // Extra custom plots
      onb_hists_MIPdEdx[i_pl] = new hist1D(std::string("h_ondat_")+std::string("MIPregiondEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
    }

    // Loop through on-beam data tree and fill plots
    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);
      CalcPIDvars(&onbeam_vars, false);
      std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&onbeam_vars);

      if (onbeam_vars.track_theta_x > 75 && onbeam_vars.track_theta_x < 90){
        for (size_t i_pl=0; i_pl < nplanes; i_pl++){
          for (size_t i_h = 0; i_h < nplots; i_h++){
            FillHist(onb_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),0); // 0 because there is no "true PDG" for data
          }
        }
      }

      for (size_t i_hit=0; i_hit<onbeam_vars.track_resrange_perhit_u->size(); i_hit++){
        if (onbeam_vars.track_resrange_perhit_u->at(i_hit)>100 && onbeam_vars.track_resrange_perhit_u->at(i_hit)<150){
          FillHist(onb_hists_MIPdEdx[0],onbeam_vars.track_dEdx_perhit_u->at(i_hit),0);
        }
      }

      for (size_t i_hit=0; i_hit<onbeam_vars.track_resrange_perhit_v->size(); i_hit++){
        if (onbeam_vars.track_resrange_perhit_v->at(i_hit)>100 && onbeam_vars.track_resrange_perhit_v->at(i_hit)<150){
          FillHist(onb_hists_MIPdEdx[1],onbeam_vars.track_dEdx_perhit_v->at(i_hit),0);
        }
      }

      for (size_t i_hit=0; i_hit<onbeam_vars.track_resrange_perhit_y->size(); i_hit++){
        if (onbeam_vars.track_resrange_perhit_y->at(i_hit)>100 && onbeam_vars.track_resrange_perhit_y->at(i_hit)<150){
          FillHist(onb_hists_MIPdEdx[2],onbeam_vars.track_dEdx_perhit_y->at(i_hit),0);
        }
      }
      double rr;
      double dedx;
      int slice;
      int nhits;
      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        if (i_pl==0){
          nhits=onbeam_vars.track_resrange_perhit_u->size();}
        if (i_pl==1){
          nhits=onbeam_vars.track_resrange_perhit_v->size();}
        if (i_pl==2){
          nhits=onbeam_vars.track_resrange_perhit_y->size();}
        for (size_t i_hit=0; i_hit<nhits; i_hit++){
          if (i_pl==0){
            rr = onbeam_vars.track_resrange_perhit_u->at(i_hit);
            dedx = onbeam_vars.track_dEdx_perhit_u->at(i_hit);
            slice=rr/5;
            if (slice>10){slice=10;}
          }
          if (i_pl==1){
            rr = onbeam_vars.track_resrange_perhit_v->at(i_hit);
            dedx = onbeam_vars.track_dEdx_perhit_v->at(i_hit);
            slice=rr/5;
            if (slice>10){slice=10;}
          }
          if (i_pl==2){
            rr = onbeam_vars.track_resrange_perhit_y->at(i_hit);
            dedx = onbeam_vars.track_dEdx_perhit_y->at(i_hit);
            slice=rr/5;
            if (slice>10){slice=10;}
          }
          FillHist(onb_hists_dEdxSlices[i_pl][slice],dedx,0);
        }
      }
    }
  }
    // ----------------- Off-beam data
  hist1D *offb_hists[nplanes][nplots];
  hist1D *offb_hists_MIPdEdx[nplanes];
  hist1D *offb_hists_dEdxSlices[nplanes][10];
  if (t_offbeam){
    // Make histograms to fill
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){
        offb_hists[i_pl][i_h] = new hist1D(std::string("h_offdat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      }
      for (int slice(0); slice <=10; slice++){
        offb_hists_dEdxSlices[i_pl][slice] = new hist1D("h_offb_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,7);
      }
      // Extra custom plots
      offb_hists_MIPdEdx[i_pl] = new hist1D(std::string("h_offdat_")+std::string("MIPregiondEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
    }


    // Loop through tree and fill plots
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
      t_offbeam->GetEntry(i);
      CalcPIDvars(&offbeam_vars, false);
      std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&offbeam_vars);

      if (offbeam_vars.track_theta_x >75 && offbeam_vars.track_theta_x < 90){
        for (size_t i_pl=0; i_pl < nplanes; i_pl++){
          for (size_t i_h = 0; i_h < nplots; i_h++){
            FillHist(offb_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),0.); // 0 because there is no "true PDG" for data
          }
        }
      }

      for (size_t i_hit=0; i_hit<offbeam_vars.track_resrange_perhit_u->size(); i_hit++){
        if (offbeam_vars.track_resrange_perhit_u->at(i_hit)>100 && offbeam_vars.track_resrange_perhit_u->at(i_hit)<150){
          FillHist(offb_hists_MIPdEdx[0],offbeam_vars.track_dEdx_perhit_u->at(i_hit),offbeam_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<offbeam_vars.track_resrange_perhit_v->size(); i_hit++){
        if (offbeam_vars.track_resrange_perhit_v->at(i_hit)>100 && offbeam_vars.track_resrange_perhit_v->at(i_hit)<150){
          FillHist(offb_hists_MIPdEdx[1],offbeam_vars.track_dEdx_perhit_v->at(i_hit),offbeam_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<offbeam_vars.track_resrange_perhit_y->size(); i_hit++){
        if (offbeam_vars.track_resrange_perhit_y->at(i_hit)>100 && offbeam_vars.track_resrange_perhit_y->at(i_hit)<150){
          FillHist(offb_hists_MIPdEdx[2],offbeam_vars.track_dEdx_perhit_y->at(i_hit),offbeam_vars.true_PDG);
        }
      }
      double rr;
      double dedx;
      int slice;
      int nhits;
      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        if (i_pl==0){
          nhits=offbeam_vars.track_resrange_perhit_u->size();}
        if (i_pl==1){
          nhits=offbeam_vars.track_resrange_perhit_v->size();}
        if (i_pl==2){
          nhits=offbeam_vars.track_resrange_perhit_y->size();}

        for (size_t i_hit=0; i_hit<nhits; i_hit++){
          if (i_pl==0){
            rr = offbeam_vars.track_resrange_perhit_u->at(i_hit);
            dedx = offbeam_vars.track_dEdx_perhit_u->at(i_hit);
            slice=rr/5;
            if (slice>10){slice=10;}
          }
          if (i_pl==1){
            rr = offbeam_vars.track_resrange_perhit_v->at(i_hit);
            dedx = offbeam_vars.track_dEdx_perhit_v->at(i_hit);
            slice=rr/5;
            if (slice>10){slice=10;}
          }
          if (i_pl==2){
            rr = offbeam_vars.track_resrange_perhit_y->at(i_hit);
            dedx = offbeam_vars.track_dEdx_perhit_y->at(i_hit);
            slice=rr/5;
            if (slice>10){slice=10;}
          }
          FillHist(offb_hists_dEdxSlices[i_pl][slice],dedx,0);
        }
      }
    } // end loop over entries in tree
  }


  // -------------------- Now make all the plots

  for (size_t i_pl=0; i_pl < nplanes; i_pl++){
    for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas();

      double POTscaling_tmp = POTscaling; // Reset POT scaling for the next plot

      if (templatefit){
        TemplateFit(mc_hists[i_pl][i_h], onb_hists[i_pl][i_h], offb_hists[i_pl][i_h], offbeamscaling, POTscaling_tmp,yrange.at(i_h));
        POTscaling_tmp = 1.; // We have already applied POT scaling to MC in the template fit function. Set it to 1 so it doesn't get applied again in DrawMC.
      }

      if (onminusoffbeam){
        DrawMC(mc_hists[i_pl][i_h],POTscaling_tmp,yrange.at(i_h));
        if (f_onbeam && f_offbeam){
          OverlayOnMinusOffData(c1,onb_hists[i_pl][i_h],offb_hists[i_pl][i_h],offbeamscaling,POTscaling_tmp);
          TString e_str("h_err");
          TString o_str("h_ondat_"+histnames[i_h]+"_plane"+std::to_string(i_pl)+"_all");
          OverlayChi2(c1, e_str, o_str);
        }
      }
      else{
        if (f_onbeam && f_offbeam){
          DrawMCPlusOffbeam(mc_hists[i_pl][i_h], offb_hists[i_pl][i_h], POTscaling_tmp, offbeamscaling,yrange.at(i_h));
          OverlayOnBeamData(c1, onb_hists[i_pl][i_h]);
          TString e_str("h_err");
          TString o_str("h_ondat_"+histnames[i_h]+"_plane"+std::to_string(i_pl)+"_all");
          OverlayChi2(c1, e_str, o_str);
        }
        else{
          DrawMC(mc_hists[i_pl][i_h],POTscaling_tmp,yrange.at(i_h));
        }
      }
      c1->Print(std::string(histnames[i_h]+std::string("_plane")+std::to_string(i_pl)+".png").c_str());
    } // end loop over plot variables
// AF comment this out for now
//    TCanvas *c1 = new TCanvas();
//    double POTscaling_tmp = POTscaling; // Reset POT scaling for the next plot
//
//    if (templatefit){
//      TemplateFit(mc_hists_MIPdEdx[i_pl], onb_hists_MIPdEdx[i_pl], offb_hists_MIPdEdx[i_pl], offbeamscaling, POTscaling_tmp,-999);
//      POTscaling_tmp = 1.; // We have already applied POT scaling to MC in the template fit function. Set it to 1 so it doesn't get applied again in DrawMC.
//    }
//
//    if (onminusoffbeam){
//      DrawMC(mc_hists_MIPdEdx[i_pl],POTscaling_tmp,-999);
//      if (f_onbeam && f_offbeam){
//        OverlayOnMinusOffData(c1,onb_hists_MIPdEdx[i_pl],offb_hists_MIPdEdx[i_pl],offbeamscaling,POTscaling_tmp);
//        TString e_str("h_err");
//        TString o_str("h_ondat_MIPdEdx_plane"+std::to_string(i_pl)+"_all");
//        OverlayChi2(c1, e_str, o_str);
//      }
//    }
//    else{
//      if (f_onbeam && f_offbeam){
//        DrawMCPlusOffbeam(mc_hists_MIPdEdx[i_pl], offb_hists_MIPdEdx[i_pl], POTscaling_tmp, offbeamscaling,-999);
//        OverlayOnBeamData(c1, onb_hists_MIPdEdx[i_pl]);
//        TString e_str("h_err");
//        TString o_str("h_ondat_MIPdEdx_plane"+std::to_string(i_pl)+"_all");
//        OverlayChi2(c1, e_str, o_str);
//      }
//      else{
//        DrawMC(mc_hists_MIPdEdx[i_pl],POTscaling_tmp,-999);
//      }
//    }
//    c1->Print(std::string(std::string("h_MIPregiondEdx_plane")+std::to_string(i_pl)+".png").c_str());
  } // end loop over planes
  
  // Make the dE/dx plots in RR slices
  TCanvas *c1;
  for (size_t i_pl=0; i_pl<3; i_pl++){
    for (size_t slice(0); slice<10; slice++){
      delete c1;
      c1 = new TCanvas();
      DrawMCPlusOffbeam(mc_hists_dEdxSlices[i_pl][slice], offb_hists_dEdxSlices[i_pl][slice], POTscaling, offbeamscaling,-999);
      OverlayOnBeamData(c1, onb_hists_dEdxSlices[i_pl][slice]);
      c1->Print(std::string(std::string("h_dEdx_plane_")+std::to_string(i_pl)+"_slice_"+std::to_string(slice)+".png").c_str());
      c1->Print(std::string(std::string("h_dEdx_plane_")+std::to_string(i_pl)+"_slice_"+std::to_string(slice)+".eps").c_str());
    }
  }

}

#include "/uboone/app/users/afurmans/CCInclusive_PID_Jan/srcs/ubobj/ubobj/UBXSec/UBXSecEvent.h"
#include "plotFromTreeHeader.h"
#include "/uboone/app/users/piphamil/CCInclusive_PID_March11/srcs/ubana/ubana/ParticleID/Algorithms/Theory_dEdx_resrange.h"
#include "/uboone/app/users/piphamil/CCInclusive_PID_March11/srcs/ubana/ubana/ParticleID/Algorithms/Theory_dEdx_resrange.cxx"

//dEdx curves
particleid::Theory_dEdx_resrange particle_knower;

TGraph *muon_graph = particle_knower.g_ThdEdxRR_Muon;
TGraph *proton_graph = particle_knower.g_ThdEdxRR_Proton;

double MuonFunction(double *x, double *) {return muon_graph->Eval(x[0]);}
double ProtonFunction(double *x, double *) {return proton_graph->Eval(x[0]);}

TF1 *muon_curve = new TF1("muoncurve", MuonFunction, 0, 100, 0);
TF1 *proton_curve = new TF1("protoncurve", ProtonFunction, 0, 100, 0);

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
      //vars->track_dEdx_mean_atstart->at(i),
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
                    // {50,0,10}, // track_dEdx_mean_atstart
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
                    // ";Truncated Mean dE/dx at start of track;",
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
                  // "dEdx_truncmean_atstart",
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
                  // -999, // track_dEdx_mean_atstart
                  -999 // track_shift_bestL
                };

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void plotDataMCfromTree_Vandalised(std::string treename, std::string mcfile, double POTscaling=0., std::string onbeamdatafile="", std::string offbeamdatafile="", double offbeamscaling=0., bool onminusoffbeam=true, bool templatefit=false, int containment_switch=0, int muon_switch=0){

  //void plotDataMCfromTree(){ // in case you want to hard-code things, you can uncomment this block
  
 // std::string treename="pidvalidcaliSCE/pidTree";
////  std::string mcfile="/uboone/data/users/sfehlber/Mar/Mar11/mc/mc.root"; // MC pure
 // std::string mcfile="/uboone/data/users/renlu23/MCC9/April_v12/overlay_75k_v12.root"; // Overlay
 // std::string onbeamdatafile="/uboone/data/users/renlu23/MCC9/April_v12/bnb_v12.root"; // onbeam
 // std::string offbeamdatafile="/uboone/data/users/renlu23/MCC9/April_v12/extbnb_v12.root"; // onbeam
 // double offbeamscaling=0.587;
 // double POTscaling =0.547; //MC pure
 // //double POTscaling =6.202; Overlay
 // bool onminusoffbeam=false;
 // bool templatefit=false;
 // //PIP'S CHANGE: switches for plot splitting
 // int containment_switch = 0;
 // int muon_switch = 0;

  std::string dir0 = std::string("plots");
  std::string delimiter = "/";
  std::string dir1 = treename.substr(0, treename.find(delimiter));
  std::string dir2;
  if (muon_switch == 0) dir2 = std::string("all_tracks");
  if (muon_switch == 1) dir2 = std::string("muon_candidates");
  if (muon_switch == 2) dir2 = std::string("non_muon_candidates");
  std::string dir3;
  if (containment_switch == 0) dir3 = std::string("all");
  if (containment_switch == 1) dir3 = std::string("contained");
  if (containment_switch == 2) dir3 = std::string("uncontained");
  std::string output_dir = dir0 + std::string("/") + dir1 + std::string("/") + dir2 + std::string("/") + dir3 + std::string("/");


  //gStyle->SetTitleX(0.1f);
  //gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0.);
  gStyle->SetOptStat(0);

  gROOT->SetBatch(); // fill with 0 to turn off batch mode, if you want your screen full of plots

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
//  if (! f_bnbcos->GetListOfKeys()->Contains(treename.c_str())){
//    std::cout << "cannot find tree with name " << treename << " in file " << mcfile << std::endl;
//    return -1;
//  }
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get(treename.c_str());
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  TFile *f_onbeam=nullptr;
  TTree *t_onbeam=nullptr;
  treevars onbeam_vars;
  if (onbeamdatafile!=""){
    std::cout << "Making data-MC comparisons" << std::endl;
    f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
//    if (! f_onbeam->GetListOfKeys()->Contains(treename.c_str())){
//      std::cout << "cannot find tree with name " << treename << " in file " << onbeamdatafile << std::endl;
//      return -1;
//    }
    t_onbeam = (TTree*)f_onbeam->Get(treename.c_str());
    settreevars(t_onbeam,&onbeam_vars);
  }

  TFile *f_offbeam=nullptr;
  TTree *t_offbeam=nullptr;
  treevars offbeam_vars;
  if (offbeamdatafile!=""){
    f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
//    if (! f_offbeam->GetListOfKeys()->Contains(treename.c_str())){
//      std::cout << "cannot find tree with name " << treename << " in file " << offbeamdatafile << std::endl;
//      return -1;
//    }
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
  hist1D *mc_hists_alldEdx[nplanes];
  hist1D *mc_hists_MIPdEdx[nplanes];
  hist1D *mc_hists_dEdxSlices[nplanes][10];
  
  hist1D *mc_hists_trackShowerScore = new hist1D("h_mc_hists_trackShowerScore","",50,0,1);
  hist1D *mc_hists_trackShowerScore_less10cm = new hist1D("h_mc_hists_trackShowerScore_less10cm","",50,0,1);
  hist1D *mc_hists_llr_showers[nplanes];
  
  hist1D *mc_hists_dEdxslicedphi[nplanes][10];
  hist1D *mc_hists_tmeandEdx_slicedphi[nplanes][10];
  TH2F *mc_hists_dEdxvsresrange[nplanes];
  TH2F *mc_hists_dEdxvsthetaxz[nplanes];
  TH2F *mc_hists_dEdxvsthetayz[nplanes];
  TH2F *mc_hists_dEdxvsphi[nplanes];
  TH2F *mc_hists_dEdxvscostheta[nplanes];

  for (int i_pl=0; i_pl<nplanes; i_pl++){
    for (int i_h=0; i_h<nplots; i_h++){

      //PLANE LABELLING HACK - VERY SOPHISTICATED
      int i_plane = i_pl;
      if (i_h > 5 && i_h < 10)
      {
        if (i_pl == 0) i_plane = 2;
        if (i_pl == 2) i_plane = 0;
      }

      mc_hists[i_pl][i_h] = new hist1D(std::string("h_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_plane)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }
    mc_hists_llr_showers[i_pl] = new hist1D(std::string("mc_hists_llr_showers_")+std::to_string(i_pl),"",50,-10,10);
    for (int slice(0); slice < 10; slice++){

      mc_hists_dEdxslicedphi[i_pl][slice] = new hist1D("h_mc_hists_dEdxSlicedPhi_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);
      mc_hists_tmeandEdx_slicedphi[i_pl][slice] = new hist1D("h_mc_hists_truncmeandEdx_SlicedPhi_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);

      if (slice<2){
        mc_hists_dEdxSlices[i_pl][slice] = new hist1D("h_mc_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);}
      else{
        mc_hists_dEdxSlices[i_pl][slice] = new hist1D("h_mc_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,7);}

    }
    // Extra custom plots
    mc_hists_MIPdEdx[i_pl] = new hist1D(std::string("h_")+std::string("MIPregiondEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
    mc_hists_alldEdx[i_pl] = new hist1D(std::string("h_")+std::string("alldEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (all hits);"),100,0,10);
    mc_hists_dEdxvsresrange[i_pl] = new TH2F((std::string("h_")+std::string("dEdx_vs_resrange")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs residual range")).c_str(), 100, 0, 30, 100, 0, 10);
    mc_hists_dEdxvsthetaxz[i_pl] = new TH2F((std::string("h_")+std::string("dEdx_vs_thetaxz")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs theta_(xz)")).c_str(), 100, 0, 3.14, 100, 0, 10);
    mc_hists_dEdxvsthetayz[i_pl] = new TH2F((std::string("h_")+std::string("dEdx_vs_thetayz")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs theta_(yz)")).c_str(), 100, 0, 3.14, 100, 0, 10);
    mc_hists_dEdxvsphi[i_pl] = new TH2F((std::string("h_")+std::string("dEdx_vs_phi")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs phi")).c_str(), 100, -3.14, 3.14, 100, 0, 25);
    mc_hists_dEdxvscostheta[i_pl] = new TH2F((std::string("h_")+std::string("dEdx_vs_costheta")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs cos(theta)")).c_str(), 100, -1, 1, 100, 0, 10);

  }


  // Loop through MC tree and fill plots
  for (int i = 0; i < t_bnbcos->GetEntries(); i++){
    t_bnbcos->GetEntry(i);
    CalcPIDvars(&mc_vars, true);
    std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&mc_vars);

    //PIP'S CHANGE: switch splitting
    //--------------------------------
    if (containment_switch == 1 && !IsContained(mc_vars)) continue;
    if (containment_switch == 2 && IsContained(mc_vars)) continue;
    if (muon_switch == 1 && !mc_vars.track_ismuoncandidate) continue;
    if (muon_switch == 2 && mc_vars.track_ismuoncandidate) continue;
    //--------------------------------

    // if (mc_vars.track_theta_x > 75 && mc_vars.track_theta_x < 90){
      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        for (size_t i_h = 0; i_h < nplots; i_h++){
          FillHist(mc_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),mc_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<mc_vars.track_resrange_perhit_u->size(); i_hit++){
        FillHist(mc_hists_alldEdx[0],mc_vars.track_dEdx_perhit_u->at(i_hit),mc_vars.true_PDG);

        if (mc_vars.track_resrange_perhit_u->at(i_hit) > 0.3) mc_hists_dEdxvsresrange[0]->Fill(mc_vars.track_resrange_perhit_u->at(i_hit), mc_vars.track_dEdx_perhit_u->at(i_hit));
        mc_hists_dEdxvsthetaxz[0] ->Fill(ThetaXZ(mc_vars)                          , mc_vars.track_dEdx_perhit_u->at(i_hit));
        mc_hists_dEdxvsthetayz[0] ->Fill(ThetaYZ(mc_vars)                          , mc_vars.track_dEdx_perhit_u->at(i_hit));
        mc_hists_dEdxvsphi[0]     ->Fill(mc_vars.track_phi                         , mc_vars.track_dEdx_perhit_u->at(i_hit));
        mc_hists_dEdxvscostheta[0]->Fill(TMath::Cos(mc_vars.track_theta)           , mc_vars.track_dEdx_perhit_u->at(i_hit));

        if (mc_vars.track_resrange_perhit_u->at(i_hit)>100 && mc_vars.track_resrange_perhit_u->at(i_hit)<150){
          FillHist(mc_hists_MIPdEdx[0],mc_vars.track_dEdx_perhit_u->at(i_hit),mc_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<mc_vars.track_resrange_perhit_v->size(); i_hit++){
        FillHist(mc_hists_alldEdx[1],mc_vars.track_dEdx_perhit_v->at(i_hit),mc_vars.true_PDG);
        
        if (mc_vars.track_resrange_perhit_v->at(i_hit) > 0.3) mc_hists_dEdxvsresrange[1]->Fill(mc_vars.track_resrange_perhit_v->at(i_hit), mc_vars.track_dEdx_perhit_v->at(i_hit));
        mc_hists_dEdxvsthetaxz[1] ->Fill(ThetaXZ(mc_vars)                          , mc_vars.track_dEdx_perhit_v->at(i_hit));
        mc_hists_dEdxvsthetayz[1] ->Fill(ThetaYZ(mc_vars)                          , mc_vars.track_dEdx_perhit_v->at(i_hit));
        mc_hists_dEdxvsphi[1]     ->Fill(mc_vars.track_phi                         , mc_vars.track_dEdx_perhit_v->at(i_hit));
        mc_hists_dEdxvscostheta[1]->Fill(TMath::Cos(mc_vars.track_theta)           , mc_vars.track_dEdx_perhit_v->at(i_hit));

        if (mc_vars.track_resrange_perhit_v->at(i_hit)>100 && mc_vars.track_resrange_perhit_v->at(i_hit)<150){
          FillHist(mc_hists_MIPdEdx[1],mc_vars.track_dEdx_perhit_v->at(i_hit),mc_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<mc_vars.track_resrange_perhit_y->size(); i_hit++){
        FillHist(mc_hists_alldEdx[2],mc_vars.track_dEdx_perhit_y->at(i_hit),mc_vars.true_PDG);

        if (mc_vars.track_resrange_perhit_y->at(i_hit) > 0.3) mc_hists_dEdxvsresrange[2]->Fill(mc_vars.track_resrange_perhit_y->at(i_hit), mc_vars.track_dEdx_perhit_y->at(i_hit));
        mc_hists_dEdxvsthetaxz[2] ->Fill(ThetaXZ(mc_vars)                          , mc_vars.track_dEdx_perhit_y->at(i_hit));
        mc_hists_dEdxvsthetayz[2] ->Fill(ThetaYZ(mc_vars)                          , mc_vars.track_dEdx_perhit_y->at(i_hit));
        mc_hists_dEdxvsphi[2]     ->Fill(mc_vars.track_phi                         , mc_vars.track_dEdx_perhit_y->at(i_hit));
        mc_hists_dEdxvscostheta[2]->Fill(TMath::Cos(mc_vars.track_theta)           , mc_vars.track_dEdx_perhit_y->at(i_hit));

        if (mc_vars.track_resrange_perhit_y->at(i_hit)>100 && mc_vars.track_resrange_perhit_y->at(i_hit)<150){
          FillHist(mc_hists_MIPdEdx[2],mc_vars.track_dEdx_perhit_y->at(i_hit),mc_vars.true_PDG);
        }
      }

    double rr;
    double dedx;
    double phi;
    int slice;
    int phislice;
    int nhits;
    //if (mc_vars.track_ismuoncandidate==1){continue;}
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      if (i_pl==0){
        nhits=mc_vars.track_resrange_perhit_u->size();}
      if (i_pl==1){
        nhits=mc_vars.track_resrange_perhit_v->size();}
      if (i_pl==2){
        nhits=mc_vars.track_resrange_perhit_y->size();}

      std::vector<float> dEdx_vec;

      for (size_t i_hit=0; i_hit<nhits; i_hit++){
        if (i_pl==0){
          rr = mc_vars.track_resrange_perhit_u->at(i_hit);
          dedx = mc_vars.track_dEdx_perhit_u->at(i_hit);
          dEdx_vec.push_back((float) dedx);
          phi = mc_vars.track_phi; 
          slice=rr/5;
          phislice=5+phi/0.628;
        }
        if (i_pl==1){
          rr = mc_vars.track_resrange_perhit_v->at(i_hit);
          dedx = mc_vars.track_dEdx_perhit_v->at(i_hit);
          dEdx_vec.push_back((float) dedx);
          phi = mc_vars.track_phi;
          slice=rr/5;
          phislice=5+phi/0.628;
        }
        if (i_pl==2){
          rr = mc_vars.track_resrange_perhit_y->at(i_hit);
          dedx = mc_vars.track_dEdx_perhit_y->at(i_hit);
          dEdx_vec.push_back((float) dedx);
          phi = mc_vars.track_phi;
          slice=rr/5;
          phislice=5+phi/0.628;
        }
        if (slice>10){slice=10;}
        if (slice<0){continue;}
        FillHist(mc_hists_dEdxSlices[i_pl][slice],dedx,mc_vars.true_PDG);

        if (phislice>10){phislice=10;}
        if (phislice<0){continue;}
        FillHist(mc_hists_dEdxslicedphi[i_pl][phislice],dedx,mc_vars.true_PDG);
      }
      if (dEdx_vec.size() > 0)
      {
        double truncmean = TruncMeandEdx(dEdx_vec, 0);
        FillHist(mc_hists_tmeandEdx_slicedphi[i_pl][phislice], truncmean, mc_vars.true_PDG);
      }
      if (mc_vars.track_shower_score<0.5 and i_pl < 3){
        FillHist(mc_hists_llr_showers[i_pl], mc_vars.track_lnlikelihood_mipoverp->at(i_pl), mc_vars.true_PDG);
      }
    }
    if (mc_vars.track_lnlikelihood_mipoverp->at(2) < -1){
      FillHist(mc_hists_trackShowerScore, mc_vars.track_shower_score, mc_vars.true_PDG);
      if (mc_vars.track_length < 10){
        FillHist(mc_hists_trackShowerScore_less10cm, mc_vars.track_shower_score, mc_vars.true_PDG);
      }
    }
    //} // end theta_x cut
  } // end loop over entries in tree

  std::cout << "Done MC, now on-beam" << std::endl;
  // ----------------- On-beam data
  hist1D *onb_hists[nplanes][nplots];
  hist1D *onb_hists_alldEdx[nplanes];
  hist1D *onb_hists_MIPdEdx[nplanes];
  hist1D *onb_hists_dEdxSlices[nplanes][10];
  hist1D *onb_hists_trackShowerScore = new hist1D("h_onb_hists_trackShowerScore","",50,0,1);
  hist1D *onb_hists_trackShowerScore_less10cm = new hist1D("h_onb_hists_trackShowerScore_less10cm","",50,0,1);
  hist1D *onb_hists_llr_showers[nplanes];
  hist1D *onb_hists_dEdxslicedphi[nplanes][10];
  hist1D *onb_hists_tmeandEdx_slicedphi[nplanes][10];
  TH2F *onb_hists_dEdxvsresrange[nplanes];
  TH2F *onb_hists_dEdxvsthetaxz[nplanes];
  TH2F *onb_hists_dEdxvsthetayz[nplanes];
  TH2F *onb_hists_dEdxvsphi[nplanes];
  TH2F *onb_hists_dEdxvscostheta[nplanes];

  if (t_onbeam){
    // Make histograms to fill
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){

        //PLANE LABELLING HACK - VERY SOPHISTICATED
        int i_plane = i_pl;
        if (i_h > 5 && i_h < 10)
        {
          if (i_pl == 0) i_plane = 2;
          if (i_pl == 2) i_plane = 0;
        }

        onb_hists[i_pl][i_h] = new hist1D(std::string("h_ondat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_plane)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      }
      onb_hists_llr_showers[i_pl] = new hist1D(std::string("onb_hists_llr_showers_")+std::to_string(i_pl),"",50,-10,10);
      for (int slice(0); slice <10; slice++){

        onb_hists_dEdxslicedphi[i_pl][slice] = new hist1D("h_onb_hists_dEdxSlicedPhi_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);
        onb_hists_tmeandEdx_slicedphi[i_pl][slice] = new hist1D("h_onb_hists_truncmean_dEdxSlicedPhi_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);

        if (slice<2){
          onb_hists_dEdxSlices[i_pl][slice] = new hist1D("h_onb_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);}
        else{
          onb_hists_dEdxSlices[i_pl][slice] = new hist1D("h_onb_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,7);}
      }
      // Extra custom plots
      onb_hists_MIPdEdx[i_pl] = new hist1D(std::string("h_ondat_")+std::string("MIPregiondEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
      onb_hists_alldEdx[i_pl] = new hist1D(std::string("h_ondat_")+std::string("alldEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
      onb_hists_dEdxvsresrange[i_pl] = new TH2F((std::string("h_ondat_")+std::string("dEdx_vs_resrange")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs residual range")).c_str(), 100, 0, 30, 100, 0, 10);
      onb_hists_dEdxvsthetaxz[i_pl] = new TH2F((std::string("h_ondat__")+std::string("dEdx_vs_thetaxz")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs theta_(xz)")).c_str(), 100, 0, 3.14, 100, 0, 10);
      onb_hists_dEdxvsthetayz[i_pl] = new TH2F((std::string("h_ondat_")+std::string("dEdx_vs_thetayz")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs theta_(yz)")).c_str(), 100, 0, 3.14, 100, 0, 10);
      onb_hists_dEdxvsphi[i_pl] = new TH2F((std::string("h_ondat_")+std::string("dEdx_vs_phi")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs phi")).c_str(), 100, -3.14, 3.14, 100, 0, 25);
      onb_hists_dEdxvscostheta[i_pl] = new TH2F((std::string("h_ondat_")+std::string("dEdx_vs_costheta")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs cos(theta)")).c_str(), 100, -1, 1, 100, 0, 10);

    }

    // Loop through on-beam data tree and fill plots
    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);
      CalcPIDvars(&onbeam_vars, false);
      std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&onbeam_vars);

      //PIP'S CHANGE: switch splitting
      //--------------------------------
      if (containment_switch == 1 && !IsContained(onbeam_vars)) continue;
      if (containment_switch == 2 && IsContained(onbeam_vars)) continue;
      if (muon_switch == 1 && !onbeam_vars.track_ismuoncandidate) continue;
      if (muon_switch == 2 && onbeam_vars.track_ismuoncandidate) continue;
      //--------------------------------

//      if (onbeam_vars.track_theta_x > 75 && onbeam_vars.track_theta_x < 90){
        for (size_t i_pl=0; i_pl < nplanes; i_pl++){
          for (size_t i_h = 0; i_h < nplots; i_h++){
            FillHist(onb_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),0); // 0 because there is no "true PDG" for data
          }
        }
//      }

      for (size_t i_hit=0; i_hit<onbeam_vars.track_resrange_perhit_u->size(); i_hit++){
        FillHist(onb_hists_alldEdx[0],onbeam_vars.track_dEdx_perhit_u->at(i_hit),onbeam_vars.true_PDG);

        if(onbeam_vars.track_resrange_perhit_u->at(i_hit) > 0.3) onb_hists_dEdxvsresrange[0]->Fill(onbeam_vars.track_resrange_perhit_u->at(i_hit), onbeam_vars.track_dEdx_perhit_u->at(i_hit));
        onb_hists_dEdxvsthetaxz[0] ->Fill(ThetaXZ(onbeam_vars)                          , onbeam_vars.track_dEdx_perhit_u->at(i_hit));
        onb_hists_dEdxvsthetayz[0] ->Fill(ThetaYZ(onbeam_vars)                          , onbeam_vars.track_dEdx_perhit_u->at(i_hit));
        onb_hists_dEdxvsphi[0]     ->Fill(onbeam_vars.track_phi                         , onbeam_vars.track_dEdx_perhit_u->at(i_hit));
        onb_hists_dEdxvscostheta[0]->Fill(TMath::Cos(onbeam_vars.track_theta)           , onbeam_vars.track_dEdx_perhit_u->at(i_hit));

        if (onbeam_vars.track_resrange_perhit_u->at(i_hit)>100 && onbeam_vars.track_resrange_perhit_u->at(i_hit)<150){
          FillHist(onb_hists_MIPdEdx[0],onbeam_vars.track_dEdx_perhit_u->at(i_hit),0);
        }
      }

      for (size_t i_hit=0; i_hit<onbeam_vars.track_resrange_perhit_v->size(); i_hit++){
        FillHist(onb_hists_alldEdx[1],onbeam_vars.track_dEdx_perhit_v->at(i_hit),onbeam_vars.true_PDG);

        if(onbeam_vars.track_resrange_perhit_v->at(i_hit) > 0.3) onb_hists_dEdxvsresrange[1]->Fill(onbeam_vars.track_resrange_perhit_v->at(i_hit), onbeam_vars.track_dEdx_perhit_v->at(i_hit));
        onb_hists_dEdxvsthetaxz[1] ->Fill(ThetaXZ(onbeam_vars)                          , onbeam_vars.track_dEdx_perhit_v->at(i_hit));
        onb_hists_dEdxvsthetayz[1] ->Fill(ThetaYZ(onbeam_vars)                          , onbeam_vars.track_dEdx_perhit_v->at(i_hit));
        onb_hists_dEdxvsphi[1]     ->Fill(onbeam_vars.track_phi                         , onbeam_vars.track_dEdx_perhit_v->at(i_hit));
        onb_hists_dEdxvscostheta[1]->Fill(TMath::Cos(onbeam_vars.track_theta)           , onbeam_vars.track_dEdx_perhit_v->at(i_hit));

        if (onbeam_vars.track_resrange_perhit_v->at(i_hit)>100 && onbeam_vars.track_resrange_perhit_v->at(i_hit)<150){
          FillHist(onb_hists_MIPdEdx[1],onbeam_vars.track_dEdx_perhit_v->at(i_hit),0);
        }
      }

      for (size_t i_hit=0; i_hit<onbeam_vars.track_resrange_perhit_y->size(); i_hit++){
        FillHist(onb_hists_alldEdx[2],onbeam_vars.track_dEdx_perhit_y->at(i_hit),onbeam_vars.true_PDG);

        if (onbeam_vars.track_resrange_perhit_y->at(i_hit) > 0.3) onb_hists_dEdxvsresrange[2]->Fill(onbeam_vars.track_resrange_perhit_y->at(i_hit), onbeam_vars.track_dEdx_perhit_y->at(i_hit));
        onb_hists_dEdxvsthetaxz[2] ->Fill(ThetaXZ(onbeam_vars)                          , onbeam_vars.track_dEdx_perhit_y->at(i_hit));
        onb_hists_dEdxvsthetayz[2] ->Fill(ThetaYZ(onbeam_vars)                          , onbeam_vars.track_dEdx_perhit_y->at(i_hit));
        onb_hists_dEdxvsphi[2]     ->Fill(onbeam_vars.track_phi                         , onbeam_vars.track_dEdx_perhit_y->at(i_hit));
        onb_hists_dEdxvscostheta[2]->Fill(TMath::Cos(onbeam_vars.track_theta)           , onbeam_vars.track_dEdx_perhit_y->at(i_hit));

        if (onbeam_vars.track_resrange_perhit_y->at(i_hit)>100 && onbeam_vars.track_resrange_perhit_y->at(i_hit)<150){
          FillHist(onb_hists_MIPdEdx[2],onbeam_vars.track_dEdx_perhit_y->at(i_hit),0);
        }
      }
      double rr;
      double dedx;
      double phi;
      int slice;
      int phislice;
      int nhits;
      //if (onbeam_vars.track_ismuoncandidate==1){continue;}
      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      
        std::vector<float> dEdx_vec;
       
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
            dEdx_vec.push_back((float) dedx);
            phi = onbeam_vars.track_phi;
            slice=rr/5;
            phislice=5+phi/0.628;
            if (slice>10){slice=10;}
          }
          if (i_pl==1){
            rr = onbeam_vars.track_resrange_perhit_v->at(i_hit);
            dedx = onbeam_vars.track_dEdx_perhit_v->at(i_hit);
            dEdx_vec.push_back((float) dedx);
            phi = onbeam_vars.track_phi;
            slice=rr/5;
            phislice=5+phi/0.628;
            if (slice>10){slice=10;}
            if (phislice>10){phislice=10;}
          }
          if (i_pl==2){
            rr = onbeam_vars.track_resrange_perhit_y->at(i_hit);
            dedx = onbeam_vars.track_dEdx_perhit_y->at(i_hit);
            dEdx_vec.push_back((float) dedx);
            phi = onbeam_vars.track_phi;
            slice=rr/5;
            phislice=5+phi/0.628;
            if (slice>10){slice=10;}
            if (phislice>10){phislice=10;}
          }
          if (slice<0){continue;}
          FillHist(onb_hists_dEdxSlices[i_pl][slice],dedx,0);
          if (phislice<0){continue;}
          FillHist(onb_hists_dEdxslicedphi[i_pl][phislice],dedx,0);
        }
        if (dEdx_vec.size() > 0) FillHist(onb_hists_tmeandEdx_slicedphi[i_pl][phislice],TruncMeandEdx(dEdx_vec,0),0);
        if (onbeam_vars.track_shower_score<0.5 and i_pl < 3){
          FillHist(onb_hists_llr_showers[i_pl], onbeam_vars.track_lnlikelihood_mipoverp->at(i_pl), 0);
        }
      }
      if (onbeam_vars.track_lnlikelihood_mipoverp->at(2) < -1){
        FillHist(onb_hists_trackShowerScore, onbeam_vars.track_shower_score, 0);
        if (onbeam_vars.track_length<10){
          FillHist(onb_hists_trackShowerScore_less10cm, onbeam_vars.track_shower_score, 0);
        }
      }
    }
  }
    // ----------------- Off-beam data
  hist1D *offb_hists[nplanes][nplots];
  hist1D *offb_hists_alldEdx[nplanes];
  hist1D *offb_hists_MIPdEdx[nplanes];
  hist1D *offb_hists_dEdxSlices[nplanes][10];
  hist1D *offb_hists_trackShowerScore = new hist1D("h_offb_hists_trackShowerScore","",50,0,1);
  hist1D *offb_hists_trackShowerScore_less10cm = new hist1D("h_offb_hists_trackShowerScore_less10cm","",50,0,1);
  hist1D *offb_hists_llr_showers[nplanes];
  hist1D *offb_hists_dEdxslicedphi[nplanes][10];
  hist1D *offb_hists_tmeandEdx_slicedphi[nplanes][10];
  TH2F *offb_hists_dEdxvsresrange[nplanes];
  TH2F *offb_hists_dEdxvsthetaxz[nplanes];
  TH2F *offb_hists_dEdxvsthetayz[nplanes];
  TH2F *offb_hists_dEdxvsphi[nplanes];
  TH2F *offb_hists_dEdxvscostheta[nplanes];


  std::cout << "Off-beam time!" << std::endl;
  if (t_offbeam){
    // Make histograms to fill
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){

        //PLANE LABELLING HACK - VERY SOPHISTICATED
        int i_plane = i_pl;
        if (i_h > 5 && i_h < 10)
        {
          if (i_pl == 0) i_plane = 2;
          if (i_pl == 2) i_plane = 0;
        }

        offb_hists[i_pl][i_h] = new hist1D(std::string("h_offdat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_plane)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));

      }
      offb_hists_llr_showers[i_pl] = new hist1D(std::string("offb_hists_llr_showers_")+std::to_string(i_pl),"",50,-10,10);

      for (int slice(0); slice <10; slice++){

        offb_hists_dEdxslicedphi[i_pl][slice] = new hist1D("h_offb_hists_dEdxSlicedPhi_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);
        offb_hists_tmeandEdx_slicedphi[i_pl][slice] = new hist1D("h_offb_hists_truncmeandEdx_SlicedPhi_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);

        if (slice<2){
          offb_hists_dEdxSlices[i_pl][slice] = new hist1D("h_offb_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,10);}
        else{
          offb_hists_dEdxSlices[i_pl][slice] = new hist1D("h_offb_hists_dEdxSlices_"+std::to_string(i_pl)+"_"+std::to_string(slice),"",50,0,7);}
      }
      // Extra custom plots
      offb_hists_MIPdEdx[i_pl] = new hist1D(std::string("h_offdat_")+std::string("MIPregiondEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
      offb_hists_alldEdx[i_pl] = new hist1D(std::string("h_offdat_")+std::string("alldEdx")+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx per hit (res. range 100-150 cm);"),100,0,10);
      offb_hists_dEdxvsresrange[i_pl] = new TH2F((std::string("h_offdat_")+std::string("dEdx_vs_resrange")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs residual range")).c_str(), 100, 0, 30, 100, 0, 10);
      offb_hists_dEdxvsthetaxz[i_pl] = new TH2F((std::string("h_offdat__")+std::string("dEdx_vs_thetaxz")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs theta_(xz)")).c_str(), 100, 0, 3.14, 100, 0, 10);
      offb_hists_dEdxvsthetayz[i_pl] = new TH2F((std::string("h_offdat_")+std::string("dEdx_vs_thetayz")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs theta_(yz)")).c_str(), 100, 0, 3.14, 100, 0, 10);
      offb_hists_dEdxvsphi[i_pl] = new TH2F((std::string("h_offdat_")+std::string("dEdx_vs_phi")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs phi")).c_str(), 100, -3.14, 3.14, 100, 0, 25);
      offb_hists_dEdxvscostheta[i_pl] = new TH2F((std::string("h_offdat_")+std::string("dEdx_vs_costheta")+std::string("_plane")+std::to_string(i_pl)).c_str(), (std::string("Plane ")+std::to_string(i_pl)+std::string(";dE/dx vs cos(theta)")).c_str(), 100, -1, 1, 100, 0, 10);


    }


    // Loop through tree and fill plots
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
      t_offbeam->GetEntry(i);
      CalcPIDvars(&offbeam_vars, false);
      std::vector<std::vector<double>> PIDvarstoplot = GetPIDvarstoplot(&offbeam_vars);

      //PIP'S CHANGE: switch splitting
      //--------------------------------
      if (containment_switch == 1 && !IsContained(offbeam_vars)) continue;
      if (containment_switch == 2 && IsContained(offbeam_vars)) continue;
      if (muon_switch == 1 && !offbeam_vars.track_ismuoncandidate) continue;
      if (muon_switch == 2 && offbeam_vars.track_ismuoncandidate) continue;
      //--------------------------------

      //if (offbeam_vars.track_theta_x >75 && offbeam_vars.track_theta_x < 90){
        for (size_t i_pl=0; i_pl < nplanes; i_pl++){
          for (size_t i_h = 0; i_h < nplots; i_h++){
            FillHist(offb_hists[i_pl][i_h],PIDvarstoplot.at(i_pl).at(i_h),0.); // 0 because there is no "true PDG" for data
          }
        }
      //}

      for (size_t i_hit=0; i_hit<offbeam_vars.track_resrange_perhit_u->size(); i_hit++){
        FillHist(offb_hists_alldEdx[0],offbeam_vars.track_dEdx_perhit_u->at(i_hit),offbeam_vars.true_PDG);

        if (offbeam_vars.track_resrange_perhit_u->at(i_hit) > 0.3) offb_hists_dEdxvsresrange[0]->Fill(offbeam_vars.track_resrange_perhit_u->at(i_hit), offbeam_vars.track_dEdx_perhit_u->at(i_hit));
        offb_hists_dEdxvsthetaxz[0] ->Fill(ThetaXZ(offbeam_vars)                          , offbeam_vars.track_dEdx_perhit_u->at(i_hit));
        offb_hists_dEdxvsthetayz[0] ->Fill(ThetaYZ(offbeam_vars)                          , offbeam_vars.track_dEdx_perhit_u->at(i_hit));
        offb_hists_dEdxvsphi[0]     ->Fill(offbeam_vars.track_phi                         , offbeam_vars.track_dEdx_perhit_u->at(i_hit));
        offb_hists_dEdxvscostheta[0]->Fill(TMath::Cos(offbeam_vars.track_theta)           , offbeam_vars.track_dEdx_perhit_u->at(i_hit));

        if (offbeam_vars.track_resrange_perhit_u->at(i_hit)>100 && offbeam_vars.track_resrange_perhit_u->at(i_hit)<150){
          FillHist(offb_hists_MIPdEdx[0],offbeam_vars.track_dEdx_perhit_u->at(i_hit),offbeam_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<offbeam_vars.track_resrange_perhit_v->size(); i_hit++){
        FillHist(offb_hists_alldEdx[1],offbeam_vars.track_dEdx_perhit_v->at(i_hit),offbeam_vars.true_PDG);

        if (offbeam_vars.track_resrange_perhit_v->at(i_hit) > 0.3) offb_hists_dEdxvsresrange[1]->Fill(offbeam_vars.track_resrange_perhit_v->at(i_hit), offbeam_vars.track_dEdx_perhit_v->at(i_hit));
        offb_hists_dEdxvsthetaxz[1] ->Fill(ThetaXZ(offbeam_vars)                          , offbeam_vars.track_dEdx_perhit_v->at(i_hit));
        offb_hists_dEdxvsthetayz[1] ->Fill(ThetaYZ(offbeam_vars)                          , offbeam_vars.track_dEdx_perhit_v->at(i_hit));
        offb_hists_dEdxvsphi[1]     ->Fill(offbeam_vars.track_phi                         , offbeam_vars.track_dEdx_perhit_v->at(i_hit));
        offb_hists_dEdxvscostheta[1]->Fill(TMath::Cos(offbeam_vars.track_theta)           , offbeam_vars.track_dEdx_perhit_v->at(i_hit));

        if (offbeam_vars.track_resrange_perhit_v->at(i_hit)>100 && offbeam_vars.track_resrange_perhit_v->at(i_hit)<150){
          FillHist(offb_hists_MIPdEdx[1],offbeam_vars.track_dEdx_perhit_v->at(i_hit),offbeam_vars.true_PDG);
        }
      }

      for (size_t i_hit=0; i_hit<offbeam_vars.track_resrange_perhit_y->size(); i_hit++){
        FillHist(offb_hists_alldEdx[2],offbeam_vars.track_dEdx_perhit_y->at(i_hit),offbeam_vars.true_PDG);

        if (offbeam_vars.track_resrange_perhit_y->at(i_hit) > 0.3) offb_hists_dEdxvsresrange[2]->Fill(offbeam_vars.track_resrange_perhit_y->at(i_hit), offbeam_vars.track_dEdx_perhit_y->at(i_hit));
        offb_hists_dEdxvsthetaxz[2] ->Fill(ThetaXZ(offbeam_vars)                          , offbeam_vars.track_dEdx_perhit_y->at(i_hit));
        offb_hists_dEdxvsthetayz[2] ->Fill(ThetaYZ(offbeam_vars)                          , offbeam_vars.track_dEdx_perhit_y->at(i_hit));
        offb_hists_dEdxvsphi[2]     ->Fill(offbeam_vars.track_phi                         , offbeam_vars.track_dEdx_perhit_y->at(i_hit));
        offb_hists_dEdxvscostheta[2]->Fill(TMath::Cos(offbeam_vars.track_theta)           , offbeam_vars.track_dEdx_perhit_y->at(i_hit));


        if (offbeam_vars.track_resrange_perhit_y->at(i_hit)>100 && offbeam_vars.track_resrange_perhit_y->at(i_hit)<150){
          FillHist(offb_hists_MIPdEdx[2],offbeam_vars.track_dEdx_perhit_y->at(i_hit),offbeam_vars.true_PDG);
        }
      }
      double rr;
      double dedx;
      double phi;
      int slice;
      int phislice;
      int nhits;
      //if (offbeam_vars.track_ismuoncandidate==1){continue;}
      for (size_t i_pl=0; i_pl < nplanes; i_pl++){

        std::vector<float> dEdx_vec;

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
            dEdx_vec.push_back((float) dedx);
            phi = offbeam_vars.track_phi;
            slice=rr/5;
            phislice=5+phi/0.628;
          }
          if (i_pl==1){
            rr = offbeam_vars.track_resrange_perhit_v->at(i_hit);
            dedx = offbeam_vars.track_dEdx_perhit_v->at(i_hit);
            dEdx_vec.push_back((float) dedx);
            phi = offbeam_vars.track_phi;
            slice=rr/5;
            phislice=5+phi/0.628;
          }
          if (i_pl==2){
            rr = offbeam_vars.track_resrange_perhit_y->at(i_hit);
            dedx = offbeam_vars.track_dEdx_perhit_y->at(i_hit);
            dEdx_vec.push_back((float) dedx);
            phi = offbeam_vars.track_phi;
            slice=rr/5;
            phislice=5+phi/0.628;
          }
          if (slice>10){slice=10;}
          if (slice<0){continue;}
          FillHist(offb_hists_dEdxSlices[i_pl][slice],dedx,0);

          if (phislice>10){phislice=10;}
          if (phislice<0){continue;}
          FillHist(offb_hists_dEdxslicedphi[i_pl][phislice],dedx,0);

        }
        if (dEdx_vec.size() > 0) FillHist(offb_hists_tmeandEdx_slicedphi[i_pl][phislice], TruncMeandEdx(dEdx_vec,0),0);
        if (offbeam_vars.track_shower_score<0.5 and i_pl < 3){
          FillHist(offb_hists_llr_showers[i_pl], offbeam_vars.track_lnlikelihood_mipoverp->at(i_pl), 0);
        }
      }
      if (offbeam_vars.track_lnlikelihood_mipoverp->at(2) < -1){
        FillHist(offb_hists_trackShowerScore, offbeam_vars.track_shower_score, 0);
        if (offbeam_vars.track_length<10){
          FillHist(offb_hists_trackShowerScore_less10cm, offbeam_vars.track_shower_score, 0);
        }
      }
    } // end loop over entries in tree
  }

  // -------------------- Now make all the plots
  std::cout << "TIME TO MAKE SOME PLOTS" << std::endl;
  for (size_t i_pl=0; i_pl < nplanes; i_pl++){
    if (i_pl == 3) continue;
    
    for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas();

      //PLANE LABELLING HACK - VERY SOPHISTICATED
      int i_plane = i_pl;
      if (i_h > 5 && i_h < 10)
      {
        if (i_pl == 0) i_plane = 2;
        if (i_pl == 2) i_plane = 0;
      }

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
      c1->Print(std::string(output_dir+histnames[i_h]+std::string("_plane")+std::to_string(i_plane)+".png").c_str());
    } // end loop over plot variables

   TCanvas *c1 = new TCanvas();
   double POTscaling_tmp = POTscaling; // Reset POT scaling for the next plot

   if (templatefit){
     TemplateFit(mc_hists_MIPdEdx[i_pl], onb_hists_MIPdEdx[i_pl], offb_hists_MIPdEdx[i_pl], offbeamscaling, POTscaling_tmp,-999);
     POTscaling_tmp = 1.; // We have already applied POT scaling to MC in the template fit function. Set it to 1 so it doesn't get applied again in DrawMC.
   }


  if (onminusoffbeam){
    DrawMC(mc_hists_MIPdEdx[i_pl],POTscaling_tmp,-999);
    if (f_onbeam && f_offbeam){
      OverlayOnMinusOffData(c1,onb_hists_MIPdEdx[i_pl],offb_hists_MIPdEdx[i_pl],offbeamscaling,POTscaling_tmp);
      TString e_str("h_err");
      TString o_str("h_ondat_MIPregiondEdx_plane"+std::to_string(i_pl)+"_all");
      OverlayChi2(c1, e_str, o_str);
    }
  }
  else{
    if (f_onbeam && f_offbeam){
      DrawMCPlusOffbeam(mc_hists_MIPdEdx[i_pl], offb_hists_MIPdEdx[i_pl], POTscaling_tmp, offbeamscaling,-999);
      OverlayOnBeamData(c1, onb_hists_MIPdEdx[i_pl]);
      TString e_str("h_err");
      TString o_str("h_ondat_MIPregiondEdx_plane"+std::to_string(i_pl)+"_all");
      OverlayChi2(c1, e_str, o_str);
    }
    else{
      DrawMC(mc_hists_MIPdEdx[i_pl],POTscaling_tmp,-999);
    }
  }
  c1->Print(std::string(output_dir+std::string("h_MIPregiondEdx_plane")+std::to_string(i_pl)+".png").c_str());

  POTscaling_tmp = POTscaling; // Reset POT scaling for the next plot

  if (templatefit){
    TemplateFit(mc_hists_alldEdx[i_pl], onb_hists_alldEdx[i_pl], offb_hists_alldEdx[i_pl], offbeamscaling, POTscaling_tmp,-999);
    POTscaling_tmp = 1.; // We have already applied POT scaling to MC in the template fit function. Set it to 1 so it doesn't get applied again in DrawMC.
  }

  if (onminusoffbeam){
    DrawMC(mc_hists_alldEdx[i_pl],POTscaling_tmp,-999);
    if (f_onbeam && f_offbeam){
      OverlayOnMinusOffData(c1,onb_hists_alldEdx[i_pl],offb_hists_alldEdx[i_pl],offbeamscaling,POTscaling_tmp);
      TString e_str("h_err");
      TString o_str("h_ondat_alldEdx_plane"+std::to_string(i_pl)+"_all");
      OverlayChi2(c1, e_str, o_str);
    }
  }
  else{
    if (f_onbeam && f_offbeam){
      DrawMCPlusOffbeam(mc_hists_alldEdx[i_pl], offb_hists_alldEdx[i_pl], POTscaling_tmp, offbeamscaling,-999);
      OverlayOnBeamData(c1, onb_hists_alldEdx[i_pl]);
      TString e_str("h_err");
      TString o_str("h_ondat_alldEdx_plane"+std::to_string(i_pl)+"_all");
      OverlayChi2(c1, e_str, o_str);
    }
    else{
      DrawMC(mc_hists_alldEdx[i_pl],POTscaling_tmp,-999);
    }
  }
  c1->Print(std::string(output_dir+std::string("h_alldEdx_plane")+std::to_string(i_pl)+".png").c_str());

  // Make plots of dE/dx vs things
 
  mc_hists_dEdxvsresrange[i_pl]->Scale(POTscaling);
  offb_hists_dEdxvsresrange[i_pl]->Scale(offbeamscaling);
  mc_hists_dEdxvsresrange[i_pl]->Add(offb_hists_dEdxvsresrange[i_pl]);
  mc_hists_dEdxvsresrange[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  mc_hists_dEdxvsresrange[i_pl]->GetXaxis()->SetTitle("Residual range (cm)");
  mc_hists_dEdxvsresrange[i_pl]->Draw("colz");
  muon_curve->Draw("same");
  proton_curve->Draw("same");
  c1->Print(std::string(output_dir+std::string("h_mc_dEdxvsresrange_plane")+std::to_string(i_pl)+".png").c_str());

  mc_hists_dEdxvsthetaxz[i_pl]->Scale(POTscaling);
  offb_hists_dEdxvsthetaxz[i_pl]->Scale(offbeamscaling);
  mc_hists_dEdxvsthetaxz[i_pl]->Add(offb_hists_dEdxvsthetaxz[i_pl]);
  mc_hists_dEdxvsthetaxz[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  mc_hists_dEdxvsthetaxz[i_pl]->GetXaxis()->SetTitle("#theta_{xz}");
  mc_hists_dEdxvsthetaxz[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_mc_dEdxvsthetaxz_plane")+std::to_string(i_pl)+".png").c_str());

  mc_hists_dEdxvsthetayz[i_pl]->Scale(POTscaling);
  offb_hists_dEdxvsthetayz[i_pl]->Scale(offbeamscaling);
  mc_hists_dEdxvsthetayz[i_pl]->Add(offb_hists_dEdxvsthetayz[i_pl]);
  mc_hists_dEdxvsthetayz[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  mc_hists_dEdxvsthetayz[i_pl]->GetXaxis()->SetTitle("#Theta_{yz}");
  mc_hists_dEdxvsthetayz[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_mc_dEdxvsthetayz_plane")+std::to_string(i_pl)+".png").c_str());

  mc_hists_dEdxvsphi[i_pl]->Scale(POTscaling);
  offb_hists_dEdxvsphi[i_pl]->Scale(offbeamscaling);
  mc_hists_dEdxvsphi[i_pl]->Add(offb_hists_dEdxvsphi[i_pl]);
  mc_hists_dEdxvsphi[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  mc_hists_dEdxvsphi[i_pl]->GetXaxis()->SetTitle("#phi");
  mc_hists_dEdxvsphi[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_mc_dEdxvsphi_plane")+std::to_string(i_pl)+".png").c_str());

  mc_hists_dEdxvscostheta[i_pl]->Scale(POTscaling);
  offb_hists_dEdxvscostheta[i_pl]->Scale(offbeamscaling);
  mc_hists_dEdxvscostheta[i_pl]->Add(offb_hists_dEdxvscostheta[i_pl]);
  mc_hists_dEdxvscostheta[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  mc_hists_dEdxvscostheta[i_pl]->GetXaxis()->SetTitle("cos(#theta)");
  mc_hists_dEdxvscostheta[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_mc_dEdxvscostheta_plane")+std::to_string(i_pl)+".png").c_str());

  onb_hists_dEdxvsresrange[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  onb_hists_dEdxvsresrange[i_pl]->GetXaxis()->SetTitle("Residual range (cm)");
  onb_hists_dEdxvsresrange[i_pl]->Draw("colz");
  muon_curve->Draw("same");
  proton_curve->Draw("same");
  c1->Print(std::string(output_dir+std::string("h_data_dEdxvsresrange_plane")+std::to_string(i_pl)+".png").c_str());

  onb_hists_dEdxvsthetaxz[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  onb_hists_dEdxvsthetaxz[i_pl]->GetXaxis()->SetTitle("#theta_{xz}");
  onb_hists_dEdxvsthetaxz[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_data_dEdxvsthetaxz_plane")+std::to_string(i_pl)+".png").c_str());

  onb_hists_dEdxvsthetayz[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  onb_hists_dEdxvsthetayz[i_pl]->GetXaxis()->SetTitle("#Theta_{yz}");
  onb_hists_dEdxvsthetayz[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_data_dEdxvsthetayz_plane")+std::to_string(i_pl)+".png").c_str());

  onb_hists_dEdxvsphi[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  onb_hists_dEdxvsphi[i_pl]->GetXaxis()->SetTitle("#phi");
  onb_hists_dEdxvsphi[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_data_dEdxvsphi_plane")+std::to_string(i_pl)+".png").c_str());

  onb_hists_dEdxvscostheta[i_pl]->GetYaxis()->SetTitle("dE/dx per hit (all hits)");
  onb_hists_dEdxvscostheta[i_pl]->GetXaxis()->SetTitle("cos(#theta)");
  onb_hists_dEdxvscostheta[i_pl]->Draw("colz");
  c1->Print(std::string(output_dir+std::string("h_data_dEdxvscostheta_plane")+std::to_string(i_pl)+".png").c_str());


  // Make the dE/dx plots in RR slices
  for (size_t slice(0); slice<10; slice++){
    c1->cd();
    FormatYourPlotsAndy(mc_hists_dEdxSlices[i_pl][slice], slice);
    FormatYourPlotsAndy(offb_hists_dEdxSlices[i_pl][slice], slice);
    FormatYourPlotsAndy(onb_hists_dEdxSlices[i_pl][slice], slice);

    //double proton_exp = 17*pow((slice*5+2.5),(-0.42));
    double lineheight = onb_hists_dEdxSlices[i_pl][slice]->h_all->GetMaximum();
    double proton_exp = particle_knower.g_ThdEdxRR_Proton->Eval(2.5+slice*5, 0, "S");
    double muon_exp   = particle_knower.g_ThdEdxRR_Muon->Eval(2.5+slice*5, 0, "S");

    TLine *protonline = new TLine(proton_exp, 0, proton_exp, lineheight);
    protonline->SetLineWidth(3);
    protonline->SetLineColor(kRed);
    TLine *muonline = new TLine(muon_exp, 0, muon_exp, lineheight);
    muonline->SetLineWidth(3);
    muonline->SetLineColor(kViolet);

    DrawMCPlusOffbeam(mc_hists_dEdxSlices[i_pl][slice], offb_hists_dEdxSlices[i_pl][slice], POTscaling, offbeamscaling,-999);
    OverlayOnBeamData(c1, onb_hists_dEdxSlices[i_pl][slice]);
    protonline->Draw("same");
    muonline->Draw("same");
    c1->Print(std::string(output_dir+std::string("h_dEdx_plane_")+std::to_string(i_pl)+"_slice_"+std::to_string(slice)+".png").c_str());
    //c1->Print(std::string(output_dir+std::string("h_dEdx_plane_")+std::to_string(i_pl)+"_slice_"+std::to_string(slice)+".eps").c_str());
    
    //Make the dE/dx plots in phi slices

    FormatPhiSlices(mc_hists_dEdxslicedphi[i_pl][slice], slice);
    FormatPhiSlices(offb_hists_dEdxslicedphi[i_pl][slice], slice);
    FormatPhiSlices(onb_hists_dEdxslicedphi[i_pl][slice], slice);

    DrawMCPlusOffbeam(mc_hists_dEdxslicedphi[i_pl][slice], offb_hists_dEdxslicedphi[i_pl][slice], POTscaling, offbeamscaling,-999);
    OverlayOnBeamData(c1, onb_hists_dEdxslicedphi[i_pl][slice]);
    c1->Print(std::string(output_dir+std::string("h_dEdx_plane_")+std::to_string(i_pl)+"_phislice_"+std::to_string(slice)+".png").c_str());

    FormatPhiSlices(mc_hists_tmeandEdx_slicedphi[i_pl][slice], slice);
    FormatPhiSlices(offb_hists_tmeandEdx_slicedphi[i_pl][slice], slice);
    FormatPhiSlices(onb_hists_tmeandEdx_slicedphi[i_pl][slice], slice);

    mc_hists_tmeandEdx_slicedphi[i_pl][slice]->h_all->GetXaxis()->SetTitle("Truncated mean dE/dx");
    offb_hists_tmeandEdx_slicedphi[i_pl][slice]->h_all->GetXaxis()->SetTitle("Truncated mean dE/dx");
    onb_hists_tmeandEdx_slicedphi[i_pl][slice]->h_all->GetXaxis()->SetTitle("Truncated mean dE/dx");

    DrawMCPlusOffbeam(mc_hists_tmeandEdx_slicedphi[i_pl][slice], offb_hists_tmeandEdx_slicedphi[i_pl][slice], POTscaling, offbeamscaling,-999);
    OverlayOnBeamData(c1, onb_hists_tmeandEdx_slicedphi[i_pl][slice]);
    c1->Print(std::string(output_dir+std::string("h_truncmeandEdx_plane_")+std::to_string(i_pl)+"_phislice_"+std::to_string(slice)+".png").c_str());


  }

  DrawMCPlusOffbeam(mc_hists_llr_showers[i_pl], offb_hists_llr_showers[i_pl], POTscaling, offbeamscaling, -999);
  OverlayOnBeamData(c1, onb_hists_llr_showers[i_pl]);
  c1->Print(std::string(output_dir+std::string("h_llr_showers_")+std::to_string(i_pl)+".eps").c_str());
  

  } // end loop over planes

  // Make track-shower score plot for protons
  TCanvas *c1 = new TCanvas();
  DrawMCPlusOffbeam(mc_hists_trackShowerScore, offb_hists_trackShowerScore, POTscaling, offbeamscaling, -999);
  OverlayOnBeamData(c1, onb_hists_trackShowerScore);
  c1->Print(std::string(output_dir+std::string("h_trackShowerScore_protonLike.eps")).c_str());
  c1->Print(std::string(output_dir+std::string("h_trackShowerScore_protonLike.png")).c_str());

  // Make track-shower score plot for protons
  //TCanvas *c1 = new TCanvas();
  c1->Clear();
  DrawMCPlusOffbeam(mc_hists_trackShowerScore_less10cm, offb_hists_trackShowerScore_less10cm, POTscaling, offbeamscaling, -999);
  OverlayOnBeamData(c1, onb_hists_trackShowerScore_less10cm);
  c1->Print(std::string(output_dir+std::string("h_trackShowerScore_protonLike_less10cm.eps")).c_str());
  c1->Print(std::string(output_dir+std::string("h_trackShowerScore_protonLike_less10cm.png")).c_str());


}

void plotDataMCfromTree()
{
  muon_curve->SetLineWidth(3);
  muon_curve->SetLineColor(kViolet);
  
  proton_curve->SetLineWidth(3);
  proton_curve->SetLineColor(kBlue);

  //std::string mcfile="/uboone/data/users/renlu23/MCC9/April_v12/overlay_75k_v12.root"; // Overlay
  //std::string onbeamdatafile="/uboone/data/users/renlu23/MCC9/April_v12/bnb_v12.root"; // onbeam
  //std::string offbeamdatafile="/uboone/data/users/renlu23/MCC9/April_v12/extbnb_v12.root"; // offbeam
  std::string mcfile="/uboone/data/users/afurmans/MCC9validation/ubxsec_trackshower_output_data_overlay.root"; // Overlay
  std::string onbeamdatafile="/uboone/data/users/afurmans/MCC9validation/ubxsec_trackshower_output_data_bnb.root"; // onbeam
  std::string offbeamdatafile="/uboone/data/users/afurmans/MCC9validation/ubxsec_trackshower_output_data_ext.root"; // offbeam
  double offbeamscaling=0.4569;
  double POTscaling =0.277; //MC pure

  //plotDataMCfromTree_Vandalised(std::string("pidvalidcaliSCEtracksFromShower/pidTree"), mcfile, POTscaling, onbeamdatafile, offbeamdatafile, offbeamscaling, false, false, 1, 2);

  std::string treename;
  for (int i = 0; i < 4; i++)
  {
    if (i == 0) treename = std::string("pidvalidcalo/pidTree");
    if (i == 1) treename = std::string("pidvalidcaloSCE/pidTree");
    if (i == 2) treename = std::string("pidvalidcali/pidTree");
    if (i == 3) treename = std::string("pidvalidcaliSCE/pidTree");
    if (i == 3) treename = std::string("pidvalidcaliSCEtracksFromShower/pidTree");

    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        plotDataMCfromTree_Vandalised(treename, mcfile, POTscaling, onbeamdatafile, offbeamdatafile, offbeamscaling, false, false, j, k);
      }
    }
  }

  return;
}



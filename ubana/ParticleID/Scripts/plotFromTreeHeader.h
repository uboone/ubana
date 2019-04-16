// #include "../../../../larana/larana/TruncatedMean/Algorithm/TruncMean.cxx"

struct treevars{
  // These are the variables that are filled directly from the tree
  int true_PDG=-9999;
  std::vector<double> *track_likelihood_fwd_p=nullptr;
  std::vector<double> *track_likelihood_fwd_mu=nullptr;
  std::vector<double> *track_likelihood_fwd_pi=nullptr;
  std::vector<double> *track_likelihood_fwd_k=nullptr;
  std::vector<double> *track_likelihood_fwd_other=nullptr;
  std::vector<double> *track_likelihood_fwd_mip=nullptr;
  std::vector<double> *track_likelihood_bwd_p=nullptr;
  std::vector<double> *track_likelihood_bwd_mu=nullptr;
  std::vector<double> *track_likelihood_bwd_pi=nullptr;
  std::vector<double> *track_likelihood_bwd_k=nullptr;
  std::vector<double> *track_likelihood_bwd_other=nullptr;
  std::vector<double> *track_PIDA_mean=nullptr;
  std::vector<double> *track_PIDA_kde=nullptr;
  std::vector<double> *track_PIDA_median=nullptr;
  std::vector<double> *track_depE=nullptr;
  std::vector<double> *track_likelihood_shift_fwd_mu=nullptr;
  std::vector<double> *track_likelihood_shift_fwd_p=nullptr;
  std::vector<double> *track_likelihood_shift_fwd_pi=nullptr;
  std::vector<double> *track_likelihood_shift_fwd_k=nullptr;
  std::vector<double> *track_likelihood_shift_bwd_mu=nullptr;
  std::vector<double> *track_likelihood_shift_bwd_p=nullptr;
  std::vector<double> *track_likelihood_shift_bwd_pi=nullptr;
  std::vector<double> *track_likelihood_shift_bwd_k=nullptr;
  std::vector<double> *track_dEdx_perhit_u=nullptr;
  std::vector<double> *track_dEdx_perhit_v=nullptr;
  std::vector<double> *track_dEdx_perhit_y=nullptr;
  std::vector<double> *track_resrange_perhit_u=nullptr;
  std::vector<double> *track_resrange_perhit_v=nullptr;
  std::vector<double> *track_resrange_perhit_y=nullptr;
  std::vector<std::vector<double>> *track_Lmip_perhit=nullptr;
  double track_chi2mu_plane2=-9999;
  double track_chi2p_plane2=-9999;
  double track_chi2pi_plane2=-9999;
  double track_chi2k_plane2=-9999;
  double track_rangeE_mu=-9999;
  double track_rangeE_p=-9999;
  double track_length = -9999;
  double track_theta = -9999;
  double track_phi = -9999;
  double track_theta_x = -9999;
  double true_purity = -9999;
  double true_start_momentum = -9999;
  double true_start_x = -9999;
  double true_start_y = -9999;
  double true_start_z = -9999;
  double true_end_momentum = -9999;
  double true_end_x = -9999;
  double true_end_y = -9999;
  double true_end_z = -9999;
  int track_id = -9999;
  double track_start_x = -9999;
  double track_start_y = -9999;
  double track_start_z = -9999;
  double track_end_x = -9999;
  double track_end_y = -9999;
  double track_end_z = -9999;
  double track_shower_score = -9999;
  bool track_ismuoncandidate = -1;
  std::vector<double> *track_nhits = nullptr;

  // Make the chi2 variables std::vector<doubles> so we can handle them in the same way as the other variables
  // This is just a cheat - we only have chi2 variables for collection plane right now, so set other values to 0 by hand. Fix this in the future!
  std::vector<double> *track_chi2mu = nullptr;
  std::vector<double> *track_chi2p = nullptr;
  std::vector<double> *track_chi2pi = nullptr;
  std::vector<double> *track_chi2k = nullptr;

  // These are derived quantities - derived from the values above in CalcPIDvars
  std::vector<double> *track_likelihood_p;
  std::vector<double> *track_likelihood_mu;
  std::vector<double> *track_likelihood_pi;
  std::vector<double> *track_likelihood_k;
  std::vector<double> *track_likelihood_mip;
  std::vector<double> *track_likelihood_maxmumip;
  std::vector<double> *track_likelihood_muoverp;
  std::vector<double> *track_likelihood_mipoverp;
  std::vector<double> *track_lnlikelihood_mipoverp;
  std::vector<double> *track_likelihood_maxmumipoverp;
  std::vector<double> *track_Lmu_0to1;
  std::vector<double> *track_Lmip_0to1;
  std::vector<double> *track_Lpi_0to1;
  std::vector<double> *track_Lp_0to1;
  std::vector<double> *track_Lk_0to1;
  std::vector<double> *track_Lmumip_0to1;
  std::vector<double> *track_Lmumippi_0to1;
  std::vector<double> *track_Lmumip_0to1_nopionkaon;
  std::vector<double> *track_Lmuovermip;
  std::vector<double> *track_Lmumipoverpi;
  std::vector<double> *track_depE_minus_rangeE_mu;
  std::vector<double> *track_depE_minus_rangeE_p;
  std::vector<double> *track_chi2_muminusp;
  std::vector<double> *track_shift_bestL;
  std::vector<double> *track_Lmip_atstart;
  std::vector<double> *track_lnLmip_atstart;
  std::vector<double> *track_dEdx_mean_atstart;
};

void settreevars(TTree *intree, treevars *varstoset){
  intree->SetBranchAddress("true_PDG"              , &(varstoset->true_PDG));
  intree->SetBranchAddress("track_likelihood_fwd_p"   , &(varstoset->track_likelihood_fwd_p));
  intree->SetBranchAddress("track_likelihood_fwd_mu"  , &(varstoset->track_likelihood_fwd_mu));
  intree->SetBranchAddress("track_likelihood_fwd_pi"  , &(varstoset->track_likelihood_fwd_pi));
  intree->SetBranchAddress("track_likelihood_fwd_k"   , &(varstoset->track_likelihood_fwd_k));
  intree->SetBranchAddress("track_likelihood_fwd_mip" , &(varstoset->track_likelihood_fwd_mip));
  intree->SetBranchAddress("track_likelihood_bwd_p"   , &(varstoset->track_likelihood_bwd_p));
  intree->SetBranchAddress("track_likelihood_bwd_mu"  , &(varstoset->track_likelihood_bwd_mu));
  intree->SetBranchAddress("track_likelihood_bwd_pi"  , &(varstoset->track_likelihood_bwd_pi));
  intree->SetBranchAddress("track_likelihood_bwd_k"   , &(varstoset->track_likelihood_bwd_k));
  intree->SetBranchAddress("track_likelihood_shift_fwd_mu"   , &(varstoset->track_likelihood_shift_fwd_mu));
  intree->SetBranchAddress("track_likelihood_shift_fwd_p"   , &(varstoset->track_likelihood_shift_fwd_p));
  intree->SetBranchAddress("track_likelihood_shift_fwd_pi"   , &(varstoset->track_likelihood_shift_fwd_pi));
  intree->SetBranchAddress("track_likelihood_shift_fwd_k"   , &(varstoset->track_likelihood_shift_fwd_k));
  intree->SetBranchAddress("track_likelihood_shift_bwd_mu"   , &(varstoset->track_likelihood_shift_bwd_mu));
  intree->SetBranchAddress("track_likelihood_shift_bwd_p"   , &(varstoset->track_likelihood_shift_bwd_p));
  intree->SetBranchAddress("track_likelihood_shift_bwd_pi"   , &(varstoset->track_likelihood_shift_bwd_pi));
  intree->SetBranchAddress("track_likelihood_shift_bwd_k"   , &(varstoset->track_likelihood_shift_bwd_k));
  intree->SetBranchAddress("track_PIDA_mean"       , &(varstoset->track_PIDA_mean));
  intree->SetBranchAddress("track_PIDA_kde"        , &(varstoset->track_PIDA_kde));
  intree->SetBranchAddress("track_PIDA_median"     ,&(varstoset->track_PIDA_median));
  intree->SetBranchAddress("track_Chi2Muon", &(varstoset->track_chi2mu));
  intree->SetBranchAddress("track_Chi2Proton", &(varstoset->track_chi2p));
  intree->SetBranchAddress("track_Chi2Pion", &(varstoset->track_chi2pi));
  intree->SetBranchAddress("track_Chi2Kaon", &(varstoset->track_chi2k));
  intree->SetBranchAddress("track_depE", &(varstoset->track_depE));
  intree->SetBranchAddress("track_rangeE_mu", &(varstoset->track_rangeE_mu));
  intree->SetBranchAddress("track_rangeE_p", &(varstoset->track_rangeE_p));
  intree->SetBranchAddress("track_length", &(varstoset->track_length));
  intree->SetBranchAddress("track_Lmip_perhit", &(varstoset->track_Lmip_perhit));
  intree->SetBranchAddress("track_dEdx_perhit_u", &(varstoset->track_dEdx_perhit_u));
  intree->SetBranchAddress("track_dEdx_perhit_v", &(varstoset->track_dEdx_perhit_v));
  intree->SetBranchAddress("track_dEdx_perhit_y", &(varstoset->track_dEdx_perhit_y));
  intree->SetBranchAddress("track_resrange_perhit_u", &(varstoset->track_resrange_perhit_u));
  intree->SetBranchAddress("track_resrange_perhit_v", &(varstoset->track_resrange_perhit_v));
  intree->SetBranchAddress("track_resrange_perhit_y", &(varstoset->track_resrange_perhit_y));
  intree->SetBranchAddress("track_theta", &(varstoset->track_theta));
  intree->SetBranchAddress("track_phi", &(varstoset->track_phi));
  intree->SetBranchAddress("true_purity", &(varstoset->true_purity));
  intree->SetBranchAddress("true_start_momentum", &(varstoset->true_start_momentum));
  intree->SetBranchAddress("true_start_x", &(varstoset->true_start_x));
  intree->SetBranchAddress("true_start_y", &(varstoset->true_start_y));
  intree->SetBranchAddress("true_start_z", &(varstoset->true_start_z));
  intree->SetBranchAddress("true_end_momentum", &(varstoset->true_end_momentum));
  intree->SetBranchAddress("true_end_x", &(varstoset->true_end_x));
  intree->SetBranchAddress("true_end_y", &(varstoset->true_end_y));
  intree->SetBranchAddress("true_end_z", &(varstoset->true_end_z));
  intree->SetBranchAddress("track_id", &(varstoset->track_id));
  intree->SetBranchAddress("track_start_x", &(varstoset->track_start_x));
  intree->SetBranchAddress("track_start_y", &(varstoset->track_start_y));
  intree->SetBranchAddress("track_start_z", &(varstoset->track_start_z));
  intree->SetBranchAddress("track_end_x", &(varstoset->track_end_x));
  intree->SetBranchAddress("track_end_y", &(varstoset->track_end_y));
  intree->SetBranchAddress("track_end_z", &(varstoset->track_end_z));
  intree->SetBranchAddress("track_shower_score", &(varstoset->track_shower_score));
  intree->SetBranchAddress("track_ismuoncandidate", &(varstoset->track_ismuoncandidate));
  intree->SetBranchAddress("track_nhits", &(varstoset->track_nhits));

  intree->GetEntry(0);
  size_t nplanes = varstoset->track_likelihood_fwd_p->size();

  varstoset->track_chi2mu = new std::vector<double>(nplanes);
  varstoset->track_chi2p = new std::vector<double>(nplanes);
  varstoset->track_chi2pi = new std::vector<double>(nplanes);
  varstoset->track_chi2k = new std::vector<double>(nplanes);

  varstoset->track_likelihood_p = new std::vector<double>(nplanes);
  varstoset->track_likelihood_mu = new std::vector<double>(nplanes);
  varstoset->track_likelihood_pi = new std::vector<double>(nplanes);
  varstoset->track_likelihood_k = new std::vector<double>(nplanes);
  varstoset->track_likelihood_mip = new std::vector<double>(nplanes);
  varstoset->track_likelihood_maxmumip = new std::vector<double>(nplanes);
  varstoset->track_likelihood_muoverp = new std::vector<double>(nplanes);
  varstoset->track_likelihood_mipoverp = new std::vector<double>(nplanes);
  varstoset->track_lnlikelihood_mipoverp = new std::vector<double>(nplanes);
  varstoset->track_likelihood_maxmumipoverp = new std::vector<double>(nplanes);
  varstoset->track_Lmu_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lmip_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lpi_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lp_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lk_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lmumip_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lmumippi_0to1 = new std::vector<double>(nplanes);
  varstoset->track_Lmumip_0to1_nopionkaon = new std::vector<double>(nplanes);
  varstoset->track_Lmuovermip = new std::vector<double>(nplanes);
  varstoset->track_Lmumipoverpi = new std::vector<double>(nplanes);
  varstoset->track_depE_minus_rangeE_mu = new std::vector<double>(nplanes);
  varstoset->track_depE_minus_rangeE_p = new std::vector<double>(nplanes);
  varstoset->track_chi2_muminusp = new std::vector<double>(nplanes);
  varstoset->track_shift_bestL = new std::vector<double>(nplanes);
  varstoset->track_Lmip_atstart = new std::vector<double>(nplanes);
  varstoset->track_lnLmip_atstart = new std::vector<double>(nplanes);
  varstoset->track_dEdx_mean_atstart = new std::vector<double>(nplanes);
}

std::pair<double,double> GetChi2(TH1D *o, TH1D* e){

  double chi2 = 0;
  double dof = 0;

  for (int i = 0; i < o->GetNbinsX(); i++){

    if (e->GetBinContent(i) > 0){

      chi2 += std::pow(o->GetBinContent(i) - e->GetBinContent(i),2)/e->GetBinContent(i);

    }
    if (e->GetBinContent(i) != 0 || o->GetBinContent(i) !=0){

      dof++;

    }

  }

  // std::cout << "Chi2 from our calculation: " << chi2 << std::endl;
  // std::cout << "Chi2 from TH1: " << std::endl;
  // double chi2_th1 = o->Chi2Test(e,"NORM UU P CHI2");

  std::pair<double,double> returner;
  returner.first = chi2;
  returner.second = dof;
  return returner;

}

void OverlayChi2(TCanvas *c1, TString e_str, TString o_str){

  TH1D* h_e = (TH1D*)c1->GetPrimitive(e_str.Data());
  TH1D* h_o = (TH1D*)c1->GetPrimitive(o_str.Data());

  std::pair<double,double> chi2 = GetChi2(h_o, h_e);

  TPaveText *pt = new TPaveText(0.15, 0.87, 0.45, 0.92, "NDC");
  TString chi2string = Form("Chi2/NDF: %.2f/%g", chi2.first, chi2.second);
  pt->SetFillColor(kWhite);
  pt->AddText(chi2string.Data());
  pt->SetTextSize(0.05);
  pt->Draw("same");

}

void CalcPIDvars(treevars *vars, bool isScale){
  //std::cout << "Calculating PID variables for " << vars->track_likelihood_fwd_p->size() << " planes" << std::endl;
  for (size_t i_pl=0; i_pl < vars->track_likelihood_fwd_p->size(); i_pl++){
    /*    if (i_pl==0 || i_pl==1){
          vars->track_chi2_muminusp->at(i_pl) = 0;
          vars->track_chi2mu->at(i_pl) = 0;
          vars->track_chi2p->at(i_pl) = 0;
          vars->track_chi2k->at(i_pl) = 0;
          vars->track_chi2pi->at(i_pl) = 0;
          }
          else{
          */
    // vars->track_chi2mu->at(i_pl) = vars->track_chi2mu->at(i_pl);
    // vars->track_chi2p->at(i_pl) = vars->track_chi2p->at(i_pl);
    // vars->track_chi2k->at(i_pl) = vars->track_chi2k->at(i_pl);
    // vars->track_chi2pi->at(i_pl) = vars->track_chi2pi->at(i_pl);
    if (vars->track_chi2mu->at(i_pl)>=0 && vars->track_chi2p->at(i_pl)>=0){
      vars->track_chi2_muminusp->at(i_pl) = vars->track_chi2mu->at(i_pl) - vars->track_chi2p->at(i_pl);
    }
    //}

    //double scalefactor_mu = 0.821;
    //double scalefactor_p  = 0.913;
    double scalefactor_mu = 1.0;
    double scalefactor_p  = 1.0;

    if (!isScale) {
      scalefactor_mu = 1;
      scalefactor_p = 1;
    }

    vars->track_likelihood_p->at(i_pl) = std::max(vars->track_likelihood_fwd_p->at(i_pl)     , vars->track_likelihood_bwd_p->at(i_pl)) * scalefactor_p;
    vars->track_likelihood_mu->at(i_pl) = std::max(vars->track_likelihood_fwd_mu->at(i_pl)    , vars->track_likelihood_bwd_mu->at(i_pl)) * scalefactor_mu;
    vars->track_likelihood_pi->at(i_pl) = std::max(vars->track_likelihood_fwd_pi->at(i_pl)    , vars->track_likelihood_bwd_pi->at(i_pl)) * scalefactor_mu;
    vars->track_likelihood_k->at(i_pl) = std::max(vars->track_likelihood_fwd_k->at(i_pl)     , vars->track_likelihood_bwd_k->at(i_pl)) * scalefactor_mu;
    vars->track_likelihood_mip->at(i_pl) = vars->track_likelihood_fwd_mip->at(i_pl) * scalefactor_mu;
    vars->track_likelihood_maxmumip->at(i_pl) = std::max(vars->track_likelihood_mu->at(i_pl), vars->track_likelihood_mip->at(i_pl));

    // Skip tracks with invalid likelihood scores
    if (vars->track_likelihood_p->at(i_pl)<0 || vars->track_likelihood_mu->at(i_pl)<0 || vars->track_likelihood_k->at(i_pl)<0 || vars->track_likelihood_pi->at(i_pl)<0 || vars->track_likelihood_mip->at(i_pl)<0){
      // Do nothing
    }
    else{
      vars->track_likelihood_muoverp->at(i_pl) = vars->track_likelihood_mu->at(i_pl) / vars->track_likelihood_p->at(i_pl);
      vars->track_likelihood_mipoverp->at(i_pl) = vars->track_likelihood_mip->at(i_pl) / vars->track_likelihood_p->at(i_pl);
      vars->track_likelihood_maxmumipoverp->at(i_pl) = vars->track_likelihood_maxmumip->at(i_pl) / vars->track_likelihood_p->at(i_pl);
      
      vars->track_lnlikelihood_mipoverp->at(i_pl) = TMath::Log(vars->track_likelihood_mipoverp->at(i_pl));

      double denom = (vars->track_likelihood_p->at(i_pl)+vars->track_likelihood_mu->at(i_pl)+vars->track_likelihood_k->at(i_pl)+vars->track_likelihood_pi->at(i_pl)+vars->track_likelihood_mip->at(i_pl));

      vars->track_Lmu_0to1->at(i_pl) = vars->track_likelihood_mu->at(i_pl)/denom;
      vars->track_Lmip_0to1->at(i_pl) = vars->track_likelihood_mip->at(i_pl)/denom;
      vars->track_Lpi_0to1->at(i_pl) = vars->track_likelihood_pi->at(i_pl)/denom;
      vars->track_Lp_0to1->at(i_pl) = vars->track_likelihood_p->at(i_pl)/denom;
      vars->track_Lk_0to1->at(i_pl) = vars->track_likelihood_k->at(i_pl)/denom;
      vars->track_Lmumip_0to1->at(i_pl) = (vars->track_likelihood_mu->at(i_pl)+vars->track_likelihood_mip->at(i_pl))/denom;
      vars->track_Lmumippi_0to1->at(i_pl) = (vars->track_likelihood_mu->at(i_pl)+vars->track_likelihood_mip->at(i_pl)+vars->track_likelihood_pi->at(i_pl))/denom;

      vars->track_Lmumip_0to1_nopionkaon->at(i_pl) = (vars->track_likelihood_mu->at(i_pl)+vars->track_likelihood_mip->at(i_pl))/(vars->track_likelihood_mu->at(i_pl)+vars->track_likelihood_mip->at(i_pl)+vars->track_likelihood_p->at(i_pl));
    }

    vars->track_depE_minus_rangeE_mu->at(i_pl) = vars->track_depE->at(i_pl) - vars->track_rangeE_mu;
    vars->track_depE_minus_rangeE_p->at(i_pl) = vars->track_depE->at(i_pl) - vars->track_rangeE_p;

    vars->track_Lmuovermip->at(i_pl) = vars->track_likelihood_mu->at(i_pl) / vars->track_likelihood_mip->at(i_pl);
    vars->track_Lmumipoverpi->at(i_pl) = std::max(vars->track_likelihood_mu->at(i_pl), vars->track_likelihood_mip->at(i_pl)) / vars->track_likelihood_pi->at(i_pl);

    double maxl = std::max({vars->track_likelihood_fwd_p->at(i_pl),vars->track_likelihood_bwd_p->at(i_pl),vars->track_likelihood_fwd_mu->at(i_pl),vars->track_likelihood_bwd_mu->at(i_pl),vars->track_likelihood_fwd_pi->at(i_pl),vars->track_likelihood_bwd_pi->at(i_pl),vars->track_likelihood_fwd_k->at(i_pl),vars->track_likelihood_bwd_k->at(i_pl)});

    if (maxl == vars->track_likelihood_fwd_p->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_fwd_p->at(i_pl);
    else if (maxl == vars->track_likelihood_bwd_p->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_bwd_p->at(i_pl);
    else if (maxl == vars->track_likelihood_fwd_mu->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_fwd_mu->at(i_pl);
    else if (maxl == vars->track_likelihood_bwd_mu->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_bwd_mu->at(i_pl);
    else if (maxl == vars->track_likelihood_fwd_pi->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_fwd_pi->at(i_pl);
    else if (maxl == vars->track_likelihood_bwd_pi->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_bwd_pi->at(i_pl);
    else if (maxl == vars->track_likelihood_fwd_k->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_fwd_k->at(i_pl);
    else if (maxl == vars->track_likelihood_bwd_k->at(i_pl)) vars->track_shift_bestL->at(i_pl) = vars->track_likelihood_shift_bwd_k->at(i_pl);

    int nhits_start = 10;
    double Lmip_start_mean = 0.;
    double dEdx_start_mean = 0.;
     // std::cout << i_pl << ": " << vars->track_Lmip_perhit->at(i_pl).size() << ", " << vars->track_dEdx_perhit_u->size() << ", " << vars->track_dEdx_perhit_v->size() << ", " << vars->track_dEdx_perhit_y->size() << std::endl;
    std::vector<float> dEdx_float;
    //std::cout << "---" << std::endl;
    for (int i=1; i<nhits_start; i++){
      // Skip first hit (start from i=1) and last hit
      if (i>=vars->track_Lmip_perhit->at(i_pl).size()) continue;
      Lmip_start_mean += vars->track_Lmip_perhit->at(i_pl).at(i);
      if (i_pl==0){
        size_t perhit_size = vars->track_dEdx_perhit_u->size();
        if (i>=perhit_size) continue;
        // dEdx_start_mean += vars->track_dEdx_perhit_u->at(i);
        int index = i;
        if (vars->track_resrange_perhit_u->at(0)<vars->track_resrange_perhit_u->at(perhit_size-1)){ // start of vector is end of track
          index = perhit_size-i;
        }
        dEdx_float.push_back((float)(vars->track_dEdx_perhit_u->at(index)));
        // std::cout << "Pushing back residual range " << vars->track_resrange_perhit_u->at(index) << " instead of " << vars->track_resrange_perhit_u->at(i) << std::endl;
      }
      else if (i_pl==1){
        size_t perhit_size = vars->track_dEdx_perhit_v->size();
        if (i>=perhit_size) continue;
        // dEdx_start_mean += vars->track_dEdx_perhit_v->at(i);
        int index = i;
        if (vars->track_resrange_perhit_v->at(0)<vars->track_resrange_perhit_v->at(perhit_size-1)){ // start of vector is end of track
          index = perhit_size-i;
        }
        dEdx_float.push_back((float)(vars->track_dEdx_perhit_v->at(index)));
        // std::cout << "Pushing back residual range " << vars->track_resrange_perhit_v->at(index) << " instead of " << vars->track_resrange_perhit_v->at(i) << std::endl;
      }
      else if (i_pl==2){
        size_t perhit_size = vars->track_dEdx_perhit_y->size();
        if (i>=perhit_size) continue;
        // dEdx_start_mean += vars->track_dEdx_perhit_y->at(i);
        int index = i;
        if (vars->track_resrange_perhit_y->at(0)<vars->track_resrange_perhit_y->at(perhit_size-1)){ // start of vector is end of track
          index = perhit_size-i;
          // std::cout << "Pushing back residual range " << vars->track_resrange_perhit_y->at(index) << " instead of " << vars->track_resrange_perhit_y->at(i) << std::endl;
        }
        else{
          // std::cout << "Pushing back residual range " << vars->track_resrange_perhit_y->at(index) << " instead of " << vars->track_resrange_perhit_y->at(perhit_size-i) << std::endl;
        }
        dEdx_float.push_back((float)(vars->track_dEdx_perhit_y->at(index)));
      }
    }
    Lmip_start_mean /= (double)nhits_start;
    vars->track_Lmip_atstart->at(i_pl) = Lmip_start_mean;
    vars->track_lnLmip_atstart->at(i_pl) = TMath::Log(Lmip_start_mean);
    // dEdx_start_mean /= (double)nhits_start;
    // vars->track_dEdx_mean_atstart->at(i_pl) = dEdx_start_mean;
    // TruncMean trm;
    // if (dEdx_float.size()>0) vars->track_dEdx_mean_atstart->at(i_pl) = (double)trm.CalcIterativeTruncMean(dEdx_float, 1, 1, 0, 1, 0.1, 1.0);
  }

  // Theta_x: angle w.r.t. x axis (drift direction), between 0 and 90 degrees
  TVector3 vec;
  vec.SetMagThetaPhi(1,vars->track_theta,vars->track_phi);
  TVector3 x_axis(1,0,0);
  double theta_x = vec.Angle(x_axis)*180.0/TMath::Pi();
  if (theta_x>90) theta_x = 180-theta_x;
  vars->track_theta_x = theta_x;

}


// --------------------------------------------------- //
// This struct contains the histograms and all functions related to them

struct hist1D{
  TH1D *h_mu;
  TH1D *h_p;
  TH1D *h_pi;
  TH1D *h_k;
  TH1D *h_other;
  TH1D *h_all;

  TLegend *l;

  // Constructor for this struct of hists
  hist1D(std::string name, std::string title, double nbins, double binlow, double binhigh){
    h_mu = new TH1D(std::string(name+"_mu").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_p = new TH1D(std::string(name+"_p").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_pi = new TH1D(std::string(name+"_pi").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_k = new TH1D(std::string(name+"_k").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_other = new TH1D(std::string(name+"_other").c_str(),title.c_str(),nbins,binlow,binhigh);

    h_all = new TH1D(std::string(name+"_all").c_str(),title.c_str(),nbins,binlow,binhigh);

    h_mu->SetFillColor(TColor::GetColor(8,64,129));
    h_p->SetFillColor(TColor::GetColor(215, 48, 39));
    h_pi->SetFillColor(TColor::GetColor(166,217,106));
    h_k->SetFillColor(TColor::GetColor(133,1,98));
    h_other->SetFillColor(TColor::GetColor(197,197,197));
    // "All" styling for MC (data styling is set in DrawData)
    h_all->SetFillStyle(3345);
    h_all->SetFillColor(kGray+2);
    h_all->SetMarkerSize(0.); // bad hack because root keeps drawing markers and I can't make it stop

    l = new TLegend(0.59,0.64,0.81,0.87);
    l->AddEntry(h_p,"True proton","f");
    l->AddEntry(h_mu,"True muon","f");
    l->AddEntry(h_pi,"True pion","f");
    l->AddEntry(h_k,"True kaon","f");
    l->AddEntry(h_other,"True other","f");
  }
};

struct hist2D{
  TH3D *h2D_mu;
  TH3D *h2D_p;
  TH3D *h2D_pi;
  TH3D *h2D_k;
  TH3D *h2D_other;
  TH3D *h2D_all;

  hist2D(std::string name                                 , std::string title , double nbinsx , double binlowx , double binhighx , double nbinsy , double binlowy , double binhighy, double nbinsz, double binlowz, double binhighz){
    h2D_mu = new TH3D(std::string(name+"_mu").c_str()       , title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy,nbinsz,binlowz,binhighz);
    h2D_p = new TH3D(std::string(name+"_p").c_str()         , title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy,nbinsz,binlowz,binhighz);
    h2D_pi = new TH3D(std::string(name+"_pi").c_str()       , title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy,nbinsz,binlowz,binhighz);
    h2D_k = new TH3D(std::string(name+"_k").c_str()         , title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy,nbinsz,binlowz,binhighz);
    h2D_other = new TH3D(std::string(name+"_other").c_str() , title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy,nbinsz,binlowz,binhighz);
    h2D_all = new TH3D(std::string(name+"_all").c_str(),title.c_str(),nbinsx,binlowx,binhighx,nbinsy, binlowy,binhighy,nbinsz,binlowz,binhighz);

  }

};

void FillHist(hist1D *hists, double value, int pdg){
  // Fill "all" histogram for every entry
  hists->h_all->Fill(value);

  // Now fill histograms by particle type
  if (TMath::Abs(pdg)==13){ // muon
    hists->h_mu->Fill(value);
  }
  else if (TMath::Abs(pdg)==2212){ // proton
    hists->h_p->Fill(value);
  }
  else if (TMath::Abs(pdg)==211){ // pion
    hists->h_pi->Fill(value);
  }
  else if (TMath::Abs(pdg)==321){ // kaon
    hists->h_k->Fill(value);
  }
  else{ // other
    hists->h_other->Fill(value);
  }
}

void Fill2DHist(hist2D *hists, double valuex, double valuey, double valuez, int pdg){

  // Fill "all" histogram for every entry
  hists->h2D_all->Fill(valuex, valuey, valuez);

  // Now fill histograms by particle type
  if (TMath::Abs(pdg)==13){ // muon
    hists->h2D_mu->Fill(valuex, valuey, valuez);
  }
  else if (TMath::Abs(pdg)==2212){ // proton
    hists->h2D_p->Fill(valuex, valuey, valuez);
  }
  else if (TMath::Abs(pdg)==211){ // pion
    hists->h2D_pi->Fill(valuex, valuey, valuez);
  }
  else if (TMath::Abs(pdg)==321){ // kaon
    hists->h2D_k->Fill(valuex, valuey, valuez);
  }
  else{ // other
    hists->h2D_other->Fill(valuex, valuey, valuez);
  }
}

void DrawMC(hist1D *hists, double POTScaling, double yrange){
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);

  if (POTScaling == 0.){ // area normalise
    POTScaling = 1./hists->h_all->Integral();
    hists->h_all->GetYaxis()->SetTitle("No. tracks (area normalised)");
  }
  else hists->h_all->GetYaxis()->SetTitle("No. tracks (POT normalised)");

  hists->h_mu->Sumw2();
  hists->h_p->Sumw2();
  hists->h_pi->Sumw2();
  hists->h_k->Sumw2();
  hists->h_other->Sumw2();
  hists->h_all->Sumw2();

  hists->h_mu->Scale(POTScaling);
  hists->h_p->Scale(POTScaling);
  hists->h_pi->Scale(POTScaling);
  hists->h_k->Scale(POTScaling);
  hists->h_other->Scale(POTScaling);
  hists->h_all->Scale(POTScaling);

  std::cout << "h_all MC->Integral() = " << hists->h_all->Integral() << std::endl;

  THStack *hs = new THStack("hs","hs");
  hs->Add(hists->h_p);
  hs->Add(hists->h_mu);
  hs->Add(hists->h_pi);
  hs->Add(hists->h_k);
  hs->Add(hists->h_other);

  if (yrange == -999){
    hists->h_all->SetMaximum((hists->h_all->GetMaximum())*1.2);
  }
  else{
    hists->h_all->SetMaximum(yrange);
  }
  hists->h_all->SetMinimum(0);
  hists->h_all->Draw("hist"); // Draw this one first because it knows about the axis titles
  hs->Draw("same hist");
  TH1D *h_err = (TH1D*)hists->h_all->Clone("h_err");
  h_err->Draw("same E2"); // Draw it again so errors are on top

  bool legendleft=false;
  int nbinsx = hists->h_all->GetXaxis()->GetNbins();
  for (int i_bin=(int)(nbinsx/2); i_bin < nbinsx+1; i_bin++){
    if (hists->h_all->GetBinContent(i_bin)>0.6*hists->h_all->GetMaximum()){
      legendleft=true;
    }
  }
  if (legendleft){
    hists->l->SetX1(0.19);
    hists->l->SetX2(0.41);
  }
  hists->l->Draw();
}

void DrawMC2D(hist2D *hists, double POTScaling, double yrange){
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);

  if (POTScaling == 0.){ // area normalise
    POTScaling = 1./hists->h2D_all->Integral();
    //hists->h2D_all->GetYaxis()->SetTitle("No. tracks (area normalised)");
  }
 // else hists->h2D_all->GetYaxis()->SetTitle("No. tracks (POT normalised)");

  //hists->h2D_mu->Sumw2();
  //hists->h2D_p->Sumw2();
  //hists->h2D_pi->Sumw2();
  //hists->h2D_k->Sumw2();
  //hists->h2D_other->Sumw2();
  hists->h2D_all->Sumw2();

  //hists->h2D_mu->Scale(POTScaling);
  //hists->h2D_p->Scale(POTScaling);
  //hists->h2D_pi->Scale(POTScaling);
  //hists->h2D_k->Scale(POTScaling);
  //hists->h2D_other->Scale(POTScaling);
  hists->h2D_all->Scale(POTScaling);

  std::cout << "h2D_all MC->Integral() = " << hists->h2D_all->Integral() << std::endl;

  //THStack *hs = new THStack("hs","hs");
  //hs->Add(hists->h_p);
  //hs->Add(hists->h_mu);
  //hs->Add(hists->h_pi);
  //hs->Add(hists->h_k);
  //hs->Add(hists->h_other);

  //if (yrange == -999){
  //  hists->h_all->SetMaximum((hists->h_all->GetMaximum())*1.2);
  //}
  //else{
  //  hists->h_all->SetMaximum(yrange);
  //}
  //hists->h_all->SetMinimum(0);

  hists->h2D_all->Draw("scat"); // Draw this one first because it knows about the axis titles
  //hists->h2D_other->Draw("same");
  //hists->h2D_k->Draw("same");
  //hists->h2D_pi->Draw("same");
  //hists->h2D_p->Draw("same");
  //hists->h2D_mu->Draw("same");

  //bool legendleft=false;
  //int nbinsx = hists->h2D_all->GetXaxis()->GetNbins();
  //for (int i_bin=(int)(nbinsx/2); i_bin < nbinsx+1; i_bin++){
  //  if (hists->h2D_all->GetBinContent(i_bin)>0.6*hists->h2D_all->GetMaximum()){
  //    legendleft=true;
  //  }
  //}
  //if (legendleft){
  //  hists->l->SetX1(0.19);
  //  hists->l->SetX2(0.41);
  //}
  //hists->l->Draw();
}


void DrawMCPlusOffbeam(hist1D *hists, hist1D *offbeam, double POTScaling, double OffBeamScaling, double yrange){
  // Note that there are no area-normalised options here because I'm not sure that makes sense
  hists->h_all->GetYaxis()->SetTitle("No. tracks (POT normalised)");

  hists->h_mu->Sumw2();
  hists->h_p->Sumw2();
  hists->h_pi->Sumw2();
  hists->h_k->Sumw2();
  hists->h_other->Sumw2();
  hists->h_all->Sumw2();
  offbeam->h_all->Sumw2();

  hists->h_mu->Scale(POTScaling);
  hists->h_p->Scale(POTScaling);
  hists->h_pi->Scale(POTScaling);
  hists->h_k->Scale(POTScaling);
  hists->h_other->Scale(POTScaling);
  hists->h_all->Scale(POTScaling);

  offbeam->h_all->Scale(OffBeamScaling);
  offbeam->h_all->SetFillColor(kBlack);
  offbeam->h_all->SetFillStyle(3345);
  offbeam->h_all->SetLineColor(kBlack);

  THStack *hs = new THStack("hs","hs");
  hs->Add(offbeam->h_all);
  hs->Add(hists->h_p);
  hs->Add(hists->h_mu);
  hs->Add(hists->h_pi);
  hs->Add(hists->h_k);
  hs->Add(hists->h_other);

  TH1D *h_err = (TH1D*)offbeam->h_all->Clone("h_err");
  h_err->Add(hists->h_p);
  h_err->Add(hists->h_mu);
  h_err->Add(hists->h_pi);
  h_err->Add(hists->h_k);
  h_err->Add(hists->h_other);

  h_err->SetFillColor(kBlack);
  h_err->SetFillStyle(3345);

  if (yrange == -999){
    hists->h_all->SetMaximum((hists->h_all->GetMaximum()+offbeam->h_all->GetMaximum())*1.2);
  }
  else{
    hists->h_all->SetMaximum(yrange);
  }
  hists->h_all->SetMinimum(0);
  hists->h_all->Draw("hist"); // Draw this one first because it knows about the axis titles
  hs->Draw("same hist");
  h_err->Draw("same E2"); // Draw it again so errors are on top

  bool legendleft=false;
  int nbinsx = hists->h_all->GetXaxis()->GetNbins();
  for (int i_bin=(int)(nbinsx/2); i_bin < nbinsx+1; i_bin++){
    if (hists->h_all->GetBinContent(i_bin)>0.6*hists->h_all->GetMaximum()){
      legendleft=true;
    }
  }
  if (legendleft){
    //std::cout << "Putting legend on the left" << std::endl;
    hists->l->SetX1(0.19);
    hists->l->SetX2(0.41);
  }
  hists->l->AddEntry(offbeam->h_all,"Data (off-beam)","f");
  hists->l->Draw();
}

void OverlayOnMinusOffData(TCanvas *c, hist1D *onbeam, hist1D *offbeam, double OffBeamScaling, double POTScaling){
  TH1D *h_onminusoff = (TH1D*)onbeam->h_all->Clone();

  h_onminusoff->Sumw2();
  offbeam->h_all->Sumw2();

  h_onminusoff->Add(offbeam->h_all,-1.0*OffBeamScaling);
  if (POTScaling==0){
    h_onminusoff->Scale(1.0/h_onminusoff->Integral());
  }

  std::cout << "h_onminusoff->Integral() = " << h_onminusoff->Integral() << std::endl;

  h_onminusoff->SetMarkerStyle(20);
  h_onminusoff->SetMarkerSize(0.6);

  c->cd();
  h_onminusoff->Draw("same p E1");

  TLegend *l = (TLegend*)c->GetPrimitive("TPave");
  l->AddEntry(h_onminusoff,"Data (on-off beam)","lp");
}

void OverlayOnBeamData(TCanvas *c, hist1D *onbeam){

  onbeam->h_all->SetMarkerStyle(20);
  onbeam->h_all->SetMarkerSize(0.6);

  c->cd();
  onbeam->h_all->Draw("same p E1");

  TLegend *l = (TLegend*)c->GetPrimitive("TPave");
  l->AddEntry(onbeam->h_all,"Data (on-beam)","lp");
}

/**
 * Function to perform template fit and draw result to canvas
 */
void TemplateFit(hist1D* mchists, hist1D* onb_hists, hist1D* offb_hists, double offbeamscaling, double POTscaling, double yrange, bool onminusoffbeam=false){

  /** Scale offbeam to onbeam */
  offb_hists->h_all->Sumw2();
  // offb_hists->h_all->Scale(offbeamscaling);

  /** create on-beam minus off-beam */
  /** do this for the template fit, whether you want the final plot to show on-off beam or not */
  TH1D* data = (TH1D*)onb_hists->h_all->Clone("data");
  data->Add(offb_hists->h_all, -1.*offbeamscaling);

  /** check kaons exist in sample */
  bool isK = true;
  if (mchists->h_k->Integral() == 0) isK = false;

  /** get current fractions */
  double mufrac = mchists->h_mu->Integral()/mchists->h_all->Integral();
  double pfrac = mchists->h_p->Integral()/mchists->h_all->Integral();
  double pifrac = mchists->h_pi->Integral()/mchists->h_all->Integral();
  double kfrac = 0;
  if (isK) kfrac = mchists->h_k->Integral()/mchists->h_all->Integral();
  double otherfrac = mchists->h_other->Integral()/mchists->h_all->Integral();
  double fracfloat = 0.5;

  /** scale mchists to POT normally */
  mchists->h_mu->Sumw2();
  mchists->h_p->Sumw2();
  mchists->h_pi->Sumw2();
  if (isK) mchists->h_k->Sumw2();
  mchists->h_other->Sumw2();

  mchists->h_mu->Scale(POTscaling);
  mchists->h_p->Scale(POTscaling);
  mchists->h_pi->Scale(POTscaling);
  if (isK) mchists->h_k->Scale(POTscaling);
  mchists->h_other->Scale(POTscaling);

  /** put mchists objects in to TObjArray */
  TObjArray *mc = new TObjArray(4);
  mc->Add(mchists->h_mu);
  mc->Add(mchists->h_p);
  mc->Add(mchists->h_pi);
  mc->Add(mchists->h_other);
  if (isK) mc->Add(mchists->h_k);

  /** call fitting function */
  TFractionFitter* fit = new TFractionFitter(data, mc);

  fit->Constrain(0, mufrac*fracfloat, mufrac*(1+fracfloat));
  fit->Constrain(1, pfrac*fracfloat, pfrac*(1+fracfloat));
  fit->Constrain(2, pifrac*fracfloat, pifrac*(1+fracfloat));
  fit->Constrain(3, otherfrac*fracfloat, otherfrac*(1+fracfloat));
  if (isK) fit->Constrain(4, kfrac*fracfloat, kfrac*(1+fracfloat));

  int stat = fit->Fit();

  double result_mu = 0;
  double result_mu_err = 0;
  fit->GetResult(0,result_mu, result_mu_err);

  double result_p = 0;
  double result_p_err = 0;
  fit->GetResult(1,result_p, result_p_err);

  double result_pi = 1;
  double result_pi_err = 0;
  fit->GetResult(2,result_pi, result_pi_err);

  double result_other = 0;
  double result_other_err = 0;
  fit->GetResult(3,result_other, result_other_err);

  double result_k = 0;
  double result_k_err = 0;
  if (isK) fit->GetResult(4,result_k, result_k_err);

  std::cout << "[TemplateFit] mu    : " << result_mu << "+/-" << result_mu_err << std::endl;
  std::cout << "[TemplateFit] p     : " << result_p << "+/-" << result_p_err << std::endl;
  std::cout << "[TemplateFit] pi    : " << result_pi << "+/-" << result_pi_err << std::endl;
  std::cout << "[TemplateFit] k     : " << result_k << "+/-" << result_k_err << std::endl;
  std::cout << "[TemplateFit] other : " << result_other << "+/-" << result_other_err << std::endl;
  // std::cout << "[TemplateFit] chi2 = " << fit->GetChisquare() << std::endl;

  /**
   * Apply scaling directly to histograms
   * Modify mchists histograms directly so that they can be drawn using DrawMC or DrawMCPlusOffbeam
   */

  mchists->h_p->Sumw2();
  mchists->h_mu->Sumw2();
  mchists->h_pi->Sumw2();
  if (isK) mchists->h_k->Sumw2();
  mchists->h_other->Sumw2();
  mchists->h_all->Sumw2();

  mchists->h_mu->Scale(result_mu*data->Integral()/mchists->h_mu->Integral());
  mchists->h_p->Scale(result_p* (data->Integral()/mchists->h_p->Integral()));
  mchists->h_pi->Scale(result_pi* (data->Integral()/mchists->h_pi->Integral()));
  if (isK) mchists->h_k->Scale(result_k* (data->Integral()/mchists->h_k->Integral()));
  mchists->h_other->Scale(result_other* (data->Integral()/mchists->h_other->Integral()));

  mchists->h_all->Reset();
  mchists->h_all->Add(mchists->h_mu);
  mchists->h_all->Add(mchists->h_p);
  mchists->h_all->Add(mchists->h_pi);
  mchists->h_all->Add(mchists->h_k);
  mchists->h_all->Add(mchists->h_other);

  /*TH1D *hTotal = (TH1D*)mchists->h_mu->Clone("hTotal");
  hTotal->Add(mchists->h_p);
  hTotal->Add(mchists->h_pi);
  hTotal->Add(mchists->h_k);
  hTotal->Add(mchists->h_other);*/


  /**
   * draw result: because data histograms modified directly just
   * draw using the function rather than drawing here (so this is commented out)
   */

/*
  hTotal->SetFillColor(kBlack);
  hTotal->SetFillStyle(3345);
  hTotal->SetMarkerSize(0);

  THStack* hs = new THStack();
  hs->Add(mchists->h_p);
  hs->Add(mchists->h_mu);
  hs->Add(mchists->h_pi);
  if (isK) hs->Add(mchists->h_k);
  hs->Add(mchists->h_other);

  if (yrange == -999){
    mchists->h_all->SetMaximum(data->GetMaximum()*1.3);
  }
  else{
    mchists->h_all->SetMaximum(yrange);
  }
  mchists->h_all->SetFillColor(kWhite);
  mchists->h_all->SetLineColor(kWhite);
  mchists->h_all->SetMarkerSize(0.);
  mchists->h_all->Draw();
  hs->Draw("histsame");
  hTotal->Draw("e2same");
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.6);
  data->Draw("same");

  std::pair<double,double> chi2 = GetChi2(data, hTotal);

  TPaveText *pt = new TPaveText(0.15, 0.87, 0.45, 0.92, "NDC");
  TString chi2string = Form("Chi2/NDF: %.2f/%g", chi2.first, chi2.second);
  pt->SetFillColor(kWhite);
  pt->AddText(chi2string.Data());
  pt->SetTextSize(0.05);
  pt->Draw("same");

  bool legendleft=false;
  int nbinsx = hTotal->GetXaxis()->GetNbins();
  for (int i_bin=(int)(nbinsx/2); i_bin < nbinsx+1; i_bin++){
    if (hTotal->GetBinContent(i_bin)>0.6*hTotal->GetMaximum()){
      legendleft=true;
    }
  }
  if (legendleft){
    mchists->l->SetX1(0.19);
    mchists->l->SetX2(0.41);
  }
  mchists->l->AddEntry(data,"Data (on beam - off beam)","lp");
  mchists->l->Draw();*/
}

void DrawMCEffPur2D(TCanvas *c_eff, TCanvas *c_pur, TCanvas *c_effpur, hist2D *hists, bool MIPlow, double cut_value, double twodvar1, double twodvar2, std::vector<string> var_names, TFile *fout = nullptr){

  if (cut_value == -999) return;

  c_eff->Divide(2,2,0.0005,0.0005);
  c_pur->Divide(2,2,0.0005,0.0005);
  c_effpur->Divide(2,2,0.0005,0.0005);

  c_eff->SetLogz();
  c_pur->SetLogz();
  c_effpur->SetLogz();

  std::vector<TH3D*> histstoeval = {
    hists->h2D_mu,
    hists->h2D_pi,
    hists->h2D_p
  };

  std::vector<std::string> histtitles = {
    "True muons",
    "True pions",
    "True protons"
  };


  for (int i_h=0; i_h<histstoeval.size(); i_h++){
    TH2D *heff = nullptr;
    TH2D *hpur = nullptr;
    TH2D *heffpur = nullptr;

    //heff->SetTitle(histtitles.at(i_h).c_str());

    int binVal = -999;
    double temp = 1000;
    for (int i_bin=1; i_bin <= histstoeval.at(i_h)->GetZaxis()->GetNbins(); i_bin++){

      if (std::abs(histstoeval.at(i_h)->GetZaxis()->GetBinLowEdge(i_bin) - cut_value) < temp){
        temp = std::abs(histstoeval.at(i_h)->GetZaxis()->GetBinLowEdge(i_bin) - cut_value);
        binVal = i_bin;
      }

    }

    double xlow = histstoeval.at(i_h)->GetXaxis()->GetBinLowEdge(1);
    double xhigh = histstoeval.at(i_h)->GetXaxis()->GetBinLowEdge(histstoeval.at(i_h)->GetNbinsX()+1);
    double ylow = histstoeval.at(i_h)->GetYaxis()->GetBinLowEdge(1);
    double yhigh = histstoeval.at(i_h)->GetYaxis()->GetBinLowEdge(histstoeval.at(i_h)->GetNbinsY()+1);
    double zlow = histstoeval.at(i_h)->GetZaxis()->GetBinLowEdge(1);
    double zhigh = histstoeval.at(i_h)->GetZaxis()->GetBinLowEdge(histstoeval.at(i_h)->GetNbinsZ()+1);

    TH2D* eff_denom = nullptr;
    TH2D* eff_num   = nullptr;
    TH2D* pur_denom = nullptr;
    TH2D* pur_num   = nullptr;

    double zlow_all = hists->h2D_all->GetZaxis()->GetBinLowEdge(0);
    double zhigh_all = hists->h2D_all->GetZaxis()->GetBinLowEdge(hists->h2D_all->GetNbinsZ()+1);

   if ((MIPlow && i_h < 2) || (!MIPlow && i_h ==2)){ // integrate from the bottom

      histstoeval.at(i_h)->GetZaxis()->SetRangeUser(zlow, zhigh);
      eff_denom = (TH2D*)histstoeval.at(i_h)->Project3D("xy")->Clone("eff_denom");

      histstoeval.at(i_h)->GetZaxis()->SetRangeUser(zlow, cut_value);
      eff_num = (TH2D*)histstoeval.at(i_h)->Project3D("xy")->Clone("eff_num");
      pur_num = (TH2D*)histstoeval.at(i_h)->Project3D("xy")->Clone("pur_num");

      hists->h2D_all->GetZaxis()->SetRangeUser(zlow_all, cut_value);
      pur_denom = (TH2D*)hists->h2D_all->Project3D("xy")->Clone("pur_denom");

      hists->h2D_all->GetZaxis()->SetRangeUser(zlow_all, zhigh_all);

    }
    else{ // integrate up to the top

      histstoeval.at(i_h)->GetZaxis()->SetRangeUser(zlow, zhigh);
      eff_denom = (TH2D*)histstoeval.at(i_h)->Project3D("xy")->Clone("eff_denom");

      histstoeval.at(i_h)->GetZaxis()->SetRangeUser(cut_value, zhigh);
      eff_num = (TH2D*)histstoeval.at(i_h)->Project3D("xy")->Clone("eff_num");
      pur_num = (TH2D*)histstoeval.at(i_h)->Project3D("xy")->Clone("pur_num");

      hists->h2D_all->GetZaxis()->SetRangeUser(cut_value, zhigh_all);
      pur_denom = (TH2D*)hists->h2D_all->Project3D("xy")->Clone("pur_denom");

      hists->h2D_all->GetZaxis()->SetRangeUser(zlow_all, zhigh_all);

    }

    TH2D* eff_s = (TH2D*)eff_num->Clone("eff_s");
    eff_s->Divide(eff_denom);

    TH2D* pur_s = (TH2D*)pur_num->Clone("pur_s");
    pur_s->Divide(pur_denom);

    TH2D* effpur_s = (TH2D*)eff_s->Clone("effpur_s");
    effpur_s->Multiply(pur_s);

    c_eff->cd(i_h+1);
    eff_s->SetTitle(std::string(histtitles.at(i_h)+";"+var_names.at(1)+";"+var_names.at(0)).c_str());
    eff_s->GetZaxis()->SetRangeUser(0,1);
    eff_s->Draw("colz");

    c_pur->cd(i_h+1);
    pur_s->SetTitle(std::string(histtitles.at(i_h)+";"+var_names.at(1)+";"+var_names.at(0)).c_str());
    pur_s->GetZaxis()->SetRangeUser(0,1);
    pur_s->Draw("colz");

    c_effpur->cd(i_h+1);
    effpur_s->SetTitle(std::string(histtitles.at(i_h)+";"+var_names.at(1)+";"+var_names.at(0)).c_str());
    effpur_s->Draw("colz");

  }

}


void DrawMCEffPur(TCanvas *c, hist1D *hists, bool MIPlow, double cutval, TFile *fout = nullptr){
  std::vector<TH1D*> histstoeval = {
    hists->h_mu,
    hists->h_pi,
    hists->h_p,
    hists->h_k
  };

  std::vector<std::string> histtitles = {
    "True muons",
    "True pions",
    "True protons",
    "True kaons"
  };

  gStyle->SetOptStat(0);

  c->Divide(2,2,0.0005,0.0005);

  TLegend *l;

  // Print out efficiency and purity at a given cut value to a txt file
  ofstream outtxtfile;
  outtxtfile.open("effpuratgivencutvals.txt",std::ios::app);
  if (!outtxtfile.is_open()) std::cout << "[WARNING] Unable to open txt file to save muon/proton efficiency and purity at given cut values." << std::endl;

  for (int i_h=0; i_h<histstoeval.size(); i_h++){
    TH1D *heff = (TH1D*)hists->h_all->Clone("heff");
    TH1D *hpur = (TH1D*)hists->h_all->Clone("hpur");
    TH1D *heffpur = (TH1D*)hists->h_all->Clone("heffpur");

    heff->Clear();
    hpur->Clear();
    heffpur->Clear();

    heff->SetTitle(histtitles.at(i_h).c_str());

    double legendy=0.;
    double legendpts=0.;

    for (int i_bin=1; i_bin <= histstoeval.at(i_h)->GetXaxis()->GetNbins(); i_bin++){

      double eff, pur;//, efferr, purerr;
      if ((MIPlow && i_h != 2) || (!MIPlow && i_h ==2)){ // integrate from the bottom
      // Assume cut is placed at the high edge of the bin (so includes this bin if integrating down)
        double selected_i = histstoeval.at(i_h)->Integral(0,i_bin);
        double total_i = histstoeval.at(i_h)->Integral(0,histstoeval.at(i_h)->GetXaxis()->GetNbins()+1);
        double selected_all = hists->h_all->Integral(0,i_bin);

        eff = selected_i/total_i;
        pur = selected_i/selected_all;

        if (selected_i==0 && selected_all==0) pur = 0;
        if (total_i==0) eff = 0;
      }
      else{ // integrate up to the top
      // Assume cut is placed at the low edge of the bin (so doesn't include this bin if integrating up)
        double selected_i = histstoeval.at(i_h)->Integral(i_bin+1,histstoeval.at(i_h)->GetXaxis()->GetNbins()+1);
        double total_i = histstoeval.at(i_h)->Integral(0,histstoeval.at(i_h)->GetXaxis()->GetNbins()+1);
        double selected_all = hists->h_all->Integral(i_bin+1,histstoeval.at(i_h)->GetXaxis()->GetNbins()+1);

        eff = selected_i/total_i;
        pur = selected_i/selected_all;

        if (selected_i==0 && selected_all==0) pur = 0;
        if (total_i==0) eff = 0;
      }

      double effpur = eff*pur;
      heff->SetBinContent(i_bin,eff);
      hpur->SetBinContent(i_bin,pur);
      heffpur->SetBinContent(i_bin,effpur);

      if (heff->GetXaxis()->GetBinCenter(i_bin)>0.7*heff->GetXaxis()->GetBinUpEdge(heff->GetXaxis()->GetLast())){
        legendy+=eff;
        legendpts++;
      }
    }

    heff->SetLineColor(kRed);
    heff->SetMarkerColor(kRed);
    heff->SetMarkerStyle(20);
    heff->SetMarkerSize(.3);

    hpur->SetLineColor(kBlue);
    hpur->SetMarkerColor(kBlue);
    hpur->SetMarkerStyle(20);
    hpur->SetMarkerSize(.3);

    heffpur->SetLineColor(kBlack);
    heffpur->SetMarkerColor(kBlack);
    heffpur->SetMarkerStyle(20);
    heffpur->SetMarkerSize(.3);

    heff->GetYaxis()->SetRangeUser(0,1.1);

    c->cd(i_h+1);
    heff->Draw("p");
    hpur->Draw("same p");
    heffpur->Draw("same p");

    if (i_h==0){
      // if (MIPlow) l = new TLegend(0.5,0.3,0.88,0.6);
      // else l = new TLegend(0.5,0.6,0.88,0.88);
      legendy = legendy/legendpts;
      if (legendy<0.6) l = new TLegend(0.5,0.6,0.88,0.88);
      else l = new TLegend(0.5,0.3,0.88,0.6);
      //if (legendy>0.6) l = new TLegend(0.5,(legendy*0.9)-0.3,0.88,legendy*0.9);
      //else l = new TLegend(0.5,legendy*1.1,0.88,legendy*1.1+0.3);
      l->SetTextFont(132);
      l->SetLineColor(kWhite);
      l->SetFillColor(kWhite);
      l->AddEntry(heff,"Efficiency","p");
      l->AddEntry(hpur,"Purity","p");
      l->AddEntry(heffpur,"Efficiency #times Purity","p");
    }

    c->cd(2);
    l->Draw();

    if (fout){
      fout->cd();
      TString name = TString(histstoeval.at(i_h)->GetName());
      name.Remove(0,9);
      heff->Write(TString::Format("heff_%s",name.Data()).Data());
      hpur->Write(TString::Format("hpur_%s",name.Data()).Data());
      heffpur->Write(TString::Format("heffpur_%s",name.Data()).Data());
    }
    // Print out efficiency and purity at a given cut value
    if (cutval!=-999){
      if (histtitles.at(i_h)=="True muons"){
        outtxtfile << "\n" << hists->h_mu->GetName() << "\n";
        outtxtfile << "Cut value = " << cutval << "\n";
        outtxtfile << "    Muon Efficiency = " << heff->Interpolate(cutval) << "\n";
        outtxtfile << "         Purity = " << hpur->Interpolate(cutval) << "\n";
        outtxtfile << "         Eff*Pur = " << heffpur->Interpolate(cutval) << "\n";
      }
      else if (histtitles.at(i_h)=="True protons"){
        outtxtfile << "    Proton Efficiency = " << heff->Interpolate(cutval) << "\n";
        outtxtfile << "           Purity = " << hpur->Interpolate(cutval) << "\n";
        outtxtfile << "           Eff*Pur = " << heffpur->Interpolate(cutval) << "\n";
      }
    }

  }
  outtxtfile.close();
}


//PIP'S CHANGES
//---------------------------------
bool IsContained(treevars vars)
{

  //TPC boundary + FV cut
  double FVxmin = 10;
  double FVxmax = 245;
  double FVymin = -100;
  double FVymax = 100;
  double FVzmin = 10;
  double FVzmax = 1030;
  
  ////Cut out additional section of mostly dead wires in z
  double deadzmin = 675.1;
  double deadzmax = 775.1;

  bool contained = true;
  if ((vars.track_start_x > FVxmax) || (vars.track_start_x < FVxmin)) contained = false;
  if ((vars.track_start_y > FVymax) || (vars.track_start_y < FVymin)) contained = false;
  if ((vars.track_start_z > FVzmax) || (vars.track_start_z < FVzmin)) contained = false;
 
  if ((vars.track_end_x > FVxmax) || (vars.track_end_x < FVxmin)) contained = false;
  if ((vars.track_end_y > FVymax) || (vars.track_end_y < FVymin)) contained = false;
  if ((vars.track_end_z > FVzmax) || (vars.track_end_z < FVzmin)) contained = false;

  //if ((vars.track_end_z < deadzmax) && (vars.track_end_z > deadzmin)) contained = false;

  return contained;
}

double ThetaXZ(treevars vars)
{
  double x = vars.track_end_x - vars.track_start_x;
  double z = vars.track_end_z - vars.track_start_z;

  TVector3 trackvec(x, 0, z);
  trackvec = trackvec.Unit();
  TVector3 zaxis(0, 0, 1);
  double costhetaxz = trackvec.Dot(zaxis);
  double thetaxz = TMath::ACos(costhetaxz);

  trackvec.Delete();
  zaxis.Delete();
  return thetaxz;
}

double ThetaYZ(treevars vars)
{
  double y = vars.track_end_y - vars.track_start_y;
  double z = vars.track_end_z - vars.track_start_z;

  TVector3 trackvec(0, y, z);
  trackvec = trackvec.Unit();
  TVector3 zaxis(0, 0, 1);
  double costhetayz = trackvec.Dot(zaxis);
  double thetayz = TMath::ACos(costhetayz);

  trackvec.Delete();
  zaxis.Delete();
  return thetayz;
}

void FormatYourPlotsAndy(hist1D *sliceplot, int slice)
{
  int rr_floor = 0 + slice*5;
  int rr_ceil  = (slice+1)*5;
  std::string title = std::string(std::to_string(rr_floor) + std::string(" < Residual Range < ") + std::to_string(rr_ceil) + std::string(" cm"));
  sliceplot->h_all->SetTitle(title.c_str());
  sliceplot->h_all->GetXaxis()->SetTitle("dE/dx per hit (all hits)");
}

void FormatPhiSlices(hist1D *sliceplot, int slice)
{
  double phi_floor = 0 + slice*0.628;
  double phi_ceil  = (slice+1)*0.628;
  std::string title = std::string(std::to_string(phi_floor) + std::string(" < #phi < ") + std::to_string(phi_ceil) + std::string(" rad"));
  sliceplot->h_all->SetTitle(title.c_str());
  sliceplot->h_all->GetXaxis()->SetTitle("dE/dx per hit (all hits)");
}

double FunctionFromGraph(TGraph *graph, double *x, double *)
{
  return graph->Eval(x[0]);
}

float Mean(std::vector<float> v)
{
  float mean = 0.;
  for (auto const& n : v) mean += n;
  mean /= v.size();
  
  return mean;
}

float Median(std::vector<float> v)
{
  if (v.size() == 1) return v[0];
  std::vector<float> vcpy = v;
  std::sort(vcpy.begin(), vcpy.end());
  float median = vcpy[ vcpy.size() / 2 ];

  return median;
}

float RMS(std::vector<float> v)
{
  float avg = 0.;
  for (auto const& val : v) avg += val;
  avg /= v.size();
  float rms = 0.;
  for (auto const& val : v) rms += (val-avg)*(val-avg);
  rms = sqrt( rms / ( v.size() -  1 ) );
 
  return rms;
}

static const float  kINVALID_FLOAT = std::numeric_limits<float>::max();

double TruncMeandEdx(std::vector<float> v, const size_t & currentiteration, const float & oldmed = kINVALID_FLOAT)
{
  const size_t nmin = 1;
  const size_t nmax = 1;
  const size_t lmin = 1;
  const float convergencelimit = 0.1;
  const float nsigma = 1;


  float const& mean = Mean(v);
  float const& med  = Median(v);
  float const& rms  = RMS(v);
  
  // if the vector length is below the lower limit -> return
  if (v.size() < lmin) 
  {
    return mean;
  }
  
  // if we have passed the maximum number of iterations -> return
  if (currentiteration >= nmax) 
  {
    return mean;
  }
      
  // if we passed the minimum number of iterations and the mean is close enough to the old value
  float fracdiff = fabs(med-oldmed) / oldmed;
  if ( (currentiteration >= nmin) && (fracdiff < convergencelimit) ) 
  {
    return mean;
  }
   
  // if reached here it means we have to go on for another iteration
  // cutoff tails of distribution surrounding the mean
  // use erase-remove : https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
  // https://stackoverflow.com/questions/17270837/stdvector-removing-elements-which-fulfill-some-conditions
  
  v.erase( std::remove_if( v.begin(), v.end(), [med,nsigma,rms](const float& x) { return ( (x < (med-nsigma*rms)) || (x > (med+nsigma*rms)) ); }), v.end());

  return TruncMeandEdx(v, currentiteration+1, med);
}

//---------------------------------


// void SaveHists(hist1D *hists, TFile *fout){
//   fout->cd();
//
//
// }

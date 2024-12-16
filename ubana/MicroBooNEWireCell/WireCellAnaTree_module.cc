////////////////////////////////////////////////////////////////////////
// Class:       WireCellAnaTree
// Plugin Type: analyzer (art v3_01_02)
// File:        WireCellAnaTree_module.cc
//
// Generated at Fri Oct  4 15:24:42 2019 by Hanyu Wei using cetskelgen
// from cetlib version v3_05_01.
//
// 12.03.2020 modified by Wenqiang Gu (wgu@bnl.gov)
// 03.15.2021 modified by Haiwang Yu (hyu@bnl.gov)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"
#include "ubobj/WcpPort/NuSelectionKINE.h"

#include "TRandom3.h"
#include "dk2nu/tree/dk2nu.h"

#include "lardataobj/MCBase/MCShower.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"

#include <iostream>
#include <fstream>


class WireCellAnaTree;


class WireCellAnaTree : public art::EDAnalyzer {
public:
  explicit WireCellAnaTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireCellAnaTree(WireCellAnaTree const&) = delete;
  WireCellAnaTree(WireCellAnaTree&&) = delete;
  WireCellAnaTree& operator=(WireCellAnaTree const&) = delete;
  WireCellAnaTree& operator=(WireCellAnaTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void endSubRun(art::SubRun const& sr) override;

  // user defined
  void reconfigure(fhicl::ParameterSet const& pset);
  void initOutput();
  void resetOutput();
  void ShowerID(int trackId);
  void MuonID(int trackId);
  void ProtonID(int trackId);
  void save_weights(art::Event const& e);
  void save_LEEweights(art::Event const& e);
  void ReadBDTvar(nsm::NuSelectionBDT const& bdt);
  void ReadSSMBDTvar(nsm::NuSelectionBDT const& bdt);
  void ReadKINEvar(nsm::NuSelectionKINE const& kine);
  void nsbeamtiming(art::Event const& e);
  void getPMTwf(art::Event const& e, double maxP[32],double timeP[32],bool Sat[32]);
  double getBeamWF(art::Event const& e);
  std::tuple< std::vector<float>*,std::vector<float>*,std::vector<float>*,std::vector<float>* > get_extrapolated_times(art::Ptr<simb::MCParticle> particle, double mother_time);
  double get_dE_dx_range(double R, int pdg);
  void CalculateCPDF(std::vector<double> bi);
  double TimeOffset();

private:
  //debugging
  bool flag_bad;

  // Declare member data here.

  // fcl config
  std::string fContainmentLabel;
  std::string fChargeLabel;
  std::string fTruthLabel;
  std::string fMatchLabel;
  std::string fSTMLabel;
  std::string fFileType;
  std::string fWeightLabel;
  std::string fWeightLeeLabel;
  bool fMC;
  bool f_wirecellPF;
  bool fSaveWeights;
  bool fSaveLeeWeights;
  bool f_ssmBDT;
  bool f_BDTvars;
  bool f_KINEvars;
  bool fIsNuMI;
  // handling for reweighting to new NuMI flux (4.10.4) for files generated using old flux(4.9.2)
  Bool_t fNuMIOldReweight;

  bool fPFValidation; // switch of particle flow validation
  std::string fPFInputTag; // inputTag -- label:instance:process
  std::string fPFtruthInputTag; // inputTag -- label:instance:process
  float fthreshold_showerKE;
  std::vector<int> fPrimaryID;
  std::vector<int> fMuonID;
  std::vector<int> fProtonID;
  std::vector<int> fShowerID;
  std::map<int, simb::MCParticle> fParticleMap; // map from trackId to particle instance

  bool f_PFDump;
  bool f_save_track_position;
  float f_PFDump_min_truth_energy;

  bool f_savesps;

  bool f_savepmt;

  bool f_ns_time_useSSMvtx; 
  bool f_ns_time_usePID; 
  bool f_ns_time_no_photon;
  float fsol; 

  float f_shiftoffset;
  bool f_isrun3;
  float f_ccnd1_a;
  float f_ccnd1_b;
  float f_ccnd2_a;
  float f_ccnd2_b;
  float f_ccnd3_a;
  float f_ccnd3_b;
  float f_ccnd3_c;
  float f_ccnd3_d;
  float f_ccnd4_a;
  float f_ccnd4_b;

  bool f_get_redk2nu_time;
  // output
  /// PF validation
  /// when fPFValidation is true
  TTree* fPFeval;
  Int_t		f_neutrino_type;
  Float_t 	f_nuvtx_diff;
  Float_t	f_showervtx_diff;
  Float_t	f_muonvtx_diff;
  Float_t	f_reco_nuvtxX;
  Float_t	f_reco_nuvtxY;
  Float_t	f_reco_nuvtxZ;
  std::vector<float>	*f_reco_vec_showervtxX; // all showers
  std::vector<float>	*f_reco_vec_showervtxY;
  std::vector<float>	*f_reco_vec_showervtxZ;
  std::vector<float> 	*f_reco_vec_showerKE;
  Float_t	f_reco_showervtxX; // primary shower [highest energy electron & not pi0 daughters]
  Float_t	f_reco_showervtxY;
  Float_t	f_reco_showervtxZ;
  Float_t 	f_reco_showerKE;
  Float_t	f_reco_showerMomentum[4];
  Float_t	f_reco_muonvtxX; //
  Float_t	f_reco_muonvtxY;
  Float_t	f_reco_muonvtxZ;
  Float_t	f_reco_muonMomentum[4];
  Int_t		f_reco_Nproton;
  Float_t	f_reco_protonvtxX; //
  Float_t	f_reco_protonvtxY;
  Float_t	f_reco_protonvtxZ;
  Float_t	f_reco_protonMomentum[4];

  Float_t f_evtDeltaTimeNS;
  Float_t f_evtTimeNS;
  Float_t f_evtTimeNS_redk2nu;
  Float_t f_Ph_Tot;
  double  calib[32];

  Int_t		f_mcflux_run;
  Int_t		f_mcflux_evtno;
  Int_t		f_mcflux_ndecay;
  Int_t		f_mcflux_ntype;
  Int_t		f_mcflux_ptype;
  Int_t		f_mcflux_tptype;
  Float_t	f_mcflux_nuEnergy;
  Float_t	f_mcflux_vx;
  Float_t	f_mcflux_vy;
  Float_t	f_mcflux_vz;
  Float_t	f_mcflux_genx; // origin of ray from flux generator
  Float_t	f_mcflux_geny;
  Float_t	f_mcflux_genz;
  Float_t	f_mcflux_dk2gen; // distance from decay to ray origin
  Float_t	f_mcflux_gen2vtx; // distance from ray origin to event vtx

  Float_t	f_truth_corr_nuvtxX; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
  Float_t	f_truth_corr_nuvtxY;
  Float_t	f_truth_corr_nuvtxZ;
  Float_t	f_truth_corr_showervtxX; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
  Float_t	f_truth_corr_showervtxY;
  Float_t	f_truth_corr_showervtxZ;
  Float_t f_truth_showerKE;
  Float_t	f_truth_showerMomentum[4];
  Int_t   f_truth_showerPdg;
  Int_t   f_truth_showerMother;
  Float_t	f_truth_corr_muonvtxX; // primary muon, vtx = nu vtx
  Float_t	f_truth_corr_muonvtxY;
  Float_t	f_truth_corr_muonvtxZ;
  Float_t	f_truth_muonvtxX;
  Float_t	f_truth_muonvtxY;
  Float_t	f_truth_muonvtxZ;
  Float_t	f_truth_muonendX; // may not in TPC active
  Float_t	f_truth_muonendY;
  Float_t	f_truth_muonendZ;
  Float_t	f_truth_muonMomentum[4];
  Int_t		f_truth_nuIntType;
  Int_t		f_truth_nuScatType;
  Int_t   f_truth_Npi0;
  Int_t         f_truth_NprimPio;
  Float_t       f_truth_pio_energy_1;
  Float_t       f_truth_pio_energy_2;
  Float_t       f_truth_pio_angle;
  Int_t 	f_truth_NCDelta; // Radiative Delta label (both CC and NC)
  Int_t   f_truth_single_photon;
  Float_t f_truth_photon_angle;
  Float_t f_truth_photon_dis;
  Float_t	f_truth_nu_pos[4]; // X,Y,Z,T
  Float_t	f_truth_nu_momentum[4]; // Px,Py,Pz,E
  Float_t       f_redk2nu_time;
  Float_t       f_redk2nu_time_nospill;
  Float_t       f_redk2nu_deltatime;
  /// other truth info as follows save in this tree

  /// BDT input vars
  TTree* fBDT;

  //single track kdar tagger
  float ssm_flag_st_kdar;
  float ssm_Nsm;
  float ssm_Nsm_wivtx;

  float ssm_dq_dx_fwd_1;
  float ssm_dq_dx_fwd_2;
  float ssm_dq_dx_fwd_3;
  float ssm_dq_dx_fwd_4;
  float ssm_dq_dx_fwd_5;
  float ssm_dq_dx_bck_1;
  float ssm_dq_dx_bck_2;
  float ssm_dq_dx_bck_3;
  float ssm_dq_dx_bck_4;
  float ssm_dq_dx_bck_5;
  float ssm_d_dq_dx_fwd_12;
  float ssm_d_dq_dx_fwd_23;
  float ssm_d_dq_dx_fwd_34;
  float ssm_d_dq_dx_fwd_45;
  float ssm_d_dq_dx_bck_12;
  float ssm_d_dq_dx_bck_23;
  float ssm_d_dq_dx_bck_34;
  float ssm_d_dq_dx_bck_45;
  float ssm_max_dq_dx_fwd_3;
  float ssm_max_dq_dx_fwd_5;
  float ssm_max_dq_dx_bck_3;
  float ssm_max_dq_dx_bck_5;
  float ssm_max_d_dq_dx_fwd_3;
  float ssm_max_d_dq_dx_fwd_5;
  float ssm_max_d_dq_dx_bck_3;
  float ssm_max_d_dq_dx_bck_5;
  float ssm_medium_dq_dx;
  float ssm_medium_dq_dx_bp;
      //angluar info
  float ssm_angle_to_z;
  float ssm_angle_to_target;
  float ssm_angle_to_absorber;
  float ssm_angle_to_vertical;
      //directional info
  float ssm_x_dir;
  float ssm_y_dir;
  float ssm_z_dir;
      //energy info
  float ssm_kine_energy;
  float ssm_kine_energy_reduced;
      //general properties
  float ssm_vtx_activity;
  float ssm_pdg;
  float ssm_dQ_dx_cut;
  float ssm_score_mu_fwd;
  float ssm_score_p_fwd;
  float ssm_score_e_fwd;
  float ssm_score_mu_bck;
  float ssm_score_p_bck;
  float ssm_score_e_bck;
  float ssm_score_mu_fwd_bp;
  float ssm_score_p_fwd_bp;
  float ssm_score_e_fwd_bp;
      //track "straighness"
  float ssm_length;
  float ssm_direct_length;
  float ssm_length_ratio;
  float ssm_max_dev;
    //number of other particles
  float ssm_n_prim_tracks_1;
  float ssm_n_prim_tracks_3;
  float ssm_n_prim_tracks_5;
  float ssm_n_prim_tracks_8;
  float ssm_n_prim_tracks_11;
  float ssm_n_all_tracks_1;
  float ssm_n_all_tracks_3;
  float ssm_n_all_tracks_5;
  float ssm_n_all_tracks_8;
  float ssm_n_all_tracks_11;
  float ssm_n_daughter_tracks_1;
  float ssm_n_daughter_tracks_3;
  float ssm_n_daughter_tracks_5;
  float ssm_n_daughter_tracks_8;
  float ssm_n_daughter_tracks_11;
  float ssm_n_daughter_all_1;
  float ssm_n_daughter_all_3;
  float ssm_n_daughter_all_5;
  float ssm_n_daughter_all_8;
  float ssm_n_daughter_all_11;
    //properties of leading other primary track
  float ssm_prim_track1_pdg;
  float ssm_prim_track1_score_mu_fwd;
  float ssm_prim_track1_score_p_fwd;
  float ssm_prim_track1_score_e_fwd;
  float ssm_prim_track1_score_mu_bck;
  float ssm_prim_track1_score_p_bck;
  float ssm_prim_track1_score_e_bck;
  float ssm_prim_track1_length;
  float ssm_prim_track1_direct_length;
  float ssm_prim_track1_length_ratio;
  float ssm_prim_track1_max_dev;
  float ssm_prim_track1_kine_energy_range;
  float ssm_prim_track1_kine_energy_range_mu;
  float ssm_prim_track1_kine_energy_range_p;
  float ssm_prim_track1_kine_energy_range_e;
  float ssm_prim_track1_kine_energy_cal;
  float ssm_prim_track1_medium_dq_dx;
  float ssm_prim_track1_x_dir;
  float ssm_prim_track1_y_dir;
  float ssm_prim_track1_z_dir;
  float ssm_prim_track1_add_daught_track_counts_1;
  float ssm_prim_track1_add_daught_all_counts_1;
  float ssm_prim_track1_add_daught_track_counts_5;
  float ssm_prim_track1_add_daught_all_counts_5;
  float ssm_prim_track1_add_daught_track_counts_11;
  float ssm_prim_track1_add_daught_all_counts_11;
  //properties of sub-leading other primary track
  float ssm_prim_track2_pdg;
  float ssm_prim_track2_score_mu_fwd;
  float ssm_prim_track2_score_p_fwd;
  float ssm_prim_track2_score_e_fwd;
  float ssm_prim_track2_score_mu_bck;
  float ssm_prim_track2_score_p_bck;
  float ssm_prim_track2_score_e_bck;
  float ssm_prim_track2_length;
  float ssm_prim_track2_direct_length;
  float ssm_prim_track2_length_ratio;
  float ssm_prim_track2_max_dev;
  float ssm_prim_track2_kine_energy_range;
  float ssm_prim_track2_kine_energy_range_mu;
  float ssm_prim_track2_kine_energy_range_p;
  float ssm_prim_track2_kine_energy_range_e;
  float ssm_prim_track2_kine_energy_cal;
  float ssm_prim_track2_medium_dq_dx;
  float ssm_prim_track2_x_dir;
  float ssm_prim_track2_y_dir;
  float ssm_prim_track2_z_dir;
  float ssm_prim_track2_add_daught_track_counts_1;
  float ssm_prim_track2_add_daught_all_counts_1;
  float ssm_prim_track2_add_daught_track_counts_5;
  float ssm_prim_track2_add_daught_all_counts_5;
  float ssm_prim_track2_add_daught_track_counts_11;
  float ssm_prim_track2_add_daught_all_counts_11;
    //properties of leading daughter track
  float ssm_daught_track1_pdg;
  float ssm_daught_track1_score_mu_fwd;
  float ssm_daught_track1_score_p_fwd;
  float ssm_daught_track1_score_e_fwd;
  float ssm_daught_track1_score_mu_bck;
  float ssm_daught_track1_score_p_bck;
  float ssm_daught_track1_score_e_bck;
  float ssm_daught_track1_length;
  float ssm_daught_track1_direct_length;
  float ssm_daught_track1_length_ratio;
  float ssm_daught_track1_max_dev;
  float ssm_daught_track1_kine_energy_range;
  float ssm_daught_track1_kine_energy_range_mu;
  float ssm_daught_track1_kine_energy_range_p;
  float ssm_daught_track1_kine_energy_range_e;
  float ssm_daught_track1_kine_energy_cal;
  float ssm_daught_track1_medium_dq_dx;
  float ssm_daught_track1_x_dir;
  float ssm_daught_track1_y_dir;
  float ssm_daught_track1_z_dir;
  float ssm_daught_track1_add_daught_track_counts_1;
  float ssm_daught_track1_add_daught_all_counts_1;
  float ssm_daught_track1_add_daught_track_counts_5;
  float ssm_daught_track1_add_daught_all_counts_5;
  float ssm_daught_track1_add_daught_track_counts_11;
  float ssm_daught_track1_add_daught_all_counts_11;
  //properties of sub-leading daughter track
  float ssm_daught_track2_pdg;
  float ssm_daught_track2_score_mu_fwd;
  float ssm_daught_track2_score_p_fwd;
  float ssm_daught_track2_score_e_fwd;
  float ssm_daught_track2_score_mu_bck;
  float ssm_daught_track2_score_p_bck;
  float ssm_daught_track2_score_e_bck;
  float ssm_daught_track2_length;
  float ssm_daught_track2_direct_length;
  float ssm_daught_track2_length_ratio;
  float ssm_daught_track2_max_dev;
  float ssm_daught_track2_kine_energy_range;
  float ssm_daught_track2_kine_energy_range_mu;
  float ssm_daught_track2_kine_energy_range_p;
  float ssm_daught_track2_kine_energy_range_e;
  float ssm_daught_track2_kine_energy_cal;
  float ssm_daught_track2_medium_dq_dx;
  float ssm_daught_track2_x_dir;
  float ssm_daught_track2_y_dir;
  float ssm_daught_track2_z_dir;
  float ssm_daught_track2_add_daught_track_counts_1;
  float ssm_daught_track2_add_daught_all_counts_1;
  float ssm_daught_track2_add_daught_track_counts_5;
  float ssm_daught_track2_add_daught_all_counts_5;
  float ssm_daught_track2_add_daught_track_counts_11;
  float ssm_daught_track2_add_daught_all_counts_11;
    //properties of leading other primary shower
  float ssm_prim_shw1_pdg;
  float ssm_prim_shw1_score_mu_fwd;
  float ssm_prim_shw1_score_p_fwd;
  float ssm_prim_shw1_score_e_fwd;
  float ssm_prim_shw1_score_mu_bck;
  float ssm_prim_shw1_score_p_bck;
  float ssm_prim_shw1_score_e_bck;
  float ssm_prim_shw1_length;
  float ssm_prim_shw1_direct_length;
  float ssm_prim_shw1_length_ratio;
  float ssm_prim_shw1_max_dev;
  float ssm_prim_shw1_kine_energy_range;
  float ssm_prim_shw1_kine_energy_range_mu;
  float ssm_prim_shw1_kine_energy_range_p;
  float ssm_prim_shw1_kine_energy_range_e;
  float ssm_prim_shw1_kine_energy_cal;
  float ssm_prim_shw1_kine_energy_best;
  float ssm_prim_shw1_medium_dq_dx;
  float ssm_prim_shw1_x_dir;
  float ssm_prim_shw1_y_dir;
  float ssm_prim_shw1_z_dir;
  float ssm_prim_shw1_add_daught_track_counts_1;
  float ssm_prim_shw1_add_daught_all_counts_1;
  float ssm_prim_shw1_add_daught_track_counts_5;
  float ssm_prim_shw1_add_daught_all_counts_5;
  float ssm_prim_shw1_add_daught_track_counts_11;
  float ssm_prim_shw1_add_daught_all_counts_11;
  //properties of sub-leading other primary shower
  float ssm_prim_shw2_pdg;
  float ssm_prim_shw2_score_mu_fwd;
  float ssm_prim_shw2_score_p_fwd;
  float ssm_prim_shw2_score_e_fwd;
  float ssm_prim_shw2_score_mu_bck;
  float ssm_prim_shw2_score_p_bck;
  float ssm_prim_shw2_score_e_bck;
  float ssm_prim_shw2_length;
  float ssm_prim_shw2_direct_length;
  float ssm_prim_shw2_length_ratio;
  float ssm_prim_shw2_max_dev;
  float ssm_prim_shw2_kine_energy_range;
  float ssm_prim_shw2_kine_energy_range_mu;
  float ssm_prim_shw2_kine_energy_range_p;
  float ssm_prim_shw2_kine_energy_range_e;
  float ssm_prim_shw2_kine_energy_cal;
  float ssm_prim_shw2_kine_energy_best;
  float ssm_prim_shw2_medium_dq_dx;
  float ssm_prim_shw2_x_dir;
  float ssm_prim_shw2_y_dir;
  float ssm_prim_shw2_z_dir;
  float ssm_prim_shw2_add_daught_track_counts_1;
  float ssm_prim_shw2_add_daught_all_counts_1;
  float ssm_prim_shw2_add_daught_track_counts_5;
  float ssm_prim_shw2_add_daught_all_counts_5;
  float ssm_prim_shw2_add_daught_track_counts_11;
  float ssm_prim_shw2_add_daught_all_counts_11;
    //properties of leading daughter shower
  float ssm_daught_shw1_pdg;
  float ssm_daught_shw1_score_mu_fwd;
  float ssm_daught_shw1_score_p_fwd;
  float ssm_daught_shw1_score_e_fwd;
  float ssm_daught_shw1_score_mu_bck;
  float ssm_daught_shw1_score_p_bck;
  float ssm_daught_shw1_score_e_bck;
  float ssm_daught_shw1_length;
  float ssm_daught_shw1_direct_length;
  float ssm_daught_shw1_length_ratio;
  float ssm_daught_shw1_max_dev;
  float ssm_daught_shw1_kine_energy_range;
  float ssm_daught_shw1_kine_energy_range_mu;
  float ssm_daught_shw1_kine_energy_range_p;
  float ssm_daught_shw1_kine_energy_range_e;
  float ssm_daught_shw1_kine_energy_cal;
  float ssm_daught_shw1_kine_energy_best;
  float ssm_daught_shw1_medium_dq_dx;
  float ssm_daught_shw1_x_dir;
  float ssm_daught_shw1_y_dir;
  float ssm_daught_shw1_z_dir;
  float ssm_daught_shw1_add_daught_track_counts_1;
  float ssm_daught_shw1_add_daught_all_counts_1;
  float ssm_daught_shw1_add_daught_track_counts_5;
  float ssm_daught_shw1_add_daught_all_counts_5;
  float ssm_daught_shw1_add_daught_track_counts_11;
  float ssm_daught_shw1_add_daught_all_counts_11;
    //properties of sub-leading daughter shower
  float ssm_daught_shw2_pdg;
  float ssm_daught_shw2_score_mu_fwd;
  float ssm_daught_shw2_score_p_fwd;
  float ssm_daught_shw2_score_e_fwd;
  float ssm_daught_shw2_score_mu_bck;
  float ssm_daught_shw2_score_p_bck;
  float ssm_daught_shw2_score_e_bck;
  float ssm_daught_shw2_length;
  float ssm_daught_shw2_direct_length;
  float ssm_daught_shw2_length_ratio;
  float ssm_daught_shw2_max_dev;
  float ssm_daught_shw2_kine_energy_range;
  float ssm_daught_shw2_kine_energy_range_mu;
  float ssm_daught_shw2_kine_energy_range_p;
  float ssm_daught_shw2_kine_energy_range_e;
  float ssm_daught_shw2_kine_energy_cal;
  float ssm_daught_shw2_kine_energy_best;
  float ssm_daught_shw2_medium_dq_dx;
  float ssm_daught_shw2_x_dir;
  float ssm_daught_shw2_y_dir;
  float ssm_daught_shw2_z_dir;
  float ssm_daught_shw2_add_daught_track_counts_1;
  float ssm_daught_shw2_add_daught_all_counts_1;
  float ssm_daught_shw2_add_daught_track_counts_5;
  float ssm_daught_shw2_add_daught_all_counts_5;
  float ssm_daught_shw2_add_daught_track_counts_11;
  float ssm_daught_shw2_add_daught_all_counts_11;
    //event level properties
  float ssm_nu_angle_z;
  float ssm_nu_angle_target;
  float ssm_nu_angle_absorber;
  float ssm_nu_angle_vertical;
  float ssm_con_nu_angle_z;
  float ssm_con_nu_angle_target;
  float ssm_con_nu_angle_absorber;
  float ssm_con_nu_angle_vertical;
  float ssm_prim_nu_angle_z;
  float ssm_prim_nu_angle_target;
  float ssm_prim_nu_angle_absorber;
  float ssm_prim_nu_angle_vertical;
  float ssm_track_angle_z;
  float ssm_track_angle_target;
  float ssm_track_angle_absorber;
  float ssm_track_angle_vertical;
  float ssm_vtxX;
  float ssm_vtxY;
  float ssm_vtxZ;
    //off vertex stuff
  float ssm_offvtx_length;
  float ssm_offvtx_energy;
  float ssm_n_offvtx_tracks_1;
  float ssm_n_offvtx_tracks_3;
  float ssm_n_offvtx_tracks_5;
  float ssm_n_offvtx_tracks_8;
  float ssm_n_offvtx_tracks_11;
  float ssm_n_offvtx_showers_1;
  float ssm_n_offvtx_showers_3;
  float ssm_n_offvtx_showers_5;
  float ssm_n_offvtx_showers_8;
  float ssm_n_offvtx_showers_11;
    //properties of leading off vertex track
  float ssm_offvtx_track1_pdg;
  float ssm_offvtx_track1_score_mu_fwd;
  float ssm_offvtx_track1_score_p_fwd;
  float ssm_offvtx_track1_score_e_fwd;
  float ssm_offvtx_track1_score_mu_bck;
  float ssm_offvtx_track1_score_p_bck;
  float ssm_offvtx_track1_score_e_bck;
  float ssm_offvtx_track1_length;
  float ssm_offvtx_track1_direct_length;
  float ssm_offvtx_track1_max_dev;
  float ssm_offvtx_track1_kine_energy_range;
  float ssm_offvtx_track1_kine_energy_range_mu;
  float ssm_offvtx_track1_kine_energy_range_p;
  float ssm_offvtx_track1_kine_energy_range_e;
  float ssm_offvtx_track1_kine_energy_cal;
  float ssm_offvtx_track1_medium_dq_dx;
  float ssm_offvtx_track1_x_dir;
  float ssm_offvtx_track1_y_dir;
  float ssm_offvtx_track1_z_dir;
  float ssm_offvtx_track1_dist_mainvtx;
    //properties of leading off vertex shower
  float ssm_offvtx_shw1_pdg_offvtx;
  float ssm_offvtx_shw1_score_mu_fwd;
  float ssm_offvtx_shw1_score_p_fwd;
  float ssm_offvtx_shw1_score_e_fwd;
  float ssm_offvtx_shw1_score_mu_bck;
  float ssm_offvtx_shw1_score_p_bck;
  float ssm_offvtx_shw1_score_e_bck;
  float ssm_offvtx_shw1_length;
  float ssm_offvtx_shw1_direct_length;
  float ssm_offvtx_shw1_max_dev;
  float ssm_offvtx_shw1_kine_energy_best;
  float ssm_offvtx_shw1_kine_energy_range;
  float ssm_offvtx_shw1_kine_energy_range_mu;
  float ssm_offvtx_shw1_kine_energy_range_p;
  float ssm_offvtx_shw1_kine_energy_range_e;
  float ssm_offvtx_shw1_kine_energy_cal;
  float ssm_offvtx_shw1_medium_dq_dx;
  float ssm_offvtx_shw1_x_dir;
  float ssm_offvtx_shw1_y_dir;
  float ssm_offvtx_shw1_z_dir;
  float ssm_offvtx_shw1_dist_mainvtx;
    //Spacepoints
  int ssmsp_Ntrack;
  std::vector<int> *ssmsp_Nsp= new std::vector<int>;
  int ssmsp_Nsp_tot;
  std::vector<int> *ssmsp_pdg= new std::vector<int>;
  std::vector<int> *ssmsp_id= new std::vector<int>;
  std::vector<int> *ssmsp_mother= new std::vector<int>;
  std::vector<float> *ssmsp_x= new std::vector<float>;
  std::vector<float> *ssmsp_y= new std::vector<float>;
  std::vector<float> *ssmsp_z= new std::vector<float>;
  std::vector<float> *ssmsp_dx= new std::vector<float>;
  std::vector<float> *ssmsp_dQ= new std::vector<float>;
  std::vector<float> *ssmsp_KE= new std::vector<float>;
  std::vector<float> *ssmsp_containing_shower_id= new std::vector<float>;
  std::vector<float> *ssmsp_containing_shower_ke= new std::vector<float>;
  std::vector<float> *ssmsp_containing_shower_flag= new std::vector<float>;
    // KINE and other BDT variables saved seperatly for KDAR 
  float ssm_kine_reco_Enu;
  float ssm_kine_reco_add_energy;
  std::vector<float> *ssm_kine_energy_particle = new std::vector<float>;
  std::vector<int> *ssm_kine_energy_info = new std::vector<int>;
  std::vector<int> *ssm_kine_particle_type = new std::vector<int>;
  std::vector<int> *ssm_kine_energy_included = new std::vector<int>;
  float ssm_kine_pio_mass;
  int   ssm_kine_pio_flag;
  float ssm_kine_pio_vtx_dis;
  float ssm_kine_pio_energy_1;
  float ssm_kine_pio_theta_1;
  float ssm_kine_pio_phi_1;
  float ssm_kine_pio_dis_1;
  float ssm_kine_pio_energy_2;
  float ssm_kine_pio_theta_2;
  float ssm_kine_pio_phi_2;
  float ssm_kine_pio_dis_2;
  float ssm_kine_pio_angle;
  float ssm_numu_cc_flag;
  float ssm_cosmict_flag_1; // fiducial volume vertex
  float ssm_cosmict_flag_2;  // single muon
  float ssm_cosmict_flag_3;  // single muon (long)
  float ssm_cosmict_flag_4;  // kinematics muon
  float ssm_cosmict_flag_5; // kinematics muon (long)
  float ssm_cosmict_flag_6; // special ...
  float ssm_cosmict_flag_7;  // muon+ michel
  float ssm_cosmict_flag_8;  // muon + michel + special
  float ssm_cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
  std::vector<float> *ssm_cosmict_flag_10 = new std::vector<float>;  // front upstream (dirt)
  float ssm_cosmict_flag;

  //single photon vars
  float shw_sp_flag;
  float shw_sp_filled;
  float shw_sp_num_mip_tracks;
  float shw_sp_num_muons;
  float shw_sp_num_pions;
  float shw_sp_num_protons;
  float shw_sp_proton_length_1;
  float shw_sp_proton_dqdx_1;
  float shw_sp_proton_energy_1;
  float shw_sp_proton_length_2;
  float shw_sp_proton_dqdx_2;
  float shw_sp_proton_energy_2;
  float shw_sp_n_good_showers;
  float shw_sp_n_20mev_showers;
  float shw_sp_n_br1_showers;
  float shw_sp_n_br2_showers;
  float shw_sp_n_br3_showers;
  float shw_sp_n_br4_showers;
  float shw_sp_n_20br1_showers;
  std::vector<int> *shw_sp_20mev_showers;
  std::vector<int> *shw_sp_br1_showers;
  std::vector<int> *shw_sp_br2_showers;
  std::vector<int> *shw_sp_br3_showers;
  std::vector<int> *shw_sp_br4_showers;
  float shw_sp_shw_vtx_dis;
  float shw_sp_max_shw_dis;
  float shw_sp_energy;
  float shw_sp_vec_dQ_dx_0;
  float shw_sp_vec_dQ_dx_1;
  float shw_sp_max_dQ_dx_sample;
  float shw_sp_n_below_threshold;
  float shw_sp_n_below_zero;
  float shw_sp_n_lowest;
  float shw_sp_n_highest;
  float shw_sp_lowest_dQ_dx;
  float shw_sp_highest_dQ_dx;
  float shw_sp_medium_dQ_dx;
  float shw_sp_stem_length;
  float shw_sp_length_main;
  float shw_sp_length_total;
  float shw_sp_angle_beam;
  float shw_sp_iso_angle;
  float shw_sp_n_vertex;
  float shw_sp_n_good_tracks;
  float shw_sp_E_indirect_max_energy;
  float shw_sp_flag_all_above;
  float shw_sp_min_dQ_dx_5;
  float shw_sp_n_other_vertex;
  float shw_sp_n_stem_size;
  float shw_sp_flag_stem_trajectory;
  float shw_sp_min_dis;

  // extra
  float shw_sp_vec_dQ_dx_2;
  float shw_sp_vec_dQ_dx_3;
  float shw_sp_vec_dQ_dx_4;
  float shw_sp_vec_dQ_dx_5;
  float shw_sp_vec_dQ_dx_6;
  float shw_sp_vec_dQ_dx_7;
  float shw_sp_vec_dQ_dx_8;
  float shw_sp_vec_dQ_dx_9;
  float shw_sp_vec_dQ_dx_10;
  float shw_sp_vec_dQ_dx_11;
  float shw_sp_vec_dQ_dx_12;
  float shw_sp_vec_dQ_dx_13;
  float shw_sp_vec_dQ_dx_14;
  float shw_sp_vec_dQ_dx_15;
  float shw_sp_vec_dQ_dx_16;
  float shw_sp_vec_dQ_dx_17;
  float shw_sp_vec_dQ_dx_18;
  float shw_sp_vec_dQ_dx_19;
  float shw_sp_vec_median_dedx;
  float shw_sp_vec_mean_dedx;

  // photon shower pi0 identification
  float shw_sp_pio_flag;
  float shw_sp_pio_mip_id;
  float shw_sp_pio_filled;
  float shw_sp_pio_flag_pio;
  // first part of tagger ...
  float shw_sp_pio_1_flag;
  float shw_sp_pio_1_mass;
  float shw_sp_pio_1_pio_type;
  float shw_sp_pio_1_energy_1;
  float shw_sp_pio_1_energy_2;
  float shw_sp_pio_1_dis_1;
  float shw_sp_pio_1_dis_2;
  // second part of tagger
  std::vector<float> *shw_sp_pio_2_v_dis2;
  std::vector<float> *shw_sp_pio_2_v_angle2;
  std::vector<float> *shw_sp_pio_2_v_acc_length;
  std::vector<float> *shw_sp_pio_2_v_flag;

  // single photon low-energy michel
  float shw_sp_lem_shower_total_length;
  float shw_sp_lem_shower_main_length;
  float shw_sp_lem_n_3seg;
  float shw_sp_lem_e_charge;
  float shw_sp_lem_e_dQdx;
  float shw_sp_lem_shower_num_segs;
  float shw_sp_lem_shower_num_main_segs;
  float shw_sp_lem_flag;

  // bad reconstruction_1
  float shw_sp_br_filled;
  float shw_sp_br1_flag;

  //bad reconstruction 1_1
  float shw_sp_br1_1_flag;
  float shw_sp_br1_1_shower_type;
  float shw_sp_br1_1_vtx_n_segs;
  float shw_sp_br1_1_energy;
  float shw_sp_br1_1_n_segs;
  float shw_sp_br1_1_flag_sg_topology;
  float shw_sp_br1_1_flag_sg_trajectory;
  float shw_sp_br1_1_sg_length;

  // bad reconstruction 1_2
  float shw_sp_br1_2_flag;
  float shw_sp_br1_2_energy;
  float shw_sp_br1_2_n_connected;
  float shw_sp_br1_2_max_length;
  float shw_sp_br1_2_n_connected_1;
  float shw_sp_br1_2_vtx_n_segs;
  float shw_sp_br1_2_n_shower_segs;
  float shw_sp_br1_2_max_length_ratio;
  float shw_sp_br1_2_shower_length;

  // bad_reconstruction 1_3
  float shw_sp_br1_3_flag;
  float shw_sp_br1_3_energy;
  float shw_sp_br1_3_n_connected_p;
  float shw_sp_br1_3_max_length_p;
  float shw_sp_br1_3_n_shower_segs;
  float shw_sp_br1_3_flag_sg_topology;
  float shw_sp_br1_3_flag_sg_trajectory;
  float shw_sp_br1_3_n_shower_main_segs;
  float shw_sp_br1_3_sg_length;


  // bad reconstruction 2
  float shw_sp_br2_flag;
  float shw_sp_br2_flag_single_shower;
  float shw_sp_br2_num_valid_tracks;
  float shw_sp_br2_energy;
  float shw_sp_br2_angle1;
  float shw_sp_br2_angle2;
  float shw_sp_br2_angle;
  float shw_sp_br2_angle3;
  float shw_sp_br2_n_shower_main_segs;
  float shw_sp_br2_max_angle;
  float shw_sp_br2_sg_length;
  float shw_sp_br2_flag_sg_trajectory;


   //bad reconstruction 3
  float shw_sp_br3_1_energy;
  float shw_sp_br3_1_n_shower_segments;
  float shw_sp_br3_1_sg_flag_trajectory;
  float shw_sp_br3_1_sg_direct_length;
  float shw_sp_br3_1_sg_length;
  float shw_sp_br3_1_total_main_length;
  float shw_sp_br3_1_total_length;
  float shw_sp_br3_1_iso_angle;
  float shw_sp_br3_1_sg_flag_topology;
  float shw_sp_br3_1_flag;

  float shw_sp_br3_2_n_ele;
  float shw_sp_br3_2_n_other;
  float shw_sp_br3_2_energy;
  float shw_sp_br3_2_total_main_length;
  float shw_sp_br3_2_total_length;
  float shw_sp_br3_2_other_fid;
  float shw_sp_br3_2_flag;

  std::vector<float> *shw_sp_br3_3_v_energy;
  std::vector<float> *shw_sp_br3_3_v_angle;
  std::vector<float> *shw_sp_br3_3_v_dir_length;
  std::vector<float> *shw_sp_br3_3_v_length;
  std::vector<float> *shw_sp_br3_3_v_flag;

  float shw_sp_br3_4_acc_length;
  float shw_sp_br3_4_total_length;
  float shw_sp_br3_4_energy;
  float shw_sp_br3_4_flag;

  std::vector<float> *shw_sp_br3_5_v_dir_length;
  std::vector<float> *shw_sp_br3_5_v_total_length;
  std::vector<float> *shw_sp_br3_5_v_flag_avoid_muon_check;
  std::vector<float> *shw_sp_br3_5_v_n_seg;
  std::vector<float> *shw_sp_br3_5_v_angle;
  std::vector<float> *shw_sp_br3_5_v_sg_length;
  std::vector<float> *shw_sp_br3_5_v_energy;
  std::vector<float> *shw_sp_br3_5_v_n_main_segs;
  std::vector<float> *shw_sp_br3_5_v_n_segs;
  std::vector<float> *shw_sp_br3_5_v_shower_main_length;
  std::vector<float> *shw_sp_br3_5_v_shower_total_length;
  std::vector<float> *shw_sp_br3_5_v_flag;

  std::vector<float> *shw_sp_br3_6_v_angle;
  std::vector<float> *shw_sp_br3_6_v_angle1;
  std::vector<float> *shw_sp_br3_6_v_flag_shower_trajectory;
  std::vector<float> *shw_sp_br3_6_v_direct_length;
  std::vector<float> *shw_sp_br3_6_v_length;
  std::vector<float> *shw_sp_br3_6_v_n_other_vtx_segs;
  std::vector<float> *shw_sp_br3_6_v_energy;
  std::vector<float> *shw_sp_br3_6_v_flag;

  float shw_sp_br3_7_energy;
  float shw_sp_br3_7_min_angle;
  float shw_sp_br3_7_sg_length;
  float shw_sp_br3_7_shower_main_length;
  float shw_sp_br3_7_flag;

  float shw_sp_br3_8_max_dQ_dx;
  float shw_sp_br3_8_energy;
  float shw_sp_br3_8_n_main_segs;
  float shw_sp_br3_8_shower_main_length;
  float shw_sp_br3_8_shower_length;
  float shw_sp_br3_8_flag;

  float shw_sp_br3_flag;

  // BR 4
  float shw_sp_br4_1_shower_main_length;
  float shw_sp_br4_1_shower_total_length;
  float shw_sp_br4_1_min_dis;
  float shw_sp_br4_1_energy;
  float shw_sp_br4_1_flag_avoid_muon_check;
  float shw_sp_br4_1_n_vtx_segs;
  float shw_sp_br4_1_n_main_segs;
  float shw_sp_br4_1_flag;

  float shw_sp_br4_2_ratio_45;
  float shw_sp_br4_2_ratio_35;
  float shw_sp_br4_2_ratio_25;
  float shw_sp_br4_2_ratio_15;
  float shw_sp_br4_2_energy;
  float shw_sp_br4_2_ratio1_45;
  float shw_sp_br4_2_ratio1_35;
  float shw_sp_br4_2_ratio1_25;
  float shw_sp_br4_2_ratio1_15;
  float shw_sp_br4_2_iso_angle;
  float shw_sp_br4_2_iso_angle1;
  float shw_sp_br4_2_angle;
  float shw_sp_br4_2_flag;

  float shw_sp_br4_flag;

  // high energy overlap
  float shw_sp_hol_1_n_valid_tracks;
  float shw_sp_hol_1_min_angle;
  float shw_sp_hol_1_energy;
  float shw_sp_hol_1_flag_all_shower;
  float shw_sp_hol_1_min_length;
  float shw_sp_hol_1_flag;

  float shw_sp_hol_2_min_angle;
  float shw_sp_hol_2_medium_dQ_dx;
  float shw_sp_hol_2_ncount;
  float shw_sp_hol_2_energy;
  float shw_sp_hol_2_flag;

  float shw_sp_hol_flag;

  // low-energy overlap
  float shw_sp_lol_flag;

  std::vector<float> *shw_sp_lol_1_v_energy;
  std::vector<float> *shw_sp_lol_1_v_vtx_n_segs;
  std::vector<float> *shw_sp_lol_1_v_nseg;
  std::vector<float> *shw_sp_lol_1_v_angle;
  std::vector<float> *shw_sp_lol_1_v_flag;

  std::vector<float> *shw_sp_lol_2_v_length;
  std::vector<float> *shw_sp_lol_2_v_angle;
  std::vector<float> *shw_sp_lol_2_v_type;
  std::vector<float> *shw_sp_lol_2_v_vtx_n_segs;
  std::vector<float> *shw_sp_lol_2_v_energy;
  std::vector<float> *shw_sp_lol_2_v_shower_main_length;
  std::vector<float> *shw_sp_lol_2_v_flag_dir_weak;
  std::vector<float> *shw_sp_lol_2_v_flag;

  float shw_sp_lol_3_angle_beam;
  float shw_sp_lol_3_n_valid_tracks;
  float shw_sp_lol_3_min_angle;
  float shw_sp_lol_3_vtx_n_segs;
  float shw_sp_lol_3_energy;
  float shw_sp_lol_3_shower_main_length;
  float shw_sp_lol_3_n_out;
  float shw_sp_lol_3_n_sum;
  float shw_sp_lol_3_flag;
  //end single photon vars

  float cosmic_filled;
  float cosmic_flag;
  float cosmic_n_solid_tracks;
  float cosmic_energy_main_showers;
  float cosmic_energy_direct_showers;
  float cosmic_energy_indirect_showers;
  float cosmic_n_direct_showers;
  float cosmic_n_indirect_showers;
  float cosmic_n_main_showers;

  float gap_filled;
  float gap_flag;
  float gap_flag_prolong_u;
  float gap_flag_prolong_v;
  float gap_flag_prolong_w;
  float gap_flag_parallel;
  float gap_n_points;
  float gap_n_bad;
  float gap_energy;
  float gap_num_valid_tracks;
  float gap_flag_single_shower;

  float mip_quality_filled;
  float mip_quality_flag;
  float mip_quality_energy;
  float mip_quality_overlap;
  float mip_quality_n_showers;
  float mip_quality_n_tracks;
  float mip_quality_flag_inside_pi0;
  float mip_quality_n_pi0_showers;
  float mip_quality_shortest_length;
  float mip_quality_acc_length;
  float mip_quality_shortest_angle;
  float mip_quality_flag_proton;

  float mip_filled;
  float mip_flag;
  float mip_energy;
  float mip_n_end_reduction;
  float mip_n_first_mip;
  float mip_n_first_non_mip;
  float mip_n_first_non_mip_1;
  float mip_n_first_non_mip_2;
  float mip_vec_dQ_dx_0;
  float mip_vec_dQ_dx_1;
  float mip_max_dQ_dx_sample;
  float mip_n_below_threshold;
  float mip_n_below_zero;
  float mip_n_lowest;
  float mip_n_highest;
  float mip_lowest_dQ_dx;
  float mip_highest_dQ_dx;
  float mip_medium_dQ_dx;
  float mip_stem_length;
  float mip_length_main;
  float mip_length_total;
  float mip_angle_beam;
  float mip_iso_angle;
  float mip_n_vertex;
  float mip_n_good_tracks;
  float mip_E_indirect_max_energy;
  float mip_flag_all_above;
  float mip_min_dQ_dx_5;
  float mip_n_other_vertex;
  float mip_n_stem_size;
  float mip_flag_stem_trajectory;
  float mip_min_dis;

  float mip_vec_dQ_dx_2;
  float mip_vec_dQ_dx_3;
  float mip_vec_dQ_dx_4;
  float mip_vec_dQ_dx_5;
  float mip_vec_dQ_dx_6;
  float mip_vec_dQ_dx_7;
  float mip_vec_dQ_dx_8;
  float mip_vec_dQ_dx_9;
  float mip_vec_dQ_dx_10;
  float mip_vec_dQ_dx_11;
  float mip_vec_dQ_dx_12;
  float mip_vec_dQ_dx_13;
  float mip_vec_dQ_dx_14;
  float mip_vec_dQ_dx_15;
  float mip_vec_dQ_dx_16;
  float mip_vec_dQ_dx_17;
  float mip_vec_dQ_dx_18;
  float mip_vec_dQ_dx_19;

  float pio_filled;
  float pio_flag;
  float pio_mip_id;
  float pio_flag_pio;
  float pio_1_flag;
  float pio_1_mass;
  float pio_1_pio_type;
  float pio_1_energy_1;
  float pio_1_energy_2;
  float pio_1_dis_1;
  float pio_1_dis_2;
  std::vector<float> *pio_2_v_flag = new std::vector<float>;
  std::vector<float> *pio_2_v_dis2 = new std::vector<float>;
  std::vector<float> *pio_2_v_angle2 = new std::vector<float>;
  std::vector<float> *pio_2_v_acc_length = new std::vector<float>;

  float sig_flag;
  std::vector<float> *sig_1_v_flag= new std::vector<float>;
  std::vector<float> *sig_1_v_angle = new std::vector<float>;
  std::vector<float> *sig_1_v_flag_single_shower= new std::vector<float>;
  std::vector<float> *sig_1_v_energy= new std::vector<float>;
  std::vector<float> *sig_1_v_energy_1= new std::vector<float>;
  std::vector<float> *sig_2_v_flag= new std::vector<float>;
  std::vector<float> *sig_2_v_energy= new std::vector<float>;
  std::vector<float> *sig_2_v_shower_angle= new std::vector<float>;
  std::vector<float> *sig_2_v_flag_single_shower= new std::vector<float>;
  std::vector<float> *sig_2_v_medium_dQ_dx= new std::vector<float>;
  std::vector<float> *sig_2_v_start_dQ_dx= new std::vector<float>;

  float mgo_flag;
  float mgo_energy;
  float mgo_max_energy;
  float mgo_total_energy;
  float mgo_n_showers;
  float mgo_max_energy_1;
  float mgo_max_energy_2;
  float mgo_total_other_energy;
  float mgo_n_total_showers;
  float mgo_total_other_energy_1;

  float mgt_flag;
  float mgt_flag_single_shower;
  float mgt_max_energy;
  float mgt_energy;
  float mgt_total_other_energy;
  float mgt_max_energy_1;
  float mgt_e_indirect_max_energy;
  float mgt_e_direct_max_energy;
  float mgt_n_direct_showers;
  float mgt_e_direct_total_energy;
  float mgt_flag_indirect_max_pio;
  float mgt_e_indirect_total_energy;

  float stw_flag;
  float stw_1_flag;
  float stw_1_energy;
  float stw_1_dis;
  float stw_1_dQ_dx;
  float stw_1_flag_single_shower;
  float stw_1_n_pi0;
  float stw_1_num_valid_tracks;
  std::vector<float> *stw_2_v_flag = new std::vector<float>;
  std::vector<float> *stw_2_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *stw_2_v_energy = new std::vector<float>;
  std::vector<float> *stw_2_v_angle = new std::vector<float>;
  std::vector<float> *stw_2_v_dir_length = new std::vector<float>;
  std::vector<float> *stw_2_v_max_dQ_dx = new std::vector<float>;
  std::vector<float> *stw_3_v_flag = new std::vector<float>;
  std::vector<float> *stw_3_v_angle = new std::vector<float>;
  std::vector<float> *stw_3_v_dir_length = new std::vector<float>;
  std::vector<float> *stw_3_v_energy = new std::vector<float>;
  std::vector<float> *stw_3_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *stw_4_v_flag = new std::vector<float>;
  std::vector<float> *stw_4_v_angle = new std::vector<float>;
  std::vector<float> *stw_4_v_dis = new std::vector<float>;
  std::vector<float> *stw_4_v_energy = new std::vector<float>;

  float spt_flag;
  float spt_flag_single_shower;
  float spt_energy;
  float spt_shower_main_length;
  float spt_shower_total_length;
  float spt_angle_beam;
  float spt_angle_vertical;
  float spt_max_dQ_dx;
  float spt_angle_beam_1;
  float spt_angle_drift;
  float spt_angle_drift_1;
  float spt_num_valid_tracks;
  float spt_n_vtx_segs;
  float spt_max_length;

  float stem_len_flag;
  float stem_len_energy;
  float stem_len_length;
  float stem_len_flag_avoid_muon_check;
  float stem_len_num_daughters;
  float stem_len_daughter_length;

  float lem_flag;
  float lem_shower_total_length;
  float lem_shower_main_length;
  float lem_n_3seg;
  float lem_e_charge;
  float lem_e_dQdx;
  float lem_shower_num_segs;
  float lem_shower_num_main_segs;

  float brm_flag;
  float brm_n_mu_segs;
  float brm_Ep;
  float brm_energy;
  float brm_acc_length;
  float brm_shower_total_length;
  float brm_connected_length;
  float brm_n_size;
  float brm_acc_direct_length;
  float brm_n_shower_main_segs;
  float brm_n_mu_main;

  float cme_flag;
  float cme_mu_energy;
  float cme_energy;
  float cme_mu_length;
  float cme_length;
  float cme_angle_beam;

  float anc_flag;
  float anc_energy;
  float anc_angle;
  float anc_max_angle;
  float anc_max_length;
  float anc_acc_forward_length;
  float anc_acc_backward_length;
  float anc_acc_forward_length1;
  float anc_shower_main_length;
  float anc_shower_total_length;
  float anc_flag_main_outside;

  float stem_dir_filled;
  float stem_dir_flag;
  float stem_dir_flag_single_shower;
  float stem_dir_angle;
  float stem_dir_energy;
  float stem_dir_angle1;
  float stem_dir_angle2;
  float stem_dir_angle3;
  float stem_dir_ratio;

  float vis_flag;
  float vis_1_filled;
  float vis_1_flag;
  float vis_1_n_vtx_segs;
  float vis_1_energy;
  float vis_1_num_good_tracks;
  float vis_1_max_angle;
  float vis_1_max_shower_angle;
  float vis_1_tmp_length1;
  float vis_1_tmp_length2;
  float vis_1_particle_type;
  float vis_2_filled;
  float vis_2_flag;
  float vis_2_n_vtx_segs;
  float vis_2_min_angle;
  float vis_2_min_weak_track;
  float vis_2_angle_beam;
  float vis_2_min_angle1;
  float vis_2_iso_angle1;
  float vis_2_min_medium_dQ_dx;
  float vis_2_min_length;
  float vis_2_sg_length;
  float vis_2_max_angle;
  float vis_2_max_weak_track;

  float br_filled;
  float br1_flag;
  float br1_1_flag;
  float br1_1_shower_type;
  float br1_1_vtx_n_segs;
  float br1_1_energy;
  float br1_1_n_segs;
  float br1_1_flag_sg_topology;
  float br1_1_flag_sg_trajectory;
  float br1_1_sg_length;
  float br1_2_flag;
  float br1_2_energy;
  float br1_2_n_connected;
  float br1_2_max_length;
  float br1_2_n_connected_1;
  float br1_2_vtx_n_segs;
  float br1_2_n_shower_segs;
  float br1_2_max_length_ratio;
  float br1_2_shower_length;
  float br1_3_flag;
  float br1_3_energy;
  float br1_3_n_connected_p;
  float br1_3_max_length_p;
  float br1_3_n_shower_segs;
  float br1_3_flag_sg_topology;
  float br1_3_flag_sg_trajectory;
  float br1_3_n_shower_main_segs;
  float br1_3_sg_length;

  float br2_flag;
  float br2_flag_single_shower;
  float br2_num_valid_tracks;
  float br2_energy;
  float br2_angle1;
  float br2_angle2;
  float br2_angle;
  float br2_angle3;
  float br2_n_shower_main_segs;
  float br2_max_angle;
  float br2_sg_length;
  float br2_flag_sg_trajectory;

  float br3_flag;
  float br3_1_flag;
  float br3_1_energy;
  float br3_1_n_shower_segments;
  float br3_1_sg_flag_trajectory;
  float br3_1_sg_direct_length;
  float br3_1_sg_length;
  float br3_1_total_main_length;
  float br3_1_total_length;
  float br3_1_iso_angle;
  float br3_1_sg_flag_topology;
  float br3_2_flag;
  float br3_2_n_ele;
  float br3_2_n_other;
  float br3_2_energy;
  float br3_2_total_main_length;
  float br3_2_total_length;
  float br3_2_other_fid;
  std::vector<float> *br3_3_v_flag = new std::vector<float>;
  std::vector<float> *br3_3_v_energy = new std::vector<float>;
  std::vector<float> *br3_3_v_angle = new std::vector<float>;
  std::vector<float> *br3_3_v_dir_length = new std::vector<float>;
  std::vector<float> *br3_3_v_length = new std::vector<float>;
  float br3_4_flag;
  float br3_4_acc_length;
  float br3_4_total_length;
  float br3_4_energy;
  std::vector<float> *br3_5_v_flag = new std::vector<float>;
  std::vector<float> *br3_5_v_dir_length = new std::vector<float>;
  std::vector<float> *br3_5_v_total_length = new std::vector<float>;
  std::vector<float> *br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  std::vector<float> *br3_5_v_n_seg = new std::vector<float>;
  std::vector<float> *br3_5_v_angle = new std::vector<float>;
  std::vector<float> *br3_5_v_sg_length = new std::vector<float>;
  std::vector<float> *br3_5_v_energy = new std::vector<float>;
  std::vector<float> *br3_5_v_n_main_segs = new std::vector<float>;
  std::vector<float> *br3_5_v_n_segs = new std::vector<float>;
  std::vector<float> *br3_5_v_shower_main_length = new std::vector<float>;
  std::vector<float> *br3_5_v_shower_total_length = new std::vector<float>;
  std::vector<float> *br3_6_v_flag = new std::vector<float>;
  std::vector<float> *br3_6_v_angle = new std::vector<float>;
  std::vector<float> *br3_6_v_angle1 = new std::vector<float>;
  std::vector<float> *br3_6_v_flag_shower_trajectory = new std::vector<float>;
  std::vector<float> *br3_6_v_direct_length = new std::vector<float>;
  std::vector<float> *br3_6_v_length = new std::vector<float>;
  std::vector<float> *br3_6_v_n_other_vtx_segs = new std::vector<float>;
  std::vector<float> *br3_6_v_energy = new std::vector<float>;
  float br3_7_flag;
  float br3_7_energy;
  float br3_7_min_angle;
  float br3_7_sg_length;
  float br3_7_shower_main_length;
  float br3_8_flag;
  float br3_8_max_dQ_dx;
  float br3_8_energy;
  float br3_8_n_main_segs;
  float br3_8_shower_main_length;
  float br3_8_shower_length;

  float br4_flag;
  float br4_1_flag;
  float br4_1_shower_main_length;
  float br4_1_shower_total_length;
  float br4_1_min_dis;
  float br4_1_energy;
  float br4_1_flag_avoid_muon_check;
  float br4_1_n_vtx_segs;
  float br4_1_n_main_segs;
  float br4_2_flag;
  float br4_2_ratio_45;
  float br4_2_ratio_35;
  float br4_2_ratio_25;
  float br4_2_ratio_15;
  float br4_2_energy;
  float br4_2_ratio1_45;
  float br4_2_ratio1_35;
  float br4_2_ratio1_25;
  float br4_2_ratio1_15;
  float br4_2_iso_angle;
  float br4_2_iso_angle1;
  float br4_2_angle;

  float tro_flag;
  std::vector<float> *tro_1_v_flag= new std::vector<float>;
  std::vector<float> *tro_1_v_particle_type= new std::vector<float>;
  std::vector<float> *tro_1_v_flag_dir_weak= new std::vector<float>;
  std::vector<float> *tro_1_v_min_dis = new std::vector<float>;
  std::vector<float> *tro_1_v_sg1_length = new std::vector<float>;
  std::vector<float> *tro_1_v_shower_main_length = new std::vector<float>;
  std::vector<float> *tro_1_v_max_n_vtx_segs= new std::vector<float>;
  std::vector<float> *tro_1_v_tmp_length = new std::vector<float>;
  std::vector<float> *tro_1_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *tro_1_v_dQ_dx_cut = new std::vector<float>;
  std::vector<float> *tro_1_v_flag_shower_topology= new std::vector<float>;
  std::vector<float> *tro_2_v_flag= new std::vector<float>;
  std::vector<float> *tro_2_v_energy = new std::vector<float>;
  std::vector<float> *tro_2_v_stem_length = new std::vector<float>;
  std::vector<float> *tro_2_v_iso_angle = new std::vector<float>;
  std::vector<float> *tro_2_v_max_length = new std::vector<float>;
  std::vector<float> *tro_2_v_angle = new std::vector<float>;
  float tro_3_flag;
  float tro_3_stem_length;
  float tro_3_n_muon_segs;
  float tro_3_energy;
  std::vector<float> *tro_4_v_flag= new std::vector<float>;
  std::vector<float> *tro_4_v_dir2_mag = new std::vector<float>;
  std::vector<float> *tro_4_v_angle = new std::vector<float>;
  std::vector<float> *tro_4_v_angle1 = new std::vector<float>;
  std::vector<float> *tro_4_v_angle2 = new std::vector<float>;
  std::vector<float> *tro_4_v_length = new std::vector<float>;
  std::vector<float> *tro_4_v_length1 = new std::vector<float>;
  std::vector<float> *tro_4_v_medium_dQ_dx = new std::vector<float>;
  std::vector<float> *tro_4_v_end_dQ_dx = new std::vector<float>;
  std::vector<float> *tro_4_v_energy = new std::vector<float>;
  std::vector<float> *tro_4_v_shower_main_length = new std::vector<float>;
  std::vector<float> *tro_4_v_flag_shower_trajectory= new std::vector<float>;
  std::vector<float> *tro_5_v_flag = new std::vector<float>;
  std::vector<float> *tro_5_v_max_angle = new std::vector<float>;
  std::vector<float> *tro_5_v_min_angle = new std::vector<float>;
  std::vector<float> *tro_5_v_max_length = new std::vector<float>;
  std::vector<float> *tro_5_v_iso_angle = new std::vector<float>;
  std::vector<float> *tro_5_v_n_vtx_segs= new std::vector<float>;
  std::vector<float> *tro_5_v_min_count= new std::vector<float>;
  std::vector<float> *tro_5_v_max_count= new std::vector<float>;
  std::vector<float> *tro_5_v_energy = new std::vector<float>;

  float hol_flag;
  float hol_1_flag;
  float hol_1_n_valid_tracks;
  float hol_1_min_angle;
  float hol_1_energy;
  float hol_1_flag_all_shower;
  float hol_1_min_length;
  float hol_2_flag;
  float hol_2_min_angle;
  float hol_2_medium_dQ_dx;
  float hol_2_ncount;
  float hol_2_energy;

  float lol_flag;
  std::vector<float> *lol_1_v_flag= new std::vector<float>;
  std::vector<float> *lol_1_v_energy = new std::vector<float>;
  std::vector<float> *lol_1_v_vtx_n_segs= new std::vector<float>;
  std::vector<float> *lol_1_v_nseg= new std::vector<float>;
  std::vector<float> *lol_1_v_angle= new std::vector<float>;
  std::vector<float> *lol_2_v_flag = new std::vector<float>;
  std::vector<float> *lol_2_v_length= new std::vector<float>;
  std::vector<float> *lol_2_v_angle= new std::vector<float>;
  std::vector<float> *lol_2_v_type= new std::vector<float>;
  std::vector<float> *lol_2_v_vtx_n_segs= new std::vector<float>;
  std::vector<float> *lol_2_v_energy= new std::vector<float>;
  std::vector<float> *lol_2_v_shower_main_length= new std::vector<float>;
  std::vector<float> *lol_2_v_flag_dir_weak= new std::vector<float>;
  float lol_3_flag;
  float lol_3_angle_beam;
  float lol_3_n_valid_tracks;
  float lol_3_min_angle;
  float lol_3_vtx_n_segs;
  float lol_3_energy;
  float lol_3_shower_main_length;
  float lol_3_n_out;
  float lol_3_n_sum;

  float cosmict_flag_1; // fiducial volume vertex
  float cosmict_flag_2;  // single muon
  float cosmict_flag_3;  // single muon (long)
  float cosmict_flag_4;  // kinematics muon
  float cosmict_flag_5; // kinematics muon (long)
  float cosmict_flag_6; // special ...
  float cosmict_flag_7;  // muon+ michel
  float cosmict_flag_8;  // muon + michel + special
  float cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
  std::vector<float> *cosmict_flag_10 = new std::vector<float>;  // front upstream (dirt)
  float cosmict_flag;
  float cosmict_2_filled;
  float cosmict_2_particle_type;
  float cosmict_2_n_muon_tracks;
  float cosmict_2_total_shower_length;
  float cosmict_2_flag_inside;
  float cosmict_2_angle_beam;
  float cosmict_2_flag_dir_weak;
  float cosmict_2_dQ_dx_end;
  float cosmict_2_dQ_dx_front;
  float cosmict_2_theta;
  float cosmict_2_phi;
  float cosmict_2_valid_tracks;
  float cosmict_3_filled;
  float cosmict_3_flag_inside;
  float cosmict_3_angle_beam;
  float cosmict_3_flag_dir_weak;
  float cosmict_3_dQ_dx_end;
  float cosmict_3_dQ_dx_front;
  float cosmict_3_theta;
  float cosmict_3_phi;
  float cosmict_3_valid_tracks;
  float cosmict_4_filled;
  float cosmict_4_flag_inside;
  float cosmict_4_angle_beam;
  float cosmict_4_connected_showers;  // need to be careful about the nueCC ...
  float cosmict_5_filled;
  float cosmict_5_flag_inside;
  float cosmict_5_angle_beam;
  float cosmict_5_connected_showers;
  float cosmict_6_filled;
  float cosmict_6_flag_dir_weak;
  float cosmict_6_flag_inside;
  float cosmict_6_angle;
  float cosmict_7_filled;
  float cosmict_7_flag_sec;
  float cosmict_7_n_muon_tracks;
  float cosmict_7_total_shower_length;
  float cosmict_7_flag_inside;
  float cosmict_7_angle_beam;
  float cosmict_7_flag_dir_weak;
  float cosmict_7_dQ_dx_end;
  float cosmict_7_dQ_dx_front;
  float cosmict_7_theta;
  float cosmict_7_phi;
  float cosmict_8_filled;
  float cosmict_8_flag_out;
  float cosmict_8_muon_length;
  float cosmict_8_acc_length;
  std::vector<float> *cosmict_10_flag_inside= new std::vector<float>;
  std::vector<float> *cosmict_10_vtx_z= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_shower= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_dir_weak= new std::vector<float>;
  std::vector<float> *cosmict_10_angle_beam= new std::vector<float>;
  std::vector<float> *cosmict_10_length = new std::vector<float>;

  float numu_cc_flag;
  std::vector<float> *numu_cc_flag_1= new std::vector<float>;
  std::vector<float> *numu_cc_1_particle_type= new std::vector<float>;
  std::vector<float> *numu_cc_1_length= new std::vector<float>;
  std::vector<float> *numu_cc_1_medium_dQ_dx= new std::vector<float>;
  std::vector<float> *numu_cc_1_dQ_dx_cut= new std::vector<float>;
  std::vector<float> *numu_cc_1_direct_length= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_tracks= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_all= new std::vector<float>;
  std::vector<float> *numu_cc_flag_2= new std::vector<float>;
  std::vector<float> *numu_cc_2_length= new std::vector<float>;
  std::vector<float> *numu_cc_2_total_length= new std::vector<float>;
  std::vector<float> *numu_cc_2_n_daughter_tracks= new std::vector<float>;
  std::vector<float> *numu_cc_2_n_daughter_all = new std::vector<float>;
  float numu_cc_flag_3;
  float numu_cc_3_particle_type;
  float numu_cc_3_max_length;
  float numu_cc_3_acc_track_length;
  float numu_cc_3_max_length_all;
  float numu_cc_3_max_muon_length;
  float numu_cc_3_n_daughter_tracks;
  float numu_cc_3_n_daughter_all;

  float cosmict_2_4_score;
  float cosmict_3_5_score;
  float cosmict_6_score;
  float cosmict_7_score;
  float cosmict_8_score;
  float cosmict_10_score;
  float numu_1_score;
  float numu_2_score;
  float numu_3_score;
  float cosmict_score;
  float numu_score;
  float mipid_score;
  float gap_score;
  float hol_lol_score;
  float cme_anc_score;
  float mgo_mgt_score;
  float br1_score;
  float br3_score;
  float br3_3_score;
  float br3_5_score;
  float br3_6_score;
  float stemdir_br2_score;
  float trimuon_score;
  float br4_tro_score;
  float mipquality_score;
  float pio_1_score;
  float pio_2_score;
  float stw_spt_score;
  float vis_1_score;
  float vis_2_score;
  float stw_2_score;
  float stw_3_score;
  float stw_4_score;
  float sig_1_score;
  float sig_2_score;
  float lol_1_score;
  float lol_2_score;
  float tro_1_score;
  float tro_2_score;
  float tro_4_score;
  float tro_5_score;
  float nue_score;

  /// kinematic variables
  TTree* fKINE;
  float kine_reco_Enu;
  float kine_reco_add_energy;
  std::vector<float> *kine_energy_particle = new std::vector<float>;
  std::vector<int> *kine_energy_info = new std::vector<int>;
  std::vector<int> *kine_particle_type = new std::vector<int>;
  std::vector<int> *kine_energy_included = new std::vector<int>;
  float kine_pio_mass;
  int 	kine_pio_flag;
  float	kine_pio_vtx_dis;
  float kine_pio_energy_1;
  float kine_pio_theta_1;
  float kine_pio_phi_1;
  float kine_pio_dis_1;
  float kine_pio_energy_2;
  float kine_pio_theta_2;
  float kine_pio_phi_2;
  float kine_pio_dis_2;
  float kine_pio_angle;

  ///

  TTree* fTreeEval;
  Int_t           f_run;
  Int_t           f_subRun;
  Int_t           f_event;
  UInt_t          f_trigger_bits;
  Bool_t          f_flash_found;
  Int_t           f_flash_found_asInt;
  Float_t         f_flash_time;
  Float_t         f_flash_measPe;
  Float_t         f_flash_predPe;
  Bool_t          f_match_found;
  Int_t           f_match_found_asInt;
  UInt_t          f_match_type;
  Bool_t          f_match_isFC;
  Bool_t          f_match_isTgm;
  Bool_t          f_match_notFC_FV;
  Bool_t          f_match_notFC_SP;
  Bool_t          f_match_notFC_DC;
  Float_t         f_match_charge; // main flag collection plane charge
  Float_t         f_match_energy;
  Float_t         f_lm_cluster_length;
  Bool_t          f_image_fail;
  Float_t	  f_match_chargeU;
  Float_t	  f_match_chargeV;
  Float_t	  f_match_chargeY;
  Float_t	  f_match_energyY;
  Bool_t	  f_lightmismatch;

  Float_t         f_truth_nuEnergy;
  Float_t         f_truth_energyInside;
  Float_t         f_truth_electronInside;
  Int_t           f_truth_nuPdg;
  Bool_t          f_truth_isCC;
  Bool_t          f_truth_isEligible;
  Bool_t          f_NC_truth_isEligible;
  Bool_t          f_truth_isFC;
  Bool_t          f_truth_vtxInside;
  Float_t         f_truth_vtxX;
  Float_t         f_truth_vtxY;
  Float_t         f_truth_vtxZ;
  Float_t         f_truth_nuTime;
  Float_t         f_match_completeness;
  Float_t         f_match_completeness_energy;
  Float_t         f_match_purity;
  Float_t         f_match_purity_xz;
  Float_t         f_match_purity_xy;

  Float_t 	  f_weight_spline;
  Float_t	  f_weight_cv;
  Float_t	  f_weight_lee;

  Int_t		  f_stm_eventtype;
  Int_t		  f_stm_lowenergy;
  Int_t		  f_stm_LM;
  Int_t		  f_stm_TGM;
  Int_t		  f_stm_STM;
  Int_t		  f_stm_FullDead;
  Float_t	  f_stm_clusterlength;

  TTree* fTreePot;
  std::string fPOT_inputTag;
  bool fPOT_counting=false;
 // beamdata:bnbETOR875
 // EXT?
 // MC [label: generator] : [instance: ]

  //int frun;
  //int fsubRun;
  double fpot_tor875;
  double fpot_tor875good;
  double fspill_tor875;
  double fspill_tor875good;

  enum LIMITS {
      MAX_TRACKS = 30000,
  };

  int truth_Ntrack;  // number of tracks in MC
  int truth_id[MAX_TRACKS];  // track id; size == truth_Ntrack
  int truth_pdg[MAX_TRACKS];  // track particle pdg; size == truth_Ntrack
  std::vector<std::string > *truth_process;
  int truth_mother[MAX_TRACKS];  // mother id of this track; size == truth_Ntrack
  float truth_startXYZT[MAX_TRACKS][4];  // start position of this track; size == truth_Ntrack
  float truth_endXYZT[MAX_TRACKS][4];  // end position of this track; size == truth_Ntrack
  float truth_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == truth_Ntrack
  float truth_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == truth_Ntrack
  std::vector<std::vector<int> > *truth_daughters;  // daughters id of this track; vector
  TObjArray *fMC_trackPosition;

  int reco_Ntrack;  // number of tracks in MC
  int reco_id[MAX_TRACKS];  // track id; size == reco_Ntrack
  int reco_pdg[MAX_TRACKS];  // track particle pdg; size == reco_Ntrack
  std::vector<std::string > *reco_process;
  int reco_mother[MAX_TRACKS];  // mother id of this track; size == reco_Ntrack
  float reco_startXYZT[MAX_TRACKS][4];  // start position of this track; size == reco_Ntrack
  float reco_endXYZT[MAX_TRACKS][4];  // end position of this track; size == reco_Ntrack
  float reco_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == reco_Ntrack
  float reco_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == reco_Ntrack
  std::vector<std::vector<int> > *reco_daughters;  // daughters id of this track; vector

  std::vector<float> *f_sps_x;
  std::vector<float> *f_sps_y;
  std::vector<float> *f_sps_z;
  std::vector<float> *f_sps_e;
  std::vector<float> *f_sps_pdg;
  std::vector<float> *f_sps_id;

  std::vector<int> *f_PMT_ID;
  std::vector<float> *f_PMT_Time;
  std::vector<float> *f_PMT_Amp;
  std::vector<float> *f_PMT_TimeProp;
  std::vector<float> *f_PMT_TimeDP;
  std::vector<float> *f_PMT_TimeDL;
  std::vector<bool> *f_PMT_Sat;
  float f_RWM_Time;

  int mc_isnu; // is neutrino interaction
  int mc_nGeniePrimaries; // number of Genie primaries
  int mc_nu_pdg; // pdg code of neutrino
  int mc_nu_ccnc; // cc or nc
  int mc_nu_mode; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
  int mc_nu_intType; // interaction type
  int mc_nu_target; // target interaction
  int mc_hitnuc; // hit nucleon
  int mc_hitquark; // hit quark
  double mc_nu_Q2; // Q^2
  double mc_nu_W; // W
  double mc_nu_X; // X
  double mc_nu_Y; // Y
  double mc_nu_Pt; // Pt
  double mc_nu_Theta; // angle relative to lepton
  float mc_nu_pos[4];  // interaction position of nu
  float mc_nu_mom[4];  // interaction momentum of nu

  //beam timing plots
  TH1F *H_time;
  TH1F *H_maxH;
  TH1F *H_t0_Beam;
  TH2F *H_TimeVsPh;
  TH1F *ns_time;

  //for redk2nu resimulating beam structure
  double fTimeBetweenBuckets;  //time between buckets
  double fBucketTimeSigma;  //how wide is distribution in bucket
  double fNBucketsPerBatch;   
  double fNFilledBucketsPerBatch;
  //double fDisallowedBatchMask; //disallow individual batches
  double fGlobalOffset;  //always displaced by this (in ns)
  std::vector<double> fbi;
  std::vector<int>    fDisallowedBatchMask;  ///< disallow individual batches
  std::vector<double> fCummulativeBatchPDF;  ///< summed prob for batches
  TRandom* fRndmGen;
};


WireCellAnaTree::WireCellAnaTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  // fcl config
  reconfigure(p);

  // T_eval / event
  // Histograms / event
  // T_pot / subrun
  initOutput();

}

void WireCellAnaTree::reconfigure(fhicl::ParameterSet const& pset)
{
  std::cout<<"------------ WireCellAnaTree::reconfigure ----------"<<std::endl;

  fContainmentLabel = pset.get<std::string>("ContainmentLabel");
  fChargeLabel = pset.get<std::string>("ChargeLabel");
  fTruthLabel = pset.get<std::string>("TruthLabel");
  fMatchLabel = pset.get<std::string>("MatchLabel");
  fMC = pset.get<bool>("MC"); // overlay and full mc
  f_wirecellPF = pset.get<bool>("wirecellPF", false);
  f_ssmBDT = pset.get<bool>("ssmBDT", false);
  f_BDTvars = pset.get<bool>("BDTvars", false);
  f_KINEvars = pset.get<bool>("KINEvars", false);
  fSaveWeights = pset.get<bool>("SaveWeights", false); // GENIE weights
  fSaveLeeWeights = pset.get<bool>("SaveLeeWeights", false); // LEE weights
  fIsNuMI = pset.get<bool>("IsNuMI", false); // is true, convert to BNB style
  fNuMIOldReweight = pset.get<bool>("NuMIOldReweight", false); // if true, converts underlying numi flux distribution from old (4.9.2) to new (4.10.4)
  fSTMLabel = pset.get<std::string>("STMLabel");
  fFileType = pset.get<std::string>("FileType", "empty");
  fWeightLabel = pset.get<std::string>("WeightLabel", "");
  fWeightLeeLabel = pset.get<std::string>("WeightLeeLabel","");

  fPOT_inputTag = pset.get<std::string>("POT_inputTag");
  fPOT_counting = pset.get<bool>("POT_counting");

  fPFValidation = pset.get<bool>("PF_validation");
  fPFInputTag = pset.get<std::string>("PF_inputtag");
  fPFtruthInputTag = pset.get<std::string>("PFtruth_inputtag");
  fthreshold_showerKE = pset.get<float>("Threshold_showerKE"); // GeV

  f_PFDump = pset.get<bool>("PFDump", false);
  f_save_track_position = pset.get<bool>("save_track_position", false);
  f_PFDump_min_truth_energy = pset.get<float>("PFDump_min_truth_energy", 0);

  f_savesps = pset.get<bool>("SaveSPS", false);

  f_savepmt  = pset.get<bool>("SavePMT", false);

  f_ns_time_useSSMvtx = pset.get<bool>("ns_time_useSSMvtx", false);
  f_ns_time_usePID = pset.get<bool>("ns_time_usePID", false);
  f_ns_time_no_photon = pset.get<bool>("ns_time_no_photon", false);
  fsol = pset.get<float>("ns_time_sol", 0.033356);

  f_shiftoffset = pset.get<float>("ShiftOffset", 0);
  f_isrun3 = pset.get<bool>("isRun3", false);
  f_ccnd1_a = pset.get<float>("ccnd1_a", 0.529594);
  f_ccnd1_b = pset.get<float>("ccnd1_b", 7.13804);
  f_ccnd2_a = pset.get<float>("ccnd2_a", 0.068752);
  f_ccnd2_b = pset.get<float>("ccnd2_b", 2.32023);
  f_ccnd3_a = pset.get<float>("ccnd3_a", 0.4697);
  f_ccnd3_b = pset.get<float>("ccnd3_b", 0.004233);
  f_ccnd3_c = pset.get<float>("ccnd3_c", 0.000001006);
  f_ccnd3_d = pset.get<float>("ccnd3_d", -0.195);
  f_ccnd4_a = pset.get<float>("ccnd4_a", 0);
  f_ccnd4_b = pset.get<float>("ccnd4_b", 0);

  f_get_redk2nu_time = pset.get<bool>("get_redk2nu_time", false);
  fTimeBetweenBuckets     = pset.get<float>("TimeBetweenBuckets",1e9/53.103e6);
  fBucketTimeSigma        = pset.get<float>("BucketTimeSigma",0.750);
  fNBucketsPerBatch       = pset.get<int>("NBucketsPerBatch",84);
  fNFilledBucketsPerBatch = pset.get<int>("NFilledBucketsPerBatch",81);
  fGlobalOffset           = pset.get<double>("GlobalOffset",0);
  fbi = pset.get<std::vector<double>>("BatchIntensities",{1,1,1,1,1,1}); // all ones would be equal batches, set last to lower to emulate slip stacking
  if(fbi.size()==0){fbi = {1,1,1,1,1,1};}
  if(f_get_redk2nu_time){
    fRndmGen = new TRandom3(0);
    CalculateCPDF(fbi);
  }
}

void WireCellAnaTree::initOutput()
{
  std::cout<<"------------ WireCellAnaTree::initOutput ----------"<<std::endl;

  art::ServiceHandle<art::TFileService> tfs;

  if(!fMC && !fIsNuMI){
    //ns beam timing plots for validation
    H_time= tfs->make<TH1F>("H_time","Time PMT",500, 0,6000);
    H_maxH= tfs->make<TH1F>("H_maxH","Max amplitude",800,2000,2100);
    H_t0_Beam= tfs->make<TH1F>("H_t0_Beam","T_0 beam",800,3000,8450);
    H_TimeVsPh= tfs->make<TH2F>("H_TimeVsPh","H_TimeVsPh",  100, -50,50,  100, 0,500);
  }else if(fIsNuMI){
    H_time= tfs->make<TH1F>("H_time","Time PMT",2000, 0,20000);
    H_maxH= tfs->make<TH1F>("H_maxH","Max amplitude",2400,1800,2100);
    H_t0_Beam= tfs->make<TH1F>("H_t0_Beam","T_0 beam",300,0,150);
    H_TimeVsPh= tfs->make<TH2F>("H_TimeVsPh","H_TimeVsPh",  100, -50,50,  100, 0,500);
 }

  fTreeEval = tfs->make<TTree>("T_eval", "T_eval");

  fTreeEval->Branch("run", 			&f_run);
  fTreeEval->Branch("subrun", 			&f_subRun);
  fTreeEval->Branch("event", 			&f_event);
  fTreeEval->Branch("file_type", 		&fFileType);
  fTreeEval->Branch("trigger_bits", 		&f_trigger_bits);
  fTreeEval->Branch("flash_found", 		&f_flash_found);
  fTreeEval->Branch("flash_found_asInt", 	&f_flash_found_asInt);
  fTreeEval->Branch("flash_time", 		&f_flash_time);
  fTreeEval->Branch("flash_measPe", 		&f_flash_measPe);
  fTreeEval->Branch("flash_predPe", 		&f_flash_predPe);
  fTreeEval->Branch("match_found", 		&f_match_found);
  fTreeEval->Branch("match_found_asInt", 	&f_match_found_asInt);
  fTreeEval->Branch("match_type", 		&f_match_type);
  fTreeEval->Branch("match_isFC", 		&f_match_isFC);
  fTreeEval->Branch("match_isTgm", 		&f_match_isTgm);
  fTreeEval->Branch("match_notFC_FV", 		&f_match_notFC_FV);
  fTreeEval->Branch("match_notFC_SP", 		&f_match_notFC_SP);
  fTreeEval->Branch("match_notFC_DC", 		&f_match_notFC_DC);
  fTreeEval->Branch("match_chargeU", 		&f_match_chargeU);
  fTreeEval->Branch("match_chargeV", 		&f_match_chargeV);
  fTreeEval->Branch("match_chargeY", 		&f_match_chargeY);
  fTreeEval->Branch("match_energyY", 		&f_match_energyY);
  fTreeEval->Branch("lm_cluster_length",        &f_lm_cluster_length);
  fTreeEval->Branch("image_fail",               &f_image_fail);
  fTreeEval->Branch("light_mismatch", 		&f_lightmismatch);
  fTreeEval->Branch("match_charge", 		&f_match_charge);
  fTreeEval->Branch("match_energy", 		&f_match_energy);
  fTreeEval->Branch("stm_eventtype",		&f_stm_eventtype);
  fTreeEval->Branch("stm_lowenergy",		&f_stm_lowenergy);
  fTreeEval->Branch("stm_LM",			&f_stm_LM);
  fTreeEval->Branch("stm_TGM",			&f_stm_TGM);
  fTreeEval->Branch("stm_STM",			&f_stm_STM);
  fTreeEval->Branch("stm_FullDead",		&f_stm_FullDead);
  fTreeEval->Branch("stm_clusterlength",	&f_stm_clusterlength);

  if( fMC==true ){
  fTreeEval->Branch("truth_nuEnergy", 		&f_truth_nuEnergy);
  fTreeEval->Branch("truth_energyInside", 	&f_truth_energyInside);
  fTreeEval->Branch("truth_electronInside", 	&f_truth_electronInside);
  fTreeEval->Branch("truth_nuPdg", 		&f_truth_nuPdg);
  fTreeEval->Branch("truth_isCC", 		&f_truth_isCC);
  fTreeEval->Branch("truth_isEligible", 	&f_truth_isEligible);
  fTreeEval->Branch("truth_NCisEligible", 	&f_NC_truth_isEligible);
  fTreeEval->Branch("truth_isFC", 		&f_truth_isFC);
  fTreeEval->Branch("truth_vtxInside", 		&f_truth_vtxInside);
  fTreeEval->Branch("truth_vtxX", 		&f_truth_vtxX);
  fTreeEval->Branch("truth_vtxY", 		&f_truth_vtxY);
  fTreeEval->Branch("truth_vtxZ", 		&f_truth_vtxZ);
  fTreeEval->Branch("truth_nuTime", 		&f_truth_nuTime);
  fTreeEval->Branch("match_completeness", 	&f_match_completeness);
  fTreeEval->Branch("match_completeness_energy",&f_match_completeness_energy);
  fTreeEval->Branch("match_purity", 		&f_match_purity);
  fTreeEval->Branch("match_purity_xz", 		&f_match_purity_xz);
  fTreeEval->Branch("match_purity_xy", 		&f_match_purity_xy);

  fTreeEval->Branch("weight_spline", 		&f_weight_spline); //MicroBooNE GENIE tune on top of weight_CV; weight_spline*weight_cv = weight
  fTreeEval->Branch("weight_cv",		&f_weight_cv); //MicroBooNE MCC9 untuned GENIE v3
  fTreeEval->Branch("weight_lee",		&f_weight_lee); //MicroBooNE MCC9 LEE weight

  }
  fTreePot = tfs->make<TTree>("T_pot", "T_pot");
  fTreePot->Branch("runNo", &f_run);
  fTreePot->Branch("subRunNo", &f_subRun);
  fTreePot->Branch("pot_tor875", &fpot_tor875);
  fTreePot->Branch("pot_tor875good", &fpot_tor875good);
  fTreePot->Branch("spill_tor875", &fspill_tor875);
  fTreePot->Branch("spill_tor875good", &fspill_tor875good);

  /// PF validation
  fPFeval = tfs->make<TTree>("T_PFeval", "T_PFeval");
  fPFeval->Branch("run", 			&f_run);
  fPFeval->Branch("subrun", 			&f_subRun);
  fPFeval->Branch("event", 			&f_event);
  fPFeval->Branch("neutrino_type", 		&f_neutrino_type);
  fPFeval->Branch("reco_nuvtxX", 		&f_reco_nuvtxX);
  fPFeval->Branch("reco_nuvtxY", 		&f_reco_nuvtxY);
  fPFeval->Branch("reco_nuvtxZ", 		&f_reco_nuvtxZ);
  fPFeval->Branch("reco_vec_showervtxX", 		&f_reco_vec_showervtxX);
  fPFeval->Branch("reco_vec_showervtxY", 		&f_reco_vec_showervtxY);
  fPFeval->Branch("reco_vec_showervtxZ", 		&f_reco_vec_showervtxZ);
  fPFeval->Branch("reco_vec_showerKE", 		&f_reco_vec_showerKE);
  fPFeval->Branch("reco_showervtxX", 		&f_reco_showervtxX);
  fPFeval->Branch("reco_showervtxY", 		&f_reco_showervtxY);
  fPFeval->Branch("reco_showervtxZ", 		&f_reco_showervtxZ);
  fPFeval->Branch("reco_showerKE", 		&f_reco_showerKE);
  fPFeval->Branch("reco_showerMomentum", 	&f_reco_showerMomentum, "reco_showerMomentum[4]/F");
  fPFeval->Branch("reco_muonvtxX", 		&f_reco_muonvtxX);
  fPFeval->Branch("reco_muonvtxY", 		&f_reco_muonvtxY);
  fPFeval->Branch("reco_muonvtxZ", 		&f_reco_muonvtxZ);
  fPFeval->Branch("reco_muonMomentum", 		&f_reco_muonMomentum, "reco_muonMomentum[4]/F");
  fPFeval->Branch("reco_Nproton", 		&f_reco_Nproton);
  fPFeval->Branch("reco_protonvtxX", 		&f_reco_protonvtxX);
  fPFeval->Branch("reco_protonvtxY", 		&f_reco_protonvtxY);
  fPFeval->Branch("reco_protonvtxZ", 		&f_reco_protonvtxZ);
  fPFeval->Branch("reco_protonMomentum", 	&f_reco_protonMomentum, "reco_protonMomentum[4]/F");

  if(!fMC || fIsNuMI) {
    fPFeval->Branch("evtDeltaTimeNS", &f_evtDeltaTimeNS);
    fPFeval->Branch("evtTimeNS", &f_evtTimeNS);
    fPFeval->Branch("Ph_Tot", &f_Ph_Tot);
    if(f_get_redk2nu_time){
      fPFeval->Branch("evtTimeNS_redk2nu", &f_evtTimeNS_redk2nu);
    }
  }

  if( fMC==true ){
  fPFeval->Branch("nuvtx_diff",			&f_nuvtx_diff);
  fPFeval->Branch("showervtx_diff",		&f_showervtx_diff);
  fPFeval->Branch("muonvtx_diff",		&f_muonvtx_diff);
  fPFeval->Branch("mcflux_run", 		&f_mcflux_run);
  fPFeval->Branch("mcflux_evtno", 		&f_mcflux_evtno);
  fPFeval->Branch("mcflux_ndecay", 		&f_mcflux_ndecay);
  fPFeval->Branch("mcflux_ntype", 		&f_mcflux_ntype);
  fPFeval->Branch("mcflux_ptype", 		&f_mcflux_ptype);
  fPFeval->Branch("mcflux_tptype", 		&f_mcflux_tptype);
  fPFeval->Branch("mcflux_nuEnergy", 		&f_mcflux_nuEnergy);
  fPFeval->Branch("mcflux_vx", 			&f_mcflux_vx);
  fPFeval->Branch("mcflux_vy", 			&f_mcflux_vy);
  fPFeval->Branch("mcflux_vz", 			&f_mcflux_vz);
  fPFeval->Branch("mcflux_genx", 		&f_mcflux_genx);
  fPFeval->Branch("mcflux_geny", 		&f_mcflux_geny);
  fPFeval->Branch("mcflux_genz", 		&f_mcflux_genz);
  fPFeval->Branch("mcflux_dk2gen", 		&f_mcflux_dk2gen);
  fPFeval->Branch("mcflux_gen2vtx", 		&f_mcflux_gen2vtx);

  fPFeval->Branch("truth_corr_nuvtxX", 		&f_truth_corr_nuvtxX);
  fPFeval->Branch("truth_corr_nuvtxY", 		&f_truth_corr_nuvtxY);
  fPFeval->Branch("truth_corr_nuvtxZ", 		&f_truth_corr_nuvtxZ);
  fPFeval->Branch("truth_corr_showervtxX", 	&f_truth_corr_showervtxX);
  fPFeval->Branch("truth_corr_showervtxY", 	&f_truth_corr_showervtxY);
  fPFeval->Branch("truth_corr_showervtxZ", 	&f_truth_corr_showervtxZ);
  fPFeval->Branch("truth_showerKE", 		&f_truth_showerKE);
  fPFeval->Branch("truth_showerMomentum", 	&f_truth_showerMomentum, "truth_showerMomentum[4]/F");
  fPFeval->Branch("truth_showerPdg", 		&f_truth_showerPdg);
  fPFeval->Branch("truth_showerMother", 		&f_truth_showerMother);
  fPFeval->Branch("truth_corr_muonvtxX",	&f_truth_corr_muonvtxX);
  fPFeval->Branch("truth_corr_muonvtxY", 	&f_truth_corr_muonvtxY);
  fPFeval->Branch("truth_corr_muonvtxZ", 	&f_truth_corr_muonvtxZ);
  fPFeval->Branch("truth_muonvtxX",		&f_truth_muonvtxX);
  fPFeval->Branch("truth_muonvtxY", 		&f_truth_muonvtxY);
  fPFeval->Branch("truth_muonvtxZ", 		&f_truth_muonvtxZ);
  fPFeval->Branch("truth_muonendX",		&f_truth_muonendX);
  fPFeval->Branch("truth_muonendY", 		&f_truth_muonendY);
  fPFeval->Branch("truth_muonendZ", 		&f_truth_muonendZ);
  fPFeval->Branch("truth_muonMomentum", 	&f_truth_muonMomentum, "truth_muonMomentum[4]/F");
  fPFeval->Branch("truth_nuEnergy", 		&f_truth_nuEnergy);
  fPFeval->Branch("truth_energyInside", 	&f_truth_energyInside);
  fPFeval->Branch("truth_electronInside", 	&f_truth_electronInside);
  fPFeval->Branch("truth_nuPdg", 		&f_truth_nuPdg);
  fPFeval->Branch("truth_isCC", 		&f_truth_isCC);
  fPFeval->Branch("truth_vtxX", 		&f_truth_vtxX);
  fPFeval->Branch("truth_vtxY", 		&f_truth_vtxY);
  fPFeval->Branch("truth_vtxZ", 		&f_truth_vtxZ);
  fPFeval->Branch("truth_nuTime", 		&f_truth_nuTime);
  fPFeval->Branch("truth_nuIntType",		&f_truth_nuIntType);
  fPFeval->Branch("truth_nuScatType",		&f_truth_nuScatType);
  fPFeval->Branch("truth_Npi0",		&f_truth_Npi0);
  fPFeval->Branch("truth_NprimPio",		&f_truth_NprimPio);
  fPFeval->Branch("truth_pio_energy_1",         &f_truth_pio_energy_1);
  fPFeval->Branch("truth_pio_energy_2",         &f_truth_pio_energy_2);
  fPFeval->Branch("truth_pio_angle",            &f_truth_pio_angle);
  fPFeval->Branch("truth_NCDelta",            	&f_truth_NCDelta);
  fPFeval->Branch("truth_single_photon",        &f_truth_single_photon);
  fPFeval->Branch("truth_photon_angle", &f_truth_photon_angle);
  fPFeval->Branch("truth_photon_dis", &f_truth_photon_dis);
  fPFeval->Branch("truth_nu_pos",         	&f_truth_nu_pos, "truth_nu_pos[4]/F");
  fPFeval->Branch("truth_nu_momentum",         	&f_truth_nu_momentum, "truth_nu_momentum[4]/F");
    if(f_get_redk2nu_time){
      fPFeval->Branch("redk2nu_time",          &f_redk2nu_time);
      fPFeval->Branch("redk2nu_time_nospill",          &f_redk2nu_time_nospill);
      fPFeval->Branch("redk2nu_deltatime",          &f_redk2nu_deltatime);
    }
  }

  // PFDump
  if(f_PFDump) {

      if( fMC==true ){
        fPFeval->Branch("truth_Ntrack", &truth_Ntrack);
        fPFeval->Branch("truth_id", &truth_id, "truth_id[truth_Ntrack]/I");
        fPFeval->Branch("truth_pdg", &truth_pdg, "truth_pdg[truth_Ntrack]/I");
        fPFeval->Branch("truth_process", &truth_process);
        fPFeval->Branch("truth_mother", &truth_mother, "truth_mother[truth_Ntrack]/I");
        fPFeval->Branch("truth_startXYZT", &truth_startXYZT, "truth_startXYZT[truth_Ntrack][4]/F");
        fPFeval->Branch("truth_endXYZT", &truth_endXYZT, "truth_endXYZT[truth_Ntrack][4]/F");
        fPFeval->Branch("truth_startMomentum", &truth_startMomentum, "truth_startMomentum[truth_Ntrack][4]/F");
        fPFeval->Branch("truth_endMomentum", &truth_endMomentum, "truth_endMomentum[truth_Ntrack][4]/F");
        fPFeval->Branch("truth_daughters", &truth_daughters);
        fMC_trackPosition = new TObjArray();
        fMC_trackPosition->SetOwner(kTRUE);
        fPFeval->Branch("fMC_trackPosition", &fMC_trackPosition);

        fPFeval->Branch("mc_isnu", &mc_isnu);
        fPFeval->Branch("mc_nGeniePrimaries", &mc_nGeniePrimaries);
        fPFeval->Branch("mc_nu_pdg", &mc_nu_pdg);
        fPFeval->Branch("mc_nu_ccnc", &mc_nu_ccnc);
        fPFeval->Branch("mc_nu_mode", &mc_nu_mode);
        fPFeval->Branch("mc_nu_intType", &mc_nu_intType);
        fPFeval->Branch("mc_nu_target", &mc_nu_target);
        fPFeval->Branch("mc_hitnuc", &mc_hitnuc);
        fPFeval->Branch("mc_hitquark", &mc_hitquark);
        fPFeval->Branch("mc_nu_Q2", &mc_nu_Q2);
        fPFeval->Branch("mc_nu_W", &mc_nu_W);
        fPFeval->Branch("mc_nu_X", &mc_nu_X);
        fPFeval->Branch("mc_nu_Y", &mc_nu_Y);
        fPFeval->Branch("mc_nu_Pt", &mc_nu_Pt);
        fPFeval->Branch("mc_nu_Theta", &mc_nu_Theta);
        fPFeval->Branch("mc_nu_pos", &mc_nu_pos, "mc_nu_pos[4]/F");
        fPFeval->Branch("mc_nu_mom", &mc_nu_mom, "mc_nu_mom[4]/F");
      }

      fPFeval->Branch("reco_Ntrack", &reco_Ntrack);
      fPFeval->Branch("reco_id", &reco_id, "reco_id[reco_Ntrack]/I");
      fPFeval->Branch("reco_pdg", &reco_pdg, "reco_pdg[reco_Ntrack]/I");
      fPFeval->Branch("reco_process", &reco_process);
      fPFeval->Branch("reco_mother", &reco_mother, "reco_mother[reco_Ntrack]/I");
      fPFeval->Branch("reco_startXYZT", &reco_startXYZT, "reco_startXYZT[reco_Ntrack][4]/F");
      fPFeval->Branch("reco_endXYZT", &reco_endXYZT, "reco_endXYZT[reco_Ntrack][4]/F");
      fPFeval->Branch("reco_startMomentum", &reco_startMomentum, "reco_startMomentum[reco_Ntrack][4]/F");
      fPFeval->Branch("reco_endMomentum", &reco_endMomentum, "reco_endMomentum[reco_Ntrack][4]/F");
      fPFeval->Branch("reco_daughters", &reco_daughters);

  }

  if (f_savesps){
    fPFeval->Branch("reco_sps_x", &f_sps_x);
    fPFeval->Branch("reco_sps_y", &f_sps_y);
    fPFeval->Branch("reco_sps_z", &f_sps_z);
    fPFeval->Branch("reco_sps_e", &f_sps_e);
    fPFeval->Branch("reco_sps_pdg", &f_sps_pdg);
    fPFeval->Branch("reco_sps_id", &f_sps_id);
  }

  if (f_savepmt){
    fPFeval->Branch("PMT_ID", &f_PMT_ID);
    fPFeval->Branch("PMT_Time", &f_PMT_Time);
    fPFeval->Branch("PMT_Amp", &f_PMT_Amp);
    fPFeval->Branch("PMT_TimeProp", &f_PMT_TimeProp);
    fPFeval->Branch("PMT_TimeDP", &f_PMT_TimeDP);
    fPFeval->Branch("PMT_TimeDL", &f_PMT_TimeDL);
    fPFeval->Branch("PMT_Sat", &f_PMT_Sat);
    fPFeval->Branch("RWM_Time", &f_RWM_Time);

  }

  fBDT = tfs->make<TTree>("T_BDTvars", "T_BDTvars");
  if(f_BDTvars){
  fBDT->Branch("run",				&f_run);
  fBDT->Branch("subrun", 			&f_subRun);
  fBDT->Branch("event", 			&f_event);
  fBDT->Branch("nuvtx_diff",			&f_nuvtx_diff);
  fBDT->Branch("showervtx_diff",		&f_showervtx_diff);
  fBDT->Branch("muonvtx_diff",			&f_muonvtx_diff);

  if(f_ssmBDT){
    //single track kdar tagger
    fBDT->Branch("ssm_flag_st_kdar",&ssm_flag_st_kdar,"ssm_flag_st_kdar/F");
    fBDT->Branch("ssm_Nsm",&ssm_Nsm,"ssm_Nsm/F");
    fBDT->Branch("ssm_Nsm_wivtx",&ssm_Nsm_wivtx,"ssm_Nsm_wivtx/F");

    //only filled if there is one ssm
    //properties of the ssm
    //dq/dx info
      fBDT->Branch("ssm_dq_dx_fwd_1", &ssm_dq_dx_fwd_1, "ssm_dq_dx_fwd_1/F");
      fBDT->Branch("ssm_dq_dx_fwd_2", &ssm_dq_dx_fwd_2, "ssm_dq_dx_fwd_2/F");
      fBDT->Branch("ssm_dq_dx_fwd_3", &ssm_dq_dx_fwd_3, "ssm_dq_dx_fwd_3/F");
      fBDT->Branch("ssm_dq_dx_fwd_4", &ssm_dq_dx_fwd_4, "ssm_dq_dx_fwd_4/F");
      fBDT->Branch("ssm_dq_dx_fwd_5", &ssm_dq_dx_fwd_5, "ssm_dq_dx_fwd_5/F");
      fBDT->Branch("ssm_dq_dx_bck_1", &ssm_dq_dx_bck_1, "ssm_dq_dx_bck_1/F");
      fBDT->Branch("ssm_dq_dx_bck_2", &ssm_dq_dx_bck_2, "ssm_dq_dx_bck_2/F");
      fBDT->Branch("ssm_dq_dx_bck_3", &ssm_dq_dx_bck_3, "ssm_dq_dx_bck_3/F");
      fBDT->Branch("ssm_dq_dx_bck_4", &ssm_dq_dx_bck_4, "ssm_dq_dx_bck_4/F");
      fBDT->Branch("ssm_dq_dx_bck_5", &ssm_dq_dx_bck_5, "ssm_dq_dx_bck_5/F");
      fBDT->Branch("ssm_d_dq_dx_fwd_12", &ssm_d_dq_dx_fwd_12, "ssm_d_dq_dx_fwd_12/F");
      fBDT->Branch("ssm_d_dq_dx_fwd_23", &ssm_d_dq_dx_fwd_23, "ssm_d_dq_dx_fwd_23/F");
      fBDT->Branch("ssm_d_dq_dx_fwd_34", &ssm_d_dq_dx_fwd_34, "ssm_d_dq_dx_fwd_34/F");
      fBDT->Branch("ssm_d_dq_dx_fwd_45", &ssm_d_dq_dx_fwd_45, "ssm_d_dq_dx_fwd_45/F");
      fBDT->Branch("ssm_d_dq_dx_bck_12", &ssm_d_dq_dx_bck_12, "ssm_d_dq_dx_bck_12/F");
      fBDT->Branch("ssm_d_dq_dx_bck_23", &ssm_d_dq_dx_bck_23, "ssm_d_dq_dx_bck_23/F");
      fBDT->Branch("ssm_d_dq_dx_bck_34", &ssm_d_dq_dx_bck_34, "ssm_d_dq_dx_bck_34/F");
      fBDT->Branch("ssm_d_dq_dx_bck_45", &ssm_d_dq_dx_bck_45, "ssm_d_dq_dx_bck_45/F");
      fBDT->Branch("ssm_max_dq_dx_fwd_3", &ssm_max_dq_dx_fwd_3, "ssm_max_dq_dx_fwd_3/F");
      fBDT->Branch("ssm_max_dq_dx_fwd_5", &ssm_max_dq_dx_fwd_5, "ssm_max_dq_dx_fwd_5/F");
      fBDT->Branch("ssm_max_dq_dx_bck_3", &ssm_max_dq_dx_bck_3, "ssm_max_dq_dx_bck_3/F");
      fBDT->Branch("ssm_max_dq_dx_bck_5", &ssm_max_dq_dx_bck_5, "ssm_max_dq_dx_bck_5/F");
      fBDT->Branch("ssm_max_d_dq_dx_fwd_3", &ssm_max_d_dq_dx_fwd_3, "ssm_max_d_dq_dx_fwd_3/F");
      fBDT->Branch("ssm_max_d_dq_dx_fwd_5", &ssm_max_d_dq_dx_fwd_5, "ssm_max_d_dq_dx_fwd_5/F");
      fBDT->Branch("ssm_max_d_dq_dx_bck_3", &ssm_max_d_dq_dx_bck_3, "ssm_max_d_dq_dx_bck_3/F");
      fBDT->Branch("ssm_max_d_dq_dx_bck_5", &ssm_max_d_dq_dx_bck_5, "ssm_max_d_dq_dx_bck_5/F");
      fBDT->Branch("ssm_medium_dq_dx", &ssm_medium_dq_dx, "ssm_medium_dq_dx/F");
      fBDT->Branch("ssm_medium_dq_dx_bp", &ssm_medium_dq_dx_bp, "ssm_medium_dq_dx_bp/F");
    //angluar info
      fBDT->Branch("ssm_angle_to_z", &ssm_angle_to_z, "ssm_angle_to_z/F");
      fBDT->Branch("ssm_angle_to_target", &ssm_angle_to_target, "ssm_angle_to_target/F");
      fBDT->Branch("ssm_angle_to_absorber", &ssm_angle_to_absorber, "ssm_angle_to_absorber/F");
      fBDT->Branch("ssm_angle_to_vertical", &ssm_angle_to_vertical, "ssm_angle_to_vertical/F");
    //directional info
      fBDT->Branch("ssm_x_dir", &ssm_x_dir, "ssm_x_dir/F");
      fBDT->Branch("ssm_y_dir", &ssm_y_dir, "ssm_y_dir/F");
      fBDT->Branch("ssm_z_dir", &ssm_z_dir, "ssm_z_dir/F");
    //energy info
      fBDT->Branch("ssm_kine_energy", &ssm_kine_energy, "ssm_kine_energy/F");
      fBDT->Branch("ssm_kine_energy_reduced", &ssm_kine_energy_reduced, "ssm_kine_energy_reduced/F");
    //general properties
      fBDT->Branch("ssm_vtx_activity", &ssm_vtx_activity, "ssm_vtx_activity/F");
      fBDT->Branch("ssm_pdg", &ssm_pdg, "ssm_pdg/F");
      fBDT->Branch("ssm_dQ_dx_cut", &ssm_dQ_dx_cut, "ssm_dQ_dx_cut/F");
      fBDT->Branch("ssm_score_mu_fwd", &ssm_score_mu_fwd, "ssm_score_mu_fwd/F");
      fBDT->Branch("ssm_score_p_fwd", &ssm_score_p_fwd, "ssm_score_p_fwd/F");
      fBDT->Branch("ssm_score_e_fwd", &ssm_score_e_fwd, "ssm_score_e_fwd/F");
      fBDT->Branch("ssm_score_mu_bck", &ssm_score_mu_bck, "ssm_score_mu_bck/F");
      fBDT->Branch("ssm_score_p_bck", &ssm_score_p_bck, "ssm_score_p_bck/F");
      fBDT->Branch("ssm_score_e_bck", &ssm_score_e_bck, "ssm_score_e_bck/F");
      fBDT->Branch("ssm_score_mu_fwd_bp", &ssm_score_mu_fwd_bp, "ssm_score_mu_fwd_bp/F");
      fBDT->Branch("ssm_score_p_fwd_bp", &ssm_score_p_fwd_bp, "ssm_score_p_fwd_bp/F");
      fBDT->Branch("ssm_score_e_fwd_bp", &ssm_score_e_fwd_bp, "ssm_score_e_fwd_bp/F");
    //track "straighness"
      fBDT->Branch("ssm_length", &ssm_length, "ssm_length/F");
      fBDT->Branch("ssm_direct_length", &ssm_direct_length, "ssm_direct_length/F");
      fBDT->Branch("ssm_length_ratio", &ssm_length_ratio, "ssm_length_ratio/F");
      fBDT->Branch("ssm_max_dev", &ssm_max_dev, "ssm_max_dev/F");
    //number of other particles
      fBDT->Branch("ssm_n_prim_tracks_1", &ssm_n_prim_tracks_1, "ssm_n_prim_tracks_1/F");
      fBDT->Branch("ssm_n_prim_tracks_3", &ssm_n_prim_tracks_3, "ssm_n_prim_tracks_3/F");
      fBDT->Branch("ssm_n_prim_tracks_5", &ssm_n_prim_tracks_5, "ssm_n_prim_tracks_5/F");
      fBDT->Branch("ssm_n_prim_tracks_8", &ssm_n_prim_tracks_8, "ssm_n_prim_tracks_8/F");
      fBDT->Branch("ssm_n_prim_tracks_11", &ssm_n_prim_tracks_11, "ssm_n_prim_tracks_11/F");
      fBDT->Branch("ssm_n_all_tracks_1", &ssm_n_all_tracks_1, "ssm_n_all_tracks_1/F");
      fBDT->Branch("ssm_n_all_tracks_3", &ssm_n_all_tracks_3, "ssm_n_all_tracks_3/F");
      fBDT->Branch("ssm_n_all_tracks_5", &ssm_n_all_tracks_5, "ssm_n_all_tracks_5/F");
      fBDT->Branch("ssm_n_all_tracks_8", &ssm_n_all_tracks_8, "ssm_n_all_tracks_8/F");
      fBDT->Branch("ssm_n_all_tracks_11", &ssm_n_all_tracks_11, "ssm_n_all_tracks_11/F");
      fBDT->Branch("ssm_n_daughter_tracks_1", &ssm_n_daughter_tracks_1, "ssm_n_daughter_tracks_1/F");
      fBDT->Branch("ssm_n_daughter_tracks_3", &ssm_n_daughter_tracks_3, "ssm_n_daughter_tracks_3/F");
      fBDT->Branch("ssm_n_daughter_tracks_5", &ssm_n_daughter_tracks_5, "ssm_n_daughter_tracks_5/F");
      fBDT->Branch("ssm_n_daughter_tracks_8", &ssm_n_daughter_tracks_8, "ssm_n_daughter_tracks_8/F");
      fBDT->Branch("ssm_n_daughter_tracks_11", &ssm_n_daughter_tracks_11, "ssm_n_daughter_tracks_11/F");
      fBDT->Branch("ssm_n_daughter_all_1", &ssm_n_daughter_all_1, "ssm_n_daughter_all_1/F");
      fBDT->Branch("ssm_n_daughter_all_3", &ssm_n_daughter_all_3, "ssm_n_daughter_all_3/F");
      fBDT->Branch("ssm_n_daughter_all_5", &ssm_n_daughter_all_5, "ssm_n_daughter_all_5/F");
      fBDT->Branch("ssm_n_daughter_all_8", &ssm_n_daughter_all_8, "ssm_n_daughter_all_8/F");
      fBDT->Branch("ssm_n_daughter_all_11", &ssm_n_daughter_all_11, "ssm_n_daughter_all_11/F");
    //properties of leading other primary track
      fBDT->Branch("ssm_prim_track1_pdg", &ssm_prim_track1_pdg, "ssm_prim_track1_pdg/F");
      fBDT->Branch("ssm_prim_track1_score_mu_fwd", &ssm_prim_track1_score_mu_fwd, "ssm_prim_track1_score_mu_fwd/F");
      fBDT->Branch("ssm_prim_track1_score_p_fwd", &ssm_prim_track1_score_p_fwd, "ssm_prim_track1_score_p_fwd/F");
      fBDT->Branch("ssm_prim_track1_score_e_fwd", &ssm_prim_track1_score_e_fwd, "ssm_prim_track1_score_e_fwd/F");
      fBDT->Branch("ssm_prim_track1_score_mu_bck", &ssm_prim_track1_score_mu_bck, "ssm_prim_track1_score_mu_bck/F");
      fBDT->Branch("ssm_prim_track1_score_p_bck", &ssm_prim_track1_score_p_bck, "ssm_prim_track1_score_p_bck/F");
      fBDT->Branch("ssm_prim_track1_score_e_bck", &ssm_prim_track1_score_e_bck, "ssm_prim_track1_score_e_bck/F");
      fBDT->Branch("ssm_prim_track1_length", &ssm_prim_track1_length, "ssm_prim_track1_length/F");
      fBDT->Branch("ssm_prim_track1_direct_length", &ssm_prim_track1_direct_length, "ssm_prim_track1_direct_length/F");
      fBDT->Branch("ssm_prim_track1_length_ratio", &ssm_prim_track1_length_ratio, "ssm_prim_track1_length_ratio/F");
      fBDT->Branch("ssm_prim_track1_max_dev", &ssm_prim_track1_max_dev, "ssm_prim_track1_max_dev/F");
      fBDT->Branch("ssm_prim_track1_kine_energy_range", &ssm_prim_track1_kine_energy_range, "ssm_prim_track1_kine_energy_range/F");
      fBDT->Branch("ssm_prim_track1_kine_energy_range_mu", &ssm_prim_track1_kine_energy_range_mu, "ssm_prim_track1_kine_energy_range_mu/F");
      fBDT->Branch("ssm_prim_track1_kine_energy_range_p", &ssm_prim_track1_kine_energy_range_p, "ssm_prim_track1_kine_energy_range_p/F");
      fBDT->Branch("ssm_prim_track1_kine_energy_range_e", &ssm_prim_track1_kine_energy_range_e, "ssm_prim_track1_kine_energy_range_e/F");
      fBDT->Branch("ssm_prim_track1_kine_energy_cal", &ssm_prim_track1_kine_energy_cal, "ssm_prim_track1_kine_energy_cal/F");
      fBDT->Branch("ssm_prim_track1_medium_dq_dx", &ssm_prim_track1_medium_dq_dx, "ssm_prim_track1_medium_dq_dx/F");
      fBDT->Branch("ssm_prim_track1_x_dir", &ssm_prim_track1_x_dir, "ssm_prim_track1_x_dir/F");
      fBDT->Branch("ssm_prim_track1_y_dir", &ssm_prim_track1_y_dir, "ssm_prim_track1_y_dir/F");
      fBDT->Branch("ssm_prim_track1_z_dir", &ssm_prim_track1_z_dir, "ssm_prim_track1_z_dir/F");
      fBDT->Branch("ssm_prim_track1_add_daught_track_counts_1", &ssm_prim_track1_add_daught_track_counts_1, "ssm_prim_track1_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_prim_track1_add_daught_all_counts_1", &ssm_prim_track1_add_daught_all_counts_1, "ssm_prim_track1_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_prim_track1_add_daught_track_counts_5", &ssm_prim_track1_add_daught_track_counts_5, "ssm_prim_track1_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_prim_track1_add_daught_all_counts_5", &ssm_prim_track1_add_daught_all_counts_5, "ssm_prim_track1_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_prim_track1_add_daught_track_counts_11", &ssm_prim_track1_add_daught_track_counts_11, "ssm_prim_track1_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_prim_track1_add_daught_all_counts_11", &ssm_prim_track1_add_daught_all_counts_11, "ssm_prim_track1_add_daught_all_counts_11/F");
    //properties of sub-leading other primary track
      fBDT->Branch("ssm_prim_track2_pdg", &ssm_prim_track2_pdg, "ssm_prim_track2_pdg/F");
      fBDT->Branch("ssm_prim_track2_score_mu_fwd", &ssm_prim_track2_score_mu_fwd, "ssm_prim_track2_score_mu_fwd/F");
      fBDT->Branch("ssm_prim_track2_score_p_fwd", &ssm_prim_track2_score_p_fwd, "ssm_prim_track2_score_p_fwd/F");
      fBDT->Branch("ssm_prim_track2_score_e_fwd", &ssm_prim_track2_score_e_fwd, "ssm_prim_track2_score_e_fwd/F");
      fBDT->Branch("ssm_prim_track2_score_mu_bck", &ssm_prim_track2_score_mu_bck, "ssm_prim_track2_score_mu_bck/F");
      fBDT->Branch("ssm_prim_track2_score_p_bck", &ssm_prim_track2_score_p_bck, "ssm_prim_track2_score_p_bck/F");
      fBDT->Branch("ssm_prim_track2_score_e_bck", &ssm_prim_track2_score_e_bck, "ssm_prim_track2_score_e_bck/F");
      fBDT->Branch("ssm_prim_track2_length", &ssm_prim_track2_length, "ssm_prim_track2_length/F");
      fBDT->Branch("ssm_prim_track2_direct_length", &ssm_prim_track2_direct_length, "ssm_prim_track2_direct_length/F");
      fBDT->Branch("ssm_prim_track2_length_ratio", &ssm_prim_track2_length_ratio, "ssm_prim_track2_length_ratio/F");
      fBDT->Branch("ssm_prim_track2_max_dev", &ssm_prim_track2_max_dev, "ssm_prim_track2_max_dev/F");
      fBDT->Branch("ssm_prim_track2_kine_energy_range", &ssm_prim_track2_kine_energy_range, "ssm_prim_track2_kine_energy_range/F");
      fBDT->Branch("ssm_prim_track2_kine_energy_range_mu", &ssm_prim_track2_kine_energy_range_mu, "ssm_prim_track2_kine_energy_range_mu/F");
      fBDT->Branch("ssm_prim_track2_kine_energy_range_p", &ssm_prim_track2_kine_energy_range_p, "ssm_prim_track2_kine_energy_range_p/F");
      fBDT->Branch("ssm_prim_track2_kine_energy_range_e", &ssm_prim_track2_kine_energy_range_e, "ssm_prim_track2_kine_energy_range_e/F");
      fBDT->Branch("ssm_prim_track2_kine_energy_cal", &ssm_prim_track2_kine_energy_cal, "ssm_prim_track2_kine_energy_cal/F");
      fBDT->Branch("ssm_prim_track2_medium_dq_dx", &ssm_prim_track2_medium_dq_dx, "ssm_prim_track2_medium_dq_dx/F");
      fBDT->Branch("ssm_prim_track2_x_dir", &ssm_prim_track2_x_dir, "ssm_prim_track2_x_dir/F");
      fBDT->Branch("ssm_prim_track2_y_dir", &ssm_prim_track2_y_dir, "ssm_prim_track2_y_dir/F");
      fBDT->Branch("ssm_prim_track2_z_dir", &ssm_prim_track2_z_dir, "ssm_prim_track2_z_dir/F");
      fBDT->Branch("ssm_prim_track2_add_daught_track_counts_1", &ssm_prim_track2_add_daught_track_counts_1, "ssm_prim_track2_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_prim_track2_add_daught_all_counts_1", &ssm_prim_track2_add_daught_all_counts_1, "ssm_prim_track2_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_prim_track2_add_daught_track_counts_5", &ssm_prim_track2_add_daught_track_counts_5, "ssm_prim_track2_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_prim_track2_add_daught_all_counts_5", &ssm_prim_track2_add_daught_all_counts_5, "ssm_prim_track2_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_prim_track2_add_daught_track_counts_11", &ssm_prim_track2_add_daught_track_counts_11, "ssm_prim_track2_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_prim_track2_add_daught_all_counts_11", &ssm_prim_track2_add_daught_all_counts_11, "ssm_prim_track2_add_daught_all_counts_11/F");
    //properties of leading daughter track
      fBDT->Branch("ssm_daught_track1_pdg", &ssm_daught_track1_pdg, "ssm_daught_track1_pdg/F");
      fBDT->Branch("ssm_daught_track1_score_mu_fwd", &ssm_daught_track1_score_mu_fwd, "ssm_daught_track1_score_mu_fwd/F");
      fBDT->Branch("ssm_daught_track1_score_p_fwd", &ssm_daught_track1_score_p_fwd, "ssm_daught_track1_score_p_fwd/F");
      fBDT->Branch("ssm_daught_track1_score_e_fwd", &ssm_daught_track1_score_e_fwd, "ssm_daught_track1_score_e_fwd/F");
      fBDT->Branch("ssm_daught_track1_score_mu_bck", &ssm_daught_track1_score_mu_bck, "ssm_daught_track1_score_mu_bck/F");
      fBDT->Branch("ssm_daught_track1_score_p_bck", &ssm_daught_track1_score_p_bck, "ssm_daught_track1_score_p_bck/F");
      fBDT->Branch("ssm_daught_track1_score_e_bck", &ssm_daught_track1_score_e_bck, "ssm_daught_track1_score_e_bck/F");
      fBDT->Branch("ssm_daught_track1_length", &ssm_daught_track1_length, "ssm_daught_track1_length/F");
      fBDT->Branch("ssm_daught_track1_direct_length", &ssm_daught_track1_direct_length, "ssm_daught_track1_direct_length/F");
      fBDT->Branch("ssm_daught_track1_length_ratio", &ssm_daught_track1_length_ratio, "ssm_daught_track1_length_ratio/F");
      fBDT->Branch("ssm_daught_track1_max_dev", &ssm_daught_track1_max_dev, "ssm_daught_track1_max_dev/F");
      fBDT->Branch("ssm_daught_track1_kine_energy_range", &ssm_daught_track1_kine_energy_range, "ssm_daught_track1_kine_energy_range/F");
      fBDT->Branch("ssm_daught_track1_kine_energy_range_mu", &ssm_daught_track1_kine_energy_range_mu, "ssm_daught_track1_kine_energy_range_mu/F");
      fBDT->Branch("ssm_daught_track1_kine_energy_range_p", &ssm_daught_track1_kine_energy_range_p, "ssm_daught_track1_kine_energy_range_p/F");
      fBDT->Branch("ssm_daught_track1_kine_energy_range_e", &ssm_daught_track1_kine_energy_range_e, "ssm_daught_track1_kine_energy_range_e/F");
      fBDT->Branch("ssm_daught_track1_kine_energy_cal", &ssm_daught_track1_kine_energy_cal, "ssm_daught_track1_kine_energy_cal/F");
      fBDT->Branch("ssm_daught_track1_medium_dq_dx", &ssm_daught_track1_medium_dq_dx, "ssm_daught_track1_medium_dq_dx/F");
      fBDT->Branch("ssm_daught_track1_x_dir", &ssm_daught_track1_x_dir, "ssm_daught_track1_x_dir/F");
      fBDT->Branch("ssm_daught_track1_y_dir", &ssm_daught_track1_y_dir, "ssm_daught_track1_y_dir/F");
      fBDT->Branch("ssm_daught_track1_z_dir", &ssm_daught_track1_z_dir, "ssm_daught_track1_z_dir/F");
      fBDT->Branch("ssm_daught_track1_add_daught_track_counts_1", &ssm_daught_track1_add_daught_track_counts_1, "ssm_daught_track1_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_daught_track1_add_daught_all_counts_1", &ssm_daught_track1_add_daught_all_counts_1, "ssm_daught_track1_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_daught_track1_add_daught_track_counts_5", &ssm_daught_track1_add_daught_track_counts_5, "ssm_daught_track1_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_daught_track1_add_daught_all_counts_5", &ssm_daught_track1_add_daught_all_counts_5, "ssm_daught_track1_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_daught_track1_add_daught_track_counts_11", &ssm_daught_track1_add_daught_track_counts_11, "ssm_daught_track1_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_daught_track1_add_daught_all_counts_11", &ssm_daught_track1_add_daught_all_counts_11, "ssm_daught_track1_add_daught_all_counts_11/F");
    //properties of sub-leading daughter track
      fBDT->Branch("ssm_daught_track2_pdg", &ssm_daught_track2_pdg, "ssm_daught_track2_pdg/F");
      fBDT->Branch("ssm_daught_track2_score_mu_fwd", &ssm_daught_track2_score_mu_fwd, "ssm_daught_track2_score_mu_fwd/F");
      fBDT->Branch("ssm_daught_track2_score_p_fwd", &ssm_daught_track2_score_p_fwd, "ssm_daught_track2_score_p_fwd/F");
      fBDT->Branch("ssm_daught_track2_score_e_fwd", &ssm_daught_track2_score_e_fwd, "ssm_daught_track2_score_e_fwd/F");
      fBDT->Branch("ssm_daught_track2_score_mu_bck", &ssm_daught_track2_score_mu_bck, "ssm_daught_track2_score_mu_bck/F");
      fBDT->Branch("ssm_daught_track2_score_p_bck", &ssm_daught_track2_score_p_bck, "ssm_daught_track2_score_p_bck/F");
      fBDT->Branch("ssm_daught_track2_score_e_bck", &ssm_daught_track2_score_e_bck, "ssm_daught_track2_score_e_bck/F");
      fBDT->Branch("ssm_daught_track2_length", &ssm_daught_track2_length, "ssm_daught_track2_length/F");
      fBDT->Branch("ssm_daught_track2_direct_length", &ssm_daught_track2_direct_length, "ssm_daught_track2_direct_length/F");
      fBDT->Branch("ssm_daught_track2_length_ratio", &ssm_daught_track2_length_ratio, "ssm_daught_track2_length_ratio/F");
      fBDT->Branch("ssm_daught_track2_max_dev", &ssm_daught_track2_max_dev, "ssm_daught_track2_max_dev/F");
      fBDT->Branch("ssm_daught_track2_kine_energy_range", &ssm_daught_track2_kine_energy_range, "ssm_daught_track2_kine_energy_range/F");
      fBDT->Branch("ssm_daught_track2_kine_energy_range_mu", &ssm_daught_track2_kine_energy_range_mu, "ssm_daught_track2_kine_energy_range_mu/F");
      fBDT->Branch("ssm_daught_track2_kine_energy_range_p", &ssm_daught_track2_kine_energy_range_p, "ssm_daught_track2_kine_energy_range_p/F");
      fBDT->Branch("ssm_daught_track2_kine_energy_range_e", &ssm_daught_track2_kine_energy_range_e, "ssm_daught_track2_kine_energy_range_e/F");
      fBDT->Branch("ssm_daught_track2_kine_energy_cal", &ssm_daught_track2_kine_energy_cal, "ssm_daught_track2_kine_energy_cal/F");
      fBDT->Branch("ssm_daught_track2_medium_dq_dx", &ssm_daught_track2_medium_dq_dx, "ssm_daught_track2_medium_dq_dx/F");
      fBDT->Branch("ssm_daught_track2_x_dir", &ssm_daught_track2_x_dir, "ssm_daught_track2_x_dir/F");
      fBDT->Branch("ssm_daught_track2_y_dir", &ssm_daught_track2_y_dir, "ssm_daught_track2_y_dir/F");
      fBDT->Branch("ssm_daught_track2_z_dir", &ssm_daught_track2_z_dir, "ssm_daught_track2_z_dir/F");
      fBDT->Branch("ssm_daught_track2_add_daught_track_counts_1", &ssm_daught_track2_add_daught_track_counts_1, "ssm_daught_track2_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_daught_track2_add_daught_all_counts_1", &ssm_daught_track2_add_daught_all_counts_1, "ssm_daught_track2_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_daught_track2_add_daught_track_counts_5", &ssm_daught_track2_add_daught_track_counts_5, "ssm_daught_track2_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_daught_track2_add_daught_all_counts_5", &ssm_daught_track2_add_daught_all_counts_5, "ssm_daught_track2_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_daught_track2_add_daught_track_counts_11", &ssm_daught_track2_add_daught_track_counts_11, "ssm_daught_track2_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_daught_track2_add_daught_all_counts_11", &ssm_daught_track2_add_daught_all_counts_11, "ssm_daught_track2_add_daught_all_counts_11/F");
    //properties of leading other primary shower
      fBDT->Branch("ssm_prim_shw1_pdg", &ssm_prim_shw1_pdg, "ssm_prim_shw1_pdg/F");
      fBDT->Branch("ssm_prim_shw1_score_mu_fwd", &ssm_prim_shw1_score_mu_fwd, "ssm_prim_shw1_score_mu_fwd/F");
      fBDT->Branch("ssm_prim_shw1_score_p_fwd", &ssm_prim_shw1_score_p_fwd, "ssm_prim_shw1_score_p_fwd/F");
      fBDT->Branch("ssm_prim_shw1_score_e_fwd", &ssm_prim_shw1_score_e_fwd, "ssm_prim_shw1_score_e_fwd/F");
      fBDT->Branch("ssm_prim_shw1_score_mu_bck", &ssm_prim_shw1_score_mu_bck, "ssm_prim_shw1_score_mu_bck/F");
      fBDT->Branch("ssm_prim_shw1_score_p_bck", &ssm_prim_shw1_score_p_bck, "ssm_prim_shw1_score_p_bck/F");
      fBDT->Branch("ssm_prim_shw1_score_e_bck", &ssm_prim_shw1_score_e_bck, "ssm_prim_shw1_score_e_bck/F");
      fBDT->Branch("ssm_prim_shw1_length", &ssm_prim_shw1_length, "ssm_prim_shw1_length/F");
      fBDT->Branch("ssm_prim_shw1_direct_length", &ssm_prim_shw1_direct_length, "ssm_prim_shw1_direct_length/F");
      fBDT->Branch("ssm_prim_shw1_length_ratio", &ssm_prim_shw1_length_ratio, "ssm_prim_shw1_length_ratio/F");
      fBDT->Branch("ssm_prim_shw1_max_dev", &ssm_prim_shw1_max_dev, "ssm_prim_shw1_max_dev/F");
      fBDT->Branch("ssm_prim_shw1_kine_energy_range", &ssm_prim_shw1_kine_energy_range, "ssm_prim_shw1_kine_energy_range/F");
      fBDT->Branch("ssm_prim_shw1_kine_energy_range_mu", &ssm_prim_shw1_kine_energy_range_mu, "ssm_prim_shw1_kine_energy_range_mu/F");
      fBDT->Branch("ssm_prim_shw1_kine_energy_range_p", &ssm_prim_shw1_kine_energy_range_p, "ssm_prim_shw1_kine_energy_range_p/F");
      fBDT->Branch("ssm_prim_shw1_kine_energy_range_e", &ssm_prim_shw1_kine_energy_range_e, "ssm_prim_shw1_kine_energy_range_e/F");
      fBDT->Branch("ssm_prim_shw1_kine_energy_cal", &ssm_prim_shw1_kine_energy_cal, "ssm_prim_shw1_kine_energy_cal/F");
      fBDT->Branch("ssm_prim_shw1_kine_energy_best", &ssm_prim_shw1_kine_energy_best, "ssm_prim_shw1_kine_energy_best/F");
      fBDT->Branch("ssm_prim_shw1_medium_dq_dx", &ssm_prim_shw1_medium_dq_dx, "ssm_prim_shw1_medium_dq_dx/F");
      fBDT->Branch("ssm_prim_shw1_x_dir", &ssm_prim_shw1_x_dir, "ssm_prim_shw1_x_dir/F");
      fBDT->Branch("ssm_prim_shw1_y_dir", &ssm_prim_shw1_y_dir, "ssm_prim_shw1_y_dir/F");
      fBDT->Branch("ssm_prim_shw1_z_dir", &ssm_prim_shw1_z_dir, "ssm_prim_shw1_z_dir/F");
      fBDT->Branch("ssm_prim_shw1_add_daught_track_counts_1", &ssm_prim_shw1_add_daught_track_counts_1, "ssm_prim_shw1_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_prim_shw1_add_daught_all_counts_1", &ssm_prim_shw1_add_daught_all_counts_1, "ssm_prim_shw1_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_prim_shw1_add_daught_track_counts_5", &ssm_prim_shw1_add_daught_track_counts_5, "ssm_prim_shw1_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_prim_shw1_add_daught_all_counts_5", &ssm_prim_shw1_add_daught_all_counts_5, "ssm_prim_shw1_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_prim_shw1_add_daught_track_counts_11", &ssm_prim_shw1_add_daught_track_counts_11, "ssm_prim_shw1_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_prim_shw1_add_daught_all_counts_11", &ssm_prim_shw1_add_daught_all_counts_11, "ssm_prim_shw1_add_daught_all_counts_11/F");
    //properties of sub-leading other primary shower
      fBDT->Branch("ssm_prim_shw2_pdg", &ssm_prim_shw2_pdg, "ssm_prim_shw2_pdg/F");
      fBDT->Branch("ssm_prim_shw2_score_mu_fwd", &ssm_prim_shw2_score_mu_fwd, "ssm_prim_shw2_score_mu_fwd/F");
      fBDT->Branch("ssm_prim_shw2_score_p_fwd", &ssm_prim_shw2_score_p_fwd, "ssm_prim_shw2_score_p_fwd/F");
      fBDT->Branch("ssm_prim_shw2_score_e_fwd", &ssm_prim_shw2_score_e_fwd, "ssm_prim_shw2_score_e_fwd/F");
      fBDT->Branch("ssm_prim_shw2_score_mu_bck", &ssm_prim_shw2_score_mu_bck, "ssm_prim_shw2_score_mu_bck/F");
      fBDT->Branch("ssm_prim_shw2_score_p_bck", &ssm_prim_shw2_score_p_bck, "ssm_prim_shw2_score_p_bck/F");
      fBDT->Branch("ssm_prim_shw2_score_e_bck", &ssm_prim_shw2_score_e_bck, "ssm_prim_shw2_score_e_bck/F");
      fBDT->Branch("ssm_prim_shw2_length", &ssm_prim_shw2_length, "ssm_prim_shw2_length/F");
      fBDT->Branch("ssm_prim_shw2_direct_length", &ssm_prim_shw2_direct_length, "ssm_prim_shw2_direct_length/F");
      fBDT->Branch("ssm_prim_shw2_length_ratio", &ssm_prim_shw2_length_ratio, "ssm_prim_shw2_length_ratio/F");
      fBDT->Branch("ssm_prim_shw2_max_dev", &ssm_prim_shw2_max_dev, "ssm_prim_shw2_max_dev/F");
      fBDT->Branch("ssm_prim_shw2_kine_energy_range", &ssm_prim_shw2_kine_energy_range, "ssm_prim_shw2_kine_energy_range/F");
      fBDT->Branch("ssm_prim_shw2_kine_energy_range_mu", &ssm_prim_shw2_kine_energy_range_mu, "ssm_prim_shw2_kine_energy_range_mu/F");
      fBDT->Branch("ssm_prim_shw2_kine_energy_range_p", &ssm_prim_shw2_kine_energy_range_p, "ssm_prim_shw2_kine_energy_range_p/F");
      fBDT->Branch("ssm_prim_shw2_kine_energy_range_e", &ssm_prim_shw2_kine_energy_range_e, "ssm_prim_shw2_kine_energy_range_e/F");
      fBDT->Branch("ssm_prim_shw2_kine_energy_cal", &ssm_prim_shw2_kine_energy_cal, "ssm_prim_shw2_kine_energy_cal/F");
      fBDT->Branch("ssm_prim_shw2_kine_energy_best", &ssm_prim_shw2_kine_energy_best, "ssm_prim_shw2_kine_energy_best/F");
      fBDT->Branch("ssm_prim_shw2_medium_dq_dx", &ssm_prim_shw2_medium_dq_dx, "ssm_prim_shw2_medium_dq_dx/F");
      fBDT->Branch("ssm_prim_shw2_x_dir", &ssm_prim_shw2_x_dir, "ssm_prim_shw2_x_dir/F");
      fBDT->Branch("ssm_prim_shw2_y_dir", &ssm_prim_shw2_y_dir, "ssm_prim_shw2_y_dir/F");
      fBDT->Branch("ssm_prim_shw2_z_dir", &ssm_prim_shw2_z_dir, "ssm_prim_shw2_z_dir/F");
      fBDT->Branch("ssm_prim_shw2_add_daught_track_counts_1", &ssm_prim_shw2_add_daught_track_counts_1, "ssm_prim_shw2_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_prim_shw2_add_daught_all_counts_1", &ssm_prim_shw2_add_daught_all_counts_1, "ssm_prim_shw2_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_prim_shw2_add_daught_track_counts_5", &ssm_prim_shw2_add_daught_track_counts_5, "ssm_prim_shw2_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_prim_shw2_add_daught_all_counts_5", &ssm_prim_shw2_add_daught_all_counts_5, "ssm_prim_shw2_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_prim_shw2_add_daught_track_counts_11", &ssm_prim_shw2_add_daught_track_counts_11, "ssm_prim_shw2_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_prim_shw2_add_daught_all_counts_11", &ssm_prim_shw2_add_daught_all_counts_11, "ssm_prim_shw2_add_daught_all_counts_11/F");
    //properties of leading daughter shower
      fBDT->Branch("ssm_daught_shw1_pdg", &ssm_daught_shw1_pdg, "ssm_daught_shw1_pdg/F");
      fBDT->Branch("ssm_daught_shw1_score_mu_fwd", &ssm_daught_shw1_score_mu_fwd, "ssm_daught_shw1_score_mu_fwd/F");
      fBDT->Branch("ssm_daught_shw1_score_p_fwd", &ssm_daught_shw1_score_p_fwd, "ssm_daught_shw1_score_p_fwd/F");
      fBDT->Branch("ssm_daught_shw1_score_e_fwd", &ssm_daught_shw1_score_e_fwd, "ssm_daught_shw1_score_e_fwd/F");
      fBDT->Branch("ssm_daught_shw1_score_mu_bck", &ssm_daught_shw1_score_mu_bck, "ssm_daught_shw1_score_mu_bck/F");
      fBDT->Branch("ssm_daught_shw1_score_p_bck", &ssm_daught_shw1_score_p_bck, "ssm_daught_shw1_score_p_bck/F");
      fBDT->Branch("ssm_daught_shw1_score_e_bck", &ssm_daught_shw1_score_e_bck, "ssm_daught_shw1_score_e_bck/F");
      fBDT->Branch("ssm_daught_shw1_length", &ssm_daught_shw1_length, "ssm_daught_shw1_length/F");
      fBDT->Branch("ssm_daught_shw1_direct_length", &ssm_daught_shw1_direct_length, "ssm_daught_shw1_direct_length/F");
      fBDT->Branch("ssm_daught_shw1_length_ratio", &ssm_daught_shw1_length_ratio, "ssm_daught_shw1_length_ratio/F");
      fBDT->Branch("ssm_daught_shw1_max_dev", &ssm_daught_shw1_max_dev, "ssm_daught_shw1_max_dev/F");
      fBDT->Branch("ssm_daught_shw1_kine_energy_range", &ssm_daught_shw1_kine_energy_range, "ssm_daught_shw1_kine_energy_range/F");
      fBDT->Branch("ssm_daught_shw1_kine_energy_range_mu", &ssm_daught_shw1_kine_energy_range_mu, "ssm_daught_shw1_kine_energy_range_mu/F");
      fBDT->Branch("ssm_daught_shw1_kine_energy_range_p", &ssm_daught_shw1_kine_energy_range_p, "ssm_daught_shw1_kine_energy_range_p/F");
      fBDT->Branch("ssm_daught_shw1_kine_energy_range_e", &ssm_daught_shw1_kine_energy_range_e, "ssm_daught_shw1_kine_energy_range_e/F");
      fBDT->Branch("ssm_daught_shw1_kine_energy_cal", &ssm_daught_shw1_kine_energy_cal, "ssm_daught_shw1_kine_energy_cal/F");
      fBDT->Branch("ssm_daught_shw1_kine_energy_best", &ssm_daught_shw1_kine_energy_best, "ssm_daught_shw1_kine_energy_best/F");
      fBDT->Branch("ssm_daught_shw1_medium_dq_dx", &ssm_daught_shw1_medium_dq_dx, "ssm_daught_shw1_medium_dq_dx/F");
      fBDT->Branch("ssm_daught_shw1_x_dir", &ssm_daught_shw1_x_dir, "ssm_daught_shw1_x_dir/F");
      fBDT->Branch("ssm_daught_shw1_y_dir", &ssm_daught_shw1_y_dir, "ssm_daught_shw1_y_dir/F");
      fBDT->Branch("ssm_daught_shw1_z_dir", &ssm_daught_shw1_z_dir, "ssm_daught_shw1_z_dir/F");
      fBDT->Branch("ssm_daught_shw1_add_daught_track_counts_1", &ssm_daught_shw1_add_daught_track_counts_1, "ssm_daught_shw1_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_daught_shw1_add_daught_all_counts_1", &ssm_daught_shw1_add_daught_all_counts_1, "ssm_daught_shw1_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_daught_shw1_add_daught_track_counts_5", &ssm_daught_shw1_add_daught_track_counts_5, "ssm_daught_shw1_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_daught_shw1_add_daught_all_counts_5", &ssm_daught_shw1_add_daught_all_counts_5, "ssm_daught_shw1_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_daught_shw1_add_daught_track_counts_11", &ssm_daught_shw1_add_daught_track_counts_11, "ssm_daught_shw1_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_daught_shw1_add_daught_all_counts_11", &ssm_daught_shw1_add_daught_all_counts_11, "ssm_daught_shw1_add_daught_all_counts_11/F");
    //properties of sub-leading daughter shower
      fBDT->Branch("ssm_daught_shw2_pdg", &ssm_daught_shw2_pdg, "ssm_daught_shw2_pdg/F");
      fBDT->Branch("ssm_daught_shw2_score_mu_fwd", &ssm_daught_shw2_score_mu_fwd, "ssm_daught_shw2_score_mu_fwd/F");
      fBDT->Branch("ssm_daught_shw2_score_p_fwd", &ssm_daught_shw2_score_p_fwd, "ssm_daught_shw2_score_p_fwd/F");
      fBDT->Branch("ssm_daught_shw2_score_e_fwd", &ssm_daught_shw2_score_e_fwd, "ssm_daught_shw2_score_e_fwd/F");
      fBDT->Branch("ssm_daught_shw2_score_mu_bck", &ssm_daught_shw2_score_mu_bck, "ssm_daught_shw2_score_mu_bck/F");
      fBDT->Branch("ssm_daught_shw2_score_p_bck", &ssm_daught_shw2_score_p_bck, "ssm_daught_shw2_score_p_bck/F");
      fBDT->Branch("ssm_daught_shw2_score_e_bck", &ssm_daught_shw2_score_e_bck, "ssm_daught_shw2_score_e_bck/F");
      fBDT->Branch("ssm_daught_shw2_length", &ssm_daught_shw2_length, "ssm_daught_shw2_length/F");
      fBDT->Branch("ssm_daught_shw2_direct_length", &ssm_daught_shw2_direct_length, "ssm_daught_shw2_direct_length/F");
      fBDT->Branch("ssm_daught_shw2_length_ratio", &ssm_daught_shw2_length_ratio, "ssm_daught_shw2_length_ratio/F");
      fBDT->Branch("ssm_daught_shw2_max_dev", &ssm_daught_shw2_max_dev, "ssm_daught_shw2_max_dev/F");
      fBDT->Branch("ssm_daught_shw2_kine_energy_range", &ssm_daught_shw2_kine_energy_range, "ssm_daught_shw2_kine_energy_range/F");
      fBDT->Branch("ssm_daught_shw2_kine_energy_range_mu", &ssm_daught_shw2_kine_energy_range_mu, "ssm_daught_shw2_kine_energy_range_mu/F");
      fBDT->Branch("ssm_daught_shw2_kine_energy_range_p", &ssm_daught_shw2_kine_energy_range_p, "ssm_daught_shw2_kine_energy_range_p/F");
      fBDT->Branch("ssm_daught_shw2_kine_energy_range_e", &ssm_daught_shw2_kine_energy_range_e, "ssm_daught_shw2_kine_energy_range_e/F");
      fBDT->Branch("ssm_daught_shw2_kine_energy_cal", &ssm_daught_shw2_kine_energy_cal, "ssm_daught_shw2_kine_energy_cal/F");
      fBDT->Branch("ssm_daught_shw2_kine_energy_best", &ssm_daught_shw2_kine_energy_best, "ssm_daught_shw2_kine_energy_best/F");
      fBDT->Branch("ssm_daught_shw2_medium_dq_dx", &ssm_daught_shw2_medium_dq_dx, "ssm_daught_shw2_medium_dq_dx/F");
      fBDT->Branch("ssm_daught_shw2_x_dir", &ssm_daught_shw2_x_dir, "ssm_daught_shw2_x_dir/F");
      fBDT->Branch("ssm_daught_shw2_y_dir", &ssm_daught_shw2_y_dir, "ssm_daught_shw2_y_dir/F");
      fBDT->Branch("ssm_daught_shw2_z_dir", &ssm_daught_shw2_z_dir, "ssm_daught_shw2_z_dir/F");
      fBDT->Branch("ssm_daught_shw2_add_daught_track_counts_1", &ssm_daught_shw2_add_daught_track_counts_1, "ssm_daught_shw2_add_daught_track_counts_1/F");
      fBDT->Branch("ssm_daught_shw2_add_daught_all_counts_1", &ssm_daught_shw2_add_daught_all_counts_1, "ssm_daught_shw2_add_daught_all_counts_1/F");
      fBDT->Branch("ssm_daught_shw2_add_daught_track_counts_5", &ssm_daught_shw2_add_daught_track_counts_5, "ssm_daught_shw2_add_daught_track_counts_5/F");
      fBDT->Branch("ssm_daught_shw2_add_daught_all_counts_5", &ssm_daught_shw2_add_daught_all_counts_5, "ssm_daught_shw2_add_daught_all_counts_5/F");
      fBDT->Branch("ssm_daught_shw2_add_daught_track_counts_11", &ssm_daught_shw2_add_daught_track_counts_11, "ssm_daught_shw2_add_daught_track_counts_11/F");
      fBDT->Branch("ssm_daught_shw2_add_daught_all_counts_11", &ssm_daught_shw2_add_daught_all_counts_11, "ssm_daught_shw2_add_daught_all_counts_11/F");
    //event level properties
      fBDT->Branch("ssm_nu_angle_z", &ssm_nu_angle_z, "ssm_nu_angle_z/F");
      fBDT->Branch("ssm_nu_angle_target", &ssm_nu_angle_target, "ssm_nu_angle_target/F");
      fBDT->Branch("ssm_nu_angle_absorber", &ssm_nu_angle_absorber, "ssm_nu_angle_absorber/F");
      fBDT->Branch("ssm_nu_angle_vertical", &ssm_nu_angle_vertical, "ssm_nu_angle_vertical/F");
      fBDT->Branch("ssm_con_nu_angle_z", &ssm_con_nu_angle_z, "ssm_con_nu_angle_z/F");
      fBDT->Branch("ssm_con_nu_angle_target", &ssm_con_nu_angle_target, "ssm_con_nu_angle_target/F");
      fBDT->Branch("ssm_con_nu_angle_absorber", &ssm_con_nu_angle_absorber, "ssm_con_nu_angle_absorber/F");
      fBDT->Branch("ssm_con_nu_angle_vertical", &ssm_con_nu_angle_vertical, "ssm_con_nu_angle_vertical/F");
      fBDT->Branch("ssm_prim_nu_angle_z", &ssm_prim_nu_angle_z, "ssm_prim_nu_angle_z/F");
      fBDT->Branch("ssm_prim_nu_angle_target", &ssm_prim_nu_angle_target, "ssm_prim_nu_angle_target/F");
      fBDT->Branch("ssm_prim_nu_angle_absorber", &ssm_prim_nu_angle_absorber, "ssm_prim_nu_angle_absorber/F");
      fBDT->Branch("ssm_prim_nu_angle_vertical", &ssm_prim_nu_angle_vertical, "ssm_prim_nu_angle_vertical/F");
      fBDT->Branch("ssm_track_angle_z", &ssm_track_angle_z, "ssm_track_angle_z/F");
      fBDT->Branch("ssm_track_angle_target", &ssm_track_angle_target, "ssm_track_angle_target/F");
      fBDT->Branch("ssm_track_angle_absorber", &ssm_track_angle_absorber, "ssm_track_angle_absorber/F");
      fBDT->Branch("ssm_track_angle_vertical", &ssm_track_angle_vertical, "ssm_track_angle_vertical/F");
      fBDT->Branch("ssm_vtxX", &ssm_vtxX, "ssm_vtxX/F");
      fBDT->Branch("ssm_vtxY", &ssm_vtxY, "ssm_vtxY/F");
      fBDT->Branch("ssm_vtxZ", &ssm_vtxZ, "ssm_vtxZ/F");
    //off vertex stuff
      fBDT->Branch("ssm_offvtx_length",&ssm_offvtx_length,"ssm_offvtx_length/F");
      fBDT->Branch("ssm_offvtx_energy",&ssm_offvtx_energy,"ssm_offvtx_energy/F");
      fBDT->Branch("ssm_n_offvtx_tracks_1",&ssm_n_offvtx_tracks_1,"ssm_n_offvtx_tracks_1/F");
      fBDT->Branch("ssm_n_offvtx_tracks_3",&ssm_n_offvtx_tracks_3,"ssm_n_offvtx_tracks_3/F");
      fBDT->Branch("ssm_n_offvtx_tracks_5",&ssm_n_offvtx_tracks_5,"ssm_n_offvtx_tracks_5/F");
      fBDT->Branch("ssm_n_offvtx_tracks_8",&ssm_n_offvtx_tracks_8,"ssm_n_offvtx_tracks_8/F");
      fBDT->Branch("ssm_n_offvtx_tracks_11",&ssm_n_offvtx_tracks_11,"ssm_n_offvtx_tracks_11/F");
      fBDT->Branch("ssm_n_offvtx_showers_1",&ssm_n_offvtx_showers_1,"ssm_n_offvtx_showers_1/F");
      fBDT->Branch("ssm_n_offvtx_showers_3",&ssm_n_offvtx_showers_3,"ssm_n_offvtx_showers_3/F");
      fBDT->Branch("ssm_n_offvtx_showers_5",&ssm_n_offvtx_showers_5,"ssm_n_offvtx_showers_5/F");
      fBDT->Branch("ssm_n_offvtx_showers_8",&ssm_n_offvtx_showers_8,"ssm_n_offvtx_showers_8/F");
      fBDT->Branch("ssm_n_offvtx_showers_11",&ssm_n_offvtx_showers_11,"ssm_n_offvtx_showers_11/F");
    //properties of leading off vertex track
      fBDT->Branch("ssm_offvtx_track1_pdg",&ssm_offvtx_track1_pdg,"ssm_offvtx_track1_pdg/F");
      fBDT->Branch("ssm_offvtx_track1_score_mu_fwd",&ssm_offvtx_track1_score_mu_fwd,"ssm_offvtx_track1_score_mu_fwd/F");
      fBDT->Branch("ssm_offvtx_track1_score_p_fwd",&ssm_offvtx_track1_score_p_fwd,"ssm_offvtx_track1_score_p_fwd/F");
      fBDT->Branch("ssm_offvtx_track1_score_e_fwd",&ssm_offvtx_track1_score_e_fwd,"ssm_offvtx_track1_score_e_fwd/F");
      fBDT->Branch("ssm_offvtx_track1_score_mu_bck",&ssm_offvtx_track1_score_mu_bck,"ssm_offvtx_track1_score_mu_bck/F");
      fBDT->Branch("ssm_offvtx_track1_score_p_bck",&ssm_offvtx_track1_score_p_bck,"ssm_offvtx_track1_score_p_bck/F");
      fBDT->Branch("ssm_offvtx_track1_score_e_bck",&ssm_offvtx_track1_score_e_bck,"ssm_offvtx_track1_score_e_bck/F");
      fBDT->Branch("ssm_offvtx_track1_length",&ssm_offvtx_track1_length,"ssm_offvtx_track1_length/F");
      fBDT->Branch("ssm_offvtx_track1_direct_length",&ssm_offvtx_track1_direct_length,"ssm_offvtx_track1_direct_length/F");
      fBDT->Branch("ssm_offvtx_track1_max_dev",&ssm_offvtx_track1_max_dev,"ssm_offvtx_track1_max_dev/F");
      fBDT->Branch("ssm_offvtx_track1_kine_energy_range",&ssm_offvtx_track1_kine_energy_range,"ssm_offvtx_track1_kine_energy_range/F");
      fBDT->Branch("ssm_offvtx_track1_kine_energy_range_mu",&ssm_offvtx_track1_kine_energy_range_mu,"ssm_offvtx_track1_kine_energy_range_mu/F");
      fBDT->Branch("ssm_offvtx_track1_kine_energy_range_p",&ssm_offvtx_track1_kine_energy_range_p,"ssm_offvtx_track1_kine_energy_range_p/F");
      fBDT->Branch("ssm_offvtx_track1_kine_energy_range_e",&ssm_offvtx_track1_kine_energy_range_e,"ssm_offvtx_track1_kine_energy_range_e/F");
      fBDT->Branch("ssm_offvtx_track1_kine_energy_cal",&ssm_offvtx_track1_kine_energy_cal,"ssm_offvtx_track1_kine_energy_cal/F");
      fBDT->Branch("ssm_offvtx_track1_medium_dq_dx",&ssm_offvtx_track1_medium_dq_dx,"ssm_offvtx_track1_medium_dq_dx/F");
      fBDT->Branch("ssm_offvtx_track1_x_dir",&ssm_offvtx_track1_x_dir,"ssm_offvtx_track1_x_dir/F");
      fBDT->Branch("ssm_offvtx_track1_y_dir",&ssm_offvtx_track1_y_dir,"ssm_offvtx_track1_y_dir/F");
      fBDT->Branch("ssm_offvtx_track1_z_dir",&ssm_offvtx_track1_z_dir,"ssm_offvtx_track1_z_dir/F");
      fBDT->Branch("ssm_offvtx_track1_dist_mainvtx",&ssm_offvtx_track1_dist_mainvtx,"ssm_offvtx_track1_dist_mainvtx/F");
    //properties of leading off vertex shower
      fBDT->Branch("ssm_offvtx_shw1_pdg_offvtx",&ssm_offvtx_shw1_pdg_offvtx,"ssm_offvtx_shw1_pdg_offvtx/F");
      fBDT->Branch("ssm_offvtx_shw1_score_mu_fwd",&ssm_offvtx_shw1_score_mu_fwd,"ssm_offvtx_shw1_score_mu_fwd/F");
      fBDT->Branch("ssm_offvtx_shw1_score_p_fwd",&ssm_offvtx_shw1_score_p_fwd,"ssm_offvtx_shw1_score_p_fwd/F");
      fBDT->Branch("ssm_offvtx_shw1_score_e_fwd",&ssm_offvtx_shw1_score_e_fwd,"ssm_offvtx_shw1_score_e_fwd/F");
      fBDT->Branch("ssm_offvtx_shw1_score_mu_bck",&ssm_offvtx_shw1_score_mu_bck,"ssm_offvtx_shw1_score_mu_bck/F");
      fBDT->Branch("ssm_offvtx_shw1_score_p_bck",&ssm_offvtx_shw1_score_p_bck,"ssm_offvtx_shw1_score_p_bck/F");
      fBDT->Branch("ssm_offvtx_shw1_score_e_bck",&ssm_offvtx_shw1_score_e_bck,"ssm_offvtx_shw1_score_e_bck/F");
      fBDT->Branch("ssm_offvtx_shw1_length",&ssm_offvtx_shw1_length,"ssm_offvtx_shw1_length/F");
      fBDT->Branch("ssm_offvtx_shw1_direct_length",&ssm_offvtx_shw1_direct_length,"ssm_offvtx_shw1_direct_length/F");
      fBDT->Branch("ssm_offvtx_shw1_max_dev",&ssm_offvtx_shw1_max_dev,"ssm_offvtx_shw1_max_dev/F");
      fBDT->Branch("ssm_offvtx_shw1_kine_energy_best",&ssm_offvtx_shw1_kine_energy_best,"ssm_offvtx_shw1_kine_energy_best/F");
      fBDT->Branch("ssm_offvtx_shw1_kine_energy_range",&ssm_offvtx_shw1_kine_energy_range,"ssm_offvtx_shw1_kine_energy_range/F");
      fBDT->Branch("ssm_offvtx_shw1_kine_energy_range_mu",&ssm_offvtx_shw1_kine_energy_range_mu,"ssm_offvtx_shw1_kine_energy_range_mu/F");
      fBDT->Branch("ssm_offvtx_shw1_kine_energy_range_p",&ssm_offvtx_shw1_kine_energy_range_p,"ssm_offvtx_shw1_kine_energy_range_p/F");
      fBDT->Branch("ssm_offvtx_shw1_kine_energy_range_e",&ssm_offvtx_shw1_kine_energy_range_e,"ssm_offvtx_shw1_kine_energy_range_e/F");
      fBDT->Branch("ssm_offvtx_shw1_kine_energy_cal",&ssm_offvtx_shw1_kine_energy_cal,"ssm_offvtx_shw1_kine_energy_cal/F");
      fBDT->Branch("ssm_offvtx_shw1_medium_dq_dx",&ssm_offvtx_shw1_medium_dq_dx,"ssm_offvtx_shw1_medium_dq_dx/F");
      fBDT->Branch("ssm_offvtx_shw1_x_dir",&ssm_offvtx_shw1_x_dir,"ssm_offvtx_shw1_x_dir/F");
      fBDT->Branch("ssm_offvtx_shw1_y_dir",&ssm_offvtx_shw1_y_dir,"ssm_offvtx_shw1_y_dir/F");
      fBDT->Branch("ssm_offvtx_shw1_z_dir",&ssm_offvtx_shw1_z_dir,"ssm_offvtx_shw1_z_dir/F");
      fBDT->Branch("ssm_offvtx_shw1_dist_mainvtx",&ssm_offvtx_shw1_dist_mainvtx,"ssm_offvtx_shw1_dist_mainvtx/F");
    // Sapcepoints
      fBDT->Branch("ssmsp_Ntrack", &ssmsp_Ntrack, "ssmsp_Ntrack/I");
      fBDT->Branch("ssmsp_Nsp", &ssmsp_Nsp);
      fBDT->Branch("ssmsp_Nsp_tot", &ssmsp_Nsp_tot, "ssmsp_Nsp_tot/I");
      fBDT->Branch("ssmsp_pdg", &ssmsp_pdg);
      fBDT->Branch("ssmsp_id", &ssmsp_id);
      fBDT->Branch("ssmsp_mother", &ssmsp_mother);
      fBDT->Branch("ssmsp_x", &ssmsp_x);
      fBDT->Branch("ssmsp_y", &ssmsp_y);
      fBDT->Branch("ssmsp_z", &ssmsp_z);
      fBDT->Branch("ssmsp_dx", &ssmsp_dx);
      fBDT->Branch("ssmsp_dQ", &ssmsp_dQ);
      fBDT->Branch("ssmsp_KE", &ssmsp_KE);
      fBDT->Branch("ssmsp_containing_shower_id", &ssmsp_containing_shower_id);
      fBDT->Branch("ssmsp_containing_shower_ke", &ssmsp_containing_shower_ke);
      fBDT->Branch("ssmsp_containing_shower_flag", &ssmsp_containing_shower_flag);
   // Kine vars
      if(f_KINEvars){
        fBDT->Branch("ssm_kine_reco_Enu",&ssm_kine_reco_Enu,"ssm_kine_reco_Enu/F");
        fBDT->Branch("ssm_kine_reco_add_energy",&ssm_kine_reco_add_energy,"ssm_kine_reco_add_energy/F");
        fBDT->Branch("ssm_kine_energy_particle",&ssm_kine_energy_particle);
        fBDT->Branch("ssm_kine_energy_info",&ssm_kine_energy_info);
        fBDT->Branch("ssm_kine_particle_type",&ssm_kine_particle_type);
        fBDT->Branch("ssm_kine_energy_included",&ssm_kine_energy_included);
        fBDT->Branch("ssm_kine_pio_mass",&ssm_kine_pio_mass,"ssm_kine_pio_mass/F");
        fBDT->Branch("ssm_kine_pio_flag",&ssm_kine_pio_flag,"ssm_kine_pio_flag/I");
        fBDT->Branch("ssm_kine_pio_vtx_dis",&ssm_kine_pio_vtx_dis,"ssm_kine_pio_vtx_dis/F");
        fBDT->Branch("ssm_kine_pio_energy_1",&ssm_kine_pio_energy_1,"ssm_kine_pio_energy_1/F");
        fBDT->Branch("ssm_kine_pio_theta_1",&ssm_kine_pio_theta_1,"ssm_kine_pio_theta_1/F");
        fBDT->Branch("ssm_kine_pio_phi_1",&ssm_kine_pio_phi_1,"ssm_kine_pio_phi_1/F");
        fBDT->Branch("ssm_kine_pio_dis_1",&ssm_kine_pio_dis_1,"ssm_kine_pio_dis_1/F");
        fBDT->Branch("ssm_kine_pio_energy_2",&ssm_kine_pio_energy_2,"ssm_kine_pio_energy_2/F");
        fBDT->Branch("ssm_kine_pio_theta_2",&ssm_kine_pio_theta_2,"ssm_kine_pio_theta_2/F");
        fBDT->Branch("ssm_kine_pio_phi_2",&ssm_kine_pio_phi_2,"ssm_kine_pio_phi_2/F");
        fBDT->Branch("ssm_kine_pio_dis_2",&ssm_kine_pio_dis_2,"ssm_kine_pio_dis_2/F");
        fBDT->Branch("ssm_kine_pio_angle",&ssm_kine_pio_angle,"ssm_kine_pio_angle/F");
      }
      fBDT->Branch("ssm_numu_cc_flag",&ssm_numu_cc_flag);
      fBDT->Branch("ssm_cosmict_flag_1",&ssm_cosmict_flag_1,"ssm_cosmict_flag_1/F");
      fBDT->Branch("ssm_cosmict_flag_2",&ssm_cosmict_flag_2,"ssm_cosmict_flag_2/F");
      fBDT->Branch("ssm_cosmict_flag_3",&ssm_cosmict_flag_3,"ssm_cosmict_flag_3/F");
      fBDT->Branch("ssm_cosmict_flag_4",&ssm_cosmict_flag_4,"ssm_cosmict_flag_4/F");
      fBDT->Branch("ssm_cosmict_flag_5",&ssm_cosmict_flag_5,"ssm_cosmict_flag_5/F");
      fBDT->Branch("ssm_cosmict_flag_6",&ssm_cosmict_flag_6,"ssm_cosmict_flag_6/F");
      fBDT->Branch("ssm_cosmict_flag_7",&ssm_cosmict_flag_7,"ssm_cosmict_flag_7/F");
      fBDT->Branch("ssm_cosmict_flag_8",&ssm_cosmict_flag_8,"ssm_cosmict_flag_8/F");
      fBDT->Branch("ssm_cosmict_flag_9",&ssm_cosmict_flag_9,"ssm_cosmict_flag_9/F");
      fBDT->Branch("ssm_cosmict_flag_10",&ssm_cosmict_flag_10);
      fBDT->Branch("ssm_cosmict_flag",&ssm_cosmict_flag,"ssm_cosmict_flag/F");
  }

  //single photon shower
  fBDT->Branch("shw_sp_flag",&shw_sp_flag,"shw_sp_flag/F");
  fBDT->Branch("shw_sp_num_mip_tracks",&shw_sp_num_mip_tracks,"shw_sp_num_mip_tracks/F");
  fBDT->Branch("shw_sp_num_muons",&shw_sp_num_muons,"shw_sp_num_muons/F");
  fBDT->Branch("shw_sp_num_pions",&shw_sp_num_pions,"shw_sp_num_pions/F");
  fBDT->Branch("shw_sp_num_protons",&shw_sp_num_protons,"shw_sp_num_protons/F");
  fBDT->Branch("shw_sp_proton_length_1",&shw_sp_proton_length_1,"shw_sp_proton_length_1/F");
  fBDT->Branch("shw_sp_proton_dqdx_1",&shw_sp_proton_dqdx_1,"shw_sp_proton_dqdx_1/F");
  fBDT->Branch("shw_sp_proton_energy_1",&shw_sp_proton_energy_1,"shw_sp_proton_energy_1/F");
  fBDT->Branch("shw_sp_proton_length_2",&shw_sp_proton_length_2,"shw_sp_proton_length_2/F");
  fBDT->Branch("shw_sp_proton_dqdx_2",&shw_sp_proton_dqdx_2,"shw_sp_proton_dqdx_2/F");
  fBDT->Branch("shw_sp_proton_energy_2",&shw_sp_proton_energy_2,"shw_sp_proton_energy_2/F");
  fBDT->Branch("shw_sp_n_good_showers",&shw_sp_n_good_showers,"shw_sp_n_good_showers/F");
  fBDT->Branch("shw_sp_n_20mev_showers",&shw_sp_n_20mev_showers,"shw_sp_n_20mev_showers/F");
  fBDT->Branch("shw_sp_n_br1_showers",&shw_sp_n_br1_showers,"shw_sp_n_br1_showers/F");
  fBDT->Branch("shw_sp_n_br2_showers",&shw_sp_n_br2_showers,"shw_sp_n_br2_showers/F");
  fBDT->Branch("shw_sp_n_br3_showers",&shw_sp_n_br3_showers,"shw_sp_n_br3_showers/F");
  fBDT->Branch("shw_sp_n_br4_showers",&shw_sp_n_br4_showers,"shw_sp_n_br4_showers/F");
  fBDT->Branch("shw_sp_n_20br1_showers",&shw_sp_n_20br1_showers,"shw_sp_n_20br1_showers/F");
  fBDT->Branch("shw_sp_20mev_showers",&shw_sp_20mev_showers);
  fBDT->Branch("shw_sp_br1_showers",&shw_sp_br1_showers);
  fBDT->Branch("shw_sp_br2_showers",&shw_sp_br2_showers);
  fBDT->Branch("shw_sp_br3_showers",&shw_sp_br3_showers);
  fBDT->Branch("shw_sp_br4_showers",&shw_sp_br4_showers);
  fBDT->Branch("shw_sp_shw_vtx_dis",&shw_sp_shw_vtx_dis,"shw_sp_shw_vtx_dis/F");
  fBDT->Branch("shw_sp_max_shw_dis",&shw_sp_max_shw_dis,"shw_sp_max_shw_dis/F");
  fBDT->Branch("shw_sp_energy",&shw_sp_energy,"shw_sp_energy/F");

  fBDT->Branch("shw_sp_vec_median_dedx",&shw_sp_vec_median_dedx,"shw_sp_vec_median_dedx/F");
  fBDT->Branch("shw_sp_vec_mean_dedx",&shw_sp_vec_mean_dedx,"shw_sp_vec_mean_dedx/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_0",&shw_sp_vec_dQ_dx_0,"shw_sp_vec_dQ_dx_0/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_1",&shw_sp_vec_dQ_dx_1,"shw_sp_vec_dQ_dx_1/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_2",&shw_sp_vec_dQ_dx_2,"shw_sp_vec_dQ_dx_2/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_3",&shw_sp_vec_dQ_dx_3,"shw_sp_vec_dQ_dx_3/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_4",&shw_sp_vec_dQ_dx_4,"shw_sp_vec_dQ_dx_4/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_5",&shw_sp_vec_dQ_dx_5,"shw_sp_vec_dQ_dx_5/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_6",&shw_sp_vec_dQ_dx_6,"shw_sp_vec_dQ_dx_6/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_7",&shw_sp_vec_dQ_dx_7,"shw_sp_vec_dQ_dx_7/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_8",&shw_sp_vec_dQ_dx_8,"shw_sp_vec_dQ_dx_8/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_9",&shw_sp_vec_dQ_dx_9,"shw_sp_vec_dQ_dx_9/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_10",&shw_sp_vec_dQ_dx_10,"shw_sp_vec_dQ_dx_10/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_11",&shw_sp_vec_dQ_dx_11,"shw_sp_vec_dQ_dx_11/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_12",&shw_sp_vec_dQ_dx_12,"shw_sp_vec_dQ_dx_12/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_13",&shw_sp_vec_dQ_dx_13,"shw_sp_vec_dQ_dx_13/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_14",&shw_sp_vec_dQ_dx_14,"shw_sp_vec_dQ_dx_14/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_15",&shw_sp_vec_dQ_dx_15,"shw_sp_vec_dQ_dx_15/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_16",&shw_sp_vec_dQ_dx_16,"shw_sp_vec_dQ_dx_16/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_17",&shw_sp_vec_dQ_dx_17,"shw_sp_vec_dQ_dx_17/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_18",&shw_sp_vec_dQ_dx_18,"shw_sp_vec_dQ_dx_18/F");
  fBDT->Branch("shw_sp_vec_dQ_dx_19",&shw_sp_vec_dQ_dx_19,"shw_sp_vec_dQ_dx_19/F");

  fBDT->Branch("shw_sp_max_dQ_dx_sample",&shw_sp_max_dQ_dx_sample,"shw_sp_max_dQ_dx_sample/F");
  fBDT->Branch("shw_sp_n_below_threshold",&shw_sp_n_below_threshold,"shw_sp_n_below_threshold/F");
  fBDT->Branch("shw_sp_n_below_zero",&shw_sp_n_below_zero,"shw_sp_n_below_zero/F");
  fBDT->Branch("shw_sp_n_lowest",&shw_sp_n_lowest,"shw_sp_n_lowest/F");
  fBDT->Branch("shw_sp_n_highest",&shw_sp_n_highest,"shw_sp_n_highest/F");

  fBDT->Branch("shw_sp_lowest_dQ_dx",&shw_sp_lowest_dQ_dx,"shw_sp_lowest_dQ_dx/F");
  fBDT->Branch("shw_sp_highest_dQ_dx",&shw_sp_highest_dQ_dx,"shw_sp_highest_dQ_dx/F");
  fBDT->Branch("shw_sp_medium_dQ_dx",&shw_sp_medium_dQ_dx,"shw_sp_medium_dQ_dx/F");
  fBDT->Branch("shw_sp_stem_length",&shw_sp_stem_length,"shw_sp_stem_length/F");
  fBDT->Branch("shw_sp_length_main",&shw_sp_length_main,"shw_sp_length_main/F");
  fBDT->Branch("shw_sp_length_total",&shw_sp_length_total,"shw_sp_length_total/F");
  fBDT->Branch("shw_sp_angle_beam",&shw_sp_angle_beam,"shw_sp_angle_beam/F");
  fBDT->Branch("shw_sp_iso_angle",&shw_sp_iso_angle,"shw_sp_iso_angle/F");

  fBDT->Branch("shw_sp_n_vertex",&shw_sp_n_vertex,"shw_sp_n_vertex/F");
  fBDT->Branch("shw_sp_n_good_tracks",&shw_sp_n_good_tracks,"shw_sp_n_good_tracks/F");
  fBDT->Branch("shw_sp_E_indirect_max_energy",&shw_sp_E_indirect_max_energy,"shw_sp_E_indirect_max_energy/F");
  fBDT->Branch("shw_sp_flag_all_above",&shw_sp_flag_all_above,"shw_sp_flag_all_above/F");
  fBDT->Branch("shw_sp_min_dQ_dx_5",&shw_sp_min_dQ_dx_5,"shw_sp_min_dQ_dx_5/F");
  fBDT->Branch("shw_sp_n_other_vertex",&shw_sp_n_other_vertex,"shw_sp_n_other_vertex/F");
  fBDT->Branch("shw_sp_n_stem_size",&shw_sp_n_stem_size,"shw_sp_n_stem_size/F");
  fBDT->Branch("shw_sp_flag_stem_trajectory",&shw_sp_flag_stem_trajectory,"shw_sp_flag_stem_trajectory/F");
  fBDT->Branch("shw_sp_min_dis",&shw_sp_min_dis,"shw_sp_min_dis/F");
  fBDT->Branch("shw_sp_filled",&shw_sp_filled,"shw_sp_filled/F");
  // pio ...
  fBDT->Branch("shw_sp_pio_flag",&shw_sp_pio_flag,"shw_sp_pio_flag/F");
  fBDT->Branch("shw_sp_pio_mip_id",&shw_sp_pio_mip_id,"shw_sp_pio_mip_id/F");
  fBDT->Branch("shw_sp_pio_filled",&shw_sp_pio_filled,"shw_sp_pio_filled/F");
  fBDT->Branch("shw_sp_pio_flag_pio",&shw_sp_pio_flag_pio,"shw_sp_pio_flag_pio/F");

  fBDT->Branch("shw_sp_pio_1_flag",&shw_sp_pio_1_flag,"shw_sp_pio_1_flag/F");
  fBDT->Branch("shw_sp_pio_1_mass",&shw_sp_pio_1_mass,"shw_sp_pio_1_mass/F");
  fBDT->Branch("shw_sp_pio_1_pio_type",&shw_sp_pio_1_pio_type,"shw_sp_pio_1_pio_type/F");
  fBDT->Branch("shw_sp_pio_1_energy_1",&shw_sp_pio_1_energy_1,"shw_sp_pio_1_energy_1/F");
  fBDT->Branch("shw_sp_pio_1_energy_2",&shw_sp_pio_1_energy_2,"shw_sp_pio_1_energy_2/F");
  fBDT->Branch("shw_sp_pio_1_dis_1",&shw_sp_pio_1_dis_1,"shw_sp_pio_1_dis_1/F");
  fBDT->Branch("shw_sp_pio_1_dis_2",&shw_sp_pio_1_dis_2,"shw_sp_pio_1_dis_2/F");

  fBDT->Branch("shw_sp_pio_2_v_dis2",&shw_sp_pio_2_v_dis2);
  fBDT->Branch("shw_sp_pio_2_v_angle2",&shw_sp_pio_2_v_angle2);
  fBDT->Branch("shw_sp_pio_2_v_acc_length",&shw_sp_pio_2_v_acc_length);
  fBDT->Branch("shw_sp_pio_2_v_flag",&shw_sp_pio_2_v_flag);

  fBDT->Branch("shw_sp_br_filled",&shw_sp_br_filled,"shw_sp_br_filled/F");

  fBDT->Branch("shw_sp_br1_flag",&shw_sp_br1_flag,"shw_sp_br1_flag/F");

  fBDT->Branch("shw_sp_br1_1_flag",&shw_sp_br1_1_flag,"shw_sp_br1_1_flag/F");
  fBDT->Branch("shw_sp_br1_1_shower_type",&shw_sp_br1_1_shower_type,"shw_sp_br1_1_shower_type/F");
  fBDT->Branch("shw_sp_br1_1_vtx_n_segs",&shw_sp_br1_1_vtx_n_segs,"shw_sp_br1_1_vtx_n_segs/F");
  fBDT->Branch("shw_sp_br1_1_energy",&shw_sp_br1_1_energy,"shw_sp_br1_1_energy/F");
  fBDT->Branch("shw_sp_br1_1_n_segs",&shw_sp_br1_1_n_segs,"shw_sp_br1_1_n_segs/F");
  fBDT->Branch("shw_sp_br1_1_flag_sg_topology",&shw_sp_br1_1_flag_sg_topology,"shw_sp_br1_1_flag_sg_topology/F");
  fBDT->Branch("shw_sp_br1_1_flag_sg_trajectory",&shw_sp_br1_1_flag_sg_trajectory,"shw_sp_br1_1_flag_sg_trajectory/F");
  fBDT->Branch("shw_sp_br1_1_sg_length",&shw_sp_br1_1_sg_length,"shw_sp_br1_1_sg_length/F");

  fBDT->Branch("shw_sp_br1_2_flag",&shw_sp_br1_2_flag,"shw_sp_br1_2_flag/F");
  fBDT->Branch("shw_sp_br1_2_energy",&shw_sp_br1_2_energy,"shw_sp_br1_2_energy/F");
  fBDT->Branch("shw_sp_br1_2_n_connected",&shw_sp_br1_2_n_connected,"shw_sp_br1_2_n_connected/F");
  fBDT->Branch("shw_sp_br1_2_max_length",&shw_sp_br1_2_max_length,"shw_sp_br1_2_max_length/F");
  fBDT->Branch("shw_sp_br1_2_n_connected_1",&shw_sp_br1_2_n_connected_1,"shw_sp_br1_2_n_connected_1/F");
  fBDT->Branch("shw_sp_br1_2_vtx_n_segs",&shw_sp_br1_2_vtx_n_segs,"shw_sp_br1_2_vtx_n_segs/F");
  fBDT->Branch("shw_sp_br1_2_n_shower_segs",&shw_sp_br1_2_n_shower_segs,"shw_sp_br1_2_n_shower_segs/F");
  fBDT->Branch("shw_sp_br1_2_max_length_ratio",&shw_sp_br1_2_max_length_ratio,"shw_sp_br1_2_max_length_ratio/F");
  fBDT->Branch("shw_sp_br1_2_shower_length",&shw_sp_br1_2_shower_length,"shw_sp_br1_2_shower_length/F");

  fBDT->Branch("shw_sp_br1_3_flag",&shw_sp_br1_3_flag,"shw_sp_br1_3_flag/F");
  fBDT->Branch("shw_sp_br1_3_energy",&shw_sp_br1_3_energy,"shw_sp_br1_3_energy/F");
  fBDT->Branch("shw_sp_br1_3_n_connected_p",&shw_sp_br1_3_n_connected_p,"shw_sp_br1_3_n_connected_p/F");
  fBDT->Branch("shw_sp_br1_3_max_length_p",&shw_sp_br1_3_max_length_p,"shw_sp_br1_3_max_length_p/F");
  fBDT->Branch("shw_sp_br1_3_n_shower_segs",&shw_sp_br1_3_n_shower_segs,"shw_sp_br1_3_n_shower_segs/F");
  fBDT->Branch("shw_sp_br1_3_flag_sg_topology",&shw_sp_br1_3_flag_sg_topology,"shw_sp_br1_3_flag_sg_topology/F");
  fBDT->Branch("shw_sp_br1_3_flag_sg_trajectory",&shw_sp_br1_3_flag_sg_trajectory,"shw_sp_br1_3_flag_sg_trajectory/F");
  fBDT->Branch("shw_sp_br1_3_n_shower_main_segs",&shw_sp_br1_3_n_shower_main_segs,"shw_sp_br1_3_n_shower_main_segs/F");
  fBDT->Branch("shw_sp_br1_3_sg_length",&shw_sp_br1_3_sg_length,"shw_sp_br1_3_sg_length/F");

  fBDT->Branch("shw_sp_br2_flag",&shw_sp_br2_flag,"shw_sp_br2_flag/F");
  fBDT->Branch("shw_sp_br2_flag_single_shower",&shw_sp_br2_flag_single_shower,"shw_sp_br2_flag_single_shower/F");
  fBDT->Branch("shw_sp_br2_num_valid_tracks",&shw_sp_br2_num_valid_tracks,"shw_sp_br2_num_valid_tracks/F");
  fBDT->Branch("shw_sp_br2_energy",&shw_sp_br2_energy,"shw_sp_br2_energy/F");
  fBDT->Branch("shw_sp_br2_angle1",&shw_sp_br2_angle1,"shw_sp_br2_angle1/F");
  fBDT->Branch("shw_sp_br2_angle2",&shw_sp_br2_angle2,"shw_sp_br2_angle2/F");
  fBDT->Branch("shw_sp_br2_angle",&shw_sp_br2_angle,"shw_sp_br2_angle/F");
  fBDT->Branch("shw_sp_br2_angle3",&shw_sp_br2_angle3,"shw_sp_br2_angle3/F");
  fBDT->Branch("shw_sp_br2_n_shower_main_segs",&shw_sp_br2_n_shower_main_segs,"shw_sp_br2_n_shower_main_segs/F");
  fBDT->Branch("shw_sp_br2_max_angle",&shw_sp_br2_max_angle,"shw_sp_br2_max_angle/F");
  fBDT->Branch("shw_sp_br2_sg_length",&shw_sp_br2_sg_length,"shw_sp_br2_sg_length/F");
  fBDT->Branch("shw_sp_br2_flag_sg_trajectory",&shw_sp_br2_flag_sg_trajectory,"shw_sp_br2_flag_sg_trajectory/F");


  fBDT->Branch("shw_sp_lol_flag",&shw_sp_lol_flag,"shw_sp_lol_flag/F");

  fBDT->Branch("shw_sp_lol_1_v_energy",&shw_sp_lol_1_v_energy);
  fBDT->Branch("shw_sp_lol_1_v_vtx_n_segs",&shw_sp_lol_1_v_vtx_n_segs);
  fBDT->Branch("shw_sp_lol_1_v_nseg",&shw_sp_lol_1_v_nseg);
  fBDT->Branch("shw_sp_lol_1_v_angle",&shw_sp_lol_1_v_angle);
  fBDT->Branch("shw_sp_lol_1_v_flag",&shw_sp_lol_1_v_flag);

  fBDT->Branch("shw_sp_lol_2_v_length",&shw_sp_lol_2_v_length);
  fBDT->Branch("shw_sp_lol_2_v_angle",&shw_sp_lol_2_v_angle);
  fBDT->Branch("shw_sp_lol_2_v_type",&shw_sp_lol_2_v_type);
  fBDT->Branch("shw_sp_lol_2_v_vtx_n_segs",&shw_sp_lol_2_v_vtx_n_segs);
  fBDT->Branch("shw_sp_lol_2_v_energy",&shw_sp_lol_2_v_energy);
  fBDT->Branch("shw_sp_lol_2_v_shower_main_length",&shw_sp_lol_2_v_shower_main_length);
  fBDT->Branch("shw_sp_lol_2_v_flag_dir_weak",&shw_sp_lol_2_v_flag_dir_weak);
  fBDT->Branch("shw_sp_lol_2_v_flag",&shw_sp_lol_2_v_flag);

  fBDT->Branch("shw_sp_lol_3_angle_beam",&shw_sp_lol_3_angle_beam,"shw_sp_lol_3_angle_beam/F");
  fBDT->Branch("shw_sp_lol_3_n_valid_tracks",&shw_sp_lol_3_n_valid_tracks,"shw_sp_lol_3_n_valid_tracks/F");
  fBDT->Branch("shw_sp_lol_3_min_angle",&shw_sp_lol_3_min_angle,"shw_sp_lol_3_min_angle/F");
  fBDT->Branch("shw_sp_lol_3_vtx_n_segs",&shw_sp_lol_3_vtx_n_segs,"shw_sp_lol_3_vtx_n_segs/F");
  fBDT->Branch("shw_sp_lol_3_energy",&shw_sp_lol_3_energy,"shw_sp_lol_3_energy/F");
  fBDT->Branch("shw_sp_lol_3_shower_main_length",&shw_sp_lol_3_shower_main_length,"shw_sp_lol_3_shower_main_length/F");
  fBDT->Branch("shw_sp_lol_3_n_out",&shw_sp_lol_3_n_out,"shw_sp_lol_3_n_out/F");
  fBDT->Branch("shw_sp_lol_3_n_sum",&shw_sp_lol_3_n_sum,"shw_sp_lol_3_n_sum/F");
  fBDT->Branch("shw_sp_lol_3_flag",&shw_sp_lol_3_flag,"shw_sp_lol_3_flag/F");

  fBDT->Branch("shw_sp_br3_1_energy",&shw_sp_br3_1_energy,"shw_sp_br3_1_energy/F");
  fBDT->Branch("shw_sp_br3_1_n_shower_segments",&shw_sp_br3_1_n_shower_segments,"shw_sp_br3_1_n_shower_segments/F");
  fBDT->Branch("shw_sp_br3_1_sg_flag_trajectory",&shw_sp_br3_1_sg_flag_trajectory,"shw_sp_br3_1_sg_flag_trajectory/F");
  fBDT->Branch("shw_sp_br3_1_sg_direct_length",&shw_sp_br3_1_sg_direct_length,"shw_sp_br3_1_sg_direct_length/F");
  fBDT->Branch("shw_sp_br3_1_sg_length",&shw_sp_br3_1_sg_length,"shw_sp_br3_1_sg_length/F");
  fBDT->Branch("shw_sp_br3_1_total_main_length",&shw_sp_br3_1_total_main_length,"shw_sp_br3_1_total_main_length/F");
  fBDT->Branch("shw_sp_br3_1_total_length",&shw_sp_br3_1_total_length,"shw_sp_br3_1_total_length/F");
  fBDT->Branch("shw_sp_br3_1_iso_angle",&shw_sp_br3_1_iso_angle,"shw_sp_br3_1_iso_angle/F");
  fBDT->Branch("shw_sp_br3_1_sg_flag_topology",&shw_sp_br3_1_sg_flag_topology,"shw_sp_br3_1_sg_flag_topology/F");
  fBDT->Branch("shw_sp_br3_1_flag",&shw_sp_br3_1_flag,"shw_sp_br3_1_flag/F");

  fBDT->Branch("shw_sp_br3_2_n_ele",&shw_sp_br3_2_n_ele,"shw_sp_br3_2_n_ele/F");
  fBDT->Branch("shw_sp_br3_2_n_other",&shw_sp_br3_2_n_other,"shw_sp_br3_2_n_other/F");
  fBDT->Branch("shw_sp_br3_2_energy",&shw_sp_br3_2_energy,"shw_sp_br3_2_energy/F");
  fBDT->Branch("shw_sp_br3_2_total_main_length",&shw_sp_br3_2_total_main_length,"shw_sp_br3_2_total_main_length/F");
  fBDT->Branch("shw_sp_br3_2_total_length",&shw_sp_br3_2_total_length,"shw_sp_br3_2_total_length/F");
  fBDT->Branch("shw_sp_br3_2_other_fid",&shw_sp_br3_2_other_fid,"shw_sp_br3_2_other_fid/F");
  fBDT->Branch("shw_sp_br3_2_flag",&shw_sp_br3_2_flag,"shw_sp_br3_2_flag/F");

  fBDT->Branch("shw_sp_br3_3_v_energy",&shw_sp_br3_3_v_energy);
  fBDT->Branch("shw_sp_br3_3_v_angle",&shw_sp_br3_3_v_angle);
  fBDT->Branch("shw_sp_br3_3_v_dir_length",&shw_sp_br3_3_v_dir_length);
  fBDT->Branch("shw_sp_br3_3_v_length",&shw_sp_br3_3_v_length);
  fBDT->Branch("shw_sp_br3_3_v_flag",&shw_sp_br3_3_v_flag);

  fBDT->Branch("shw_sp_br3_4_acc_length", &shw_sp_br3_4_acc_length, "shw_sp_br3_4_acc_length/F");
  fBDT->Branch("shw_sp_br3_4_total_length", &shw_sp_br3_4_total_length, "shw_sp_br3_4_total_length/F");
  fBDT->Branch("shw_sp_br3_4_energy", &shw_sp_br3_4_energy, "shw_sp_br3_4_energy/F");
  fBDT->Branch("shw_sp_br3_4_flag", &shw_sp_br3_4_flag, "shw_sp_br3_4_flag/F");

  fBDT->Branch("shw_sp_br3_5_v_dir_length", &shw_sp_br3_5_v_dir_length);
  fBDT->Branch("shw_sp_br3_5_v_total_length", &shw_sp_br3_5_v_total_length);
  fBDT->Branch("shw_sp_br3_5_v_flag_avoid_muon_check", &shw_sp_br3_5_v_flag_avoid_muon_check);
  fBDT->Branch("shw_sp_br3_5_v_n_seg", &shw_sp_br3_5_v_n_seg);
  fBDT->Branch("shw_sp_br3_5_v_angle", &shw_sp_br3_5_v_angle);
  fBDT->Branch("shw_sp_br3_5_v_sg_length", &shw_sp_br3_5_v_sg_length);
  fBDT->Branch("shw_sp_br3_5_v_energy", &shw_sp_br3_5_v_energy);
  fBDT->Branch("shw_sp_br3_5_v_n_main_segs", &shw_sp_br3_5_v_n_main_segs);
  fBDT->Branch("shw_sp_br3_5_v_n_segs", &shw_sp_br3_5_v_n_segs);
  fBDT->Branch("shw_sp_br3_5_v_shower_main_length", &shw_sp_br3_5_v_shower_main_length);
  fBDT->Branch("shw_sp_br3_5_v_shower_total_length", &shw_sp_br3_5_v_shower_total_length);
  fBDT->Branch("shw_sp_br3_5_v_flag", &shw_sp_br3_5_v_flag);

  fBDT->Branch("shw_sp_br3_6_v_angle",&shw_sp_br3_6_v_angle);
  fBDT->Branch("shw_sp_br3_6_v_angle1",&shw_sp_br3_6_v_angle1);
  fBDT->Branch("shw_sp_br3_6_v_flag_shower_trajectory",&shw_sp_br3_6_v_flag_shower_trajectory);
  fBDT->Branch("shw_sp_br3_6_v_direct_length",&shw_sp_br3_6_v_direct_length);
  fBDT->Branch("shw_sp_br3_6_v_length",&shw_sp_br3_6_v_length);
  fBDT->Branch("shw_sp_br3_6_v_n_other_vtx_segs",&shw_sp_br3_6_v_n_other_vtx_segs);
  fBDT->Branch("shw_sp_br3_6_v_energy",&shw_sp_br3_6_v_energy);
  fBDT->Branch("shw_sp_br3_6_v_flag",&shw_sp_br3_6_v_flag);

  fBDT->Branch("shw_sp_br3_7_energy",&shw_sp_br3_7_energy,"shw_sp_br3_7_energy/F");
  fBDT->Branch("shw_sp_br3_7_min_angle",&shw_sp_br3_7_min_angle,"shw_sp_br3_7_min_angle/F");
  fBDT->Branch("shw_sp_br3_7_sg_length",&shw_sp_br3_7_sg_length,"shw_sp_br3_7_sg_length/F");
  fBDT->Branch("shw_sp_br3_7_main_length",&shw_sp_br3_7_shower_main_length,"shw_sp_br3_7_shower_main_length/F");
  fBDT->Branch("shw_sp_br3_7_flag",&shw_sp_br3_7_flag,"shw_sp_br3_7_flag/F");

  fBDT->Branch("shw_sp_br3_8_max_dQ_dx",&shw_sp_br3_8_max_dQ_dx,"shw_sp_br3_8_max_dQ_dx/F");
  fBDT->Branch("shw_sp_br3_8_energy",&shw_sp_br3_8_energy,"shw_sp_br3_8_energy/F");
  fBDT->Branch("shw_sp_br3_8_n_main_segs",&shw_sp_br3_8_n_main_segs,"shw_sp_br3_8_n_main_segs/F");
  fBDT->Branch("shw_sp_br3_8_shower_main_length",&shw_sp_br3_8_shower_main_length,"shw_sp_br3_8_shower_main_length/F");
  fBDT->Branch("shw_sp_br3_8_shower_length",&shw_sp_br3_8_shower_length,"shw_sp_br3_8_shower_length/F");
  fBDT->Branch("shw_sp_br3_8_flag",&shw_sp_br3_8_flag,"shw_sp_br3_8_flag/F");

  fBDT->Branch("shw_sp_br3_flag",&shw_sp_br3_flag,"shw_sp_br3_flag/F");


  fBDT->Branch("shw_sp_br4_1_shower_main_length", &shw_sp_br4_1_shower_main_length, "shw_sp_br4_1_shower_main_length/F");
  fBDT->Branch("shw_sp_br4_1_shower_total_length", &shw_sp_br4_1_shower_total_length, "shw_sp_br4_1_shower_total_length/F");
  fBDT->Branch("shw_sp_br4_1_min_dis", &shw_sp_br4_1_min_dis, "shw_sp_br4_1_min_dis/F");
  fBDT->Branch("shw_sp_br4_1_energy", &shw_sp_br4_1_energy, "shw_sp_br4_1_energy/F");
  fBDT->Branch("shw_sp_br4_1_flag_avoid_muon_check", &shw_sp_br4_1_flag_avoid_muon_check, "shw_sp_br4_1_flag_avoid_muon_check/F");
  fBDT->Branch("shw_sp_br4_1_n_vtx_segs", &shw_sp_br4_1_n_vtx_segs, "shw_sp_br4_1_n_vtx_segs/F");
  fBDT->Branch("shw_sp_br4_1_n_main_segs", &shw_sp_br4_1_n_main_segs, "shw_sp_br4_1_n_main_segs/F");
  fBDT->Branch("shw_sp_br4_1_flag", &shw_sp_br4_1_flag, "shw_sp_br4_1_flag/F");

  fBDT->Branch("shw_sp_br4_2_ratio_45", &shw_sp_br4_2_ratio_45, "shw_sp_br4_2_ratio_45/F");
  fBDT->Branch("shw_sp_br4_2_ratio_35", &shw_sp_br4_2_ratio_35, "shw_sp_br4_2_ratio_35/F");
  fBDT->Branch("shw_sp_br4_2_ratio_25", &shw_sp_br4_2_ratio_25, "shw_sp_br4_2_ratio_25/F");
  fBDT->Branch("shw_sp_br4_2_ratio_15", &shw_sp_br4_2_ratio_15, "shw_sp_br4_2_ratio_15/F");
  fBDT->Branch("shw_sp_br4_2_energy",   &shw_sp_br4_2_energy, "shw_sp_br4_2_energy/F");
  fBDT->Branch("shw_sp_br4_2_ratio1_45", &shw_sp_br4_2_ratio1_45, "shw_sp_br4_2_ratio1_45/F");
  fBDT->Branch("shw_sp_br4_2_ratio1_35", &shw_sp_br4_2_ratio1_35, "shw_sp_br4_2_ratio1_35/F");
  fBDT->Branch("shw_sp_br4_2_ratio1_25", &shw_sp_br4_2_ratio1_25, "shw_sp_br4_2_ratio1_25/F");
  fBDT->Branch("shw_sp_br4_2_ratio1_15", &shw_sp_br4_2_ratio1_15, "shw_sp_br4_2_ratio1_15/F");
  fBDT->Branch("shw_sp_br4_2_iso_angle", &shw_sp_br4_2_iso_angle, "shw_sp_br4_2_iso_angle/F");
  fBDT->Branch("shw_sp_br4_2_iso_angle1", &shw_sp_br4_2_iso_angle1, "shw_sp_br4_2_iso_angle1/F");
  fBDT->Branch("shw_sp_br4_2_angle", &shw_sp_br4_2_angle, "shw_sp_br4_2_angle/F");
  fBDT->Branch("shw_sp_br4_2_flag", &shw_sp_br4_2_flag, "shw_sp_br4_2_flag/F");

  fBDT->Branch("shw_sp_br4_flag", &shw_sp_br4_flag, "shw_sp_br4_flag/F");


  fBDT->Branch("shw_sp_hol_1_n_valid_tracks", &shw_sp_hol_1_n_valid_tracks,"shw_sp_hol_1_n_valid_tracks/F");
  fBDT->Branch("shw_sp_hol_1_min_angle", &shw_sp_hol_1_min_angle,"shw_sp_hol_1_min_angle/F");
  fBDT->Branch("shw_sp_hol_1_energy", &shw_sp_hol_1_energy,"shw_sp_hol_1_energy/F");
  fBDT->Branch("shw_sp_hol_1_flag_all_shower", &shw_sp_hol_1_flag_all_shower,"shw_sp_hol_1_flag_all_shower/F");
  fBDT->Branch("shw_sp_hol_1_min_length", &shw_sp_hol_1_min_length,"shw_sp_hol_1_min_length/F");
  fBDT->Branch("shw_sp_hol_1_flag", &shw_sp_hol_1_flag,"shw_sp_hol_1_flag/F");

  fBDT->Branch("shw_sp_hol_2_min_angle", &shw_sp_hol_2_min_angle,"shw_sp_hol_2_min_angle/F");
  fBDT->Branch("shw_sp_hol_2_medium_dQ_dx", &shw_sp_hol_2_medium_dQ_dx,"shw_sp_hol_2_medium_dQ_dx/F");
  fBDT->Branch("shw_sp_hol_2_ncount", &shw_sp_hol_2_ncount,"shw_sp_hol_2_ncount/F");
  fBDT->Branch("shw_sp_hol_2_energy", &shw_sp_hol_2_energy,"shw_sp_hol_2_energy/F");
  fBDT->Branch("shw_sp_hol_2_flag", &shw_sp_hol_2_flag,"shw_sp_hol_2_flag/F");

  fBDT->Branch("shw_sp_hol_flag", &shw_sp_hol_flag,"shw_sp_hol_flag/F");

  fBDT->Branch("shw_sp_lem_shower_total_length",&shw_sp_lem_shower_total_length,"shw_sp_lem_shower_total_length/F");
  fBDT->Branch("shw_sp_lem_shower_main_length",&shw_sp_lem_shower_main_length,"shw_sp_lem_shower_main_length/F");
  fBDT->Branch("shw_sp_lem_n_3seg",&shw_sp_lem_n_3seg,"shw_sp_lem_n_3seg/F");
  fBDT->Branch("shw_sp_lem_e_charge",&shw_sp_lem_e_charge,"shw_sp_lem_e_charge/F");
  fBDT->Branch("shw_sp_lem_e_dQdx",&shw_sp_lem_e_dQdx,"shw_sp_lem_e_dQdx/F");
  fBDT->Branch("shw_sp_lem_shower_num_segs",&shw_sp_lem_shower_num_segs,"shw_sp_lem_shower_num_segs/F");
  fBDT->Branch("shw_sp_lem_shower_num_main_segs",&shw_sp_lem_shower_num_main_segs,"shw_sp_lem_shower_num_main_segs/F");
  fBDT->Branch("shw_sp_lem_flag",&shw_sp_lem_flag,"shw_sp_lem_flag/F");


  fBDT->Branch("cosmic_flag", &cosmic_flag, "cosmic_flag/F");
  fBDT->Branch("cosmic_n_solid_tracks",&cosmic_n_solid_tracks,"cosmic_n_solid_tracks/F");
  fBDT->Branch("cosmic_energy_main_showers",&cosmic_energy_main_showers,"cosmic_energy_main_showers/F");
  fBDT->Branch("cosmic_energy_direct_showers",&cosmic_energy_direct_showers,"cosmic_energy_direct_showers/F");
  fBDT->Branch("cosmic_energy_indirect_showers",&cosmic_energy_indirect_showers,"cosmic_energy_indirect_showers/F");
  fBDT->Branch("cosmic_n_direct_showers",&cosmic_n_direct_showers,"cosmic_n_direct_showers/F");
  fBDT->Branch("cosmic_n_indirect_showers",&cosmic_n_indirect_showers,"cosmic_n_indirect_showers/F");
  fBDT->Branch("cosmic_n_main_showers",&cosmic_n_main_showers,"cosmic_n_main_showers/F");
  fBDT->Branch("cosmic_filled",&cosmic_filled,"cosmic_filled/F");

  fBDT->Branch("gap_flag",&gap_flag,"gap_flag/F");
  fBDT->Branch("gap_flag_prolong_u",&gap_flag_prolong_u,"gap_flag_prolong_u/F");
  fBDT->Branch("gap_flag_prolong_v",&gap_flag_prolong_v,"gap_flag_prolong_v/F");
  fBDT->Branch("gap_flag_prolong_w",&gap_flag_prolong_w,"gap_flag_prolong_w/F");
  fBDT->Branch("gap_flag_parallel",&gap_flag_parallel,"gap_flag_parallel/F");
  fBDT->Branch("gap_n_points",&gap_n_points,"gap_n_points/F");
  fBDT->Branch("gap_n_bad",&gap_n_bad,"gap_n_bad/F");
  fBDT->Branch("gap_energy",&gap_energy,"gap_energy/F");
  fBDT->Branch("gap_num_valid_tracks",&gap_num_valid_tracks,"gap_num_valid_tracks/F");
  fBDT->Branch("gap_flag_single_shower",&gap_flag_single_shower,"gap_flag_single_shower/F");
  fBDT->Branch("gap_filled",&gap_filled,"gap_filled/F");

  fBDT->Branch("mip_quality_flag",&mip_quality_flag,"mip_quality_flag/F");
  fBDT->Branch("mip_quality_energy",&mip_quality_energy,"mip_quality_energy/F");
  fBDT->Branch("mip_quality_overlap",&mip_quality_overlap,"mip_quality_overlap/F");
  fBDT->Branch("mip_quality_n_showers",&mip_quality_n_showers,"mip_quality_n_showers/F");
  fBDT->Branch("mip_quality_n_tracks",&mip_quality_n_tracks,"mip_quality_n_tracks/F");
  fBDT->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0,"mip_quality_flag_inside_pi0/F");
  fBDT->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers,"mip_quality_n_pi0_showers/F");
  fBDT->Branch("mip_quality_shortest_length",&mip_quality_shortest_length,"mip_quality_shortest_length/F");
  fBDT->Branch("mip_quality_acc_length",&mip_quality_acc_length,"mip_quality_acc_length/F");
  fBDT->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle,"mip_quality_shortest_angle/F");
  fBDT->Branch("mip_quality_flag_proton",&mip_quality_flag_proton,"mip_quality_flag_proton/F");
  fBDT->Branch("mip_quality_filled",&mip_quality_flag,"mip_quality_filled/F");

  fBDT->Branch("mip_flag",&mip_flag,"mip_flag/F");
  fBDT->Branch("mip_energy",&mip_energy,"mip_energy/F");
  fBDT->Branch("mip_n_end_reduction",&mip_n_end_reduction,"mip_n_end_reduction/F");
  fBDT->Branch("mip_n_first_mip",&mip_n_first_mip,"mip_n_first_mip/F");
  fBDT->Branch("mip_n_first_non_mip",&mip_n_first_non_mip,"mip_n_first_non_mip/F");
  fBDT->Branch("mip_n_first_non_mip_1",&mip_n_first_non_mip_1,"mip_n_first_non_mip_1/F");
  fBDT->Branch("mip_n_first_non_mip_2",&mip_n_first_non_mip_2,"mip_n_first_non_mip_2/F");
  fBDT->Branch("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0,"mip_vec_dQ_dx_0/F");
  fBDT->Branch("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1,"mip_vec_dQ_dx_1/F");
  fBDT->Branch("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample,"mip_max_dQ_dx_sample/F");
  fBDT->Branch("mip_n_below_threshold",&mip_n_below_threshold,"mip_n_below_threshold/F");
  fBDT->Branch("mip_n_below_zero",&mip_n_below_zero,"mip_n_below_zero/F");
  fBDT->Branch("mip_n_lowest",&mip_n_lowest,"mip_n_lowest/F");
  fBDT->Branch("mip_n_highest",&mip_n_highest,"mip_n_highest/F");
  fBDT->Branch("mip_lowest_dQ_dx",&mip_lowest_dQ_dx,"mip_lowest_dQ_dx/F");
  fBDT->Branch("mip_highest_dQ_dx",&mip_highest_dQ_dx,"mip_highest_dQ_dx/F");
  fBDT->Branch("mip_medium_dQ_dx",&mip_medium_dQ_dx,"mip_medium_dQ_dx/F");
  fBDT->Branch("mip_stem_length",&mip_stem_length,"mip_stem_length/F");
  fBDT->Branch("mip_length_main",&mip_length_main,"mip_length_main/F");
  fBDT->Branch("mip_length_total",&mip_length_total,"mip_length_total/F");
  fBDT->Branch("mip_angle_beam",&mip_angle_beam,"mip_angle_beam/F");
  fBDT->Branch("mip_iso_angle",&mip_iso_angle,"mip_iso_angle/F");
  fBDT->Branch("mip_n_vertex",&mip_n_vertex,"mip_n_vertex/F");
  fBDT->Branch("mip_n_good_tracks",&mip_n_good_tracks,"mip_n_good_tracks/F");
  fBDT->Branch("mip_E_indirect_max_energy",&mip_E_indirect_max_energy,"mip_E_indirect_max_energy/F");
  fBDT->Branch("mip_flag_all_above",&mip_flag_all_above,"mip_flag_all_above/F");
  fBDT->Branch("mip_min_dQ_dx_5",&mip_min_dQ_dx_5,"mip_min_dQ_dx_5/F");
  fBDT->Branch("mip_n_other_vertex",&mip_n_other_vertex,"mip_n_other_vertex/F");
  fBDT->Branch("mip_n_stem_size",&mip_n_stem_size,"mip_n_stem_size/F");
  fBDT->Branch("mip_flag_stem_trajectory",&mip_flag_stem_trajectory,"mip_flag_stem_trajectory/F");
  fBDT->Branch("mip_min_dis",&mip_min_dis,"mip_min_dis/F");
  fBDT->Branch("mip_filled",&mip_filled,"mip_filled/F");

  fBDT->Branch("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2,"mip_vec_dQ_dx_2/F");
  fBDT->Branch("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3,"mip_vec_dQ_dx_3/F");
  fBDT->Branch("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4,"mip_vec_dQ_dx_4/F");
  fBDT->Branch("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5,"mip_vec_dQ_dx_5/F");
  fBDT->Branch("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6,"mip_vec_dQ_dx_6/F");
  fBDT->Branch("mip_vec_dQ_dx_7",&mip_vec_dQ_dx_7,"mip_vec_dQ_dx_7/F");
  fBDT->Branch("mip_vec_dQ_dx_8",&mip_vec_dQ_dx_8,"mip_vec_dQ_dx_8/F");
  fBDT->Branch("mip_vec_dQ_dx_9",&mip_vec_dQ_dx_9,"mip_vec_dQ_dx_9/F");
  fBDT->Branch("mip_vec_dQ_dx_10",&mip_vec_dQ_dx_10,"mip_vec_dQ_dx_10/F");
  fBDT->Branch("mip_vec_dQ_dx_11",&mip_vec_dQ_dx_11,"mip_vec_dQ_dx_11/F");
  fBDT->Branch("mip_vec_dQ_dx_12",&mip_vec_dQ_dx_12,"mip_vec_dQ_dx_12/F");
  fBDT->Branch("mip_vec_dQ_dx_13",&mip_vec_dQ_dx_13,"mip_vec_dQ_dx_13/F");
  fBDT->Branch("mip_vec_dQ_dx_14",&mip_vec_dQ_dx_14,"mip_vec_dQ_dx_14/F");
  fBDT->Branch("mip_vec_dQ_dx_15",&mip_vec_dQ_dx_15,"mip_vec_dQ_dx_15/F");
  fBDT->Branch("mip_vec_dQ_dx_16",&mip_vec_dQ_dx_16,"mip_vec_dQ_dx_16/F");
  fBDT->Branch("mip_vec_dQ_dx_17",&mip_vec_dQ_dx_17,"mip_vec_dQ_dx_17/F");
  fBDT->Branch("mip_vec_dQ_dx_18",&mip_vec_dQ_dx_18,"mip_vec_dQ_dx_18/F");
  fBDT->Branch("mip_vec_dQ_dx_19",&mip_vec_dQ_dx_19,"mip_vec_dQ_dx_19/F");

  fBDT->Branch("pio_flag",&pio_flag,"pio_flag/F");
  fBDT->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  fBDT->Branch("pio_filled",&pio_filled,"pio_filled/F");
  fBDT->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");
  fBDT->Branch("pio_1_flag",&pio_1_flag,"pio_1_flag/F");
  fBDT->Branch("pio_1_mass",&pio_1_mass,"pio_1_mass/F");
  fBDT->Branch("pio_1_pio_type",&pio_1_pio_type,"pio_1_pio_type/F");
  fBDT->Branch("pio_1_energy_1",&pio_1_energy_1,"pio_1_energy_1/F");
  fBDT->Branch("pio_1_energy_2",&pio_1_energy_2,"pio_1_energy_2/F");
  fBDT->Branch("pio_1_dis_1",&pio_1_dis_1,"pio_1_dis_1/F");
  fBDT->Branch("pio_1_dis_2",&pio_1_dis_2,"pio_1_dis_2/F");
  fBDT->Branch("pio_2_v_dis2",&pio_2_v_dis2);
  fBDT->Branch("pio_2_v_angle2",&pio_2_v_angle2);
  fBDT->Branch("pio_2_v_acc_length",&pio_2_v_acc_length);
  fBDT->Branch("pio_2_v_flag",&pio_2_v_flag);

  fBDT->Branch("sig_1_v_angle",&sig_1_v_angle);
  fBDT->Branch("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
  fBDT->Branch("sig_1_v_energy",&sig_1_v_energy);
  fBDT->Branch("sig_1_v_energy_1",&sig_1_v_energy_1);
  fBDT->Branch("sig_1_v_flag",&sig_1_v_flag);
  fBDT->Branch("sig_2_v_energy",&sig_2_v_energy);
  fBDT->Branch("sig_2_v_shower_angle",&sig_2_v_shower_angle);
  fBDT->Branch("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
  fBDT->Branch("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);
  fBDT->Branch("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);
  fBDT->Branch("sig_2_v_flag",&sig_2_v_flag);
  fBDT->Branch("sig_flag",&sig_flag, "sig_flag/F");

  fBDT->Branch("mgo_energy",&mgo_energy,"mgo_energy/F");
  fBDT->Branch("mgo_max_energy",&mgo_max_energy,"mgo_max_energy/F");
  fBDT->Branch("mgo_total_energy",&mgo_total_energy,"mgo_total_energy/F");
  fBDT->Branch("mgo_n_showers",&mgo_n_showers,"mgo_n_showers/F");
  fBDT->Branch("mgo_max_energy_1",&mgo_max_energy_1,"mgo_max_energy_1/F");
  fBDT->Branch("mgo_max_energy_2",&mgo_max_energy_2,"mgo_max_energy_2/F");
  fBDT->Branch("mgo_total_other_energy",&mgo_total_other_energy,"mgo_total_other_energy/F");
  fBDT->Branch("mgo_n_total_showers",&mgo_n_total_showers,"mgo_n_total_showers/F");
  fBDT->Branch("mgo_total_other_energy_1",&mgo_total_other_energy_1,"mgo_total_other_energy_1/F");
  fBDT->Branch("mgo_flag",&mgo_flag,"mgo_flag/F");

  fBDT->Branch("mgt_flag_single_shower",&mgt_flag_single_shower,"mgt_flag_single_shower/F");
  fBDT->Branch("mgt_max_energy",&mgt_max_energy,"mgt_max_energy/F");
  fBDT->Branch("mgt_energy",&mgt_energy,"mgt_energy/F");
  fBDT->Branch("mgt_total_other_energy",&mgt_total_other_energy,"mgt_total_other_energy/F");
  fBDT->Branch("mgt_max_energy_1",&mgt_max_energy_1,"mgt_max_energy_1/F");
  fBDT->Branch("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy,"mgt_e_indirect_max_energy/F");
  fBDT->Branch("mgt_e_direct_max_energy",&mgt_e_direct_max_energy,"mgt_e_direct_max_energy/F");
  fBDT->Branch("mgt_n_direct_showers",&mgt_n_direct_showers,"mgt_n_direct_showers/F");
  fBDT->Branch("mgt_e_direct_total_energy",&mgt_e_direct_total_energy,"mgt_e_direct_total_energy/F");
  fBDT->Branch("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio,"mgt_flag_indirect_max_pio/F");
  fBDT->Branch("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy,"mgt_e_indirect_total_energy/F");
  fBDT->Branch("mgt_flag",&mgt_flag,"mgt_flag/F");

  fBDT->Branch("stw_1_energy",&stw_1_energy,"stw_1_energy/F");
  fBDT->Branch("stw_1_dis",&stw_1_dis,"stw_1_dis/F");
  fBDT->Branch("stw_1_dQ_dx",&stw_1_dQ_dx,"stw_1_dQ_dx/F");
  fBDT->Branch("stw_1_flag_single_shower",&stw_1_flag_single_shower,"stw_1_flag_single_shower/F");
  fBDT->Branch("stw_1_n_pi0",&stw_1_n_pi0,"stw_1_n_pi0/F");
  fBDT->Branch("stw_1_num_valid_tracks",&stw_1_num_valid_tracks,"stw_1_num_valid_tracks/F");
  fBDT->Branch("stw_1_flag",&stw_1_flag,"stw_1_flag/F");
  fBDT->Branch("stw_2_v_medium_dQ_dx", &stw_2_v_medium_dQ_dx);
  fBDT->Branch("stw_2_v_energy", &stw_2_v_energy);
  fBDT->Branch("stw_2_v_angle", &stw_2_v_angle);
  fBDT->Branch("stw_2_v_dir_length", &stw_2_v_dir_length);
  fBDT->Branch("stw_2_v_max_dQ_dx", &stw_2_v_max_dQ_dx);
  fBDT->Branch("stw_2_v_flag", &stw_2_v_flag);
  fBDT->Branch("stw_3_v_angle",&stw_3_v_angle);
  fBDT->Branch("stw_3_v_dir_length",&stw_3_v_dir_length);
  fBDT->Branch("stw_3_v_energy",&stw_3_v_energy);
  fBDT->Branch("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);
  fBDT->Branch("stw_3_v_flag",&stw_3_v_flag);
  fBDT->Branch("stw_4_v_angle",&stw_4_v_angle);
  fBDT->Branch("stw_4_v_dis",&stw_4_v_dis);
  fBDT->Branch("stw_4_v_energy",&stw_4_v_energy);
  fBDT->Branch("stw_4_v_flag",&stw_4_v_flag);
  fBDT->Branch("stw_flag", &stw_flag,"stw_flag/F");

  fBDT->Branch("spt_flag_single_shower", &spt_flag_single_shower, "spt_flag_single_shower/F");
  fBDT->Branch("spt_energy", &spt_energy, "spt_energy/F");
  fBDT->Branch("spt_shower_main_length", &spt_shower_main_length, "spt_shower_main_length/F");
  fBDT->Branch("spt_shower_total_length", &spt_shower_total_length, "spt_shower_total_length/F");
  fBDT->Branch("spt_angle_beam", &spt_angle_beam, "spt_angle_beam/F");
  fBDT->Branch("spt_angle_vertical", &spt_angle_vertical, "spt_angle_vertical/F");
  fBDT->Branch("spt_max_dQ_dx", &spt_max_dQ_dx, "spt_max_dQ_dx/F");
  fBDT->Branch("spt_angle_beam_1", &spt_angle_beam_1, "spt_angle_beam_1/F");
  fBDT->Branch("spt_angle_drift", &spt_angle_drift, "spt_angle_drift/F");
  fBDT->Branch("spt_angle_drift_1", &spt_angle_drift_1, "spt_angle_drift_1/F");
  fBDT->Branch("spt_num_valid_tracks", &spt_num_valid_tracks, "spt_num_valid_tracks/F");
  fBDT->Branch("spt_n_vtx_segs", &spt_n_vtx_segs, "spt_n_vtx_segs/F");
  fBDT->Branch("spt_max_length", &spt_max_length, "spt_max_length/F");
  fBDT->Branch("spt_flag", &spt_flag, "spt_flag/F");

  fBDT->Branch("stem_len_energy", &stem_len_energy, "stem_len_energy/F");
  fBDT->Branch("stem_len_length", &stem_len_length, "stem_len_length/F");
  fBDT->Branch("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check, "stem_len_flag_avoid_muon_check/F");
  fBDT->Branch("stem_len_num_daughters", &stem_len_num_daughters, "stem_len_num_daughters/F");
  fBDT->Branch("stem_len_daughter_length", &stem_len_daughter_length, "stem_len_daughter_length/F");
  fBDT->Branch("stem_len_flag", &stem_len_flag, "stem_len_flag/F");

  fBDT->Branch("lem_shower_total_length",&lem_shower_total_length,"lem_shower_total_length/F");
  fBDT->Branch("lem_shower_main_length",&lem_shower_main_length,"lem_shower_main_length/F");
  fBDT->Branch("lem_n_3seg",&lem_n_3seg,"lem_n_3seg/F");
  fBDT->Branch("lem_e_charge",&lem_e_charge,"lem_e_charge/F");
  fBDT->Branch("lem_e_dQdx",&lem_e_dQdx,"lem_e_dQdx/F");
  fBDT->Branch("lem_shower_num_segs",&lem_shower_num_segs,"lem_shower_num_segs/F");
  fBDT->Branch("lem_shower_num_main_segs",&lem_shower_num_main_segs,"lem_shower_num_main_segs/F");
  fBDT->Branch("lem_flag",&lem_flag,"lem_flag/F");

  fBDT->Branch("brm_n_mu_segs",&brm_n_mu_segs,"brm_n_mu_segs/F");
  fBDT->Branch("brm_Ep",&brm_Ep,"brm_Ep/F");
  fBDT->Branch("brm_energy",&brm_energy,"brm_energy/F");
  fBDT->Branch("brm_acc_length",&brm_acc_length,"brm_acc_length/F");
  fBDT->Branch("brm_shower_total_length",&brm_shower_total_length,"brm_shower_total_length/F");
  fBDT->Branch("brm_connected_length",&brm_connected_length,"brm_connected_length/F");
  fBDT->Branch("brm_n_size",&brm_n_size,"brm_n_size/F");
  fBDT->Branch("brm_acc_direct_length",&brm_acc_direct_length,"brm_acc_direct_length/F");
  fBDT->Branch("brm_n_shower_main_segs",&brm_n_shower_main_segs,"brm_n_shower_main_segs/F");
  fBDT->Branch("brm_n_mu_main",&brm_n_mu_main,"brm_n_mu_main/F");
  fBDT->Branch("brm_flag",&brm_flag,"brm_flag/F");

  fBDT->Branch("cme_mu_energy",&cme_mu_energy,"cme_mu_energy/F");
  fBDT->Branch("cme_energy",&cme_energy,"cme_energy/F");
  fBDT->Branch("cme_mu_length",&cme_mu_length,"cme_mu_length/F");
  fBDT->Branch("cme_length",&cme_length,"cme_length/F");
  fBDT->Branch("cme_angle_beam",&cme_angle_beam,"cme_angle_beam/F");
  fBDT->Branch("cme_flag",&cme_flag,"cme_flag/F");

  fBDT->Branch("anc_energy",&anc_energy,"anc_energy/F");
  fBDT->Branch("anc_angle",&anc_angle,"anc_angle/F");
  fBDT->Branch("anc_max_angle",&anc_max_angle,"anc_max_angle/F");
  fBDT->Branch("anc_max_length",&anc_max_length,"anc_max_length/F");
  fBDT->Branch("anc_acc_forward_length",&anc_acc_forward_length,"anc_acc_forward_length/F");
  fBDT->Branch("anc_acc_backward_length",&anc_acc_backward_length,"anc_acc_backward_length/F");
  fBDT->Branch("anc_acc_forward_length1",&anc_acc_forward_length1,"anc_acc_forward_length1/F");
  fBDT->Branch("anc_shower_main_length",&anc_shower_main_length,"anc_shower_main_length/F");
  fBDT->Branch("anc_shower_total_length",&anc_shower_total_length,"anc_shower_total_length/F");
  fBDT->Branch("anc_flag_main_outside",&anc_flag_main_outside,"anc_flag_main_outside/F");
  fBDT->Branch("anc_flag",&anc_flag,"anc_flag/F");

  fBDT->Branch("stem_dir_flag",&stem_dir_flag,"stem_dir_flag/F");
  fBDT->Branch("stem_dir_flag_single_shower",&stem_dir_flag_single_shower,"stem_dir_flag_single_shower/F");
  fBDT->Branch("stem_dir_filled",&stem_dir_filled,"stem_dir_filled/F");
  fBDT->Branch("stem_dir_angle",&stem_dir_angle,"stem_dir_angle/F");
  fBDT->Branch("stem_dir_energy",&stem_dir_energy,"stem_dir_energy/F");
  fBDT->Branch("stem_dir_angle1",&stem_dir_angle1,"stem_dir_angle1/F");
  fBDT->Branch("stem_dir_angle2",&stem_dir_angle2,"stem_dir_angle2/F");
  fBDT->Branch("stem_dir_angle3",&stem_dir_angle3,"stem_dir_angle3/F");
  fBDT->Branch("stem_dir_ratio",&stem_dir_ratio,"stem_dir_ratio/F");

  fBDT->Branch("vis_1_filled",&vis_1_filled,"vis_1_filled/F");
  fBDT->Branch("vis_1_n_vtx_segs",&vis_1_n_vtx_segs,"vis_1_n_vtx_segs/F");
  fBDT->Branch("vis_1_energy",&vis_1_energy,"vis_1_energy/F");
  fBDT->Branch("vis_1_num_good_tracks",&vis_1_num_good_tracks,"vis_1_num_good_tracks/F");
  fBDT->Branch("vis_1_max_angle",&vis_1_max_angle,"vis_1_max_angle/F");
  fBDT->Branch("vis_1_max_shower_angle",&vis_1_max_shower_angle,"vis_1_max_shower_angle/F");
  fBDT->Branch("vis_1_tmp_length1",&vis_1_tmp_length1,"vis_1_tmp_length1/F");
  fBDT->Branch("vis_1_tmp_length2",&vis_1_tmp_length2,"vis_1_tmp_length2/F");
  fBDT->Branch("vis_1_particle_type",&vis_1_particle_type,"vis_1_particle_type/F");
  fBDT->Branch("vis_1_flag",&vis_1_flag,"vis_1_flag/F");
  fBDT->Branch("vis_2_filled",&vis_2_filled,"vis_2_filled/F");
  fBDT->Branch("vis_2_n_vtx_segs",&vis_2_n_vtx_segs,"vis_2_n_vtx_segs/F");
  fBDT->Branch("vis_2_min_angle",&vis_2_min_angle,"vis_2_min_angle/F");
  fBDT->Branch("vis_2_min_weak_track",&vis_2_min_weak_track,"vis_2_min_weak_track/F");
  fBDT->Branch("vis_2_angle_beam",&vis_2_angle_beam,"vis_2_angle_beam/F");
  fBDT->Branch("vis_2_min_angle1",&vis_2_min_angle1,"vis_2_min_angle1/F");
  fBDT->Branch("vis_2_iso_angle1",&vis_2_iso_angle1,"vis_2_iso_angle1/F");
  fBDT->Branch("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx,"vis_2_min_medium_dQ_dx/F");
  fBDT->Branch("vis_2_min_length",&vis_2_min_length,"vis_2_min_length/F");
  fBDT->Branch("vis_2_sg_length",&vis_2_sg_length,"vis_2_sg_length/F");
  fBDT->Branch("vis_2_max_angle",&vis_2_max_angle,"vis_2_max_angle/F");
  fBDT->Branch("vis_2_max_weak_track",&vis_2_max_weak_track,"vis_2_max_weak_track/F");
  fBDT->Branch("vis_2_flag",&vis_2_flag,"vis_2_flag/F");
  fBDT->Branch("vis_flag",&vis_flag,"vis_flag/F");

  fBDT->Branch("br_filled",&br_filled,"br_filled/F");
  fBDT->Branch("br1_flag",&br1_flag,"br1_flag/F");
  fBDT->Branch("br1_1_flag",&br1_1_flag,"br1_1_flag/F");
  fBDT->Branch("br1_1_shower_type",&br1_1_shower_type,"br1_1_shower_type/F");
  fBDT->Branch("br1_1_vtx_n_segs",&br1_1_vtx_n_segs,"br1_1_vtx_n_segs/F");
  fBDT->Branch("br1_1_energy",&br1_1_energy,"br1_1_energy/F");
  fBDT->Branch("br1_1_n_segs",&br1_1_n_segs,"br1_1_n_segs/F");
  fBDT->Branch("br1_1_flag_sg_topology",&br1_1_flag_sg_topology,"br1_1_flag_sg_topology/F");
  fBDT->Branch("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory,"br1_1_flag_sg_trajectory/F");
  fBDT->Branch("br1_1_sg_length",&br1_1_sg_length,"br1_1_sg_length/F");
  fBDT->Branch("br1_2_flag",&br1_2_flag,"br1_2_flag/F");
  fBDT->Branch("br1_2_energy",&br1_2_energy,"br1_2_energy/F");
  fBDT->Branch("br1_2_n_connected",&br1_2_n_connected,"br1_2_n_connected/F");
  fBDT->Branch("br1_2_max_length",&br1_2_max_length,"br1_2_max_length/F");
  fBDT->Branch("br1_2_n_connected_1",&br1_2_n_connected_1,"br1_2_n_connected_1/F");
  fBDT->Branch("br1_2_vtx_n_segs",&br1_2_vtx_n_segs,"br1_2_vtx_n_segs/F");
  fBDT->Branch("br1_2_n_shower_segs",&br1_2_n_shower_segs,"br1_2_n_shower_segs/F");
  fBDT->Branch("br1_2_max_length_ratio",&br1_2_max_length_ratio,"br1_2_max_length_ratio/F");
  fBDT->Branch("br1_2_shower_length",&br1_2_shower_length,"br1_2_shower_length/F");
  fBDT->Branch("br1_3_flag",&br1_3_flag,"br1_3_flag/F");
  fBDT->Branch("br1_3_energy",&br1_3_energy,"br1_3_energy/F");
  fBDT->Branch("br1_3_n_connected_p",&br1_3_n_connected_p,"br1_3_n_connected_p/F");
  fBDT->Branch("br1_3_max_length_p",&br1_3_max_length_p,"br1_3_max_length_p/F");
  fBDT->Branch("br1_3_n_shower_segs",&br1_3_n_shower_segs,"br1_3_n_shower_segs/F");
  fBDT->Branch("br1_3_flag_sg_topology",&br1_3_flag_sg_topology,"br1_3_flag_sg_topology/F");
  fBDT->Branch("br1_3_flag_sg_trajectory",&br1_3_flag_sg_trajectory,"br1_3_flag_sg_trajectory/F");
  fBDT->Branch("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs,"br1_3_n_shower_main_segs/F");
  fBDT->Branch("br1_3_sg_length",&br1_3_sg_length,"br1_3_sg_length/F");

  //fBDT->Branch("br_filled",&br_filled,"br_filled/F");
  fBDT->Branch("br2_flag",&br2_flag,"br2_flag/F");
  fBDT->Branch("br2_flag_single_shower",&br2_flag_single_shower,"br2_flag_single_shower/F");
  fBDT->Branch("br2_num_valid_tracks",&br2_num_valid_tracks,"br2_num_valid_tracks/F");
  fBDT->Branch("br2_energy",&br2_energy,"br2_energy/F");
  fBDT->Branch("br2_angle1",&br2_angle1,"br2_angle1/F");
  fBDT->Branch("br2_angle2",&br2_angle2,"br2_angle2/F");
  fBDT->Branch("br2_angle",&br2_angle,"br2_angle/F");
  fBDT->Branch("br2_angle3",&br2_angle3,"br2_angle3/F");
  fBDT->Branch("br2_n_shower_main_segs",&br2_n_shower_main_segs,"br2_n_shower_main_segs/F");
  fBDT->Branch("br2_max_angle",&br2_max_angle,"br2_max_angle/F");
  fBDT->Branch("br2_sg_length",&br2_sg_length,"br2_sg_length/F");
  fBDT->Branch("br2_flag_sg_trajectory",&br2_flag_sg_trajectory,"br2_flag_sg_trajectory/F");

  //fBDT->Branch("br_filled",&br_filled,"br_filled/F");
  fBDT->Branch("br3_1_energy",&br3_1_energy,"br3_1_energy/F");
  fBDT->Branch("br3_1_n_shower_segments",&br3_1_n_shower_segments,"br3_1_n_shower_segments/F");
  fBDT->Branch("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory,"br3_1_sg_flag_trajectory/F");
  fBDT->Branch("br3_1_sg_direct_length",&br3_1_sg_direct_length,"br3_1_sg_direct_length/F");
  fBDT->Branch("br3_1_sg_length",&br3_1_sg_length,"br3_1_sg_length/F");
  fBDT->Branch("br3_1_total_main_length",&br3_1_total_main_length,"br3_1_total_main_length/F");
  fBDT->Branch("br3_1_total_length",&br3_1_total_length,"br3_1_total_length/F");
  fBDT->Branch("br3_1_iso_angle",&br3_1_iso_angle,"br3_1_iso_angle/F");
  fBDT->Branch("br3_1_sg_flag_topology",&br3_1_sg_flag_topology,"br3_1_sg_flag_topology/F");
  fBDT->Branch("br3_1_flag",&br3_1_flag,"br3_1_flag/F");
  fBDT->Branch("br3_2_n_ele",&br3_2_n_ele,"br3_2_n_ele/F");
  fBDT->Branch("br3_2_n_other",&br3_2_n_other,"br3_2_n_other/F");
  fBDT->Branch("br3_2_energy",&br3_2_energy,"br3_2_energy/F");
  fBDT->Branch("br3_2_total_main_length",&br3_2_total_main_length,"br3_2_total_main_length/F");
  fBDT->Branch("br3_2_total_length",&br3_2_total_length,"br3_2_total_length/F");
  fBDT->Branch("br3_2_other_fid",&br3_2_other_fid,"br3_2_other_fid/F");
  fBDT->Branch("br3_2_flag",&br3_2_flag,"br3_2_flag/F");
  fBDT->Branch("br3_3_v_energy",&br3_3_v_energy);
  fBDT->Branch("br3_3_v_angle",&br3_3_v_angle);
  fBDT->Branch("br3_3_v_dir_length",&br3_3_v_dir_length);
  fBDT->Branch("br3_3_v_length",&br3_3_v_length);
  fBDT->Branch("br3_3_v_flag",&br3_3_v_flag);
  fBDT->Branch("br3_4_acc_length", &br3_4_acc_length, "br3_4_acc_length/F");
  fBDT->Branch("br3_4_total_length", &br3_4_total_length, "br3_4_total_length/F");
  fBDT->Branch("br3_4_energy", &br3_4_energy, "br3_4_energy/F");
  fBDT->Branch("br3_4_flag", &br3_4_flag, "br3_4_flag/F");
  fBDT->Branch("br3_5_v_dir_length", &br3_5_v_dir_length);
  fBDT->Branch("br3_5_v_total_length", &br3_5_v_total_length);
  fBDT->Branch("br3_5_v_flag_avoid_muon_check", &br3_5_v_flag_avoid_muon_check);
  fBDT->Branch("br3_5_v_n_seg", &br3_5_v_n_seg);
  fBDT->Branch("br3_5_v_angle", &br3_5_v_angle);
  fBDT->Branch("br3_5_v_sg_length", &br3_5_v_sg_length);
  fBDT->Branch("br3_5_v_energy", &br3_5_v_energy);
  fBDT->Branch("br3_5_v_n_main_segs", &br3_5_v_n_main_segs);
  fBDT->Branch("br3_5_v_n_segs", &br3_5_v_n_segs);
  fBDT->Branch("br3_5_v_shower_main_length", &br3_5_v_shower_main_length);
  fBDT->Branch("br3_5_v_shower_total_length", &br3_5_v_shower_total_length);
  fBDT->Branch("br3_5_v_flag", &br3_5_v_flag);
  fBDT->Branch("br3_6_v_angle",&br3_6_v_angle);
  fBDT->Branch("br3_6_v_angle1",&br3_6_v_angle1);
  fBDT->Branch("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
  fBDT->Branch("br3_6_v_direct_length",&br3_6_v_direct_length);
  fBDT->Branch("br3_6_v_length",&br3_6_v_length);
  fBDT->Branch("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
  fBDT->Branch("br3_6_v_energy",&br3_6_v_energy);
  fBDT->Branch("br3_6_v_flag",&br3_6_v_flag);
  fBDT->Branch("br3_7_energy",&br3_7_energy,"br3_7_energy/F");
  fBDT->Branch("br3_7_min_angle",&br3_7_min_angle,"br3_7_min_angle/F");
  fBDT->Branch("br3_7_sg_length",&br3_7_sg_length,"br3_7_sg_length/F");
  fBDT->Branch("br3_7_main_length",&br3_7_shower_main_length,"br3_7_shower_main_length/F");
  fBDT->Branch("br3_7_flag",&br3_7_flag,"br3_7_flag/F");
  fBDT->Branch("br3_8_max_dQ_dx",&br3_8_max_dQ_dx,"br3_8_max_dQ_dx/F");
  fBDT->Branch("br3_8_energy",&br3_8_energy,"br3_8_energy/F");
  fBDT->Branch("br3_8_n_main_segs",&br3_8_n_main_segs,"br3_8_n_main_segs/F");
  fBDT->Branch("br3_8_shower_main_length",&br3_8_shower_main_length,"br3_8_shower_main_length/F");
  fBDT->Branch("br3_8_shower_length",&br3_8_shower_length,"br3_8_shower_length/F");
  fBDT->Branch("br3_8_flag",&br3_8_flag,"br3_8_flag/F");
  fBDT->Branch("br3_flag",&br3_flag,"br3_flag/F");

  //fBDT->Branch("br_filled",&br_filled,"br_filled/F");
  fBDT->Branch("br4_1_shower_main_length", &br4_1_shower_main_length, "br4_1_shower_main_length/F");
  fBDT->Branch("br4_1_shower_total_length", &br4_1_shower_total_length, "br4_1_shower_total_length/F");
  fBDT->Branch("br4_1_min_dis", &br4_1_min_dis, "br4_1_min_dis/F");
  fBDT->Branch("br4_1_energy", &br4_1_energy, "br4_1_energy/F");
  fBDT->Branch("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check, "br4_1_flag_avoid_muon_check/F");
  fBDT->Branch("br4_1_n_vtx_segs", &br4_1_n_vtx_segs, "br4_1_n_vtx_segs/F");
  fBDT->Branch("br4_1_n_main_segs", &br4_1_n_main_segs, "br4_1_n_main_segs/F");
  fBDT->Branch("br4_1_flag", &br4_1_flag, "br4_1_flag/F");
  fBDT->Branch("br4_2_ratio_45", &br4_2_ratio_45, "br4_2_ratio_45/F");
  fBDT->Branch("br4_2_ratio_35", &br4_2_ratio_35, "br4_2_ratio_35/F");
  fBDT->Branch("br4_2_ratio_25", &br4_2_ratio_25, "br4_2_ratio_25/F");
  fBDT->Branch("br4_2_ratio_15", &br4_2_ratio_15, "br4_2_ratio_15/F");
  fBDT->Branch("br4_2_energy",   &br4_2_energy, "br4_2_energy/F");
  fBDT->Branch("br4_2_ratio1_45", &br4_2_ratio1_45, "br4_2_ratio1_45/F");
  fBDT->Branch("br4_2_ratio1_35", &br4_2_ratio1_35, "br4_2_ratio1_35/F");
  fBDT->Branch("br4_2_ratio1_25", &br4_2_ratio1_25, "br4_2_ratio1_25/F");
  fBDT->Branch("br4_2_ratio1_15", &br4_2_ratio1_15, "br4_2_ratio1_15/F");
  fBDT->Branch("br4_2_iso_angle", &br4_2_iso_angle, "br4_2_iso_angle/F");
  fBDT->Branch("br4_2_iso_angle1", &br4_2_iso_angle1, "br4_2_iso_angle1/F");
  fBDT->Branch("br4_2_angle", &br4_2_angle, "br4_2_angle/F");
  fBDT->Branch("br4_2_flag", &br4_2_flag, "br4_2_flag/F");
  fBDT->Branch("br4_flag", &br4_flag, "br4_flag/F");

  fBDT->Branch("tro_1_v_particle_type",&tro_1_v_particle_type);
  fBDT->Branch("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
  fBDT->Branch("tro_1_v_min_dis",&tro_1_v_min_dis);
  fBDT->Branch("tro_1_v_sg1_length",&tro_1_v_sg1_length);
  fBDT->Branch("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
  fBDT->Branch("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
  fBDT->Branch("tro_1_v_tmp_length",&tro_1_v_tmp_length);
  fBDT->Branch("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
  fBDT->Branch("tro_1_v_dQ_dx_cut",&tro_1_v_dQ_dx_cut);
  fBDT->Branch("tro_1_v_flag_shower_topology",&tro_1_v_flag_shower_topology);
  fBDT->Branch("tro_1_v_flag",&tro_1_v_flag);
  fBDT->Branch("tro_2_v_energy",&tro_2_v_energy);
  fBDT->Branch("tro_2_v_stem_length",&tro_2_v_stem_length);
  fBDT->Branch("tro_2_v_iso_angle",&tro_2_v_iso_angle);
  fBDT->Branch("tro_2_v_max_length",&tro_2_v_max_length);
  fBDT->Branch("tro_2_v_angle",&tro_2_v_angle);
  fBDT->Branch("tro_2_v_flag",&tro_2_v_flag);
  fBDT->Branch("tro_3_stem_length",&tro_3_stem_length,"tro_3_stem_length/F");
  fBDT->Branch("tro_3_n_muon_segs",&tro_3_n_muon_segs,"tro_3_n_muon_segs/F");
  fBDT->Branch("tro_3_energy",&tro_3_energy,"tro_3_energy/F");
  fBDT->Branch("tro_3_flag",&tro_3_flag,"tro_3_flag/F");
  fBDT->Branch("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
  fBDT->Branch("tro_4_v_angle",&tro_4_v_angle);
  fBDT->Branch("tro_4_v_angle1",&tro_4_v_angle1);
  fBDT->Branch("tro_4_v_angle2",&tro_4_v_angle2);
  fBDT->Branch("tro_4_v_length",&tro_4_v_length);
  fBDT->Branch("tro_4_v_length1",&tro_4_v_length1);
  fBDT->Branch("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
  fBDT->Branch("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
  fBDT->Branch("tro_4_v_energy",&tro_4_v_energy);
  fBDT->Branch("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
  fBDT->Branch("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);
  fBDT->Branch("tro_4_v_flag",&tro_4_v_flag);
  fBDT->Branch("tro_5_v_max_angle",&tro_5_v_max_angle);
  fBDT->Branch("tro_5_v_min_angle",&tro_5_v_min_angle);
  fBDT->Branch("tro_5_v_max_length",&tro_5_v_max_length);
  fBDT->Branch("tro_5_v_iso_angle",&tro_5_v_iso_angle);
  fBDT->Branch("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
  fBDT->Branch("tro_5_v_min_count",&tro_5_v_min_count);
  fBDT->Branch("tro_5_v_max_count",&tro_5_v_max_count);
  fBDT->Branch("tro_5_v_energy",&tro_5_v_energy);
  fBDT->Branch("tro_5_v_flag",&tro_5_v_flag);
  fBDT->Branch("tro_flag",&tro_flag,"tro_flag/F");

  fBDT->Branch("hol_1_n_valid_tracks", &hol_1_n_valid_tracks,"hol_1_n_valid_tracks/F");
  fBDT->Branch("hol_1_min_angle", &hol_1_min_angle,"hol_1_min_angle/F");
  fBDT->Branch("hol_1_energy", &hol_1_energy,"hol_1_energy/F");
  fBDT->Branch("hol_1_flag_all_shower", &hol_1_flag_all_shower,"hol_1_flag_all_shower/F");
  fBDT->Branch("hol_1_min_length", &hol_1_min_length,"hol_1_min_length/F");
  fBDT->Branch("hol_1_flag", &hol_1_flag,"hol_1_flag/F");
  fBDT->Branch("hol_2_min_angle", &hol_2_min_angle,"hol_2_min_angle/F");
  fBDT->Branch("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx,"hol_2_medium_dQ_dx/F");
  fBDT->Branch("hol_2_ncount", &hol_2_ncount,"hol_2_ncount/F");
  fBDT->Branch("hol_2_energy", &hol_2_energy,"hol_2_energy/F");
  fBDT->Branch("hol_2_flag", &hol_2_flag,"hol_2_flag/F");
  fBDT->Branch("hol_flag", &hol_flag,"hol_flag/F");

  fBDT->Branch("lol_flag",&lol_flag,"lol_flag/F");
  fBDT->Branch("lol_1_v_energy",&lol_1_v_energy);
  fBDT->Branch("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
  fBDT->Branch("lol_1_v_nseg",&lol_1_v_nseg);
  fBDT->Branch("lol_1_v_angle",&lol_1_v_angle);
  fBDT->Branch("lol_1_v_flag",&lol_1_v_flag);
  fBDT->Branch("lol_2_v_length",&lol_2_v_length);
  fBDT->Branch("lol_2_v_angle",&lol_2_v_angle);
  fBDT->Branch("lol_2_v_type",&lol_2_v_type);
  fBDT->Branch("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
  fBDT->Branch("lol_2_v_energy",&lol_2_v_energy);
  fBDT->Branch("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
  fBDT->Branch("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);
  fBDT->Branch("lol_2_v_flag",&lol_2_v_flag);
  fBDT->Branch("lol_3_angle_beam",&lol_3_angle_beam,"lol_3_angle_beam/F");
  fBDT->Branch("lol_3_n_valid_tracks",&lol_3_n_valid_tracks,"lol_3_n_valid_tracks/F");
  fBDT->Branch("lol_3_min_angle",&lol_3_min_angle,"lol_3_min_angle/F");
  fBDT->Branch("lol_3_vtx_n_segs",&lol_3_vtx_n_segs,"lol_3_vtx_n_segs/F");
  fBDT->Branch("lol_3_energy",&lol_3_energy,"lol_3_energy/F");
  fBDT->Branch("lol_3_shower_main_length",&lol_3_shower_main_length,"lol_3_shower_main_length/F");
  fBDT->Branch("lol_3_n_out",&lol_3_n_out,"lol_3_n_out/F");
  fBDT->Branch("lol_3_n_sum",&lol_3_n_sum,"lol_3_n_sum/F");
  fBDT->Branch("lol_3_flag",&lol_3_flag,"lol_3_flag/F");

  fBDT->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");
  fBDT->Branch("cosmict_flag_2",&cosmict_flag_2,"cosmict_flag_2/F");
  fBDT->Branch("cosmict_flag_3",&cosmict_flag_3,"cosmict_flag_3/F");
  fBDT->Branch("cosmict_flag_4",&cosmict_flag_4,"cosmict_flag_4/F");
  fBDT->Branch("cosmict_flag_5",&cosmict_flag_5,"cosmict_flag_5/F");
  fBDT->Branch("cosmict_flag_6",&cosmict_flag_6,"cosmict_flag_6/F");
  fBDT->Branch("cosmict_flag_7",&cosmict_flag_7,"cosmict_flag_7/F");
  fBDT->Branch("cosmict_flag_8",&cosmict_flag_8,"cosmict_flag_8/F");
  fBDT->Branch("cosmict_flag_9",&cosmict_flag_9,"cosmict_flag_9/F");
  fBDT->Branch("cosmict_flag_10",&cosmict_flag_10);
  fBDT->Branch("cosmict_flag",&cosmict_flag,"cosmict_flag/F");
  fBDT->Branch("cosmict_2_filled",&cosmict_2_filled,"cosmict_2_filled/F");
  fBDT->Branch("cosmict_2_particle_type",&cosmict_2_particle_type,"cosmict_2_particle_type/F");
  fBDT->Branch("cosmict_2_n_muon_tracks",&cosmict_2_n_muon_tracks,"cosmict_2_n_muon_tracks/F");
  fBDT->Branch("cosmict_2_total_shower_length",&cosmict_2_total_shower_length,"cosmict_2_total_shower_length/F");
  fBDT->Branch("cosmict_2_flag_inside",&cosmict_2_flag_inside,"cosmict_2_flag_inside/F");
  fBDT->Branch("cosmict_2_angle_beam",&cosmict_2_angle_beam,"cosmict_2_angle_beam/F");
  fBDT->Branch("cosmict_2_flag_dir_weak",&cosmict_2_flag_dir_weak,"cosmict_2_flag_dir_weak/F");
  fBDT->Branch("cosmict_2_dQ_dx_end",&cosmict_2_dQ_dx_end,"cosmict_2_dQ_dx_end/F");
  fBDT->Branch("cosmict_2_dQ_dx_front",&cosmict_2_dQ_dx_front,"cosmict_2_dQ_dx_front/F");
  fBDT->Branch("cosmict_2_theta",&cosmict_2_theta,"cosmict_2_theta/F");
  fBDT->Branch("cosmict_2_phi",&cosmict_2_phi,"cosmict_2_phi/F");
  fBDT->Branch("cosmict_2_valid_tracks",&cosmict_2_valid_tracks,"cosmict_2_valid_tracks/F");
  fBDT->Branch("cosmict_3_filled",&cosmict_3_filled,"cosmict_3_filled/F");
  fBDT->Branch("cosmict_3_flag_inside",&cosmict_3_flag_inside,"cosmict_3_flag_inside/F");
  fBDT->Branch("cosmict_3_angle_beam",&cosmict_3_angle_beam,"cosmict_3_angle_beam/F");
  fBDT->Branch("cosmict_3_flag_dir_weak",&cosmict_3_flag_dir_weak,"cosmict_3_flag_dir_weak/F");
  fBDT->Branch("cosmict_3_dQ_dx_end",&cosmict_3_dQ_dx_end,"cosmict_3_dQ_dx_end/F");
  fBDT->Branch("cosmict_3_dQ_dx_front",&cosmict_3_dQ_dx_front,"cosmict_3_dQ_dx_front/F");
  fBDT->Branch("cosmict_3_theta",&cosmict_3_theta,"cosmict_3_theta/F");
  fBDT->Branch("cosmict_3_phi",&cosmict_3_phi,"cosmict_3_phi/F");
  fBDT->Branch("cosmict_3_valid_tracks",&cosmict_3_valid_tracks,"cosmict_3_valid_tracks/F");
  fBDT->Branch("cosmict_4_filled",&cosmict_4_filled,"cosmict_4_filled/F");
  fBDT->Branch("cosmict_4_flag_inside",&cosmict_4_flag_inside,"cosmict_4_flag_inside/F");
  fBDT->Branch("cosmict_4_angle_beam",&cosmict_4_angle_beam,"cosmict_4_angle_beam/F");
  fBDT->Branch("cosmict_4_connected_showers",&cosmict_4_connected_showers,"cosmict_4_connected_showers/F");
  fBDT->Branch("cosmict_5_filled",&cosmict_5_filled,"cosmict_5_filled/F");
  fBDT->Branch("cosmict_5_flag_inside",&cosmict_5_flag_inside,"cosmict_5_flag_inside/F");
  fBDT->Branch("cosmict_5_angle_beam",&cosmict_5_angle_beam,"cosmict_5_angle_beam/F");
  fBDT->Branch("cosmict_5_connected_showers",&cosmict_5_connected_showers,"cosmict_5_connected_showers/F");
  fBDT->Branch("cosmict_6_filled",&cosmict_6_filled,"cosmict_6_filled/F");
  fBDT->Branch("cosmict_6_flag_dir_weak",&cosmict_6_flag_dir_weak,"cosmict_6_flag_dir_weak/F");
  fBDT->Branch("cosmict_6_flag_inside",&cosmict_6_flag_inside,"cosmict_6_flag_inside/F");
  fBDT->Branch("cosmict_6_angle",&cosmict_6_angle,"cosmict_6_angle/F");
  fBDT->Branch("cosmict_7_filled",&cosmict_7_filled,"cosmict_7_filled/F");
  fBDT->Branch("cosmict_7_flag_sec",&cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
  fBDT->Branch("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
  fBDT->Branch("cosmict_7_total_shower_length",&cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
  fBDT->Branch("cosmict_7_flag_inside",&cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
  fBDT->Branch("cosmict_7_angle_beam",&cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
  fBDT->Branch("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
  fBDT->Branch("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
  fBDT->Branch("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
  fBDT->Branch("cosmict_7_theta",&cosmict_7_theta,"cosmict_7_theta/F");
  fBDT->Branch("cosmict_7_phi",&cosmict_7_phi,"cosmict_7_phi/F");
  fBDT->Branch("cosmict_8_filled",&cosmict_8_filled,"cosmict_8_filled/F");
  fBDT->Branch("cosmict_8_flag_out",&cosmict_8_flag_out,"cosmict_8_flag_out/F");
  fBDT->Branch("cosmict_8_muon_length",&cosmict_8_muon_length,"cosmict_8_muon_length/F");
  fBDT->Branch("cosmict_8_acc_length",&cosmict_8_acc_length,"cosmict_8_acc_length/F");
  fBDT->Branch("cosmict_10_flag_inside",&cosmict_10_flag_inside);
  fBDT->Branch("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  fBDT->Branch("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  fBDT->Branch("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  fBDT->Branch("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  fBDT->Branch("cosmict_10_length",&cosmict_10_length);

  fBDT->Branch("numu_cc_flag",&numu_cc_flag,"numu_cc_flag/F");
  fBDT->Branch("numu_cc_flag_1",&numu_cc_flag_1);
  fBDT->Branch("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  fBDT->Branch("numu_cc_1_length",&numu_cc_1_length);
  fBDT->Branch("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  fBDT->Branch("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  fBDT->Branch("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  fBDT->Branch("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  fBDT->Branch("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
  fBDT->Branch("numu_cc_flag_2",&numu_cc_flag_2);
  fBDT->Branch("numu_cc_2_length",&numu_cc_2_length);
  fBDT->Branch("numu_cc_2_total_length",&numu_cc_2_total_length);
  fBDT->Branch("numu_cc_2_n_daughter_tracks",&numu_cc_2_n_daughter_tracks);
  fBDT->Branch("numu_cc_2_n_daughter_all",&numu_cc_2_n_daughter_all);
  fBDT->Branch("numu_cc_flag_3",&numu_cc_flag_3,"numu_cc_flag_3/F");
  fBDT->Branch("numu_cc_3_particle_type",&numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
  fBDT->Branch("numu_cc_3_max_length",&numu_cc_3_max_length,"numu_cc_3_max_length/F");
  fBDT->Branch("numu_cc_3_track_length",&numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
  fBDT->Branch("numu_cc_3_max_length_all",&numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
  fBDT->Branch("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
  fBDT->Branch("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
  fBDT->Branch("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");

  fBDT->Branch("cosmict_2_4_score",&cosmict_2_4_score,"cosmict_2_4_score/F");
  fBDT->Branch("cosmict_3_5_score",&cosmict_3_5_score,"cosmict_3_5_score/F");
  fBDT->Branch("cosmict_6_score",&cosmict_6_score,"cosmict_6_score/F");
  fBDT->Branch("cosmict_7_score",&cosmict_7_score,"cosmict_7_score/F");
  fBDT->Branch("cosmict_8_score",&cosmict_8_score,"cosmict_8_score/F");
  fBDT->Branch("cosmict_10_score",&cosmict_10_score,"cosmict_10_score/F");
  fBDT->Branch("numu_1_score",&numu_1_score,"numu_1_score/F");
  fBDT->Branch("numu_2_score",&numu_2_score,"numu_2_score/F");
  fBDT->Branch("numu_3_score",&numu_3_score,"numu_3_score/F");
  fBDT->Branch("cosmict_score",&cosmict_score,"cosmict_score/F");
  fBDT->Branch("numu_score",&numu_score,"numu_score/F");
  fBDT->Branch("mipid_score",&mipid_score,"mipid_score/F");
  fBDT->Branch("gap_score",&gap_score,"gap_score/F");
  fBDT->Branch("hol_lol_score",&hol_lol_score,"hol_lol_score/F");
  fBDT->Branch("cme_anc_score",&cme_anc_score,"cme_anc_score/F");
  fBDT->Branch("mgo_mgt_score",&mgo_mgt_score,"mgo_mgt_score/F");
  fBDT->Branch("br1_score",&br1_score,"br1_score/F");
  fBDT->Branch("br3_score",&br3_score,"br3_score/F");
  fBDT->Branch("br3_3_score",&br3_3_score,"br3_3_score/F");
  fBDT->Branch("br3_5_score",&br3_5_score,"br3_5_score/F");
  fBDT->Branch("br3_6_score",&br3_6_score,"br3_6_score/F");
  fBDT->Branch("stemdir_br2_score",&stemdir_br2_score,"stemdir_br2_score/F");
  fBDT->Branch("trimuon_score",&trimuon_score,"trimuon_score/F");
  fBDT->Branch("br4_tro_score",&br4_tro_score,"br4_tro_score/F");
  fBDT->Branch("mipquality_score",&mipquality_score,"mipquality_score/F");
  fBDT->Branch("pio_1_score",&pio_1_score,"pio_1_score/F");
  fBDT->Branch("pio_2_score",&pio_2_score,"pio_2_score/F");
  fBDT->Branch("stw_spt_score",&stw_spt_score,"stw_spt_score/F");
  fBDT->Branch("vis_1_score",&vis_1_score,"vis_1_score/F");
  fBDT->Branch("vis_2_score",&vis_2_score,"vis_2_score/F");
  fBDT->Branch("stw_2_score",&stw_2_score,"stw_2_score/F");
  fBDT->Branch("stw_3_score",&stw_3_score,"stw_3_score/F");
  fBDT->Branch("stw_4_score",&stw_4_score,"stw_4_score/F");
  fBDT->Branch("sig_1_score",&sig_1_score,"sig_1_score/F");
  fBDT->Branch("sig_2_score",&sig_2_score,"sig_2_score/F");
  fBDT->Branch("lol_1_score",&lol_1_score,"lol_1_score/F");
  fBDT->Branch("lol_2_score",&lol_2_score,"lol_2_score/F");
  fBDT->Branch("tro_1_score",&tro_1_score,"tro_1_score/F");
  fBDT->Branch("tro_2_score",&tro_2_score,"tro_2_score/F");
  fBDT->Branch("tro_4_score",&tro_4_score,"tro_4_score/F");
  fBDT->Branch("tro_5_score",&tro_5_score,"tro_5_score/F");
  fBDT->Branch("nue_score",&nue_score,"nue_score/F");
  }

  fKINE = tfs->make<TTree>("T_KINEvars", "T_KINEvars");
  if(f_KINEvars){
  fKINE->Branch("kine_reco_Enu",&kine_reco_Enu,"kine_reco_Enu/F");
  fKINE->Branch("kine_reco_add_energy",&kine_reco_add_energy,"kine_reco_add_energy/F");
  fKINE->Branch("kine_energy_particle",&kine_energy_particle);
  fKINE->Branch("kine_energy_info",&kine_energy_info);
  fKINE->Branch("kine_particle_type",&kine_particle_type);
  fKINE->Branch("kine_energy_included",&kine_energy_included);
  fKINE->Branch("kine_pio_mass",&kine_pio_mass,"kine_pio_mass/F");
  fKINE->Branch("kine_pio_flag",&kine_pio_flag,"kine_pio_flag/I");
  fKINE->Branch("kine_pio_vtx_dis",&kine_pio_vtx_dis,"kine_pio_vtx_dis/F");
  fKINE->Branch("kine_pio_energy_1",&kine_pio_energy_1,"kine_pio_energy_1/F");
  fKINE->Branch("kine_pio_theta_1",&kine_pio_theta_1,"kine_pio_theta_1/F");
  fKINE->Branch("kine_pio_phi_1",&kine_pio_phi_1,"kine_pio_phi_1/F");
  fKINE->Branch("kine_pio_dis_1",&kine_pio_dis_1,"kine_pio_dis_1/F");
  fKINE->Branch("kine_pio_energy_2",&kine_pio_energy_2,"kine_pio_energy_2/F");
  fKINE->Branch("kine_pio_theta_2",&kine_pio_theta_2,"kine_pio_theta_2/F");
  fKINE->Branch("kine_pio_phi_2",&kine_pio_phi_2,"kine_pio_phi_2/F");
  fKINE->Branch("kine_pio_dis_2",&kine_pio_dis_2,"kine_pio_dis_2/F");
  fKINE->Branch("kine_pio_angle",&kine_pio_angle,"kine_pio_angle/F");
  }

}

void WireCellAnaTree::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // reset output at the beginning of each event
  resetOutput();

  // read NuMetrics
	std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;
	f_run = e.run();
 	f_subRun = e.subRun();
	f_event = e.id().event();

        auto triggerH = e.getHandle<std::vector<raw::Trigger>>("daq");
        if(! triggerH){
          std::cout << "WARNING: no raw::Trigger label: daq"  << std::endl;
          f_trigger_bits = 0;
          //return;
        }else{
          std::vector<art::Ptr<raw::Trigger> > t_v;
          art::fill_ptr_vector(t_v, triggerH);
          if(t_v.size()){
            f_trigger_bits = t_v[0]->TriggerBits();
          }
          else{
            f_trigger_bits = 0;
          }
       }

        auto const& containment_vec = e.getProduct<std::vector<nsm::NuSelectionContainment>>(fContainmentLabel);
	std::cout<<"--- NuSelectionContainment ---"<<std::endl;
	if(containment_vec.size()>1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	}
	if(containment_vec.size()<1) {
		f_flash_found = false;
		f_flash_found_asInt = -1;
		f_flash_time = -1;
		f_flash_measPe = -1;
		f_flash_predPe = -1;
		f_match_found = false;
		f_match_found_asInt = -1;
		f_match_type = -1;
		f_match_isFC = false;
		f_match_isTgm = false;
		f_match_notFC_FV = false;
		f_match_notFC_SP = false;
		f_match_notFC_DC = false;
		f_match_charge = -1;
		f_match_energy = -1;
                f_lm_cluster_length = -1;
                f_image_fail = false;
	}
        for(nsm::NuSelectionContainment const& c : containment_vec) {
                f_flash_found = c.GetFlashFound();
		f_flash_found_asInt = f_flash_found? 1:0;
                f_flash_time = c.GetFlashTime();
                f_flash_measPe = c.GetFlashMeasPe();
                f_flash_predPe = c.GetFlashPredPe();
                f_match_found = c.GetMatchFound();
		f_match_found_asInt = f_match_found? 1:0;
                f_match_type = c.GetMatchType();
                f_match_isFC = c.GetIsFC();
                f_match_isTgm = c.GetIsTGM();
                f_match_notFC_FV = c.GetNotFCFV();
                f_match_notFC_SP = c.GetNotFCSP();
                f_match_notFC_DC = c.GetNotFCDC();
                f_match_charge = c.GetCharge();
                f_match_energy = c.GetEnergy();
                f_lm_cluster_length = c.GetLength();
                f_image_fail = c.GetImageFail();
	}

        auto const& charge_vec = e.getProduct<std::vector<nsm::NuSelectionCharge>>(fChargeLabel);
	std::cout<<"--- NuSelectionCharge  ---"<<std::endl;
	if(charge_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	}
        for(nsm::NuSelectionCharge const& c : charge_vec) {
                f_match_chargeU = c.GetChargeU();
                f_match_chargeV = c.GetChargeV();
                f_match_chargeY = c.GetChargeY();
	}

        auto const& stm_vec = e.getProduct<std::vector<nsm::NuSelectionSTM>>(fSTMLabel);
	std::cout<<"--- NuSelectionSTM ---"<<std::endl;
	if(stm_vec.size()>1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	}
	if(stm_vec.size()<1) {
		f_stm_eventtype = -1;
		f_stm_lowenergy = -1;
		f_stm_LM = -1;
		f_stm_TGM = -1;
		f_stm_STM = -1;
		f_stm_FullDead = -1;
		f_stm_clusterlength = -1.0;
	}
        for(nsm::NuSelectionSTM const& s : stm_vec) {
                f_stm_eventtype = s.GetEventType();
                f_stm_lowenergy = s.GetLowEnergy();
                f_stm_LM = s.GetLM();
                f_stm_TGM = s.GetTGM();
                f_stm_STM = s.GetSTM();
                f_stm_FullDead = s.GetFullDead();
                f_stm_clusterlength = s.GetClusterLength();
	}

	if(fMC==true){

        auto const& truth_vec = e.getProduct<std::vector<nsm::NuSelectionTruth>>(fTruthLabel);
	std::cout<<"--- NuSelectionTruth  ---"<<std::endl;
	if(truth_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	}
        for(nsm::NuSelectionTruth const& t : truth_vec) {
                f_truth_nuEnergy = t.GetNuEnergy();
                f_truth_energyInside = t.GetEnergyInside();
                f_truth_electronInside = t.GetElectronInside();
                f_truth_nuPdg = t.GetNuPdg();
                f_truth_isCC = t.GetIsCC();
                f_truth_isEligible = t.GetIsEligible();
                f_truth_isFC = t.GetIsFC();
                f_truth_vtxInside = t.GetIsVtxInside();
                f_truth_vtxX = t.GetVtxX();
                f_truth_vtxY = t.GetVtxY();
                f_truth_vtxZ = t.GetVtxZ();
                f_truth_nuTime = t.GetTime();
	}

        auto const& match_vec = e.getProduct<std::vector<nsm::NuSelectionMatch>>(fMatchLabel);
	std::cout<<"--- NuSelectionMatch  ---"<<std::endl;
	if(match_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	}
        for(nsm::NuSelectionMatch const& m : match_vec) {
                f_match_completeness = m.GetCompleteness();
                f_match_completeness_energy = m.GetCompletenessEnergy();
                f_match_purity = m.GetPurity();
                f_match_purity_xz = m.GetPurityXZ();
                f_match_purity_xy = m.GetPurityXY();
	}

	/// save GENIE weights
	if(fSaveWeights) save_weights(e);
	if(fSaveLeeWeights) save_LEEweights(e);

	if ( !f_truth_isCC && f_truth_vtxInside && f_truth_nuPdg==14 ) f_NC_truth_isEligible = true;
	else f_NC_truth_isEligible = false;

	}

	if ( fMC==false ) f_truth_isEligible = false;
	if ( fMC==false ) f_NC_truth_isEligible = false;
	f_match_energyY = f_match_chargeY*23.6*1e-6/0.55;
	if( f_match_type&1U || (f_match_type>>1)&1U ) f_lightmismatch = true;
	else f_lightmismatch = false;

	/// PF validation starts
	// reco start [nested loop]
      if( f_wirecellPF ){
        auto particleHandle = e.getHandle<std::vector<simb::MCParticle>>(fPFInputTag);
        if (! particleHandle ) return;
        std::cout << "particles.size(): " << particleHandle->size() << std::endl;


        for (auto const& particle: *particleHandle){

        if(f_PFDump) {
            ++reco_Ntrack;
            reco_id[reco_Ntrack-1] = particle.TrackId();
            reco_pdg[reco_Ntrack-1] = particle.PdgCode();
            reco_process->push_back(particle.Process());
            reco_mother[reco_Ntrack-1] = particle.Mother();
            auto start_pos = particle.Position();
            reco_startXYZT[reco_Ntrack-1][0] = start_pos.X();
            reco_startXYZT[reco_Ntrack-1][1] = start_pos.Y();
            reco_startXYZT[reco_Ntrack-1][2] = start_pos.Z();
            reco_startXYZT[reco_Ntrack-1][3] = start_pos.T();
            auto end_pos = particle.EndPosition();
            reco_endXYZT[reco_Ntrack-1][0] = end_pos.X();
            reco_endXYZT[reco_Ntrack-1][1] = end_pos.Y();
            reco_endXYZT[reco_Ntrack-1][2] = end_pos.Z();
            reco_endXYZT[reco_Ntrack-1][3] = end_pos.T();
            auto start_mom = particle.Momentum();
            reco_startMomentum[reco_Ntrack-1][0] = start_mom.Px();
            reco_startMomentum[reco_Ntrack-1][1] = start_mom.Py();
            reco_startMomentum[reco_Ntrack-1][2] = start_mom.Pz();
            reco_startMomentum[reco_Ntrack-1][3] = start_mom.E();
            auto end_mom = particle.EndMomentum();
            reco_endMomentum[reco_Ntrack-1][0] = end_mom.Px();
            reco_endMomentum[reco_Ntrack-1][1] = end_mom.Py();
            reco_endMomentum[reco_Ntrack-1][2] = end_mom.Pz();
            reco_endMomentum[reco_Ntrack-1][3] = end_mom.E();
            reco_daughters->push_back(std::vector<int>());
            for (int i=0; i<particle.NumberDaughters(); ++i) {
                reco_daughters->back().push_back(particle.Daughter(i));
            }
        }

        if(f_savesps){
          //std::cout << "Num Trajectory Points: "<< particle->NumberTrajectoryPoints() << std::endl;
          for (uint i_pos=0; i_pos<particle.NumberTrajectoryPoints(); i_pos++){
              const TLorentzVector& pos = particle.Position(i_pos);
              auto mom = particle.Momentum(i_pos);
              f_sps_x->push_back(pos.X());
              f_sps_y->push_back(pos.Y());
              f_sps_z->push_back(pos.Z());
              f_sps_e->push_back(mom.E());
	      f_sps_pdg->push_back(particle.PdgCode());
	      f_sps_id->push_back(particle.TrackId());
            }
        }

		int trkID = particle.TrackId();
		fParticleMap[trkID] = particle;
		if(particle.Mother() == 0){
			if(fPrimaryID.size()<1){ // fill once
                        const TLorentzVector& position = particle.Position(0);
			f_reco_nuvtxX = position.X(); // units: cm inherit from larsoft
			f_reco_nuvtxY = position.Y(); // units: cm inherit from larsoft
			f_reco_nuvtxZ = position.Z(); // units: cm inherit from larsoft
                        f_neutrino_type = particle.StatusCode(); // neutrino type
			}
			fPrimaryID.push_back(trkID);
		}
		// TEST on mc_included information
		// For neutrino energy reconstruction
		// 0: for gamma or pi0, not included: their converted electrons will be included
		// 1: particles should be included
		// 3: low-energy gamma or distant activity < 80 cm to main cluster
		// 4: same as 3, but distance > 80 cm, not included

                /*std::cout<<"DEBUG -- mc_included information: "<<particle.Mass()<<std::endl;*/

		// not an actual mass reconstruction since PID tells us the mass if you believe
		// END
	}
	for (size_t i=0; i<fPrimaryID.size(); i++){
		//std::cout<<"Primary particle:  "<< fPrimaryID.at(i) <<std::endl;
		MuonID(fPrimaryID.at(i));
		ProtonID(fPrimaryID.at(i));
		ShowerID(fPrimaryID.at(i)); // nested loop to dump descendant (not one generation) shower IDs with additional constraint
	}
	//muon
	float MuonKE = 0; //temp shower kinetic energy
	for (size_t j=0; j<fMuonID.size(); j++){
		auto const& p = fParticleMap[fMuonID.at(j)];
		const TLorentzVector& pos = p.Position(0);
		const TLorentzVector& momentum = p.Momentum(0);
		float tempKE = momentum.E() - momentum.M();
		if( MuonKE>=tempKE ) continue;
		MuonKE = tempKE;
		f_reco_muonvtxX = pos.X(); // cm
		f_reco_muonvtxY = pos.Y();
		f_reco_muonvtxZ = pos.Z();
		f_reco_muonMomentum[0] = momentum.Px();
		f_reco_muonMomentum[1] = momentum.Py();
		f_reco_muonMomentum[2] = momentum.Pz();
		f_reco_muonMomentum[3] = momentum.E(); // GeV
	}
	// proton
	float ProtonKE = 0; //temp proton kinetic energy
	f_reco_Nproton = 0;
	for (size_t j=0; j<fProtonID.size(); j++){
		auto const& p = fParticleMap[fProtonID.at(j)];
		const TLorentzVector& pos = p.Position(0);
		const TLorentzVector& momentum = p.Momentum(0);
		float tempKE = momentum.E() - momentum.M();
		if(tempKE > 0.035) f_reco_Nproton ++; // 35 MeV threshold
		if( ProtonKE>=tempKE ) continue;
		ProtonKE = tempKE;
		f_reco_protonvtxX = pos.X(); // cm
		f_reco_protonvtxY = pos.Y();
		f_reco_protonvtxZ = pos.Z();
		f_reco_protonMomentum[0] = momentum.Px();
		f_reco_protonMomentum[1] = momentum.Py();
		f_reco_protonMomentum[2] = momentum.Pz();
		f_reco_protonMomentum[3] = momentum.E(); // GeV
	}

 	//shower
	float ShowerKE = 0; //temp shower kinetic energy
	bool ShowerPrimary = false; //temp shower primary particle?
	for (size_t j=0; j<fShowerID.size(); j++){
		//std::cout<<"Shower: "<< fShowerID.at(j) <<std::endl;
		auto const& p = fParticleMap[fShowerID.at(j)];
		bool IsPrimary = false;
		// find primary shower (electron): primary particle > others; high energy > low energyi
		const TLorentzVector& pos = p.Position(0);
		TVector3 vnu(f_reco_nuvtxX, f_reco_nuvtxY, f_reco_nuvtxZ);
		TVector3 vshower(pos.X(), pos.Y(), pos.Z());
		// here the "primary" may be from a pi0 but connected to reco nu vertex
		if ( p.Mother()==0 || (vnu-vshower).Mag()<1 ) IsPrimary = true;
		const TLorentzVector& momentum = p.Momentum(0);
		float tempKE = momentum.E() - momentum.M();
    f_reco_vec_showervtxX->push_back(pos.X()); // units: cm inherit from larsoft
		f_reco_vec_showervtxY->push_back(pos.Y()); // units: cm inherit from larsoft
		f_reco_vec_showervtxZ->push_back(pos.Z()); // units: cm inherit from larsoft
		f_reco_vec_showerKE->push_back(tempKE);
		if( (ShowerKE>=tempKE || ShowerPrimary) && !IsPrimary ) continue; // logic bug: ignore coincidence equal shower KE particle ordering
		if( IsPrimary && ShowerPrimary && ShowerKE>=tempKE) continue;
		if( IsPrimary ) ShowerPrimary = true; // a primary shower recorded
		ShowerKE = tempKE;
		f_reco_showervtxX = pos.X(); // units: cm inherit from larsoft
		f_reco_showervtxY = pos.Y(); // units: cm inherit from larsoft
		f_reco_showervtxZ = pos.Z(); // units: cm inherit from larsoft
		f_reco_showerKE = ShowerKE;
		f_reco_showerMomentum[0] = momentum.Px();
		f_reco_showerMomentum[1] = momentum.Py();
		f_reco_showerMomentum[2] = momentum.Pz();
		f_reco_showerMomentum[3] = momentum.E(); // GeV
	}
	//std::cout<<"Primary shower KE: "<< ShowerKE << std::endl;
	// reco end

	if(fMC == true){
	/// truth start
        auto particleHandle2 = e.getHandle<std::vector<simb::MCParticle>>(fPFtruthInputTag);
        if (! particleHandle2 ) return;
	// Get space charge correction
	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

        /**
         *  Construct track-id to MCParticle index mapping
         */
         // load MCShowers
         // auto const &mcs_h = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");

         // load MCParticles
         auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");

         // map mcp track-id to index
         std::map<int,int> MCPmap;
         // map mcp track-id to daughter track-ids
         // std::map<int,std::vector<int>> MCPdaughterMap;
         for (size_t p=0; p < mcp_h->size(); p++) {
           auto mcp = mcp_h->at(p);
           MCPmap[mcp.TrackId()] = p;
           // size_t ndaughters = mcp.NumberDaughters();
           // std::vector<int> daughter_v;
           // for (size_t d=0; d < ndaughters; d++) {
           //   daughter_v.push_back( mcp.Daughter(d) );
           // }// for all daughter particles
           // MCPdaughterMap[mcp.TrackId()] = daughter_v;
         }
         // end mapping
         // find two gammas
         // simb::MCParticle d0; // first daughter
         // simb::MCParticle d1; // second daughter
         // for (size_t i=0; i < mcs_h->size(); i++) {
         //   auto const& mcs = mcs_h->at(i);

         //   // get track-ID of shower
         //   auto MCSID = mcs.TrackID();
         //   // link this back to MCParticle from which it comes
         //   auto gammaMCP = mcp_h->at( MCPmap[MCSID] );
         //   // find daughters of gamma
         //   auto GAMMAdaughters = MCPdaughterMap[ MCSID ];
         //   //std::cout << "[Pi0TruthAnalysis] photon has " << GAMMAdaughters.size() << " daughter particles" << std::endl;
         //   if (GAMMAdaughters.size() == 2) {
         //     d0 = mcp_h->at( MCPmap[ GAMMAdaughters[0] ] );
         //     d1 = mcp_h->at( MCPmap[ GAMMAdaughters[1] ] );
         //     //std::cout << "[Pi0TruthAnalysis] \t gamma [ " << etot << " MeV ] -> electron (" << d0.PdgCode() << ") of [ " << d0.E(0) * 1000. << " MeV ] and "
         //     //        << " electron (" << d1.PdgCode() << ") of [ " << d1.E(0) * 1000. << " MeV ]" << std::endl;
         //   }

         // }
         //

         int true_photons = 0;
         bool truth_dalitz = false;

         f_truth_photon_angle = -1;
         f_truth_photon_dis = -1;
         TVector3 dir_1;
         TVector3 start_1;
         int mother_2 = -1;

         // Generator Info
         auto const& mclist = e.getProduct<std::vector<simb::MCTruth>>("generator");
         if (mclist.size()>0) {
             simb::MCTruth const& mctruth = mclist[0];
             if (mctruth.NeutrinoSet()) {
                 simb::MCNeutrino nu = mctruth.GetNeutrino();
                 mc_isnu = 1;
                 mc_nGeniePrimaries = mctruth.NParticles();
                 mc_nu_pdg = nu.Nu().PdgCode();
                 mc_nu_ccnc = nu.CCNC();
                 mc_nu_mode = nu.Mode();
                 mc_nu_intType = nu.InteractionType();
                 mc_nu_target = nu.Target();
                 mc_hitnuc = nu.HitNuc();
                 mc_hitquark = nu.HitQuark();
                 mc_nu_Q2 = nu.QSqr();
                 mc_nu_W = nu.W();
                 mc_nu_X = nu.X();
                 mc_nu_Y = nu.Y();
                 mc_nu_Pt = nu.Pt();
                 mc_nu_Theta = nu.Theta();
                 const TLorentzVector& position = nu.Nu().Position(0);
                 const TLorentzVector& momentum = nu.Nu().Momentum(0);
                 position.GetXYZT(mc_nu_pos);
                 momentum.GetXYZT(mc_nu_mom);
                 //cout << "nu: " << mc_nu_pdg << ", nPrim: " << mc_nGeniePrimaries
                 //     << ", ccnc: " << mc_nu_ccnc << endl;
                 //for (int i=0; i<mc_nGeniePrimaries; i++) {
                 //    simb::MCParticle particle = mctruth.GetParticle(i);
                 //    cout << "id: " << particle.TrackId()
                 //         << ", pdg: " << particle.PdgCode()
                 //         << endl;
                 //}
             }
         }

        for (auto const& particle: *particleHandle2){

        if(f_PFDump) {
            auto start_mom = particle.Momentum();
            //std::cout << particle.Mother()
            //<< ", " << start_mom.E()
            //<< ", " << f_PFDump_min_truth_energy
            //<< std::endl;
            if (!(particle.Mother()==0) && !(start_mom.E()>f_PFDump_min_truth_energy)) {
                continue;
            }

            ++truth_Ntrack;
            truth_id[truth_Ntrack-1] = particle.TrackId();
            truth_pdg[truth_Ntrack-1] = particle.PdgCode();
            truth_process->push_back(particle.Process());
            truth_mother[truth_Ntrack-1] = particle.Mother();
            truth_startMomentum[truth_Ntrack-1][0] = start_mom.Px();
            truth_startMomentum[truth_Ntrack-1][1] = start_mom.Py();
            truth_startMomentum[truth_Ntrack-1][2] = start_mom.Pz();
            truth_startMomentum[truth_Ntrack-1][3] = start_mom.E();
            auto end_mom = particle.EndMomentum();
            truth_endMomentum[truth_Ntrack-1][0] = end_mom.Px();
            truth_endMomentum[truth_Ntrack-1][1] = end_mom.Py();
            truth_endMomentum[truth_Ntrack-1][2] = end_mom.Pz();
            truth_endMomentum[truth_Ntrack-1][3] = end_mom.E();
            auto start_pos = particle.Position();
            auto start_sce_offset = SCE->GetPosOffsets(geo::Point_t(start_pos.X(), start_pos.Y(), start_pos.Z()));
            truth_startXYZT[truth_Ntrack-1][0] = start_pos.X() - start_sce_offset.X();
            truth_startXYZT[truth_Ntrack-1][1] = start_pos.Y() + start_sce_offset.Y();
            truth_startXYZT[truth_Ntrack-1][2] = start_pos.Z() + start_sce_offset.Z();
            truth_startXYZT[truth_Ntrack-1][3] = start_pos.T();
            truth_startXYZT[truth_Ntrack-1][0] = (truth_startXYZT[truth_Ntrack-1][0] + 0.6)*1.101/1.098 + start_pos.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
            auto end_pos = particle.EndPosition();
            auto end_sce_offset = SCE->GetPosOffsets(geo::Point_t(end_pos.X(), end_pos.Y(), end_pos.Z()));
            truth_endXYZT[truth_Ntrack-1][0] = end_pos.X() - end_sce_offset.X();
            truth_endXYZT[truth_Ntrack-1][1] = end_pos.Y() + end_sce_offset.Y();
            truth_endXYZT[truth_Ntrack-1][2] = end_pos.Z() + end_sce_offset.Z();
            truth_endXYZT[truth_Ntrack-1][3] = end_pos.T();
            truth_endXYZT[truth_Ntrack-1][0] = (truth_endXYZT[truth_Ntrack-1][0] + 0.6)*1.101/1.098 + end_pos.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
            truth_daughters->push_back(std::vector<int>());
            for (int i=0; i<particle.NumberDaughters(); ++i) {
                truth_daughters->back().push_back(particle.Daughter(i));
            }
            if(f_save_track_position) {
                size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
                TClonesArray *Lposition = new TClonesArray("TLorentzVector", numberTrajectoryPoints);
                // Read the position and momentum along this particle track
                for(unsigned int j=0; j<numberTrajectoryPoints; j++) {
                    new ((*Lposition)[j]) TLorentzVector(particle.Position(j));
                }
                fMC_trackPosition->Add(Lposition);
            }
        }

                if( particle.Mother() == 0 && (particle.PdgCode() == 11 || particle.PdgCode() == -11) ){
      f_truth_showerPdg = particle.PdgCode();
                        const TLorentzVector& position = particle.Position(0);
			//f_truth_corr_showervtxX = position.X(); // units: cm inherit from larsoft
			//f_truth_corr_showervtxY = position.Y(); // units: cm inherit from larsoft
			//f_truth_corr_showervtxZ = position.Z(); // units: cm inherit from larsoft
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_showervtxX = position.X() - sce_offset.X();
			f_truth_corr_showervtxY = position.Y() + sce_offset.Y();
			f_truth_corr_showervtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_showervtxX = (f_truth_corr_showervtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
			std::cout<<"Shower info: "<<position.X() <<", "<<position.Y() <<", "<<position.Z()<<", "<<position.T()<<" ns"<<std::endl;
			std::cout<<"Shower vertex SCE offset: "<<sce_offset.X() <<", "<<sce_offset.Y() <<", "<<sce_offset.Z()<<std::endl;
                        const TLorentzVector& showerMom = particle.Momentum(0);
			f_truth_showerKE = showerMom.E() - showerMom.M();
			f_truth_showerMomentum[0] = showerMom.Px();
			f_truth_showerMomentum[1] = showerMom.Py();
			f_truth_showerMomentum[2] = showerMom.Pz();
			f_truth_showerMomentum[3] = showerMom.E();

		}
                if( particle.Mother()==0 && (particle.PdgCode() == 13 || particle.PdgCode() == -13) ){
                        const TLorentzVector& position = particle.Position(0);
			f_truth_muonvtxX = position.X();
			f_truth_muonvtxY = position.Y();
			f_truth_muonvtxZ = position.Z();
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_muonvtxX = position.X() - sce_offset.X();
			f_truth_corr_muonvtxY = position.Y() + sce_offset.Y();
			f_truth_corr_muonvtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_muonvtxX = (f_truth_corr_muonvtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us

                        const TLorentzVector& endposition = particle.EndPosition();
			f_truth_muonendX = endposition.X();
			f_truth_muonendY = endposition.Y();
			f_truth_muonendZ = endposition.Z();

                        const TLorentzVector& momentum = particle.Momentum(0);
			f_truth_muonMomentum[0] = momentum.Px();
			f_truth_muonMomentum[1] = momentum.Py();
			f_truth_muonMomentum[2] = momentum.Pz();
			f_truth_muonMomentum[3] = momentum.E();
		}

    if( particle.PdgCode() == 22 && particle.Process()!="eBrem" && particle.Process()!="annihil" && particle.StatusCode() == 1 &&
        !(f_truth_isCC==1 && f_truth_nuPdg==12)){
      f_truth_showerPdg = particle.PdgCode();
                        const TLorentzVector& position = particle.Position(0);
			//f_truth_corr_showervtxX = position.X(); // units: cm inherit from larsoft
			//f_truth_corr_showervtxY = position.Y(); // units: cm inherit from larsoft
			//f_truth_corr_showervtxZ = position.Z(); // units: cm inherit from larsoft
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_showervtxX = position.X() - sce_offset.X();
			f_truth_corr_showervtxY = position.Y() + sce_offset.Y();
			f_truth_corr_showervtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_showervtxX = (f_truth_corr_showervtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
			std::cout<<"Shower info: "<<position.X() <<", "<<position.Y() <<", "<<position.Z()<<", "<<position.T()<<" ns"<<std::endl;
			std::cout<<"Shower vertex SCE offset: "<<sce_offset.X() <<", "<<sce_offset.Y() <<", "<<sce_offset.Z()<<std::endl;
                        const TLorentzVector& showerMom = particle.Momentum(0);
			f_truth_showerKE = showerMom.E() - showerMom.M();
			f_truth_showerMomentum[0] = showerMom.Px();
			f_truth_showerMomentum[1] = showerMom.Py();
			f_truth_showerMomentum[2] = showerMom.Pz();
			f_truth_showerMomentum[3] = showerMom.E();

      float ex = particle.EndX();
      float ey = particle.EndY();
      float ez = particle.EndZ();
      auto end_sce_offset = SCE->GetPosOffsets(geo::Point_t(ex, ey, ez));
      ex -= end_sce_offset.X();
      ey += end_sce_offset.Y();
      ez += end_sce_offset.Z();
      if (f_truth_showerKE  > 0.02 &&
          /*f_truth_corr_showervtxX > 3.0 && f_truth_corr_showervtxX < 253.0 &&
          f_truth_corr_showervtxY > -113.0 && f_truth_corr_showervtxY < 114.0 &&
          f_truth_corr_showervtxZ > 3.0 && f_truth_corr_showervtxZ < 1034.0 &&*/
          ex > 3.0 && ex < 253.0 && ey > -113.0 && ey < 114.0 && ez > 3.0 && ez < 1034.0){
            true_photons+=1;
            if (true_photons==1){
              TVector3 dir(showerMom.Px(),showerMom.Py(),showerMom.Pz());
              TVector3 dis(ex,ey,ez);
              dir_1 = dir;
              start_1 = dis;
              f_truth_photon_dis = start_1.Mag();
              mother_2 = particle.Mother();
            }

            if (true_photons==2 && particle.Mother()==mother_2){
              TVector3 dir_2(showerMom.Px(),showerMom.Py(),showerMom.Pz());
              TVector3 start_2(ex,ey,ez);
              f_truth_photon_angle = (dir_1.Angle(dir_2))*(180./TMath::Pi());
              f_truth_photon_dis = (start_1-start_2).Mag();
            }

            if (particle.Mother()!=-1){
              simb::MCParticle mother = mcp_h->at( MCPmap[particle.Mother()] );
              if (mother.PdgCode()==22){
                if (mother.Mother()!=-1.){
                  simb::MCParticle mothermother = mcp_h->at( MCPmap[mother.Mother()] );
                  f_truth_showerMother = mothermother.PdgCode();
                }
              }else{
                f_truth_showerMother = mother.PdgCode();
              }
            }
      }

		}



                int ndaughters = particle.NumberDaughters();
                if( particle.Mother()==0 && particle.PdgCode() == 111 && ndaughters == 2) {

                        int gamma_0_trkId = particle.Daughter(0);
                        int gamma_1_trkId = particle.Daughter(1);
                        simb::MCParticle d0 = mcp_h->at( MCPmap[gamma_0_trkId] ); // 1st gamma
                        simb::MCParticle d1 = mcp_h->at( MCPmap[gamma_1_trkId] ); // 2nd gamma
                        f_truth_NprimPio ++;
                        f_truth_pio_energy_1 = d0.E(0);
                        f_truth_pio_energy_2 = d1.E(0);
                        auto dir0 = d0.Momentum().Vect().Unit();
                        auto dir1 = d1.Momentum().Vect().Unit();
                        f_truth_pio_angle = dir0.Angle(dir1) / 3.1416 * 180;
                        std::cout << "mass: " << particle.Mass() << " , gamma energys: " << d0.E(0) << " " << d1.E(0) << " , angle: " << f_truth_pio_angle << std::endl;
                }

      if ( particle.PdgCode() == 111 ) {
        f_truth_Npi0++;
      for (int i_d=0; i_d<ndaughters; i_d++){
        simb::MCParticle daughter = mcp_h->at( MCPmap[particle.Daughter(i_d)] );
        if (ndaughters>2 && abs(daughter.PdgCode())==11){ truth_dalitz = true;}
      }
  	}

	}
	auto sce_offset = SCE->GetPosOffsets(geo::Point_t(f_truth_vtxX, f_truth_vtxY, f_truth_vtxZ));
	f_truth_corr_nuvtxX = f_truth_vtxX - sce_offset.X();
	f_truth_corr_nuvtxY = f_truth_vtxY + sce_offset.Y();
	f_truth_corr_nuvtxZ = f_truth_vtxZ + sce_offset.Z();
	f_truth_corr_nuvtxX = (f_truth_corr_nuvtxX + 0.6)*1.101/1.098 + f_truth_nuTime*1.101*0.1; //nuTime: us; 1.101 mm/us
	std::cout<<"Neutrino info: "<<f_truth_vtxX <<", "<<f_truth_vtxY <<", "<<f_truth_vtxZ<<", "<<f_truth_nuTime<<" us"<<std::endl;
	std::cout<<"Neutrino vertex SCE offset: "<<sce_offset.X() <<", "<<sce_offset.Y() <<", "<<sce_offset.Z()<<std::endl;

	// neutrino interaction type. Integer, see MCNeutrino.h for more details.
	//art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
	//e.getByLabel("generator",mctruthListHandle);
	//std::vector<art::Ptr<simb::MCTruth> > mclist;
	//art::fill_ptr_vector(mclist, mctruthListHandle);
	//art::Ptr<simb::MCTruth> mctruth;

	if (mclist.size()>0) {
                simb::MCTruth const& mctruth = mclist[0];
                if (mctruth.NeutrinoSet()) {
                        simb::MCNeutrino nu = mctruth.GetNeutrino();
			f_truth_nuIntType = nu.InteractionType();
			// one can access more neutrino GENIE info

			// NC Delta radiate events
			// Delta (top 2 BRs: Delta+ and Delta0)
			// Type 1: Delta --> Gamma
			// Type 2: Delta --> Gamma (internal) --> Gamma (final state)
			f_truth_NCDelta = 0;
			for(int k=0; k<mctruth.NParticles(); k++){
				simb::MCParticle const & mcp = mctruth.GetParticle(k);
        if(mcp.PdgCode()==22 && (mcp.StatusCode()==1 || mcp.StatusCode()==14)){
					const simb::MCParticle mother = mctruth.GetParticle(mcp.Mother());
          f_truth_showerMother = mother.PdgCode();
          if (mother.PdgCode()==22){
            const simb::MCParticle mothermother = mctruth.GetParticle(mother.Mother());
            f_truth_showerMother = mothermother.PdgCode();
          }
					//if(abs(mother.PdgCode()) == 2114 || abs(mother.PdgCode()) == 2214 || abs(mother.PdgCode()) == 2224 || abs(mother.PdgCode()) == 1114) {
          if(abs(f_truth_showerMother) == 2114 || abs(f_truth_showerMother) == 2214 || abs(f_truth_showerMother) == 2224 || abs(f_truth_showerMother) == 1114) {
            f_truth_NCDelta = 1;
					  break;
					}
				}
			}


                        const TLorentzVector& position = nu.Nu().Position(0);
                        const TLorentzVector& momentum = nu.Nu().Momentum(0);
                        f_truth_nu_pos[0] = position.X(); // cm
                        f_truth_nu_pos[1] = position.Y();
                        f_truth_nu_pos[2] = position.Z();
                        f_truth_nu_pos[3] = position.T();
                        f_truth_nu_momentum[0] = momentum.Px();
                        f_truth_nu_momentum[1] = momentum.Py();
                        f_truth_nu_momentum[2] = momentum.Pz();
                        f_truth_nu_momentum[3] = momentum.E(); // GeV

		}

	}

  f_truth_single_photon = 0;
  //f_truth_photon_angle = true_angle;
  if (true_photons==1 && !truth_dalitz){
    f_truth_single_photon = 1;
  }else if (true_photons==2 && !truth_dalitz && f_truth_photon_angle>-1 && f_truth_photon_angle<20){
    f_truth_single_photon = 1;
  }

	// neutrino scattering code. Integer, see GTruth.h for more details.
        if(auto gtruthListHandle = e.getHandle<std::vector<simb::GTruth>>("generator")){
          if (gtruthListHandle->size()>0) {
                f_truth_nuScatType = gtruthListHandle->front().fGscatter;
	  }
        }
        // by default, f_truth_nuScatType = -1

	// flux truth, see MCFlux.h for more details
	// and the description at http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/v19/output_gnumi.html
        if (auto mcfluxListHandle = e.getHandle<std::vector<simb::MCFlux>>("generator")){
          if (mcfluxListHandle->size()>0) {
                simb::MCFlux const& mcflux = mcfluxListHandle->front();
                f_mcflux_run = mcflux.frun;
                f_mcflux_evtno = mcflux.fevtno;
                f_mcflux_ndecay = mcflux.fndecay;
                f_mcflux_ntype = mcflux.fntype;
                f_mcflux_ptype = mcflux.fptype;
                f_mcflux_tptype = mcflux.ftptype;
                f_mcflux_nuEnergy = mcflux.fnenergyn; // neutrino energy for a decay
                f_mcflux_vx = mcflux.fvx; // vertex of hadron decay
                f_mcflux_vy = mcflux.fvy;
                f_mcflux_vz = mcflux.fvz;
                f_mcflux_genx = mcflux.fgenx;
                f_mcflux_geny = mcflux.fgeny;
                f_mcflux_genz = mcflux.fgenz;
                f_mcflux_dk2gen = mcflux.fdk2gen; // distance from decay to ray origin
                f_mcflux_gen2vtx = mcflux.fgen2vtx; // distance from ray origin to event vtx
	  }
        }

	}
	//

        if(f_get_redk2nu_time){
          double spill_time = TimeOffset();
	  double propegation_time = (f_mcflux_dk2gen+f_mcflux_gen2vtx)*100*0.033356;
	  art::Handle<std::vector<bsim::Dk2Nu>> dk2nu_flux_handle;
	  std::vector<art::Ptr<bsim::Dk2Nu> > dk2nu_flux;
          e.getByLabel("generator",dk2nu_flux_handle);
          art::fill_ptr_vector(dk2nu_flux,dk2nu_flux_handle);
	  if(dk2nu_flux.size()==0){
            std::cout<<"no redk2nu found"<<std::endl;
	  }
	  else{
            art::Ptr<bsim::Dk2Nu> this_dk2nu_flux = dk2nu_flux.at(0);
	    std::vector<bsim::Ancestor> ancestors = this_dk2nu_flux->ancestor;
            double decay_time = ancestors.back().startt;
            f_redk2nu_time = propegation_time + decay_time + spill_time;
	    f_redk2nu_time_nospill = propegation_time + decay_time;
	    f_redk2nu_deltatime = f_redk2nu_time - f_truth_nu_pos[3];
            std::cout<<"redk2nu: "<<"og time "<<f_truth_nu_pos[3]<<"  redk2nu_time "<<f_redk2nu_time<<"  f_redk2nu_deltatime "<<f_redk2nu_deltatime<<"  spill_time "<<spill_time<<" decay_time "<<decay_time<<" propegation_time "<<propegation_time<<"  time to window "<<decay_time+(f_mcflux_dk2gen)*100*0.033356<<"  time to vtx "<<decay_time+(f_mcflux_dk2gen+f_mcflux_gen2vtx)*100*0.033356<<"  ancestors.front().startt "<<ancestors.front().startt<<"  seed "<<fRndmGen->GetSeed()<<std::endl;
          }
  	}
	
	/// truth end
	std::cout<<"Corrected Truth Neutrino vertex: ("<<f_truth_corr_nuvtxX<<", "<<f_truth_corr_nuvtxY<<", "<<f_truth_corr_nuvtxZ<<")"<<"\n";
	std::cout<<"Corrected Truth Shower vertex: ("<<f_truth_corr_showervtxX<<", "<<f_truth_corr_showervtxY<<", "<<f_truth_corr_showervtxZ<<")"<<"\n";
	std::cout<<"Reco neutrino vertex: ("<<f_reco_nuvtxX<<", "<<f_reco_nuvtxY<<", "<<f_reco_nuvtxZ<<")"<<"\n";
	std::cout<<"Reco shower vertex: ("<<f_reco_showervtxX<<", "<<f_reco_showervtxY<<", "<<f_reco_showervtxZ<<")"<<"\n";
	std::cout<<"Reco/True shower KE: "<<f_reco_showerKE<<" / "<<f_truth_showerKE<<"\n";
	/// PF validation ends

	/// vertex distance [diff]
	TVector3 vr(f_reco_nuvtxX, f_reco_nuvtxY, f_reco_nuvtxZ);
	TVector3 vt(f_truth_corr_nuvtxX, f_truth_corr_nuvtxY, f_truth_corr_nuvtxZ);
	TVector3 vrshower(f_reco_showervtxX, f_reco_showervtxY, f_reco_showervtxZ);
	TVector3 vtshower(f_truth_corr_showervtxX, f_truth_corr_showervtxY, f_truth_corr_showervtxZ);
	TVector3 vrmuon(f_reco_muonvtxX, f_reco_muonvtxY, f_reco_muonvtxZ);
	TVector3 vtmuon(f_truth_corr_muonvtxX, f_truth_corr_muonvtxY, f_truth_corr_muonvtxZ);
	f_nuvtx_diff = (vr-vt).Mag();
	f_showervtx_diff = (vrshower-vtshower).Mag();
	f_muonvtx_diff = (vrmuon-vtmuon).Mag();
      }

      /// BDT input variables
      if(f_BDTvars){
        auto bdthandle = e.getHandle<std::vector<nsm::NuSelectionBDT>>(fPFInputTag);
        if (! bdthandle) return;
	std::cout<<"--- NuSelectionBDT ---"<<std::endl;
        if(bdthandle->size()>1) {
		std::cout<<"WARNING: >1 set of BDT input variables" << std::endl;
		return;
	}
        for(nsm::NuSelectionBDT const& bdt : *bdthandle) {
                //Always fill, we get nue files when we pass KDAR gen sel 
                ReadSSMBDTvar(bdt);
                //Only fill the rest if we passed the regular generic nu selection
                if(f_lm_cluster_length>10){
                  ReadBDTvar(bdt);
                }
	std::cout<<"BDT input vars check: \n"<<
	"Cosmic Tagger: "<<
        bdt.GetCosmicTagger().cosmic_filled<<" "<<
        bdt.GetCosmicTagger().cosmic_flag<<" "<<
        bdt.GetCosmicTagger().cosmic_n_solid_tracks<<" "<<
        bdt.GetCosmicTagger().cosmic_energy_main_showers<<" "<<
        bdt.GetCosmicTagger().cosmic_energy_direct_showers<<" "<<
        bdt.GetCosmicTagger().cosmic_energy_indirect_showers<<" "<<
        bdt.GetCosmicTagger().cosmic_n_direct_showers<<" "<<
        bdt.GetCosmicTagger().cosmic_n_indirect_showers<<" "<<
        bdt.GetCosmicTagger().cosmic_n_main_showers<<"\n";
	}
      }

      if(f_KINEvars){
        auto kinehandle = e.getHandle<std::vector<nsm::NuSelectionKINE>>(fPFInputTag);
        if (! kinehandle) return;
	std::cout<<"--- NuSelectionKINE ---"<<std::endl;
        if(kinehandle->size()>1) {
		std::cout<<"WARNING: >1 set of KINE input variables" << std::endl;
		return;
	}
        for(nsm::NuSelectionKINE const& kine : *kinehandle) {
                if(f_lm_cluster_length>10){
                  ReadKINEvar(kine);
                }
	std::cout<<"KINE input vars check: \n"<<

	kine.GetKineInfo().kine_reco_Enu<<" "<<
	kine.GetKineInfo().kine_reco_add_energy<<" ";
        if(kine.GetKineInfo().kine_energy_particle->empty()){
          std::cout<<"0"<<" "<<
          "0"<<" "<<
          "0"<<" "<<
          "0"<<" ";
        }
        else{
	  std::cout<<kine.GetKineInfo().kine_energy_particle->at(0)<<" "<<
	  kine.GetKineInfo().kine_energy_info->at(0)<<" "<<
	  kine.GetKineInfo().kine_particle_type->at(0)<<" "<<
	  kine.GetKineInfo().kine_energy_included->at(0)<<" ";
        }
	std::cout<<kine.GetKineInfo().kine_pio_mass<<" "<<
	kine.GetKineInfo().kine_pio_flag<<" "<<
	kine.GetKineInfo().kine_pio_vtx_dis<<" "<<
	kine.GetKineInfo().kine_pio_energy_1<<" "<<
	kine.GetKineInfo().kine_pio_theta_1<<" "<<
	kine.GetKineInfo().kine_pio_phi_1<<" "<<
	kine.GetKineInfo().kine_pio_dis_1<<" "<<
	kine.GetKineInfo().kine_pio_energy_2<<" "<<
	kine.GetKineInfo().kine_pio_theta_2<<" "<<
	kine.GetKineInfo().kine_pio_phi_2<<" "<<
	kine.GetKineInfo().kine_pio_dis_2<<" "<<
	kine.GetKineInfo().kine_pio_angle<<"\n";
	}
      }
  //if ( (!fMC || fIsNuMI) && kine_reco_Enu>-1.) {nsbeamtiming(e);}
  //if ( (!fMC || fIsNuMI) && numu_cc_flag!=-1) {nsbeamtiming(e);}
  if ( (!fMC || fIsNuMI) && (numu_cc_flag!=-1 || (ssm_numu_cc_flag!=-1 && f_ssmBDT) ) ) {nsbeamtiming(e);}
        fTreeEval->Fill();
	fPFeval->Fill();
	fBDT->Fill();
	fKINE->Fill();
}

void WireCellAnaTree::endSubRun(art::SubRun const& sr)
{
  // Implementation of optional member function here.

  // POT counting
  if( fPOT_counting==true ){
	art::Handle<sumdata::POTSummary> pots;
	if(sr.getByLabel(fPOT_inputTag, pots)){
	  sumdata::POTSummary const& p1(*pots);
	  fpot_tor875 = p1.totpot;
	  fpot_tor875good = p1.totgoodpot;
	  fspill_tor875 = p1.totspills;
	  fspill_tor875good = p1.goodspills;
	}
        else{
	  std::cout << "WARNING:  no sumdata::POTSummary inputTag " << fPOT_inputTag << std::endl;
        }
  }
  fTreePot->Fill();
  // check if MC POT is cumulative for each file as MCC8 overlay
  // sumdata::POTSummary used for Overlay and ?MC
  // Zarko's script used for data (bnb, ext)
  // If fPOT_counting is false, we still need run and subrun info
}

void WireCellAnaTree::ShowerID(int trackId)
{
  //std::cout<<"TrackID: "<<trackId<<std::endl;
  // nested loop to fill ShowerID into a vector
  // initial trackId should be looped over primary particles
  auto const& p = fParticleMap[trackId];
  // if (p.PdgCode() == 111) return; // not pi0 induced --> deferred to tagger
  const TLorentzVector& mom = p.Momentum(0);
  if (mom.E()-mom.M()<fthreshold_showerKE) return; //energy threshold
  if (p.PdgCode() == 11 || p.PdgCode() == -11) fShowerID.push_back(trackId); // key: fill shower ID

  // keep looping
  int Ndaughters =  p.NumberDaughters();
  if (Ndaughters == 0) return;
  for(int i=0; i<Ndaughters; i++){
    //std::cout<<"Daughter ID: "<<p.Daughter(i)<<std::endl;
    ShowerID(p.Daughter(i));
  }
}

void WireCellAnaTree::MuonID(int trackId)
{
  auto const& p = fParticleMap[trackId];
  if (p.Mother()==0 && (p.PdgCode() == 13 || p.PdgCode() == -13)) fMuonID.push_back(trackId);
}

void WireCellAnaTree::ProtonID(int trackId)
{
  auto const& p = fParticleMap[trackId];
  if (p.Mother()==0 && p.PdgCode() == 2212) fProtonID.push_back(trackId);
}


void WireCellAnaTree::resetOutput()
{
	// live period within each event
	// maybe redundant here
	if(f_BDTvars){
		f_neutrino_type = -1;
		f_nuvtx_diff = -1;
		f_showervtx_diff = -1;
		f_muonvtx_diff = -1;

  //single track kdar
  ssm_flag_st_kdar = -999;
  ssm_Nsm = -999;
  ssm_Nsm_wivtx = -999;

  ssm_dq_dx_fwd_1 = -999;
  ssm_dq_dx_fwd_2 = -999;
  ssm_dq_dx_fwd_3 = -999;
  ssm_dq_dx_fwd_4 = -999;
  ssm_dq_dx_fwd_5 = -999;
  ssm_dq_dx_bck_1 = -999;
  ssm_dq_dx_bck_2 = -999;
  ssm_dq_dx_bck_3 = -999;
  ssm_dq_dx_bck_4 = -999;
  ssm_dq_dx_bck_5 = -999;
  ssm_d_dq_dx_fwd_12 = -999;
  ssm_d_dq_dx_fwd_23 = -999;
  ssm_d_dq_dx_fwd_34 = -999;
  ssm_d_dq_dx_fwd_45 = -999;
  ssm_d_dq_dx_bck_12 = -999;
  ssm_d_dq_dx_bck_23 = -999;
  ssm_d_dq_dx_bck_34 = -999;
  ssm_d_dq_dx_bck_45 = -999;
  ssm_max_dq_dx_fwd_3 = -999;
  ssm_max_dq_dx_fwd_5 = -999;
  ssm_max_dq_dx_bck_3 = -999;
  ssm_max_dq_dx_bck_5 = -999;
  ssm_max_d_dq_dx_fwd_3 = -999;
  ssm_max_d_dq_dx_fwd_5 = -999;
  ssm_max_d_dq_dx_bck_3 = -999;
  ssm_max_d_dq_dx_bck_5 = -999;
  ssm_medium_dq_dx = -999;
  ssm_medium_dq_dx_bp = -999;
      //angluar info
  ssm_angle_to_z = -999;
  ssm_angle_to_target = -999;
  ssm_angle_to_absorber = -999;
  ssm_angle_to_vertical = -999;
      //directional info
  ssm_x_dir = -999;
  ssm_y_dir = -999;
  ssm_z_dir = -999;
      //energy info
  ssm_kine_energy = -999;
  ssm_kine_energy_reduced = -999;
      //general properties
  ssm_vtx_activity = -999;
  ssm_pdg = -999;
  ssm_dQ_dx_cut = -999;
  ssm_score_mu_fwd = -999;
  ssm_score_p_fwd = -999;
  ssm_score_e_fwd = -999;
  ssm_score_mu_bck = -999;
  ssm_score_p_bck = -999;
  ssm_score_e_bck = -999;
  ssm_score_mu_fwd_bp = -999;
  ssm_score_p_fwd_bp = -999;
  ssm_score_e_fwd_bp = -999;
      //track "straighness"
  ssm_length = -999;
  ssm_direct_length = -999;
  ssm_length_ratio = -999;
  ssm_max_dev = -999;
    //number of other particles
  ssm_n_prim_tracks_1 = -999;
  ssm_n_prim_tracks_3 = -999;
  ssm_n_prim_tracks_5 = -999;
  ssm_n_prim_tracks_8 = -999;
  ssm_n_prim_tracks_11 = -999;
  ssm_n_all_tracks_1 = -999;
  ssm_n_all_tracks_3 = -999;
  ssm_n_all_tracks_5 = -999;
  ssm_n_all_tracks_8 = -999;
  ssm_n_all_tracks_11 = -999;
  ssm_n_daughter_tracks_1 = -999;
  ssm_n_daughter_tracks_3 = -999;
  ssm_n_daughter_tracks_5 = -999;
  ssm_n_daughter_tracks_8 = -999;
  ssm_n_daughter_tracks_11 = -999;
  ssm_n_daughter_all_1 = -999;
  ssm_n_daughter_all_3 = -999;
  ssm_n_daughter_all_5 = -999;
  ssm_n_daughter_all_8 = -999;
  ssm_n_daughter_all_11 = -999;
    //properties of leading other primary track
  ssm_prim_track1_pdg = -999;
  ssm_prim_track1_score_mu_fwd = -999;
  ssm_prim_track1_score_p_fwd = -999;
  ssm_prim_track1_score_e_fwd = -999;
  ssm_prim_track1_score_mu_bck = -999;
  ssm_prim_track1_score_p_bck = -999;
  ssm_prim_track1_score_e_bck = -999;
  ssm_prim_track1_length = -999;
  ssm_prim_track1_direct_length = -999;
  ssm_prim_track1_length_ratio = -999;
  ssm_prim_track1_max_dev = -999;
  ssm_prim_track1_kine_energy_range = -999;
  ssm_prim_track1_kine_energy_range_mu = -999;
  ssm_prim_track1_kine_energy_range_p = -999;
  ssm_prim_track1_kine_energy_range_e = -999;
  ssm_prim_track1_kine_energy_cal = -999;
  ssm_prim_track1_medium_dq_dx = -999;
  ssm_prim_track1_x_dir = -999;
  ssm_prim_track1_y_dir = -999;
  ssm_prim_track1_z_dir = -999;
  ssm_prim_track1_add_daught_track_counts_1 = -999;
  ssm_prim_track1_add_daught_all_counts_1 = -999;
  ssm_prim_track1_add_daught_track_counts_5 = -999;
  ssm_prim_track1_add_daught_all_counts_5 = -999;
  ssm_prim_track1_add_daught_track_counts_11 = -999;
  ssm_prim_track1_add_daught_all_counts_11 = -999;
    //properties of sub-leading other primary track
  ssm_prim_track2_pdg = -999;
  ssm_prim_track2_score_mu_fwd = -999;
  ssm_prim_track2_score_p_fwd = -999;
  ssm_prim_track2_score_e_fwd = -999;
  ssm_prim_track2_score_mu_bck = -999;
  ssm_prim_track2_score_p_bck = -999;
  ssm_prim_track2_score_e_bck = -999;
  ssm_prim_track2_length = -999;
  ssm_prim_track2_direct_length = -999;
  ssm_prim_track2_length_ratio = -999;
  ssm_prim_track2_max_dev = -999;
  ssm_prim_track2_kine_energy_range = -999;
  ssm_prim_track2_kine_energy_range_mu = -999;
  ssm_prim_track2_kine_energy_range_p = -999;
  ssm_prim_track2_kine_energy_range_e = -999;
  ssm_prim_track2_kine_energy_cal = -999;
  ssm_prim_track2_medium_dq_dx = -999;
  ssm_prim_track2_x_dir = -999;
  ssm_prim_track2_y_dir = -999;
  ssm_prim_track2_z_dir = -999;
  ssm_prim_track2_add_daught_track_counts_1 = -999;
  ssm_prim_track2_add_daught_all_counts_1 = -999;
  ssm_prim_track2_add_daught_track_counts_5 = -999;
  ssm_prim_track2_add_daught_all_counts_5 = -999;
  ssm_prim_track2_add_daught_track_counts_11 = -999;
  ssm_prim_track2_add_daught_all_counts_11 = -999;
    //properties of leading daughter track
  ssm_daught_track1_pdg = -999;
  ssm_daught_track1_score_mu_fwd = -999;
  ssm_daught_track1_score_p_fwd = -999;
  ssm_daught_track1_score_e_fwd = -999;
  ssm_daught_track1_score_mu_bck = -999;
  ssm_daught_track1_score_p_bck = -999;
  ssm_daught_track1_score_e_bck = -999;
  ssm_daught_track1_length = -999;
  ssm_daught_track1_direct_length = -999;
  ssm_daught_track1_length_ratio = -999;
  ssm_daught_track1_max_dev = -999;
  ssm_daught_track1_kine_energy_range = -999;
  ssm_daught_track1_kine_energy_range_mu = -999;
  ssm_daught_track1_kine_energy_range_p = -999;
  ssm_daught_track1_kine_energy_range_e = -999;
  ssm_daught_track1_kine_energy_cal = -999;
  ssm_daught_track1_medium_dq_dx = -999;
  ssm_daught_track1_x_dir = -999;
  ssm_daught_track1_y_dir = -999;
  ssm_daught_track1_z_dir = -999;
  ssm_daught_track1_add_daught_track_counts_1 = -999;
  ssm_daught_track1_add_daught_all_counts_1 = -999;
  ssm_daught_track1_add_daught_track_counts_5 = -999;
  ssm_daught_track1_add_daught_all_counts_5 = -999;
  ssm_daught_track1_add_daught_track_counts_11 = -999;
  ssm_daught_track1_add_daught_all_counts_11 = -999;
    //properties of sub-leading daughter track
  ssm_daught_track2_pdg = -999;
  ssm_daught_track2_score_mu_fwd = -999;
  ssm_daught_track2_score_p_fwd = -999;
  ssm_daught_track2_score_e_fwd = -999;
  ssm_daught_track2_score_mu_bck = -999;
  ssm_daught_track2_score_p_bck = -999;
  ssm_daught_track2_score_e_bck = -999;
  ssm_daught_track2_length = -999;
  ssm_daught_track2_direct_length = -999;
  ssm_daught_track2_length_ratio = -999;
  ssm_daught_track2_max_dev = -999;
  ssm_daught_track2_kine_energy_range = -999;
  ssm_daught_track2_kine_energy_range_mu = -999;
  ssm_daught_track2_kine_energy_range_p = -999;
  ssm_daught_track2_kine_energy_range_e = -999;
  ssm_daught_track2_kine_energy_cal = -999;
  ssm_daught_track2_medium_dq_dx = -999;
  ssm_daught_track2_x_dir = -999;
  ssm_daught_track2_y_dir = -999;
  ssm_daught_track2_z_dir = -999;
  ssm_daught_track2_add_daught_track_counts_1 = -999;
  ssm_daught_track2_add_daught_all_counts_1 = -999;
  ssm_daught_track2_add_daught_track_counts_5 = -999;
  ssm_daught_track2_add_daught_all_counts_5 = -999;
  ssm_daught_track2_add_daught_track_counts_11 = -999;
  ssm_daught_track2_add_daught_all_counts_11 = -999;
    //properties of leading other primary shower
  ssm_prim_shw1_pdg = -999;
  ssm_prim_shw1_score_mu_fwd = -999;
  ssm_prim_shw1_score_p_fwd = -999;
  ssm_prim_shw1_score_e_fwd = -999;
  ssm_prim_shw1_score_mu_bck = -999;
  ssm_prim_shw1_score_p_bck = -999;
  ssm_prim_shw1_score_e_bck = -999;
  ssm_prim_shw1_length = -999;
  ssm_prim_shw1_direct_length = -999;
  ssm_prim_shw1_length_ratio = -999;
  ssm_prim_shw1_max_dev = -999;
  ssm_prim_shw1_kine_energy_range = -999;
  ssm_prim_shw1_kine_energy_range_mu = -999;
  ssm_prim_shw1_kine_energy_range_p = -999;
  ssm_prim_shw1_kine_energy_range_e = -999;
  ssm_prim_shw1_kine_energy_cal = -999;
  ssm_prim_shw1_kine_energy_best = -999;
  ssm_prim_shw1_medium_dq_dx = -999;
  ssm_prim_shw1_x_dir = -999;
  ssm_prim_shw1_y_dir = -999;
  ssm_prim_shw1_z_dir = -999;
  ssm_prim_shw1_add_daught_track_counts_1 = -999;
  ssm_prim_shw1_add_daught_all_counts_1 = -999;
  ssm_prim_shw1_add_daught_track_counts_5 = -999;
  ssm_prim_shw1_add_daught_all_counts_5 = -999;
  ssm_prim_shw1_add_daught_track_counts_11 = -999;
  ssm_prim_shw1_add_daught_all_counts_11 = -999;
    //properties of sub-leading other primary shower
  ssm_prim_shw2_pdg = -999;
  ssm_prim_shw2_score_mu_fwd = -999;
  ssm_prim_shw2_score_p_fwd = -999;
  ssm_prim_shw2_score_e_fwd = -999;
  ssm_prim_shw2_score_mu_bck = -999;
  ssm_prim_shw2_score_p_bck = -999;
  ssm_prim_shw2_score_e_bck = -999;
  ssm_prim_shw2_length = -999;
  ssm_prim_shw2_direct_length = -999;
  ssm_prim_shw2_length_ratio = -999;
  ssm_prim_shw2_max_dev = -999;
  ssm_prim_shw2_kine_energy_range = -999;
  ssm_prim_shw2_kine_energy_range_mu = -999;
  ssm_prim_shw2_kine_energy_range_p = -999;
  ssm_prim_shw2_kine_energy_range_e = -999;
  ssm_prim_shw2_kine_energy_cal = -999;
  ssm_prim_shw2_kine_energy_best = -999;
  ssm_prim_shw2_medium_dq_dx = -999;
  ssm_prim_shw2_x_dir = -999;
  ssm_prim_shw2_y_dir = -999;
  ssm_prim_shw2_z_dir = -999;
  ssm_prim_shw2_add_daught_track_counts_1 = -999;
  ssm_prim_shw2_add_daught_all_counts_1 = -999;
  ssm_prim_shw2_add_daught_track_counts_5 = -999;
  ssm_prim_shw2_add_daught_all_counts_5 = -999;
  ssm_prim_shw2_add_daught_track_counts_11 = -999;
  ssm_prim_shw2_add_daught_all_counts_11 = -999;
    //properties of leading daughter shower
  ssm_daught_shw1_pdg = -999;
  ssm_daught_shw1_score_mu_fwd = -999;
  ssm_daught_shw1_score_p_fwd = -999;
  ssm_daught_shw1_score_e_fwd = -999;
  ssm_daught_shw1_score_mu_bck = -999;
  ssm_daught_shw1_score_p_bck = -999;
  ssm_daught_shw1_score_e_bck = -999;
  ssm_daught_shw1_length = -999;
  ssm_daught_shw1_direct_length = -999;
  ssm_daught_shw1_length_ratio = -999;
  ssm_daught_shw1_max_dev = -999;
  ssm_daught_shw1_kine_energy_range = -999;
  ssm_daught_shw1_kine_energy_range_mu = -999;
  ssm_daught_shw1_kine_energy_range_p = -999;
  ssm_daught_shw1_kine_energy_range_e = -999;
  ssm_daught_shw1_kine_energy_cal = -999;
  ssm_daught_shw1_kine_energy_best = -999;
  ssm_daught_shw1_medium_dq_dx = -999;
  ssm_daught_shw1_x_dir = -999;
  ssm_daught_shw1_y_dir = -999;
  ssm_daught_shw1_z_dir = -999;
  ssm_daught_shw1_add_daught_track_counts_1 = -999;
  ssm_daught_shw1_add_daught_all_counts_1 = -999;
  ssm_daught_shw1_add_daught_track_counts_5 = -999;
  ssm_daught_shw1_add_daught_all_counts_5 = -999;
  ssm_daught_shw1_add_daught_track_counts_11 = -999;
  ssm_daught_shw1_add_daught_all_counts_11 = -999;
    //properties of sub-leading daughter shower
  ssm_daught_shw2_pdg = -999;
  ssm_daught_shw2_score_mu_fwd = -999;
  ssm_daught_shw2_score_p_fwd = -999;
  ssm_daught_shw2_score_e_fwd = -999;
  ssm_daught_shw2_score_mu_bck = -999;
  ssm_daught_shw2_score_p_bck = -999;
  ssm_daught_shw2_score_e_bck = -999;
  ssm_daught_shw2_length = -999;
  ssm_daught_shw2_direct_length = -999;
  ssm_daught_shw2_length_ratio = -999;
  ssm_daught_shw2_max_dev = -999;
  ssm_daught_shw2_kine_energy_range = -999;
  ssm_daught_shw2_kine_energy_range_mu = -999;
  ssm_daught_shw2_kine_energy_range_p = -999;
  ssm_daught_shw2_kine_energy_range_e = -999;
  ssm_daught_shw2_kine_energy_cal = -999;
  ssm_daught_shw2_kine_energy_best = -999;
  ssm_daught_shw2_medium_dq_dx = -999;
  ssm_daught_shw2_x_dir = -999;
  ssm_daught_shw2_y_dir = -999;
  ssm_daught_shw2_z_dir = -999;
  ssm_daught_shw2_add_daught_track_counts_1 = -999;
  ssm_daught_shw2_add_daught_all_counts_1 = -999;
  ssm_daught_shw2_add_daught_track_counts_5 = -999;
  ssm_daught_shw2_add_daught_all_counts_5 = -999;
  ssm_daught_shw2_add_daught_track_counts_11 = -999;
  ssm_daught_shw2_add_daught_all_counts_11 = -999;
    //event level properties
  ssm_nu_angle_z = -999;
  ssm_nu_angle_target = -999;
  ssm_nu_angle_absorber = -999;
  ssm_nu_angle_vertical = -999;
  ssm_con_nu_angle_z = -999;
  ssm_con_nu_angle_target = -999;
  ssm_con_nu_angle_absorber = -999;
  ssm_con_nu_angle_vertical = -999;
  ssm_prim_nu_angle_z = -999;
  ssm_prim_nu_angle_target = -999;
  ssm_prim_nu_angle_absorber = -999;
  ssm_prim_nu_angle_vertical = -999;
  ssm_track_angle_z = -999;
  ssm_track_angle_target = -999;
  ssm_track_angle_absorber = -999;
  ssm_track_angle_vertical = -999;
  ssm_vtxX = -999;
  ssm_vtxY = -999;
  ssm_vtxZ = -999;
    //off vertex stuff
  ssm_offvtx_length  =  -999;
  ssm_offvtx_energy  =  -999;
  ssm_n_offvtx_tracks_1  =  -999;
  ssm_n_offvtx_tracks_3  =  -999;
  ssm_n_offvtx_tracks_5  =  -999;
  ssm_n_offvtx_tracks_8  =  -999;
  ssm_n_offvtx_tracks_11  =  -999;
  ssm_n_offvtx_showers_1  =  -999;
  ssm_n_offvtx_showers_3  =  -999;
  ssm_n_offvtx_showers_5  =  -999;
  ssm_n_offvtx_showers_8  =  -999;
  ssm_n_offvtx_showers_11  =  -999;
    //properties of leading off vertex track
  ssm_offvtx_track1_pdg  =  -999;
  ssm_offvtx_track1_score_mu_fwd  =  -999;
  ssm_offvtx_track1_score_p_fwd  =  -999;
  ssm_offvtx_track1_score_e_fwd  =  -999;
  ssm_offvtx_track1_score_mu_bck  =  -999;
  ssm_offvtx_track1_score_p_bck  =  -999;
  ssm_offvtx_track1_score_e_bck  =  -999;
  ssm_offvtx_track1_length  =  -999;
  ssm_offvtx_track1_direct_length  =  -999;
  ssm_offvtx_track1_max_dev  =  -999;
  ssm_offvtx_track1_kine_energy_range  =  -999;
  ssm_offvtx_track1_kine_energy_range_mu  = -999;
  ssm_offvtx_track1_kine_energy_range_p  =  -999;
  ssm_offvtx_track1_kine_energy_range_e  =  -999;
  ssm_offvtx_track1_kine_energy_cal  =  -999;
  ssm_offvtx_track1_medium_dq_dx  =  -999;
  ssm_offvtx_track1_x_dir  =  -999;
  ssm_offvtx_track1_y_dir  =  -999;
  ssm_offvtx_track1_z_dir  =  -999;
  ssm_offvtx_track1_dist_mainvtx  =  -999;
    //properties of leading off vertex shower
  ssm_offvtx_shw1_pdg_offvtx  =  -999;
  ssm_offvtx_shw1_score_mu_fwd  = -999;
  ssm_offvtx_shw1_score_p_fwd  =  -999;
  ssm_offvtx_shw1_score_e_fwd  =  -999;
  ssm_offvtx_shw1_score_mu_bck  =  -999;
  ssm_offvtx_shw1_score_p_bck  =  -999;
  ssm_offvtx_shw1_score_e_bck  =  -999;
  ssm_offvtx_shw1_length  =  -999;
  ssm_offvtx_shw1_direct_length  =  -999;
  ssm_offvtx_shw1_max_dev  =  -999;
  ssm_offvtx_shw1_kine_energy_best  =  -999;
  ssm_offvtx_shw1_kine_energy_range  =  -999;
  ssm_offvtx_shw1_kine_energy_range_mu  =  -999;
  ssm_offvtx_shw1_kine_energy_range_p  =  -999;
  ssm_offvtx_shw1_kine_energy_range_e  =  -999;
  ssm_offvtx_shw1_kine_energy_cal  =  -999;
  ssm_offvtx_shw1_medium_dq_dx  =  -999;
  ssm_offvtx_shw1_x_dir  =  -999;
  ssm_offvtx_shw1_y_dir  =  -999;
  ssm_offvtx_shw1_z_dir  =  -999;
  ssm_offvtx_shw1_dist_mainvtx  =  -999;
    // Spacepoints
  ssmsp_Ntrack = 0;
  ssmsp_Nsp=nullptr;
  ssmsp_Nsp_tot = 0;
  ssmsp_pdg=nullptr;
  ssmsp_id=nullptr;
  ssmsp_mother=nullptr;
  ssmsp_x=nullptr;
  ssmsp_y=nullptr;
  ssmsp_z=nullptr;
  ssmsp_dx=nullptr;
  ssmsp_dQ=nullptr;
  ssmsp_KE=nullptr;
  ssmsp_containing_shower_id=nullptr;
  ssmsp_containing_shower_ke=nullptr;
  ssmsp_containing_shower_flag=nullptr;
    //Kine vars
  ssm_kine_reco_Enu=-1;
  ssm_kine_reco_add_energy=-1;
  ssm_kine_energy_particle=nullptr;
  ssm_kine_energy_info=nullptr;
  ssm_kine_particle_type=nullptr;
  ssm_kine_energy_included=nullptr;
  ssm_kine_pio_mass=-1;
  ssm_kine_pio_flag=-1;
  ssm_kine_pio_vtx_dis=-1;
  ssm_kine_pio_energy_1=-1;
  ssm_kine_pio_theta_1=-1;
  ssm_kine_pio_phi_1=-1;
  ssm_kine_pio_dis_1=-1;
  ssm_kine_pio_energy_2=-1;
  ssm_kine_pio_theta_2=-1;
  ssm_kine_pio_phi_2=-1;
  ssm_kine_pio_dis_2=-1;
  ssm_kine_pio_angle=-1;
  ssm_numu_cc_flag = -1;
  ssm_cosmict_flag_1=-1; // fiducial volume vertex
  ssm_cosmict_flag_2=-1;  // single muon
  ssm_cosmict_flag_3=-1;  // single muon (long)
  ssm_cosmict_flag_4=-1;  // kinematics muon
  ssm_cosmict_flag_5=-1; // kinematics muon (long)
  ssm_cosmict_flag_6=-1; // special ...
  ssm_cosmict_flag_7=-1;  // muon+ michel
  ssm_cosmict_flag_8=-1;  // muon + michel + special
  ssm_cosmict_flag_9=-1;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
  ssm_cosmict_flag_10=nullptr;  // front upstream (dirt)
  ssm_cosmict_flag=-1;

    // single photon shower identification
  shw_sp_flag = -1;
  shw_sp_num_mip_tracks = -1;
  shw_sp_num_muons = -1;
  shw_sp_num_pions = -1;
  shw_sp_num_protons = -1;
  shw_sp_proton_length_1 = -1;
  shw_sp_proton_dqdx_1 = -1;
  shw_sp_proton_energy_1 = -1;
  shw_sp_proton_length_2 = -1;
  shw_sp_proton_dqdx_2 = -1;
  shw_sp_proton_energy_2 = -1;
  shw_sp_n_good_showers = -1;
  shw_sp_n_20mev_showers = -1;
  shw_sp_n_br1_showers = -1;
  shw_sp_n_br2_showers = -1;
  shw_sp_n_br3_showers = -1;
  shw_sp_n_br4_showers = -1;
  shw_sp_n_20mev_showers = -1;
  shw_sp_n_br1_showers = -1;
  shw_sp_n_br2_showers = -1;
  shw_sp_n_br3_showers = -1;
  shw_sp_n_br4_showers = -1;
  shw_sp_n_20br1_showers = -1;
  shw_sp_shw_vtx_dis = -1;
  shw_sp_max_shw_dis = -1;
  shw_sp_energy = -1;
  shw_sp_max_dQ_dx_sample = -1;
  shw_sp_vec_dQ_dx_0 = -1;
  shw_sp_vec_dQ_dx_1 = -1;
  shw_sp_n_below_threshold = -1;
  shw_sp_n_good_tracks = -1;
  shw_sp_n_vertex = -1;
  shw_sp_angle_beam = -1;
  shw_sp_flag_all_above = false;
  shw_sp_length_main = -1;
  shw_sp_length_total = -1;
  shw_sp_min_dQ_dx_5 = -1;
  shw_sp_lowest_dQ_dx = -1;
  shw_sp_iso_angle = -1;
  shw_sp_n_below_zero = -1;
  shw_sp_highest_dQ_dx = -1;
  shw_sp_n_lowest = -1;
  shw_sp_n_highest = -1;
  shw_sp_stem_length = -1;
  shw_sp_E_indirect_max_energy = -1;
  shw_sp_flag_stem_trajectory = -1;
  shw_sp_min_dis = -1;
  shw_sp_n_other_vertex = 2;
  shw_sp_n_stem_size = 20;
  shw_sp_medium_dQ_dx = -1;
  shw_sp_filled = -1;

  shw_sp_vec_dQ_dx_2 = -1;
  shw_sp_vec_dQ_dx_3 = -1;
  shw_sp_vec_dQ_dx_4 = -1;
  shw_sp_vec_dQ_dx_5 = -1;
  shw_sp_vec_dQ_dx_6 = -1;
  shw_sp_vec_dQ_dx_7 = -1;
  shw_sp_vec_dQ_dx_8 = -1;
  shw_sp_vec_dQ_dx_9 = -1;
  shw_sp_vec_dQ_dx_10 = -1;
  shw_sp_vec_dQ_dx_11 = -1;
  shw_sp_vec_dQ_dx_12 = -1;
  shw_sp_vec_dQ_dx_13 = -1;
  shw_sp_vec_dQ_dx_14 = -1;
  shw_sp_vec_dQ_dx_15 = -1;
  shw_sp_vec_dQ_dx_16 = -1;
  shw_sp_vec_dQ_dx_17 = -1;
  shw_sp_vec_dQ_dx_18 = -1;
  shw_sp_vec_dQ_dx_19 = -1;
  shw_sp_vec_median_dedx = -1;
  shw_sp_vec_mean_dedx = -1;

  // shower pi0 identification
  shw_sp_pio_flag = -1;
  shw_sp_pio_mip_id = -1;
  shw_sp_pio_filled = -1;
  shw_sp_pio_flag_pio = -1;

  shw_sp_pio_1_flag = -1;
  shw_sp_pio_1_mass = -1;
  shw_sp_pio_1_pio_type = -1;
  shw_sp_pio_1_energy_1 = -1;
  shw_sp_pio_1_energy_2 = -1;
  shw_sp_pio_1_dis_1 = -1;
  shw_sp_pio_1_dis_2 = -1;

  // bad reconstruction
  shw_sp_br_filled = -1;

  shw_sp_br1_flag = -1;
  // br1_1
  shw_sp_br1_1_flag = -1;
  shw_sp_br1_1_shower_type = -1;
  shw_sp_br1_1_vtx_n_segs = -1;
  shw_sp_br1_1_energy = -1;
  shw_sp_br1_1_n_segs = -1;
  shw_sp_br1_1_flag_sg_topology = -1;
  shw_sp_br1_1_flag_sg_trajectory = -1;
  shw_sp_br1_1_sg_length = -1;

  //br1_2
  shw_sp_br1_2_flag = -1;
  shw_sp_br1_2_energy = -1;
  shw_sp_br1_2_n_connected = -1;
  shw_sp_br1_2_max_length = -1;
  shw_sp_br1_2_n_connected_1 = -1;
  shw_sp_br1_2_vtx_n_segs = -1;
  shw_sp_br1_2_n_shower_segs = -1;
  shw_sp_br1_2_max_length_ratio = -1;
  shw_sp_br1_2_shower_length = -1;

  //br1_3
  shw_sp_br1_3_flag = -1;
  shw_sp_br1_3_energy = -1;
  shw_sp_br1_3_n_connected_p = -1;
  shw_sp_br1_3_max_length_p = -1;
  shw_sp_br1_3_n_shower_segs = -1;
  shw_sp_br1_3_flag_sg_topology = -1;
  shw_sp_br1_3_flag_sg_trajectory = -1;
  shw_sp_br1_3_n_shower_main_segs = -1;
  shw_sp_br1_3_sg_length = -1;

  // br2
  shw_sp_br2_flag = -1;
  shw_sp_br2_flag_single_shower = -1;
  shw_sp_br2_num_valid_tracks = -1;
  shw_sp_br2_energy = -1;
  shw_sp_br2_angle1 = -1;
  shw_sp_br2_angle2 = -1;
  shw_sp_br2_angle = -1;
  shw_sp_br2_angle3 = -1;
  shw_sp_br2_n_shower_main_segs = -1;
  shw_sp_br2_max_angle = -1;
  shw_sp_br2_sg_length = -1;
  shw_sp_br2_flag_sg_trajectory = -1;

  // low-energy overlap
  shw_sp_lol_flag = -1;
  shw_sp_lol_3_flag = -1;
  shw_sp_lol_3_angle_beam = -1;
  shw_sp_lol_3_min_angle = -1;
  shw_sp_lol_3_n_valid_tracks = -1;
  shw_sp_lol_3_vtx_n_segs = -1;
  shw_sp_lol_3_energy = -1;
  shw_sp_lol_3_shower_main_length = -1;
  shw_sp_lol_3_n_sum = -1;
  shw_sp_lol_3_n_out = -1;


  // br3
  shw_sp_br3_1_energy = -1;
  shw_sp_br3_1_n_shower_segments = -1;
  shw_sp_br3_1_sg_flag_trajectory = -1;
  shw_sp_br3_1_sg_direct_length = -1;
  shw_sp_br3_1_sg_length = -1;
  shw_sp_br3_1_total_main_length = -1;
  shw_sp_br3_1_total_length = -1;
  shw_sp_br3_1_iso_angle = -1;
  shw_sp_br3_1_sg_flag_topology = -1;
  shw_sp_br3_1_flag = -1;

  shw_sp_br3_2_n_ele = -1;
  shw_sp_br3_2_n_other = -1;
  shw_sp_br3_2_energy = -1;
  shw_sp_br3_2_total_main_length = -1;
  shw_sp_br3_2_total_length = -1;
  shw_sp_br3_2_other_fid = -1;
  shw_sp_br3_2_flag = -1;

  shw_sp_br3_4_acc_length = -1;
  shw_sp_br3_4_total_length = -1;
  shw_sp_br3_4_energy = -1;
  shw_sp_br3_4_flag = -1;

  shw_sp_br3_7_energy = -1;
  shw_sp_br3_7_min_angle = -1;
  shw_sp_br3_7_sg_length = -1;
  shw_sp_br3_7_shower_main_length = -1;
  shw_sp_br3_7_flag = -1;

  shw_sp_br3_8_max_dQ_dx = -1;
  shw_sp_br3_8_energy = -1;
  shw_sp_br3_8_n_main_segs = -1;
  shw_sp_br3_8_shower_main_length = -1;
  shw_sp_br3_8_shower_length = -1;
  shw_sp_br3_8_flag = -1;

  shw_sp_br3_flag = -1;


  shw_sp_br4_1_shower_main_length = -1;
  shw_sp_br4_1_shower_total_length = -1;
  shw_sp_br4_1_min_dis = -1;
  shw_sp_br4_1_energy = -1;
  shw_sp_br4_1_flag_avoid_muon_check = -1;
  shw_sp_br4_1_n_vtx_segs = -1;
  shw_sp_br4_1_n_main_segs = -1;
  shw_sp_br4_1_flag = -1;

  shw_sp_br4_2_ratio_45 = -1;
  shw_sp_br4_2_ratio_35 = -1;
  shw_sp_br4_2_ratio_25 = -1;
  shw_sp_br4_2_ratio_15 = -1;
  shw_sp_br4_2_energy = -1;
  shw_sp_br4_2_ratio1_45 = -1;
  shw_sp_br4_2_ratio1_35 = -1;
  shw_sp_br4_2_ratio1_25 = -1;
  shw_sp_br4_2_ratio1_15 = -1;
  shw_sp_br4_2_iso_angle = -1;
  shw_sp_br4_2_iso_angle1 = -1;
  shw_sp_br4_2_angle = -1;
  shw_sp_br4_2_flag = -1;

  shw_sp_br4_flag = -1;


  shw_sp_hol_1_n_valid_tracks = -1;
  shw_sp_hol_1_min_angle = -1;
  shw_sp_hol_1_energy = -1;
  shw_sp_hol_1_flag_all_shower = -1;
  shw_sp_hol_1_min_length = -1;
  shw_sp_hol_1_flag = -1;

  shw_sp_hol_2_min_angle = -1;
  shw_sp_hol_2_medium_dQ_dx = -1;
  shw_sp_hol_2_ncount = -1;
  shw_sp_hol_2_energy = -1;
  shw_sp_hol_2_flag = -1;

  shw_sp_hol_flag = -1;

  shw_sp_lem_shower_total_length = -1;
  shw_sp_lem_shower_main_length = -1;
  shw_sp_lem_n_3seg = -1;
  shw_sp_lem_e_charge = -1;
  shw_sp_lem_e_dQdx = -1;
  shw_sp_lem_shower_num_segs = -1;
  shw_sp_lem_shower_num_main_segs = -1;
  shw_sp_lem_flag = -1;

		cosmic_filled=-1;
		cosmic_flag=-1;
		cosmic_n_solid_tracks=-1;
		cosmic_energy_main_showers=-1;
		cosmic_energy_direct_showers=-1;
		cosmic_energy_indirect_showers=-1;
		cosmic_n_direct_showers=-1;
		cosmic_n_indirect_showers=-1;
		cosmic_n_main_showers=-1;
		gap_filled=-1;
		gap_flag=-1;
		gap_flag_prolong_u=-1;
		gap_flag_prolong_v=-1;
		gap_flag_prolong_w=-1;
		gap_flag_parallel=-1;
		gap_n_points=-1;
		gap_n_bad=-1;
		gap_energy=-1;
		gap_num_valid_tracks=-1;
		gap_flag_single_shower=-1;
		mip_quality_filled=-1;
		mip_quality_flag=-1;
		mip_quality_energy=-1;
		mip_quality_overlap=-1;
		mip_quality_n_showers=-1;
		mip_quality_n_tracks=-1;
		mip_quality_flag_inside_pi0=-1;
		mip_quality_n_pi0_showers=-1;
		mip_quality_shortest_length=-1;
		mip_quality_acc_length=-1;
		mip_quality_shortest_angle=-1;
		mip_quality_flag_proton=-1;
		mip_filled=-1;
		mip_flag=-1;
		mip_energy=-1;
		mip_n_end_reduction=-1;
		mip_n_first_mip=-1;
		mip_n_first_non_mip=-1;
		mip_n_first_non_mip_1=-1;
		mip_n_first_non_mip_2=-1;
		mip_vec_dQ_dx_0=-1;
		mip_vec_dQ_dx_1=-1;
		mip_max_dQ_dx_sample=-1;
		mip_n_below_threshold=-1;
		mip_n_below_zero=-1;
		mip_n_lowest=-1;
		mip_n_highest=-1;
		mip_lowest_dQ_dx=-1;
		mip_highest_dQ_dx=-1;
		mip_medium_dQ_dx=-1;
		mip_stem_length=-1;
		mip_length_main=-1;
		mip_length_total=-1;
		mip_angle_beam=-1;
		mip_iso_angle=-1;
		mip_n_vertex=-1;
		mip_n_good_tracks=-1;
		mip_E_indirect_max_energy=-1;
		mip_flag_all_above=-1;
		mip_min_dQ_dx_5=-1;
		mip_n_other_vertex=-1;
		mip_n_stem_size=-1;
		mip_flag_stem_trajectory=-1;
		mip_min_dis=-1;
		mip_vec_dQ_dx_2=-1;
		mip_vec_dQ_dx_3=-1;
		mip_vec_dQ_dx_4=-1;
		mip_vec_dQ_dx_5=-1;
		mip_vec_dQ_dx_6=-1;
		mip_vec_dQ_dx_7=-1;
		mip_vec_dQ_dx_8=-1;
		mip_vec_dQ_dx_9=-1;
		mip_vec_dQ_dx_10=-1;
		mip_vec_dQ_dx_11=-1;
		mip_vec_dQ_dx_12=-1;
		mip_vec_dQ_dx_13=-1;
		mip_vec_dQ_dx_14=-1;
		mip_vec_dQ_dx_15=-1;
		mip_vec_dQ_dx_16=-1;
		mip_vec_dQ_dx_17=-1;
		mip_vec_dQ_dx_18=-1;
		mip_vec_dQ_dx_19=-1;
		pio_filled=-1;
		pio_flag=-1;
		pio_mip_id=-1;
		pio_flag_pio=-1;
		pio_1_flag=-1;
		pio_1_mass=-1;
		pio_1_pio_type=-1;
		pio_1_energy_1=-1;
		pio_1_energy_2=-1;
		pio_1_dis_1=-1;
		pio_1_dis_2=-1;
		pio_2_v_flag=nullptr;
		pio_2_v_dis2=nullptr;
		pio_2_v_angle2=nullptr;
		pio_2_v_acc_length=nullptr;
		sig_flag=-1;
		sig_1_v_flag=nullptr;
		sig_1_v_angle=nullptr;
		sig_1_v_flag_single_shower=nullptr;
		sig_1_v_energy=nullptr;
		sig_1_v_energy_1=nullptr;
		sig_2_v_flag=nullptr;
		sig_2_v_energy=nullptr;
		sig_2_v_shower_angle=nullptr;
		sig_2_v_flag_single_shower=nullptr;
		sig_2_v_medium_dQ_dx=nullptr;
		sig_2_v_start_dQ_dx=nullptr;
		mgo_flag=-1;
		mgo_energy=-1;
		mgo_max_energy=-1;
		mgo_total_energy=-1;
		mgo_n_showers=-1;
		mgo_max_energy_1=-1;
		mgo_max_energy_2=-1;
		mgo_total_other_energy=-1;
		mgo_n_total_showers=-1;
		mgo_total_other_energy_1=-1;
		mgt_flag=-1;
		mgt_flag_single_shower=-1;
		mgt_max_energy=-1;
		mgt_energy=-1;
		mgt_total_other_energy=-1;
		mgt_max_energy_1=-1;
		mgt_e_indirect_max_energy=-1;
		mgt_e_direct_max_energy=-1;
		mgt_n_direct_showers=-1;
		mgt_e_direct_total_energy=-1;
		mgt_flag_indirect_max_pio=-1;
		mgt_e_indirect_total_energy=-1;
		stw_flag=-1;
		stw_1_flag=-1;
		stw_1_energy=-1;
		stw_1_dis=-1;
		stw_1_dQ_dx=-1;
		stw_1_flag_single_shower=-1;
		stw_1_n_pi0=-1;
		stw_1_num_valid_tracks=-1;
		stw_2_v_flag=nullptr;
		stw_2_v_medium_dQ_dx=nullptr;
		stw_2_v_energy=nullptr;
		stw_2_v_angle=nullptr;
		stw_2_v_dir_length=nullptr;
		stw_2_v_max_dQ_dx=nullptr;
		stw_3_v_flag=nullptr;
		stw_3_v_angle=nullptr;
		stw_3_v_dir_length=nullptr;
		stw_3_v_energy=nullptr;
		stw_3_v_medium_dQ_dx=nullptr;
		stw_4_v_flag=nullptr;
		stw_4_v_angle=nullptr;
		stw_4_v_dis=nullptr;
		stw_4_v_energy=nullptr;
		spt_flag=-1;
		spt_flag_single_shower=-1;
		spt_energy=-1;
		spt_shower_main_length=-1;
		spt_shower_total_length=-1;
		spt_angle_beam=-1;
		spt_angle_vertical=-1;
		spt_max_dQ_dx=-1;
		spt_angle_beam_1=-1;
		spt_angle_drift=-1;
		spt_angle_drift_1=-1;
		spt_num_valid_tracks=-1;
		spt_n_vtx_segs=-1;
		spt_max_length=-1;
		stem_len_flag=-1;
		stem_len_energy=-1;
		stem_len_length=-1;
		stem_len_flag_avoid_muon_check=-1;
		stem_len_num_daughters=-1;
		stem_len_daughter_length=-1;
		lem_flag=-1;
		lem_shower_total_length=-1;
		lem_shower_main_length=-1;
		lem_n_3seg=-1;
		lem_e_charge=-1;
		lem_e_dQdx=-1;
		lem_shower_num_segs=-1;
		lem_shower_num_main_segs=-1;
		brm_flag=-1;
		brm_n_mu_segs=-1;
		brm_Ep=-1;
		brm_energy=-1;
		brm_acc_length=-1;
		brm_shower_total_length=-1;
		brm_connected_length=-1;
		brm_n_size=-1;
		brm_acc_direct_length=-1;
		brm_n_shower_main_segs=-1;
		brm_n_mu_main=-1;
		cme_flag=-1;
		cme_mu_energy=-1;
		cme_energy=-1;
		cme_mu_length=-1;
		cme_length=-1;
		cme_angle_beam=-1;
		anc_flag=-1;
		anc_energy=-1;
		anc_angle=-1;
		anc_max_angle=-1;
		anc_max_length=-1;
		anc_acc_forward_length=-1;
		anc_acc_backward_length=-1;
		anc_acc_forward_length1=-1;
		anc_shower_main_length=-1;
		anc_shower_total_length=-1;
		anc_flag_main_outside=-1;
		stem_dir_filled=-1;
		stem_dir_flag=-1;
		stem_dir_flag_single_shower=-1;
		stem_dir_angle=-1;
		stem_dir_energy=-1;
		stem_dir_angle1=-1;
		stem_dir_angle2=-1;
		stem_dir_angle3=-1;
		stem_dir_ratio=-1;
		vis_flag=-1;
		vis_1_filled=-1;
		vis_1_flag=-1;
		vis_1_n_vtx_segs=-1;
		vis_1_energy=-1;
		vis_1_num_good_tracks=-1;
		vis_1_max_angle=-1;
		vis_1_max_shower_angle=-1;
		vis_1_tmp_length1=-1;
		vis_1_tmp_length2=-1;
		vis_1_particle_type=-1;
		vis_2_filled=-1;
		vis_2_flag=-1;
		vis_2_n_vtx_segs=-1;
		vis_2_min_angle=-1;
		vis_2_min_weak_track=-1;
		vis_2_angle_beam=-1;
		vis_2_min_angle1=-1;
		vis_2_iso_angle1=-1;
		vis_2_min_medium_dQ_dx=-1;
		vis_2_min_length=-1;
		vis_2_sg_length=-1;
		vis_2_max_angle=-1;
		vis_2_max_weak_track=-1;
		br_filled=-1;
		br1_flag=-1;
		br1_1_flag=-1;
		br1_1_shower_type=-1;
		br1_1_vtx_n_segs=-1;
		br1_1_energy=-1;
		br1_1_n_segs=-1;
		br1_1_flag_sg_topology=-1;
		br1_1_flag_sg_trajectory=-1;
		br1_1_sg_length=-1;
		br1_2_flag=-1;
		br1_2_energy=-1;
		br1_2_n_connected=-1;
		br1_2_max_length=-1;
		br1_2_n_connected_1=-1;
		br1_2_vtx_n_segs=-1;
		br1_2_n_shower_segs=-1;
		br1_2_max_length_ratio=-1;
		br1_2_shower_length=-1;
		br1_3_flag=-1;
		br1_3_energy=-1;
		br1_3_n_connected_p=-1;
		br1_3_max_length_p=-1;
		br1_3_n_shower_segs=-1;
		br1_3_flag_sg_topology=-1;
		br1_3_flag_sg_trajectory=-1;
		br1_3_n_shower_main_segs=-1;
		br1_3_sg_length=-1;
		br2_flag=-1;
		br2_flag_single_shower=-1;
		br2_num_valid_tracks=-1;
		br2_energy=-1;
		br2_angle1=-1;
		br2_angle2=-1;
		br2_angle=-1;
		br2_angle3=-1;
		br2_n_shower_main_segs=-1;
		br2_max_angle=-1;
		br2_sg_length=-1;
		br2_flag_sg_trajectory=-1;
		br3_flag=-1;
		br3_1_flag=-1;
		br3_1_energy=-1;
		br3_1_n_shower_segments=-1;
		br3_1_sg_flag_trajectory=-1;
		br3_1_sg_direct_length=-1;
		br3_1_sg_length=-1;
		br3_1_total_main_length=-1;
		br3_1_total_length=-1;
		br3_1_iso_angle=-1;
		br3_1_sg_flag_topology=-1;
		br3_2_flag=-1;
		br3_2_n_ele=-1;
		br3_2_n_other=-1;
		br3_2_energy=-1;
		br3_2_total_main_length=-1;
		br3_2_total_length=-1;
		br3_2_other_fid=-1;
		br3_3_v_flag=nullptr;
		br3_3_v_energy=nullptr;
		br3_3_v_angle=nullptr;
		br3_3_v_dir_length=nullptr;
		br3_3_v_length=nullptr;
		br3_4_flag=-1;
		br3_4_acc_length=-1;
		br3_4_total_length=-1;
		br3_4_energy=-1;
		br3_5_v_flag=nullptr;
		br3_5_v_dir_length=nullptr;
		br3_5_v_total_length=nullptr;
		br3_5_v_flag_avoid_muon_check=nullptr;
		br3_5_v_n_seg=nullptr;
		br3_5_v_angle=nullptr;
		br3_5_v_sg_length=nullptr;
		br3_5_v_energy=nullptr;
		br3_5_v_n_main_segs=nullptr;
		br3_5_v_n_segs=nullptr;
		br3_5_v_shower_main_length=nullptr;
		br3_5_v_shower_total_length=nullptr;
		br3_6_v_flag=nullptr;
		br3_6_v_angle=nullptr;
		br3_6_v_angle1=nullptr;
		br3_6_v_flag_shower_trajectory=nullptr;
		br3_6_v_direct_length=nullptr;
		br3_6_v_length=nullptr;
		br3_6_v_n_other_vtx_segs=nullptr;
		br3_6_v_energy=nullptr;
		br3_7_flag=-1;
		br3_7_energy=-1;
		br3_7_min_angle=-1;
		br3_7_sg_length=-1;
		br3_7_shower_main_length=-1;
		br3_8_flag=-1;
		br3_8_max_dQ_dx=-1;
		br3_8_energy=-1;
		br3_8_n_main_segs=-1;
		br3_8_shower_main_length=-1;
		br3_8_shower_length=-1;
		br4_flag=-1;
		br4_1_flag=-1;
		br4_1_shower_main_length=-1;
		br4_1_shower_total_length=-1;
		br4_1_min_dis=-1;
		br4_1_energy=-1;
		br4_1_flag_avoid_muon_check=-1;
		br4_1_n_vtx_segs=-1;
		br4_1_n_main_segs=-1;
		br4_2_flag=-1;
		br4_2_ratio_45=-1;
		br4_2_ratio_35=-1;
		br4_2_ratio_25=-1;
		br4_2_ratio_15=-1;
		br4_2_energy=-1;
		br4_2_ratio1_45=-1;
		br4_2_ratio1_35=-1;
		br4_2_ratio1_25=-1;
		br4_2_ratio1_15=-1;
		br4_2_iso_angle=-1;
		br4_2_iso_angle1=-1;
		br4_2_angle=-1;
		tro_flag=-1;
		tro_1_v_flag=nullptr;
		tro_1_v_particle_type=nullptr;
		tro_1_v_flag_dir_weak=nullptr;
		tro_1_v_min_dis=nullptr;
		tro_1_v_sg1_length=nullptr;
		tro_1_v_shower_main_length=nullptr;
		tro_1_v_max_n_vtx_segs=nullptr;
		tro_1_v_tmp_length=nullptr;
		tro_1_v_medium_dQ_dx=nullptr;
		tro_1_v_dQ_dx_cut=nullptr;
		tro_1_v_flag_shower_topology=nullptr;
		tro_2_v_flag=nullptr;
		tro_2_v_energy=nullptr;
		tro_2_v_stem_length=nullptr;
		tro_2_v_iso_angle=nullptr;
		tro_2_v_max_length=nullptr;
		tro_2_v_angle=nullptr;
		tro_3_flag=-1;
		tro_3_stem_length=-1;
		tro_3_n_muon_segs=-1;
		tro_3_energy=-1;
		tro_4_v_flag=nullptr;
		tro_4_v_dir2_mag=nullptr;
		tro_4_v_angle=nullptr;
		tro_4_v_angle1=nullptr;
		tro_4_v_angle2=nullptr;
		tro_4_v_length=nullptr;
		tro_4_v_length1=nullptr;
		tro_4_v_medium_dQ_dx=nullptr;
		tro_4_v_end_dQ_dx=nullptr;
		tro_4_v_energy=nullptr;
		tro_4_v_shower_main_length=nullptr;
		tro_4_v_flag_shower_trajectory=nullptr;
		tro_5_v_flag=nullptr;
		tro_5_v_max_angle=nullptr;
		tro_5_v_min_angle=nullptr;
		tro_5_v_max_length=nullptr;
		tro_5_v_iso_angle=nullptr;
		tro_5_v_n_vtx_segs=nullptr;
		tro_5_v_min_count=nullptr;
		tro_5_v_max_count=nullptr;
		tro_5_v_energy=nullptr;
		hol_flag=-1;
		hol_1_flag=-1;
		hol_1_n_valid_tracks=-1;
		hol_1_min_angle=-1;
		hol_1_energy=-1;
		hol_1_flag_all_shower=-1;
		hol_1_min_length=-1;
		hol_2_flag=-1;
		hol_2_min_angle=-1;
		hol_2_medium_dQ_dx=-1;
		hol_2_ncount=-1;
		hol_2_energy=-1;
		lol_flag=-1;
		lol_1_v_flag=nullptr;
		lol_1_v_energy=nullptr;
		lol_1_v_vtx_n_segs=nullptr;
		lol_1_v_nseg=nullptr;
		lol_1_v_angle=nullptr;
		lol_2_v_flag=nullptr;
		lol_2_v_length=nullptr;
		lol_2_v_angle=nullptr;
		lol_2_v_type=nullptr;
		lol_2_v_vtx_n_segs=nullptr;
		lol_2_v_energy=nullptr;
		lol_2_v_shower_main_length=nullptr;
		lol_2_v_flag_dir_weak=nullptr;
		lol_3_flag=-1;
		lol_3_angle_beam=-1;
		lol_3_n_valid_tracks=-1;
		lol_3_min_angle=-1;
		lol_3_vtx_n_segs=-1;
		lol_3_energy=-1;
		lol_3_shower_main_length=-1;
		lol_3_n_out=-1;
		lol_3_n_sum=-1;
		cosmict_flag_1=-1; // fiducial volume vertex
		cosmict_flag_2=-1;  // single muon
		cosmict_flag_3=-1;  // single muon (long)
		cosmict_flag_4=-1;  // kinematics muon
		cosmict_flag_5=-1; // kinematics muon (long)
		cosmict_flag_6=-1; // special ...
		cosmict_flag_7=-1;  // muon+ michel
		cosmict_flag_8=-1;  // muon + michel + special
		cosmict_flag_9=-1;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
		cosmict_flag_10=nullptr;  // front upstream (dirt)
		cosmict_flag=-1;
		cosmict_2_filled=-1;
		cosmict_2_particle_type=-1;
		cosmict_2_n_muon_tracks=-1;
		cosmict_2_total_shower_length=-1;
		cosmict_2_flag_inside=-1;
		cosmict_2_angle_beam=-1;
		cosmict_2_flag_dir_weak=-1;
		cosmict_2_dQ_dx_end=-1;
		cosmict_2_dQ_dx_front=-1;
		cosmict_2_theta=-1;
		cosmict_2_phi=-1;
		cosmict_2_valid_tracks=-1;
		cosmict_3_filled=-1;
		cosmict_3_flag_inside=-1;
		cosmict_3_angle_beam=-1;
		cosmict_3_flag_dir_weak=-1;
		cosmict_3_dQ_dx_end=-1;
		cosmict_3_dQ_dx_front=-1;
		cosmict_3_theta=-1;
		cosmict_3_phi=-1;
		cosmict_3_valid_tracks=-1;
		cosmict_4_filled=-1;
		cosmict_4_flag_inside=-1;
		cosmict_4_angle_beam=-1;
		cosmict_4_connected_showers=-1;  // need to be careful about the nueCC ...
		cosmict_5_filled=-1;
		cosmict_5_flag_inside=-1;
		cosmict_5_angle_beam=-1;
		cosmict_5_connected_showers=-1;
		cosmict_6_filled=-1;
		cosmict_6_flag_dir_weak=-1;
		cosmict_6_flag_inside=-1;
		cosmict_6_angle=-1;
		cosmict_7_filled=-1;
		cosmict_7_flag_sec=-1;
		cosmict_7_n_muon_tracks=-1;
		cosmict_7_total_shower_length=-1;
		cosmict_7_flag_inside=-1;
		cosmict_7_angle_beam=-1;
		cosmict_7_flag_dir_weak=-1;
		cosmict_7_dQ_dx_end=-1;
		cosmict_7_dQ_dx_front=-1;
		cosmict_7_theta=-1;
		cosmict_7_phi=-1;
		cosmict_8_filled=-1;
		cosmict_8_flag_out=-1;
		cosmict_8_muon_length=-1;
		cosmict_8_acc_length=-1;
		cosmict_10_flag_inside=nullptr;
		cosmict_10_vtx_z=nullptr;
		cosmict_10_flag_shower=nullptr;
		cosmict_10_flag_dir_weak=nullptr;
		cosmict_10_angle_beam=nullptr;
		cosmict_10_length=nullptr;
		numu_cc_flag=-1;
		numu_cc_flag_1=nullptr;
		numu_cc_1_particle_type=nullptr;
		numu_cc_1_length=nullptr;
		numu_cc_1_medium_dQ_dx=nullptr;
		numu_cc_1_dQ_dx_cut=nullptr;
		numu_cc_1_direct_length=nullptr;
		numu_cc_1_n_daughter_tracks=nullptr;
		numu_cc_1_n_daughter_all=nullptr;
		numu_cc_flag_2=nullptr;
		numu_cc_2_length=nullptr;
		numu_cc_2_total_length=nullptr;
		numu_cc_2_n_daughter_tracks=nullptr;
		numu_cc_2_n_daughter_all=nullptr;
		numu_cc_flag_3=-1;
		numu_cc_3_particle_type=-1;
		numu_cc_3_max_length=-1;
		numu_cc_3_acc_track_length=-1;
		numu_cc_3_max_length_all=-1;
		numu_cc_3_max_muon_length=-1;
		numu_cc_3_n_daughter_tracks=-1;
		numu_cc_3_n_daughter_all=-1;
		cosmict_2_4_score=-1;
		cosmict_3_5_score=-1;
		cosmict_6_score=-1;
		cosmict_7_score=-1;
		cosmict_8_score=-1;
		cosmict_10_score=-1;
		numu_1_score=-1;
		numu_2_score=-1;
		numu_3_score=-1;
		cosmict_score=-1;
		numu_score=-1;
		mipid_score=-1;
		gap_score=-1;
		hol_lol_score=-1;
		cme_anc_score=-1;
		mgo_mgt_score=-1;
		br1_score=-1;
		br3_score=-1;
		br3_3_score=-1;
		br3_5_score=-1;
		br3_6_score=-1;
		stemdir_br2_score=-1;
		trimuon_score=-1;
		br4_tro_score=-1;
		mipquality_score=-1;
		pio_1_score=-1;
		pio_2_score=-1;
		stw_spt_score=-1;
		vis_1_score=-1;
		vis_2_score=-1;
		stw_2_score=-1;
		stw_3_score=-1;
		stw_4_score=-1;
		sig_1_score=-1;
		sig_2_score=-1;
		lol_1_score=-1;
		lol_2_score=-1;
		tro_1_score=-1;
		tro_2_score=-1;
		tro_4_score=-1;
		tro_5_score=-1;
		nue_score=-1;
	}

	if(f_KINEvars){
		  kine_reco_Enu=-1;
		  kine_reco_add_energy=-1;
		  kine_energy_particle=nullptr;
		  kine_energy_info=nullptr;
		  kine_particle_type=nullptr;
		  kine_energy_included=nullptr;
		  kine_pio_mass=-1;
		  kine_pio_flag=-1;
		  kine_pio_vtx_dis=-1;
		  kine_pio_energy_1=-1;
		  kine_pio_theta_1=-1;
		  kine_pio_phi_1=-1;
		  kine_pio_dis_1=-1;
		  kine_pio_energy_2=-1;
		  kine_pio_theta_2=-1;
		  kine_pio_phi_2=-1;
		  kine_pio_dis_2=-1;
		  kine_pio_angle=-1;
	}

	if(f_PFDump){
          truth_Ntrack = 0;
          if (truth_process) truth_process->clear(); // not nullptr
          if (truth_daughters) truth_daughters->clear();
          reco_Ntrack = 0;
          reco_process->clear();
          reco_daughters->clear();
          if (fMC_trackPosition) fMC_trackPosition->Clear();

          mc_isnu = 0;
          mc_nGeniePrimaries = -1;
          mc_nu_pdg = -1;
          mc_nu_ccnc = -1;
          mc_nu_mode = -1;
          mc_nu_intType = -1;
          mc_nu_target = -1;
          mc_hitnuc = -1;
          mc_hitquark = -1;
          mc_nu_Q2 = -1;
          mc_nu_W = -1;
          mc_nu_X = -1;
          mc_nu_Y = -1;
          mc_nu_Pt = -1;
          mc_nu_Theta = -1;
          for (int i=0; i<4; i++) {
              mc_nu_pos[i] = 0;
              mc_nu_mom[i] = 0;
          }
	}

  if (f_savesps){
    f_sps_x->clear();
    f_sps_y->clear();
    f_sps_z->clear();
    f_sps_e->clear();
    f_sps_pdg->clear();
    f_sps_id->clear();
  }

  if (f_savepmt){
    f_PMT_ID->clear();
    f_PMT_Time->clear();
    f_PMT_Amp->clear();
    f_PMT_TimeProp->clear();
    f_PMT_TimeDP->clear();
    f_PMT_TimeDL->clear();
    f_PMT_Sat->clear();
    f_RWM_Time = -999;
  }

	f_neutrino_type = -1;
	f_nuvtx_diff = -1;
	f_showervtx_diff = -1;
	f_muonvtx_diff = -1;
	f_reco_nuvtxX = -1;
	f_reco_nuvtxY = -1;
	f_reco_nuvtxZ = -1;
  f_reco_vec_showervtxX->clear(); // all reco showers
	f_reco_vec_showervtxY->clear();
	f_reco_vec_showervtxZ->clear();
	f_reco_vec_showerKE->clear();
	f_reco_showervtxX = -1; // primary shower [highest energy electron & not pi0 daughters]
	f_reco_showervtxY = -1;
	f_reco_showervtxZ = -1;
	f_reco_showerKE = -1;
	f_reco_muonvtxX = -1; //
	f_reco_muonvtxY = -1;
	f_reco_muonvtxZ = -1;
	f_reco_Nproton = -1;
	f_reco_protonvtxX = -1; //
	f_reco_protonvtxY = -1;
	f_reco_protonvtxZ = -1;
	f_reco_showerMomentum[0] = -1;
	f_reco_showerMomentum[1] = -1;
	f_reco_showerMomentum[2] = -1;
	f_reco_showerMomentum[3] = -1;
	f_reco_muonMomentum[0] = -1;
	f_reco_muonMomentum[1] = -1;
	f_reco_muonMomentum[2] = -1;
	f_reco_muonMomentum[3] = -1;
	f_reco_protonMomentum[0] = -1;
	f_reco_protonMomentum[1] = -1;
	f_reco_protonMomentum[2] = -1;
	f_reco_protonMomentum[3] = -1;

  f_evtDeltaTimeNS = -99999.;
  f_evtTimeNS = -99999.;
  f_evtTimeNS_redk2nu = -99999.;
  f_Ph_Tot = -99999.;
  for(int i=0;i<32;i++){calib[i]=0;}

	f_mcflux_run = -1;
	f_mcflux_evtno = -1;
	f_mcflux_ndecay = -1;
	f_mcflux_ntype = -1;
	f_mcflux_ptype = -1;
	f_mcflux_tptype = -1;
	f_mcflux_nuEnergy = -1;
	f_mcflux_vx = -1;
	f_mcflux_vy = -1;
	f_mcflux_vz = -1;
	f_mcflux_genx = -1;
	f_mcflux_geny = -1;
	f_mcflux_genz = -1;
	f_mcflux_dk2gen = -1;
	f_mcflux_gen2vtx = -1;

	f_truth_corr_nuvtxX = -1; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
	f_truth_corr_nuvtxY = -1;
	f_truth_corr_nuvtxZ = -1;
	f_truth_corr_showervtxX = -1; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
	f_truth_corr_showervtxY = -1;
	f_truth_corr_showervtxZ = -1;
	f_truth_showerKE = -1;
	f_truth_showerMomentum[0] = -1;
	f_truth_showerMomentum[1] = -1;
	f_truth_showerMomentum[2] = -1;
	f_truth_showerMomentum[3] = -1;
  f_truth_showerPdg = -1;
  f_truth_showerMother = -1.;
	f_truth_corr_muonvtxX = -1; //
	f_truth_corr_muonvtxY = -1;
	f_truth_corr_muonvtxZ = -1;
	f_truth_muonvtxX = -1; //
	f_truth_muonvtxY = -1;
	f_truth_muonvtxZ = -1;
	f_truth_muonendX = -1; //
	f_truth_muonendY = -1;
	f_truth_muonendZ = -1;
	f_truth_muonMomentum[0] = -1;
	f_truth_muonMomentum[1] = -1;
	f_truth_muonMomentum[2] = -1;
	f_truth_muonMomentum[3] = -1;
	f_truth_nuIntType = -1;
	f_truth_nuScatType = -1;
  f_truth_Npi0 = 0;
        f_truth_NprimPio = 0;
        f_truth_pio_energy_1 = -1;
        f_truth_pio_energy_2 = -1;
        f_truth_pio_angle = -1;
        f_truth_NCDelta = -1;
  f_truth_single_photon = -1;
  f_truth_photon_angle = -1;
  f_truth_photon_dis = -1;
	f_truth_nu_pos[0] = -1;
	f_truth_nu_pos[1] = -1;
	f_truth_nu_pos[2] = -1;
	f_truth_nu_pos[3] = -1;
	f_truth_nu_momentum[0] = -1;
	f_truth_nu_momentum[1] = -1;
	f_truth_nu_momentum[2] = -1;
	f_truth_nu_momentum[3] = -1;

	f_redk2nu_time = -1;
	f_redk2nu_time_nospill = -1;
	f_redk2nu_deltatime = 0;

	fPrimaryID.clear();
	fShowerID.clear();
	fMuonID.clear();
	fProtonID.clear();
	fPrimaryID.shrink_to_fit();
	fShowerID.shrink_to_fit();
	fMuonID.shrink_to_fit();
	fProtonID.shrink_to_fit();
	fParticleMap.clear();

 	f_weight_spline = -1.0;
	f_weight_cv = -1.0;
	f_weight_lee = -1.0;

	f_stm_eventtype = -1;
	f_stm_lowenergy = -1;
	f_stm_LM = -1;
	f_stm_TGM = -1;
	f_stm_STM = -1;
	f_stm_FullDead = -1;
	f_stm_clusterlength = -1;
}

void WireCellAnaTree::save_weights(art::Event const& e)
{
  double ppfx_cv_UBPPFXCV = 1.0; // for NuMI

  // Use the EventWeight producer label here
  art::Handle<std::vector<evwgh::MCEventWeight> > weightsHandle;
  // e.getByLabel("eventweight", weightsHandle);
  e.getByLabel(fWeightLabel, weightsHandle);

  // Loop through these objects for each neutrino vertex in the event
  for(size_t i=0; i<weightsHandle->size(); i++){
    const evwgh::MCEventWeight& mc_weights = weightsHandle->at(i);
    // Loop over all of the weights in the MCEventWeight object
    for ( const auto& pair : mc_weights.fWeight ) {
      std::string knob_name = pair.first;
      std::vector<double> weights = pair.second;
      //std::cout<<"Knob name: "<<knob_name<<std::endl;
      //std::cout<<"Weight size: "<<weights.size()<<std::endl;

      if( knob_name == "TunedCentralValue_UBGenie"){
          f_weight_cv = weights.at(0);
          if (std::isnan(f_weight_cv) or std::isinf(f_weight_cv)) {
            f_weight_cv = 1.0;
          }
      }
      if (knob_name == "splines_general_Spline"){
          f_weight_spline = weights.at(0);
          if (std::isnan(f_weight_spline) or std::isinf(f_weight_spline)) {
            f_weight_spline = 1.0;
          }
      }
      if (knob_name == "ppfx_cv_UBPPFXCV" and fIsNuMI){
          double value = weights.at(0);
          if (not std::isnan(value) and not std::isinf(value)) {
            ppfx_cv_UBPPFXCV = value;
          }
      }
      if (knob_name == "ppfx_oldrw_cv_UBOLDPPFXCV" and fIsNuMI and fNuMIOldReweight){
          double value = weights.at(0);
          if (not std::isnan(value) and not std::isinf(value)) {
            ppfx_cv_UBPPFXCV = value;
          }
      }
    }
  }

  if (fIsNuMI) {
    f_weight_spline *= ppfx_cv_UBPPFXCV; // absorb NuMI's cv correction into spline
    // std::cout << "weight_spline *= ppfx_cv_UBPPFXCV, where ppfx_cv_UBPPFXCV= " << ppfx_cv_UBPPFXCV << std::endl;
  }

  //std::cout<<"cv weight: "<<f_weight_cv<<std::endl;
  //std::cout<<"spline weight: "<<f_weight_spline<<std::endl;

}

void WireCellAnaTree::save_LEEweights(art::Event const& e)
{
  // Use the EventWeight producer label here
  art::Handle<std::vector<evwgh::MCEventWeight> > weightsHandle;
  // e.getByLabel("eventweightLEE", "", "EventWeightLEE", weightsHandle); // producer, instance, process
  e.getByLabel(fWeightLeeLabel, weightsHandle); // producer, instance, process

  // Loop through these objects for each neutrino vertex in the event
  for(size_t i=0; i<weightsHandle->size(); i++){
    const evwgh::MCEventWeight& mc_weights = weightsHandle->at(i);
    // Loop over all of the weights in the MCEventWeight object
    for ( const auto& pair : mc_weights.fWeight ) {
      std::string knob_name = pair.first;
      std::vector<double> weights = pair.second;
      //std::cout<<"Knob name: "<<knob_name<<std::endl;
      //std::cout<<"Weight size: "<<weights.size()<<std::endl;
      if( knob_name == "eLEE_Combined_Oct2018_LEESignalElectron" ){
          f_weight_lee = weights.at(0);
      }
    }
  }
}


void WireCellAnaTree::ReadSSMBDTvar(nsm::NuSelectionBDT const& bdt)
{
  
  ssm_flag_st_kdar = bdt.Getstkdar().ssm_flag_st_kdar;
  ssm_Nsm = bdt.Getstkdar().ssm_Nsm;
  ssm_Nsm_wivtx = bdt.Getstkdar().ssm_Nsm_wivtx;

  ssm_dq_dx_fwd_1 = bdt.Getstkdar().ssm_dq_dx_fwd_1;
  ssm_dq_dx_fwd_2 = bdt.Getstkdar().ssm_dq_dx_fwd_2;
  ssm_dq_dx_fwd_3 = bdt.Getstkdar().ssm_dq_dx_fwd_3;
  ssm_dq_dx_fwd_4 = bdt.Getstkdar().ssm_dq_dx_fwd_4;
  ssm_dq_dx_fwd_5 = bdt.Getstkdar().ssm_dq_dx_fwd_5;
  ssm_dq_dx_bck_1 = bdt.Getstkdar().ssm_dq_dx_bck_1;
  ssm_dq_dx_bck_2 = bdt.Getstkdar().ssm_dq_dx_bck_2;
  ssm_dq_dx_bck_3 = bdt.Getstkdar().ssm_dq_dx_bck_3;
  ssm_dq_dx_bck_4 = bdt.Getstkdar().ssm_dq_dx_bck_4;
  ssm_dq_dx_bck_5 = bdt.Getstkdar().ssm_dq_dx_bck_5;
  ssm_d_dq_dx_fwd_12 = bdt.Getstkdar().ssm_d_dq_dx_fwd_12;
  ssm_d_dq_dx_fwd_23 = bdt.Getstkdar().ssm_d_dq_dx_fwd_23;
  ssm_d_dq_dx_fwd_34 = bdt.Getstkdar().ssm_d_dq_dx_fwd_34;
  ssm_d_dq_dx_fwd_45 = bdt.Getstkdar().ssm_d_dq_dx_fwd_45;
  ssm_d_dq_dx_bck_12 = bdt.Getstkdar().ssm_d_dq_dx_bck_12;
  ssm_d_dq_dx_bck_23 = bdt.Getstkdar().ssm_d_dq_dx_bck_23;
  ssm_d_dq_dx_bck_34 = bdt.Getstkdar().ssm_d_dq_dx_bck_34;
  ssm_d_dq_dx_bck_45 = bdt.Getstkdar().ssm_d_dq_dx_bck_45;
  ssm_max_dq_dx_fwd_3 = bdt.Getstkdar().ssm_max_dq_dx_fwd_3;
  ssm_max_dq_dx_fwd_5 = bdt.Getstkdar().ssm_max_dq_dx_fwd_5;
  ssm_max_dq_dx_bck_3 = bdt.Getstkdar().ssm_max_dq_dx_bck_3;
  ssm_max_dq_dx_bck_5 = bdt.Getstkdar().ssm_max_dq_dx_bck_5;
  ssm_max_d_dq_dx_fwd_3 = bdt.Getstkdar().ssm_max_d_dq_dx_fwd_3;
  ssm_max_d_dq_dx_fwd_5 = bdt.Getstkdar().ssm_max_d_dq_dx_fwd_5;
  ssm_max_d_dq_dx_bck_3 = bdt.Getstkdar().ssm_max_d_dq_dx_bck_3;
  ssm_max_d_dq_dx_bck_5 = bdt.Getstkdar().ssm_max_d_dq_dx_bck_5;
  ssm_medium_dq_dx = bdt.Getstkdar().ssm_medium_dq_dx;
  ssm_medium_dq_dx_bp = bdt.Getstkdar().ssm_medium_dq_dx_bp;
      //angluar info
  ssm_angle_to_z = bdt.Getstkdar().ssm_angle_to_z;
  ssm_angle_to_target = bdt.Getstkdar().ssm_angle_to_target;
  ssm_angle_to_absorber = bdt.Getstkdar().ssm_angle_to_absorber;
  ssm_angle_to_vertical = bdt.Getstkdar().ssm_angle_to_vertical;
      //directional info
  ssm_x_dir = bdt.Getstkdar().ssm_x_dir;
  ssm_y_dir = bdt.Getstkdar().ssm_y_dir;
  ssm_z_dir = bdt.Getstkdar().ssm_z_dir;
      //energy info
  ssm_kine_energy = bdt.Getstkdar().ssm_kine_energy;
  ssm_kine_energy_reduced = bdt.Getstkdar().ssm_kine_energy_reduced;
      //general properties
  ssm_vtx_activity = bdt.Getstkdar().ssm_vtx_activity;
  ssm_pdg = bdt.Getstkdar().ssm_pdg;
  ssm_dQ_dx_cut = bdt.Getstkdar().ssm_dQ_dx_cut;
  ssm_score_mu_fwd = bdt.Getstkdar().ssm_score_mu_fwd;
  ssm_score_p_fwd = bdt.Getstkdar().ssm_score_p_fwd;
  ssm_score_e_fwd = bdt.Getstkdar().ssm_score_e_fwd;
  ssm_score_mu_bck = bdt.Getstkdar().ssm_score_mu_bck;
  ssm_score_p_bck = bdt.Getstkdar().ssm_score_p_bck;
  ssm_score_e_bck = bdt.Getstkdar().ssm_score_e_bck;
  ssm_score_mu_fwd_bp = bdt.Getstkdar().ssm_score_mu_fwd_bp;
  ssm_score_p_fwd_bp = bdt.Getstkdar().ssm_score_p_fwd_bp;
  ssm_score_e_fwd_bp = bdt.Getstkdar().ssm_score_e_fwd_bp;
      //track "straighness"
  ssm_length = bdt.Getstkdar().ssm_length;
  ssm_direct_length = bdt.Getstkdar().ssm_direct_length;
  ssm_length_ratio = bdt.Getstkdar().ssm_length_ratio;
  ssm_max_dev = bdt.Getstkdar().ssm_max_dev;
    //number of other particles
  ssm_n_prim_tracks_1 = bdt.Getstkdar().ssm_n_prim_tracks_1;
  ssm_n_prim_tracks_3 = bdt.Getstkdar().ssm_n_prim_tracks_3;
  ssm_n_prim_tracks_5 = bdt.Getstkdar().ssm_n_prim_tracks_5;
  ssm_n_prim_tracks_8 = bdt.Getstkdar().ssm_n_prim_tracks_8;
  ssm_n_prim_tracks_11 = bdt.Getstkdar().ssm_n_prim_tracks_11;
  ssm_n_all_tracks_1 = bdt.Getstkdar().ssm_n_all_tracks_1;
  ssm_n_all_tracks_3 = bdt.Getstkdar().ssm_n_all_tracks_3;
  ssm_n_all_tracks_5 = bdt.Getstkdar().ssm_n_all_tracks_5;
  ssm_n_all_tracks_8 = bdt.Getstkdar().ssm_n_all_tracks_8;
  ssm_n_all_tracks_11 = bdt.Getstkdar().ssm_n_all_tracks_11;
  ssm_n_daughter_tracks_1 = bdt.Getstkdar().ssm_n_daughter_tracks_1;
  ssm_n_daughter_tracks_3 = bdt.Getstkdar().ssm_n_daughter_tracks_3;
  ssm_n_daughter_tracks_5 = bdt.Getstkdar().ssm_n_daughter_tracks_5;
  ssm_n_daughter_tracks_8 = bdt.Getstkdar().ssm_n_daughter_tracks_8;
  ssm_n_daughter_tracks_11 = bdt.Getstkdar().ssm_n_daughter_tracks_11;
  ssm_n_daughter_all_1 = bdt.Getstkdar().ssm_n_daughter_all_1;
  ssm_n_daughter_all_3 = bdt.Getstkdar().ssm_n_daughter_all_3;
  ssm_n_daughter_all_5 = bdt.Getstkdar().ssm_n_daughter_all_5;
  ssm_n_daughter_all_8 = bdt.Getstkdar().ssm_n_daughter_all_8;
  ssm_n_daughter_all_11 = bdt.Getstkdar().ssm_n_daughter_all_11;
    //properties of leading other primary track
  ssm_prim_track1_pdg = bdt.Getstkdar().ssm_prim_track1_pdg;
  ssm_prim_track1_score_mu_fwd = bdt.Getstkdar().ssm_prim_track1_score_mu_fwd;
  ssm_prim_track1_score_p_fwd = bdt.Getstkdar().ssm_prim_track1_score_p_fwd;
  ssm_prim_track1_score_e_fwd = bdt.Getstkdar().ssm_prim_track1_score_e_fwd;
  ssm_prim_track1_score_mu_bck = bdt.Getstkdar().ssm_prim_track1_score_mu_bck;
  ssm_prim_track1_score_p_bck = bdt.Getstkdar().ssm_prim_track1_score_p_bck;
  ssm_prim_track1_score_e_bck = bdt.Getstkdar().ssm_prim_track1_score_e_bck;
  ssm_prim_track1_length = bdt.Getstkdar().ssm_prim_track1_length;
  ssm_prim_track1_direct_length = bdt.Getstkdar().ssm_prim_track1_direct_length;
  ssm_prim_track1_length_ratio = bdt.Getstkdar().ssm_prim_track1_length_ratio;
  ssm_prim_track1_max_dev = bdt.Getstkdar().ssm_prim_track1_max_dev;
  ssm_prim_track1_kine_energy_range = bdt.Getstkdar().ssm_prim_track1_kine_energy_range;
  ssm_prim_track1_kine_energy_range_mu = bdt.Getstkdar().ssm_prim_track1_kine_energy_range_mu;
  ssm_prim_track1_kine_energy_range_p = bdt.Getstkdar().ssm_prim_track1_kine_energy_range_p;
  ssm_prim_track1_kine_energy_range_e = bdt.Getstkdar().ssm_prim_track1_kine_energy_range_e;
  ssm_prim_track1_kine_energy_cal = bdt.Getstkdar().ssm_prim_track1_kine_energy_cal;
  ssm_prim_track1_medium_dq_dx = bdt.Getstkdar().ssm_prim_track1_medium_dq_dx;
  ssm_prim_track1_x_dir = bdt.Getstkdar().ssm_prim_track1_x_dir;
  ssm_prim_track1_y_dir = bdt.Getstkdar().ssm_prim_track1_y_dir;
  ssm_prim_track1_z_dir = bdt.Getstkdar().ssm_prim_track1_z_dir;
  ssm_prim_track1_add_daught_track_counts_1 = bdt.Getstkdar().ssm_prim_track1_add_daught_track_counts_1;
  ssm_prim_track1_add_daught_all_counts_1 = bdt.Getstkdar().ssm_prim_track1_add_daught_all_counts_1;
  ssm_prim_track1_add_daught_track_counts_5 = bdt.Getstkdar().ssm_prim_track1_add_daught_track_counts_5;
  ssm_prim_track1_add_daught_all_counts_5 = bdt.Getstkdar().ssm_prim_track1_add_daught_all_counts_5;
  ssm_prim_track1_add_daught_track_counts_11 = bdt.Getstkdar().ssm_prim_track1_add_daught_track_counts_11;
  ssm_prim_track1_add_daught_all_counts_11 = bdt.Getstkdar().ssm_prim_track1_add_daught_all_counts_11;
  //properties of sub-leading other primary track
  ssm_prim_track2_pdg = bdt.Getstkdar().ssm_prim_track2_pdg;
  ssm_prim_track2_score_mu_fwd = bdt.Getstkdar().ssm_prim_track2_score_mu_fwd;
  ssm_prim_track2_score_p_fwd = bdt.Getstkdar().ssm_prim_track2_score_p_fwd;
  ssm_prim_track2_score_e_fwd = bdt.Getstkdar().ssm_prim_track2_score_e_fwd;
  ssm_prim_track2_score_mu_bck = bdt.Getstkdar().ssm_prim_track2_score_mu_bck;
  ssm_prim_track2_score_p_bck = bdt.Getstkdar().ssm_prim_track2_score_p_bck;
  ssm_prim_track2_score_e_bck = bdt.Getstkdar().ssm_prim_track2_score_e_bck;
  ssm_prim_track2_length = bdt.Getstkdar().ssm_prim_track2_length;
  ssm_prim_track2_direct_length = bdt.Getstkdar().ssm_prim_track2_direct_length;
  ssm_prim_track2_length_ratio = bdt.Getstkdar().ssm_prim_track2_length_ratio;
  ssm_prim_track2_max_dev = bdt.Getstkdar().ssm_prim_track2_max_dev;
  ssm_prim_track2_kine_energy_range = bdt.Getstkdar().ssm_prim_track2_kine_energy_range;
  ssm_prim_track2_kine_energy_range_mu = bdt.Getstkdar().ssm_prim_track2_kine_energy_range_mu;
  ssm_prim_track2_kine_energy_range_p = bdt.Getstkdar().ssm_prim_track2_kine_energy_range_p;
  ssm_prim_track2_kine_energy_range_e = bdt.Getstkdar().ssm_prim_track2_kine_energy_range_e;
  ssm_prim_track2_kine_energy_cal = bdt.Getstkdar().ssm_prim_track2_kine_energy_cal;
  ssm_prim_track2_medium_dq_dx = bdt.Getstkdar().ssm_prim_track2_medium_dq_dx;
  ssm_prim_track2_x_dir = bdt.Getstkdar().ssm_prim_track2_x_dir;
  ssm_prim_track2_y_dir = bdt.Getstkdar().ssm_prim_track2_y_dir;
  ssm_prim_track2_z_dir = bdt.Getstkdar().ssm_prim_track2_z_dir;
  ssm_prim_track2_add_daught_track_counts_1 = bdt.Getstkdar().ssm_prim_track2_add_daught_track_counts_1;
  ssm_prim_track2_add_daught_all_counts_1 = bdt.Getstkdar().ssm_prim_track2_add_daught_all_counts_1;
  ssm_prim_track2_add_daught_track_counts_5 = bdt.Getstkdar().ssm_prim_track2_add_daught_track_counts_5;
  ssm_prim_track2_add_daught_all_counts_5 = bdt.Getstkdar().ssm_prim_track2_add_daught_all_counts_5;
  ssm_prim_track2_add_daught_track_counts_11 = bdt.Getstkdar().ssm_prim_track2_add_daught_track_counts_11;
  ssm_prim_track2_add_daught_all_counts_11 = bdt.Getstkdar().ssm_prim_track2_add_daught_all_counts_11;
  //properties of leading daughter track
  ssm_daught_track1_pdg = bdt.Getstkdar().ssm_daught_track1_pdg;
  ssm_daught_track1_score_mu_fwd = bdt.Getstkdar().ssm_daught_track1_score_mu_fwd;
  ssm_daught_track1_score_p_fwd = bdt.Getstkdar().ssm_daught_track1_score_p_fwd;
  ssm_daught_track1_score_e_fwd = bdt.Getstkdar().ssm_daught_track1_score_e_fwd;
  ssm_daught_track1_score_mu_bck = bdt.Getstkdar().ssm_daught_track1_score_mu_bck;
  ssm_daught_track1_score_p_bck = bdt.Getstkdar().ssm_daught_track1_score_p_bck;
  ssm_daught_track1_score_e_bck = bdt.Getstkdar().ssm_daught_track1_score_e_bck;
  ssm_daught_track1_length = bdt.Getstkdar().ssm_daught_track1_length;
  ssm_daught_track1_direct_length = bdt.Getstkdar().ssm_daught_track1_direct_length;
  ssm_daught_track1_length_ratio = bdt.Getstkdar().ssm_daught_track1_length_ratio;
  ssm_daught_track1_max_dev = bdt.Getstkdar().ssm_daught_track1_max_dev;
  ssm_daught_track1_kine_energy_range = bdt.Getstkdar().ssm_daught_track1_kine_energy_range;
  ssm_daught_track1_kine_energy_range_mu = bdt.Getstkdar().ssm_daught_track1_kine_energy_range_mu;
  ssm_daught_track1_kine_energy_range_p = bdt.Getstkdar().ssm_daught_track1_kine_energy_range_p;
  ssm_daught_track1_kine_energy_range_e = bdt.Getstkdar().ssm_daught_track1_kine_energy_range_e;
  ssm_daught_track1_kine_energy_cal = bdt.Getstkdar().ssm_daught_track1_kine_energy_cal;
  ssm_daught_track1_medium_dq_dx = bdt.Getstkdar().ssm_daught_track1_medium_dq_dx;
  ssm_daught_track1_x_dir = bdt.Getstkdar().ssm_daught_track1_x_dir;
  ssm_daught_track1_y_dir = bdt.Getstkdar().ssm_daught_track1_y_dir;
  ssm_daught_track1_z_dir = bdt.Getstkdar().ssm_daught_track1_z_dir;
  ssm_daught_track1_add_daught_track_counts_1 = bdt.Getstkdar().ssm_daught_track1_add_daught_track_counts_1;
  ssm_daught_track1_add_daught_all_counts_1 = bdt.Getstkdar().ssm_daught_track1_add_daught_all_counts_1;
  ssm_daught_track1_add_daught_track_counts_5 = bdt.Getstkdar().ssm_daught_track1_add_daught_track_counts_5;
  ssm_daught_track1_add_daught_all_counts_5 = bdt.Getstkdar().ssm_daught_track1_add_daught_all_counts_5;
  ssm_daught_track1_add_daught_track_counts_11 = bdt.Getstkdar().ssm_daught_track1_add_daught_track_counts_11;
  ssm_daught_track1_add_daught_all_counts_11 = bdt.Getstkdar().ssm_daught_track1_add_daught_all_counts_11;
  //properties of sub-leading daughter track
  ssm_daught_track2_pdg = bdt.Getstkdar().ssm_daught_track2_pdg;
  ssm_daught_track2_score_mu_fwd = bdt.Getstkdar().ssm_daught_track2_score_mu_fwd;
  ssm_daught_track2_score_p_fwd = bdt.Getstkdar().ssm_daught_track2_score_p_fwd;
  ssm_daught_track2_score_e_fwd = bdt.Getstkdar().ssm_daught_track2_score_e_fwd;
  ssm_daught_track2_score_mu_bck = bdt.Getstkdar().ssm_daught_track2_score_mu_bck;
  ssm_daught_track2_score_p_bck = bdt.Getstkdar().ssm_daught_track2_score_p_bck;
  ssm_daught_track2_score_e_bck = bdt.Getstkdar().ssm_daught_track2_score_e_bck;
  ssm_daught_track2_length = bdt.Getstkdar().ssm_daught_track2_length;
  ssm_daught_track2_direct_length = bdt.Getstkdar().ssm_daught_track2_direct_length;
  ssm_daught_track2_length_ratio = bdt.Getstkdar().ssm_daught_track2_length_ratio;
  ssm_daught_track2_max_dev = bdt.Getstkdar().ssm_daught_track2_max_dev;
  ssm_daught_track2_kine_energy_range = bdt.Getstkdar().ssm_daught_track2_kine_energy_range;
  ssm_daught_track2_kine_energy_range_mu = bdt.Getstkdar().ssm_daught_track2_kine_energy_range_mu;
  ssm_daught_track2_kine_energy_range_p = bdt.Getstkdar().ssm_daught_track2_kine_energy_range_p;
  ssm_daught_track2_kine_energy_range_e = bdt.Getstkdar().ssm_daught_track2_kine_energy_range_e;
  ssm_daught_track2_kine_energy_cal = bdt.Getstkdar().ssm_daught_track2_kine_energy_cal;
  ssm_daught_track2_medium_dq_dx = bdt.Getstkdar().ssm_daught_track2_medium_dq_dx;
  ssm_daught_track2_x_dir = bdt.Getstkdar().ssm_daught_track2_x_dir;
  ssm_daught_track2_y_dir = bdt.Getstkdar().ssm_daught_track2_y_dir;
  ssm_daught_track2_z_dir = bdt.Getstkdar().ssm_daught_track2_z_dir;
  ssm_daught_track2_add_daught_track_counts_1 = bdt.Getstkdar().ssm_daught_track2_add_daught_track_counts_1;
  ssm_daught_track2_add_daught_all_counts_1 = bdt.Getstkdar().ssm_daught_track2_add_daught_all_counts_1;
  ssm_daught_track2_add_daught_track_counts_5 = bdt.Getstkdar().ssm_daught_track2_add_daught_track_counts_5;
  ssm_daught_track2_add_daught_all_counts_5 = bdt.Getstkdar().ssm_daught_track2_add_daught_all_counts_5;
  ssm_daught_track2_add_daught_track_counts_11 = bdt.Getstkdar().ssm_daught_track2_add_daught_track_counts_11;
  ssm_daught_track2_add_daught_all_counts_11 = bdt.Getstkdar().ssm_daught_track2_add_daught_all_counts_11;
  //properties of leading other primary shower
  ssm_prim_shw1_pdg = bdt.Getstkdar().ssm_prim_shw1_pdg;
  ssm_prim_shw1_score_mu_fwd = bdt.Getstkdar().ssm_prim_shw1_score_mu_fwd;
  ssm_prim_shw1_score_p_fwd = bdt.Getstkdar().ssm_prim_shw1_score_p_fwd;
  ssm_prim_shw1_score_e_fwd = bdt.Getstkdar().ssm_prim_shw1_score_e_fwd;
  ssm_prim_shw1_score_mu_bck = bdt.Getstkdar().ssm_prim_shw1_score_mu_bck;
  ssm_prim_shw1_score_p_bck = bdt.Getstkdar().ssm_prim_shw1_score_p_bck;
  ssm_prim_shw1_score_e_bck = bdt.Getstkdar().ssm_prim_shw1_score_e_bck;
  ssm_prim_shw1_length = bdt.Getstkdar().ssm_prim_shw1_length;
  ssm_prim_shw1_direct_length = bdt.Getstkdar().ssm_prim_shw1_direct_length;
  ssm_prim_shw1_length_ratio = bdt.Getstkdar().ssm_prim_shw1_length_ratio;
  ssm_prim_shw1_max_dev = bdt.Getstkdar().ssm_prim_shw1_max_dev;
  ssm_prim_shw1_kine_energy_range = bdt.Getstkdar().ssm_prim_shw1_kine_energy_range;
  ssm_prim_shw1_kine_energy_range_mu = bdt.Getstkdar().ssm_prim_shw1_kine_energy_range_mu;
  ssm_prim_shw1_kine_energy_range_p = bdt.Getstkdar().ssm_prim_shw1_kine_energy_range_p;
  ssm_prim_shw1_kine_energy_range_e = bdt.Getstkdar().ssm_prim_shw1_kine_energy_range_e;
  ssm_prim_shw1_kine_energy_cal = bdt.Getstkdar().ssm_prim_shw1_kine_energy_cal;
  ssm_prim_shw1_kine_energy_best = bdt.Getstkdar().ssm_prim_shw1_kine_energy_best;
  ssm_prim_shw1_medium_dq_dx = bdt.Getstkdar().ssm_prim_shw1_medium_dq_dx;
  ssm_prim_shw1_x_dir = bdt.Getstkdar().ssm_prim_shw1_x_dir;
  ssm_prim_shw1_y_dir = bdt.Getstkdar().ssm_prim_shw1_y_dir;
  ssm_prim_shw1_z_dir = bdt.Getstkdar().ssm_prim_shw1_z_dir;
  ssm_prim_shw1_add_daught_track_counts_1 = bdt.Getstkdar().ssm_prim_shw1_add_daught_track_counts_1;
  ssm_prim_shw1_add_daught_all_counts_1 = bdt.Getstkdar().ssm_prim_shw1_add_daught_all_counts_1;
  ssm_prim_shw1_add_daught_track_counts_5 = bdt.Getstkdar().ssm_prim_shw1_add_daught_track_counts_5;
  ssm_prim_shw1_add_daught_all_counts_5 = bdt.Getstkdar().ssm_prim_shw1_add_daught_all_counts_5;
  ssm_prim_shw1_add_daught_track_counts_11 = bdt.Getstkdar().ssm_prim_shw1_add_daught_track_counts_11;
  ssm_prim_shw1_add_daught_all_counts_11 = bdt.Getstkdar().ssm_prim_shw1_add_daught_all_counts_11;
  //properties of sub-leading other primary shower
  ssm_prim_shw2_pdg = bdt.Getstkdar().ssm_prim_shw2_pdg;
  ssm_prim_shw2_score_mu_fwd = bdt.Getstkdar().ssm_prim_shw2_score_mu_fwd;
  ssm_prim_shw2_score_p_fwd = bdt.Getstkdar().ssm_prim_shw2_score_p_fwd;
  ssm_prim_shw2_score_e_fwd = bdt.Getstkdar().ssm_prim_shw2_score_e_fwd;
  ssm_prim_shw2_score_mu_bck = bdt.Getstkdar().ssm_prim_shw2_score_mu_bck;
  ssm_prim_shw2_score_p_bck = bdt.Getstkdar().ssm_prim_shw2_score_p_bck;
  ssm_prim_shw2_score_e_bck = bdt.Getstkdar().ssm_prim_shw2_score_e_bck;
  ssm_prim_shw2_length = bdt.Getstkdar().ssm_prim_shw2_length;
  ssm_prim_shw2_direct_length = bdt.Getstkdar().ssm_prim_shw2_direct_length;
  ssm_prim_shw2_length_ratio = bdt.Getstkdar().ssm_prim_shw2_length_ratio;
  ssm_prim_shw2_max_dev = bdt.Getstkdar().ssm_prim_shw2_max_dev;
  ssm_prim_shw2_kine_energy_range = bdt.Getstkdar().ssm_prim_shw2_kine_energy_range;
  ssm_prim_shw2_kine_energy_range_mu = bdt.Getstkdar().ssm_prim_shw2_kine_energy_range_mu;
  ssm_prim_shw2_kine_energy_range_p = bdt.Getstkdar().ssm_prim_shw2_kine_energy_range_p;
  ssm_prim_shw2_kine_energy_range_e = bdt.Getstkdar().ssm_prim_shw2_kine_energy_range_e;
  ssm_prim_shw2_kine_energy_cal = bdt.Getstkdar().ssm_prim_shw2_kine_energy_cal;
  ssm_prim_shw2_kine_energy_best = bdt.Getstkdar().ssm_prim_shw2_kine_energy_best;
  ssm_prim_shw2_medium_dq_dx = bdt.Getstkdar().ssm_prim_shw2_medium_dq_dx;
  ssm_prim_shw2_x_dir = bdt.Getstkdar().ssm_prim_shw2_x_dir;
  ssm_prim_shw2_y_dir = bdt.Getstkdar().ssm_prim_shw2_y_dir;
  ssm_prim_shw2_z_dir = bdt.Getstkdar().ssm_prim_shw2_z_dir;
  ssm_prim_shw2_add_daught_track_counts_1 = bdt.Getstkdar().ssm_prim_shw2_add_daught_track_counts_1;
  ssm_prim_shw2_add_daught_all_counts_1 = bdt.Getstkdar().ssm_prim_shw2_add_daught_all_counts_1;
  ssm_prim_shw2_add_daught_track_counts_5 = bdt.Getstkdar().ssm_prim_shw2_add_daught_track_counts_5;
  ssm_prim_shw2_add_daught_all_counts_5 = bdt.Getstkdar().ssm_prim_shw2_add_daught_all_counts_5;
  ssm_prim_shw2_add_daught_track_counts_11 = bdt.Getstkdar().ssm_prim_shw2_add_daught_track_counts_11;
  ssm_prim_shw2_add_daught_all_counts_11 = bdt.Getstkdar().ssm_prim_shw2_add_daught_all_counts_11;
  //properties of leading daughter shower
  ssm_daught_shw1_pdg = bdt.Getstkdar().ssm_daught_shw1_pdg;
  ssm_daught_shw1_score_mu_fwd = bdt.Getstkdar().ssm_daught_shw1_score_mu_fwd;
  ssm_daught_shw1_score_p_fwd = bdt.Getstkdar().ssm_daught_shw1_score_p_fwd;
  ssm_daught_shw1_score_e_fwd = bdt.Getstkdar().ssm_daught_shw1_score_e_fwd;
  ssm_daught_shw1_score_mu_bck = bdt.Getstkdar().ssm_daught_shw1_score_mu_bck;
  ssm_daught_shw1_score_p_bck = bdt.Getstkdar().ssm_daught_shw1_score_p_bck;
  ssm_daught_shw1_score_e_bck = bdt.Getstkdar().ssm_daught_shw1_score_e_bck;
  ssm_daught_shw1_length = bdt.Getstkdar().ssm_daught_shw1_length;
  ssm_daught_shw1_direct_length = bdt.Getstkdar().ssm_daught_shw1_direct_length;
  ssm_daught_shw1_length_ratio = bdt.Getstkdar().ssm_daught_shw1_length_ratio;
  ssm_daught_shw1_max_dev = bdt.Getstkdar().ssm_daught_shw1_max_dev;
  ssm_daught_shw1_kine_energy_range = bdt.Getstkdar().ssm_daught_shw1_kine_energy_range;
  ssm_daught_shw1_kine_energy_range_mu = bdt.Getstkdar().ssm_daught_shw1_kine_energy_range_mu;
  ssm_daught_shw1_kine_energy_range_p = bdt.Getstkdar().ssm_daught_shw1_kine_energy_range_p;
  ssm_daught_shw1_kine_energy_range_e = bdt.Getstkdar().ssm_daught_shw1_kine_energy_range_e;
  ssm_daught_shw1_kine_energy_cal = bdt.Getstkdar().ssm_daught_shw1_kine_energy_cal;
  ssm_daught_shw1_kine_energy_best = bdt.Getstkdar().ssm_daught_shw1_kine_energy_best;
  ssm_daught_shw1_medium_dq_dx = bdt.Getstkdar().ssm_daught_shw1_medium_dq_dx;
  ssm_daught_shw1_x_dir = bdt.Getstkdar().ssm_daught_shw1_x_dir;
  ssm_daught_shw1_y_dir = bdt.Getstkdar().ssm_daught_shw1_y_dir;
  ssm_daught_shw1_z_dir = bdt.Getstkdar().ssm_daught_shw1_z_dir;
  ssm_daught_shw1_add_daught_track_counts_1 = bdt.Getstkdar().ssm_daught_shw1_add_daught_track_counts_1;
  ssm_daught_shw1_add_daught_all_counts_1 = bdt.Getstkdar().ssm_daught_shw1_add_daught_all_counts_1;
  ssm_daught_shw1_add_daught_track_counts_5 = bdt.Getstkdar().ssm_daught_shw1_add_daught_track_counts_5;
  ssm_daught_shw1_add_daught_all_counts_5 = bdt.Getstkdar().ssm_daught_shw1_add_daught_all_counts_5;
  ssm_daught_shw1_add_daught_track_counts_11 = bdt.Getstkdar().ssm_daught_shw1_add_daught_track_counts_11;
  ssm_daught_shw1_add_daught_all_counts_11 = bdt.Getstkdar().ssm_daught_shw1_add_daught_all_counts_11;
  //properties of sub-leading daughter shower
  ssm_daught_shw2_pdg = bdt.Getstkdar().ssm_daught_shw2_pdg;
  ssm_daught_shw2_score_mu_fwd = bdt.Getstkdar().ssm_daught_shw2_score_mu_fwd;
  ssm_daught_shw2_score_p_fwd = bdt.Getstkdar().ssm_daught_shw2_score_p_fwd;
  ssm_daught_shw2_score_e_fwd = bdt.Getstkdar().ssm_daught_shw2_score_e_fwd;
  ssm_daught_shw2_score_mu_bck = bdt.Getstkdar().ssm_daught_shw2_score_mu_bck;
  ssm_daught_shw2_score_p_bck = bdt.Getstkdar().ssm_daught_shw2_score_p_bck;
  ssm_daught_shw2_score_e_bck = bdt.Getstkdar().ssm_daught_shw2_score_e_bck;
  ssm_daught_shw2_length = bdt.Getstkdar().ssm_daught_shw2_length;
  ssm_daught_shw2_direct_length = bdt.Getstkdar().ssm_daught_shw2_direct_length;
  ssm_daught_shw2_length_ratio = bdt.Getstkdar().ssm_daught_shw2_length_ratio;
  ssm_daught_shw2_max_dev = bdt.Getstkdar().ssm_daught_shw2_max_dev;
  ssm_daught_shw2_kine_energy_range = bdt.Getstkdar().ssm_daught_shw2_kine_energy_range;
  ssm_daught_shw2_kine_energy_range_mu = bdt.Getstkdar().ssm_daught_shw2_kine_energy_range_mu;
  ssm_daught_shw2_kine_energy_range_p = bdt.Getstkdar().ssm_daught_shw2_kine_energy_range_p;
  ssm_daught_shw2_kine_energy_range_e = bdt.Getstkdar().ssm_daught_shw2_kine_energy_range_e;
  ssm_daught_shw2_kine_energy_cal = bdt.Getstkdar().ssm_daught_shw2_kine_energy_cal;
  ssm_daught_shw2_kine_energy_best = bdt.Getstkdar().ssm_daught_shw2_kine_energy_best;
  ssm_daught_shw2_medium_dq_dx = bdt.Getstkdar().ssm_daught_shw2_medium_dq_dx;
  ssm_daught_shw2_x_dir = bdt.Getstkdar().ssm_daught_shw2_x_dir;
  ssm_daught_shw2_y_dir = bdt.Getstkdar().ssm_daught_shw2_y_dir;
  ssm_daught_shw2_z_dir = bdt.Getstkdar().ssm_daught_shw2_z_dir;
  ssm_daught_shw2_add_daught_track_counts_1 = bdt.Getstkdar().ssm_daught_shw2_add_daught_track_counts_1;
  ssm_daught_shw2_add_daught_all_counts_1 = bdt.Getstkdar().ssm_daught_shw2_add_daught_all_counts_1;
  ssm_daught_shw2_add_daught_track_counts_5 = bdt.Getstkdar().ssm_daught_shw2_add_daught_track_counts_5;
  ssm_daught_shw2_add_daught_all_counts_5 = bdt.Getstkdar().ssm_daught_shw2_add_daught_all_counts_5;
  ssm_daught_shw2_add_daught_track_counts_11 = bdt.Getstkdar().ssm_daught_shw2_add_daught_track_counts_11;
  ssm_daught_shw2_add_daught_all_counts_11 = bdt.Getstkdar().ssm_daught_shw2_add_daught_all_counts_11;
  //event level properties
  ssm_nu_angle_z = bdt.Getstkdar().ssm_nu_angle_z;
  ssm_nu_angle_target = bdt.Getstkdar().ssm_nu_angle_target;
  ssm_nu_angle_absorber = bdt.Getstkdar().ssm_nu_angle_absorber;
  ssm_nu_angle_vertical = bdt.Getstkdar().ssm_nu_angle_vertical;
  ssm_con_nu_angle_z = bdt.Getstkdar().ssm_con_nu_angle_z;
  ssm_con_nu_angle_target = bdt.Getstkdar().ssm_con_nu_angle_target;
  ssm_con_nu_angle_absorber = bdt.Getstkdar().ssm_con_nu_angle_absorber;
  ssm_con_nu_angle_vertical = bdt.Getstkdar().ssm_con_nu_angle_vertical;
  ssm_prim_nu_angle_z = bdt.Getstkdar().ssm_prim_nu_angle_z;
  ssm_prim_nu_angle_target = bdt.Getstkdar().ssm_prim_nu_angle_target;
  ssm_prim_nu_angle_absorber = bdt.Getstkdar().ssm_prim_nu_angle_absorber;
  ssm_prim_nu_angle_vertical = bdt.Getstkdar().ssm_prim_nu_angle_vertical;
  ssm_track_angle_z = bdt.Getstkdar().ssm_track_angle_z;
  ssm_track_angle_target = bdt.Getstkdar().ssm_track_angle_target;
  ssm_track_angle_absorber = bdt.Getstkdar().ssm_track_angle_absorber;
  ssm_track_angle_vertical = bdt.Getstkdar().ssm_track_angle_vertical;
  ssm_vtxX = bdt.Getstkdar().ssm_vtxX;
  ssm_vtxY = bdt.Getstkdar().ssm_vtxY;
  ssm_vtxZ = bdt.Getstkdar().ssm_vtxZ;
  //off vertex stuff
  ssm_offvtx_length = bdt.Getstkdar().ssm_offvtx_length;
  ssm_offvtx_energy = bdt.Getstkdar().ssm_offvtx_energy;
  ssm_n_offvtx_tracks_1 = bdt.Getstkdar().ssm_n_offvtx_tracks_1;
  ssm_n_offvtx_tracks_3 = bdt.Getstkdar().ssm_n_offvtx_tracks_3;
  ssm_n_offvtx_tracks_5 = bdt.Getstkdar().ssm_n_offvtx_tracks_5;
  ssm_n_offvtx_tracks_8 = bdt.Getstkdar().ssm_n_offvtx_tracks_8;
  ssm_n_offvtx_tracks_11 = bdt.Getstkdar().ssm_n_offvtx_tracks_11;
  ssm_n_offvtx_showers_1 = bdt.Getstkdar().ssm_n_offvtx_showers_1;
  ssm_n_offvtx_showers_3 = bdt.Getstkdar().ssm_n_offvtx_showers_3;
  ssm_n_offvtx_showers_5 = bdt.Getstkdar().ssm_n_offvtx_showers_5;
  ssm_n_offvtx_showers_8 = bdt.Getstkdar().ssm_n_offvtx_showers_8;
  ssm_n_offvtx_showers_11 = bdt.Getstkdar().ssm_n_offvtx_showers_11;
  //properties of leading off vertex track
  ssm_offvtx_track1_pdg = bdt.Getstkdar().ssm_offvtx_track1_pdg;
  ssm_offvtx_track1_score_mu_fwd = bdt.Getstkdar().ssm_offvtx_track1_score_mu_fwd;
  ssm_offvtx_track1_score_p_fwd = bdt.Getstkdar().ssm_offvtx_track1_score_p_fwd;
  ssm_offvtx_track1_score_e_fwd = bdt.Getstkdar().ssm_offvtx_track1_score_e_fwd;
  ssm_offvtx_track1_score_mu_bck = bdt.Getstkdar().ssm_offvtx_track1_score_mu_bck;
  ssm_offvtx_track1_score_p_bck = bdt.Getstkdar().ssm_offvtx_track1_score_p_bck;
  ssm_offvtx_track1_score_e_bck = bdt.Getstkdar().ssm_offvtx_track1_score_e_bck;
  ssm_offvtx_track1_length = bdt.Getstkdar().ssm_offvtx_track1_length;
  ssm_offvtx_track1_direct_length = bdt.Getstkdar().ssm_offvtx_track1_direct_length;
  ssm_offvtx_track1_max_dev = bdt.Getstkdar().ssm_offvtx_track1_max_dev;
  ssm_offvtx_track1_kine_energy_range = bdt.Getstkdar().ssm_offvtx_track1_kine_energy_range;
  ssm_offvtx_track1_kine_energy_range_mu = bdt.Getstkdar().ssm_offvtx_track1_kine_energy_range_mu;
  ssm_offvtx_track1_kine_energy_range_p = bdt.Getstkdar().ssm_offvtx_track1_kine_energy_range_p;
  ssm_offvtx_track1_kine_energy_range_e = bdt.Getstkdar().ssm_offvtx_track1_kine_energy_range_e;
  ssm_offvtx_track1_kine_energy_cal = bdt.Getstkdar().ssm_offvtx_track1_kine_energy_cal;
  ssm_offvtx_track1_medium_dq_dx = bdt.Getstkdar().ssm_offvtx_track1_medium_dq_dx;
  ssm_offvtx_track1_x_dir = bdt.Getstkdar().ssm_offvtx_track1_x_dir;
  ssm_offvtx_track1_y_dir = bdt.Getstkdar().ssm_offvtx_track1_y_dir;
  ssm_offvtx_track1_z_dir = bdt.Getstkdar().ssm_offvtx_track1_z_dir;
  ssm_offvtx_track1_dist_mainvtx = bdt.Getstkdar().ssm_offvtx_track1_dist_mainvtx;
  //properties of leading off vertex shower
  ssm_offvtx_shw1_pdg_offvtx = bdt.Getstkdar().ssm_offvtx_shw1_pdg_offvtx;
  ssm_offvtx_shw1_score_mu_fwd = bdt.Getstkdar().ssm_offvtx_shw1_score_mu_fwd;
  ssm_offvtx_shw1_score_p_fwd = bdt.Getstkdar().ssm_offvtx_shw1_score_p_fwd;
  ssm_offvtx_shw1_score_e_fwd = bdt.Getstkdar().ssm_offvtx_shw1_score_e_fwd;
  ssm_offvtx_shw1_score_mu_bck = bdt.Getstkdar().ssm_offvtx_shw1_score_mu_bck;
  ssm_offvtx_shw1_score_p_bck = bdt.Getstkdar().ssm_offvtx_shw1_score_p_bck;
  ssm_offvtx_shw1_score_e_bck = bdt.Getstkdar().ssm_offvtx_shw1_score_e_bck;
  ssm_offvtx_shw1_length = bdt.Getstkdar().ssm_offvtx_shw1_length;
  ssm_offvtx_shw1_direct_length = bdt.Getstkdar().ssm_offvtx_shw1_direct_length;
  ssm_offvtx_shw1_max_dev = bdt.Getstkdar().ssm_offvtx_shw1_max_dev;
  ssm_offvtx_shw1_kine_energy_best = bdt.Getstkdar().ssm_offvtx_shw1_kine_energy_best;
  ssm_offvtx_shw1_kine_energy_range = bdt.Getstkdar().ssm_offvtx_shw1_kine_energy_range;
  ssm_offvtx_shw1_kine_energy_range_mu = bdt.Getstkdar().ssm_offvtx_shw1_kine_energy_range_mu;
  ssm_offvtx_shw1_kine_energy_range_p = bdt.Getstkdar().ssm_offvtx_shw1_kine_energy_range_p;
  ssm_offvtx_shw1_kine_energy_range_e = bdt.Getstkdar().ssm_offvtx_shw1_kine_energy_range_e;
  ssm_offvtx_shw1_kine_energy_cal = bdt.Getstkdar().ssm_offvtx_shw1_kine_energy_cal;
  ssm_offvtx_shw1_medium_dq_dx = bdt.Getstkdar().ssm_offvtx_shw1_medium_dq_dx;
  ssm_offvtx_shw1_x_dir = bdt.Getstkdar().ssm_offvtx_shw1_x_dir;
  ssm_offvtx_shw1_y_dir = bdt.Getstkdar().ssm_offvtx_shw1_y_dir;
  ssm_offvtx_shw1_z_dir = bdt.Getstkdar().ssm_offvtx_shw1_z_dir;
  ssm_offvtx_shw1_dist_mainvtx = bdt.Getstkdar().ssm_offvtx_shw1_dist_mainvtx;
  // Spacepoints
  ssmsp_Ntrack = bdt.Getstkdar().ssmsp_Ntrack;
  ssmsp_Nsp = bdt.Getstkdar().ssmsp_Nsp;
  ssmsp_Nsp_tot = bdt.Getstkdar().ssmsp_Nsp_tot;
  ssmsp_pdg = bdt.Getstkdar().ssmsp_pdg;
  ssmsp_id = bdt.Getstkdar().ssmsp_id;
  ssmsp_mother = bdt.Getstkdar().ssmsp_mother;
  ssmsp_x = bdt.Getstkdar().ssmsp_x;
  ssmsp_y = bdt.Getstkdar().ssmsp_y;
  ssmsp_z = bdt.Getstkdar().ssmsp_z;
  ssmsp_dx = bdt.Getstkdar().ssmsp_dx;
  ssmsp_dQ = bdt.Getstkdar().ssmsp_dQ;
  ssmsp_KE = bdt.Getstkdar().ssmsp_KE;
  ssmsp_containing_shower_id = bdt.Getstkdar().ssmsp_containing_shower_id;
  ssmsp_containing_shower_ke = bdt.Getstkdar().ssmsp_containing_shower_ke;
  ssmsp_containing_shower_flag = bdt.Getstkdar().ssmsp_containing_shower_flag;
  // Kine Vars
  ssm_kine_reco_Enu = bdt.Getstkdar().ssm_kine_reco_Enu;
  ssm_kine_reco_add_energy = bdt.Getstkdar().ssm_kine_reco_add_energy;
  ssm_kine_energy_particle = bdt.Getstkdar().ssm_kine_energy_particle;
  ssm_kine_energy_info = bdt.Getstkdar().ssm_kine_energy_info;
  ssm_kine_particle_type = bdt.Getstkdar().ssm_kine_particle_type;
  ssm_kine_energy_included = bdt.Getstkdar().ssm_kine_energy_included;
  ssm_kine_pio_mass = bdt.Getstkdar().ssm_kine_pio_mass;
  ssm_kine_pio_flag = bdt.Getstkdar().ssm_kine_pio_flag;
  ssm_kine_pio_vtx_dis = bdt.Getstkdar().ssm_kine_pio_vtx_dis;
  ssm_kine_pio_energy_1 = bdt.Getstkdar().ssm_kine_pio_energy_1;
  ssm_kine_pio_theta_1 = bdt.Getstkdar().ssm_kine_pio_theta_1;
  ssm_kine_pio_phi_1 = bdt.Getstkdar().ssm_kine_pio_phi_1;
  ssm_kine_pio_dis_1 = bdt.Getstkdar().ssm_kine_pio_dis_1;
  ssm_kine_pio_energy_2 = bdt.Getstkdar().ssm_kine_pio_energy_2;
  ssm_kine_pio_theta_2 = bdt.Getstkdar().ssm_kine_pio_theta_2;
  ssm_kine_pio_phi_2 = bdt.Getstkdar().ssm_kine_pio_phi_2;
  ssm_kine_pio_dis_2 = bdt.Getstkdar().ssm_kine_pio_dis_2;
  ssm_kine_pio_angle = bdt.Getstkdar().ssm_kine_pio_angle;
  ssm_numu_cc_flag = bdt.Getstkdar().ssm_numu_cc_flag;
  ssm_cosmict_flag_1 = bdt.Getstkdar().ssm_cosmict_flag_1; // fiducial volume vertex
  ssm_cosmict_flag_2 = bdt.Getstkdar().ssm_cosmict_flag_2;  // single muon
  ssm_cosmict_flag_3 = bdt.Getstkdar().ssm_cosmict_flag_3;  // single muon (long)
  ssm_cosmict_flag_4 = bdt.Getstkdar().ssm_cosmict_flag_4;  // kinematics muon
  ssm_cosmict_flag_5 = bdt.Getstkdar().ssm_cosmict_flag_5; // kinematics muon (long)
  ssm_cosmict_flag_6 = bdt.Getstkdar().ssm_cosmict_flag_6; // special ...
  ssm_cosmict_flag_7 = bdt.Getstkdar().ssm_cosmict_flag_7;  // muon+ michel
  ssm_cosmict_flag_8 = bdt.Getstkdar().ssm_cosmict_flag_8;  // muon + michel + special
  ssm_cosmict_flag_9 = bdt.Getstkdar().ssm_cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
  ssm_cosmict_flag_10 = bdt.Getstkdar().ssm_cosmict_flag_10;  // front upstream (dirt)
  ssm_cosmict_flag = bdt.Getstkdar().ssm_cosmict_flag;
}

void WireCellAnaTree::ReadBDTvar(nsm::NuSelectionBDT const& bdt)
{
  shw_sp_num_mip_tracks = bdt.GetSPID().shw_sp_num_mip_tracks;
  shw_sp_num_muons = bdt.GetSPID().shw_sp_num_muons;
  shw_sp_num_pions = bdt.GetSPID().shw_sp_num_pions;
  shw_sp_num_protons = bdt.GetSPID().shw_sp_num_protons;
  shw_sp_proton_length_1 = bdt.GetSPID().shw_sp_proton_length_1;
  shw_sp_proton_dqdx_1 = bdt.GetSPID().shw_sp_proton_dqdx_1;
  shw_sp_proton_energy_1 = bdt.GetSPID().shw_sp_proton_energy_1;
  shw_sp_proton_length_2 = bdt.GetSPID().shw_sp_proton_length_2;
  shw_sp_proton_dqdx_2 = bdt.GetSPID().shw_sp_proton_dqdx_2;
  shw_sp_proton_energy_2 = bdt.GetSPID().shw_sp_proton_energy_2;
  shw_sp_n_good_showers = bdt.GetSPID().shw_sp_n_good_showers;
  shw_sp_n_20mev_showers = bdt.GetSPID().shw_sp_n_20mev_showers;
  shw_sp_n_br1_showers = bdt.GetSPID().shw_sp_n_br1_showers;
  shw_sp_n_br2_showers = bdt.GetSPID().shw_sp_n_br2_showers;
  shw_sp_n_br3_showers = bdt.GetSPID().shw_sp_n_br3_showers;
  shw_sp_n_br4_showers = bdt.GetSPID().shw_sp_n_br4_showers;
  shw_sp_n_20br1_showers = bdt.GetSPID().shw_sp_n_20br1_showers;
  shw_sp_20mev_showers = bdt.GetSPID().shw_sp_20mev_showers;
  shw_sp_br1_showers = bdt.GetSPID().shw_sp_br1_showers;
  shw_sp_br2_showers = bdt.GetSPID().shw_sp_br2_showers;
  shw_sp_br3_showers = bdt.GetSPID().shw_sp_br3_showers;
  shw_sp_br4_showers = bdt.GetSPID().shw_sp_br4_showers;
  shw_sp_shw_vtx_dis = bdt.GetSPID().shw_sp_shw_vtx_dis;
  shw_sp_max_shw_dis = bdt.GetSPID().shw_sp_max_shw_dis;
  shw_sp_filled = bdt.GetSPSHWID1().shw_sp_filled;
        shw_sp_flag = bdt.GetSPSHWID1().shw_sp_flag;
        shw_sp_energy = bdt.GetSPSHWID1().shw_sp_energy;
        shw_sp_vec_dQ_dx_0 = bdt.GetSPSHWID1().shw_sp_vec_dQ_dx_0;
        shw_sp_vec_dQ_dx_1 = bdt.GetSPSHWID1().shw_sp_vec_dQ_dx_1;
        shw_sp_max_dQ_dx_sample = bdt.GetSPSHWID1().shw_sp_max_dQ_dx_sample;
        shw_sp_n_below_threshold = bdt.GetSPSHWID1().shw_sp_n_below_threshold;
        shw_sp_n_below_zero = bdt.GetSPSHWID1().shw_sp_n_below_zero;
        shw_sp_n_lowest = bdt.GetSPSHWID1().shw_sp_n_lowest;
        shw_sp_n_highest = bdt.GetSPSHWID1().shw_sp_n_highest;
        shw_sp_lowest_dQ_dx = bdt.GetSPSHWID1().shw_sp_lowest_dQ_dx;
        shw_sp_highest_dQ_dx = bdt.GetSPSHWID1().shw_sp_highest_dQ_dx;
        shw_sp_medium_dQ_dx = bdt.GetSPSHWID1().shw_sp_medium_dQ_dx;
        shw_sp_stem_length = bdt.GetSPSHWID1().shw_sp_stem_length;
        shw_sp_length_main = bdt.GetSPSHWID1().shw_sp_length_main;
        shw_sp_length_total = bdt.GetSPSHWID1().shw_sp_length_total;
        shw_sp_angle_beam = bdt.GetSPSHWID1().shw_sp_angle_beam;
        shw_sp_iso_angle = bdt.GetSPSHWID1().shw_sp_iso_angle;
        shw_sp_n_vertex = bdt.GetSPSHWID1().shw_sp_n_vertex;
        shw_sp_n_good_tracks = bdt.GetSPSHWID1().shw_sp_n_good_tracks;
        shw_sp_E_indirect_max_energy = bdt.GetSPSHWID1().shw_sp_E_indirect_max_energy;
        shw_sp_flag_all_above = bdt.GetSPSHWID1().shw_sp_flag_all_above;
        shw_sp_min_dQ_dx_5 = bdt.GetSPSHWID1().shw_sp_min_dQ_dx_5;
        shw_sp_n_other_vertex = bdt.GetSPSHWID1().shw_sp_n_other_vertex;
        shw_sp_n_stem_size = bdt.GetSPSHWID1().shw_sp_n_stem_size;
        shw_sp_flag_stem_trajectory = bdt.GetSPSHWID1().shw_sp_flag_stem_trajectory;
        shw_sp_min_dis = bdt.GetSPSHWID1().shw_sp_min_dis;
  shw_sp_vec_mean_dedx = bdt.GetSPSHWID2().shw_sp_vec_mean_dedx;
        shw_sp_vec_median_dedx = bdt.GetSPSHWID2().shw_sp_vec_median_dedx;
        shw_sp_vec_dQ_dx_2 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_2;
        shw_sp_vec_dQ_dx_3 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_3;
        shw_sp_vec_dQ_dx_4 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_4;
        shw_sp_vec_dQ_dx_5 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_5;
        shw_sp_vec_dQ_dx_6 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_6;
        shw_sp_vec_dQ_dx_7 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_7;
        shw_sp_vec_dQ_dx_8 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_8;
        shw_sp_vec_dQ_dx_9 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_9;
        shw_sp_vec_dQ_dx_10 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_10;
        shw_sp_vec_dQ_dx_11 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_11;
        shw_sp_vec_dQ_dx_12 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_12;
        shw_sp_vec_dQ_dx_13 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_13;
        shw_sp_vec_dQ_dx_14 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_14;
        shw_sp_vec_dQ_dx_15 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_15;
        shw_sp_vec_dQ_dx_16 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_16;
        shw_sp_vec_dQ_dx_17 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_17;
        shw_sp_vec_dQ_dx_18 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_18;
        shw_sp_vec_dQ_dx_19 = bdt.GetSPSHWID2().shw_sp_vec_dQ_dx_19;
        shw_sp_pio_filled = bdt.GetSPPi0Tagger1().shw_sp_pio_filled;
        shw_sp_pio_flag = bdt.GetSPPi0Tagger1().shw_sp_pio_flag;
        shw_sp_pio_mip_id = bdt.GetSPPi0Tagger1().shw_sp_pio_mip_id;
        shw_sp_pio_flag_pio = bdt.GetSPPi0Tagger1().shw_sp_pio_flag_pio;
        shw_sp_pio_1_flag = bdt.GetSPPi0Tagger1().shw_sp_pio_1_flag;
        shw_sp_pio_1_mass = bdt.GetSPPi0Tagger1().shw_sp_pio_1_mass;
        shw_sp_pio_1_pio_type = bdt.GetSPPi0Tagger1().shw_sp_pio_1_pio_type;
        shw_sp_pio_1_energy_1 = bdt.GetSPPi0Tagger1().shw_sp_pio_1_energy_1;
        shw_sp_pio_1_energy_2 = bdt.GetSPPi0Tagger1().shw_sp_pio_1_energy_2;
        shw_sp_pio_1_dis_1 = bdt.GetSPPi0Tagger1().shw_sp_pio_1_dis_1;
        shw_sp_pio_1_dis_2 = bdt.GetSPPi0Tagger1().shw_sp_pio_1_dis_2;
        shw_sp_pio_2_v_flag = bdt.GetSPPi0Tagger1().shw_sp_pio_2_v_flag;
        shw_sp_pio_2_v_dis2 = bdt.GetSPPi0Tagger1().shw_sp_pio_2_v_dis2;
        shw_sp_pio_2_v_angle2 = bdt.GetSPPi0Tagger1().shw_sp_pio_2_v_angle2;
        shw_sp_pio_2_v_acc_length = bdt.GetSPPi0Tagger1().shw_sp_pio_2_v_acc_length;
        shw_sp_lem_flag = bdt.GetSPLowEMichel().shw_sp_lem_flag;
        shw_sp_lem_shower_total_length = bdt.GetSPLowEMichel().shw_sp_lem_shower_total_length;
        shw_sp_lem_shower_main_length = bdt.GetSPLowEMichel().shw_sp_lem_shower_main_length;
        shw_sp_lem_n_3seg = bdt.GetSPLowEMichel().shw_sp_lem_n_3seg;
        shw_sp_lem_e_charge = bdt.GetSPLowEMichel().shw_sp_lem_e_charge;
        shw_sp_lem_e_dQdx = bdt.GetSPLowEMichel().shw_sp_lem_e_dQdx;
        shw_sp_lem_shower_num_segs = bdt.GetSPLowEMichel().shw_sp_lem_shower_num_segs;
        shw_sp_lem_shower_num_main_segs = bdt.GetSPLowEMichel().shw_sp_lem_shower_num_main_segs;
        shw_sp_br_filled = bdt.GetSPBadReco1().shw_sp_br_filled;
        shw_sp_br1_flag = bdt.GetSPBadReco1().shw_sp_br1_flag;
        shw_sp_br1_1_flag = bdt.GetSPBadReco1().shw_sp_br1_1_flag;
        shw_sp_br1_1_shower_type = bdt.GetSPBadReco1().shw_sp_br1_1_shower_type;
        shw_sp_br1_1_vtx_n_segs = bdt.GetSPBadReco1().shw_sp_br1_1_vtx_n_segs;
        shw_sp_br1_1_energy = bdt.GetSPBadReco1().shw_sp_br1_1_energy;
        shw_sp_br1_1_n_segs = bdt.GetSPBadReco1().shw_sp_br1_1_n_segs;
        shw_sp_br1_1_flag_sg_topology = bdt.GetSPBadReco1().shw_sp_br1_1_flag_sg_topology;
        shw_sp_br1_1_flag_sg_trajectory = bdt.GetSPBadReco1().shw_sp_br1_1_flag_sg_trajectory;
        shw_sp_br1_1_sg_length = bdt.GetSPBadReco1().shw_sp_br1_1_sg_length;
        shw_sp_br1_2_flag = bdt.GetSPBadReco1().shw_sp_br1_2_flag;
        shw_sp_br1_2_energy = bdt.GetSPBadReco1().shw_sp_br1_2_energy;
        shw_sp_br1_2_n_connected = bdt.GetSPBadReco1().shw_sp_br1_2_n_connected;
        shw_sp_br1_2_max_length = bdt.GetSPBadReco1().shw_sp_br1_2_max_length;
        shw_sp_br1_2_n_connected_1 = bdt.GetSPBadReco1().shw_sp_br1_2_n_connected_1;
        shw_sp_br1_2_vtx_n_segs = bdt.GetSPBadReco1().shw_sp_br1_2_vtx_n_segs;
        shw_sp_br1_2_n_shower_segs = bdt.GetSPBadReco1().shw_sp_br1_2_n_shower_segs;
        shw_sp_br1_2_max_length_ratio = bdt.GetSPBadReco1().shw_sp_br1_2_max_length_ratio;
        shw_sp_br1_2_shower_length = bdt.GetSPBadReco1().shw_sp_br1_2_shower_length;
        shw_sp_br1_3_flag = bdt.GetSPBadReco1().shw_sp_br1_3_flag;
        shw_sp_br1_3_energy = bdt.GetSPBadReco1().shw_sp_br1_3_energy;
        shw_sp_br1_3_n_connected_p = bdt.GetSPBadReco1().shw_sp_br1_3_n_connected_p;
        shw_sp_br1_3_max_length_p = bdt.GetSPBadReco1().shw_sp_br1_3_max_length_p;
        shw_sp_br1_3_n_shower_segs = bdt.GetSPBadReco1().shw_sp_br1_3_n_shower_segs;
        shw_sp_br1_3_flag_sg_topology = bdt.GetSPBadReco1().shw_sp_br1_3_flag_sg_topology;
        shw_sp_br1_3_flag_sg_trajectory = bdt.GetSPBadReco1().shw_sp_br1_3_flag_sg_trajectory;
        shw_sp_br1_3_n_shower_main_segs = bdt.GetSPBadReco1().shw_sp_br1_3_n_shower_main_segs;
        shw_sp_br1_3_sg_length = bdt.GetSPBadReco1().shw_sp_br1_3_sg_length;
        shw_sp_br_filled = bdt.GetSPBadReco2().shw_sp_br_filled;
        shw_sp_br2_flag = bdt.GetSPBadReco2().shw_sp_br2_flag;
        shw_sp_br2_flag_single_shower = bdt.GetSPBadReco2().shw_sp_br2_flag_single_shower;
        shw_sp_br2_num_valid_tracks = bdt.GetSPBadReco2().shw_sp_br2_num_valid_tracks;
        shw_sp_br2_energy = bdt.GetSPBadReco2().shw_sp_br2_energy;
        shw_sp_br2_angle1 = bdt.GetSPBadReco2().shw_sp_br2_angle1;
        shw_sp_br2_angle2 = bdt.GetSPBadReco2().shw_sp_br2_angle2;
        shw_sp_br2_angle = bdt.GetSPBadReco2().shw_sp_br2_angle;
        shw_sp_br2_angle3 = bdt.GetSPBadReco2().shw_sp_br2_angle3;
        shw_sp_br2_n_shower_main_segs = bdt.GetSPBadReco2().shw_sp_br2_n_shower_main_segs;
        shw_sp_br2_max_angle = bdt.GetSPBadReco2().shw_sp_br2_max_angle;
        shw_sp_br2_sg_length = bdt.GetSPBadReco2().shw_sp_br2_sg_length;
        shw_sp_br2_flag_sg_trajectory = bdt.GetSPBadReco2().shw_sp_br2_flag_sg_trajectory;
        shw_sp_br_filled = bdt.GetSPBadReco3().shw_sp_br_filled;
        shw_sp_br3_flag = bdt.GetSPBadReco3().shw_sp_br3_flag;
        shw_sp_br3_1_flag = bdt.GetSPBadReco3().shw_sp_br3_1_flag;
        shw_sp_br3_1_energy = bdt.GetSPBadReco3().shw_sp_br3_1_energy;
        shw_sp_br3_1_n_shower_segments = bdt.GetSPBadReco3().shw_sp_br3_1_n_shower_segments;
        shw_sp_br3_1_sg_flag_trajectory = bdt.GetSPBadReco3().shw_sp_br3_1_sg_flag_trajectory;
        shw_sp_br3_1_sg_direct_length = bdt.GetSPBadReco3().shw_sp_br3_1_sg_direct_length;
        shw_sp_br3_1_sg_length = bdt.GetSPBadReco3().shw_sp_br3_1_sg_length;
        shw_sp_br3_1_total_main_length = bdt.GetSPBadReco3().shw_sp_br3_1_total_main_length;
        shw_sp_br3_1_total_length = bdt.GetSPBadReco3().shw_sp_br3_1_total_length;
        shw_sp_br3_1_iso_angle = bdt.GetSPBadReco3().shw_sp_br3_1_iso_angle;
        shw_sp_br3_1_sg_flag_topology = bdt.GetSPBadReco3().shw_sp_br3_1_sg_flag_topology;
        shw_sp_br3_2_flag = bdt.GetSPBadReco3().shw_sp_br3_2_flag;
        shw_sp_br3_2_n_ele = bdt.GetSPBadReco3().shw_sp_br3_2_n_ele;
        shw_sp_br3_2_n_other = bdt.GetSPBadReco3().shw_sp_br3_2_n_other;
        shw_sp_br3_2_energy = bdt.GetSPBadReco3().shw_sp_br3_2_energy;
        shw_sp_br3_2_total_main_length = bdt.GetSPBadReco3().shw_sp_br3_2_total_main_length;
        shw_sp_br3_2_total_length = bdt.GetSPBadReco3().shw_sp_br3_2_total_length;
        shw_sp_br3_2_other_fid = bdt.GetSPBadReco3().shw_sp_br3_2_other_fid;
        shw_sp_br3_3_v_flag = bdt.GetSPBadReco3().shw_sp_br3_3_v_flag;
        shw_sp_br3_3_v_energy = bdt.GetSPBadReco3().shw_sp_br3_3_v_energy;
        shw_sp_br3_3_v_angle = bdt.GetSPBadReco3().shw_sp_br3_3_v_angle;
        shw_sp_br3_3_v_dir_length = bdt.GetSPBadReco3().shw_sp_br3_3_v_dir_length;
        shw_sp_br3_3_v_length = bdt.GetSPBadReco3().shw_sp_br3_3_v_length;
        shw_sp_br3_4_flag = bdt.GetSPBadReco3().shw_sp_br3_4_flag;
        shw_sp_br3_4_acc_length = bdt.GetSPBadReco3().shw_sp_br3_4_acc_length;
        shw_sp_br3_4_total_length = bdt.GetSPBadReco3().shw_sp_br3_4_total_length;
        shw_sp_br3_4_energy = bdt.GetSPBadReco3().shw_sp_br3_4_energy;
        shw_sp_br3_5_v_flag = bdt.GetSPBadReco3().shw_sp_br3_5_v_flag;
        shw_sp_br3_5_v_dir_length = bdt.GetSPBadReco3().shw_sp_br3_5_v_dir_length;
        shw_sp_br3_5_v_total_length = bdt.GetSPBadReco3().shw_sp_br3_5_v_total_length;
        shw_sp_br3_5_v_flag_avoid_muon_check = bdt.GetSPBadReco3().shw_sp_br3_5_v_flag_avoid_muon_check;
        shw_sp_br3_5_v_n_seg = bdt.GetSPBadReco3().shw_sp_br3_5_v_n_seg;
        shw_sp_br3_5_v_angle = bdt.GetSPBadReco3().shw_sp_br3_5_v_angle;
        shw_sp_br3_5_v_sg_length = bdt.GetSPBadReco3().shw_sp_br3_5_v_sg_length;
        shw_sp_br3_5_v_energy = bdt.GetSPBadReco3().shw_sp_br3_5_v_energy;
        shw_sp_br3_5_v_n_main_segs = bdt.GetSPBadReco3().shw_sp_br3_5_v_n_main_segs;
        shw_sp_br3_5_v_n_segs = bdt.GetSPBadReco3().shw_sp_br3_5_v_n_segs;
        shw_sp_br3_5_v_shower_main_length = bdt.GetSPBadReco3().shw_sp_br3_5_v_shower_main_length;
        shw_sp_br3_5_v_shower_total_length = bdt.GetSPBadReco3().shw_sp_br3_5_v_shower_total_length;
        shw_sp_br3_6_v_flag = bdt.GetSPBadReco3().shw_sp_br3_6_v_flag;
        shw_sp_br3_6_v_angle = bdt.GetSPBadReco3().shw_sp_br3_6_v_angle;
        shw_sp_br3_6_v_angle1 = bdt.GetSPBadReco3().shw_sp_br3_6_v_angle1;
        shw_sp_br3_6_v_flag_shower_trajectory = bdt.GetSPBadReco3().shw_sp_br3_6_v_flag_shower_trajectory;
        shw_sp_br3_6_v_direct_length = bdt.GetSPBadReco3().shw_sp_br3_6_v_direct_length;
        shw_sp_br3_6_v_length = bdt.GetSPBadReco3().shw_sp_br3_6_v_length;
        shw_sp_br3_6_v_n_other_vtx_segs = bdt.GetSPBadReco3().shw_sp_br3_6_v_n_other_vtx_segs;
        shw_sp_br3_6_v_energy = bdt.GetSPBadReco3().shw_sp_br3_6_v_energy;
        shw_sp_br3_7_flag = bdt.GetSPBadReco3().shw_sp_br3_7_flag;
        shw_sp_br3_7_energy = bdt.GetSPBadReco3().shw_sp_br3_7_energy;
        shw_sp_br3_7_min_angle = bdt.GetSPBadReco3().shw_sp_br3_7_min_angle;
        shw_sp_br3_7_sg_length = bdt.GetSPBadReco3().shw_sp_br3_7_sg_length;
        shw_sp_br3_7_shower_main_length = bdt.GetSPBadReco3().shw_sp_br3_7_shower_main_length;
        shw_sp_br3_8_flag = bdt.GetSPBadReco3().shw_sp_br3_8_flag;
        shw_sp_br3_8_max_dQ_dx = bdt.GetSPBadReco3().shw_sp_br3_8_max_dQ_dx;
        shw_sp_br3_8_energy = bdt.GetSPBadReco3().shw_sp_br3_8_energy;
        shw_sp_br3_8_n_main_segs = bdt.GetSPBadReco3().shw_sp_br3_8_n_main_segs;
        shw_sp_br3_8_shower_main_length = bdt.GetSPBadReco3().shw_sp_br3_8_shower_main_length;
        shw_sp_br3_8_shower_length = bdt.GetSPBadReco3().shw_sp_br3_8_shower_length;
        shw_sp_br_filled = bdt.GetSPBadReco4().shw_sp_br_filled;
        shw_sp_br4_flag = bdt.GetSPBadReco4().shw_sp_br4_flag;
        shw_sp_br4_1_flag = bdt.GetSPBadReco4().shw_sp_br4_1_flag;
        shw_sp_br4_1_shower_main_length = bdt.GetSPBadReco4().shw_sp_br4_1_shower_main_length;
        shw_sp_br4_1_shower_total_length = bdt.GetSPBadReco4().shw_sp_br4_1_shower_total_length;
        shw_sp_br4_1_min_dis = bdt.GetSPBadReco4().shw_sp_br4_1_min_dis;
        shw_sp_br4_1_energy = bdt.GetSPBadReco4().shw_sp_br4_1_energy;
        shw_sp_br4_1_flag_avoid_muon_check = bdt.GetSPBadReco4().shw_sp_br4_1_flag_avoid_muon_check;
        shw_sp_br4_1_n_vtx_segs = bdt.GetSPBadReco4().shw_sp_br4_1_n_vtx_segs;
        shw_sp_br4_1_n_main_segs = bdt.GetSPBadReco4().shw_sp_br4_1_n_main_segs;
        shw_sp_br4_2_flag = bdt.GetSPBadReco4().shw_sp_br4_2_flag;
        shw_sp_br4_2_ratio_45 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio_45;
        shw_sp_br4_2_ratio_35 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio_35;
        shw_sp_br4_2_ratio_25 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio_25;
        shw_sp_br4_2_ratio_15 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio_15;
        shw_sp_br4_2_energy = bdt.GetSPBadReco4().shw_sp_br4_2_energy;
        shw_sp_br4_2_ratio1_45 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio1_45;
        shw_sp_br4_2_ratio1_35 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio1_35;
        shw_sp_br4_2_ratio1_25 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio1_25;
        shw_sp_br4_2_ratio1_15 = bdt.GetSPBadReco4().shw_sp_br4_2_ratio1_15;
        shw_sp_br4_2_iso_angle = bdt.GetSPBadReco4().shw_sp_br4_2_iso_angle;
        shw_sp_br4_2_iso_angle1 = bdt.GetSPBadReco4().shw_sp_br4_2_iso_angle1;
        shw_sp_br4_2_angle = bdt.GetSPBadReco4().shw_sp_br4_2_angle;
        shw_sp_hol_flag = bdt.GetSPHighEoverlap().shw_sp_hol_flag;
        shw_sp_hol_1_flag = bdt.GetSPHighEoverlap().shw_sp_hol_1_flag;
        shw_sp_hol_1_n_valid_tracks = bdt.GetSPHighEoverlap().shw_sp_hol_1_n_valid_tracks;
        shw_sp_hol_1_min_angle = bdt.GetSPHighEoverlap().shw_sp_hol_1_min_angle;
        shw_sp_hol_1_energy = bdt.GetSPHighEoverlap().shw_sp_hol_1_energy;
        shw_sp_hol_1_flag_all_shower = bdt.GetSPHighEoverlap().shw_sp_hol_1_flag_all_shower;
        shw_sp_hol_1_min_length = bdt.GetSPHighEoverlap().shw_sp_hol_1_min_length;
        shw_sp_hol_2_flag = bdt.GetSPHighEoverlap().shw_sp_hol_2_flag;
        shw_sp_hol_2_min_angle = bdt.GetSPHighEoverlap().shw_sp_hol_2_min_angle;
        shw_sp_hol_2_medium_dQ_dx = bdt.GetSPHighEoverlap().shw_sp_hol_2_medium_dQ_dx;
        shw_sp_hol_2_ncount = bdt.GetSPHighEoverlap().shw_sp_hol_2_ncount;
        shw_sp_hol_2_energy = bdt.GetSPHighEoverlap().shw_sp_hol_2_energy;
        shw_sp_lol_flag = bdt.GetSPLowEoverlap().shw_sp_lol_flag;
        shw_sp_lol_1_v_flag = bdt.GetSPLowEoverlap().shw_sp_lol_1_v_flag;
        shw_sp_lol_1_v_energy = bdt.GetSPLowEoverlap().shw_sp_lol_1_v_energy;
        shw_sp_lol_1_v_vtx_n_segs = bdt.GetSPLowEoverlap().shw_sp_lol_1_v_vtx_n_segs;
        shw_sp_lol_1_v_nseg = bdt.GetSPLowEoverlap().shw_sp_lol_1_v_nseg;
        shw_sp_lol_1_v_angle = bdt.GetSPLowEoverlap().shw_sp_lol_1_v_angle;
        shw_sp_lol_2_v_flag = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_flag;
        shw_sp_lol_2_v_length = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_length;
        shw_sp_lol_2_v_angle = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_angle;
        shw_sp_lol_2_v_type = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_type;
        shw_sp_lol_2_v_vtx_n_segs = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_vtx_n_segs;
        shw_sp_lol_2_v_energy = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_energy;
        shw_sp_lol_2_v_shower_main_length = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_shower_main_length;
        shw_sp_lol_2_v_flag_dir_weak = bdt.GetSPLowEoverlap().shw_sp_lol_2_v_flag_dir_weak;
        shw_sp_lol_3_flag = bdt.GetSPLowEoverlap().shw_sp_lol_3_flag;
        shw_sp_lol_3_angle_beam = bdt.GetSPLowEoverlap().shw_sp_lol_3_angle_beam;
        shw_sp_lol_3_n_valid_tracks = bdt.GetSPLowEoverlap().shw_sp_lol_3_n_valid_tracks;
        shw_sp_lol_3_min_angle = bdt.GetSPLowEoverlap().shw_sp_lol_3_min_angle;
        shw_sp_lol_3_vtx_n_segs = bdt.GetSPLowEoverlap().shw_sp_lol_3_vtx_n_segs;
        shw_sp_lol_3_energy = bdt.GetSPLowEoverlap().shw_sp_lol_3_energy;
        shw_sp_lol_3_shower_main_length = bdt.GetSPLowEoverlap().shw_sp_lol_3_shower_main_length;
        shw_sp_lol_3_n_out = bdt.GetSPLowEoverlap().shw_sp_lol_3_n_out;
        shw_sp_lol_3_n_sum = bdt.GetSPLowEoverlap().shw_sp_lol_3_n_sum;

        cosmic_filled = bdt.GetCosmicTagger().cosmic_filled;
        cosmic_flag = bdt.GetCosmicTagger().cosmic_flag;
        cosmic_n_solid_tracks = bdt.GetCosmicTagger().cosmic_n_solid_tracks;
        cosmic_energy_main_showers = bdt.GetCosmicTagger().cosmic_energy_main_showers;
        cosmic_energy_direct_showers = bdt.GetCosmicTagger().cosmic_energy_direct_showers;
        cosmic_energy_indirect_showers = bdt.GetCosmicTagger().cosmic_energy_indirect_showers;
        cosmic_n_direct_showers = bdt.GetCosmicTagger().cosmic_n_direct_showers;
        cosmic_n_indirect_showers = bdt.GetCosmicTagger().cosmic_n_indirect_showers;
        cosmic_n_main_showers = bdt.GetCosmicTagger().cosmic_n_main_showers;
        gap_filled = bdt.GetGapID().gap_filled;
        gap_flag = bdt.GetGapID().gap_flag;
        gap_flag_prolong_u = bdt.GetGapID().gap_flag_prolong_u;
        gap_flag_prolong_v = bdt.GetGapID().gap_flag_prolong_v;
        gap_flag_prolong_w = bdt.GetGapID().gap_flag_prolong_w;
        gap_flag_parallel = bdt.GetGapID().gap_flag_parallel;
        gap_n_points = bdt.GetGapID().gap_n_points;
        gap_n_bad = bdt.GetGapID().gap_n_bad;
        gap_energy = bdt.GetGapID().gap_energy;
        gap_num_valid_tracks = bdt.GetGapID().gap_num_valid_tracks;
        gap_flag_single_shower = bdt.GetGapID().gap_flag_single_shower;
        mip_quality_filled = bdt.GetMipCheck().mip_quality_filled;
        mip_quality_flag = bdt.GetMipCheck().mip_quality_flag;
        mip_quality_energy = bdt.GetMipCheck().mip_quality_energy;
        mip_quality_overlap = bdt.GetMipCheck().mip_quality_overlap;
        mip_quality_n_showers = bdt.GetMipCheck().mip_quality_n_showers;
        mip_quality_n_tracks = bdt.GetMipCheck().mip_quality_n_tracks;
        mip_quality_flag_inside_pi0 = bdt.GetMipCheck().mip_quality_flag_inside_pi0;
        mip_quality_n_pi0_showers = bdt.GetMipCheck().mip_quality_n_pi0_showers;
        mip_quality_shortest_length = bdt.GetMipCheck().mip_quality_shortest_length;
        mip_quality_acc_length = bdt.GetMipCheck().mip_quality_acc_length;
        mip_quality_shortest_angle = bdt.GetMipCheck().mip_quality_shortest_angle;
        mip_quality_flag_proton = bdt.GetMipCheck().mip_quality_flag_proton;
        mip_filled = bdt.GetMipID1().mip_filled;
        mip_flag = bdt.GetMipID1().mip_flag;
        mip_energy = bdt.GetMipID1().mip_energy;
        mip_n_end_reduction = bdt.GetMipID1().mip_n_end_reduction;
        mip_n_first_mip = bdt.GetMipID1().mip_n_first_mip;
        mip_n_first_non_mip = bdt.GetMipID1().mip_n_first_non_mip;
        mip_n_first_non_mip_1 = bdt.GetMipID1().mip_n_first_non_mip_1;
        mip_n_first_non_mip_2 = bdt.GetMipID1().mip_n_first_non_mip_2;
        mip_vec_dQ_dx_0 = bdt.GetMipID1().mip_vec_dQ_dx_0;
        mip_vec_dQ_dx_1 = bdt.GetMipID1().mip_vec_dQ_dx_1;
        mip_max_dQ_dx_sample = bdt.GetMipID1().mip_max_dQ_dx_sample;
        mip_n_below_threshold = bdt.GetMipID1().mip_n_below_threshold;
        mip_n_below_zero = bdt.GetMipID1().mip_n_below_zero;
        mip_n_lowest = bdt.GetMipID1().mip_n_lowest;
        mip_n_highest = bdt.GetMipID1().mip_n_highest;
        mip_lowest_dQ_dx = bdt.GetMipID1().mip_lowest_dQ_dx;
        mip_highest_dQ_dx = bdt.GetMipID1().mip_highest_dQ_dx;
        mip_medium_dQ_dx = bdt.GetMipID1().mip_medium_dQ_dx;
        mip_stem_length = bdt.GetMipID1().mip_stem_length;
        mip_length_main = bdt.GetMipID1().mip_length_main;
        mip_length_total = bdt.GetMipID1().mip_length_total;
        mip_angle_beam = bdt.GetMipID1().mip_angle_beam;
        mip_iso_angle = bdt.GetMipID1().mip_iso_angle;
        mip_n_vertex = bdt.GetMipID1().mip_n_vertex;
        mip_n_good_tracks = bdt.GetMipID1().mip_n_good_tracks;
        mip_E_indirect_max_energy = bdt.GetMipID1().mip_E_indirect_max_energy;
        mip_flag_all_above = bdt.GetMipID1().mip_flag_all_above;
        mip_min_dQ_dx_5 = bdt.GetMipID1().mip_min_dQ_dx_5;
        mip_n_other_vertex = bdt.GetMipID1().mip_n_other_vertex;
        mip_n_stem_size = bdt.GetMipID1().mip_n_stem_size;
        mip_flag_stem_trajectory = bdt.GetMipID1().mip_flag_stem_trajectory;
        mip_min_dis = bdt.GetMipID1().mip_min_dis;
        mip_vec_dQ_dx_2 = bdt.GetMipID2().mip_vec_dQ_dx_2;
        mip_vec_dQ_dx_3 = bdt.GetMipID2().mip_vec_dQ_dx_3;
        mip_vec_dQ_dx_4 = bdt.GetMipID2().mip_vec_dQ_dx_4;
        mip_vec_dQ_dx_5 = bdt.GetMipID2().mip_vec_dQ_dx_5;
        mip_vec_dQ_dx_6 = bdt.GetMipID2().mip_vec_dQ_dx_6;
        mip_vec_dQ_dx_7 = bdt.GetMipID2().mip_vec_dQ_dx_7;
        mip_vec_dQ_dx_8 = bdt.GetMipID2().mip_vec_dQ_dx_8;
        mip_vec_dQ_dx_9 = bdt.GetMipID2().mip_vec_dQ_dx_9;
        mip_vec_dQ_dx_10 = bdt.GetMipID2().mip_vec_dQ_dx_10;
        mip_vec_dQ_dx_11 = bdt.GetMipID2().mip_vec_dQ_dx_11;
        mip_vec_dQ_dx_12 = bdt.GetMipID2().mip_vec_dQ_dx_12;
        mip_vec_dQ_dx_13 = bdt.GetMipID2().mip_vec_dQ_dx_13;
        mip_vec_dQ_dx_14 = bdt.GetMipID2().mip_vec_dQ_dx_14;
        mip_vec_dQ_dx_15 = bdt.GetMipID2().mip_vec_dQ_dx_15;
        mip_vec_dQ_dx_16 = bdt.GetMipID2().mip_vec_dQ_dx_16;
        mip_vec_dQ_dx_17 = bdt.GetMipID2().mip_vec_dQ_dx_17;
        mip_vec_dQ_dx_18 = bdt.GetMipID2().mip_vec_dQ_dx_18;
        mip_vec_dQ_dx_19 = bdt.GetMipID2().mip_vec_dQ_dx_19;
        pio_filled = bdt.GetPi0Tagger1().pio_filled;
        pio_flag = bdt.GetPi0Tagger1().pio_flag;
        pio_mip_id = bdt.GetPi0Tagger1().pio_mip_id;
        pio_flag_pio = bdt.GetPi0Tagger1().pio_flag_pio;
        pio_1_flag = bdt.GetPi0Tagger1().pio_1_flag;
        pio_1_mass = bdt.GetPi0Tagger1().pio_1_mass;
        pio_1_pio_type = bdt.GetPi0Tagger1().pio_1_pio_type;
        pio_1_energy_1 = bdt.GetPi0Tagger1().pio_1_energy_1;
        pio_1_energy_2 = bdt.GetPi0Tagger1().pio_1_energy_2;
        pio_1_dis_1 = bdt.GetPi0Tagger1().pio_1_dis_1;
        pio_1_dis_2 = bdt.GetPi0Tagger1().pio_1_dis_2;
        pio_2_v_flag = bdt.GetPi0Tagger1().pio_2_v_flag;
        pio_2_v_dis2 = bdt.GetPi0Tagger1().pio_2_v_dis2;
        pio_2_v_angle2 = bdt.GetPi0Tagger1().pio_2_v_angle2;
        pio_2_v_acc_length = bdt.GetPi0Tagger1().pio_2_v_acc_length;
        sig_flag = bdt.GetPi0Tagger2().sig_flag;
        sig_1_v_flag = bdt.GetPi0Tagger2().sig_1_v_flag;
        sig_1_v_angle = bdt.GetPi0Tagger2().sig_1_v_angle;
        sig_1_v_flag_single_shower = bdt.GetPi0Tagger2().sig_1_v_flag_single_shower;
        sig_1_v_energy = bdt.GetPi0Tagger2().sig_1_v_energy;
        sig_1_v_energy_1 = bdt.GetPi0Tagger2().sig_1_v_energy_1;
        sig_2_v_flag = bdt.GetPi0Tagger2().sig_2_v_flag;
        sig_2_v_energy = bdt.GetPi0Tagger2().sig_2_v_energy;
        sig_2_v_shower_angle = bdt.GetPi0Tagger2().sig_2_v_shower_angle;
        sig_2_v_flag_single_shower = bdt.GetPi0Tagger2().sig_2_v_flag_single_shower;
        sig_2_v_medium_dQ_dx = bdt.GetPi0Tagger2().sig_2_v_medium_dQ_dx;
        sig_2_v_start_dQ_dx = bdt.GetPi0Tagger2().sig_2_v_start_dQ_dx;
        mgo_flag = bdt.GetMultiGamma1().mgo_flag;
        mgo_energy = bdt.GetMultiGamma1().mgo_energy;
        mgo_max_energy = bdt.GetMultiGamma1().mgo_max_energy;
        mgo_total_energy = bdt.GetMultiGamma1().mgo_total_energy;
        mgo_n_showers = bdt.GetMultiGamma1().mgo_n_showers;
        mgo_max_energy_1 = bdt.GetMultiGamma1().mgo_max_energy_1;
        mgo_max_energy_2 = bdt.GetMultiGamma1().mgo_max_energy_2;
        mgo_total_other_energy = bdt.GetMultiGamma1().mgo_total_other_energy;
        mgo_n_total_showers = bdt.GetMultiGamma1().mgo_n_total_showers;
        mgo_total_other_energy_1 = bdt.GetMultiGamma1().mgo_total_other_energy_1;
        mgt_flag = bdt.GetMultiGamma2().mgt_flag;
        mgt_flag_single_shower = bdt.GetMultiGamma2().mgt_flag_single_shower;
        mgt_max_energy = bdt.GetMultiGamma2().mgt_max_energy;
        mgt_energy = bdt.GetMultiGamma2().mgt_energy;
        mgt_total_other_energy = bdt.GetMultiGamma2().mgt_total_other_energy;
        mgt_max_energy_1 = bdt.GetMultiGamma2().mgt_max_energy_1;
        mgt_e_indirect_max_energy = bdt.GetMultiGamma2().mgt_e_indirect_max_energy;
        mgt_e_direct_max_energy = bdt.GetMultiGamma2().mgt_e_direct_max_energy;
        mgt_n_direct_showers = bdt.GetMultiGamma2().mgt_n_direct_showers;
        mgt_e_direct_total_energy = bdt.GetMultiGamma2().mgt_e_direct_total_energy;
        mgt_flag_indirect_max_pio = bdt.GetMultiGamma2().mgt_flag_indirect_max_pio;
        mgt_e_indirect_total_energy = bdt.GetMultiGamma2().mgt_e_indirect_total_energy;
        stw_flag = bdt.GetSingleGamma1().stw_flag;
        stw_1_flag = bdt.GetSingleGamma1().stw_1_flag;
        stw_1_energy = bdt.GetSingleGamma1().stw_1_energy;
        stw_1_dis = bdt.GetSingleGamma1().stw_1_dis;
        stw_1_dQ_dx = bdt.GetSingleGamma1().stw_1_dQ_dx;
        stw_1_flag_single_shower = bdt.GetSingleGamma1().stw_1_flag_single_shower;
        stw_1_n_pi0 = bdt.GetSingleGamma1().stw_1_n_pi0;
        stw_1_num_valid_tracks = bdt.GetSingleGamma1().stw_1_num_valid_tracks;
        stw_2_v_flag = bdt.GetSingleGamma1().stw_2_v_flag;
        stw_2_v_medium_dQ_dx = bdt.GetSingleGamma1().stw_2_v_medium_dQ_dx;
        stw_2_v_energy = bdt.GetSingleGamma1().stw_2_v_energy;
        stw_2_v_angle = bdt.GetSingleGamma1().stw_2_v_angle;
        stw_2_v_dir_length = bdt.GetSingleGamma1().stw_2_v_dir_length;
        stw_2_v_max_dQ_dx = bdt.GetSingleGamma1().stw_2_v_max_dQ_dx;
        stw_3_v_flag = bdt.GetSingleGamma1().stw_3_v_flag;
        stw_3_v_angle = bdt.GetSingleGamma1().stw_3_v_angle;
        stw_3_v_dir_length = bdt.GetSingleGamma1().stw_3_v_dir_length;
        stw_3_v_energy = bdt.GetSingleGamma1().stw_3_v_energy;
        stw_3_v_medium_dQ_dx = bdt.GetSingleGamma1().stw_3_v_medium_dQ_dx;
        stw_4_v_flag = bdt.GetSingleGamma1().stw_4_v_flag;
        stw_4_v_angle = bdt.GetSingleGamma1().stw_4_v_angle;
        stw_4_v_dis = bdt.GetSingleGamma1().stw_4_v_dis;
        stw_4_v_energy = bdt.GetSingleGamma1().stw_4_v_energy;
        spt_flag = bdt.GetSingleGamma2().spt_flag;
        spt_flag_single_shower = bdt.GetSingleGamma2().spt_flag_single_shower;
        spt_energy = bdt.GetSingleGamma2().spt_energy;
        spt_shower_main_length = bdt.GetSingleGamma2().spt_shower_main_length;
        spt_shower_total_length = bdt.GetSingleGamma2().spt_shower_total_length;
        spt_angle_beam = bdt.GetSingleGamma2().spt_angle_beam;
        spt_angle_vertical = bdt.GetSingleGamma2().spt_angle_vertical;
        spt_max_dQ_dx = bdt.GetSingleGamma2().spt_max_dQ_dx;
        spt_angle_beam_1 = bdt.GetSingleGamma2().spt_angle_beam_1;
        spt_angle_drift = bdt.GetSingleGamma2().spt_angle_drift;
        spt_angle_drift_1 = bdt.GetSingleGamma2().spt_angle_drift_1;
        spt_num_valid_tracks = bdt.GetSingleGamma2().spt_num_valid_tracks;
        spt_n_vtx_segs = bdt.GetSingleGamma2().spt_n_vtx_segs;
        spt_max_length = bdt.GetSingleGamma2().spt_max_length;
        stem_len_flag = bdt.GetStemLen().stem_len_flag;
        stem_len_energy = bdt.GetStemLen().stem_len_energy;
        stem_len_length = bdt.GetStemLen().stem_len_length;
        stem_len_flag_avoid_muon_check = bdt.GetStemLen().stem_len_flag_avoid_muon_check;
        stem_len_num_daughters = bdt.GetStemLen().stem_len_num_daughters;
        stem_len_daughter_length = bdt.GetStemLen().stem_len_daughter_length;
        lem_flag = bdt.GetLowEMichel().lem_flag;
        lem_shower_total_length = bdt.GetLowEMichel().lem_shower_total_length;
        lem_shower_main_length = bdt.GetLowEMichel().lem_shower_main_length;
        lem_n_3seg = bdt.GetLowEMichel().lem_n_3seg;
        lem_e_charge = bdt.GetLowEMichel().lem_e_charge;
        lem_e_dQdx = bdt.GetLowEMichel().lem_e_dQdx;
        lem_shower_num_segs = bdt.GetLowEMichel().lem_shower_num_segs;
        lem_shower_num_main_segs = bdt.GetLowEMichel().lem_shower_num_main_segs;
        brm_flag = bdt.GetBrokenMuon().brm_flag;
        brm_n_mu_segs = bdt.GetBrokenMuon().brm_n_mu_segs;
        brm_Ep = bdt.GetBrokenMuon().brm_Ep;
        brm_energy = bdt.GetBrokenMuon().brm_energy;
        brm_acc_length = bdt.GetBrokenMuon().brm_acc_length;
        brm_shower_total_length = bdt.GetBrokenMuon().brm_shower_total_length;
        brm_connected_length = bdt.GetBrokenMuon().brm_connected_length;
        brm_n_size = bdt.GetBrokenMuon().brm_n_size;
        brm_acc_direct_length = bdt.GetBrokenMuon().brm_acc_direct_length;
        brm_n_shower_main_segs = bdt.GetBrokenMuon().brm_n_shower_main_segs;
        brm_n_mu_main = bdt.GetBrokenMuon().brm_n_mu_main;
        cme_flag = bdt.GetMuEnergy().cme_flag;
        cme_mu_energy = bdt.GetMuEnergy().cme_mu_energy;
        cme_energy = bdt.GetMuEnergy().cme_energy;
        cme_mu_length = bdt.GetMuEnergy().cme_mu_length;
        cme_length = bdt.GetMuEnergy().cme_length;
        cme_angle_beam = bdt.GetMuEnergy().cme_angle_beam;
        anc_flag = bdt.GetShowerAngle().anc_flag;
        anc_energy = bdt.GetShowerAngle().anc_energy;
        anc_angle = bdt.GetShowerAngle().anc_angle;
        anc_max_angle = bdt.GetShowerAngle().anc_max_angle;
        anc_max_length = bdt.GetShowerAngle().anc_max_length;
        anc_acc_forward_length = bdt.GetShowerAngle().anc_acc_forward_length;
        anc_acc_backward_length = bdt.GetShowerAngle().anc_acc_backward_length;
        anc_acc_forward_length1 = bdt.GetShowerAngle().anc_acc_forward_length1;
        anc_shower_main_length = bdt.GetShowerAngle().anc_shower_main_length;
        anc_shower_total_length = bdt.GetShowerAngle().anc_shower_total_length;
        anc_flag_main_outside = bdt.GetShowerAngle().anc_flag_main_outside;
        stem_dir_filled = bdt.GetBadStem().stem_dir_filled;
        stem_dir_flag = bdt.GetBadStem().stem_dir_flag;
        stem_dir_flag_single_shower = bdt.GetBadStem().stem_dir_flag_single_shower;
        stem_dir_angle = bdt.GetBadStem().stem_dir_angle;
        stem_dir_energy = bdt.GetBadStem().stem_dir_energy;
        stem_dir_angle1 = bdt.GetBadStem().stem_dir_angle1;
        stem_dir_angle2 = bdt.GetBadStem().stem_dir_angle2;
        stem_dir_angle3 = bdt.GetBadStem().stem_dir_angle3;
        stem_dir_ratio = bdt.GetBadStem().stem_dir_ratio;
        vis_flag = bdt.GetVtxInShw().vis_flag;
        vis_1_filled = bdt.GetVtxInShw().vis_1_filled;
        vis_1_flag = bdt.GetVtxInShw().vis_1_flag;
        vis_1_n_vtx_segs = bdt.GetVtxInShw().vis_1_n_vtx_segs;
        vis_1_energy = bdt.GetVtxInShw().vis_1_energy;
        vis_1_num_good_tracks = bdt.GetVtxInShw().vis_1_num_good_tracks;
        vis_1_max_angle = bdt.GetVtxInShw().vis_1_max_angle;
        vis_1_max_shower_angle = bdt.GetVtxInShw().vis_1_max_shower_angle;
        vis_1_tmp_length1 = bdt.GetVtxInShw().vis_1_tmp_length1;
        vis_1_tmp_length2 = bdt.GetVtxInShw().vis_1_tmp_length2;
        vis_1_particle_type = bdt.GetVtxInShw().vis_1_particle_type;
        vis_2_filled = bdt.GetVtxInShw().vis_2_filled;
        vis_2_flag = bdt.GetVtxInShw().vis_2_flag;
        vis_2_n_vtx_segs = bdt.GetVtxInShw().vis_2_n_vtx_segs;
        vis_2_min_angle = bdt.GetVtxInShw().vis_2_min_angle;
        vis_2_min_weak_track = bdt.GetVtxInShw().vis_2_min_weak_track;
        vis_2_angle_beam = bdt.GetVtxInShw().vis_2_angle_beam;
        vis_2_min_angle1 = bdt.GetVtxInShw().vis_2_min_angle1;
        vis_2_iso_angle1 = bdt.GetVtxInShw().vis_2_iso_angle1;
        vis_2_min_medium_dQ_dx = bdt.GetVtxInShw().vis_2_min_medium_dQ_dx;
        vis_2_min_length = bdt.GetVtxInShw().vis_2_min_length;
        vis_2_sg_length = bdt.GetVtxInShw().vis_2_sg_length;
        vis_2_max_angle = bdt.GetVtxInShw().vis_2_max_angle;
        vis_2_max_weak_track = bdt.GetVtxInShw().vis_2_max_weak_track;
        br_filled = bdt.GetBadReco1().br_filled;
        br1_flag = bdt.GetBadReco1().br1_flag;
        br1_1_flag = bdt.GetBadReco1().br1_1_flag;
        br1_1_shower_type = bdt.GetBadReco1().br1_1_shower_type;
        br1_1_vtx_n_segs = bdt.GetBadReco1().br1_1_vtx_n_segs;
        br1_1_energy = bdt.GetBadReco1().br1_1_energy;
        br1_1_n_segs = bdt.GetBadReco1().br1_1_n_segs;
        br1_1_flag_sg_topology = bdt.GetBadReco1().br1_1_flag_sg_topology;
        br1_1_flag_sg_trajectory = bdt.GetBadReco1().br1_1_flag_sg_trajectory;
        br1_1_sg_length = bdt.GetBadReco1().br1_1_sg_length;
        br1_2_flag = bdt.GetBadReco1().br1_2_flag;
        br1_2_energy = bdt.GetBadReco1().br1_2_energy;
        br1_2_n_connected = bdt.GetBadReco1().br1_2_n_connected;
        br1_2_max_length = bdt.GetBadReco1().br1_2_max_length;
        br1_2_n_connected_1 = bdt.GetBadReco1().br1_2_n_connected_1;
        br1_2_vtx_n_segs = bdt.GetBadReco1().br1_2_vtx_n_segs;
        br1_2_n_shower_segs = bdt.GetBadReco1().br1_2_n_shower_segs;
        br1_2_max_length_ratio = bdt.GetBadReco1().br1_2_max_length_ratio;
        br1_2_shower_length = bdt.GetBadReco1().br1_2_shower_length;
        br1_3_flag = bdt.GetBadReco1().br1_3_flag;
        br1_3_energy = bdt.GetBadReco1().br1_3_energy;
        br1_3_n_connected_p = bdt.GetBadReco1().br1_3_n_connected_p;
        br1_3_max_length_p = bdt.GetBadReco1().br1_3_max_length_p;
        br1_3_n_shower_segs = bdt.GetBadReco1().br1_3_n_shower_segs;
        br1_3_flag_sg_topology = bdt.GetBadReco1().br1_3_flag_sg_topology;
        br1_3_flag_sg_trajectory = bdt.GetBadReco1().br1_3_flag_sg_trajectory;
        br1_3_n_shower_main_segs = bdt.GetBadReco1().br1_3_n_shower_main_segs;
        br1_3_sg_length = bdt.GetBadReco1().br1_3_sg_length;
        br_filled = bdt.GetBadReco2().br_filled;
        br2_flag = bdt.GetBadReco2().br2_flag;
        br2_flag_single_shower = bdt.GetBadReco2().br2_flag_single_shower;
        br2_num_valid_tracks = bdt.GetBadReco2().br2_num_valid_tracks;
        br2_energy = bdt.GetBadReco2().br2_energy;
        br2_angle1 = bdt.GetBadReco2().br2_angle1;
        br2_angle2 = bdt.GetBadReco2().br2_angle2;
        br2_angle = bdt.GetBadReco2().br2_angle;
        br2_angle3 = bdt.GetBadReco2().br2_angle3;
        br2_n_shower_main_segs = bdt.GetBadReco2().br2_n_shower_main_segs;
        br2_max_angle = bdt.GetBadReco2().br2_max_angle;
        br2_sg_length = bdt.GetBadReco2().br2_sg_length;
        br2_flag_sg_trajectory = bdt.GetBadReco2().br2_flag_sg_trajectory;
        br_filled = bdt.GetBadReco3().br_filled;
        br3_flag = bdt.GetBadReco3().br3_flag;
        br3_1_flag = bdt.GetBadReco3().br3_1_flag;
        br3_1_energy = bdt.GetBadReco3().br3_1_energy;
        br3_1_n_shower_segments = bdt.GetBadReco3().br3_1_n_shower_segments;
        br3_1_sg_flag_trajectory = bdt.GetBadReco3().br3_1_sg_flag_trajectory;
        br3_1_sg_direct_length = bdt.GetBadReco3().br3_1_sg_direct_length;
        br3_1_sg_length = bdt.GetBadReco3().br3_1_sg_length;
        br3_1_total_main_length = bdt.GetBadReco3().br3_1_total_main_length;
        br3_1_total_length = bdt.GetBadReco3().br3_1_total_length;
        br3_1_iso_angle = bdt.GetBadReco3().br3_1_iso_angle;
        br3_1_sg_flag_topology = bdt.GetBadReco3().br3_1_sg_flag_topology;
        br3_2_flag = bdt.GetBadReco3().br3_2_flag;
        br3_2_n_ele = bdt.GetBadReco3().br3_2_n_ele;
        br3_2_n_other = bdt.GetBadReco3().br3_2_n_other;
        br3_2_energy = bdt.GetBadReco3().br3_2_energy;
        br3_2_total_main_length = bdt.GetBadReco3().br3_2_total_main_length;
        br3_2_total_length = bdt.GetBadReco3().br3_2_total_length;
        br3_2_other_fid = bdt.GetBadReco3().br3_2_other_fid;
        br3_3_v_flag = bdt.GetBadReco3().br3_3_v_flag;
        br3_3_v_energy = bdt.GetBadReco3().br3_3_v_energy;
        br3_3_v_angle = bdt.GetBadReco3().br3_3_v_angle;
        br3_3_v_dir_length = bdt.GetBadReco3().br3_3_v_dir_length;
        br3_3_v_length = bdt.GetBadReco3().br3_3_v_length;
        br3_4_flag = bdt.GetBadReco3().br3_4_flag;
        br3_4_acc_length = bdt.GetBadReco3().br3_4_acc_length;
        br3_4_total_length = bdt.GetBadReco3().br3_4_total_length;
        br3_4_energy = bdt.GetBadReco3().br3_4_energy;
        br3_5_v_flag = bdt.GetBadReco3().br3_5_v_flag;
        br3_5_v_dir_length = bdt.GetBadReco3().br3_5_v_dir_length;
        br3_5_v_total_length = bdt.GetBadReco3().br3_5_v_total_length;
        br3_5_v_flag_avoid_muon_check = bdt.GetBadReco3().br3_5_v_flag_avoid_muon_check;
        br3_5_v_n_seg = bdt.GetBadReco3().br3_5_v_n_seg;
        br3_5_v_angle = bdt.GetBadReco3().br3_5_v_angle;
        br3_5_v_sg_length = bdt.GetBadReco3().br3_5_v_sg_length;
        br3_5_v_energy = bdt.GetBadReco3().br3_5_v_energy;
        br3_5_v_n_main_segs = bdt.GetBadReco3().br3_5_v_n_main_segs;
        br3_5_v_n_segs = bdt.GetBadReco3().br3_5_v_n_segs;
        br3_5_v_shower_main_length = bdt.GetBadReco3().br3_5_v_shower_main_length;
        br3_5_v_shower_total_length = bdt.GetBadReco3().br3_5_v_shower_total_length;
        br3_6_v_flag = bdt.GetBadReco3().br3_6_v_flag;
        br3_6_v_angle = bdt.GetBadReco3().br3_6_v_angle;
        br3_6_v_angle1 = bdt.GetBadReco3().br3_6_v_angle1;
        br3_6_v_flag_shower_trajectory = bdt.GetBadReco3().br3_6_v_flag_shower_trajectory;
        br3_6_v_direct_length = bdt.GetBadReco3().br3_6_v_direct_length;
        br3_6_v_length = bdt.GetBadReco3().br3_6_v_length;
        br3_6_v_n_other_vtx_segs = bdt.GetBadReco3().br3_6_v_n_other_vtx_segs;
        br3_6_v_energy = bdt.GetBadReco3().br3_6_v_energy;
        br3_7_flag = bdt.GetBadReco3().br3_7_flag;
        br3_7_energy = bdt.GetBadReco3().br3_7_energy;
        br3_7_min_angle = bdt.GetBadReco3().br3_7_min_angle;
        br3_7_sg_length = bdt.GetBadReco3().br3_7_sg_length;
        br3_7_shower_main_length = bdt.GetBadReco3().br3_7_shower_main_length;
        br3_8_flag = bdt.GetBadReco3().br3_8_flag;
        br3_8_max_dQ_dx = bdt.GetBadReco3().br3_8_max_dQ_dx;
        br3_8_energy = bdt.GetBadReco3().br3_8_energy;
        br3_8_n_main_segs = bdt.GetBadReco3().br3_8_n_main_segs;
        br3_8_shower_main_length = bdt.GetBadReco3().br3_8_shower_main_length;
        br3_8_shower_length = bdt.GetBadReco3().br3_8_shower_length;
        br_filled = bdt.GetBadReco4().br_filled;
        br4_flag = bdt.GetBadReco4().br4_flag;
        br4_1_flag = bdt.GetBadReco4().br4_1_flag;
        br4_1_shower_main_length = bdt.GetBadReco4().br4_1_shower_main_length;
        br4_1_shower_total_length = bdt.GetBadReco4().br4_1_shower_total_length;
        br4_1_min_dis = bdt.GetBadReco4().br4_1_min_dis;
        br4_1_energy = bdt.GetBadReco4().br4_1_energy;
        br4_1_flag_avoid_muon_check = bdt.GetBadReco4().br4_1_flag_avoid_muon_check;
        br4_1_n_vtx_segs = bdt.GetBadReco4().br4_1_n_vtx_segs;
        br4_1_n_main_segs = bdt.GetBadReco4().br4_1_n_main_segs;
        br4_2_flag = bdt.GetBadReco4().br4_2_flag;
        br4_2_ratio_45 = bdt.GetBadReco4().br4_2_ratio_45;
        br4_2_ratio_35 = bdt.GetBadReco4().br4_2_ratio_35;
        br4_2_ratio_25 = bdt.GetBadReco4().br4_2_ratio_25;
        br4_2_ratio_15 = bdt.GetBadReco4().br4_2_ratio_15;
        br4_2_energy = bdt.GetBadReco4().br4_2_energy;
        br4_2_ratio1_45 = bdt.GetBadReco4().br4_2_ratio1_45;
        br4_2_ratio1_35 = bdt.GetBadReco4().br4_2_ratio1_35;
        br4_2_ratio1_25 = bdt.GetBadReco4().br4_2_ratio1_25;
        br4_2_ratio1_15 = bdt.GetBadReco4().br4_2_ratio1_15;
        br4_2_iso_angle = bdt.GetBadReco4().br4_2_iso_angle;
        br4_2_iso_angle1 = bdt.GetBadReco4().br4_2_iso_angle1;
        br4_2_angle = bdt.GetBadReco4().br4_2_angle;
        tro_flag = bdt.GetTrackOverCluster().tro_flag;
        tro_1_v_flag = bdt.GetTrackOverCluster().tro_1_v_flag;
        tro_1_v_particle_type = bdt.GetTrackOverCluster().tro_1_v_particle_type;
        tro_1_v_flag_dir_weak = bdt.GetTrackOverCluster().tro_1_v_flag_dir_weak;
        tro_1_v_min_dis = bdt.GetTrackOverCluster().tro_1_v_min_dis;
        tro_1_v_sg1_length = bdt.GetTrackOverCluster().tro_1_v_sg1_length;
        tro_1_v_shower_main_length = bdt.GetTrackOverCluster().tro_1_v_shower_main_length;
        tro_1_v_max_n_vtx_segs = bdt.GetTrackOverCluster().tro_1_v_max_n_vtx_segs;
        tro_1_v_tmp_length = bdt.GetTrackOverCluster().tro_1_v_tmp_length;
        tro_1_v_medium_dQ_dx = bdt.GetTrackOverCluster().tro_1_v_medium_dQ_dx;
        tro_1_v_dQ_dx_cut = bdt.GetTrackOverCluster().tro_1_v_dQ_dx_cut;
        tro_1_v_flag_shower_topology = bdt.GetTrackOverCluster().tro_1_v_flag_shower_topology;
        tro_2_v_flag = bdt.GetTrackOverCluster().tro_2_v_flag;
        tro_2_v_energy = bdt.GetTrackOverCluster().tro_2_v_energy;
        tro_2_v_stem_length = bdt.GetTrackOverCluster().tro_2_v_stem_length;
        tro_2_v_iso_angle = bdt.GetTrackOverCluster().tro_2_v_iso_angle;
        tro_2_v_max_length = bdt.GetTrackOverCluster().tro_2_v_max_length;
        tro_2_v_angle = bdt.GetTrackOverCluster().tro_2_v_angle;
        tro_3_flag = bdt.GetTrackOverCluster().tro_3_flag;
        tro_3_stem_length = bdt.GetTrackOverCluster().tro_3_stem_length;
        tro_3_n_muon_segs = bdt.GetTrackOverCluster().tro_3_n_muon_segs;
        tro_3_energy = bdt.GetTrackOverCluster().tro_3_energy;
        tro_4_v_flag = bdt.GetTrackOverCluster().tro_4_v_flag;
        tro_4_v_dir2_mag = bdt.GetTrackOverCluster().tro_4_v_dir2_mag;
        tro_4_v_angle = bdt.GetTrackOverCluster().tro_4_v_angle;
        tro_4_v_angle1 = bdt.GetTrackOverCluster().tro_4_v_angle1;
        tro_4_v_angle2 = bdt.GetTrackOverCluster().tro_4_v_angle2;
        tro_4_v_length = bdt.GetTrackOverCluster().tro_4_v_length;
        tro_4_v_length1 = bdt.GetTrackOverCluster().tro_4_v_length1;
        tro_4_v_medium_dQ_dx = bdt.GetTrackOverCluster().tro_4_v_medium_dQ_dx;
        tro_4_v_end_dQ_dx = bdt.GetTrackOverCluster().tro_4_v_end_dQ_dx;
        tro_4_v_energy = bdt.GetTrackOverCluster().tro_4_v_energy;
        tro_4_v_shower_main_length = bdt.GetTrackOverCluster().tro_4_v_shower_main_length;
        tro_4_v_flag_shower_trajectory = bdt.GetTrackOverCluster().tro_4_v_flag_shower_trajectory;
        tro_5_v_flag = bdt.GetTrackOverCluster().tro_5_v_flag;
        tro_5_v_max_angle = bdt.GetTrackOverCluster().tro_5_v_max_angle;
        tro_5_v_min_angle = bdt.GetTrackOverCluster().tro_5_v_min_angle;
        tro_5_v_max_length = bdt.GetTrackOverCluster().tro_5_v_max_length;
        tro_5_v_iso_angle = bdt.GetTrackOverCluster().tro_5_v_iso_angle;
        tro_5_v_n_vtx_segs = bdt.GetTrackOverCluster().tro_5_v_n_vtx_segs;
        tro_5_v_min_count = bdt.GetTrackOverCluster().tro_5_v_min_count;
        tro_5_v_max_count = bdt.GetTrackOverCluster().tro_5_v_max_count;
        tro_5_v_energy = bdt.GetTrackOverCluster().tro_5_v_energy;
        hol_flag = bdt.GetHighEoverlap().hol_flag;
        hol_1_flag = bdt.GetHighEoverlap().hol_1_flag;
        hol_1_n_valid_tracks = bdt.GetHighEoverlap().hol_1_n_valid_tracks;
        hol_1_min_angle = bdt.GetHighEoverlap().hol_1_min_angle;
        hol_1_energy = bdt.GetHighEoverlap().hol_1_energy;
        hol_1_flag_all_shower = bdt.GetHighEoverlap().hol_1_flag_all_shower;
        hol_1_min_length = bdt.GetHighEoverlap().hol_1_min_length;
        hol_2_flag = bdt.GetHighEoverlap().hol_2_flag;
        hol_2_min_angle = bdt.GetHighEoverlap().hol_2_min_angle;
        hol_2_medium_dQ_dx = bdt.GetHighEoverlap().hol_2_medium_dQ_dx;
        hol_2_ncount = bdt.GetHighEoverlap().hol_2_ncount;
        hol_2_energy = bdt.GetHighEoverlap().hol_2_energy;
        lol_flag = bdt.GetLowEoverlap().lol_flag;
        lol_1_v_flag = bdt.GetLowEoverlap().lol_1_v_flag;
        lol_1_v_energy = bdt.GetLowEoverlap().lol_1_v_energy;
        lol_1_v_vtx_n_segs = bdt.GetLowEoverlap().lol_1_v_vtx_n_segs;
        lol_1_v_nseg = bdt.GetLowEoverlap().lol_1_v_nseg;
        lol_1_v_angle = bdt.GetLowEoverlap().lol_1_v_angle;
        lol_2_v_flag = bdt.GetLowEoverlap().lol_2_v_flag;
        lol_2_v_length = bdt.GetLowEoverlap().lol_2_v_length;
        lol_2_v_angle = bdt.GetLowEoverlap().lol_2_v_angle;
        lol_2_v_type = bdt.GetLowEoverlap().lol_2_v_type;
        lol_2_v_vtx_n_segs = bdt.GetLowEoverlap().lol_2_v_vtx_n_segs;
        lol_2_v_energy = bdt.GetLowEoverlap().lol_2_v_energy;
        lol_2_v_shower_main_length = bdt.GetLowEoverlap().lol_2_v_shower_main_length;
        lol_2_v_flag_dir_weak = bdt.GetLowEoverlap().lol_2_v_flag_dir_weak;
        lol_3_flag = bdt.GetLowEoverlap().lol_3_flag;
        lol_3_angle_beam = bdt.GetLowEoverlap().lol_3_angle_beam;
        lol_3_n_valid_tracks = bdt.GetLowEoverlap().lol_3_n_valid_tracks;
        lol_3_min_angle = bdt.GetLowEoverlap().lol_3_min_angle;
        lol_3_vtx_n_segs = bdt.GetLowEoverlap().lol_3_vtx_n_segs;
        lol_3_energy = bdt.GetLowEoverlap().lol_3_energy;
        lol_3_shower_main_length = bdt.GetLowEoverlap().lol_3_shower_main_length;
        lol_3_n_out = bdt.GetLowEoverlap().lol_3_n_out;
        lol_3_n_sum = bdt.GetLowEoverlap().lol_3_n_sum;
        cosmict_flag_1 = bdt.GetMajorCosmicTagger().cosmict_flag_1; // fiducial volume vertex
        cosmict_flag_2 = bdt.GetMajorCosmicTagger().cosmict_flag_2;  // single muon
        cosmict_flag_3 = bdt.GetMajorCosmicTagger().cosmict_flag_3;  // single muon (long)
        cosmict_flag_4 = bdt.GetMajorCosmicTagger().cosmict_flag_4;  // kinematics muon
        cosmict_flag_5 = bdt.GetMajorCosmicTagger().cosmict_flag_5; // kinematics muon (long)
        cosmict_flag_6 = bdt.GetMajorCosmicTagger().cosmict_flag_6; // special ...
        cosmict_flag_7 = bdt.GetMajorCosmicTagger().cosmict_flag_7;  // muon+ michel
        cosmict_flag_8 = bdt.GetMajorCosmicTagger().cosmict_flag_8;  // muon + michel + special
        cosmict_flag_9 = bdt.GetMajorCosmicTagger().cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
        cosmict_flag_10 = bdt.GetMajorCosmicTagger().cosmict_flag_10;  // front upstream (dirt)
        cosmict_flag = bdt.GetMajorCosmicTagger().cosmict_flag;
        cosmict_2_filled = bdt.GetMajorCosmicTagger().cosmict_2_filled;
        cosmict_2_particle_type = bdt.GetMajorCosmicTagger().cosmict_2_particle_type;
        cosmict_2_n_muon_tracks = bdt.GetMajorCosmicTagger().cosmict_2_n_muon_tracks;
        cosmict_2_total_shower_length = bdt.GetMajorCosmicTagger().cosmict_2_total_shower_length;
        cosmict_2_flag_inside = bdt.GetMajorCosmicTagger().cosmict_2_flag_inside;
        cosmict_2_angle_beam = bdt.GetMajorCosmicTagger().cosmict_2_angle_beam;
        cosmict_2_flag_dir_weak = bdt.GetMajorCosmicTagger().cosmict_2_flag_dir_weak;
        cosmict_2_dQ_dx_end = bdt.GetMajorCosmicTagger().cosmict_2_dQ_dx_end;
        cosmict_2_dQ_dx_front = bdt.GetMajorCosmicTagger().cosmict_2_dQ_dx_front;
        cosmict_2_theta = bdt.GetMajorCosmicTagger().cosmict_2_theta;
        cosmict_2_phi = bdt.GetMajorCosmicTagger().cosmict_2_phi;
        cosmict_2_valid_tracks = bdt.GetMajorCosmicTagger().cosmict_2_valid_tracks;
        cosmict_3_filled = bdt.GetMajorCosmicTagger().cosmict_3_filled;
        cosmict_3_flag_inside = bdt.GetMajorCosmicTagger().cosmict_3_flag_inside;
        cosmict_3_angle_beam = bdt.GetMajorCosmicTagger().cosmict_3_angle_beam;
        cosmict_3_flag_dir_weak = bdt.GetMajorCosmicTagger().cosmict_3_flag_dir_weak;
        cosmict_3_dQ_dx_end = bdt.GetMajorCosmicTagger().cosmict_3_dQ_dx_end;
        cosmict_3_dQ_dx_front = bdt.GetMajorCosmicTagger().cosmict_3_dQ_dx_front;
        cosmict_3_theta = bdt.GetMajorCosmicTagger().cosmict_3_theta;
        cosmict_3_phi = bdt.GetMajorCosmicTagger().cosmict_3_phi;
        cosmict_3_valid_tracks = bdt.GetMajorCosmicTagger().cosmict_3_valid_tracks;
        cosmict_4_filled = bdt.GetMajorCosmicTagger().cosmict_4_filled;
        cosmict_4_flag_inside = bdt.GetMajorCosmicTagger().cosmict_4_flag_inside;
        cosmict_4_angle_beam = bdt.GetMajorCosmicTagger().cosmict_4_angle_beam;
        cosmict_4_connected_showers = bdt.GetMajorCosmicTagger().cosmict_4_connected_showers;  // need to be careful about the nueCC ...
        cosmict_5_filled = bdt.GetMajorCosmicTagger().cosmict_5_filled;
        cosmict_5_flag_inside = bdt.GetMajorCosmicTagger().cosmict_5_flag_inside;
        cosmict_5_angle_beam = bdt.GetMajorCosmicTagger().cosmict_5_angle_beam;
        cosmict_5_connected_showers = bdt.GetMajorCosmicTagger().cosmict_5_connected_showers;
        cosmict_6_filled = bdt.GetMajorCosmicTagger().cosmict_6_filled;
        cosmict_6_flag_dir_weak = bdt.GetMajorCosmicTagger().cosmict_6_flag_dir_weak;
        cosmict_6_flag_inside = bdt.GetMajorCosmicTagger().cosmict_6_flag_inside;
        cosmict_6_angle = bdt.GetMajorCosmicTagger().cosmict_6_angle;
        cosmict_7_filled = bdt.GetMajorCosmicTagger().cosmict_7_filled;
        cosmict_7_flag_sec = bdt.GetMajorCosmicTagger().cosmict_7_flag_sec;
        cosmict_7_n_muon_tracks = bdt.GetMajorCosmicTagger().cosmict_7_n_muon_tracks;
        cosmict_7_total_shower_length = bdt.GetMajorCosmicTagger().cosmict_7_total_shower_length;
        cosmict_7_flag_inside = bdt.GetMajorCosmicTagger().cosmict_7_flag_inside;
        cosmict_7_angle_beam = bdt.GetMajorCosmicTagger().cosmict_7_angle_beam;
        cosmict_7_flag_dir_weak = bdt.GetMajorCosmicTagger().cosmict_7_flag_dir_weak;
        cosmict_7_dQ_dx_end = bdt.GetMajorCosmicTagger().cosmict_7_dQ_dx_end;
        cosmict_7_dQ_dx_front = bdt.GetMajorCosmicTagger().cosmict_7_dQ_dx_front;
        cosmict_7_theta = bdt.GetMajorCosmicTagger().cosmict_7_theta;
        cosmict_7_phi = bdt.GetMajorCosmicTagger().cosmict_7_phi;
        cosmict_8_filled = bdt.GetMajorCosmicTagger().cosmict_8_filled;
        cosmict_8_flag_out = bdt.GetMajorCosmicTagger().cosmict_8_flag_out;
        cosmict_8_muon_length = bdt.GetMajorCosmicTagger().cosmict_8_muon_length;
        cosmict_8_acc_length = bdt.GetMajorCosmicTagger().cosmict_8_acc_length;
        cosmict_10_flag_inside = bdt.GetMajorCosmicTagger().cosmict_10_flag_inside;
        cosmict_10_vtx_z = bdt.GetMajorCosmicTagger().cosmict_10_vtx_z;
        cosmict_10_flag_shower = bdt.GetMajorCosmicTagger().cosmict_10_flag_shower;
        cosmict_10_flag_dir_weak = bdt.GetMajorCosmicTagger().cosmict_10_flag_dir_weak;
        cosmict_10_angle_beam = bdt.GetMajorCosmicTagger().cosmict_10_angle_beam;
        cosmict_10_length = bdt.GetMajorCosmicTagger().cosmict_10_length;
        numu_cc_flag = bdt.GetNumuCCTagger().numu_cc_flag;
        numu_cc_flag_1 = bdt.GetNumuCCTagger().numu_cc_flag_1;
        numu_cc_1_particle_type = bdt.GetNumuCCTagger().numu_cc_1_particle_type;
        numu_cc_1_length = bdt.GetNumuCCTagger().numu_cc_1_length;
        numu_cc_1_medium_dQ_dx = bdt.GetNumuCCTagger().numu_cc_1_medium_dQ_dx;
        numu_cc_1_dQ_dx_cut = bdt.GetNumuCCTagger().numu_cc_1_dQ_dx_cut;
        numu_cc_1_direct_length = bdt.GetNumuCCTagger().numu_cc_1_direct_length;
        numu_cc_1_n_daughter_tracks = bdt.GetNumuCCTagger().numu_cc_1_n_daughter_tracks;
        numu_cc_1_n_daughter_all = bdt.GetNumuCCTagger().numu_cc_1_n_daughter_all;
        numu_cc_flag_2 = bdt.GetNumuCCTagger().numu_cc_flag_2;
        numu_cc_2_length = bdt.GetNumuCCTagger().numu_cc_2_length;
        numu_cc_2_total_length = bdt.GetNumuCCTagger().numu_cc_2_total_length;
        numu_cc_2_n_daughter_tracks = bdt.GetNumuCCTagger().numu_cc_2_n_daughter_tracks;
        numu_cc_2_n_daughter_all = bdt.GetNumuCCTagger().numu_cc_2_n_daughter_all;
        numu_cc_flag_3 = bdt.GetNumuCCTagger().numu_cc_flag_3;
        numu_cc_3_particle_type = bdt.GetNumuCCTagger().numu_cc_3_particle_type;
        numu_cc_3_max_length = bdt.GetNumuCCTagger().numu_cc_3_max_length;
        numu_cc_3_acc_track_length = bdt.GetNumuCCTagger().numu_cc_3_acc_track_length;
        numu_cc_3_max_length_all = bdt.GetNumuCCTagger().numu_cc_3_max_length_all;
        numu_cc_3_max_muon_length = bdt.GetNumuCCTagger().numu_cc_3_max_muon_length;
        numu_cc_3_n_daughter_tracks = bdt.GetNumuCCTagger().numu_cc_3_n_daughter_tracks;
        numu_cc_3_n_daughter_all = bdt.GetNumuCCTagger().numu_cc_3_n_daughter_all;
        cosmict_2_4_score = bdt.GetBDTscores().cosmict_2_4_score;
        cosmict_3_5_score = bdt.GetBDTscores().cosmict_3_5_score;
        cosmict_6_score = bdt.GetBDTscores().cosmict_6_score;
        cosmict_7_score = bdt.GetBDTscores().cosmict_7_score;
        cosmict_8_score = bdt.GetBDTscores().cosmict_8_score;
        cosmict_10_score = bdt.GetBDTscores().cosmict_10_score;
        numu_1_score = bdt.GetBDTscores().numu_1_score;
        numu_2_score = bdt.GetBDTscores().numu_2_score;
        numu_3_score = bdt.GetBDTscores().numu_3_score;
        cosmict_score = bdt.GetBDTscores().cosmict_score;
        numu_score = bdt.GetBDTscores().numu_score;
        mipid_score = bdt.GetBDTscores().mipid_score;
        gap_score = bdt.GetBDTscores().gap_score;
        hol_lol_score = bdt.GetBDTscores().hol_lol_score;
        cme_anc_score = bdt.GetBDTscores().cme_anc_score;
        mgo_mgt_score = bdt.GetBDTscores().mgo_mgt_score;
        br1_score = bdt.GetBDTscores().br1_score;
        br3_score = bdt.GetBDTscores().br3_score;
        br3_3_score = bdt.GetBDTscores().br3_3_score;
        br3_5_score = bdt.GetBDTscores().br3_5_score;
        br3_6_score = bdt.GetBDTscores().br3_6_score;
        stemdir_br2_score = bdt.GetBDTscores().stemdir_br2_score;
        trimuon_score = bdt.GetBDTscores().trimuon_score;
        br4_tro_score = bdt.GetBDTscores().br4_tro_score;
        mipquality_score = bdt.GetBDTscores().mipquality_score;
        pio_1_score = bdt.GetBDTscores().pio_1_score;
        pio_2_score = bdt.GetBDTscores().pio_2_score;
        stw_spt_score = bdt.GetBDTscores().stw_spt_score;
        vis_1_score = bdt.GetBDTscores().vis_1_score;
        vis_2_score = bdt.GetBDTscores().vis_2_score;
        stw_2_score = bdt.GetBDTscores().stw_2_score;
        stw_3_score = bdt.GetBDTscores().stw_3_score;
        stw_4_score = bdt.GetBDTscores().stw_4_score;
        sig_1_score = bdt.GetBDTscores().sig_1_score;
        sig_2_score = bdt.GetBDTscores().sig_2_score;
        lol_1_score = bdt.GetBDTscores().lol_1_score;
        lol_2_score = bdt.GetBDTscores().lol_2_score;
        tro_1_score = bdt.GetBDTscores().tro_1_score;
        tro_2_score = bdt.GetBDTscores().tro_2_score;
        tro_4_score = bdt.GetBDTscores().tro_4_score;
        tro_5_score = bdt.GetBDTscores().tro_5_score;
        nue_score = bdt.GetBDTscores().nue_score;

}

void WireCellAnaTree::ReadKINEvar(nsm::NuSelectionKINE const& kine)
{
        kine_reco_Enu = kine.GetKineInfo().kine_reco_Enu;
        kine_reco_add_energy = kine.GetKineInfo().kine_reco_add_energy;
        kine_energy_particle = kine.GetKineInfo().kine_energy_particle;
        kine_energy_info = kine.GetKineInfo().kine_energy_info;
        kine_particle_type = kine.GetKineInfo().kine_particle_type;
        kine_energy_included = kine.GetKineInfo().kine_energy_included;
        kine_pio_mass = kine.GetKineInfo().kine_pio_mass;
        kine_pio_flag = kine.GetKineInfo().kine_pio_flag;
        kine_pio_vtx_dis = kine.GetKineInfo().kine_pio_vtx_dis;
        kine_pio_energy_1 = kine.GetKineInfo().kine_pio_energy_1;
        kine_pio_theta_1 = kine.GetKineInfo().kine_pio_theta_1;
        kine_pio_phi_1 = kine.GetKineInfo().kine_pio_phi_1;
        kine_pio_dis_1 = kine.GetKineInfo().kine_pio_dis_1;
        kine_pio_energy_2 = kine.GetKineInfo().kine_pio_energy_2;
        kine_pio_theta_2 = kine.GetKineInfo().kine_pio_theta_2;
        kine_pio_phi_2 = kine.GetKineInfo().kine_pio_phi_2;
        kine_pio_dis_2 = kine.GetKineInfo().kine_pio_dis_2;
        kine_pio_angle = kine.GetKineInfo().kine_pio_angle;
}


void WireCellAnaTree::nsbeamtiming(art::Event const& e)
{
  double max[32],time[32];
  bool Sat[32];
  for(int i=0; i<32; i++){Sat[i]=false;}
  double BeamT0 = -99999.;
  std::vector<int>* pmtid = new std::vector<int>;
  getPMTwf(e,max,time,Sat);
  const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
  for (int i=0; i<32; i++){
    calib[i] = gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain");
  }
  if(!fMC) {BeamT0 = getBeamWF(e);}
  else{BeamT0 = 0;}//no RWM signal simulated for now, so just set this to 0
  Float_t x =        f_reco_nuvtxX;
  Float_t y =        f_reco_nuvtxY;
  Float_t z =        f_reco_nuvtxZ;
  if(f_ns_time_useSSMvtx){
    if(ssm_vtxX!=-999 && ssm_vtxY!=-999 && ssm_vtxZ!=-999 && ssm_kine_energy>0){
      x = ssm_vtxX;
      y = ssm_vtxY;
      z = ssm_vtxZ;
    }
  }
  std::vector<float> *sps_x = new std::vector<float>;
  std::vector<float> *sps_y = new std::vector<float>;
  std::vector<float> *sps_z = new std::vector<float>;
  std::vector<float> *sps_t = new std::vector<float>;
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  if (! e.getByLabel(fPFInputTag, particleHandle)) return;
  std::vector< art::Ptr<simb::MCParticle> > particles;
  art::fill_ptr_vector(particles, particleHandle);
  //origional method
  if(!f_ns_time_usePID){
    for (auto const& particle: particles){
      for (uint i_pos=0; i_pos<particle->NumberTrajectoryPoints(); i_pos++){
    	  const TLorentzVector& pos = particle->Position(i_pos);
          sps_x->push_back(pos.X());
          sps_y->push_back(pos.Y());
          sps_z->push_back(pos.Z());
        }
    }
  }  
  //method using linear extrapolation and PID
  else{
    std::map<int,std::tuple< std::vector<float>*,std::vector<float>*,std::vector<float>*,std::vector<float>* > >my_particle_times;
    //<id,<x,y,z,t>>

    //only do primary particles first
    for (auto const& particle: particles){
      int mother = particle->Mother();
      if(mother!=0) {continue;} 
      int id = particle->TrackId();
      auto this_id = my_particle_times.find(id);      
      if (this_id==my_particle_times.end()){
        my_particle_times[id] = get_extrapolated_times(particle, 0);
      }

      //now find the direct daughters of this primary particle
      std::vector< std::tuple<int,art::Ptr<simb::MCParticle>,double> > daughters; //<id,particle,mother_time>
      for(auto const& daughter_particle: particles){
        int mother_of_daughter = daughter_particle->Mother();
        int daughter_id = daughter_particle->TrackId();
        if (mother_of_daughter == id) {
          double mother_time = std::get<3>(my_particle_times[id])->back();
          daughters.push_back(std::make_tuple(daughter_id,daughter_particle,mother_time));
        }
      }

      //Now add this daughter and find the daughters of the daughter
      //keep going untill we exahust all daughters of daughters
      while (daughters.size() > 0){
        int daughter_id = std::get<0>(daughters.front());
        auto this_daughter_id = my_particle_times.find(daughter_id);
        //skip is we have already added it
	if (this_daughter_id!=my_particle_times.end()){
          daughters.erase(daughters.begin());
          continue;
        }
        my_particle_times[daughter_id] = get_extrapolated_times( std::get<1>(daughters.front()), std::get<2>(daughters.front()));
        for(auto const& daughter_daughter_particle: particles){
          int mother_of_daughter_daughter_id = daughter_daughter_particle->Mother();
          int daughter_daughter_id = daughter_daughter_particle->TrackId();
          auto this_daughter_daughter_id = my_particle_times.find(daughter_daughter_id);
	  if (daughter_id == mother_of_daughter_daughter_id && this_daughter_daughter_id==my_particle_times.end()) {
            //new daughter that we need to check on, its a daughter of this daughter.
	    //Set its mother time according to the time we just added for the daughter we were working on.
            daughters.push_back(std::make_tuple(daughter_daughter_id,daughter_daughter_particle,std::get<3>(my_particle_times[daughter_id])->back())); 
	  }
        } //no more daugters of the current daughter
        daughters.erase(daughters.begin());
      }//no more daughter or duaghters of daughters
    }//end loop over all primary particles

    //unpack this for the rest of the code
    for (auto my_particle = my_particle_times.begin(); my_particle != my_particle_times.end(); my_particle++) {
      std::vector<float>* _x = std::get<0>(my_particle->second);
      std::vector<float>* _y = std::get<1>(my_particle->second);
      std::vector<float>* _z = std::get<2>(my_particle->second);
      std::vector<float>* _t = std::get<3>(my_particle->second);
      for(uint point=0; point < _x->size(); point++){
        sps_x->push_back(_x->at(point));
        sps_y->push_back(_y->at(point));
        sps_z->push_back(_z->at(point));
        sps_t->push_back(_t->at(point));	 
      }
    }
  }
  
  double PMT0[3]={-11.4545, -28.625, 990.356};  double PMT1[3]={-11.4175, 27.607, 989.712};
  double PMT2[3]={-11.7755, -56.514, 951.865};  double PMT3[3]={-11.6415, 55.313, 951.861};
  double PMT4[3]={-12.0585, -56.309, 911.939};  double PMT5[3]={-11.8345, 55.822, 911.065};
  double PMT6[3]={-12.1765, -0.722, 865.599};   double PMT7[3]={-12.3045, -0.502, 796.208};
  double PMT8[3]={-12.6045, -56.284, 751.905};  double PMT9[3]={-12.5405, 55.625, 751.884};
  double PMT10[3]={-12.6125, -56.408, 711.274}; double PMT11[3]={-12.6615, 55.8, 711.073};
  double PMT12[3]={-12.6245, -0.051, 664.203};  double PMT13[3]={-12.6515, -0.549, 585.284};
  double PMT14[3]={-12.8735, 55.822, 540.929};  double PMT15[3]={-12.6205, -56.205, 540.616};
  double PMT16[3]={-12.5945, -56.323, 500.221}; double PMT17[3]={-12.9835, 55.771, 500.134};
  double PMT18[3]={-12.6185, -0.875, 453.096};  double PMT19[3]={-13.0855, -0.706, 373.839};
  double PMT20[3]={-12.6485, -57.022, 328.341}; double PMT21[3]={-13.1865, 54.693, 328.212};
  double PMT22[3]={-13.4175, 54.646, 287.976};  double PMT23[3]={-13.0075, -56.261, 287.639};
  double PMT24[3]={-13.1505, -0.829, 242.014};  double PMT25[3]={-13.4415, -0.303, 173.743};
  double PMT26[3]={-13.3965, 55.249, 128.354};  double PMT27[3]={-13.2784, -56.203, 128.18};
  double PMT28[3]={-13.2375, -56.615, 87.8695}; double PMT29[3]={-13.5415, 55.249, 87.7605};
  double PMT30[3]={-13.4345, 27.431, 51.1015};  double PMT31[3]={-13.1525, -28.576, 50.4745};
  double PMT[32][3];    for(int j=0; j<3; j++){ PMT[30][j]=PMT30[j]; PMT[31][j]=PMT31[j];
  PMT[0][j]=PMT0[j];   PMT[10][j]=PMT10[j]; PMT[20][j]=PMT20[j]; PMT[1][j]=PMT1[j];   PMT[11][j]=PMT11[j];
  PMT[21][j]=PMT21[j]; PMT[2][j]=PMT2[j];   PMT[12][j]=PMT12[j]; PMT[22][j]=PMT22[j]; PMT[3][j]=PMT3[j];
  PMT[13][j]=PMT13[j]; PMT[23][j]=PMT23[j]; PMT[4][j]=PMT4[j];   PMT[14][j]=PMT14[j]; PMT[24][j]=PMT24[j];
  PMT[5][j]=PMT5[j];   PMT[15][j]=PMT15[j]; PMT[25][j]=PMT25[j]; PMT[6][j]=PMT6[j];   PMT[16][j]=PMT16[j];
  PMT[18][j]=PMT18[j]; PMT[28][j]=PMT28[j]; PMT[9][j]=PMT9[j];   PMT[19][j]=PMT19[j]; PMT[29][j]=PMT29[j];
  PMT[26][j]=PMT26[j]; PMT[7][j]=PMT7[j];   PMT[17][j]=PMT17[j]; PMT[27][j]=PMT27[j]; PMT[8][j]=PMT8[j];}
  double offset[32]={1.03002, -5.18104, -2.11164, -5.99395, -1.25798, 0.633079, 2.87666, 2.21969, 0.885092, 2.35423,
    -1.63039, -1.83775, -0.859883, 3.4741, 1.84833, 1.58233, -2.71783, 0, 3.18776, 0.982666, 0.728438, 0.280592, -5.27068,
    -3.27857, -1.41196, 1.59643, 1.41425, -1.62682, -2.55772, 1.49136, -0.522791, 0.974533};
  if(fMC){
    for(int i=0; i<32; i++){offset[i]=0;}//no need to apply the additional pmt calibration to the MC
  }
  //================================================================================================================
  double gap=18.936;
  double MaxLim=2.5;
  std::vector<int> N_pmt;
  double ccnd1, ccnd2,ccnd3, ccnd4;
  double Ph_Tot, RWM_T, nuToF, DPh,DLh, tPhelp,tp, tDPhelp,tDP, tDLhelp,tDL;
  double Med_TT3=-9999.;
  double TT_merged = -9999.;
  nuToF=0;
  //===================================================================================================================
  //===================================================================================================================
  Ph_Tot=0.;
  N_pmt.clear();
  for(int q=0; q<32; q++){max[q]=max[q]/calib[q];
  //do not use a time cut in the NuMI case
  if((max[q]>MaxLim && q!=17 && q!=28) && ((time[q]>3000.0 && time[q]<5000.0) || fIsNuMI ) ){N_pmt.push_back(q); Ph_Tot=Ph_Tot+max[q]; pmtid->push_back(q);}}
  std::vector<double> timeProp = std::vector<double>(N_pmt.size(),0);
  std::vector<double> timeDP = std::vector<double>(N_pmt.size(),0);
  std::vector<double> timeDL = std::vector<double>(N_pmt.size(),0);
  //--------------------------------------------------------------------------------------------------------------------
  if(N_pmt.size()>2){
    RWM_T=BeamT0;
    double dist = z; //in BNB correct to front face of TPC, in NuMI correct to plane perpendicular to the beam
    if(fIsNuMI) {
      TVector3 target_dir(-0.46, -0.05, -0.885);
      double min_a = -122.86902944472968;  
      double min_b = 80.60659897339974; 
      double min_c = 59.34119182916038;
      dist = ( (min_a-x)*target_dir[0] + (min_b-y)*target_dir[1] + (min_c-z)*target_dir[2] ) / sqrt(target_dir[0]*target_dir[0] + target_dir[1]*target_dir[1] + target_dir[2]*target_dir[2] );
    }
    nuToF=dist*0.033356;
    for(uint i=0; i<N_pmt.size(); i++){
        tp=5000000000.0;
        tDL=5000000000.0;
        tDP=5000000000.0;
	for(uint j=0; j<sps_x->size(); j++){
          DPh=abs(sqrt(TMath::Power(x-sps_x->at(j),2)+TMath::Power(y-sps_y->at(j),2)+TMath::Power(z-sps_z->at(j),2)));
          DLh=abs(sqrt(TMath::Power(PMT[N_pmt.at(i)][0]-sps_x->at(j),2)+TMath::Power(PMT[N_pmt.at(i)][1]-sps_y->at(j),2)+TMath::Power(PMT[N_pmt.at(i)][2]-sps_z->at(j),2)));
          if(!f_ns_time_usePID){
            tPhelp=(DPh*fsol)+(DLh*0.0746);
            tDPhelp=(DPh*fsol);
            tDLhelp=(DLh*0.0746);
            if(f_ns_time_no_photon){
              tPhelp=(DPh*fsol);
              tDLhelp=0;
            }
	  }
	  else{
            tPhelp=sps_t->at(j)+(DLh*0.0746);
            tDPhelp=sps_t->at(j);
            tDLhelp=(DLh*0.0746);
	    if(f_ns_time_no_photon){
              tPhelp=sps_t->at(j);
              tDLhelp=0;
            }
	  }
  	  if(tPhelp<tp){
            tp=tPhelp;
            tDP=tDPhelp;
            tDL=tDLhelp;
          }
	}
        timeProp[i]=tp;
        timeDP[i]=tDP;
        timeDL[i]=tDL;
    }
    double TT3_array[32];
    //do not think we have to make this correction for NuMI, may need to revisit
    if(f_isrun3 && f_run>17200 && f_run<17400 && !fIsNuMI){if(RWM_T>5450){ f_shiftoffset=118.3;}}
    float RWM_offset = 5700.0 - f_shiftoffset;
    for(uint i=0; i<N_pmt.size(); i++){
      ccnd1= timeProp[i]*(f_ccnd1_a)-(f_ccnd1_b);
      ccnd2= max[N_pmt.at(i)]*(f_ccnd2_a)-(f_ccnd2_b);
      ccnd4= x*(f_ccnd4_a)-(f_ccnd4_b);//for x dependent correction
      if(Ph_Tot>150){ccnd3=f_ccnd3_a-f_ccnd3_b*Ph_Tot+f_ccnd3_c*Ph_Tot*Ph_Tot;}
      else{ccnd3=f_ccnd3_d;}

      
      //all the corrections
      TT3_array[i]=(time[N_pmt.at(i)])-RWM_T+RWM_offset-nuToF-timeProp[i]-offset[N_pmt.at(i)]+ccnd1+ccnd2+ccnd3+ccnd4;
      std::cout<<"TT3_array[i] "<<i<<" "<<TT3_array[i]<<std::endl;
    }
    Med_TT3=TMath::Median((Long64_t)N_pmt.size(),TT3_array);
    //Fill a 2d histogram with  TT3_array[i] vs max[N_pmt.at(i)] this is usefull to check for any errors
    for(uint i=0; i<N_pmt.size(); i++){
      H_TimeVsPh->Fill( TT3_array[i]-Med_TT3, max[N_pmt.at(i)]);
    }
  }
  f_evtTimeNS = Med_TT3;
  f_evtTimeNS_redk2nu = Med_TT3+f_redk2nu_deltatime;
  f_Ph_Tot = Ph_Tot;
  std::cout<<"f_evtTimeNS "<<f_evtTimeNS<<"  f_evtTimeNS_redk2nu "<<f_evtTimeNS_redk2nu<<std::endl;

  //Merge Peaks, shift may be incorrect for NuMi
  double Shift=3166.9;
  double TThelp=Med_TT3-Shift+gap*0.5;
  double bunches = 81;
  if(fIsNuMI){Shift=11567.87; gap=18.83; bunches=503;}
  //merge peaks
  if(TThelp>=0 && TThelp<gap*bunches){
    TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
  }
  else {TT_merged=-9999;}
  f_evtDeltaTimeNS = TT_merged;

  if(f_savepmt){
    f_PMT_ID = pmtid;
    f_RWM_Time = BeamT0;
    for(uint pmt=0; pmt<pmtid->size(); pmt++){
      int id = pmtid->at(pmt);
      f_PMT_TimeProp->push_back(timeProp[pmt]);
      f_PMT_TimeDP->push_back(timeDP[pmt]);
      f_PMT_TimeDL->push_back(timeDL[pmt]);
      f_PMT_Time->push_back(time[id]);
      f_PMT_Amp->push_back(max[id]);
      f_PMT_Sat->push_back(Sat[id]);
    }
  }

}
void WireCellAnaTree::getPMTwf(art::Event const& e, double maxP[32], double timeP[32], bool Sat[32])
{
    //set number of samples to look at, different for BNB and NuMI
    int samples = 500;
    int samples_64 = 5;
    if(fIsNuMI){samples=1500; samples_64=23;}//do we need to go this much further?
    //get waveforms
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle;
    if(!fMC) {e.getByLabel( "pmtreadout:OpdetBeamHighGain", wf_handle );}
    else{e.getByLabel( "mixer:OpdetBeamHighGain", wf_handle );}//MC in this instead
    if(!wf_handle.isValid()) {
      std::cout<<"BeamTiming: No pmtreadout:OpdetBeamHighGain!"<<std::endl;
      return;
    }
    //clear waveform
    //auto _ch = std::numeric_limits<unsigned int>::max();
    // auto wfsum = std::vector<double>(1500,0);
    auto _wf_v = std::vector< std::vector<double> >(32,std::vector<double>(1500,0));
    for (size_t i=0; i < wf_handle->size(); i++) {
      auto const wf = wf_handle->at(i);
      auto ch = wf.ChannelNumber();
      if (ch >= 32) continue;
      for (size_t n=0; n < wf.size(); n++){
        if (n < 1500){
          _wf_v[ch][n] = (wf_handle->at(i))[n];
	}
      }// for all channels
    }// for all waveforms
    //analyze waveforms
    //=======================================================================================================
    //=======================================================================================================
      double Help_wf_v[32][1500];
      double x_wf_v[1500], Raw_wf_v[1500], Base_wf_v[1500], Norm_wf_v[1500];
      double maxZ,max0,base;   int basebinmax,tick;
      double tca,tcb,tcc, TT[32], max[32];
      int TF,TB,tickF,tickB,FB, Nss,is;
      double  maxZhelp1,maxZhelp2,maxZhelp3, tickFit1,tickFit2;
      //Raw waveform saturation
      int saturation=4094;
      //Parameters for discrete saturated WF reconstruction (1st step)
      double Frac[100]={1.0, 0.951931, 0.93356, 0.838637, 0.719408, 0.701042, 0.565673, 0.44655, 0.39447, 0.352336, 0.28716, 0.245364, 0.216771, 0.194888, 0.178976, 0.16844, 0.161732, 0.157486, 0.154762, 0.153136, 0.152468, 0.152299, 0.152147, 0.151456, 0.150069, 0.148089, 0.145565, 0.142516, 0.139369, 0.136467, 0.133871, 0.13141, 0.128926, 0.126237, 0.12313, 0.119611, 0.115946, 0.112453, 0.111706, 0.109225, 0.106281, 0.103499, 0.100926, 0.0985215, 0.0961512, 0.0938219, 0.0917209, 0.0898817, 0.0882937, 0.0868338, 0.0854753, 0.0841353, 0.0827237, 0.081198, 0.0794796, 0.077565, 0.0756475, 0.0738863, 0.0722821, 0.0708473, 0.0695119, 0.0682504, 0.0672023, 0.0663549, 0.0656337, 0.0649918, 0.0645003, 0.0641535, 0.0638046, 0.0633435, 0.0627506, 0.0621379, 0.0615464, 0.0609178, 0.0601846, 0.0592098, 0.0580465, 0.0568861, 0.0559024, 0.0550731, 0.0541904, 0.0532532, 0.0524181, 0.0517606, 0.0512326, 0.0507392, 0.0502093, 0.0495968, 0.0488915, 0.0480173, 0.0470195, 0.0459744, 0.0448855, 0.0437359, 0.0425199, 0.0412832, 0.0400036, 0.038688, 0.0373173, 0.0358925};
      //double Delay[100]={0.0, 0.0686367, 0.0948574, 0.230854, 0.405159, 0.432604, 0.642806, 0.845613, 0.942555, 1.02616, 1.16775, 1.26923, 1.34523, 1.40802, 1.45675, 1.49068, 1.51306, 1.52757, 1.53702, 1.54272, 1.54507, 1.54567, 1.54621, 1.54865, 1.55359, 1.56069, 1.56984, 1.58104, 1.59279, 1.60379, 1.61377, 1.62337, 1.63319, 1.64398, 1.65665, 1.67129, 1.68688, 1.70208, 1.70537, 1.71644, 1.72981, 1.74271, 1.75486, 1.76644, 1.77806, 1.78969, 1.80036, 1.80987, 1.81819, 1.82594, 1.83324, 1.84054, 1.84831, 1.85684, 1.86658, 1.87763, 1.88891, 1.89947, 1.90926, 1.91815, 1.92656, 1.93462, 1.9414, 1.94695, 1.95171, 1.95598, 1.95928, 1.96162, 1.96398, 1.96711, 1.97117, 1.97539, 1.97951, 1.98391, 1.98909, 1.99605, 2.00448, 2.01302, 2.02038, 2.02665, 2.03342, 2.04069, 2.04726, 2.0525, 2.05674, 2.06073, 2.06505, 2.0701, 2.07597, 2.08334, 2.09188, 2.10099, 2.11065, 2.12106, 2.13232, 2.14403, 2.15646, 2.16957, 2.18362, 2.19868 };
      //Fit function parameters for saturated WF reconstruction (2nd step)
      double pLL[8]={3.93256,4.31002,2.44182,5.12491,0.830928,0.231375,50.9081,-2.69014};
    //=====================================================================================================
    //=======================================================================================================
    for(int q=0; q<32; q++){for(int i=0; i<samples; i++){Help_wf_v[q][i]=_wf_v[q][i];}}
    _wf_v.clear(); _wf_v.shrink_to_fit();
    //-----------------------------------------------------------------------------------------------------
    for(int q=0; q<32; q++){TT[q]=-9999.; max[q]=-9999.;
    maxZ=0.; max0=0.; base=0.; tick=0; tickB=0; tickF=0; TF=0; TB=0;
    //Getting raw waveform (Raw_wf_v[i]) only for i<500 since the beam window is between 3 and 5 us -> [i>3*64 && i<5*64], for NuMI, we need to go longer
    for(int i=0; i<samples; i++){x_wf_v[i]=i*1.0; Raw_wf_v[i]=Help_wf_v[q][i];}
    //Getting raw wf max amplitude and max amp tick
    for(int i=3*64; i<samples_64*64; i++){if(maxZ<Raw_wf_v[i]){maxZ=Raw_wf_v[i]; tick=i;}}
    //Baseline removal
    TH1F *basehelp= new TH1F("basehelp","basehelp",400, 1900,2200);
    basebinmax=0; for(int i=0; i<3*64; i++){basehelp->Fill(Raw_wf_v[i]);}
    basebinmax=basehelp->GetMaximumBin(); base=basehelp->GetXaxis()->GetBinCenter(basebinmax);
    basehelp->Delete();
    //Getting wf max amp after baseline removal (this is proportional to number of Photons in the rising endge)
    //getting wf baseline subtracted and wf baseline subtracted and normalized for the max amp.
    for(int i=0; i<samples; i++){max0=maxZ-base;
    Base_wf_v[i]=Raw_wf_v[i]-base; Norm_wf_v[i]=Base_wf_v[i]/max0;}
    //fitting the normalized baseline subtracted wf
    TGraph *gr = new TGraph(samples,x_wf_v,Norm_wf_v);
    TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power(([0]-x)/[1],4))",tick-10, tick);
    fit->SetParameters(tick,2,1);  gr->Fit("fit","Q","",tick-10, tick);
//std::cout<<std::endl;
//std::cout<<"tick "<<tick<<std::endl;
//std::cout<<std::endl;
//for(int i=0; i<3; i++){std::cout<<fit->GetParameter(i)<<", ";}
//std::cout<<std::endl;
//for(int i=0; i<samples; i++){std::cout<<gr->GetPointY(i)<<", ";}
//std::cout<<std::endl;
//std::cout<<std::endl;

    tca=fit->GetParameter(0);  tcb=fit->GetParameter(1);  tcc=fit->GetParameter(2);
    //timing is the risign edge half height
    TT[q]=(tca-abs(tcb*TMath::Power(-log(0.5/tcc),0.25)))/0.064; max[q]=max0;
    //----------------------------------------------
    //check for saturated wf
    if(maxZ<=saturation){TT[q]=TT[q]; max[q]=max[q]; Sat[q]=false;}
      else if(maxZ>saturation) { Sat[q]=true;
std::cout<<"Saturated PMT ch "<<q<<std::endl;
        //counting the number of ticks above the saturation, extended for NuMI
        for(int i=3*64; i<samples_64*64; i++){
//std::cout<<"Saturated PMT ch "<<q<<"  "<<Raw_wf_v[i]<<std::endl;
        if(TF==0){if(Raw_wf_v[i+1]>4094 && Raw_wf_v[i]<=4094){tickF=i; TF=1;}}
        if(TB==0){if(Raw_wf_v[i]>4094 && Raw_wf_v[i+1]<=4094){tickB=i; TB=1;}}}
        FB=tickB-tickF;  if(FB>99){FB=99;}
 	//amplitude discrete correction
        maxZhelp1=maxZ/Frac[FB]; tick=tickF; Nss=0; is=0;
        for(int i=3*64; i<samples_64*64; i++){if(Raw_wf_v[i]<4095){Nss=Nss+1;}}
        //double txSS[256],tySS[256],txSS2[256],tySS2[256];
	double txSS[1500],tySS[1500],txSS2[1500],tySS2[1500];
        for(int i=3*64; i<samples_64*64; i++){if(Raw_wf_v[i]<4095){txSS[is]=i*1.0; tySS[is]=Raw_wf_v[i]/maxZhelp1; is=is+1;}}
std::cout<<std::endl;
std::cout<<std::endl;
        for(int i=3*64; i<samples_64*64; i++){std::cout<<Raw_wf_v[i]<<", ";}
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<"FB "<<FB<<" tickB "<<tickB<<" tickF "<<tickF<<" Nss "<<Nss<<std::endl;
std::cout<<"Fit 1 from "<<tick-30<<" to "<<tick+250<<std::endl;
	TGraph *g1 = new TGraph(Nss,txSS,tySS);
        TF1 *fitS1 = new TF1("fitS1","[9]*(exp(-TMath::Power(([0]-(x-[8]))/[1],4))*0.5*(TMath::Erf(-(x-[8])-[7])+1.0)+([5]+[4]*exp(-TMath::Power(([2]-(x-[8]))/[3],2)))*exp((-(x-[8]))/[6])*0.5*(TMath::Erf([7]+(x-[8]))+1.0))",tick-30, tick+250);
	fitS1->SetParameters(pLL[0],pLL[1],pLL[2],pLL[3],pLL[4],pLL[5],pLL[6],pLL[7],tick,1.);
        for(int i=0; i<8; i++){fitS1->FixParameter(i,pLL[i]);} g1->Fit("fitS1","Q","",tick-30, tick+250);
        tickFit1=fitS1->GetParameter(8); maxZhelp2=fitS1->GetParameter(9);  maxZhelp3=maxZhelp1/maxZhelp2;
        //amplitude fit correction
        for(int i=0; i<Nss; i++){txSS2[i]=txSS[i]; tySS2[i]=tySS[i]/maxZhelp2;}
std::cout<<"Fit 2 from "<<tick-30<<" to "<<tick+250<<std::endl;
	TGraph *g2 = new TGraph(Nss,txSS2,tySS2);
        TF1 *fitS2 = new TF1("fitS2","exp(-TMath::Power(([0]-(x-[8]))/[1],4))*0.5*(TMath::Erf(-(x-[8])-[7])+1.0)+([5]+[4]*exp(-TMath::Power(([2]-(x-[8]))/[3],2)))*exp((-(x-[8]))/[6])*0.5*(TMath::Erf([7]+(x-[8]))+1.0)",tick-30, tick+250);
        fitS2->SetParameters(pLL[0],pLL[1],pLL[2],pLL[3],pLL[4],pLL[5],pLL[6],pLL[7],tickFit1);
        for(int i=0; i<8; i++){fitS2->FixParameter(i,pLL[i]);}
        g2->Fit("fitS2","Q","",tick-30, tick+250);  tickFit2=fitS2->GetParameter(8);
std::cout<<std::endl;
std::cout<<std::endl;
for(int i=0; i<10; i++){std::cout<<fitS1->GetParameter(i)<<", ";}
std::cout<<std::endl;
for(int i=tick-30; i<tick+250; i++){std::cout<<g1->GetPointY(i)<<", ";}
std::cout<<std::endl;
std::cout<<std::endl;

std::cout<<std::endl;
std::cout<<std::endl;
for(int i=0; i<9; i++){std::cout<<fitS2->GetParameter(i)<<", ";}
std::cout<<std::endl;
for(int i=tick-30; i<tick+250; i++){std::cout<<g2->GetPointY(i)<<", ";}
std::cout<<std::endl;
std::cout<<std::endl;

        TT[q]=tickFit2/0.064; max[q]=maxZhelp3;}
    //-------------------------------------------------------------------------------------------------------
    H_time->Fill(TT[q]);
    }
    //-------------------------------------------------------------------------------------------------------
    for(int q=0; q<32; q++){maxP[q]=max[q]; timeP[q]=TT[q];} //only two variables needed
}
double WireCellAnaTree::getBeamWF(art::Event const& e)
{
    //-------------------------------------------------------------------------------------------------------
    //get RWM--------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle_beam;
    e.getByLabel( "pmtreadout:UnspecifiedLogic", wf_handle_beam );
    if(!wf_handle_beam.isValid()) {
      std::cout<<"BeamTiming: No pmtreadout:UnspecifiedLogic!"<<std::endl;
      return -99999.0;
    }
    std::vector<double> *wf_w_03 = new std::vector<double>;
    for (size_t i=0; i < wf_handle_beam->size(); i++) {
      auto const wf = wf_handle_beam->at(i);
      auto ch = wf.ChannelNumber();
      std::cout<<"Channel: "<<ch<<"  size: "<<wf.size()<<std::endl;
      // BNB RWM is ch 39, NuMI RWM is ch 37
      if ( (ch != 39 && !fIsNuMI) || (ch != 37 && fIsNuMI) ) continue;
      for (size_t n=0; n < wf.size(); n++){
        if (n < 1500){
          wf_w_03->emplace_back((wf_handle_beam->at(i))[n]);
        }
      }// for all channels
    }// for all waveforms
    //=======================================================================================================
    //=======================================================================================================
    double beamBase,BBmax,wx[500],wy[500],pca,pcb,pcc, TT = -99999.0;
    int Btick,tickMax;
    //=======================================================================================================
    //-------baseline calculation-------------------------------------------
    beamBase=0.; BBmax=0.; Btick=0; TT=-9999.;
    for(int i=0; i<4*64; i++){beamBase=beamBase+wf_w_03->at(i);}
    beamBase=beamBase/(4.0*64.0);
    //baseline subtraction
    for(int i=0; i<500; i++){wx[i]=i*1.0; wy[i]=wf_w_03->at(i)-beamBase;
    //max amplitude
    if(BBmax<wy[i]){BBmax=wy[i]; Btick=i;}}
    H_maxH->Fill(BBmax);
    double BBmax_threshold_l0 = 2000;
    if(fIsNuMI) BBmax_threshold_l0 = 1800;//RWM pulse hight appears lower relative to baseline for NuMI, could possible dial this in more.
    if(BBmax>BBmax_threshold_l0 && BBmax<2100){
    //wf normalization
    for(int i=0; i<500; i++){wy[i]=wy[i]/BBmax;}
    //wf max check
    tickMax=0;
    for(int i=Btick-20; i<Btick+10; i++){if(wy[i-1]<1 && wy[i]==1){tickMax=i;}}
    //wf fit
    TGraph *gr0 = new TGraph(500,wx,wy);
    TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power((x-[0])/[1],4))",tickMax-6, tickMax);
    fit->SetParameters(tickMax,2,1);     gr0->Fit("fit","Q","",tickMax-6, tickMax);
    pca=fit->GetParameter(0); pcb=fit->GetParameter(1); pcc=fit->GetParameter(2);
    //timing is the risign edge half height
    TT=(pca-abs(pcb*TMath::Power(-log(0.5/pcc),0.25)))/0.064;
    H_t0_Beam->Fill(TT);
    }
    std::cout<<"Beam T0: "<<TT<<std::endl;
    return TT;
}


std::tuple< std::vector<float>*,std::vector<float>*,std::vector<float>*,std::vector<float>* > WireCellAnaTree::get_extrapolated_times(art::Ptr<simb::MCParticle> particle, double mother_time){
          //<x,y,z,t>
          std::vector<float> *_x = new std::vector<float>;
          std::vector<float> *_y = new std::vector<float>;
          std::vector<float> *_z = new std::vector<float>;
          std::vector<float> *_t = new std::vector<float>;

          double dx = 0.5;
          int pdg = particle->PdgCode();
          double start_x_pos = particle->Position().X();
          double start_y_pos = particle->Position().Y();
          double start_z_pos = particle->Position().Z();
          double end_x_pos = particle->EndPosition().X();
          double end_y_pos = particle->EndPosition().Y();
          double end_z_pos = particle->EndPosition().Z();
          double length = sqrt( pow(start_x_pos-end_x_pos,2) + pow(start_y_pos-end_y_pos,2) + pow(start_z_pos-end_z_pos,2) );
          // just assign the start and end points of neutrons
          if(pdg==2112){
             _x->push_back(start_x_pos);
             _y->push_back(start_y_pos);
             _z->push_back(start_z_pos);
             _t->push_back(mother_time);   
             _x->push_back(end_x_pos);
             _y->push_back(end_y_pos);
             _z->push_back(end_z_pos);
             _t->push_back(mother_time+length*fsol);        
             return std::make_tuple(_x,_y,_z,_t);
          }
          if(!isfinite(length)){length=0.3;}
          double residual_range = length;
          double x_pos = start_x_pos;
          double y_pos = start_y_pos;
          double z_pos = start_z_pos;
          double t_pos = mother_time;
          double mass = 0;
          if (pdg == 13){ mass = 0.1057;}
          if (pdg == 2212){ mass = 0.9397933;}//Don't do neutrons, KE is not assigned well so just assume c
          if (pdg == 211){ mass = 0.13982067;}
          double KE =  particle->Momentum()[3] - mass;
          double v = fsol*1/sqrt( 1-pow(mass/(mass+KE),2) );
          if(!isfinite(v)){v=fsol;}
	  if (pdg==22 || pdg==2112 || pdg==11){ v = fsol;}
          double gamma=0;//used for extrapolating the line from start to end position
          while (residual_range>=0){
            _x->push_back(x_pos);
            _y->push_back(y_pos);
            _z->push_back(z_pos);
            _t->push_back(t_pos);
            t_pos += v*dx;
            double dedx = get_dE_dx_range(residual_range,pdg)/1000;
            KE = KE-dedx*dx;
            v = fsol*1/sqrt( 1-pow(mass/(mass+KE),2) );
	    if (!isfinite(v)){v=fsol;}
            if (pdg==22 || pdg==2112 || pdg==11){ v = fsol;}
            gamma+=(dx/length);
            x_pos = start_x_pos + gamma*(end_x_pos-start_x_pos);
            y_pos = start_y_pos + gamma*(end_y_pos-start_y_pos);
            z_pos = start_z_pos + gamma*(end_z_pos-start_z_pos);
            residual_range = length - sqrt( pow(start_x_pos-x_pos,2) + pow(start_y_pos-y_pos,2) + pow(start_z_pos-z_pos,2));
            if (!isfinite(residual_range)){break;}
	  }

          return std::make_tuple(_x,_y,_z,_t);
}

double WireCellAnaTree::get_dE_dx_range(double R, int pdg){
    if (pdg==22 || pdg==11 || pdg==2112){ return 0; }
    double A = 8; 
    double b = -0.37; 
    if (pdg==2212){    
        A = 17;
        b = -0.42;
    }
    double dedx = A*pow(R,b);
    return dedx;
}

double WireCellAnaTree::TimeOffset(){
    // calculate in small to large

    // pick a time within a bucket
    double offset = fRndmGen->Gaus(0.0,fBucketTimeSigma);

    // pick a bucket within a batch
    // assume all ~ buckets constant in batch until we have another model
    offset +=  fTimeBetweenBuckets * (double)fRndmGen->Integer(fNFilledBucketsPerBatch);

    // pick a bucket
    bool   disallowed = true;
    size_t ibatch = 0;
    size_t nbatch = fCummulativeBatchPDF.size();
    double r = 2;
    while ( disallowed ) {
      r = fRndmGen->Uniform();
      for (ibatch=0; ibatch<nbatch; ++ibatch) {
        if ( r <= fCummulativeBatchPDF[ibatch] ) break;
      }
      disallowed = ( fDisallowedBatchMask[ibatch] != 0 );
    }
    offset += fTimeBetweenBuckets*(double)fNBucketsPerBatch*(double)ibatch;

    // finally the global offset
    return offset + fGlobalOffset;
}

void WireCellAnaTree::CalculateCPDF(std::vector<double> bi){
    fCummulativeBatchPDF.clear();
    double sum = 0;
    size_t nbi = bi.size();
    for (size_t i=0; i < nbi; ++i) {
      sum += bi[i];
      fCummulativeBatchPDF.push_back(sum);
    }
    // normalize to unit probability
    for (size_t i=0; i < nbi; ++i) fCummulativeBatchPDF[i] /= sum;
    // make sure the mask vector keeps up (but never make it smaller)
    // allowing all new batches
    if ( nbi > fDisallowedBatchMask.size() )
      fDisallowedBatchMask.resize(nbi,0);
}

DEFINE_ART_MODULE(WireCellAnaTree)

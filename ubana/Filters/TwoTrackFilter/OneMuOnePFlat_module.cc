////////////////////////////////////////////////////////////////////////
//// Class:       OneMuOneP
//// Plugin Type: analyzer (art v3_01_02)
//// File:        OneMuOneP_module.cc
////
//// Generated at Thu May  23 2019 by Yifan Chen using cetskelgen
//// from cetlib version v3_05_01.
//////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larsim/EventWeight/Base/WeightManager.h"

#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "TTree.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"

#include <tgmath.h>

#include "FiducialVolume.h"
#include "BackTrackerTruthMatch.h"
#include "RecoTruthMCParticle.h"
#include "Topology.h"
#include "PID.h"

class SingleMuon;


class SingleMuon : public art::EDAnalyzer {
public:
  explicit SingleMuon(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  //
  // Plugins should not be copied or assigned.
  SingleMuon(SingleMuon const&) = delete;
  SingleMuon(SingleMuon&&) = delete;
  SingleMuon& operator=(SingleMuon const&) = delete;
  SingleMuon& operator=(SingleMuon&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void endSubRun(art::SubRun const &sr) override;
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  ::ubana::FiducialVolume _fiducial_volume;

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<art::TFileService> tfs;

  spacecharge::SpaceCharge const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  detinfo::DetectorProperties const* detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();

  TTree * POTtree;
  int run, subrun;
  double POT_miss1E10;

  TTree * my_event_;

  void Initialize_event();

  int evt_run;
  int evt_subrun;
  int evt_evt;

  // -- MC Truth
  // GENIE reweight
  std::vector<double> MaCCQE;
  std::vector<double> CoulombCCQE;
  std::vector<double> MaNCEL;
  std::vector<double> EtaNCEL;

  std::vector<double> NormCCMEC;
  std::vector<double> NormNCMEC;
  std::vector<double> FracPN_CCMEC;
  std::vector<double> FracDelta_CCMEC;

  std::vector<double> MaCCRES;
  std::vector<double> MvCCRES;
  std::vector<double> MaNCRES;
  std::vector<double> MvNCRES;

  std::vector<double> NonRESBGvpCC1pi;
  std::vector<double> NonRESBGvpCC2pi;
  std::vector<double> NonRESBGvpNC1pi;
  std::vector<double> NonRESBGvpNC2pi;
  std::vector<double> NonRESBGvnCC1pi;
  std::vector<double> NonRESBGvnCC2pi;
  std::vector<double> NonRESBGvnNC1pi;
  std::vector<double> NonRESBGvnNC2pi;
  std::vector<double> NonRESBGvbarpCC1pi;
  std::vector<double> NonRESBGvbarpCC2pi;
  std::vector<double> NonRESBGvbarpNC1pi;
  std::vector<double> NonRESBGvbarpNC2pi;
  std::vector<double> NonRESBGvbarnCC1pi;
  std::vector<double> NonRESBGvbarnCC2pi;
  std::vector<double> NonRESBGvbarnNC1pi;
  std::vector<double> NonRESBGvbarnNC2pi;
  std::vector<double> AhtBY;
  std::vector<double> BhtBY;
  std::vector<double> CV1uBY;
  std::vector<double> CV2uBY;

  std::vector<double> AGKYxF1pi;
  std::vector<double> AGKYpT1pi;

  std::vector<double> MFP_pi;
  std::vector<double> MFP_N;
  std::vector<double> FrCEx_pi;
  std::vector<double> FrInel_pi;
  std::vector<double> FrAbs_pi;
  std::vector<double> FrCEx_N;
  std::vector<double> FrInel_N;
  std::vector<double> FrAbs_N;

  std::vector<double> RDecBR1gamma;
  std::vector<double> RDecBR1eta;

  std::vector<double> RootinoFix;
  std::vector<double> splines_general_Spline;

  std::vector<double> MaCCQE_tune;
  std::vector<double> RPA_CCQE_tune;
  std::vector<double> NormCCMEC_tune;
  std::vector<double> XSecShape_CCMEC_tune;

  //Geant4 weights
  std::vector<double> piplusReacLow;
  std::vector<double> piplusReacHigh;
  std::vector<double> piplusAbs;
  std::vector<double> piplusCex;
  std::vector<double> piplusDCex;
  std::vector<double> piplusPiProd;
  std::vector<double> piplusElast;

  std::vector<double> piminusReacLow;
  std::vector<double> piminusReacHigh;
  std::vector<double> piminusAbs;
  std::vector<double> piminusCex;
  std::vector<double> piminusDCex;
  std::vector<double> piminusPiProd;
  std::vector<double> piminusElast;

  std::vector<double> protonReac;
  std::vector<double> protonElast;

  std::vector<double> neutronAbs_max;
  std::vector<double> neutronElast_max;
  std::vector<double> neutronInela_max;
  std::vector<double> neutronPiprod_max;
  std::vector<double> neutronPprod_max;

  std::vector<double> neutronAbs_100;
  std::vector<double> neutronElast_100;
  std::vector<double> neutronInela_100;
  std::vector<double> neutronPiprod_100;
  std::vector<double> neutronPprod_100;

  std::vector<double> neutronAbs_200;
  std::vector<double> neutronElast_200;
  std::vector<double> neutronInela_200;
  std::vector<double> neutronPiprod_200;
  std::vector<double> neutronPprod_200;

  // flux parameters
  std::vector<double> expskin;
  std::vector<double> horncurrent;
  std::vector<double> nucleoninexsec;
  std::vector<double> nucleonqexsec;
  std::vector<double> nucleontotxsec;
  std::vector<double> pioninexsec;
  std::vector<double> pionqexsec;
  std::vector<double> piontotxsec;

  double EventWeight = 1; // Spine reweight using Steven's tool

  int Nr_MCNu = 0;

  int nu_parent_pdg = -999; //parent of the neutrino: mu-13 / pi-211 / k0-130 / k-321
  bool MC_beamNeutrino = false; // MCTruth beam origin
  bool MC_FV = false; // MCTruth vertex in FV = true, out of FV = false
  bool MC_if_in_active = false; // MCTruth vertex in active volume = true, out of active volume = false
  int MC_ccnc = -999; // MCTruth cc = 0 or nc = 1
  int MC_Q2 = -999; // MCTruth cc = 0 or nc = 1
  int MC_nupdg = -999; // MCTruth nupdg; numu = 14, nue = 12
  int MC_int_mode = -999; // https://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
  double MC_nu_E = -999; // MCTruth nu energy
  double MC_transfer_E = -999; // MCTruth nu energy
  double MC_nuVtxX = -999; // MCTruth nu vtx X
  double MC_nuVtxY = -999; // MCTruth nu vtx Y
  double MC_nuVtxZ = -999; // MCTruth nu vtx Z

  int MC_nMuon = 0; // Number of muon(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nElectron = 0; // Number of eletron(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nNeutron = 0; // Number of neutron(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_below260 = 0; // Number of proton(s) (p<260) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_above260 = 0; // Number of proton(s) (p >= 260) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPi0 = 0; // Number of pi0(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus_below80 = 0; // Number of pi plus(s) (p < 80MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus_above80 = 0; // Number of pi plus(s) (p > 80MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus_below80 = 0; // Number of pi minus(s) (p < 80MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus_above80 = 0; // Number of pi minus(s) (p > 80MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)

  std::vector<int> MC_Primary_PDG; // PDG of neutrino daughters
  std::vector<double> MC_Primary_Mom; // Momemtum of neutrino daughters

  std::vector<double> MC_muon_true_Mom;
  std::vector<double> MC_muon_true_theta;
  std::vector<double> MC_muon_true_cos_theta;
  std::vector<double> MC_muon_true_phi;
  std::vector<double> MC_muon_true_Px;
  std::vector<double> MC_muon_true_Py;
  std::vector<double> MC_muon_true_Pz;
  std::vector<double> MC_muon_true_E;

  std::vector<double> MC_proton_true_Mom;
  std::vector<double> MC_proton_true_theta;
  std::vector<double> MC_proton_true_cos_theta;
  std::vector<double> MC_proton_true_phi;
  std::vector<double> MC_proton_true_Px;
  std::vector<double> MC_proton_true_Py;
  std::vector<double> MC_proton_true_Pz;

  int TopologyType = -999;// The topology of true neutrino interaction + FSI products after Geant4

  double MC_0pi0p_muon_mom = -999;
  double MC_0pi0p_muon_theta = -999;
  double MC_0pi0p_muon_costheta = -999;
  double MC_0pi0p_muon_phi = -999;

  double MC_0pi1p_muon_mom = -999;
  double MC_0pi1p_muon_theta = -999;
  double MC_0pi1p_muon_costheta = -999;
  double MC_0pi1p_muon_phi = -999;

  double MC_0pi1p_proton_mom = -999;
  double MC_0pi1p_proton_theta = -999;
  double MC_0pi1p_proton_costheta = -999;
  double MC_0pi1p_proton_phi = -999;

  double MC_0pi1p_cos_ang_muon_proton = -999;
  double MC_0pi1p_PT = -999;
  double MC_0pi1p_PL = -999;

  std::vector<bool> if_MC_reco_match = {false, false};

  double flash_matching_chi2 = -999; //Chi2 of flash matching in each neutrino slice
  double trigger_time = -999;
  double flash_YCenter = -999;
  double flash_YWidth = -999;
  double flash_ZCenter = -999;
  double flash_ZWidth = -999;
  double flash_TotalPE = -999;
  double flash_time = -999;

  int n_pfp_nuDaughters;
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters

  std::vector<int> Ghost_PDG; // pdg code of the pfp which has no track or shower associated; No elements ideally

  bool if_2tracks = false; // if there are two primary tracks and nothing

  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime
  double only_crthit_time_inBeam = -999;
  bool if_t0_crt_time_inBeam_match = true;

  bool evt_CRTveto_70 = false; // If CRT veto, eliminate the events for contained (70PE threshold)
  bool evt_CRTveto_100 = false; // If CRT veto, eliminate the events for contained (100PE threshold)

  std::vector<int> Nr_trk_asso_crthit = {0, 0}; // Number of CRT hits associated to the track
  std::vector<double> trk_crt_time = {-999, -999};// The CRT time of the hit which matched to the track
  std::vector<bool> if_trk_CRT_out_Beam = {false, false}; // Check if a track matches with out of beam CRT hit(s)
  std::vector<bool> if_t0_trk_crt_time_match = {true, true};

  std::vector<bool> trk_OutOfTime = {false, false};

  std::vector<double> trk_start_x = {-999, -999};
  std::vector<double> trk_start_y = {-999, -999};
  std::vector<double> trk_start_z = {-999, -999};
  std::vector<double> trk_end_x = {-999, -999};
  std::vector<double> trk_end_y = {-999, -999};
  std::vector<double> trk_end_z = {-999, -999};

  std::vector<bool> trk_start_InFV = {false, false};
  std::vector<bool> trk_contained = {false, false};

  std::vector<double> trk_phi = {-999, -999};
  std::vector<double> trk_theta = {-999, -999};
  std::vector<double> trk_costheta = {-999, -999};

  std::vector<double> theta_pl2 = {-999, -999};
  std::vector<double> theta_pl1 = {-999, -999};
  std::vector<double> theta_pl0 = {-999, -999};
  std::vector<double> phi_readout = {-999, -999};

  std::vector<double> sin2_theta_pl2 = {-999, -999};
  std::vector<double> sin2_theta_pl1 = {-999, -999};
  std::vector<double> sin2_theta_pl0 = {-999, -999};
  std::vector<double> sin2_phi_readout = {-999, -999};

  std::vector<double> bestMCS = {-999, -999};
  std::vector<double> bestMCSLL = {-999, -999};
  std::vector<double> fwdMCS = {-999, -999};
  std::vector<double> fwdMCSLL = {-999, -999};
  std::vector<double> bwdMCS = {-999, -999};
  std::vector<double> bwdMCSLL = {-999, -999};

  std::vector<double> bestMCS_NoSCE = {-999, -999};
  std::vector<double> fwdMCS_NoSCE = {-999, -999};
  std::vector<double> bwdMCS_NoSCE = {-999, -999};

  std::vector<double> bestMCSLL_NoSCE = {-999, -999};
  std::vector<double> fwdMCSLL_NoSCE = {-999, -999};
  std::vector<double> bwdMCSLL_NoSCE = {-999, -999};

  std::vector<double> Trk_length = {-999, -999};
  std::vector<double> Mom_Range_mu = {-999, -999};
  std::vector<double> Mom_Range_p = {-999, -999};
  std::vector<double> Mom_Range_pi = {-999, -999};

  std::vector<float> dEdx_pl0; // dE/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<float> dEdx_pl1; // dE/dx of the selected (muon) track from plane 1
  std::vector<float> dEdx_pl2; // dE/dx of the selected (muon) track from plane 2 (collection)
  std::vector<float> resRange_pl0; // range from a hit to the end of the selected track end
  std::vector<float> resRange_pl1; // range from a hit to the end of the selected track end
  std::vector<float> resRange_pl2; // range from a hit to the end of the selected track end

  std::vector<double> dEdx_pl0_start_half = {-999, -999};
  std::vector<double> dEdx_pl1_start_half = {-999, -999};
  std::vector<double> dEdx_pl2_start_half = {-999, -999};
  std::vector<double> dEdx_pl0_end_half = {-999, -999};
  std::vector<double> dEdx_pl1_end_half = {-999, -999};
  std::vector<double> dEdx_pl2_end_half = {-999, -999};

  std::vector<double> dEdx_pl0_start1020 = {-999, -999};
  std::vector<double> dEdx_pl1_start1020 = {-999, -999};
  std::vector<double> dEdx_pl2_start1020 = {-999, -999};
  std::vector<double> dEdx_pl0_end1020 = {-999, -999};
  std::vector<double> dEdx_pl1_end1020 = {-999, -999};
  std::vector<double> dEdx_pl2_end1020 = {-999, -999};

  std::vector<double> dEdx_pl2_1020_ratio = {-999, -999};
  std::vector<double> dEdx_pl2_half_ratio = {-999, -999};

  std::vector<double> PID_Chi2Mu_3pl = {-999, -999};
  std::vector<double> PID_Chi2P_3pl = {-999, -999};
  std::vector<double> PID_Chi2Pi_3pl = {-999, -999};
  std::vector<double> PID_Chi2K_3pl = {-999, -999};

  std::vector<double> trk_cosmic_percent = {-999, -999};
  std::vector<double> trk_purity = {-999, -999};
  std::vector<double> trk_completeness = {-999, -999};

  std::vector<bool> if_cosmic = {true, true};
  std::vector<bool> if_matchMu = {false, false};
  std::vector<bool> if_matchP = {false, false};
  std::vector<bool> if_matchPrimary = {false, false};

  std::vector<double> true_mom = {-999, -999};
  std::vector<double> true_start_x = {-999, -999};
  std::vector<double> true_start_y = {-999, -999};
  std::vector<double> true_start_z = {-999, -999};
  std::vector<double> true_end_x = {-999, -999};
  std::vector<double> true_end_y = {-999, -999};
  std::vector<double> true_end_z = {-999, -999};
  std::vector<double> true_trk_phi = {-999, -999};
  std::vector<double> true_trk_theta = {-999, -999};
  std::vector<double> true_trk_costheta = {-999, -999};
  std::vector<double> true_trk_theta_yz = {-999, -999};
  std::vector<double> true_trk_costheta_yz = {-999, -999};
  std::vector<double> true_trk_theta_xz = {-999, -999};
  std::vector<double> true_trk_costheta_xz = {-999, -999};
  std::vector<double> true_trk_length = {-999, -999};
  std::vector<double> true_trk_PDG = {-999, -999};
  std::vector<bool> true_trk_ifcontained = {false, false};
  std::vector<bool> true_vtxFV = {false, false};

  int muon_idx = -999;
  int proton_idx = -999;

  double vtx_x = -999;
  double vtx_y = -999;
  double vtx_z = -999;
  bool vtx_InFV = false;

  double muon_cand_phi = -999;
  double muon_cand_theta = -999;
  double muon_cand_costheta = -999;
  double muon_cand_mom_MCS = -999;
  double muon_cand_mom_range = -999;
  double muon_cand_length = -999;

  double proton_cand_phi = -999;
  double proton_cand_theta = -999;
  double proton_cand_costheta = -999;
  double proton_cand_mom_range = -999;
  double proton_cand_length = -999;

  double cos_ang_muon_proton = -999;

  double dist_mu_p_start = -999;

  double PT = -999;
  double PL = -999;

  bool if_match_mu_p = false;
  bool if_match_mu_p_flipped = false;

  double reco_MC_dist_vtx = -999; // Distance of reco - MC vertex w/ SCE correction
  double reco_MC_dist_vtx_noSCE = -999; // Distance of reco - MC vertex w/o SCE correction

  int MCmatch_type = 0;
  unsigned int matched_MCid = 0;
  std::vector<int> matched_MCid_2trk = {-999, -999};

  bool                                IsMC;
  bool                                T0Corr;
  bool                                UsingCRT;
  bool                                if_GENIE_ReWeight;
  bool                                if_G4_ReWeight;
  bool                                if_Flux_ReWeight;
  double                              fBeamStart;
  double                              fBeamEnd;
  double                              fDTOffset;
  double                              fDTOffset_overlay;
  std::string                         m_DAQHeaderProducer;
  std::string                         m_generatorLabel;
  std::string                         m_geantLabel;
  std::string                         m_GENIEeventweightLabel;
  std::string                         m_G4eventweightLabel;
  std::string                         m_FluxeventweightLabel;
  std::string                         m_pandoraLabel;
  std::string                         m_T0ProducerLabel;
  std::string                         m_hitProducerLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_showerProducerLabel;
  std::string                         m_MCSmuProducerLabel;
  std::string                         m_MCSmuNoSCEProducerLabel;
  std::string                         m_calorimetryProducerLabel;
  std::string                         Hits_TrackAssLabel;
  std::string                         PID_TrackAssLabel;
  std::string                         CRT_TrackAssLabel;
  std::string                         m_CRTVetoLabel;
  std::string                         m_CRTHitLabel;
  std::string                         m_CRTcorrT0Label;
  std::string                         m_FlashLabel;

  std::vector<std::string>            genie_pars;
  std::vector<std::string>            flux_pars;
  std::vector<std::string>            g4_pars;

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};

SingleMuon::SingleMuon(fhicl::ParameterSet const& pset)
  :
  EDAnalyzer{pset},
  IsMC(pset.get<bool>("IsMC")),
  T0Corr(pset.get<bool>("T0Corr")),
  UsingCRT(pset.get<bool>("UsingCRT")),
  if_GENIE_ReWeight(pset.get<bool>("if_GENIE_ReWeight")),
  if_G4_ReWeight(pset.get<bool>("if_G4_ReWeight")),
  if_Flux_ReWeight(pset.get<bool>("if_Flux_ReWeight")),
  fBeamStart(pset.get<double>("BeamStart")),
  fBeamEnd(pset.get<double>("BeamEnd")),
  fDTOffset(pset.get<double>("DTOffset")),
  fDTOffset_overlay(pset.get<double>("DTOffset_overlay")),
  m_DAQHeaderProducer(pset.get<std::string>("DAQHeaderProducer")),
  m_generatorLabel(pset.get<std::string>("GeneratorLabel")),
  m_geantLabel(pset.get<std::string>("GeantLabel")),
  m_GENIEeventweightLabel(pset.get<std::string>("GENIE_EventweightLabel")),
  m_G4eventweightLabel(pset.get<std::string>("G4_EventweightLabel")),
  m_FluxeventweightLabel(pset.get<std::string>("Flux_EventweightLabel")),
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_T0ProducerLabel(pset.get<std::string>("T0ProducerLabel")),
  m_hitProducerLabel(pset.get<std::string>("HitProducerLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel")),
  m_MCSmuProducerLabel(pset.get<std::string>("MCSmuProducerLabel")),
  m_MCSmuNoSCEProducerLabel(pset.get<std::string>("MCSmuNoSCEProducerLabel")),
  m_calorimetryProducerLabel(pset.get<std::string>("calorimetryProducerLabel")),
  m_CRTVetoLabel(pset.get<std::string>("CRTVetoLabel")),
  m_CRTHitLabel(pset.get<std::string>("CRTHitProducer")),
  m_CRTcorrT0Label(pset.get<std::string>("CRTcorrT0Label")),
  m_FlashLabel(pset.get<std::string>("FlashLabel")),
  _min_track_len{pset.get<double>("MinTrackLength", 0.1)},
  _trk_mom_calculator{_min_track_len}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  if (IsMC && if_GENIE_ReWeight){
    genie_pars = pset.get< std::vector<std::string> >( "genie_parameter_list" );
  }
  if (IsMC && if_G4_ReWeight){
    g4_pars = pset.get< std::vector<std::string> >( "g4_parameter_list" );
  }
  if (IsMC && if_Flux_ReWeight){
    flux_pars = pset.get< std::vector<std::string> >( "flux_parameter_list" );
  }

  _fiducial_volume.Configure(pset.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();
  Hits_TrackAssLabel = pset.get<std::string>("HitsPerTrackAssLabel");
  PID_TrackAssLabel = pset.get<std::string>("PIDTrackAssLabel");
  CRT_TrackAssLabel = pset.get<std::string>("CRTTrackAssLabel");
}

void SingleMuon::analyze(art::Event const& evt)
{

  // run, subrun, evt
  evt_run = evt.run();
  evt_subrun = evt.subRun();
  evt_evt = evt.event();

  std::cout<<"-------------------------"<<std::endl;
  std::cout<<"Run: "<< evt_run <<" | Subrun: " << evt_subrun <<" | Event: " << evt_evt<<std::endl;
  std::cout<<"-------------------------"<<std::endl;

  // prepare X offset for position correction (SCE)
  double xtimeoffset = 0;
  if(IsMC){
    auto const& mct_h = evt.getValidHandle<std::vector<simb::MCTruth> >("generator");
    auto gen = mct_h->at(0);
    double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
    xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0) + 0.6; //6mm X=0 differeneces in wireplanes, wirecell and reconstruction
  }
  else{
    xtimeoffset = 0;
  }

  //// Get reco handles
  // Slice
  art::Handle<std::vector<recob::Slice> > Handle_Slice;
  evt.getByLabel(m_pandoraLabel, Handle_Slice);

  // Hit
  art::Handle<std::vector<recob::Hit> > Handle_Hit;
  evt.getByLabel(m_hitProducerLabel, Handle_Hit);

  //Track
  art::Handle<std::vector<recob::Track> > Handle_TPCtrack;
  evt.getByLabel(m_trackProducerLabel, Handle_TPCtrack);
  std::vector< art::Ptr<recob::Track> > AllTrackCollection;
  art::fill_ptr_vector(AllTrackCollection, Handle_TPCtrack);

  //Shower
  art::Handle<std::vector<recob::Shower> > Handle_TPCshower;
  evt.getByLabel(m_showerProducerLabel, Handle_TPCshower);
  std::vector< art::Ptr<recob::Shower> > AllShowerCollection;
  art::fill_ptr_vector(AllShowerCollection, Handle_TPCshower);

  //MCS momentum (the current fcl parameter has spatial correction)
  art::Handle<std::vector<recob::MCSFitResult> > Handle_mcsfitresult_mu;
  evt.getByLabel(m_MCSmuProducerLabel, Handle_mcsfitresult_mu);
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_v;
  art::fill_ptr_vector(mcsfitresult_mu_v, Handle_mcsfitresult_mu);

  //MCS momentum (without SCE)
  art::Handle<std::vector<recob::MCSFitResult> > Handle_mcsfitresult_mu_NoSCE;
  evt.getByLabel(m_MCSmuNoSCEProducerLabel, Handle_mcsfitresult_mu_NoSCE);
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_NoSCE_v;
  art::fill_ptr_vector(mcsfitresult_mu_NoSCE_v, Handle_mcsfitresult_mu_NoSCE);

  //Hits per track
  art::FindManyP<recob::Hit> hits_per_track(Handle_TPCtrack,evt,Hits_TrackAssLabel);

  //PFParticle collection
  art::Handle< std::vector<recob::PFParticle> > Handle_pfParticle;
  evt.getByLabel(m_pandoraLabel, Handle_pfParticle);
  std::vector<art::Ptr<recob::PFParticle>> pfParticle_v;
  art::fill_ptr_vector(pfParticle_v, Handle_pfParticle);

  // Flash Collection
  art::Handle< std::vector<recob::OpFlash> > Handle_opflash;
  evt.getByLabel(m_FlashLabel, Handle_opflash);
  std::vector<art::Ptr<recob::OpFlash>> flash_v;
  art::fill_ptr_vector(flash_v, Handle_opflash);

  // DAQ
  art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
  evt.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);

  // CRT
  art::Handle<std::vector<crt::CRTHit>> Handle_crthit;
  evt.getByLabel(m_CRTHitLabel, Handle_crthit);
  std::vector<art::Ptr<crt::CRTHit>> crthit_v;
  art::fill_ptr_vector(crthit_v, Handle_crthit);
  Nr_crthit_inBeam = 0;

  // CRT corr T0
  std::vector<art::Ptr<anab::T0>> crtT0_v;
  if(T0Corr){
    art::Handle<std::vector<anab::T0>> Handle_crtcorrT0;
    evt.getByLabel(m_CRTcorrT0Label, Handle_crtcorrT0);
    art::fill_ptr_vector(crtT0_v, Handle_crtcorrT0);
  }

  //// Get necessary MC handles
  std::vector<art::Ptr<simb::MCFlux>> MCFluxCollection;
  art::Handle< std::vector<simb::MCFlux>> Handle_MCFlux;

  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;

  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
  art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;

  std::vector<art::Ptr<evwgh::MCEventWeight> > GENIEWeightCollection;
  art::Handle< std::vector<evwgh::MCEventWeight> > Handle_GENIEWeight;

  std::vector<art::Ptr<evwgh::MCEventWeight> > G4WeightCollection;
  art::Handle< std::vector<evwgh::MCEventWeight> > Handle_G4Weight;

  std::vector<art::Ptr<evwgh::MCEventWeight> > FluxWeightCollection;
  art::Handle< std::vector<evwgh::MCEventWeight> > Handle_FluxWeight;

  // Neutrino Origin
  const simb::Origin_t Neutrino_Origin = simb::kBeamNeutrino;

  if(IsMC){
    // MC Flux
    evt.getByLabel(m_generatorLabel,Handle_MCFlux);
    art::fill_ptr_vector(MCFluxCollection, Handle_MCFlux);

    // MC Truth
    evt.getByLabel(m_generatorLabel, Handle_MCTruth);
    art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

    // MC Particle
    evt.getByLabel(m_geantLabel, Handle_MCParticle);
    art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

    // Genie Reweight
    if(if_GENIE_ReWeight){
      evt.getByLabel(m_GENIEeventweightLabel, Handle_GENIEWeight);
      art::fill_ptr_vector(GENIEWeightCollection, Handle_GENIEWeight);
    }

    // Flux Reweight
    if(if_Flux_ReWeight){
      evt.getByLabel(m_FluxeventweightLabel, Handle_FluxWeight);
      art::fill_ptr_vector(FluxWeightCollection, Handle_FluxWeight);
    }

    // Geant Reweight
    if(if_G4_ReWeight){
      evt.getByLabel(m_G4eventweightLabel, Handle_G4Weight);
      art::fill_ptr_vector(G4WeightCollection, Handle_G4Weight);
    }
  }

  //// Get association
  // Vertex - PFP association
  art::FindMany<recob::Vertex> pfpToVtxAsso(Handle_pfParticle, evt, m_pandoraLabel);

  // CRT Hit - Flash association
  art::FindMany<crt::CRTHit> CRThitFlashAsso(Handle_opflash, evt, m_CRTVetoLabel);

  // pfp t0 association (to get flash_matching chi2 score)
  art::FindMany<anab::T0> pfpToT0Asso(Handle_pfParticle, evt, m_T0ProducerLabel);

  //pfp track association
  art::FindManyP<recob::Track> pfpToTrackAsso(Handle_pfParticle, evt, m_trackProducerLabel);

  //pfp shower association
  art::FindManyP<recob::Shower> pfpToShowerAsso(Handle_pfParticle, evt, m_showerProducerLabel);

  //track calorimetry association
  art::FindManyP<anab::Calorimetry> trackToCalAsso(Handle_TPCtrack, evt, m_calorimetryProducerLabel);

  //PID
  art::FindManyP<anab::ParticleID> PIDTotrackAsso(Handle_TPCtrack,evt,PID_TrackAssLabel);

  // CRT Hit Track association
  art::FindManyP<crt::CRTHit> CRTToTrackAsso(Handle_TPCtrack,evt,CRT_TrackAssLabel);

  // Slice PFP association
  art::FindManyP<recob::Slice> SliceToPFPAsso(Handle_pfParticle, evt, m_pandoraLabel);

  // Hit Slice association
  art::FindManyP<recob::Hit> HitToSliceAsso(Handle_Slice, evt, m_pandoraLabel);

  // SpacePoint Hit association
  art::FindManyP<recob::SpacePoint> SpacePointToHitAsso(Handle_Hit, evt, m_pandoraLabel);

  // Trajectory Point Hit association
  art::FindManyP<recob::Hit, recob::TrackHitMeta> TrackHitToHitAsso(Handle_TPCtrack, evt, m_trackProducerLabel); //this has more information about hit-track association, only available in PMA for now

  // MC truth to MC particle association
  //art::FindMany<simb::MCParticle,sim::GeneratedParticleInfo> MCtToMCpAsso(Handle_MCTruth, evt, m_geantLabel);

  //////////////////////////////////////////////////////////////////////////////
  // If it is data, the MC info would be void, and the reco info would be the 1st(also the only) entryof the event
  // If is is MC, each neutrino interactions would get an entry
  // .. If the reco info matches with an MC interaction, then the reco info would be filled under the same entry as that MC info
  // .. The other potential unmatched MC interaction(s) would contain void reco info
  // .. Otherwise if the reco info does not match any MC interaction (say the selected track is a cosmic track), then the reco info would be registered under the 1st entry
  // .. And the other potential MC interaction(s) would contain void reco info
  //////////////////////////////////////////////////////////////////////////////


  //-------Reco

  // Get mapping from ID to PFParticle
  std::unordered_map<size_t, art::Ptr<recob::PFParticle> > pfParticleIdMap;
  for (unsigned int i = 0; i < Handle_pfParticle->size(); ++i){
    auto pfp = pfParticle_v[i];
    if (!pfParticleIdMap.emplace(pfp->Self(), pfp).second){
      throw cet::exception("[Numu0pi0p]")<< "Repeated PFParticles!" << std::endl;
    }
  }

  //Define the daughters of neutrino
  std::vector<art::Ptr<recob::PFParticle> > NeutrinoDaughters;
  std::vector<art::Ptr<recob::Track> > daughter_Tracks;
  std::vector<art::Ptr<recob::Shower> > daughter_Showers;
  std::vector<int> Track_PDG; // The oder follows daughter_Tracks
  std::vector<int> Shower_PDG; // The oder follows daughter_Showers

  std::vector<float> vTrk_len;

  TVector3 true_nuVtx; // useful to calculate the distance to the true vtx
  std::vector<TVector3> Trk_start(2);
  std::vector<TVector3> Trk_end(2);
  std::vector<TVector3> Trk_start_SCEcorr(2);
  std::vector<TVector3> Trk_end_SCEcorr(2);
  std::vector<TVector3> trk_startdir(2);

  std::vector<art::Ptr< simb::MCParticle >> selected_MCparticle(2);

  //-------- Get Reco neutrino (pfparticle)
  for(unsigned int i = 0; i < pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
    // Check if the pfp is leading a neutrino slice
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      n_pfp_nuDaughters = pfp->NumDaughters();

      // For CC0pi1p, we only consider the case with the number of neutrino daughters less than 6
      if(n_pfp_nuDaughters < 6){
        // Get the pointer for the daughters of the neutrino
        for(int j = 0; j < n_pfp_nuDaughters; j++){
          auto Iterator = pfParticleIdMap.find(pfp->Daughters().at(j));
          auto dau_pfp = Iterator->second;
          NeutrinoDaughters.push_back(dau_pfp);
          // Collect pfparticle associated track in a vector
          auto assoTrack = pfpToTrackAsso.at(dau_pfp.key()); // vector
          if(assoTrack.size()==1){
            daughter_Tracks.push_back(assoTrack.front());
            Track_PDG.push_back(dau_pfp->PdgCode());
          }
          if(assoTrack.size()>1){
            throw cet::exception("[Numu0pi0p]") << "PFParticle has >1 track!" << std::endl;
          }
          // Collect pfparticle associated shower in a vector
          auto assoShower = pfpToShowerAsso.at(dau_pfp.key()); // vector
          if(assoShower.size()==1){
            daughter_Showers.push_back(assoShower.front());
            Shower_PDG.push_back(dau_pfp->PdgCode());
          }
          if(assoShower.size()>1){
            throw cet::exception("[Numu0pi0p]") << "PFParticle has >1 shower!" << std::endl;
          }

          // If some pfparticles are built without associated tracks or showers
          if(assoTrack.empty() && assoShower.empty()){
            Ghost_PDG.push_back(dau_pfp->PdgCode());
          }
        } // finish looping of daughter pfps
      } // If the the nr of daughter pfps is < 6

      //number of tracks and showers
      n_dau_tracks = daughter_Tracks.size();
      n_dau_showers = daughter_Showers.size();

      // Selection and Fill in Info
      if(n_dau_tracks == 2 && n_dau_showers == 0){

        if_2tracks = true;

        //////////////////////
        // Information to be filled per slice: flash, number of CRT hits in beam
        //////////////////////
        
        // CRT related information
        double evt_timeGPS_nsec = 0.;
        if(!rawHandle_DAQHeader.isValid()) {
           throw cet::exception("[Numu0pi0p]") << "The DAQ header cannot be located." << std::endl;
         }
        raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
        art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
        evt_timeGPS_nsec = evtTimeGPS.timeLow();

        double CRTT0corr = 0;
        if(T0Corr){
          if(crtT0_v.size() == 1){
            CRTT0corr = crtT0_v.front()->Time();
          }
          else{
            CRTT0corr = 0;
          }
        }

        if(UsingCRT){
          // overlay is also basically data, using ts0
          for (unsigned int i_crt = 0; i_crt < crthit_v.size(); i_crt ++){
            // figure out what plane this hit comes from
            // 3 -> top, 0 -> bottom, 1 -> anode, 2 -> cathode
            double crt_time = -999;
            if(T0Corr){
              crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
            }
            else if(!T0Corr){
              crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
            }

            crthit_PE.push_back(crthit_v[i_crt]->peshit);
            crthit_plane.push_back(crthit_v[i_crt]->plane);
            crthit_time.push_back(crt_time);

            if(crt_time >= fBeamStart && crt_time <= fBeamEnd && crthit_v[i_crt]->peshit > 100){
              Nr_crthit_inBeam++; // The total number of CRT hits in beam should be less than the number of exiting tracks
              only_crthit_time_inBeam = crt_time;
            } // If CRT hit in Beam window

          } // CRT hit loop

          if(Nr_crthit_inBeam == 1){
            if(trigger_time != -999){
              if(abs(only_crthit_time_inBeam - trigger_time) > 1){
                if_t0_crt_time_inBeam_match = false;
              }
            }
          }

          //- For contained (Veto if there is any CRT hit within 1us of the flash which is in the beam window)
          if(flash_v.size() > 0){
            for(unsigned int i_fl = 0; i_fl < flash_v.size(); i_fl++){
              auto CRT_hit = CRThitFlashAsso.at(flash_v[i_fl].key());
              if(CRT_hit.size() > 0){
                if(CRT_hit.front()->peshit > 70) evt_CRTveto_70 = true;
                if(CRT_hit.front()->peshit > 100) evt_CRTveto_100 = true;
              } // if CRT veto
            } // loop over flash(es)
          } // if flash exists
        } // Using CRT

        //////////////////////
        // Information required per track (2 tracks)
        /////////////////////

        for(int i_trk = 0; i_trk < n_dau_tracks; i_trk++){
          //----- CRT info
          if(UsingCRT){
            // If a track is associated to a CRT hit which is outside of beam window, exclude them (For both contained and exiting)
            auto Track_CRThit = CRTToTrackAsso.at(daughter_Tracks[i_trk].key());
            // If a track is exiting, Nr_trk_asso_crthit <= 1
            // If a track is contained, Nr_trk_asso_crthit = 0
            // Else remove the slice
            if(Track_CRThit.size() > 0){
              for(unsigned int i_trk_crt = 0; i_trk_crt < Track_CRThit.size(); i_trk_crt++){
                if(Track_CRThit[i_trk_crt]->peshit > 100){

                  // w/ T0 correction
                  if(T0Corr){
                    trk_crt_time[i_trk] = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
                  }
                  // w/o T0 correction
                  else if(!T0Corr){
                    trk_crt_time[i_trk] = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
                  }

                  Nr_trk_asso_crthit[i_trk]++;
                }

                if(trk_crt_time[i_trk] < fBeamStart || trk_crt_time[i_trk] > fBeamEnd){
                  // For 2 tracks, as long as one of them has associated CRT hit out of time, reject!
                  if_trk_CRT_out_Beam[i_trk] = true;
                } // If matched CRT hit out of Beam window

              } // Loop over the associated CRT hits
            }  // If there is associated CRT hits

            // - crt time and t0 time match
            if(trigger_time != -999 && trk_crt_time[i_trk] != -999){
              if(abs(trk_crt_time[i_trk] - trigger_time) < 1){
                if_t0_trk_crt_time_match[i_trk] =  true;
              }
              else{
                if_t0_trk_crt_time_match[i_trk] = false;
              }
            }
          } // CRT info

          //------Info of the tracks

          //-- Track start and end (The track start of the muon candidate track will be the vertex)
          Trk_start[i_trk] = daughter_Tracks[i_trk]->Start<TVector3>();
          auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start[i_trk].X(), Trk_start[i_trk].Y(), Trk_start[i_trk].Z()));
          Trk_start_SCEcorr[i_trk].SetX(Trk_start[i_trk].X() - Trk_start_offset.X() + xtimeoffset);
          Trk_start_SCEcorr[i_trk].SetY(Trk_start[i_trk].Y() + Trk_start_offset.Y());
          Trk_start_SCEcorr[i_trk].SetZ(Trk_start[i_trk].Z() + Trk_start_offset.Z());

          trk_start_x[i_trk] = Trk_start_SCEcorr[i_trk].X();
          trk_start_y[i_trk] = Trk_start_SCEcorr[i_trk].Y();
          trk_start_z[i_trk] = Trk_start_SCEcorr[i_trk].Z();

          Trk_end[i_trk] = daughter_Tracks[i_trk]->End<TVector3>();
          auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end[i_trk].X(), Trk_end[i_trk].Y(), Trk_end[i_trk].Z()));
          Trk_end_SCEcorr[i_trk].SetX(Trk_end[i_trk].X() - Trk_end_offset.X() + xtimeoffset);
          Trk_end_SCEcorr[i_trk].SetY(Trk_end[i_trk].Y() + Trk_end_offset.Y());
          Trk_end_SCEcorr[i_trk].SetZ(Trk_end[i_trk].Z() + Trk_end_offset.Z());

          trk_end_x[i_trk] = Trk_end_SCEcorr[i_trk].X();
          trk_end_y[i_trk] = Trk_end_SCEcorr[i_trk].Y();
          trk_end_z[i_trk] = Trk_end_SCEcorr[i_trk].Z();

          //-- Track out of time
          if( Trk_start_SCEcorr[i_trk].X() < 0 || Trk_start_SCEcorr[i_trk].X() > 2. * geo->DetHalfWidth() || Trk_end_SCEcorr[i_trk].X() < 0 || Trk_end_SCEcorr[i_trk].X() > 2. * geo->DetHalfWidth()){
            trk_OutOfTime[i_trk] = true;
          }
          else{
            trk_OutOfTime[i_trk] = false;
          }

          //-- If contained
          trk_start_InFV[i_trk] = _fiducial_volume.VertexInFV(Trk_start_SCEcorr[i_trk]);
          trk_contained[i_trk] = _fiducial_volume.TrackContain(Trk_start_SCEcorr[i_trk], Trk_end_SCEcorr[i_trk]);

          //-- Track direction
          trk_startdir[i_trk] = daughter_Tracks[i_trk]->StartDirection<TVector3>();

          //-- Angle wrt to the detector coordinator
          trk_phi[i_trk] = daughter_Tracks[i_trk]->Phi();
          trk_theta[i_trk] = daughter_Tracks[i_trk]->Theta();
          trk_costheta[i_trk] = cos(daughter_Tracks[i_trk]->Theta());

          //-- Angle wrt to the readout; theta is to the wire; phi is to x
          auto Trk_pos = Trk_end_SCEcorr[i_trk] - Trk_start_SCEcorr[i_trk];
          theta_pl2[i_trk] = std::atan2(Trk_pos.Z(), Trk_pos.Y()); // atan2(y,x)
          theta_pl1[i_trk] = theta_pl2[i_trk] + M_PI/3; // If plan1 is -60 degree to Y, looking from outside to the TPC
          theta_pl0[i_trk] = theta_pl2[i_trk] - M_PI/3; // If plan0 is +60 degree to Y, looking from outside to the TPC
          phi_readout[i_trk] = std::atan2(Trk_pos.X(), Trk_pos.Y());

          sin2_theta_pl2[i_trk] = sin(theta_pl2[i_trk]) * sin(theta_pl2[i_trk]);
          sin2_theta_pl1[i_trk] = sin(theta_pl1[i_trk]) * sin(theta_pl1[i_trk]);
          sin2_theta_pl0[i_trk] = sin(theta_pl0[i_trk]) * sin(theta_pl0[i_trk]);
          sin2_phi_readout[i_trk] = sin(phi_readout[i_trk]) * sin(phi_readout[i_trk]);

          //-- MCS momentum
          bestMCS[i_trk] =  mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bestMomentum();
          bestMCSLL[i_trk] =  mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bestLogLikelihood();
          fwdMCS[i_trk] =  mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->fwdMomentum();
          fwdMCSLL[i_trk] =  mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->fwdLogLikelihood();
          bwdMCS[i_trk] =  mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bwdMomentum();
          bwdMCSLL[i_trk] =  mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bwdLogLikelihood();

          bestMCS_NoSCE[i_trk] = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bestMomentum();
          fwdMCS_NoSCE[i_trk] = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->fwdMomentum();
          bwdMCS_NoSCE[i_trk] = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bwdMomentum();
          bestMCSLL_NoSCE[i_trk] = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bestLogLikelihood();
          fwdMCSLL_NoSCE[i_trk] = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->fwdLogLikelihood();
          bwdMCSLL_NoSCE[i_trk] = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bwdLogLikelihood();

          //-- Track length and range momentum
          Trk_length[i_trk] = daughter_Tracks[i_trk]->Length();

          Mom_Range_mu[i_trk] = _trk_mom_calculator.GetTrackMomentum(Trk_length[i_trk], 13);
          Mom_Range_p[i_trk] = _trk_mom_calculator.GetTrackMomentum(Trk_length[i_trk], 2212);
          Mom_Range_pi[i_trk] = _trk_mom_calculator.GetTrackMomentum(Trk_length[i_trk], 211);

          //-- dE/dx info (pandoracaliSCE has E-field and spatial correction)
          auto assoCal = trackToCalAsso.at(daughter_Tracks[i_trk].key());
          if(assoCal.size()!=3){
            throw cet::exception("[Numu0pi0p]") << "Associated calorimetry does not match the three plane information!" << std::endl;
          }
          // induction 0 = 0, induction 1 = 1, collection = 2 (Check if the plane ID is correct)
          // assoCal[id_pl]->PlaneID().Plane == 2 (collection)
          // The vector is ordered by residual range from small to big (track end to track start)
          dEdx_pl0 = assoCal[0]->dEdx();
          resRange_pl0 = assoCal[0]->ResidualRange();

          dEdx_pl1 = assoCal[1]->dEdx();
          resRange_pl1 = assoCal[1]->ResidualRange();

          dEdx_pl2 = assoCal[2]->dEdx();
          resRange_pl2 = assoCal[2]->ResidualRange();

          auto hits_dEdx_size_pl0 = dEdx_pl0.size();
          auto hits_dEdx_size_pl1 = dEdx_pl1.size();
          auto hits_dEdx_size_pl2 = dEdx_pl2.size();

          auto half_size_pl0 = hits_dEdx_size_pl0 / 2;
          auto half_size_pl1 = hits_dEdx_size_pl1 / 2;
          auto half_size_pl2 = hits_dEdx_size_pl2 / 2;


          //- Average dE/dx close to the start and the end of the track
          if(half_size_pl0 != 0){
            dEdx_pl0_start_half[i_trk] = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
            dEdx_pl0_end_half[i_trk] = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
          }
          else{
            dEdx_pl0_start_half[i_trk] = 0;
            dEdx_pl0_end_half[i_trk] = 0;
          }
          if(half_size_pl1 != 0){
            dEdx_pl1_start_half[i_trk] = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
            dEdx_pl1_end_half[i_trk] = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
          }
          else{
            dEdx_pl1_start_half[i_trk] = 0;
            dEdx_pl1_end_half[i_trk] = 0;
          }
          if(half_size_pl2 != 0){
            dEdx_pl2_start_half[i_trk] = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
            dEdx_pl2_end_half[i_trk] = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;
          }
          else{
            dEdx_pl2_start_half[i_trk] = 0;
            dEdx_pl2_end_half[i_trk] = 0;
          }
          //- dEdx_1020
          if (dEdx_pl0.size()<=30) {
            dEdx_pl0_start1020[i_trk] = dEdx_pl0_start_half[i_trk];
            dEdx_pl0_end1020[i_trk] = dEdx_pl0_end_half[i_trk];
          }
          else{
            dEdx_pl0_start1020[i_trk] = std::accumulate(dEdx_pl0.end() - 20, dEdx_pl0.end() - 10, 0.) / 10.;
            dEdx_pl0_end1020[i_trk] = std::accumulate(dEdx_pl0.begin() + 10, dEdx_pl0.begin() + 20, 0.) / 10.;
          }
          if (dEdx_pl1.size()<=30) {
            dEdx_pl1_start1020[i_trk] = dEdx_pl1_start_half[i_trk];
            dEdx_pl1_end1020[i_trk] = dEdx_pl1_end_half[i_trk];
          }
          else{
            dEdx_pl1_start1020[i_trk] = std::accumulate(dEdx_pl1.end() - 20, dEdx_pl1.end() - 10, 0.) / 10.;
            dEdx_pl1_end1020[i_trk] = std::accumulate(dEdx_pl1.begin() + 10, dEdx_pl1.begin() + 20, 0.) / 10.;
          }
          if (dEdx_pl2.size()<=30) {
            dEdx_pl2_start1020[i_trk] = dEdx_pl2_start_half[i_trk];
            dEdx_pl2_end1020[i_trk] = dEdx_pl2_end_half[i_trk];
          }
          else{
            dEdx_pl2_start1020[i_trk] = std::accumulate(dEdx_pl2.end() - 20, dEdx_pl2.end() - 10, 0.) / 10.;
            dEdx_pl2_end1020[i_trk] = std::accumulate(dEdx_pl2.begin() + 10, dEdx_pl2.begin() + 20, 0.) / 10.;
          }

          if((dEdx_pl2_end_half[i_trk] + dEdx_pl2_start_half[i_trk]) != 0){
            dEdx_pl2_1020_ratio[i_trk] = dEdx_pl2_end1020[i_trk] / (dEdx_pl2_end1020[i_trk] + dEdx_pl2_start1020[i_trk]);
            dEdx_pl2_half_ratio[i_trk] = dEdx_pl2_end_half[i_trk] / (dEdx_pl2_end_half[i_trk] + dEdx_pl2_start_half[i_trk]);
          }
          else{
            dEdx_pl2_1020_ratio[i_trk] = 0;
            dEdx_pl2_half_ratio[i_trk] = 0;
          }

          //--PID
          PID3pl pid;
          pid.Chi2(PIDTotrackAsso,daughter_Tracks[i_trk], Trk_start_SCEcorr[i_trk], Trk_end_SCEcorr[i_trk],hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2);

          PID_Chi2Mu_3pl[i_trk] = pid.PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
          PID_Chi2P_3pl[i_trk] = pid.PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID
          PID_Chi2Pi_3pl[i_trk] = pid.PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID
          PID_Chi2K_3pl[i_trk] = pid.PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

          /////////////////
          // Fill True info from reco-truth matching
          ////////////////

          if(IsMC){
            std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(daughter_Tracks[i_trk].key());
            BackTrackerTruthMatching backtrackertruthmatch;
            backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
            auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
            trk_cosmic_percent[i_trk] = backtrackertruthmatch.ReturnCosmicPercent();
            trk_purity[i_trk] = backtrackertruthmatch.ReturnPurity();
            trk_completeness[i_trk] = backtrackertruthmatch.ReturnCompleteness();

            if(!MCparticle){
              if_cosmic[i_trk] = true;
              std::cout<<"MC particle does not exist!"<<std::endl;
            }
            else{
              if(trk_cosmic_percent[i_trk] > 0.8){
                if_cosmic[i_trk] = true;
              }
              else{
                if_cosmic[i_trk] = false;
                if(MCparticle->Process() == "primary" && trk_purity[i_trk] > 0.2 && trk_completeness[i_trk] > 0.2){
                  if_matchPrimary[i_trk] = true;
                  if(MCparticle->PdgCode() == 13){
                    if_matchMu[i_trk] = true;
                  }
                  if(MCparticle->PdgCode() == 2212){
                    if_matchP[i_trk] = true;
                  }
                }
                TVector3 TrueTrackPos(MCparticle->Px(), MCparticle->Py(), MCparticle->Pz());// The initial momentum represent the angle of true track
                true_mom[i_trk] = MCparticle->P();
                true_start_x[i_trk] = MCparticle->Position().X();
                true_start_y[i_trk] = MCparticle->Position().Y();
                true_start_z[i_trk] = MCparticle->Position().Z();
                true_end_x[i_trk] = MCparticle->EndPosition().X();
                true_end_y[i_trk] = MCparticle->EndPosition().Y();
                true_end_z[i_trk] = MCparticle->EndPosition().Z();
                TVector3 true_start(MCparticle->Position().X(), MCparticle->Position().Y(), MCparticle->Position().Z());
                TVector3 true_end(MCparticle->EndPosition().X(), MCparticle->EndPosition().Y(), MCparticle->EndPosition().Z());
                true_trk_ifcontained[i_trk] = _fiducial_volume.TrackContain(true_start, true_end);
                true_vtxFV[i_trk] = _fiducial_volume.VertexInFV(true_start);
                true_trk_phi[i_trk] = TrueTrackPos.Phi();
                true_trk_theta[i_trk] = TrueTrackPos.Theta();
                true_trk_costheta[i_trk] = cos(TrueTrackPos.Theta());
                true_trk_theta_yz[i_trk] = std::atan2(TrueTrackPos.Y(), TrueTrackPos.Z());
                true_trk_costheta_yz[i_trk] = cos(true_trk_theta_yz[i_trk]);
                true_trk_theta_xz[i_trk] = std::atan2(TrueTrackPos.X(), TrueTrackPos.Z());
                true_trk_costheta_xz[i_trk] = cos(true_trk_theta_xz[i_trk]);
                true_trk_length[i_trk] = (true_start - true_end).Mag(); // An estimation of true track length
                true_trk_PDG[i_trk] = MCparticle->PdgCode();

                // pass the MCparticle out to the MC truth
                selected_MCparticle[i_trk] = MCparticle;
              } // not cosmic
            } // If MC particle exists
          } // IsMC
        } // Finish looping the 2 tracks

        //////////////////////////
        // decide which track is muon and which is proton
        /////////////////////////

        //-- Track particle ID (Set the threshold to be PID_Chi2P_3pl 70 and 90)
        if(PID_Chi2P_3pl[0] < 70 && PID_Chi2P_3pl[1] > 90){
          muon_idx = 1;
          proton_idx = 0;
        }
        if(PID_Chi2P_3pl[0] > 90 && PID_Chi2P_3pl[1] < 70){
          muon_idx = 0;
          proton_idx = 1;
        }

        if(muon_idx >= 0 && proton_idx >= 0){
          //-- Vertex and Others
          vtx_x = Trk_start_SCEcorr[muon_idx].X();
          vtx_y = Trk_start_SCEcorr[muon_idx].Y();
          vtx_z = Trk_start_SCEcorr[muon_idx].Z();

          vtx_InFV = _fiducial_volume.VertexInFV(Trk_start_SCEcorr[muon_idx]);

          muon_cand_phi = daughter_Tracks[muon_idx]->Phi();
          muon_cand_theta = daughter_Tracks[muon_idx]->Theta();
          muon_cand_costheta = cos(daughter_Tracks[muon_idx]->Theta());
          muon_cand_mom_MCS = bestMCS[muon_idx];
          muon_cand_mom_range = Mom_Range_mu[muon_idx];
          muon_cand_length = Trk_length[muon_idx];

          proton_cand_phi = daughter_Tracks[proton_idx]->Phi();
          proton_cand_theta = daughter_Tracks[proton_idx]->Theta();
          proton_cand_costheta = cos(daughter_Tracks[proton_idx]->Theta());
          proton_cand_mom_range = Mom_Range_p[proton_idx];
          proton_cand_length = Trk_length[proton_idx];

          cos_ang_muon_proton = (trk_startdir[muon_idx] * trk_startdir[proton_idx]) / (trk_startdir[muon_idx].Mag() * trk_startdir[proton_idx].Mag());

          TVector3 muon_start(Trk_start_SCEcorr[muon_idx].X(), Trk_start_SCEcorr[muon_idx].Y(), Trk_start_SCEcorr[muon_idx].Z());
          TVector3 proton_start(Trk_start_SCEcorr[proton_idx].X(), Trk_start_SCEcorr[proton_idx].Y(), Trk_start_SCEcorr[proton_idx].Z());
          dist_mu_p_start = (muon_start - proton_start).Mag();

          //-- PT and PL
          float muon_PT_x = bestMCS[muon_idx] * sin(muon_cand_theta) * cos(muon_cand_phi);
          float muon_PT_y = bestMCS[muon_idx] * sin(muon_cand_theta) * sin(muon_cand_phi);
          float muon_PT_z = bestMCS[muon_idx] * cos(muon_cand_theta);

          float proton_PT_x = Mom_Range_p[proton_idx] * sin(proton_cand_theta) * cos(proton_cand_phi);
          float proton_PT_y = Mom_Range_p[proton_idx] * sin(proton_cand_theta) * sin(proton_cand_phi);
          float proton_PT_z = Mom_Range_p[proton_idx] * cos(proton_cand_theta);

          PT = sqrt((muon_PT_x + proton_PT_x) * (muon_PT_x + proton_PT_x) + (muon_PT_y + proton_PT_y) * (muon_PT_y + proton_PT_y));
          PL = abs(muon_PT_z + proton_PT_z);

          if(IsMC){
            if(if_matchMu[muon_idx] == true && if_matchP[proton_idx] == true){
              if_match_mu_p = true;
            }
            if(if_matchMu[proton_idx] == true && if_matchP[muon_idx] == true){
              if_match_mu_p_flipped = true;
            }
          }

        } // if there is a valid pair of muon and proton candidate
       
      } // if it is 2 single track in reco signature

      // In pandora, there would be only one neutrino slice maximum
      // If we already have one neutrino slice, we can ..
      break;

    } // if it is neutrino slice
  } // finish loop over all the pfp particles (looking for neutrino slice)

  if(!IsMC){
    my_event_->Fill();
  }
  //-------True
  else if(IsMC){

    // sanity check that there is MC truth info and they are reasonable
    if (MCTruthCollection.size() < 1){
      throw cet::exception("[Numu0pi0p]")<< "There is no MC truth info!" << std::endl;
    }
    if (if_GENIE_ReWeight){
      if(MCTruthCollection.size() != GENIEWeightCollection.size()){
        throw cet::exception("[Numu0pi0p]")<< "GENIE Reweight is not consistent with MCTruth!" << std::endl;
      }
    }
    if (if_G4_ReWeight){
      if(MCTruthCollection.size() != G4WeightCollection.size()){
        throw cet::exception("[Numu0pi0p]")<< "G4 Reweight is not consistent with MCTruth!" << std::endl;
      }
    }
    if (if_Flux_ReWeight){
      if(MCTruthCollection.size() != FluxWeightCollection.size()){
        throw cet::exception("[Numu0pi0p]")<< "Flux Reweight is not consistent with MCTruth!" << std::endl;
      }
    }

    // we first find the matching one if there is any, otherwise the the reco will be paired with the first MC interaction
    // set up which MC to fill with reco
    if_MC_reco_match = {false, false};
    // if there is no selected MC particle (the selected one is cosmic or so), which should be for the most of the cases
    if (if_cosmic[0] == true && if_cosmic[1] == true){
      matched_MCid = 0;
      //0. cosmic (there is MC, but no reco)
      MCmatch_type = 0;
    }
    else {
      // define MC truth to MC particle association when it is needed
      art::FindMany<simb::MCParticle,sim::GeneratedParticleInfo> MCtToMCpAsso(Handle_MCTruth, evt, m_geantLabel);
      for (unsigned int i_trk = 0; i_trk < 2; i_trk++){

        for (unsigned int i_mc = 0; i_mc < MCTruthCollection.size(); i_mc++){
          std::vector<const simb::MCParticle*> assoMC = MCtToMCpAsso.at(MCTruthCollection[i_mc].key());
          for (unsigned int i_mcp = 0; i_mcp < assoMC.size(); i_mcp++){
            if (assoMC[i_mcp] == selected_MCparticle[i_trk].get()){
              if_MC_reco_match[i_trk] = true;
              matched_MCid_2trk[i_trk] = i_mc;
              break;
            }
          }
          if (if_MC_reco_match[i_trk]) break;
        } // MCtruth loop
      } // 2trk
    } // not cosmic

    // 1. neither the tracks matches any MC particle, set it to the first MC interaction (cosmic)
    if (if_MC_reco_match[0] == false && if_MC_reco_match[1] == false){
      matched_MCid = 0;
      MCmatch_type = 1;
    }
    // 2. one of the tracks matches a MC particle, and the other not, set it to the matched one
    if (if_MC_reco_match[0] == true && if_MC_reco_match[1] == false){
      matched_MCid = matched_MCid_2trk[0];
      if_match_mu_p = false;
      MCmatch_type = 2;
    }
    if (if_MC_reco_match[0] == false && if_MC_reco_match[1] == true){
      matched_MCid = matched_MCid_2trk[1];
      if_match_mu_p = false;
      MCmatch_type = 2;
    }
    // 3. the two tracks match two different MC truth, set it to the one that corresponds to the first track (rather arbitary:()
    if (if_MC_reco_match[0] == true && if_MC_reco_match[1] == true && matched_MCid_2trk[0] != matched_MCid_2trk[1]){
      matched_MCid = matched_MCid_2trk[0];
      if_match_mu_p = false;
      MCmatch_type = 3;
    }
    // 4. the two tracks match the same MC truth, set it to that one!
    if (if_MC_reco_match[0] == true && if_MC_reco_match[1] == true && matched_MCid_2trk[0] == matched_MCid_2trk[1]){
      matched_MCid = matched_MCid_2trk[0];
      MCmatch_type = 4;
    }

    //sanity check that matched_MCid is valid
    if (matched_MCid < 0 || matched_MCid >= MCTruthCollection.size()){
      throw cet::exception("[Numu0pi0p]")<< "Something is wrong at the MC reco matching stage!" << std::endl;
    }

    ////////////////////////////////////////////////////
    // [the MC to be paired with the reco is decided] Fill a pair of reco and truth together (we already checked that there's at least one MC)
    if (if_GENIE_ReWeight){
      //-- Xsec reweight
      std::map<std::string, std::vector<double>> evtwgt_map = GENIEWeightCollection.at(matched_MCid)->fWeight;
      for (const auto& this_par: genie_pars){
        const std::vector<double> &weights = evtwgt_map.at(this_par);

        if (this_par == "MaCCQE_UBGenie") { MaCCQE = weights; }
        if (this_par == "CoulombCCQE_UBGenie") { CoulombCCQE = weights; }
        if (this_par == "MaNCEL_UBGenie") { MaNCEL = weights; }
        if (this_par == "EtaNCEL_UBGenie") { EtaNCEL = weights; }

        if (this_par == "NormCCMEC_UBGenie") { NormCCMEC = weights; }
        if (this_par == "NormNCMEC_UBGenie") { NormNCMEC = weights; }
        if (this_par == "FracPN_CCMEC_UBGenie") { FracPN_CCMEC = weights; }
        if (this_par == "FracDelta_CCMEC_UBGenie") { FracDelta_CCMEC = weights; }

        if (this_par == "MaCCRES_UBGenie") { MaCCRES = weights; }
        if (this_par == "MvCCRES_UBGenie") { MvCCRES = weights; }
        if (this_par == "MaNCRES_UBGenie") { MaNCRES = weights; }
        if (this_par == "MvNCRES_UBGenie") { MvNCRES = weights; }

        if (this_par == "NonRESBGvpCC1pi_UBGenie") { NonRESBGvpCC1pi = weights; }
        if (this_par == "NonRESBGvpCC2pi_UBGenie") { NonRESBGvpCC2pi = weights; }
        if (this_par == "NonRESBGvpNC1pi_UBGenie") { NonRESBGvpNC1pi = weights; }
        if (this_par == "NonRESBGvpNC2pi_UBGenie") { NonRESBGvpNC2pi = weights; }
        if (this_par == "NonRESBGvnCC1pi_UBGenie") { NonRESBGvnCC1pi = weights; }
        if (this_par == "NonRESBGvnCC2pi_UBGenie") { NonRESBGvnCC2pi = weights; }
        if (this_par == "NonRESBGvnNC1pi_UBGenie") { NonRESBGvnNC1pi = weights; }
        if (this_par == "NonRESBGvnNC2pi_UBGenie") { NonRESBGvnNC2pi = weights; }
        if (this_par == "NonRESBGvbarpCC1pi_UBGenie") { NonRESBGvbarpCC1pi = weights; }
        if (this_par == "NonRESBGvbarpCC2pi_UBGenie") { NonRESBGvbarpCC2pi = weights; }
        if (this_par == "NonRESBGvbarpNC1pi_UBGenie") { NonRESBGvbarpNC1pi = weights; }
        if (this_par == "NonRESBGvbarpNC2pi_UBGenie") { NonRESBGvbarpNC2pi = weights; }
        if (this_par == "NonRESBGvbarnCC1pi_UBGenie") { NonRESBGvbarnCC1pi = weights; }
        if (this_par == "NonRESBGvbarnCC2pi_UBGenie") { NonRESBGvbarnCC2pi = weights; }
        if (this_par == "NonRESBGvbarnNC1pi_UBGenie") { NonRESBGvbarnNC1pi = weights; }
        if (this_par == "NonRESBGvbarnNC2pi_UBGenie") { NonRESBGvbarnNC2pi = weights; }
        if (this_par == "AhtBY_UBGenie") { AhtBY = weights; }
        if (this_par == "BhtBY_UBGenie") { BhtBY = weights; }
        if (this_par == "CV1uBY_UBGenie") { CV1uBY = weights; }
        if (this_par == "CV2uBY_UBGenie") { CV2uBY = weights; }

        if (this_par == "AGKYxF1pi_UBGenie") { AGKYxF1pi = weights; }
        if (this_par == "AGKYpT1pi_UBGenie") { AGKYpT1pi = weights; }

        if (this_par == "MFP_pi_UBGenie") { MFP_pi = weights; }
        if (this_par == "MFP_N_UBGenie") { MFP_N = weights; }
        if (this_par == "FrCEx_pi_UBGenie") { FrCEx_pi = weights; }
        if (this_par == "FrInel_pi_UBGenie") { FrInel_pi = weights; }
        if (this_par == "FrAbs_pi_UBGenie") { FrAbs_pi = weights; }
        if (this_par == "FrCEx_N_UBGenie") { FrCEx_N = weights; }
        if (this_par == "FrInel_N_UBGenie") { FrInel_N = weights; }
        if (this_par == "FrAbs_N_UBGenie") { FrAbs_N = weights; }

        if (this_par == "RDecBR1gamma_UBGenie") { RDecBR1gamma = weights; }
        if (this_par == "RDecBR1eta_UBGenie") { RDecBR1eta = weights; }

        // additional fix
        if (this_par == "RootinoFix_UBGenie") { RootinoFix = weights; }
        if (this_par == "splines_general_Spline") { splines_general_Spline = weights; }

        // for fake data
        if (this_par == "MaCCQE_tune_UBGenie") { MaCCQE_tune = weights; }
        if (this_par == "RPA_CCQE_tune_UBGenie") { RPA_CCQE_tune = weights; }
        if (this_par == "NormCCMEC_tune_UBGenie") { NormCCMEC_tune = weights; }
        if (this_par == "XSecShape_CCMEC_tune_UBGenie") { XSecShape_CCMEC_tune = weights; }

      } // GENIE
    } // if GENIE rewight

    if (if_Flux_ReWeight){
      std::map<std::string, std::vector<double>> Fluxevtwgt_map = FluxWeightCollection.at(matched_MCid)->fWeight;
      for (const auto& this_par: flux_pars){
        const std::vector<double> &weights = Fluxevtwgt_map.at(this_par);

        if (this_par == "expskin_FluxUnisim") { expskin = weights; }
        if (this_par == "horncurrent_FluxUnisim") { horncurrent = weights; }
        if (this_par == "nucleoninexsec_FluxUnisim") { nucleoninexsec = weights; }
        if (this_par == "nucleonqexsec_FluxUnisim") { nucleonqexsec = weights; }
        if (this_par == "nucleontotxsec_FluxUnisim") { nucleontotxsec = weights; }
        if (this_par == "pioninexsec_FluxUnisim") { pioninexsec = weights; }
        if (this_par == "pionqexsec_FluxUnisim") { pionqexsec = weights; }
        if (this_par == "piontotxsec_FluxUnisim") { piontotxsec = weights; }

      } // flux
      //for (auto it : Fluxevtwgt_map){
      //  std::string func_name = it.first;
      //  std::vector<double> weight_v = it.second;
      //  std::cout<<"--------func_name: " << func_name <<"---------"<<std::endl;
      //  for (size_t i_w = 0; i_w < weight_v.size();i_w++){
      //    std::cout<<"i_w: "<< i_w <<" | weight: "<< weight_v.at(i_w)<<std::endl;
      //  }
      //}
    } // if flux reweight
    if (if_G4_ReWeight){
      std::map<std::string, std::vector<double>> G4evtwgt_map = G4WeightCollection.at(matched_MCid)->fWeight;
      //for(auto it : G4evtwgt_map){
      //  std::string func_name = it.first;
      //  std::vector<double> weight_v = it.second;
      //  std::cout<<"--------func_name: " << func_name <<"---------"<<std::endl;
      //  for (size_t i_w = 0; i_w < weight_v.size();i_w++){
      //    std::cout<<"i_w: "<< i_w <<" | weight: "<< weight_v.at(i_w)<<std::endl;
      //  }
      //} 
      for (const auto& this_par: g4_pars){
        const std::vector<double> &g4_weights = G4evtwgt_map.at(this_par);
        if (this_par == "reinteractions_piplusReacLow_Geant4") { piplusReacLow = g4_weights; }
        if (this_par == "reinteractions_piplusReacHigh_Geant4") { piplusReacHigh = g4_weights; }
        if (this_par == "reinteractions_piplusAbs_Geant4") { piplusAbs = g4_weights; }
        if (this_par == "reinteractions_piplusCex_Geant4") { piplusCex = g4_weights; }
        if (this_par == "reinteractions_piplusDCex_Geant4") { piplusDCex = g4_weights; }
        if (this_par == "reinteractions_piplusPiProd_Geant4") { piplusPiProd = g4_weights; }
        if (this_par == "reinteractions_piplusElast_Geant4") { piplusElast = g4_weights; }

        if (this_par == "reinteractions_piminusReacLow_Geant4") { piminusReacLow = g4_weights; }
        if (this_par == "reinteractions_piminusReacHigh_Geant4") { piminusReacHigh = g4_weights; }
        if (this_par == "reinteractions_piminusAbs_Geant4") { piminusAbs = g4_weights; }
        if (this_par == "reinteractions_piminusCex_Geant4") { piminusCex = g4_weights; }
        if (this_par == "reinteractions_piminusDCex_Geant4") { piminusDCex = g4_weights; }
        if (this_par == "reinteractions_piminusPiProd_Geant4") { piminusPiProd = g4_weights; }
        if (this_par == "reinteractions_piminusElast_Geant4") { piminusElast = g4_weights; }

        if (this_par == "reinteractions_protonReac_Geant4") { protonReac = g4_weights; }
        if (this_par == "reinteractions_protonElast_Geant4") { protonElast = g4_weights; }

        if (this_par == "reinteractions_neutronAbs_max_Geant4") { neutronAbs_max = g4_weights; }
        if (this_par == "reinteractions_neutronElast_max_Geant4") { neutronElast_max = g4_weights; }
        if (this_par == "reinteractions_neutronInela_max_Geant4") { neutronInela_max = g4_weights; }
        if (this_par == "reinteractions_neutronPiprod_max_Geant4") { neutronPiprod_max = g4_weights; }
        if (this_par == "reinteractions_neutronPprod_max_Geant4") { neutronPprod_max = g4_weights; }

        if (this_par == "reinteractions_neutronAbs_100_Geant4") { neutronAbs_100 = g4_weights; }
        if (this_par == "reinteractions_neutronElast_100_Geant4") { neutronElast_100 = g4_weights; }
        if (this_par == "reinteractions_neutronInela_100_Geant4") { neutronInela_100 = g4_weights; }
        if (this_par == "reinteractions_neutronPiprod_100_Geant4") { neutronPiprod_100 = g4_weights; }
        if (this_par == "reinteractions_neutronPprod_100_Geant4") { neutronPprod_100 = g4_weights; }

        if (this_par == "reinteractions_neutronAbs_200_Geant4") { neutronAbs_200 = g4_weights; }
        if (this_par == "reinteractions_neutronElast_200_Geant4") { neutronElast_200 = g4_weights; }
        if (this_par == "reinteractions_neutronInela_200_Geant4") { neutronInela_200 = g4_weights; }
        if (this_par == "reinteractions_neutronPiprod_200_Geant4") { neutronPiprod_200 = g4_weights; }
        if (this_par == "reinteractions_neutronPprod_200_Geant4") { neutronPprod_200 = g4_weights; }

      } // G4 parameters
    } // if G4 reweight

    // The neutrino info
    nu_parent_pdg = MCFluxCollection[matched_MCid]->fptype;
    if (MCTruthCollection[matched_MCid]->Origin() == Neutrino_Origin){ MC_beamNeutrino = true;}
    MC_int_mode = MCTruthCollection[matched_MCid]->GetNeutrino().Mode();
    MC_nupdg = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().PdgCode();
    MC_ccnc = MCTruthCollection[matched_MCid]->GetNeutrino().CCNC();
    MC_Q2 = MCTruthCollection[matched_MCid]->GetNeutrino().QSqr();

    MC_nu_E = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().E();
    MC_nuVtxX = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().Vx();
    MC_nuVtxY = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().Vy();
    MC_nuVtxZ = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().Vz();

    TVector3 true_nuVtx(MC_nuVtxX, MC_nuVtxY, MC_nuVtxZ);
    MC_FV = _fiducial_volume.VertexInFV(true_nuVtx);
    MC_if_in_active = _fiducial_volume.VertexInActive(true_nuVtx);

    // The final state particles
    // Loop all the MCParticles to determine the true topology (all the MCParticles are from the neutrino events in overlay)
    if (MC_ccnc == 0 && MC_nupdg == 14 && MC_beamNeutrino == true){
      for(unsigned int i_mcp = 0; i_mcp < MCParticleCollection.size(); i_mcp++){
        if(MCParticleCollection[i_mcp]->Process() == "primary"){

          // PDG and momemtum of neutrino daughters
          MC_Primary_PDG.push_back(MCParticleCollection[i_mcp]->PdgCode());
          MC_Primary_Mom.push_back(MCParticleCollection[i_mcp]->P());

          // muon
          if(MCParticleCollection[i_mcp]->PdgCode() == 13){
            MC_nMuon++;

            MC_muon_true_Mom.push_back(MCParticleCollection[i_mcp]->P());
            MC_muon_true_theta.push_back(acos(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P()));
            MC_muon_true_cos_theta.push_back(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P());
            MC_muon_true_phi.push_back(atan2(MCParticleCollection[i_mcp]->Py(), MCParticleCollection[i_mcp]->Px()));

            MC_muon_true_Px.push_back(MCParticleCollection[i_mcp]->Px());
            MC_muon_true_Py.push_back(MCParticleCollection[i_mcp]->Py());
            MC_muon_true_Pz.push_back(MCParticleCollection[i_mcp]->Pz());
            MC_muon_true_E.push_back(MCParticleCollection[i_mcp]->E());

          }
          // electron
          if(MCParticleCollection[i_mcp]->PdgCode() == 11) MC_nElectron++;
          // neutron
          if(MCParticleCollection[i_mcp]->PdgCode() == 2112) MC_nNeutron++;
          // proton
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() < 0.260) MC_nProton_below260++;
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() >= 0.260) {
            MC_nProton_above260++;

            MC_proton_true_Mom.push_back(MCParticleCollection[i_mcp]->P());
            MC_proton_true_theta.push_back(acos(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P()));
            MC_proton_true_cos_theta.push_back(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P());
            MC_proton_true_phi.push_back(atan2(MCParticleCollection[i_mcp]->Py(), MCParticleCollection[i_mcp]->Px()));

            MC_proton_true_Px.push_back(MCParticleCollection[i_mcp]->Px());
            MC_proton_true_Py.push_back(MCParticleCollection[i_mcp]->Py());
            MC_proton_true_Pz.push_back(MCParticleCollection[i_mcp]->Pz());
          }
          // pion0
          if(MCParticleCollection[i_mcp]->PdgCode() == 111) MC_nPi0++;
          // pion+
          if(MCParticleCollection[i_mcp]->PdgCode() == 211 && MCParticleCollection[i_mcp]->P() < 0.080) MC_nPiPlus_below80++;
          if(MCParticleCollection[i_mcp]->PdgCode() == 211 && MCParticleCollection[i_mcp]->P() >= 0.080) MC_nPiPlus_above80++;
          // pion-
          if(MCParticleCollection[i_mcp]->PdgCode() == -211 && MCParticleCollection[i_mcp]->P() < 0.080) MC_nPiMinus_below80++;
          if(MCParticleCollection[i_mcp]->PdgCode() == -211 && MCParticleCollection[i_mcp]->P() >= 0.080) MC_nPiMinus_above80++;
        }
      }
    }

    Topology topology;
    TopologyType = topology.TopologyLabel(MC_nMuon, MC_nElectron, MC_nPiPlus_above80, MC_nPiPlus_below80, MC_nPiMinus_above80, MC_nPiMinus_below80, MC_nPi0, MC_nProton_above260, MC_nProton_below260, MC_nupdg, MC_ccnc, MC_beamNeutrino, MC_FV);
    // [CC0pi0p] Store true kinematic information for efficiency study
    if(TopologyType == 1){

      MC_0pi0p_muon_mom = MC_muon_true_Mom.front();
      MC_0pi0p_muon_theta = MC_muon_true_theta.front();
      MC_0pi0p_muon_costheta = MC_muon_true_cos_theta.front();
      MC_0pi0p_muon_phi = MC_muon_true_phi.front();

    }

    if(TopologyType == 2){

      MC_0pi1p_muon_mom = MC_muon_true_Mom.front();
      MC_0pi1p_muon_theta = MC_muon_true_theta.front();
      MC_0pi1p_muon_costheta = MC_muon_true_cos_theta.front();
      MC_0pi1p_muon_phi = MC_muon_true_phi.front();

      MC_0pi1p_proton_mom = MC_proton_true_Mom.front();
      MC_0pi1p_proton_theta = MC_proton_true_theta.front();
      MC_0pi1p_proton_costheta = MC_proton_true_cos_theta.front();
      MC_0pi1p_proton_phi = MC_proton_true_phi.front();

      TVector3 muon_dir(MC_muon_true_Px.front(), MC_muon_true_Py.front(), MC_muon_true_Pz.front());
      TVector3 proton_dir(MC_proton_true_Px.front(), MC_proton_true_Py.front(), MC_proton_true_Pz.front());
      MC_0pi1p_cos_ang_muon_proton = (muon_dir * proton_dir) / (muon_dir.Mag() * proton_dir.Mag());

      MC_0pi1p_PT = sqrt((MC_muon_true_Px.front() + MC_proton_true_Px.front()) * (MC_muon_true_Px.front() + MC_proton_true_Px.front()) + (MC_muon_true_Py.front() + MC_proton_true_Py.front()) * (MC_muon_true_Py.front() + MC_proton_true_Py.front()));
      MC_0pi1p_PL = abs(MC_muon_true_Pz.front() + MC_proton_true_Pz.front());
    }

    if(MC_nMuon == 1){
      MC_transfer_E = MC_nu_E - MC_muon_true_E.front();
    }

    // vtx bias
    reco_MC_dist_vtx = (true_nuVtx - Trk_start_SCEcorr[muon_idx]).Mag();
    reco_MC_dist_vtx_noSCE = (true_nuVtx - Trk_start[muon_idx]).Mag();

    my_event_->Fill();

    // [wipe the reco info] if if_2tracks is false, it is not selected. Therefore, the reco info won't be considered
    if_2tracks = false;
    MCmatch_type = 0; // no reco info is equivalent to cosmic
    reco_MC_dist_vtx_noSCE = -999;
    reco_MC_dist_vtx = -999;

    // If there is more MC to fill
    if (MCTruthCollection.size() > 1){
      for (unsigned int i_mc = 0; i_mc < MCTruthCollection.size(); i_mc++){
        if (i_mc != matched_MCid){

          if (if_GENIE_ReWeight){
            //-- Xsec reweight
            std::map<std::string, std::vector<double>> evtwgt_map = GENIEWeightCollection.at(i_mc)->fWeight;
            for (const auto& this_par: genie_pars){
              const std::vector<double> &weights = evtwgt_map.at(this_par);
              if (this_par == "MaCCQE_UBGenie") { MaCCQE = weights; }
              if (this_par == "CoulombCCQE_UBGenie") { CoulombCCQE = weights; }
              if (this_par == "MaNCEL_UBGenie") { MaNCEL = weights; }
              if (this_par == "EtaNCEL_UBGenie") { EtaNCEL = weights; }

              if (this_par == "NormCCMEC_UBGenie") { NormCCMEC = weights; }
              if (this_par == "NormNCMEC_UBGenie") { NormNCMEC = weights; }
              if (this_par == "FracPN_CCMEC_UBGenie") { FracPN_CCMEC = weights; }
              if (this_par == "FracDelta_CCMEC_UBGenie") { FracDelta_CCMEC = weights; }

              if (this_par == "MaCCRES_UBGenie") { MaCCRES = weights; }
              if (this_par == "MvCCRES_UBGenie") { MvCCRES = weights; }
              if (this_par == "MaNCRES_UBGenie") { MaNCRES = weights; }
              if (this_par == "MvNCRES_UBGenie") { MvNCRES = weights; }

              if (this_par == "NonRESBGvpCC1pi_UBGenie") { NonRESBGvpCC1pi = weights; }
              if (this_par == "NonRESBGvpCC2pi_UBGenie") { NonRESBGvpCC2pi = weights; }
              if (this_par == "NonRESBGvpNC1pi_UBGenie") { NonRESBGvpNC1pi = weights; }
              if (this_par == "NonRESBGvpNC2pi_UBGenie") { NonRESBGvpNC2pi = weights; }
              if (this_par == "NonRESBGvnCC1pi_UBGenie") { NonRESBGvnCC1pi = weights; }
              if (this_par == "NonRESBGvnCC2pi_UBGenie") { NonRESBGvnCC2pi = weights; }
              if (this_par == "NonRESBGvnNC1pi_UBGenie") { NonRESBGvnNC1pi = weights; }
              if (this_par == "NonRESBGvnNC2pi_UBGenie") { NonRESBGvnNC2pi = weights; }
              if (this_par == "NonRESBGvbarpCC1pi_UBGenie") { NonRESBGvbarpCC1pi = weights; }
              if (this_par == "NonRESBGvbarpCC2pi_UBGenie") { NonRESBGvbarpCC2pi = weights; }
              if (this_par == "NonRESBGvbarpNC1pi_UBGenie") { NonRESBGvbarpNC1pi = weights; }
              if (this_par == "NonRESBGvbarpNC2pi_UBGenie") { NonRESBGvbarpNC2pi = weights; }
              if (this_par == "NonRESBGvbarnCC1pi_UBGenie") { NonRESBGvbarnCC1pi = weights; }
              if (this_par == "NonRESBGvbarnCC2pi_UBGenie") { NonRESBGvbarnCC2pi = weights; }
              if (this_par == "NonRESBGvbarnNC1pi_UBGenie") { NonRESBGvbarnNC1pi = weights; }
              if (this_par == "NonRESBGvbarnNC2pi_UBGenie") { NonRESBGvbarnNC2pi = weights; }
              if (this_par == "AhtBY_UBGenie") { AhtBY = weights; }
              if (this_par == "BhtBY_UBGenie") { BhtBY = weights; }
              if (this_par == "CV1uBY_UBGenie") { CV1uBY = weights; }
              if (this_par == "CV2uBY_UBGenie") { CV2uBY = weights; }

              if (this_par == "AGKYxF1pi_UBGenie") { AGKYxF1pi = weights; }
              if (this_par == "AGKYpT1pi_UBGenie") { AGKYpT1pi = weights; }

              if (this_par == "MFP_pi_UBGenie") { MFP_pi = weights; }
              if (this_par == "MFP_N_UBGenie") { MFP_N = weights; }
              if (this_par == "FrCEx_pi_UBGenie") { FrCEx_pi = weights; }
              if (this_par == "FrInel_pi_UBGenie") { FrInel_pi = weights; }
              if (this_par == "FrAbs_pi_UBGenie") { FrAbs_pi = weights; }
              if (this_par == "FrCEx_N_UBGenie") { FrCEx_N = weights; }
              if (this_par == "FrInel_N_UBGenie") { FrInel_N = weights; }
              if (this_par == "FrAbs_N_UBGenie") { FrAbs_N = weights; }

              if (this_par == "RDecBR1gamma_UBGenie") { RDecBR1gamma = weights; }
              if (this_par == "RDecBR1eta_UBGenie") { RDecBR1eta = weights; }

              // additional fix
              if (this_par == "RootinoFix_UBGenie") { RootinoFix = weights; }
              if (this_par == "splines_general_Spline") { splines_general_Spline = weights; }

              // for fake data
              if (this_par == "MaCCQE_tune_UBGenie") { MaCCQE_tune = weights; }
              if (this_par == "RPA_CCQE_tune_UBGenie") { RPA_CCQE_tune = weights; }
              if (this_par == "NormCCMEC_tune_UBGenie") { NormCCMEC_tune = weights; }
              if (this_par == "XSecShape_CCMEC_tune_UBGenie") { XSecShape_CCMEC_tune = weights; }
            } // GENIE
          } // if GENIE rewight

          if (if_Flux_ReWeight){
            std::map<std::string, std::vector<double>> Fluxevtwgt_map = FluxWeightCollection.at(i_mc)->fWeight;
            for (const auto& this_par: flux_pars){
              const std::vector<double> &weights = Fluxevtwgt_map.at(this_par);

              if (this_par == "expskin_FluxUnisim") { expskin = weights; }
              if (this_par == "horncurrent_FluxUnisim") { horncurrent = weights; }
              if (this_par == "nucleoninexsec_FluxUnisim") { nucleoninexsec = weights; }
              if (this_par == "nucleonqexsec_FluxUnisim") { nucleonqexsec = weights; }
              if (this_par == "nucleontotxsec_FluxUnisim") { nucleontotxsec = weights; }
              if (this_par == "pioninexsec_FluxUnisim") { pioninexsec = weights; }
              if (this_par == "pionqexsec_FluxUnisim") { pionqexsec = weights; }
              if (this_par == "piontotxsec_FluxUnisim") { piontotxsec = weights; }

            }
          }

          if (if_G4_ReWeight){
            std::map<std::string, std::vector<double>> G4evtwgt_map = G4WeightCollection.at(i_mc)->fWeight;
            for (const auto& this_par: g4_pars){
              const std::vector<double> &g4_weights = G4evtwgt_map.at(this_par);
                if (this_par == "reinteractions_piplusReacLow_Geant4") { piplusReacLow = g4_weights; }
                if (this_par == "reinteractions_piplusReacHigh_Geant4") { piplusReacHigh = g4_weights; }
                if (this_par == "reinteractions_piplusAbs_Geant4") { piplusAbs = g4_weights; }
                if (this_par == "reinteractions_piplusCex_Geant4") { piplusCex = g4_weights; }
                if (this_par == "reinteractions_piplusDCex_Geant4") { piplusDCex = g4_weights; }
                if (this_par == "reinteractions_piplusPiProd_Geant4") { piplusPiProd = g4_weights; }
                if (this_par == "reinteractions_piplusElast_Geant4") { piplusElast = g4_weights; }

                if (this_par == "reinteractions_piminusReacLow_Geant4") { piminusReacLow = g4_weights; }
                if (this_par == "reinteractions_piminusReacHigh_Geant4") { piminusReacHigh = g4_weights; }
                if (this_par == "reinteractions_piminusAbs_Geant4") { piminusAbs = g4_weights; }
                if (this_par == "reinteractions_piminusCex_Geant4") { piminusCex = g4_weights; }
                if (this_par == "reinteractions_piminusDCex_Geant4") { piminusDCex = g4_weights; }
                if (this_par == "reinteractions_piminusPiProd_Geant4") { piminusPiProd = g4_weights; }
                if (this_par == "reinteractions_piminusElast_Geant4") { piminusElast = g4_weights; }

                if (this_par == "reinteractions_protonReac_Geant4") { protonReac = g4_weights; }
                if (this_par == "reinteractions_protonElast_Geant4") { protonElast = g4_weights; }

                if (this_par == "reinteractions_neutronAbs_max_Geant4") { neutronAbs_max = g4_weights; }
                if (this_par == "reinteractions_neutronElast_max_Geant4") { neutronElast_max = g4_weights; }
                if (this_par == "reinteractions_neutronInela_max_Geant4") { neutronInela_max = g4_weights; }
                if (this_par == "reinteractions_neutronPiprod_max_Geant4") { neutronPiprod_max = g4_weights; }
                if (this_par == "reinteractions_neutronPprod_max_Geant4") { neutronPprod_max = g4_weights; }

                if (this_par == "reinteractions_neutronAbs_100_Geant4") { neutronAbs_100 = g4_weights; }
                if (this_par == "reinteractions_neutronElast_100_Geant4") { neutronElast_100 = g4_weights; }
                if (this_par == "reinteractions_neutronInela_100_Geant4") { neutronInela_100 = g4_weights; }
                if (this_par == "reinteractions_neutronPiprod_100_Geant4") { neutronPiprod_100 = g4_weights; }
                if (this_par == "reinteractions_neutronPprod_100_Geant4") { neutronPprod_100 = g4_weights; }

                if (this_par == "reinteractions_neutronAbs_200_Geant4") { neutronAbs_200 = g4_weights; }
                if (this_par == "reinteractions_neutronElast_200_Geant4") { neutronElast_200 = g4_weights; }
                if (this_par == "reinteractions_neutronInela_200_Geant4") { neutronInela_200 = g4_weights; }
                if (this_par == "reinteractions_neutronPiprod_200_Geant4") { neutronPiprod_200 = g4_weights; }
                if (this_par == "reinteractions_neutronPprod_200_Geant4") { neutronPprod_200 = g4_weights; }

            } //g4 par
          } // if g4 reweight

          // The neutrino info
          nu_parent_pdg = MCFluxCollection[i_mc]->fptype;
          if (MCTruthCollection[i_mc]->Origin() == Neutrino_Origin){ MC_beamNeutrino = true;}
          MC_int_mode = MCTruthCollection[i_mc]->GetNeutrino().Mode();
          MC_nupdg = MCTruthCollection[i_mc]->GetNeutrino().Nu().PdgCode();
          MC_ccnc = MCTruthCollection[i_mc]->GetNeutrino().CCNC();
          MC_Q2 = MCTruthCollection[i_mc]->GetNeutrino().QSqr();

          MC_nu_E = MCTruthCollection[i_mc]->GetNeutrino().Nu().E();
          MC_nuVtxX = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vx();
          MC_nuVtxY = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vy();
          MC_nuVtxZ = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vz();

          TVector3 true_nuVtx(MC_nuVtxX, MC_nuVtxY, MC_nuVtxZ);
          MC_FV = _fiducial_volume.VertexInFV(true_nuVtx);
          MC_if_in_active = _fiducial_volume.VertexInActive(true_nuVtx);

          // The final state particles
          // Loop all the MCParticles to determine the true topology (all the MCParticles are from the neutrino events in overlay)
          if (MC_ccnc == 0 && MC_nupdg == 14 && MC_beamNeutrino == true){
            for(unsigned int i_mcp = 0; i_mcp < MCParticleCollection.size(); i_mcp++){
              if(MCParticleCollection[i_mcp]->Process() == "primary"){

                // PDG and momemtum of neutrino daughters
                MC_Primary_PDG.push_back(MCParticleCollection[i_mcp]->PdgCode());
                MC_Primary_Mom.push_back(MCParticleCollection[i_mcp]->P());

                // muon
                if(MCParticleCollection[i_mcp]->PdgCode() == 13){
                  MC_nMuon++;

                  MC_muon_true_Mom.push_back(MCParticleCollection[i_mcp]->P());
                  MC_muon_true_theta.push_back(acos(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P()));
                  MC_muon_true_cos_theta.push_back(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P());
                  MC_muon_true_phi.push_back(atan2(MCParticleCollection[i_mcp]->Py(), MCParticleCollection[i_mcp]->Px()));

                  MC_muon_true_Px.push_back(MCParticleCollection[i_mcp]->Px());
                  MC_muon_true_Py.push_back(MCParticleCollection[i_mcp]->Py());
                  MC_muon_true_Pz.push_back(MCParticleCollection[i_mcp]->Pz());
                  MC_muon_true_E.push_back(MCParticleCollection[i_mcp]->E());

                }
                // electron
                if(MCParticleCollection[i_mcp]->PdgCode() == 11) MC_nElectron++;
                // neutron
                if(MCParticleCollection[i_mcp]->PdgCode() == 2112) MC_nNeutron++;
                // proton
                if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() < 0.260) MC_nProton_below260++;
                if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() >= 0.260) {
                  MC_nProton_above260++;

                  MC_proton_true_Mom.push_back(MCParticleCollection[i_mcp]->P());
                  MC_proton_true_theta.push_back(acos(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P()));
                  MC_proton_true_cos_theta.push_back(MCParticleCollection[i_mcp]->Pz() / MCParticleCollection[i_mcp]->P());
                  MC_proton_true_phi.push_back(atan2(MCParticleCollection[i_mcp]->Py(), MCParticleCollection[i_mcp]->Px()));

                  MC_proton_true_Px.push_back(MCParticleCollection[i_mcp]->Px());
                  MC_proton_true_Py.push_back(MCParticleCollection[i_mcp]->Py());
                  MC_proton_true_Pz.push_back(MCParticleCollection[i_mcp]->Pz());
                }
                // pion0
                if(MCParticleCollection[i_mcp]->PdgCode() == 111) MC_nPi0++;
                // pion+
                if(MCParticleCollection[i_mcp]->PdgCode() == 211 && MCParticleCollection[i_mcp]->P() < 0.080) MC_nPiPlus_below80++;
                if(MCParticleCollection[i_mcp]->PdgCode() == 211 && MCParticleCollection[i_mcp]->P() >= 0.080) MC_nPiPlus_above80++;
                // pion-
                if(MCParticleCollection[i_mcp]->PdgCode() == -211 && MCParticleCollection[i_mcp]->P() < 0.080) MC_nPiMinus_below80++;
                if(MCParticleCollection[i_mcp]->PdgCode() == -211 && MCParticleCollection[i_mcp]->P() >= 0.080) MC_nPiMinus_above80++;
              }
            }
          }

          Topology topology;
          TopologyType = topology.TopologyLabel(MC_nMuon, MC_nElectron, MC_nPiPlus_above80, MC_nPiPlus_below80, MC_nPiMinus_above80, MC_nPiMinus_below80, MC_nPi0, MC_nProton_above260, MC_nProton_below260, MC_nupdg, MC_ccnc, MC_beamNeutrino, MC_FV);
          // [CC0pi0p] Store true kinematic information for efficiency study
          if(TopologyType == 1){
      
            MC_0pi0p_muon_mom = MC_muon_true_Mom.front();
            MC_0pi0p_muon_theta = MC_muon_true_theta.front();
            MC_0pi0p_muon_costheta = MC_muon_true_cos_theta.front();
            MC_0pi0p_muon_phi = MC_muon_true_phi.front();
      
          }
      
          if(TopologyType == 2){
      
            MC_0pi1p_muon_mom = MC_muon_true_Mom.front();
            MC_0pi1p_muon_theta = MC_muon_true_theta.front();
            MC_0pi1p_muon_costheta = MC_muon_true_cos_theta.front();
            MC_0pi1p_muon_phi = MC_muon_true_phi.front();
      
            MC_0pi1p_proton_mom = MC_proton_true_Mom.front();
            MC_0pi1p_proton_theta = MC_proton_true_theta.front();
            MC_0pi1p_proton_costheta = MC_proton_true_cos_theta.front();
            MC_0pi1p_proton_phi = MC_proton_true_phi.front();
      
            TVector3 muon_dir(MC_muon_true_Px.front(), MC_muon_true_Py.front(), MC_muon_true_Pz.front());
            TVector3 proton_dir(MC_proton_true_Px.front(), MC_proton_true_Py.front(), MC_proton_true_Pz.front());
            MC_0pi1p_cos_ang_muon_proton = (muon_dir * proton_dir) / (muon_dir.Mag() * proton_dir.Mag());
      
            MC_0pi1p_PT = sqrt((MC_muon_true_Px.front() + MC_proton_true_Px.front()) * (MC_muon_true_Px.front() + MC_proton_true_Px.front()) + (MC_muon_true_Py.front() + MC_proton_true_Py.front()) * (MC_muon_true_Py.front() + MC_proton_true_Py.front()));
            MC_0pi1p_PL = abs(MC_muon_true_Pz.front() + MC_proton_true_Pz.front());
          }
      
          if(MC_nMuon == 1){
            MC_transfer_E = MC_nu_E - MC_muon_true_E.front();
          }
          

          my_event_->Fill();
        } // i_mc != matched_MCid (this MC hasn't been filled yet)
      } // loop MCTruth
    } // if MCTruth is > 1
  } // IsMC

/////////////////////////////////////
// Set all the variable to default here
  if(IsMC){
    MaCCQE.clear();
    CoulombCCQE.clear();
    MaNCEL.clear();
    EtaNCEL.clear();

    NormCCMEC.clear();
    NormNCMEC.clear();
    FracPN_CCMEC.clear();
    FracDelta_CCMEC.clear();

    MaCCRES.clear();
    MvCCRES.clear();
    MaNCRES.clear();
    MvNCRES.clear();

    NonRESBGvpCC1pi.clear();
    NonRESBGvpCC2pi.clear();
    NonRESBGvpNC1pi.clear();
    NonRESBGvpNC2pi.clear();
    NonRESBGvnCC1pi.clear();
    NonRESBGvnCC2pi.clear();
    NonRESBGvnNC1pi.clear();
    NonRESBGvnNC2pi.clear();
    NonRESBGvbarpCC1pi.clear();
    NonRESBGvbarpCC2pi.clear();
    NonRESBGvbarpNC1pi.clear();
    NonRESBGvbarpNC2pi.clear();
    NonRESBGvbarnCC1pi.clear();
    NonRESBGvbarnCC2pi.clear();
    NonRESBGvbarnNC1pi.clear();
    NonRESBGvbarnNC2pi.clear();
    AhtBY.clear();
    BhtBY.clear();
    CV1uBY.clear();
    CV2uBY.clear();

    AGKYxF1pi.clear();
    AGKYpT1pi.clear();

    MFP_pi.clear();
    MFP_N.clear();
    FrCEx_pi.clear();
    FrInel_pi.clear();
    FrAbs_pi.clear();
    FrCEx_N.clear();
    FrInel_N.clear();
    FrAbs_N.clear();

    RDecBR1gamma.clear();
    RDecBR1eta.clear();

    RootinoFix.clear();
    splines_general_Spline.clear();

    MaCCQE_tune.clear();
    RPA_CCQE_tune.clear();
    NormCCMEC_tune.clear();
    XSecShape_CCMEC_tune.clear();

    piplusReacLow.clear();
    piplusReacHigh.clear();
    piplusAbs.clear();
    piplusCex.clear();
    piplusDCex.clear();
    piplusPiProd.clear();
    piplusElast.clear();

    piminusReacLow.clear();
    piminusReacHigh.clear();
    piminusAbs.clear();
    piminusCex.clear();
    piminusDCex.clear();
    piminusPiProd.clear();
    piminusElast.clear();

    protonReac.clear();
    protonElast.clear();

    neutronAbs_max.clear();
    neutronElast_max.clear();
    neutronInela_max.clear();
    neutronPiprod_max.clear();
    neutronPprod_max.clear();

    neutronAbs_100.clear();
    neutronElast_100.clear();
    neutronInela_100.clear();
    neutronPiprod_100.clear();
    neutronPprod_100.clear();

    neutronAbs_200.clear();
    neutronElast_200.clear();
    neutronInela_200.clear();
    neutronPiprod_200.clear();
    neutronPprod_200.clear();

    expskin.clear();
    horncurrent.clear();
    nucleoninexsec.clear();
    nucleonqexsec.clear();
    nucleontotxsec.clear();
    pioninexsec.clear();
    pionqexsec.clear();
    piontotxsec.clear();

    nu_parent_pdg = -999;
    MC_beamNeutrino = false;
    MC_nupdg = -999;
    MC_ccnc = -999;
    MC_Q2 = -999;
    MC_FV = false;
    MC_if_in_active = false;
    MC_int_mode = -999;
    MC_nu_E = -999;
    MC_transfer_E = -999;
    MC_nuVtxX = -999;
    MC_nuVtxY = -999;
    MC_nuVtxZ = -999;

    MC_nMuon = 0;
    MC_nElectron = 0;
    MC_nNeutron = 0;
    MC_nProton_below260 = 0;
    MC_nProton_above260 = 0;
    MC_nPi0 = 0;
    MC_nPiPlus_below80 = 0;
    MC_nPiPlus_above80 = 0;
    MC_nPiMinus_below80 = 0;
    MC_nPiMinus_above80 = 0;

    MC_Primary_PDG.clear();
    MC_Primary_Mom.clear();
    Ghost_PDG.clear();

    TopologyType = -999;

    MC_muon_true_Mom.clear();
    MC_muon_true_theta.clear();
    MC_muon_true_cos_theta.clear();
    MC_muon_true_phi.clear();
    MC_muon_true_Px.clear();
    MC_muon_true_Py.clear();
    MC_muon_true_Pz.clear();
    MC_muon_true_E.clear();

    MC_proton_true_Mom.clear();
    MC_proton_true_theta.clear();
    MC_proton_true_cos_theta.clear();
    MC_proton_true_phi.clear();
    MC_proton_true_Px.clear();
    MC_proton_true_Py.clear();
    MC_proton_true_Pz.clear();

    MC_0pi0p_muon_mom = -999;
    MC_0pi0p_muon_theta = -999;
    MC_0pi0p_muon_costheta = -999;
    MC_0pi0p_muon_phi = -999;

    MC_0pi1p_muon_mom = -999;
    MC_0pi1p_muon_theta = -999;
    MC_0pi1p_muon_costheta = -999;
    MC_0pi1p_muon_phi = -999;

    MC_0pi1p_proton_mom = -999;
    MC_0pi1p_proton_theta = -999;
    MC_0pi1p_proton_costheta = -999;
    MC_0pi1p_proton_phi = -999;

    MC_0pi1p_cos_ang_muon_proton = -999;
    MC_0pi1p_PT = -999;
    MC_0pi1p_PL = -999;

    if_MC_reco_match = {false, false};
    matched_MCid = 0;
    matched_MCid_2trk = {-999, -999};
    MCmatch_type = 0;
//////////////////
    trk_cosmic_percent = {-999, -999};
    trk_purity = {-999, -999};
    trk_completeness = {-999, -999};

    if_matchMu = {false, false};
    if_matchP = {false, false};
    if_matchPrimary = {false, false};
    if_cosmic = {true, true};

    true_mom = {-999, -999};
    true_start_x = {-999, -999};
    true_start_y = {-999, -999};
    true_start_z = {-999, -999};
    true_end_x = {-999, -999};
    true_end_y = {-999, -999};
    true_end_z = {-999, -999};
    true_trk_ifcontained = {false, false};
    true_vtxFV = {false, false};
    true_trk_phi = {-999, -999};
    true_trk_theta = {-999, -999};
    true_trk_costheta = {-999, -999};
    true_trk_theta_yz = {-999, -999};
    true_trk_costheta_yz = {-999, -999};
    true_trk_theta_xz = {-999, -999};
    true_trk_costheta_xz = {-999, -999};
    true_trk_length = {-999, -999};
    true_trk_PDG = {-999, -999};
  }

  Ghost_PDG.clear();

  flash_matching_chi2 = -999;
  trigger_time = -999;
  flash_YCenter = -999;
  flash_YWidth = -999;
  flash_ZCenter = -999;
  flash_ZWidth = -999;
  flash_TotalPE = -999;
  flash_time = -999;

  if_2tracks = false;

  crthit_PE.clear();
  crthit_plane.clear();
  crthit_time.clear();
  Nr_crthit_inBeam = 0;
  only_crthit_time_inBeam = -999;
  if_t0_crt_time_inBeam_match = true;

  evt_CRTveto_70 = false;
  evt_CRTveto_100 = false;

  Nr_trk_asso_crthit = {0, 0};
  trk_crt_time = {-999, -999};
  if_trk_CRT_out_Beam = {false, false};
  if_t0_trk_crt_time_match = {true, true};

  trk_OutOfTime = {false, false};

  trk_start_x = {-999, -999};
  trk_start_y = {-999, -999};
  trk_start_z = {-999, -999};
  trk_end_x = {-999, -999};
  trk_end_y = {-999, -999};
  trk_end_z = {-999, -999};

  trk_start_InFV = {false, false};
  trk_contained = {false, false};

  trk_phi = {-999, -999};
  trk_theta = {-999, -999};
  trk_costheta = {-999, -999};

  theta_pl0 = {-999, -999};
  theta_pl1 = {-999, -999};
  theta_pl2 = {-999, -999};
  phi_readout = {-999, -999};

  sin2_theta_pl0 = {-999, -999};
  sin2_theta_pl1 = {-999, -999};
  sin2_theta_pl2 = {-999, -999};
  sin2_phi_readout = {-999, -999};

  bestMCS = {-999, -999};
  bestMCSLL = {-999, -999};
  fwdMCS = {-999, -999};
  fwdMCSLL = {-999, -999};
  bwdMCS = {-999, -999};
  bwdMCSLL = {-999, -999};

  bestMCS_NoSCE = {-999, -999};
  fwdMCS_NoSCE = {-999, -999};
  bwdMCS_NoSCE = {-999, -999};
  bestMCSLL_NoSCE = {-999, -999};
  fwdMCSLL_NoSCE = {-999, -999};
  bwdMCSLL_NoSCE = {-999, -999};

  Trk_length = {-999, -999};
  Mom_Range_mu = {-999, -999};
  Mom_Range_p = {-999, -999};
  Mom_Range_pi = {-999, -999};

  dEdx_pl0.clear();
  dEdx_pl1.clear();
  dEdx_pl2.clear();
  resRange_pl0.clear();
  resRange_pl1.clear();
  resRange_pl2.clear();

  dEdx_pl0_start_half = {-999, -999};
  dEdx_pl1_start_half = {-999, -999};
  dEdx_pl2_start_half = {-999, -999};
  dEdx_pl0_end_half = {-999, -999};
  dEdx_pl1_end_half = {-999, -999};
  dEdx_pl2_end_half = {-999, -999};

  dEdx_pl0_start1020 = {-999, -999};
  dEdx_pl1_start1020 = {-999, -999};
  dEdx_pl2_start1020 = {-999, -999};
  dEdx_pl0_end1020 = {-999, -999};
  dEdx_pl1_end1020 = {-999, -999};
  dEdx_pl2_end1020 = {-999, -999};

  PID_Chi2Mu_3pl = {-999, -999};
  PID_Chi2P_3pl = {-999, -999};
  PID_Chi2Pi_3pl = {-999, -999};
  PID_Chi2K_3pl = {-999, -999};

  muon_idx = -999;
  proton_idx = -999;

  vtx_x = -999;
  vtx_y = -999;
  vtx_z = -999;
  vtx_InFV = false;

  muon_cand_phi = -999;
  muon_cand_theta = -999;
  muon_cand_costheta = -999;
  muon_cand_mom_MCS = -999;
  muon_cand_mom_range = -999;
  muon_cand_length = -999;

  proton_cand_phi = -999;
  proton_cand_theta = -999;
  proton_cand_costheta = -999;
  proton_cand_mom_range = -999;
  proton_cand_length = -999;

  cos_ang_muon_proton = -999;

  dist_mu_p_start = -999;

  PT = -999;
  PL = -999;

  if_match_mu_p = false;
  if_match_mu_p_flipped = false;

  reco_MC_dist_vtx = -999;
  reco_MC_dist_vtx_noSCE = -999;
}

void SingleMuon::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;

  // Make a tree store run and subrun info for POT
  POTtree = tfs->make<TTree>("POTtree","POTtree");

  POTtree->Branch("run", &run);
  POTtree->Branch("subrun", &subrun);
  POTtree->Branch("POT_miss1E10", &POT_miss1E10);

  // Make a tree to store selection information
  my_event_ = tfs->make<TTree>("tree","tree");

  my_event_->Branch("evt_run", &evt_run);
  my_event_->Branch("evt_subrun", &evt_subrun);
  my_event_->Branch("evt_evt", &evt_evt);

  if(IsMC){
    my_event_->Branch("Nr_MCNu", &Nr_MCNu);

    my_event_->Branch("MaCCQE", &MaCCQE);
    my_event_->Branch("CoulombCCQE", &CoulombCCQE);
    my_event_->Branch("MaNCEL", &MaNCEL);
    my_event_->Branch("EtaNCEL", &EtaNCEL);

    my_event_->Branch("NormCCMEC", &NormCCMEC);
    my_event_->Branch("NormNCMEC", &NormNCMEC);
    my_event_->Branch("FracPN_CCMEC", &FracPN_CCMEC);
    my_event_->Branch("FracDelta_CCMEC", &FracDelta_CCMEC);

    my_event_->Branch("MaCCRES", &MaCCRES);
    my_event_->Branch("MvCCRES", &MvCCRES);
    my_event_->Branch("MaNCRES", &MaNCRES);
    my_event_->Branch("MvNCRES", &MvNCRES);

    my_event_->Branch("NonRESBGvpCC1pi", &NonRESBGvpCC1pi);
    my_event_->Branch("NonRESBGvpCC2pi", &NonRESBGvpCC2pi);
    my_event_->Branch("NonRESBGvpNC1pi", &NonRESBGvpNC1pi);
    my_event_->Branch("NonRESBGvpNC2pi", &NonRESBGvpNC2pi);
    my_event_->Branch("NonRESBGvnCC1pi", &NonRESBGvnCC1pi);
    my_event_->Branch("NonRESBGvnCC2pi", &NonRESBGvnCC2pi);
    my_event_->Branch("NonRESBGvnNC1pi", &NonRESBGvnNC1pi);
    my_event_->Branch("NonRESBGvnNC2pi", &NonRESBGvnNC2pi);
    my_event_->Branch("NonRESBGvbarpCC1pi", &NonRESBGvbarpCC1pi);
    my_event_->Branch("NonRESBGvbarpCC2pi", &NonRESBGvbarpCC2pi);
    my_event_->Branch("NonRESBGvbarpNC1pi", &NonRESBGvbarpNC1pi);
    my_event_->Branch("NonRESBGvbarpNC2pi", &NonRESBGvbarpNC2pi);
    my_event_->Branch("NonRESBGvbarnCC1pi", &NonRESBGvbarnCC1pi);
    my_event_->Branch("NonRESBGvbarnCC2pi", &NonRESBGvbarnCC2pi);
    my_event_->Branch("NonRESBGvbarnNC1pi", &NonRESBGvbarnNC1pi);
    my_event_->Branch("NonRESBGvbarnNC2pi", &NonRESBGvbarnNC2pi);
    my_event_->Branch("AhtBY", &AhtBY);
    my_event_->Branch("BhtBY", &BhtBY);
    my_event_->Branch("CV1uBY", &CV1uBY);
    my_event_->Branch("CV2uBY", &CV2uBY);

    my_event_->Branch("AGKYxF1pi", &AGKYxF1pi);
    my_event_->Branch("AGKYpT1pi", &AGKYpT1pi);

    my_event_->Branch("MFP_pi", &MFP_pi);
    my_event_->Branch("MFP_N", &MFP_N);
    my_event_->Branch("FrCEx_pi", &FrCEx_pi);
    my_event_->Branch("FrInel_pi", &FrInel_pi);
    my_event_->Branch("FrAbs_pi", &FrAbs_pi);
    my_event_->Branch("FrCEx_N", &FrCEx_N);
    my_event_->Branch("FrInel_N", &FrInel_N);
    my_event_->Branch("FrAbs_N", &FrAbs_N);

    my_event_->Branch("RDecBR1gamma", &RDecBR1gamma);
    my_event_->Branch("RDecBR1eta", &RDecBR1eta);

    my_event_->Branch("RootinoFix", &RootinoFix);
    my_event_->Branch("splines_general_Spline", &splines_general_Spline);

    my_event_->Branch("MaCCQE_tune", &MaCCQE_tune);
    my_event_->Branch("RPA_CCQE_tune", &RPA_CCQE_tune);
    my_event_->Branch("NormCCMEC_tune", &NormCCMEC_tune);
    my_event_->Branch("XSecShape_CCMEC_tune", &XSecShape_CCMEC_tune);

    my_event_->Branch("piplusReacLow", &piplusReacLow);
    my_event_->Branch("piplusReacHigh", &piplusReacHigh);
    my_event_->Branch("piplusAbs", &piplusAbs);
    my_event_->Branch("piplusCex", &piplusCex);
    my_event_->Branch("piplusDCex", &piplusDCex);
    my_event_->Branch("piplusPiProd", &piplusPiProd);
    my_event_->Branch("piplusElast", &piplusElast);

    my_event_->Branch("piminusReacLow", &piminusReacLow);
    my_event_->Branch("piminusReacHigh", &piminusReacHigh);
    my_event_->Branch("piminusAbs", &piminusAbs);
    my_event_->Branch("piminusCex", &piminusCex);
    my_event_->Branch("piminusDCex", &piminusDCex);
    my_event_->Branch("piminusPiProd", &piminusPiProd);
    my_event_->Branch("piminusElast", &piminusElast);

    my_event_->Branch("protonReac", &protonReac);
    my_event_->Branch("protonElast", &protonElast);

    my_event_->Branch("neutronAbs_max", &neutronAbs_max);
    my_event_->Branch("neutronElast_max", &neutronElast_max);
    my_event_->Branch("neutronInela_max", &neutronInela_max);
    my_event_->Branch("neutronPiprod_max", &neutronPiprod_max);
    my_event_->Branch("neutronPprod_max", &neutronPprod_max);

    my_event_->Branch("neutronAbs_100", &neutronAbs_100);
    my_event_->Branch("neutronElast_100", &neutronElast_100);
    my_event_->Branch("neutronInela_100", &neutronInela_100);
    my_event_->Branch("neutronPiprod_100", &neutronPiprod_100);
    my_event_->Branch("neutronPprod_100", &neutronPprod_100);

    my_event_->Branch("neutronAbs_200", &neutronAbs_200);
    my_event_->Branch("neutronElast_200", &neutronElast_200);
    my_event_->Branch("neutronInela_200", &neutronInela_200);
    my_event_->Branch("neutronPiprod_200", &neutronPiprod_200);
    my_event_->Branch("neutronPprod_200", &neutronPprod_200);

    my_event_->Branch("expskin" ,&expskin);
    my_event_->Branch("horncurrent", &horncurrent);
    my_event_->Branch("nucleoninexsec", &nucleoninexsec);
    my_event_->Branch("nucleonqexsec", &nucleonqexsec);
    my_event_->Branch("nucleontotxsec", &nucleontotxsec);
    my_event_->Branch("pioninexsec", &pioninexsec);
    my_event_->Branch("pionqexsec", &pionqexsec);
    my_event_->Branch("piontotxsec", &piontotxsec);

    my_event_->Branch("EventWeight", &EventWeight);
    my_event_->Branch("TopologyType", &TopologyType);
    my_event_->Branch("nu_parent_pdg", &nu_parent_pdg);
    my_event_->Branch("MC_beamNeutrino", &MC_beamNeutrino);
    my_event_->Branch("MC_int_mode", &MC_int_mode);
    my_event_->Branch("MC_nupdg", &MC_nupdg);
    my_event_->Branch("MC_ccnc", &MC_ccnc);
    my_event_->Branch("MC_Q2", &MC_Q2);
    my_event_->Branch("MC_nu_E", &MC_nu_E);
    my_event_->Branch("MC_transfer_E", &MC_transfer_E);
    my_event_->Branch("MC_nuVtxX", &MC_nuVtxX);
    my_event_->Branch("MC_nuVtxY", &MC_nuVtxY);
    my_event_->Branch("MC_nuVtxZ", &MC_nuVtxZ);
    my_event_->Branch("MC_FV", &MC_FV);
    my_event_->Branch("MC_if_in_active", &MC_if_in_active);
    my_event_->Branch("MC_nMuon", &MC_nMuon);
    my_event_->Branch("MC_nElectron", &MC_nElectron);
    my_event_->Branch("MC_nNeutron", &MC_nNeutron);
    my_event_->Branch("MC_nProton_below260", &MC_nProton_below260);
    my_event_->Branch("MC_nProton_above260", &MC_nProton_above260);
    my_event_->Branch("MC_nPi0", &MC_nPi0);
    my_event_->Branch("MC_nPiPlus_below80", &MC_nPiPlus_below80);
    my_event_->Branch("MC_nPiPlus_above80", &MC_nPiPlus_above80);
    my_event_->Branch("MC_nPiMinus_below80", &MC_nPiMinus_below80);
    my_event_->Branch("MC_nPiMinus_above80", &MC_nPiMinus_above80);
    my_event_->Branch("MC_Primary_PDG", &MC_Primary_PDG);
    my_event_->Branch("MC_Primary_Mom", &MC_Primary_Mom);

    my_event_->Branch("MC_0pi0p_muon_mom", &MC_0pi0p_muon_mom);
    my_event_->Branch("MC_0pi0p_muon_theta", &MC_0pi0p_muon_theta);
    my_event_->Branch("MC_0pi0p_muon_costheta", &MC_0pi0p_muon_costheta);
    my_event_->Branch("MC_0pi0p_muon_phi", &MC_0pi0p_muon_phi);

    my_event_->Branch("MC_0pi1p_muon_mom", &MC_0pi1p_muon_mom);
    my_event_->Branch("MC_0pi1p_muon_theta", &MC_0pi1p_muon_theta);
    my_event_->Branch("MC_0pi1p_muon_costheta", &MC_0pi1p_muon_costheta);
    my_event_->Branch("MC_0pi1p_muon_phi", &MC_0pi1p_muon_phi);
    my_event_->Branch("MC_0pi1p_proton_mom", &MC_0pi1p_proton_mom);
    my_event_->Branch("MC_0pi1p_proton_theta", &MC_0pi1p_proton_theta);
    my_event_->Branch("MC_0pi1p_proton_costheta", &MC_0pi1p_proton_costheta);
    my_event_->Branch("MC_0pi1p_proton_phi", &MC_0pi1p_proton_phi);
    my_event_->Branch("MC_0pi1p_cos_ang_muon_proton", &MC_0pi1p_cos_ang_muon_proton);
    my_event_->Branch("MC_0pi1p_PT", &MC_0pi1p_PT);
    my_event_->Branch("MC_0pi1p_PL", &MC_0pi1p_PL);

    my_event_->Branch("MCmatch_type", &MCmatch_type);

///////// selected
    my_event_->Branch("Ghost_PDG", &Ghost_PDG);

    my_event_->Branch("true_mom", &true_mom);
    my_event_->Branch("true_start_x", &true_start_x);
    my_event_->Branch("true_start_y", &true_start_y);
    my_event_->Branch("true_start_z", &true_start_z);
    my_event_->Branch("true_end_x", &true_end_x);
    my_event_->Branch("true_end_y", &true_end_y);
    my_event_->Branch("true_end_z", &true_end_z);
    my_event_->Branch("true_trk_phi", &true_trk_phi);
    my_event_->Branch("true_trk_theta", &true_trk_theta);
    my_event_->Branch("true_trk_costheta", &true_trk_costheta);
    my_event_->Branch("true_trk_theta_yz", &true_trk_theta_yz);
    my_event_->Branch("true_trk_costheta_yz", &true_trk_costheta_yz);
    my_event_->Branch("true_trk_theta_xz", &true_trk_theta_xz);
    my_event_->Branch("true_trk_costheta_xz", &true_trk_costheta_xz);
    my_event_->Branch("true_trk_length", &true_trk_length);
    my_event_->Branch("true_trk_PDG", &true_trk_PDG);
    my_event_->Branch("true_trk_ifcontained", &true_trk_ifcontained);
    my_event_->Branch("true_vtxFV", &true_vtxFV);

    my_event_->Branch("trk_cosmic_percent", &trk_cosmic_percent);
    my_event_->Branch("trk_purity", &trk_purity);
    my_event_->Branch("trk_completeness", &trk_completeness);

    my_event_->Branch("if_cosmic", &if_cosmic);
    my_event_->Branch("if_matchMu", &if_matchMu);
    my_event_->Branch("if_matchP", &if_matchP);
    my_event_->Branch("if_matchPrimary", &if_matchPrimary);

    my_event_->Branch("reco_MC_dist_vtx", &reco_MC_dist_vtx);
    my_event_->Branch("reco_MC_dist_vtx_noSCE", &reco_MC_dist_vtx_noSCE);

    my_event_->Branch("if_match_mu_p", &if_match_mu_p);
    my_event_->Branch("if_match_mu_p_flipped", &if_match_mu_p_flipped);
  }

  my_event_->Branch("flash_matching_chi2", &flash_matching_chi2);
  my_event_->Branch("trigger_time", &trigger_time);
  my_event_->Branch("flash_YCenter", &flash_YCenter);
  my_event_->Branch("flash_YWidth", &flash_YWidth);
  my_event_->Branch("flash_ZCenter", &flash_ZCenter);
  my_event_->Branch("flash_ZWidth", &flash_ZWidth);
  my_event_->Branch("flash_TotalPE", &flash_TotalPE);
  my_event_->Branch("flash_time", &flash_time);

  my_event_->Branch("n_pfp_nuDaughters", &n_pfp_nuDaughters);
  my_event_->Branch("n_dau_tracks", &n_dau_tracks);
  my_event_->Branch("n_dau_showers", &n_dau_showers);

  my_event_->Branch("if_2tracks", &if_2tracks);

  my_event_->Branch("crthit_PE", &crthit_PE);
  my_event_->Branch("crthit_plane", &crthit_plane);
  my_event_->Branch("crthit_time", &crthit_time);
  my_event_->Branch("Nr_crthit_inBeam", &Nr_crthit_inBeam);
  my_event_->Branch("only_crthit_time_inBeam", &only_crthit_time_inBeam);
  my_event_->Branch("if_t0_crt_time_inBeam_match", &if_t0_crt_time_inBeam_match);

  my_event_->Branch("evt_CRTveto_70", &evt_CRTveto_70);
  my_event_->Branch("evt_CRTveto_100", &evt_CRTveto_100);

  my_event_->Branch("Nr_trk_asso_crthit", &Nr_trk_asso_crthit);
  my_event_->Branch("trk_crt_time", &trk_crt_time);
  my_event_->Branch("if_trk_CRT_out_Beam", &if_trk_CRT_out_Beam);
  my_event_->Branch("if_t0_trk_crt_time_match", &if_t0_trk_crt_time_match);

  my_event_->Branch("trk_OutOfTime", &trk_OutOfTime);

  my_event_->Branch("trk_start_x", &trk_start_x);
  my_event_->Branch("trk_start_y", &trk_start_y);
  my_event_->Branch("trk_start_z", &trk_start_z);
  my_event_->Branch("trk_end_x", &trk_end_x);
  my_event_->Branch("trk_end_y", &trk_end_y);
  my_event_->Branch("trk_end_z", &trk_end_z);

  my_event_->Branch("trk_start_InFV", &trk_start_InFV);
  my_event_->Branch("trk_contained", &trk_contained);

  my_event_->Branch("trk_phi", &trk_phi);
  my_event_->Branch("trk_theta", &trk_theta);
  my_event_->Branch("trk_costheta", &trk_costheta);

  my_event_->Branch("sin2_theta_pl2", &sin2_theta_pl2);
  my_event_->Branch("sin2_theta_pl1", &sin2_theta_pl1);
  my_event_->Branch("sin2_theta_pl0", &sin2_theta_pl0);
  my_event_->Branch("sin2_phi_readout", &sin2_phi_readout);

  my_event_->Branch("bestMCS", &bestMCS);
  my_event_->Branch("bestMCSLL", &bestMCSLL);
  my_event_->Branch("fwdMCS", &fwdMCS);
  my_event_->Branch("fwdMCSLL", &fwdMCSLL);
  my_event_->Branch("bwdMCS", &bwdMCS);
  my_event_->Branch("bwdMCSLL", &bwdMCSLL);
  my_event_->Branch("bestMCS_NoSCE", &bestMCS_NoSCE);
  my_event_->Branch("fwdMCS_NoSCE", &fwdMCS_NoSCE);
  my_event_->Branch("bwdMCS_NoSCE", &bwdMCS_NoSCE);
  my_event_->Branch("bestMCSLL_NoSCE", &bestMCSLL_NoSCE);
  my_event_->Branch("fwdMCSLL_NoSCE", &fwdMCSLL_NoSCE);
  my_event_->Branch("bwdMCSLL_NoSCE", &bwdMCSLL_NoSCE);

  my_event_->Branch("Trk_length", &Trk_length);
  my_event_->Branch("Mom_Range_mu", &Mom_Range_mu);
  my_event_->Branch("Mom_Range_p", &Mom_Range_p);
  my_event_->Branch("Mom_Range_pi", &Mom_Range_pi);

  my_event_->Branch("dEdx_pl0_start_half", &dEdx_pl0_start_half);
  my_event_->Branch("dEdx_pl1_start_half", &dEdx_pl1_start_half);
  my_event_->Branch("dEdx_pl2_start_half", &dEdx_pl2_start_half);
  my_event_->Branch("dEdx_pl0_end_half", &dEdx_pl0_end_half);
  my_event_->Branch("dEdx_pl1_end_half", &dEdx_pl1_end_half);
  my_event_->Branch("dEdx_pl2_end_half", &dEdx_pl2_end_half);

  my_event_->Branch("dEdx_pl0_start1020", &dEdx_pl0_start1020);
  my_event_->Branch("dEdx_pl1_start1020", &dEdx_pl1_start1020);
  my_event_->Branch("dEdx_pl2_start1020", &dEdx_pl2_start1020);
  my_event_->Branch("dEdx_pl0_end1020", &dEdx_pl0_end1020);
  my_event_->Branch("dEdx_pl1_end1020", &dEdx_pl1_end1020);
  my_event_->Branch("dEdx_pl2_end1020", &dEdx_pl2_end1020);

  my_event_->Branch("dEdx_pl2_1020_ratio", &dEdx_pl2_1020_ratio);
  my_event_->Branch("dEdx_pl2_half_ratio", &dEdx_pl2_half_ratio);

  my_event_->Branch("PID_Chi2Mu_3pl", &PID_Chi2Mu_3pl);
  my_event_->Branch("PID_Chi2P_3pl", &PID_Chi2P_3pl);
  my_event_->Branch("PID_Chi2Pi_3pl", &PID_Chi2Pi_3pl);
  my_event_->Branch("PID_Chi2K_3pl", &PID_Chi2K_3pl);

  my_event_->Branch("muon_idx", &muon_idx);
  my_event_->Branch("proton_idx", &proton_idx);

  my_event_->Branch("vtx_x", &vtx_x);
  my_event_->Branch("vtx_y", &vtx_y);
  my_event_->Branch("vtx_z", &vtx_z);
  my_event_->Branch("vtx_InFV", &vtx_InFV);

  my_event_->Branch("muon_cand_phi", &muon_cand_phi);
  my_event_->Branch("muon_cand_theta", &muon_cand_theta);
  my_event_->Branch("muon_cand_costheta", &muon_cand_costheta);
  my_event_->Branch("muon_cand_mom_MCS", &muon_cand_mom_MCS);
  my_event_->Branch("muon_cand_mom_range", &muon_cand_mom_range);
  my_event_->Branch("muon_cand_length", &muon_cand_length);

  my_event_->Branch("proton_cand_phi", &proton_cand_phi);
  my_event_->Branch("proton_cand_theta", &proton_cand_theta);
  my_event_->Branch("proton_cand_costheta", &proton_cand_costheta);
  my_event_->Branch("proton_cand_mom_range", &proton_cand_mom_range);
  my_event_->Branch("proton_cand_length", &proton_cand_length);

  my_event_->Branch("cos_ang_muon_proton", &cos_ang_muon_proton);

  my_event_->Branch("dist_mu_p_start", &dist_mu_p_start);

  my_event_->Branch("PT", &PT);
  my_event_->Branch("PL", &PL);
}

void SingleMuon::endSubRun(art::SubRun const &sr){

  run       = sr.run();
  subrun    = sr.subRun();

  // For Overlay or MC
  art::Handle<sumdata::POTSummary> Handle_potsum;
  if(sr.getByLabel("generator", Handle_potsum)){
    POT_miss1E10 = Handle_potsum->totpot;
    POT_miss1E10 = POT_miss1E10 / 1E10;
  }

  POTtree->Fill();
}


void SingleMuon::beginJob()
{
  // Implementation of optional member function here.
  Initialize_event();
}

void SingleMuon::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SingleMuon)

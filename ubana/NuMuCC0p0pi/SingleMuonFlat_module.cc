////////////////////////////////////////////////////////////////////////
//// Class:       SingleMuon
//// Plugin Type: analyzer (art v3_01_02)
//// File:        SingleMuon_module.cc
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
#include "nusimdata/SimulationBase/GTruth.h"

#include "TTree.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"

#include "FiducialVolume.h"
#include "BackTrackerTruthMatch.h"
#include "RecoTruthMCParticle.h"
#include "BrokenTrack.h"
#include "Topology.h"
#include "PID.h"

class SingleMuon;

class SingleMuon : public art::EDAnalyzer {
public:
  explicit SingleMuon(fhicl::ParameterSet const& pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use:


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
  std::vector<std::vector<double>> MaCCQE;
  std::vector<std::vector<double>> CoulombCCQE;
  std::vector<std::vector<double>> MaNCEL;
  std::vector<std::vector<double>> EtaNCEL;

  std::vector<std::vector<double>> NormCCMEC;
  std::vector<std::vector<double>> NormNCMEC;
  std::vector<std::vector<double>> FracPN_CCMEC;
  std::vector<std::vector<double>> FracDelta_CCMEC;

  std::vector<std::vector<double>> MaCCRES;
  std::vector<std::vector<double>> MvCCRES;
  std::vector<std::vector<double>> MaNCRES;
  std::vector<std::vector<double>> MvNCRES;

  std::vector<std::vector<double>> NonRESBGvpCC1pi;
  std::vector<std::vector<double>> NonRESBGvpCC2pi;
  std::vector<std::vector<double>> NonRESBGvpNC1pi;
  std::vector<std::vector<double>> NonRESBGvpNC2pi;
  std::vector<std::vector<double>> NonRESBGvnCC1pi;
  std::vector<std::vector<double>> NonRESBGvnCC2pi;
  std::vector<std::vector<double>> NonRESBGvnNC1pi;
  std::vector<std::vector<double>> NonRESBGvnNC2pi;
  std::vector<std::vector<double>> NonRESBGvbarpCC1pi;
  std::vector<std::vector<double>> NonRESBGvbarpCC2pi;
  std::vector<std::vector<double>> NonRESBGvbarpNC1pi;
  std::vector<std::vector<double>> NonRESBGvbarpNC2pi;
  std::vector<std::vector<double>> NonRESBGvbarnCC1pi;
  std::vector<std::vector<double>> NonRESBGvbarnCC2pi;
  std::vector<std::vector<double>> NonRESBGvbarnNC1pi;
  std::vector<std::vector<double>> NonRESBGvbarnNC2pi;
  std::vector<std::vector<double>> AhtBY;
  std::vector<std::vector<double>> BhtBY;
  std::vector<std::vector<double>> CV1uBY;
  std::vector<std::vector<double>> CV2uBY;

  std::vector<std::vector<double>> AGKYxF1pi;
  std::vector<std::vector<double>> AGKYpT1pi;

  std::vector<std::vector<double>> MFP_pi;
  std::vector<std::vector<double>> MFP_N;
  std::vector<std::vector<double>> FrCEx_pi;
  std::vector<std::vector<double>> FrInel_pi;
  std::vector<std::vector<double>> FrAbs_pi;
  std::vector<std::vector<double>> FrCEx_N;
  std::vector<std::vector<double>> FrInel_N;
  std::vector<std::vector<double>> FrAbs_N;

  std::vector<std::vector<double>> RDecBR1gamma;
  std::vector<std::vector<double>> RDecBR1eta;

  double EventWeight = 1; // Spine reweight using Steven's tool

  int Nr_MCNu = 0;

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

  double MC_0pi0p_PT = -999;
  double MC_0pi0p_PL = -999; 

  double MC_0pi1p_muon_mom = -999;
  double MC_0pi1p_muon_theta = -999;
  double MC_0pi1p_muon_costheta = -999;
  double MC_0pi1p_muon_phi = -999;

  double MC_0pi1p_proton_mom = -999;
  double MC_0pi1p_proton_theta = -999;
  double MC_0pi1p_proton_costheta = -999;
  double MC_0pi1p_proton_phi = -999;

  double MC_0pi1p_cos_ang_muon_proton = -999;

  double MC_0pi2p_proton1_mom = -999;
  double MC_0pi2p_proton2_mom = -999;

  bool if_MC_reco_match = false;

  // -- Reco
  int n_pfp_nuDaughters = 0; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks = 0; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers = 0; // number of showers asssociated to pfp neutrino daughters

  int nr_granddau_shw = 0;
  int nr_granddau_trk = 0;
  int nr_granddau = 0;
  std::vector<int> MC_granddau_pdg;
  std::vector<double> granddau_trk_len;
  std::vector<double> granddau_shw_len;

  std::vector<int> Ghost_PDG; // pdg code of the pfp which has no track or shower associated; No elements ideally

  bool if_1track = false; // If selected based on the reco info

  double flash_matching_chi2 = -999; //Chi2 of flash matching in each neutrino slice
  double trigger_time = -999;
  double flash_YCenter = -999;
  double flash_YWidth = -999;
  double flash_ZCenter = -999;
  double flash_ZWidth = -999;
  double flash_TotalPE = -999;
  double flash_time = -999;

  double CRTT0corr = 0;
  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime
  double only_crthit_time_inBeam = -999;
  bool if_t0_crt_time_inBeam_match = true;

  int Nr_trk_asso_crthit = 0; // Number of CRT hits associated to the track
  double trk_crt_time = -999;// The CRT time of the hit which matched to the track
  bool if_trk_CRT_out_Beam = false; // Check if a track matches with out of beam CRT hit(s)
  bool evt_CRTveto_70 = false; // If CRT veto, eliminate the events for contained (70PE threshold)
  bool evt_CRTveto_100 = false; // If CRT veto, eliminate the events for contained (100PE threshold)
  bool evt_CRT_veto_homebrew = false;
  bool if_t0_trk_crt_time_match = true;

  double trk_vtx_x = -999;
  double trk_vtx_y = -999;
  double trk_vtx_z = -999;
  double trk_start_x = -999;
  double trk_start_y = -999;
  double trk_start_z = -999;
  double trk_end_x = -999;
  double trk_end_y = -999;
  double trk_end_z = -999;

  double trk_vtx_x_noSCE = -999;
  double trk_vtx_y_noSCE = -999;
  double trk_vtx_z_noSCE = -999;
  double trk_start_x_noSCE = -999;
  double trk_start_y_noSCE = -999;
  double trk_start_z_noSCE = -999;
  double trk_end_x_noSCE = -999;
  double trk_end_y_noSCE = -999;
  double trk_end_z_noSCE = -999;

  bool trk_OutOfTime = false;

  bool trk_contained = false;
  bool vtx_InFV = false;

  bool if_broken = false;
  double trk_broken_len = -999;
  double trk_broken_nr_merged = -999;
  bool trk_merged_ifcontained = false;
  bool vtx_merged_InFV = false;
  bool if_newTrkThroughGoing = false;

  double sin2_theta_pl0 = -999;
  double sin2_theta_pl1 = -999;
  double sin2_theta_pl2 = -999;
  double sin2_phi_readout = -999;

  double trk_phi = -999;
  double trk_theta = -999;
  double trk_costheta = -999;

  double trk_length = -999;

  double mom_Range_mu = -999;
  double mom_Range_p = -999;
  double mom_Range_pi = -999;

  double bestMCS = -999;
  double bestMCSLL = -999;
  double fwdMCS = -999;
  double fwdMCSLL = -999;
  double bwdMCS = -999;
  double bwdMCSLL = -999;

  double bestMCS_NoSCE = -999;
  double fwdMCS_NoSCE = -999;
  double bwdMCS_NoSCE = -999;
  double bestMCSLL_NoSCE = -999;
  double fwdMCSLL_NoSCE = -999;
  double bwdMCSLL_NoSCE = -999;

  double PT_range = -999;
  double PL_range = -999;
  double PT_MCS = -999;
  double PL_MCS = -999;

  double PID_Chi2Mu_3pl = -999;
  double PID_Chi2P_3pl = -999;
  double PID_Chi2Pi_3pl = -999;
  double PID_Chi2K_3pl = -999;

  std::vector<float> dEdx_pl0;
  std::vector<float> dQdx_pl0;
  std::vector<float> resRange_pl0;

  std::vector<float> dEdx_pl1;
  std::vector<float> dQdx_pl1;
  std::vector<float> resRange_pl1;

  std::vector<float> dEdx_pl2;
  std::vector<float> dQdx_pl2;
  std::vector<float> resRange_pl2;

  double dEdx_pl0_start_half = -999;
  double dEdx_pl1_start_half = -999;
  double dEdx_pl2_start_half = -999;
  double dEdx_pl0_end_half = -999;
  double dEdx_pl1_end_half = -999;
  double dEdx_pl2_end_half = -999;

  double dEdx_pl0_start1020 = -999;
  double dEdx_pl1_start1020 = -999;
  double dEdx_pl2_start1020 = -999;
  double dEdx_pl0_end1020 = -999;
  double dEdx_pl1_end1020 = -999;
  double dEdx_pl2_end1020 = -999;

  double dEdx_pl0_1020_ratio = -999;
  double dEdx_pl0_half_ratio = -999;
  double dEdx_pl1_1020_ratio = -999;
  double dEdx_pl1_half_ratio = -999;
  double dEdx_pl2_1020_ratio = -999;
  double dEdx_pl2_half_ratio = -999;

  double dEdx_pl0_start_half_clean = -999;
  double dEdx_pl1_start_half_clean = -999;
  double dEdx_pl2_start_half_clean = -999;
  double dEdx_pl0_end_half_clean = -999;
  double dEdx_pl1_end_half_clean = -999;
  double dEdx_pl2_end_half_clean = -999;

  double dEdx_pl0_start1020_clean = -999;
  double dEdx_pl1_start1020_clean = -999;
  double dEdx_pl2_start1020_clean = -999;
  double dEdx_pl0_end1020_clean = -999;
  double dEdx_pl1_end1020_clean = -999;
  double dEdx_pl2_end1020_clean = -999;

  double dEdx_pl0_start5_clean = -999;
  double dEdx_pl1_start5_clean = -999;
  double dEdx_pl2_start5_clean = -999;
  double dEdx_pl0_end5_clean = -999;
  double dEdx_pl1_end5_clean = -999;
  double dEdx_pl2_end5_clean = -999;

  double dEdx_pl0_5_ratio_clean = -999;
  double dEdx_pl0_1020_ratio_clean = -999;
  double dEdx_pl0_half_ratio_clean = -999;
  double dEdx_pl1_5_ratio_clean = -999;
  double dEdx_pl1_1020_ratio_clean = -999;
  double dEdx_pl1_half_ratio_clean = -999;
  double dEdx_pl2_5_ratio_clean = -999;
  double dEdx_pl2_1020_ratio_clean = -999;
  double dEdx_pl2_half_ratio_clean = -999;

  // Reco matched truth
  double trk_cosmic_percent = -999;
  double trk_purity = -999;
  double trk_completeness = -999;
  bool if_cosmic = true;
  bool if_matchPrimary = false;
  bool if_matchMu = false;

  double true_mom = -999;//True momentum of muon track in the every event
  double true_start_x = -999;//True start of muon track (X)
  double true_start_y = -999;//True start of muon track (Y)
  double true_start_z = -999;//True start of muon track (Z)
  double true_end_x = -999;//True end of muon track (X)
  double true_end_y = -999;//True end of muon track (Y)
  double true_end_z = -999;//True end of muon track (Z)
  double true_trk_phi = -999;//True phi of muon track
  double true_trk_theta = -999;//True theta of muon track
  double true_trk_costheta = -999;//True cos(theta) of muon track
  double true_trk_theta_yz = -999;
  double true_trk_costheta_yz = -999;
  double true_trk_theta_xz = -999;
  double true_trk_costheta_xz = -999;
  double true_trk_length = -999;//True track length (distance from the start to the end point)
  double true_trk_PDG = -999;//Track pdg
  bool true_trk_ifcontained = -999; // True track if contained or not
  bool true_vtxFV = -999; // True track if contained or not

  double reco_MC_dist_vtx = -999; // Distance of reco - MC vertex w/ SCE correction
  double reco_MC_dist_vtx_noSCE = -999; // Distance of reco - MC vertex w/o SCE correction

  bool                                IsMC;
  bool                                T0Corr;
  bool                                UsingCRT;
  double                              fBeamStart;
  double                              fBeamEnd;
  double                              fDTOffset;
  double                              fDTOffset_overlay;
  std::string                         m_DAQHeaderProducer;
  std::string                         m_generatorLabel;
  std::string                         m_geantLabel;
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

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};

//TODO add run subrun event info
//TODO add weight reader for secondary interactions and flux
//TODO maybe change the reco and MCtruth part into nested class?

SingleMuon::SingleMuon(fhicl::ParameterSet const& pset)
  :
  EDAnalyzer{pset},
  IsMC(pset.get<bool>("IsMC")),
  T0Corr(pset.get<bool>("T0Corr")),
  UsingCRT(pset.get<bool>("UsingCRT")),
  fBeamStart(pset.get<double>("BeamStart")),
  fBeamEnd(pset.get<double>("BeamEnd")),
  fDTOffset(pset.get<double>("DTOffset")),
  fDTOffset_overlay(pset.get<double>("DTOffset_overlay")),
  m_DAQHeaderProducer(pset.get<std::string>("DAQHeaderProducer")),
  m_generatorLabel(pset.get<std::string>("GeneratorLabel")),
  m_geantLabel(pset.get<std::string>("GeantLabel")),
  m_eventweightLabel(pset.get<std::string>("EventweightLabel")),
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
  if (IsMC){
    genie_pars = pset.get< std::vector<std::string> >( "genie_parameter_list" );
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
  evt_run = sr.run();
  evt_subrun = sr.subRun();
  evt_evt = evt.event();

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
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;

  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
  art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;

  std::vector<art::Ptr<evwgh::MCEventWeight> > WeightCollection;
  art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight;

  if(IsMC){
    // MC Truth
    evt.getByLabel(m_generatorLabel, Handle_MCTruth);
    art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

    // MC Particle
    evt.getByLabel(m_geantLabel, Handle_MCParticle);
    art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

    // Genie Reweight
    evt.getByLabel(m_eventweight, Handle_Weight);
    art::fill_ptr_vector(WeightCollection, Handle_Weight);
 
    // Flux Reweight
    // Geant Reweight

    // Neutrino Origin
    const simb::Origin_t Neutrino_Origin = simb::kBeamNeutrino;

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
  art::FindMany<simb::MCParticle,sim::GeneratedParticleInfo> MCtToMCpAsso(Handle_MCTruth, evt, m_geantLabel);


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

  //TODO Think about this again
  //Define the daughters of neutrino
  std::vector<art::Ptr<recob::PFParticle> > NeutrinoDaughters;
  std::vector<art::Ptr<recob::Track> > daughter_Tracks;
  std::vector<art::Ptr<recob::Shower> > daughter_Showers;
  std::vector<int> Track_PDG; // The oder follows daughter_Tracks
  std::vector<int> Shower_PDG; // The oder follows daughter_Showers

  std::vector<float> vTrk_len;

  TVector3 Trk_vtx;
  TVector3 Trk_start;
  TVector3 Trk_end;
  TVector3 Trk_vtx_SCEcorr;
  TVector3 Trk_start_SCEcorr;
  TVector3 Trk_end_SCEcorr;

  art::Ptr< simb::MCParticle > selected_MCparticle;
///////////// to be moved
  //// Flaten the info to be filled per neutrino interactions
  if(IsMC){
    Nr_Iter = MCTruthCollection.size();
  } 
  else {
    Nr_Iter = 1; 
  }
  
 //////////// to be moved 


  //-------- Get Reco neutrino (pfparticle)
  for(unsigned int i = 0; i < pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
    // Check if the pfp is leading a neutrino slice
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      n_pfp_nuDaughters = pfp->NumDaughters();

      // For CC0pi0p, we only consider the case with the number of neutrino daughters less than 4
      if(n_pfp_nuDaughters < 4){
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
      } // If the the nr of daughter pfps is < 4

      //number of tracks and showers
      n_dau_tracks = daughter_Tracks.size();
      n_dau_showers = daughter_Showers.size();

      // Selection and Fill in Info
      if(n_dau_tracks == 1 && n_dau_showers == 0){
        if_1track = true;

        // The daughter of the pfparticle which corresponds to 1 primary track
        // In case there are more than one primary pfparticle when 1 track + 0 shower (primary) precent
        for (unsigned int i_trk_dau = 0; i_trk_dau < NeutrinoDaughters.size(); i_trk_dau++){
          nr_granddau += NeutrinoDaughters[i_trk_dau]->NumDaughters();
          for(int i_granddau = 0; i_granddau < NeutrinoDaughters[i_trk_dau]->NumDaughters(); i_granddau++){
            auto Iterator_granddau = pfParticleIdMap.find(NeutrinoDaughters[i_trk_dau]->Daughters().at(i_granddau));
            auto pfp_granddau = Iterator_granddau->second;

            // Track Association
            auto asso_granddau_track = pfpToTrackAsso.at(pfp_granddau.key());
            if(asso_granddau_track.size()>0){
              nr_granddau_trk += asso_granddau_track.size(); // inclusive number of granddaughters as tracks
              for (unsigned int i_granddau_trk = 0; i_granddau_trk < asso_granddau_track.size(); i_granddau_trk++){
                granddau_trk_len.push_back(asso_granddau_track[i_granddau_trk]->Length());
                if(IsMC){
                  std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(asso_granddau_track[i_granddau_trk].key());
                  BackTrackerTruthMatch backtrackertruthmatch;
                  backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
                  auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
                  if(MCparticle){
                    // MC_granddau_pdg stores all granddaughters of all the primary pfparticle(s)
                    MC_granddau_pdg.push_back(MCparticle->PdgCode());
                  }
                }
              }
            }

            // Shower Association
            auto asso_granddau_shower = pfpToShowerAsso.at(pfp_granddau.key());
            if(asso_granddau_shower.size()>0){
              nr_granddau_shw += asso_granddau_shower.size();// inclusive number of granddaughters as showers
              for (unsigned int i_granddau_shw = 0; i_granddau_shw < asso_granddau_shower.size(); i_granddau_shw++){
                granddau_shw_len.push_back(asso_granddau_shower[i_granddau_shw]->Length());
              }
            }
          }
        } // finish loop over the daughter tracks (under the 1track 0shower condition)


        //-- OpFlash related information
        for(unsigned int i_flash = 0; i_flash < flash_v.size(); i_flash++){
          if(flash_v[i_flash]->Time() >= fBeamStart && flash_v[i_flash]->Time() <= fBeamEnd){
            flash_YCenter = flash_v[i_flash]->YCenter();
            flash_YWidth = flash_v[i_flash]->YWidth();
            flash_ZCenter = flash_v[i_flash]->ZCenter();
            flash_ZWidth = flash_v[i_flash]->ZWidth();
            flash_TotalPE = flash_v[i_flash]->TotalPE();
            flash_time = flash_v[i_flash]->Time();
          }
        }

        auto flash_matching_T0 = pfpToT0Asso.at(pfp.key());
        if(flash_matching_T0.size() == 1){
          flash_matching_chi2 = flash_matching_T0.front()->TriggerConfidence();
          trigger_time = flash_matching_T0.front()->Time();
        }

        //-- CRT related information
        if(UsingCRT){
          double evt_timeGPS_nsec = 0.;
          if(!rawHandle_DAQHeader.isValid()) {
             std::cout << "Could not locate DAQ header." << std::endl;
           }
          raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
          art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
          evt_timeGPS_nsec = evtTimeGPS.timeLow();

          if(T0Corr){
            if(!IsMC && crtT0_v.size() == 1){
              CRTT0corr = crtT0_v.front()->Time();
            }
            else{
              CRTT0corr = 0;
            }
          }

          //- CRT info in beam time
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

            // only count if the CRT hit has signal above 100 PE
            if(crt_time >= fBeamStart && crt_time <= fBeamEnd && crthit_v[i_crt]->peshit > 100){
              Nr_crthit_inBeam++;
              only_crthit_time_inBeam = crt_time;
            } // If CRT hit in Beam window

            //home-brew CRT veto
            if(trigger_time != -999){
              if(abs(crt_time - trigger_time) < 1 && crthit_v[i_crt]->peshit > 100){
                evt_CRT_veto_homebrew = true;
              }
            }

          } // CRT hit loop

          if(Nr_crthit_inBeam == 1){
            if(trigger_time != -999){
              if(abs(only_crthit_time_inBeam - trigger_time) > 1){
                if_t0_crt_time_inBeam_match = false;
              }
            }
          }

          //- If a track is associated to a CRT hit which is outside of beam window, exclude them (For both contained and exiting)
          auto Track_CRThit = CRTToTrackAsso.at(daughter_Tracks.front().key());
          if(Track_CRThit.size() > 0){
            for(unsigned int i_trk_crt = 0; i_trk_crt < Track_CRThit.size(); i_trk_crt++){
              if(Track_CRThit[i_trk_crt]->peshit > 70){

                // w/ T0 correction
                if(T0Corr){
                  trk_crt_time = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
                }
                // w/o T0 correction
                else if(!T0Corr){
                  trk_crt_time = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
                }

                Nr_trk_asso_crthit++;
              }

              if(trk_crt_time < fBeamStart || trk_crt_time > fBeamEnd){
                if_trk_CRT_out_Beam = true;
              } // If matched CRT hit out of Beam window
            }
          }

          //- For contained (Veto if there is any CRT hit within 1us of the flash which is in the beam window)
          if(flash_v.size() > 0){
            for(unsigned int i_fl = 0; i_fl < flash_v.size(); i_fl++){
              auto CRT_hit = CRThitFlashAsso.at(flash_v[i_fl].key());
              if(CRT_hit.size() == 1){
                if(CRT_hit.front()->peshit > 70) evt_CRTveto_70PE = true;
                if(CRT_hit.front()->peshit > 100) evt_CRTveto = true;
              } // if CRT veto
            } // loop over flash(es)
          } // if flash exists

        } // Using CRT

        // - crt time and t0 time match
        if(trigger_time != -999 && trk_crt_time != -999){
          if(abs(trk_crt_time - trigger_time) < 1){
            if_t0_trk_crt_time_match =  true;
          }
          else{
            if_t0_trk_crt_time_match = false;
          }
        }

        //-- Fill RECO track info (in the naive version this is selected)

        // Add spatial correction to the track start and end
        Trk_vtx = daughter_Tracks.front()->Vertex<TVector3>();
        auto Trk_vtx_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_vtx.X(), Trk_vtx.Y(), Trk_vtx.Z()));
        Trk_vtx_SCEcorr.SetX(Trk_vtx.X() - Trk_vtx_offset.X() + xtimeoffset);
        Trk_vtx_SCEcorr.SetY(Trk_vtx.Y() + Trk_vtx_offset.Y());
        Trk_vtx_SCEcorr.SetZ(Trk_vtx.Z() + Trk_vtx_offset.Z());

        trk_vtx_x = Trk_vtx_SCEcorr.X();
        trk_vtx_y = Trk_vtx_SCEcorr.Y();
        trk_vtx_z = Trk_vtx_SCEcorr.Z();

        trk_vtx_x_noSCE = Trk_vtx.X();
        trk_vtx_y_noSCE = Trk_vtx.Y();
        trk_vtx_z_noSCE = Trk_vtx.Z();

        Trk_start = daughter_Tracks.front()->Start<TVector3>();
        auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
        Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X() + xtimeoffset);
        Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
        Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

        trk_start_x = Trk_start_SCEcorr.X();
        trk_start_y = Trk_start_SCEcorr.Y();
        trk_start_z = Trk_start_SCEcorr.Z();

        trk_start_x_noSCE = Trk_start.X();
        trk_start_y_noSCE = Trk_start.Y();
        trk_start_z_noSCE = Trk_start.Z();

        Trk_end = daughter_Tracks.front()->End<TVector3>();
        auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
        Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X() + xtimeoffset);
        Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
        Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

        trk_end_x = Trk_end_SCEcorr.X();
        trk_end_y = Trk_end_SCEcorr.Y();
        trk_end_z = Trk_end_SCEcorr.Z();

        trk_end_x_noSCE = Trk_end.X();
        trk_end_y_noSCE = Trk_end.Y();
        trk_end_z_noSCE = Trk_end.Z();

        //-- if either of the track end is out of the time, label them
        if(trk_start_x < 0 || trk_start_x > 2. * geo->DetHalfWidth() || trk_end_x < 0 || trk_end_x > 2. * geo->DetHalfWidth()){
          trk_OutOfTime = true;
        }
        else{
          trk_OutOfTime = false;
        }

        //-- track containment
        trk_contained = _fiducial_volume.TrackContain(Trk_start_SCEcorr, Trk_end_SCEcorr);

        vtx_InFV = _fiducial_volume.VertexInFV(Trk_start_SCEcorr);

        //-- Preliminary broken track searching (For the moment, only support 2 track merging)
        BrokenTrack brokentrack;
        brokentrack.MatchTracks(daughter_Tracks.front(), AllTrackCollection);
        if(brokentrack.NewTrk()){
          if_broken = true;
          trk_broken_len = brokentrack.TrkLen();
          trk_broken_nr_merged = brokentrack.NumberMergedTracks();

          TVector3 trk_end1 = brokentrack.TrkEnd1();
          TVector3 trk_end2 = brokentrack.TrkEnd2();

          trk_merged_ifcontained = _fiducial_volume.TrackContain(trk_end1, trk_end2);

          if(!_fiducial_volume.PointContain(trk_end1) && !_fiducial_volume.PointContain(trk_end2)){
            vtx_merged_InFV = false;
            if_newTrkThroughGoing = true;
          }
          else{
            vtx_merged_InFV = vtx_InFV;
          }
        }
        else{
          trk_merged_ifcontained = trk_contained;
          vtx_merged_InFV = vtx_InFV;
        }

        //-- Angle wrt to the readout; theta is to the wire; phi is to x
        auto Trk_pos = Trk_end_SCEcorr - Trk_start_SCEcorr;
        double theta_pl2 = std::atan2(Trk_pos.Z(), Trk_pos.Y()); // atan2(y,x)
        double theta_pl1 = theta_pl2 + M_PI/3; // If plan1 is -60 degree to Y, looking from outside to the TPC
        double theta_pl0 = theta_pl2 - M_PI/3; // If plan0 is +60 degree to Y, looking from outside to the TPC
        double phi_readout = std::atan2(Trk_pos.X(), Trk_pos.Y());

        sin2_theta_pl2 = sin(theta_pl2) * sin(theta_pl2);
        sin2_theta_pl1 = sin(theta_pl1) * sin(theta_pl1);
        sin2_theta_pl0 = sin(theta_pl0) * sin(theta_pl0);
        sin2_phi_readout = sin(phi_readout) * sin(phi_readout);

        //-- track angle
        trk_phi = daughter_Tracks.front()->Phi();
        trk_theta = daughter_Tracks.front()->Theta();
        trk_costheta = cos(daughter_Tracks.front()->Theta());

        //-- track length
        trk_length = daughter_Tracks.front()->Length();

        //-- Momentum
        //- range
        mom_Range_mu = _trk_mom_calculator.GetTrackMomentum(trk_length, 13);
        mom_Range_p = _trk_mom_calculator.GetTrackMomentum(trk_length, 2212);
        mom_Range_pi = _trk_mom_calculator.GetTrackMomentum(trk_length, 211);
        
        //- MCS
        bestMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bestMomentum();
        bestMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bestLogLikelihood();
        fwdMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->fwdMomentum();
        fwdMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->fwdLogLikelihood();
        bwdMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bwdMomentum();
        bwdMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bwdLogLikelihood();

        bestMCS_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->bestMomentum();
        bestMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->bestLogLikelihood();
        fwdMCS_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->fwdMomentum();
        fwdMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->fwdLogLikelihood();
        bwdMCS_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->bwdMomentum();
        bwdMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->bwdLogLikelihood();

        //-- PT and PL
         PT_range = abs(mom_Range_mu * sin(trk_theta));
        PL_range = abs(mom_Range_mu * cos(trk_theta));
        PT_MCS = abs(bestMCS * sin(trk_theta));
        PL_MCS = abs(bestMCS * cos(trk_theta));

        //-- dE/dx
        // pandoracaliSCE has E-field and spatial correction
        auto assoCal = trackToCalAsso.at(daughter_Tracks.front().key());
        if(assoCal.size()!=3){
          throw cet::exception("[Numu0pi0p]") << "Where are the three planes for the calorimetry!" << std::endl;
        }
        // induction 0 = 0, induction 1 = 1, collection = 2 (Check if the plane ID is correct)
        // assoCal[id_pl]->PlaneID().Plane == 2 (collection)
        // The vector is ordered by residual range from small to big (track end to track start)
        dEdx_pl0 = assoCal[0]->dEdx();
        dQdx_pl0 = assoCal[0]->dQdx();
        resRange_pl0 = assoCal[0]->ResidualRange();

        dEdx_pl1 = assoCal[1]->dEdx();
        dQdx_pl1 = assoCal[1]->dQdx();
        resRange_pl1 = assoCal[1]->ResidualRange();

        dEdx_pl2 = assoCal[2]->dEdx();
        dQdx_pl2 = assoCal[2]->dQdx();
        resRange_pl2 = assoCal[2]->ResidualRange();

        auto hits_dEdx_size_pl0 = dEdx_pl0.size();
        auto hits_dEdx_size_pl1 = dEdx_pl1.size();
        auto hits_dEdx_size_pl2 = dEdx_pl2.size();

        auto half_size_pl0 = hits_dEdx_size_pl0 / 2;
        auto half_size_pl1 = hits_dEdx_size_pl1 / 2;
        auto half_size_pl2 = hits_dEdx_size_pl2 / 2;

        //- ratio
        // dEdx_half
        if(half_size_pl0 != 0){
          dEdx_pl0_start_half = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
          dEdx_pl0_end_half = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
        }
        else{
          dEdx_pl0_start_half = 0;
          dEdx_pl0_end_half = 0;
        }
        if(half_size_pl1 != 0){
          dEdx_pl1_start_half = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
          dEdx_pl1_end_half = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
        }
        else{
          dEdx_pl1_start_half = 0;
          dEdx_pl1_end_half = 0;
        }
        if(half_size_pl2 != 0){
          dEdx_pl2_start_half = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
          dEdx_pl2_end_half = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;
        }
        else{
          dEdx_pl2_start_half = 0;
          dEdx_pl2_end_half = 0;
        }
        // dEdx_1020
        if (dEdx_pl0.size()<=40) {
          dEdx_pl0_start1020 = dEdx_pl0_start_half;
          dEdx_pl0_end1020 = dEdx_pl0_end_half;
        }
        else{
          dEdx_pl0_start1020 = std::accumulate(dEdx_pl0.end() - 20, dEdx_pl0.end() - 10, 0.) / 10.;
          dEdx_pl0_end1020 = std::accumulate(dEdx_pl0.begin() + 10, dEdx_pl0.begin() + 20, 0.) / 10.;
        }
        if (dEdx_pl1.size()<=40) {
          dEdx_pl1_start1020 = dEdx_pl1_start_half;
          dEdx_pl1_end1020 = dEdx_pl1_end_half;
        }
        else{
          dEdx_pl1_start1020 = std::accumulate(dEdx_pl1.end() - 20, dEdx_pl1.end() - 10, 0.) / 10.;
          dEdx_pl1_end1020 = std::accumulate(dEdx_pl1.begin() + 10, dEdx_pl1.begin() + 20, 0.) / 10.;
        }
        if (dEdx_pl2.size()<=40) {
          dEdx_pl2_start1020 = dEdx_pl2_start_half;
          dEdx_pl2_end1020 = dEdx_pl2_end_half;
        }
        else{
          dEdx_pl2_start1020 = std::accumulate(dEdx_pl2.end() - 20, dEdx_pl2.end() - 10, 0.) / 10.;
          dEdx_pl2_end1020 = std::accumulate(dEdx_pl2.begin() + 10, dEdx_pl2.begin() + 20, 0.) / 10.;
        }

        if((dEdx_pl0_end_half + dEdx_pl0_start_half) != 0){
          dEdx_pl0_1020_ratio = dEdx_pl0_end1020 / (dEdx_pl0_end1020 + dEdx_pl0_start1020);
          dEdx_pl0_half_ratio = dEdx_pl0_end_half / (dEdx_pl0_end_half + dEdx_pl0_start_half);
        }
        else{
          dEdx_pl0_1020_ratio = 0;
          dEdx_pl0_half_ratio = 0;
        }

        if((dEdx_pl1_end_half + dEdx_pl1_start_half) != 0){
          dEdx_pl1_1020_ratio = dEdx_pl1_end1020 / (dEdx_pl1_end1020 + dEdx_pl1_start1020);
          dEdx_pl1_half_ratio = dEdx_pl1_end_half / (dEdx_pl1_end_half + dEdx_pl1_start_half);
        }
        else{
          dEdx_pl1_1020_ratio = 0;
          dEdx_pl1_half_ratio = 0;
        }

        if((dEdx_pl2_end_half + dEdx_pl2_start_half) != 0){
          dEdx_pl2_1020_ratio = dEdx_pl2_end1020 / (dEdx_pl2_end1020 + dEdx_pl2_start1020);
          dEdx_pl2_half_ratio = dEdx_pl2_end_half / (dEdx_pl2_end_half + dEdx_pl2_start_half);
        }
        else{
          dEdx_pl2_1020_ratio = 0;
          dEdx_pl2_half_ratio = 0;
        }

//TODO check if the clean dE/dx perform better
        //-- Clean dE/dx, remove dE/dx is below 0.5 and above 70 MeV (artificial)
        std::vector<float> dEdx_pl0_clean;
        std::vector<float> resRange_pl0_clean;
        std::vector<float> dEdx_pl1_clean;
        std::vector<float> resRange_pl1_clean;
        std::vector<float> dEdx_pl2_clean;
        std::vector<float> resRange_pl2_clean;

        if(hits_dEdx_size_pl0 != 0){
          for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl0; i_hit++){
            if(dEdx_pl0[i_hit] > 0.5 && dEdx_pl0[i_hit] < 70){
              dEdx_pl0_clean.push_back(dEdx_pl0[i_hit]);
              resRange_pl0_clean.push_back(resRange_pl0[i_hit]);
            }
          }
        }

        if(hits_dEdx_size_pl1 != 0){
          for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl1; i_hit++){
            if(dEdx_pl1[i_hit] > 0.5 && dEdx_pl1[i_hit] < 70){
              dEdx_pl1_clean.push_back(dEdx_pl1[i_hit]);
              resRange_pl1_clean.push_back(resRange_pl1[i_hit]);
            }
          }
        }

        if(hits_dEdx_size_pl2 != 0){
          for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl2; i_hit++){
            if(dEdx_pl2[i_hit] > 0.5 && dEdx_pl2[i_hit] < 70){
              dEdx_pl2_clean.push_back(dEdx_pl2[i_hit]);
              resRange_pl2_clean.push_back(resRange_pl2[i_hit]);
            }
          }
        }

        auto hits_dEdx_size_pl0_clean = dEdx_pl0_clean.size();
        auto hits_dEdx_size_pl1_clean = dEdx_pl1_clean.size();
        auto hits_dEdx_size_pl2_clean = dEdx_pl2_clean.size();

        auto half_size_pl0_clean = hits_dEdx_size_pl0_clean / 2;
        auto half_size_pl1_clean = hits_dEdx_size_pl1_clean / 2;
        auto half_size_pl2_clean = hits_dEdx_size_pl2_clean / 2;

        //- ratio
        // dEdx_half
        if(half_size_pl0_clean != 0){
          dEdx_pl0_start_half_clean = std::accumulate(dEdx_pl0_clean.end() - half_size_pl0_clean, dEdx_pl0_clean.end(), 0.) / half_size_pl0_clean;
          dEdx_pl0_end_half_clean = std::accumulate(dEdx_pl0_clean.begin(), dEdx_pl0_clean.begin() + half_size_pl0_clean, 0. ) / half_size_pl0_clean;
        }
        else{
          dEdx_pl0_start_half_clean = 0;
          dEdx_pl0_end_half_clean = 0;
        }
        if(half_size_pl1_clean != 0){
          dEdx_pl1_start_half_clean = std::accumulate(dEdx_pl1_clean.end() - half_size_pl1_clean, dEdx_pl1_clean.end(), 0.) / half_size_pl1_clean;
          dEdx_pl1_end_half_clean = std::accumulate(dEdx_pl1_clean.begin(), dEdx_pl1_clean.begin() + half_size_pl1_clean, 0. ) / half_size_pl1_clean;
        }
        else{
          dEdx_pl1_start_half_clean = 0;
          dEdx_pl1_end_half_clean = 0;
        }
        if(half_size_pl2_clean != 0){
          dEdx_pl2_start_half_clean = std::accumulate(dEdx_pl2_clean.end() - half_size_pl2_clean, dEdx_pl2_clean.end(), 0.) / half_size_pl2_clean;
          dEdx_pl2_end_half_clean = std::accumulate(dEdx_pl2_clean.begin(), dEdx_pl2_clean.begin() + half_size_pl2_clean, 0. ) / half_size_pl2_clean;
        }
        else{
          dEdx_pl2_start_half_clean = 0;
          dEdx_pl2_end_half_clean = 0;
        }

        // dEdx_1020
        if (dEdx_pl0_clean.size()<=40) {
          dEdx_pl0_start1020_clean = dEdx_pl0_start_half_clean;
          dEdx_pl0_end1020_clean = dEdx_pl0_end_half_clean;
        }
        else{
          dEdx_pl0_start1020_clean = std::accumulate(dEdx_pl0_clean.end() - 20, dEdx_pl0_clean.end() - 10, 0.) / 10.;
          dEdx_pl0_end1020_clean = std::accumulate(dEdx_pl0_clean.begin() + 10, dEdx_pl0_clean.begin() + 20, 0.) / 10.;
        }
        if (dEdx_pl1_clean.size()<=40) {
          dEdx_pl1_start1020_clean = dEdx_pl1_start_half_clean;
          dEdx_pl1_end1020_clean = dEdx_pl1_end_half_clean;
        }
        else{
          dEdx_pl1_start1020_clean = std::accumulate(dEdx_pl1_clean.end() - 20, dEdx_pl1_clean.end() - 10, 0.) / 10.;
          dEdx_pl1_end1020_clean = std::accumulate(dEdx_pl1_clean.begin() + 10, dEdx_pl1_clean.begin() + 20, 0.) / 10.;
        }
        if (dEdx_pl2_clean.size()<=40) {
          dEdx_pl2_start1020_clean = dEdx_pl2_start_half_clean;
          dEdx_pl2_end1020_clean = dEdx_pl2_end_half_clean;
        }
        else{
          dEdx_pl2_start1020_clean = std::accumulate(dEdx_pl2_clean.end() - 20, dEdx_pl2_clean.end() - 10, 0.) / 10.;
          dEdx_pl2_end1020_clean = std::accumulate(dEdx_pl2_clean.begin() + 10, dEdx_pl2_clean.begin() + 20, 0.) / 10.;
        }

        // dEdx_5
        if (dEdx_pl0_clean.size()<=10) {
          dEdx_pl0_start5_clean = dEdx_pl0_start_half_clean;
          dEdx_pl0_end5_clean = dEdx_pl0_end_half_clean;
        }
        else{
          dEdx_pl0_start5_clean = std::accumulate(dEdx_pl0_clean.end() - 5, dEdx_pl0_clean.end(), 0.) / 5.;
          dEdx_pl0_end5_clean = std::accumulate(dEdx_pl0_clean.begin(), dEdx_pl0_clean.begin() + 5, 0.) / 5.;
        }
        if (dEdx_pl1_clean.size()<=10) {
          dEdx_pl1_start5_clean = dEdx_pl1_start_half_clean;
          dEdx_pl1_end5_clean = dEdx_pl1_end_half_clean;
        }
        else{
          dEdx_pl1_start5_clean = std::accumulate(dEdx_pl1_clean.end() - 5, dEdx_pl1_clean.end(), 0.) / 5.;
          dEdx_pl1_end5_clean = std::accumulate(dEdx_pl1_clean.begin(), dEdx_pl1_clean.begin() + 5, 0.) / 5.;
        }
        if (dEdx_pl2_clean.size()<=10) {
          dEdx_pl2_start5_clean = dEdx_pl2_start_half_clean;
          dEdx_pl2_end5_clean = dEdx_pl2_end_half_clean;
        }
        else{
          dEdx_pl2_start5_clean = std::accumulate(dEdx_pl2_clean.end() - 5, dEdx_pl2_clean.end(), 0.) / 5.;
          dEdx_pl2_end5_clean = std::accumulate(dEdx_pl2_clean.begin(), dEdx_pl2_clean.begin() + 5, 0.) / 5.;
        }

        if((dEdx_pl0_end_half_clean + dEdx_pl0_start_half_clean) != 0){
          dEdx_pl0_5_ratio_clean = dEdx_pl0_end5_clean / (dEdx_pl0_end5_clean + dEdx_pl0_start5_clean);
          dEdx_pl0_1020_ratio_clean = dEdx_pl0_end1020_clean / (dEdx_pl0_end1020_clean + dEdx_pl0_start1020_clean);
          dEdx_pl0_half_ratio_clean = dEdx_pl0_end_half_clean / (dEdx_pl0_end_half_clean + dEdx_pl0_start_half_clean);
        }
        else{
          dEdx_pl0_5_ratio_clean = 0;
          dEdx_pl0_1020_ratio_clean = 0;
          dEdx_pl0_half_ratio_clean = 0;
        }

        if((dEdx_pl1_end_half_clean + dEdx_pl1_start_half_clean)){
          dEdx_pl1_5_ratio_clean = dEdx_pl1_end5_clean / (dEdx_pl1_end5_clean + dEdx_pl1_start5_clean);
          dEdx_pl1_1020_ratio_clean = dEdx_pl1_end1020_clean / (dEdx_pl1_end1020_clean + dEdx_pl1_start1020_clean);
          dEdx_pl1_half_ratio_clean = dEdx_pl1_end_half_clean / (dEdx_pl1_end_half_clean + dEdx_pl1_start_half_clean);
        }
        else{
          dEdx_pl1_5_ratio_clean = 0;
          dEdx_pl1_1020_ratio_clean = 0;
          dEdx_pl1_half_ratio_clean = 0;
        }

        if((dEdx_pl2_end_half_clean + dEdx_pl2_start_half_clean)){
          dEdx_pl2_5_ratio_clean = dEdx_pl2_end5_clean / (dEdx_pl2_end5_clean + dEdx_pl2_start5_clean);
          dEdx_pl2_1020_ratio_clean = dEdx_pl2_end1020_clean / (dEdx_pl2_end1020_clean + dEdx_pl2_start1020_clean);
          dEdx_pl2_half_ratio_clean = dEdx_pl2_end_half_clean / (dEdx_pl2_end_half_clean + dEdx_pl2_start_half_clean);
        }
        else{
          dEdx_pl2_5_ratio_clean = 0;
          dEdx_pl2_1020_ratio_clean = 0;
          dEdx_pl2_half_ratio_clean = 0;
        }

        //-- PID
        PID pid;
        pid.Chi2(PIDTotrackAsso,daughter_Tracks.front(), Trk_start_SCEcorr, Trk_end_SCEcorr,hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2);
        PID_Chi2Mu_3pl = pid.PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
        PID_Chi2P_3pl = pid.PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID
        PID_Chi2Pi_3pl = pid.PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID
        PID_Chi2K_3pl = pid.PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

        ///////////////
        // TRUE info corresponding to the selected track
        //////////////
        if(IsMC){
          std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(daughter_Tracks.front().key());
          BackTrackerTruthMatch backtrackertruthmatch;
          backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
          //art::Ptr< simb::MCParticle > selected_MCparticle;
          //auto selected_MCparticle = backtrackertruthmatch.ReturnMCParticle();
          auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
          trk_cosmic_percent = backtrackertruthmatch.ReturnCosmicPercent();
          trk_purity = backtrackertruthmatch.ReturnPurity();
          trk_completeness = backtrackertruthmatch.ReturnCompleteness();

          if(!MCparticle){
            if_cosmic = true;
            std::cout<<"MC particle does not exist!"<<std::endl;
          }
          else{
            if(trk_cosmic_percent > 0.8){
                if_cosmic = true;
              }
            else{
              if_cosmic = false;
              if(MCparticle->Process() == "primary" && trk_purity > 0.2 && trk_completeness > 0.2){
                if_matchPrimary = true;
                if(MCparticle->PdgCode() == 13){
                  if_matchMu = true;
                }
              }
              TVector3 TrueTrackPos(MCparticle->Px(), MCparticle->Py(), MCparticle->Pz());// The initial momentum represent the angle of true track
              true_mom = MCparticle->P();
              true_start_x = MCparticle->Position().X();
              true_start_y = MCparticle->Position().Y();
              true_start_z = MCparticle->Position().Z();
              true_end_x = MCparticle->EndPosition().X();
              true_end_y = MCparticle->EndPosition().Y();
              true_end_z = MCparticle->EndPosition().Z();
              TVector3 true_start(MCparticle->Position().X(), MCparticle->Position().Y(), MCparticle->Position().Z());
              TVector3 true_end(MCparticle->EndPosition().X(), MCparticle->EndPosition().Y(), MCparticle->EndPosition().Z());
              true_trk_ifcontained = _fiducial_volume.TrackContain(true_start, true_end);
              true_vtxFV = _fiducial_volume.VertexInFV(true_start);
              true_trk_phi = TrueTrackPos.Phi();
              true_trk_theta = TrueTrackPos.Theta();
              true_trk_costheta = cos(TrueTrackPos.Theta());
              true_trk_theta_yz = std::atan2(TrueTrackPos.Y(), TrueTrackPos.Z());
              true_trk_costheta_yz = cos(true_trk_theta_yz);
              true_trk_theta_xz = std::atan2(TrueTrackPos.X(), TrueTrackPos.Z());
              true_trk_costheta_xz = cos(true_trk_theta_xz);
              true_trk_length = (true_start - true_end).Mag(); // An estimation of true track length
              true_trk_PDG = MCparticle->PdgCode();

              // pass the MCparticle out to the MC truth
              selected_MCparticle = MCparticle;
            } // not cosmic
          } // If MC particle exists
        }        

      } // if it is 1 single track in reco signature

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
    if (if_ReWeight){
      if(MCTruthCollection.size() != WeightCollection.size()){
        throw cet::exception("[Numu0pi0p]")<< "Reweight is not consistent with MCTruth!" << std::endl;
      }
    } 

    // we first find the matching one if there is any, otherwise the the reco will be paired with the first MC interaction
    // set up which MC to fill with reco
    if_MC_reco_match = false;
    // if there is no selected MC particle (the selected one is cosmic or so), which should be for the most of the cases
    if (if_cosmic){
      matched_MCid = 0;
    }
    else if (!if_cosmic){
      for (unsigned int i_mc = 0; i_mc < MCTruthCollection.size(); i_mc++){
        auto assoMC = MCtToMCpAsso.at(MCTruthCollection[i_mc].key());
        for (unsigned int i_mcp = 0; i_mcp < assoMC.size(); i_mcp++){
          if (assoMC[i_mcp] == selected_MCparticle){
            if_MC_reco_match = true;
            matched_MCid = i_mc;
            break;
          }
        }
        if (if_MC_reco_match) break;
      }
    }

    // if nothing matches even MC particle exists, set it to the first MC interaction  
    if (!if_MC_reco_match){
      matched_MCid = 0;
    }


//TODO add if_MC_reco_match to ttree

    // Fill a pair of reco and truth together
    if (if_ReWeight){
      //-- Xsec reweight
      std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection.at(matched_MCid)->fWeight;
      for (const auto& this_par: genie_pars){
        const std::vector<double> &weights = evtwgt_map.at(this_par);

        if (this_par == "MaCCQE") { MaCCQE = weights; }
        if (this_par == "CoulombCCQE") { CoulombCCQE = weights; }
        if (this_par == "MaNCEL") { MaNCEL = weights; }
        if (this_par == "EtaNCEL") { EtaNCEL = weights; }

        if (this_par == "NormCCMEC") { NormCCMEC = weights; }
        if (this_par == "NormNCMEC") { NormNCMEC = weights; }
        if (this_par == "FracPN_CCMEC") { FracPN_CCMEC = weights; }
        if (this_par == "FracDelta_CCMEC") { FracDelta_CCMEC = weights; }

        if (this_par == "MaCCRES") { MaCCRES = weights; }
        if (this_par == "MvCCRES") { MvCCRES = weights; }
        if (this_par == "MaNCRES") { MaNCRES = weights; }
        if (this_par == "MvNCRES") { MvNCRES = weights; }

        if (this_par == "NonRESBGvpCC1pi") { NonRESBGvpCC1pi = weights; }
        if (this_par == "NonRESBGvpCC2pi") { NonRESBGvpCC2pi = weights; }
        if (this_par == "NonRESBGvpNC1pi") { NonRESBGvpNC1pi = weights; }
        if (this_par == "NonRESBGvpNC2pi") { NonRESBGvpNC2pi = weights; }
        if (this_par == "NonRESBGvnCC1pi") { NonRESBGvnCC1pi = weights; }
        if (this_par == "NonRESBGvnCC2pi") { NonRESBGvnCC2pi = weights; }
        if (this_par == "NonRESBGvnNC1pi") { NonRESBGvnNC1pi = weights; }
        if (this_par == "NonRESBGvnNC2pi") { NonRESBGvnNC2pi = weights; }
        if (this_par == "NonRESBGvbarpCC1pi") { NonRESBGvbarpCC1pi = weights; }
        if (this_par == "NonRESBGvbarpCC2pi") { NonRESBGvbarpCC2pi = weights; }
        if (this_par == "NonRESBGvbarpNC1pi") { NonRESBGvbarpNC1pi = weights; }
        if (this_par == "NonRESBGvbarpNC2pi") { NonRESBGvbarpNC2pi = weights; }
        if (this_par == "NonRESBGvbarnCC1pi") { NonRESBGvbarnCC1pi = weights; }
        if (this_par == "NonRESBGvbarnCC2pi") { NonRESBGvbarnCC2pi = weights; }
        if (this_par == "NonRESBGvbarnNC1pi") { NonRESBGvbarnNC1pi = weights; }
        if (this_par == "NonRESBGvbarnNC2pi") { NonRESBGvbarnNC2pi = weights; }
        if (this_par == "AhtBY") { AhtBY = weights; }
        if (this_par == "BhtBY") { BhtBY = weights; }
        if (this_par == "CV1uBY") { CV1uBY = weights; }
        if (this_par == "CV2uBY") { CV2uBY = weights; }

        if (this_par == "AGKYxF1pi") { AGKYxF1pi = weights; }
        if (this_par == "AGKYpT1pi") { AGKYpT1pi = weights; }

        if (this_par == "MFP_pi") { MFP_pi = weights; }
        if (this_par == "MFP_N") { MFP_N = weights; }
        if (this_par == "FrCEx_pi") { FrCEx_pi = weights; }
        if (this_par == "FrInel_pi") { FrInel_pi = weights; }
        if (this_par == "FrAbs_pi") { FrAbs_pi = weights; }
        if (this_par == "FrCEx_N") { FrCEx_N = weights; }
        if (this_par == "FrInel_N") { FrInel_N = weights; }
        if (this_par == "FrAbs_N") { FrAbs_N = weights; }

        if (this_par == "RDecBR1gamma") { RDecBR1gamma = weights; }
        if (this_par == "RDecBR1eta") { RDecBR1eta = weights; }
      } // GENIE
    } // if rewight

    // The neutrino info
    if (MCTruthCollection[matched_MCid]->Origin() == Neutrino_Origin){ MC_beamNeutrino = true;}
    MC_int_mode = MCTruthCollection[matched_MCid]->GetNeutrino().Mode();
    MC_nupdg = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().PdgCode();
    MC_ccnc = MCTruthCollection[matched_MCid]->GetNeutrino().CCNC();
    MC_Q2 = MCTruthCollection[matched_MCid]->GetNeutrino().QSqr();

    MC_nu_E = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().E();
    MC_nuVtxX = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().Vx();
    MC_nuVtxY = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().Vy();
    MC_nuVtxZ = MCTruthCollection[matched_MCid]->GetNeutrino().Nu().Vz();

    true_nuVtx.SetXYZ(MC_nuVtxX, MC_nuVtxY, MC_nuVtxZ);
    MC_FV = _fiducial_volume.VertexInFV(true_nuVtx);
    MC_if_in_active = _fiducial_volume.VertexInActive(true_nuVtx);

    // The final state particles
    // Loop all the MCParticles to determine the true topology (all the MCParticles are from the neutrino events in overlay)
    // Not necessary all the Genie particles go through the geant4 stage?
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

      MC_0pi0p_PT = sqrt(MC_muon_true_Px.front() * MC_muon_true_Px.front() + MC_muon_true_Py.front() * MC_muon_true_Py.front());
      MC_0pi0p_PL = abs(MC_muon_true_Pz.front());
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

    }

    if(TopologyType == 3){
      MC_0pi2p_proton1_mom = MC_proton_true_Mom[0];
      MC_0pi2p_proton2_mom = MC_proton_true_Mom[1];
    }

    if(MC_nMuon == 1){
      MC_transfer_E = MC_nu_E - MC_muon_true_E.front();
    }
    my_event_->Fill();

//TODO test Fill()

    // [wipe the reco info] if if_1track is false, it is not selected. Therefore, the reco info won't be considered
    if_1track = false;

    // If there is more MC to fill
    if (MCTruthCollection.size() > 1){
      for (unsigned int i_mc = 0; i_mc < MCTruthCollection.size(); i_mc++){
        if (i_mc != matched_MCid){
          if (if_ReWeight){
            //-- Xsec reweight
            std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection.at(i_mc)->fWeight;
            for (const auto& this_par: genie_pars){
              const std::vector<double> &weights = evtwgt_map.at(this_par);

              if (this_par == "MaCCQE") { MaCCQE = weights; }
              if (this_par == "CoulombCCQE") { CoulombCCQE = weights; }
              if (this_par == "MaNCEL") { MaNCEL = weights; }
              if (this_par == "EtaNCEL") { EtaNCEL = weights; }

              if (this_par == "NormCCMEC") { NormCCMEC = weights; }
              if (this_par == "NormNCMEC") { NormNCMEC = weights; }
              if (this_par == "FracPN_CCMEC") { FracPN_CCMEC = weights; }
              if (this_par == "FracDelta_CCMEC") { FracDelta_CCMEC = weights; }

              if (this_par == "MaCCRES") { MaCCRES = weights; }
              if (this_par == "MvCCRES") { MvCCRES = weights; }
              if (this_par == "MaNCRES") { MaNCRES = weights; }
              if (this_par == "MvNCRES") { MvNCRES = weights; }

              if (this_par == "NonRESBGvpCC1pi") { NonRESBGvpCC1pi = weights; }
              if (this_par == "NonRESBGvpCC2pi") { NonRESBGvpCC2pi = weights; }
              if (this_par == "NonRESBGvpNC1pi") { NonRESBGvpNC1pi = weights; }
              if (this_par == "NonRESBGvpNC2pi") { NonRESBGvpNC2pi = weights; }
              if (this_par == "NonRESBGvnCC1pi") { NonRESBGvnCC1pi = weights; }
              if (this_par == "NonRESBGvnCC2pi") { NonRESBGvnCC2pi = weights; }
              if (this_par == "NonRESBGvnNC1pi") { NonRESBGvnNC1pi = weights; }
              if (this_par == "NonRESBGvnNC2pi") { NonRESBGvnNC2pi = weights; }
              if (this_par == "NonRESBGvbarpCC1pi") { NonRESBGvbarpCC1pi = weights; }
              if (this_par == "NonRESBGvbarpCC2pi") { NonRESBGvbarpCC2pi = weights; }
              if (this_par == "NonRESBGvbarpNC1pi") { NonRESBGvbarpNC1pi = weights; }
              if (this_par == "NonRESBGvbarpNC2pi") { NonRESBGvbarpNC2pi = weights; }
              if (this_par == "NonRESBGvbarnCC1pi") { NonRESBGvbarnCC1pi = weights; }
              if (this_par == "NonRESBGvbarnCC2pi") { NonRESBGvbarnCC2pi = weights; }
              if (this_par == "NonRESBGvbarnNC1pi") { NonRESBGvbarnNC1pi = weights; }
              if (this_par == "NonRESBGvbarnNC2pi") { NonRESBGvbarnNC2pi = weights; }
              if (this_par == "AhtBY") { AhtBY = weights; }
              if (this_par == "BhtBY") { BhtBY = weights; }
              if (this_par == "CV1uBY") { CV1uBY = weights; }
              if (this_par == "CV2uBY") { CV2uBY = weights; }

              if (this_par == "AGKYxF1pi") { AGKYxF1pi = weights; }
              if (this_par == "AGKYpT1pi") { AGKYpT1pi = weights; }

              if (this_par == "MFP_pi") { MFP_pi = weights; }
              if (this_par == "MFP_N") { MFP_N = weights; }
              if (this_par == "FrCEx_pi") { FrCEx_pi = weights; }
              if (this_par == "FrInel_pi") { FrInel_pi = weights; }
              if (this_par == "FrAbs_pi") { FrAbs_pi = weights; }
              if (this_par == "FrCEx_N") { FrCEx_N = weights; }
              if (this_par == "FrInel_N") { FrInel_N = weights; }
              if (this_par == "FrAbs_N") { FrAbs_N = weights; }

              if (this_par == "RDecBR1gamma") { RDecBR1gamma = weights; }
              if (this_par == "RDecBR1eta") { RDecBR1eta = weights; }
            } // GENIE
          } // if rewight

          // The neutrino info
          if (MCTruthCollection[i_mc]->Origin() == Neutrino_Origin){ MC_beamNeutrino = true;}
          MC_int_mode = MCTruthCollection[i_mc]->GetNeutrino().Mode();
          MC_nupdg = MCTruthCollection[i_mc]->GetNeutrino().Nu().PdgCode();
          MC_ccnc = MCTruthCollection[i_mc]->GetNeutrino().CCNC();
          MC_Q2 = MCTruthCollection[i_mc]->GetNeutrino().QSqr();

          MC_nu_E = MCTruthCollection[i_mc]->GetNeutrino().Nu().E();
          MC_nuVtxX = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vx();
          MC_nuVtxY = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vy();
          MC_nuVtxZ = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vz();

          true_nuVtx.SetXYZ(MC_nuVtxX, MC_nuVtxY, MC_nuVtxZ);
          MC_FV = _fiducial_volume.VertexInFV(true_nuVtx);
          MC_if_in_active = _fiducial_volume.VertexInActive(true_nuVtx);

          // The final state particles 
          // Loop all the MCParticles to determine the true topology (all the MCParticles are from the neutrino events in overlay)
          // Not necessary all the Genie particles go through the geant4 stage?
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
  
            MC_0pi0p_PT = sqrt(MC_muon_true_Px.front() * MC_muon_true_Px.front() + MC_muon_true_Py.front() * MC_muon_true_Py.front());
            MC_0pi0p_PL = abs(MC_muon_true_Pz.front());
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
  
          }
  
          if(TopologyType == 3){
            MC_0pi2p_proton1_mom = MC_proton_true_Mom[0];
            MC_0pi2p_proton2_mom = MC_proton_true_Mom[1];
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

    MC_0pi0p_muon_mom = -999;
    MC_0pi0p_muon_theta = -999;
    MC_0pi0p_muon_costheta = -999;
    MC_0pi0p_muon_phi = -999;

    MC_0pi0p_PT = -999;
    MC_0pi0p_PL = -999;

    MC_0pi1p_muon_mom = -999;
    MC_0pi1p_muon_theta = -999;
    MC_0pi1p_muon_costheta = -999;
    MC_0pi1p_muon_phi = -999;

    MC_0pi1p_proton_mom = -999;
    MC_0pi1p_proton_theta = -999;
    MC_0pi1p_proton_costheta = -999;
    MC_0pi1p_proton_phi = -999;

    MC_0pi1p_cos_ang_muon_proton = -999;

    MC_0pi2p_proton1_mom = -999;
    MC_0pi2p_proton2_mom = -999;

    /////////
    true_mom = -999;
    true_start_x = -999;
    true_start_y = -999;
    true_start_z = -999;
    true_end_x = -999;
    true_end_y = -999;
    true_end_z = -999;
    true_trk_phi = -999;
    true_trk_theta = -999;
    true_trk_costheta = -999;
    true_trk_theta_yz = -999;
    true_trk_costheta_yz = -999;
    true_trk_theta_xz = -999;
    true_trk_costheta_xz = -999;
    true_trk_length = -999;
    true_trk_PDG = -999;
    true_trk_ifcontained = -999;
    true_vtxFV = -999;
    reco_MC_dist_vtx = -999;
    reco_MC_dist_vtx_noSCE = -999;

    trk_cosmic_percent = -999;
    trk_purity = -999;
    trk_completeness = -999;
    if_cosmic = false;
    if_matchPrimary = false;
    if_matchMu = false;

    if_MC_reco_match = false;    
  }

  daughter_Tracks.clear();
  daughter_Showers.clear();
  Track_PDG.clear();
  Shower_PDG.clear();
  Ghost_PDG.clear();

  n_pfp_nuDaughters = 0; // number of pfp which are the daughters of the neutrino
  n_dau_tracks = 0; // number of tracks asssociated to pfp neutrino daughters
  n_dau_showers = 0; // number of showers asssociated to pfp neutrino daughters

  nr_granddau_shw = 0;
  nr_granddau_trk = 0;
  nr_granddau = 0;
  MC_granddau_pdg.clear();
  granddau_trk_len.clear();
  granddau_shw_len.clear();

  if_1track = false;

  flash_matching_chi2 = -999; //Chi2 of flash matching in each neutrino slice
  trigger_time = -999;
  flash_YCenter = -999;
  flash_YWidth = -999;
  flash_ZCenter = -999;
  flash_ZWidth = -999;
  flash_TotalPE = -999;
  flash_time = -999;

  CRTT0corr = 0;
  crthit_PE.clear();
  crthit_plane.clear();
  crthit_time.clear();
  Nr_crthit_inBeam = 0;
  only_crthit_time_inBeam = -999;
  if_t0_crt_time_inBeam_match = true;

  Nr_trk_asso_crthit = 0;
  trk_crt_time = -999;
  evt_CRTveto_70 = false;
  evt_CRTveto_100 = false;
  evt_CRT_veto_homebrew = false;
  if_trk_CRT_out_Beam = false;
  if_t0_trk_crt_time_match = true;

  trk_vtx_x = -999;
  trk_vtx_y = -999;
  trk_vtx_z = -999;
  trk_start_x = -999;
  trk_start_y = -999;
  trk_start_z = -999;
  trk_end_x = -999;
  trk_end_y = -999;
  trk_end_z = -999;

  trk_vtx_x_noSCE = -999;
  trk_vtx_y_noSCE = -999;
  trk_vtx_z_noSCE = -999;
  trk_start_x_noSCE = -999;
  trk_start_y_noSCE = -999;
  trk_start_z_noSCE = -999;
  trk_end_x_noSCE = -999;
  trk_end_y_noSCE = -999;
  trk_end_z_noSCE = -999;

  trk_OutOfTime = false;

  trk_contained = false;
  vtx_InFV = false;

  if_broken = false;
  trk_broken_len = -999;
  trk_broken_nr_merged = -999;
  trk_merged_ifcontained = false;
  vtx_merged_InFV = false;
  if_newTrkThroughGoing = false;

  sin2_theta_pl0 = -999;
  sin2_theta_pl1 = -999;
  sin2_theta_pl2 = -999;
  sin2_phi_readout = -999;

  trk_phi = -999;
  trk_theta = -999;
  trk_costheta = -999;

  trk_length = -999;

  mom_Range_mu = -999;
  mom_Range_p = -999;
  mom_Range_pi = -999;

  bestMCS = -999;
  bestMCSLL = -999;
  fwdMCS = -999;
  fwdMCSLL = -999;
  bwdMCS = -999;
  bwdMCSLL = -999;

  bestMCS_NoSCE = -999;
  fwdMCS_NoSCE = -999;
  bwdMCS_NoSCE = -999;
  bestMCSLL_NoSCE = -999;
  fwdMCSLL_NoSCE = -999;
  bwdMCSLL_NoSCE = -999;

  PT_range = -999;
  PL_range = -999;
  PT_MCS = -999;
  PL_MCS = -999;

  PID_Chi2Mu_3pl = -999;
  PID_Chi2P_3pl = -999;
  PID_Chi2Pi_3pl = -999;
  PID_Chi2K_3pl = -999;

  dEdx_pl0.clear();
  dQdx_pl0.clear();
  resRange_pl0.clear();

  dEdx_pl1.clear();
  dQdx_pl1.clear();
  resRange_pl1.clear();

  dEdx_pl2.clear();
  dQdx_pl2.clear();
  resRange_pl2.clear();

  dEdx_pl0_start_half = -999;
  dEdx_pl1_start_half = -999;
  dEdx_pl2_start_half = -999;
  dEdx_pl0_end_half = -999;
  dEdx_pl1_end_half = -999;
  dEdx_pl2_end_half = -999;

  dEdx_pl0_start1020 = -999;
  dEdx_pl1_start1020 = -999;
  dEdx_pl2_start1020 = -999;
  dEdx_pl0_end1020 = -999;
  dEdx_pl1_end1020 = -999;
  dEdx_pl2_end1020 = -999;

  dEdx_pl0_1020_ratio = -999;
  dEdx_pl0_half_ratio = -999;
  dEdx_pl1_1020_ratio = -999;
  dEdx_pl1_half_ratio = -999;
  dEdx_pl2_1020_ratio = -999;
  dEdx_pl2_half_ratio = -999;

  dEdx_pl0_start_half_clean = -999;
  dEdx_pl1_start_half_clean = -999;
  dEdx_pl2_start_half_clean = -999;
  dEdx_pl0_end_half_clean = -999;
  dEdx_pl1_end_half_clean = -999;
  dEdx_pl2_end_half_clean = -999;

  dEdx_pl0_start1020_clean = -999;
  dEdx_pl1_start1020_clean = -999;
  dEdx_pl2_start1020_clean = -999;
  dEdx_pl0_end1020_clean = -999;
  dEdx_pl1_end1020_clean = -999;
  dEdx_pl2_end1020_clean = -999;

  dEdx_pl0_start5_clean = -999;
  dEdx_pl1_start5_clean = -999;
  dEdx_pl2_start5_clean = -999;
  dEdx_pl0_end5_clean = -999;
  dEdx_pl1_end5_clean = -999;
  dEdx_pl2_end5_clean = -999;

  dEdx_pl0_5_ratio_clean = -999;
  dEdx_pl0_1020_ratio_clean = -999;
  dEdx_pl0_half_ratio_clean = -999;
  dEdx_pl1_5_ratio_clean = -999;
  dEdx_pl1_1020_ratio_clean = -999;
  dEdx_pl1_half_ratio_clean = -999;
  dEdx_pl2_5_ratio_clean = -999;
  dEdx_pl2_1020_ratio_clean = -999;
  dEdx_pl2_half_ratio_clean = -999;

  if_fwd_true = true; // If fwd by the reco true vertex distance
  if_fwd_MCS = true; // If using forward MCS direction judge
  if_fwd_dEdx1020 = true; // If fwd by the reco dEdx 10 hits (should use for contained)
  if_fwd_dEdxhalf = true; // If fwd by the reco dEdx half of the hits (should use for contained)

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

    my_event_->Branch("EventWeight", &EventWeight);
    my_event_->Branch("TopologyType", &TopologyType);
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

    my_event_->Branch("MC_0pi0p_PT", &MC_0pi0p_PT);
    my_event_->Branch("MC_0pi0p_PL", &MC_0pi0p_PL);

    my_event_->Branch("MC_0pi1p_muon_mom", &MC_0pi1p_muon_mom);
    my_event_->Branch("MC_0pi1p_muon_theta", &MC_0pi1p_muon_theta);
    my_event_->Branch("MC_0pi1p_muon_costheta", &MC_0pi1p_muon_costheta);
    my_event_->Branch("MC_0pi1p_muon_phi", &MC_0pi1p_muon_phi);
    my_event_->Branch("MC_0pi1p_proton_mom", &MC_0pi1p_proton_mom);
    my_event_->Branch("MC_0pi1p_proton_theta", &MC_0pi1p_proton_theta);
    my_event_->Branch("MC_0pi1p_proton_costheta", &MC_0pi1p_proton_costheta);
    my_event_->Branch("MC_0pi1p_proton_phi", &MC_0pi1p_proton_phi);
    my_event_->Branch("MC_0pi1p_cos_ang_muon_proton", &MC_0pi1p_cos_ang_muon_proton);

    my_event_->Branch("MC_0pi2p_proton1_mom", &MC_0pi2p_proton1_mom);
    my_event_->Branch("MC_0pi2p_proton2_mom", &MC_0pi2p_proton2_mom);

//////// selected
    my_event_->Branch("MC_muon_mom_0pi0p", &MC_muon_mom_0pi0p);
    my_event_->Branch("MC_muon_theta_0pi0p", &MC_muon_theta_0pi0p);
    my_event_->Branch("MC_muon_costheta_0pi0p", &MC_muon_costheta_0pi0p);
    my_event_->Branch("MC_muon_phi_0pi0p", &MC_muon_phi_0pi0p);

    my_event_->Branch("MC_PT_0pi0p", &MC_PT_0pi0p);
    my_event_->Branch("MC_PL_0pi0p", &MC_PL_0pi0p);

    my_event_->Branch("MC_muon_mom_0pi1p", &MC_muon_mom_0pi1p);
    my_event_->Branch("MC_muon_theta_0pi1p", &MC_muon_theta_0pi1p);
    my_event_->Branch("MC_muon_costheta_0pi1p", &MC_muon_costheta_0pi1p);
    my_event_->Branch("MC_muon_phi_0pi1p", &MC_muon_phi_0pi1p);

    my_event_->Branch("MC_proton_mom_0pi1p", &MC_proton_mom_0pi1p);
    my_event_->Branch("MC_proton_theta_0pi1p", &MC_proton_theta_0pi1p);
    my_event_->Branch("MC_proton_costheta_0pi1p", &MC_proton_costheta_0pi1p);
    my_event_->Branch("MC_proton_phi_0pi1p", &MC_proton_phi_0pi1p);

    my_event_->Branch("MC_cos_ang_muon_proton_0pi1p", &MC_cos_ang_muon_proton_0pi1p);

    my_event_->Branch("MC_proton1_mom_0pi2p", &MC_proton1_mom_0pi2p);
    my_event_->Branch("MC_proton2_mom_0pi2p", &MC_proton2_mom_0pi2p);

    my_event_->Branch("MC_transfer_E_selected", &MC_transfer_E_selected);

    my_event_->Branch("TopologyType_selected", &TopologyType_selected);
///////  selected
    my_event_->Branch("Ghost_PDG", &Ghost_PDG);

    my_event_->Branch("MC_granddau_pdg", &MC_granddau_pdg);

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

    my_event_->Branch("reco_MC_dist_vtx", &reco_MC_dist_vtx);
    my_event_->Branch("reco_MC_dist_vtx_noSCE", &reco_MC_dist_vtx_noSCE);

    my_event_->Branch("if_cosmic", &if_cosmic);
    my_event_->Branch("if_matchPrimary", &if_matchPrimary);
    my_event_->Branch("if_matchMu", &if_matchMu);
    my_event_->Branch("trk_cosmic_percent", &trk_cosmic_percent);
    my_event_->Branch("trk_purity", &trk_purity);
    my_event_->Branch("trk_completeness", &trk_completeness);

  }

  my_event_->Branch("n_pfp_nuDaughters", &n_pfp_nuDaughters);
  my_event_->Branch("n_dau_tracks", &n_dau_tracks);
  my_event_->Branch("n_dau_showers", &n_dau_showers);

  my_event_->Branch("if_1track", &if_1track);

  my_event_->Branch("nr_granddau_shw", &nr_granddau_shw);
  my_event_->Branch("nr_granddau_trk", &nr_granddau_trk);
  my_event_->Branch("nr_granddau", &nr_granddau);
  my_event_->Branch("granddau_trk_len", &granddau_trk_len);
  my_event_->Branch("granddau_shw_len", &granddau_shw_len);

  my_event_->Branch("flash_matching_chi2", &flash_matching_chi2);
  my_event_->Branch("trigger_time", &trigger_time);
  my_event_->Branch("flash_YCenter", &flash_YCenter);
  my_event_->Branch("flash_YWidth", &flash_YWidth);
  my_event_->Branch("flash_ZCenter", &flash_ZCenter);
  my_event_->Branch("flash_ZWidth", &flash_ZWidth);
  my_event_->Branch("flash_TotalPE", &flash_TotalPE);
  my_event_->Branch("flash_time", &flash_time);

  my_event_->Branch("evt_CRTveto", &evt_CRTveto);
  my_event_->Branch("evt_CRTveto_100", &evt_CRTveto_100);
  my_event_->Branch("evt_CRT_veto_homebrew", &evt_CRT_veto_homebrew);
  my_event_->Branch("crthit_PE", &crthit_PE);
  my_event_->Branch("crthit_plane", &crthit_plane);
  my_event_->Branch("crthit_time", &crthit_time);
  my_event_->Branch("Nr_crthit_inBeam", &Nr_crthit_inBeam);
  my_event_->Branch("only_crthit_time_inBeam", &only_crthit_time_inBeam);
  my_event_->Branch("if_t0_crt_time_inBeam_match", &if_t0_crt_time_inBeam_match);

  my_event_->Branch("Nr_trk_asso_crthit", &Nr_trk_asso_crthit);
  my_event_->Branch("trk_crt_time", &trk_crt_time);
  my_event_->Branch("if_trk_CRT_out_Beam", &if_trk_CRT_out_Beam);
  my_event_->Branch("if_t0_trk_crt_time_match", &if_t0_trk_crt_time_match);

  my_event_->Branch("trk_vtx_x", &trk_vtx_x);
  my_event_->Branch("trk_vtx_y", &trk_vtx_y);
  my_event_->Branch("trk_vtx_z", &trk_vtx_z);
  my_event_->Branch("trk_start_x", &trk_start_x);
  my_event_->Branch("trk_start_y", &trk_start_y);
  my_event_->Branch("trk_start_z", &trk_start_z);
  my_event_->Branch("trk_end_x", &trk_end_x);
  my_event_->Branch("trk_end_y", &trk_end_y);
  my_event_->Branch("trk_end_z", &trk_end_z);

  my_event_->Branch("trk_vtx_x_noSCE", &trk_vtx_x_noSCE);
  my_event_->Branch("trk_vtx_y_noSCE", &trk_vtx_y_noSCE);
  my_event_->Branch("trk_vtx_z_noSCE", &trk_vtx_z_noSCE);
  my_event_->Branch("trk_start_x_noSCE", &trk_start_x_noSCE);
  my_event_->Branch("trk_start_y_noSCE", &trk_start_y_noSCE);
  my_event_->Branch("trk_start_z_noSCE", &trk_start_z_noSCE);
  my_event_->Branch("trk_end_x_noSCE", &trk_end_x_noSCE);
  my_event_->Branch("trk_end_y_noSCE", &trk_end_y_noSCE);
  my_event_->Branch("trk_end_z_noSCE", &trk_end_z_noSCE);

  my_event_->Branch("trk_OutOfTime", &trk_OutOfTime);

  my_event_->Branch("trk_contained", &trk_contained);
  my_event_->Branch("vtx_InFV", &vtx_InFV);

  my_event_->Branch("if_broken", &if_broken);
  my_event_->Branch("trk_broken_len", &trk_broken_len);
  my_event_->Branch("trk_broken_nr_merged", &trk_broken_nr_merged);
  my_event_->Branch("trk_merged_ifcontained", &trk_merged_ifcontained);
  my_event_->Branch("vtx_merged_InFV", &vtx_merged_InFV);
  my_event_->Branch("if_newTrkThroughGoing", &if_newTrkThroughGoing);

  my_event_->Branch("sin2_theta_pl0", &sin2_theta_pl0);
  my_event_->Branch("sin2_theta_pl1", &sin2_theta_pl1);
  my_event_->Branch("sin2_theta_pl2", &sin2_theta_pl2);
  my_event_->Branch("sin2_phi_readout", &sin2_phi_readout);

  my_event_->Branch("trk_phi", &trk_phi);
  my_event_->Branch("trk_theta", &trk_theta);
  my_event_->Branch("trk_costheta", &trk_costheta);

  my_event_->Branch("trk_length", &trk_length);

  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
  my_event_->Branch("mom_Range_p", &mom_Range_p);
  my_event_->Branch("mom_Range_pi", &mom_Range_pi);

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

  my_event_->Branch("PT_range", &PT_range);
  my_event_->Branch("PL_range", &PL_range);
  my_event_->Branch("PT_MCS", &PT_MCS);
  my_event_->Branch("PL_MCS", &PL_MCS);

  my_event_->Branch("PID_Chi2Mu_3pl", &PID_Chi2Mu_3pl);
  my_event_->Branch("PID_Chi2P_3pl", &PID_Chi2P_3pl);
  my_event_->Branch("PID_Chi2Pi_3pl", &PID_Chi2Pi_3pl);
  my_event_->Branch("PID_Chi2K_3pl", &PID_Chi2K_3pl);

  my_event_->Branch("dEdx_pl0", &dEdx_pl0);
  my_event_->Branch("dQdx_pl0", &dQdx_pl0);
  my_event_->Branch("resRange_pl0", &resRange_pl0);
  my_event_->Branch("dEdx_pl1", &dEdx_pl1);
  my_event_->Branch("dQdx_pl1", &dQdx_pl1);
  my_event_->Branch("resRange_pl1", &resRange_pl1);
  my_event_->Branch("dEdx_pl2", &dEdx_pl2);
  my_event_->Branch("dQdx_pl2", &dQdx_pl2);
  my_event_->Branch("resRange_pl2", &resRange_pl2);

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

  my_event_->Branch("dEdx_pl0_1020_ratio", &dEdx_pl0_1020_ratio);
  my_event_->Branch("dEdx_pl0_half_ratio", &dEdx_pl0_half_ratio);
  my_event_->Branch("dEdx_pl1_1020_ratio", &dEdx_pl1_1020_ratio);
  my_event_->Branch("dEdx_pl1_half_ratio", &dEdx_pl1_half_ratio);
  my_event_->Branch("dEdx_pl2_1020_ratio", &dEdx_pl2_1020_ratio);
  my_event_->Branch("dEdx_pl2_half_ratio", &dEdx_pl2_half_ratio);

  my_event_->Branch("dEdx_pl0_start_half_clean", &dEdx_pl0_start_half_clean);
  my_event_->Branch("dEdx_pl1_start_half_clean", &dEdx_pl1_start_half_clean);
  my_event_->Branch("dEdx_pl2_start_half_clean", &dEdx_pl2_start_half_clean);
  my_event_->Branch("dEdx_pl0_end_half_clean", &dEdx_pl0_end_half_clean);
  my_event_->Branch("dEdx_pl1_end_half_clean", &dEdx_pl1_end_half_clean);
  my_event_->Branch("dEdx_pl2_end_half_clean", &dEdx_pl2_end_half_clean);

  my_event_->Branch("dEdx_pl0_start1020_clean", &dEdx_pl0_start1020_clean);
  my_event_->Branch("dEdx_pl1_start1020_clean", &dEdx_pl1_start1020_clean);
  my_event_->Branch("dEdx_pl2_start1020_clean", &dEdx_pl2_start1020_clean);
  my_event_->Branch("dEdx_pl0_end1020_clean", &dEdx_pl0_end1020_clean);
  my_event_->Branch("dEdx_pl1_end1020_clean", &dEdx_pl1_end1020_clean);
  my_event_->Branch("dEdx_pl2_end1020_clean", &dEdx_pl2_end1020_clean);

  my_event_->Branch("dEdx_pl0_5_ratio_clean", &dEdx_pl0_5_ratio_clean);
  my_event_->Branch("dEdx_pl0_1020_ratio_clean", &dEdx_pl0_1020_ratio_clean);
  my_event_->Branch("dEdx_pl0_half_ratio_clean", &dEdx_pl0_half_ratio_clean);
  my_event_->Branch("dEdx_pl1_5_ratio_clean", &dEdx_pl1_5_ratio_clean);
  my_event_->Branch("dEdx_pl1_1020_ratio_clean", &dEdx_pl1_1020_ratio_clean);
  my_event_->Branch("dEdx_pl1_half_ratio_clean", &dEdx_pl1_half_ratio_clean);
  my_event_->Branch("dEdx_pl2_5_ratio_clean", &dEdx_pl2_5_ratio_clean);
  my_event_->Branch("dEdx_pl2_1020_ratio_clean", &dEdx_pl2_1020_ratio_clean);
  my_event_->Branch("dEdx_pl2_half_ratio_clean", &dEdx_pl2_half_ratio_clean);

  my_event_->Branch("if_fwd_true", &if_fwd_true);
  my_event_->Branch("if_fwd_MCS", &if_fwd_MCS);
  my_event_->Branch("if_fwd_dEdx1020", &if_fwd_dEdx1020);
  my_event_->Branch("if_fwd_dEdxhalf", &if_fwd_dEdxhalf);

}

void SingleMuon::endSubRun(art::SubRun const &sr){

  run       = sr.run();
  subrun    = sr.subRun();

  // For Overlay
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

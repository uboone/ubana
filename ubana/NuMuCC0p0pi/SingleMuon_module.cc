////////////////////////////////////////////////////////////////////////
// Class:       SingleMuon
// Plugin Type: analyzer (art v3_01_02)
// File:        SingleMuon_module.cc
//
// Generated at Thu May  23 2019 by Yifan Chen using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

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
//#include "lardataobj/AnalysisBase/PlaneIDBitsetHelperFunctions.h"

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
  explicit SingleMuon(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

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

  TTree * POTtree;
  int run, subrun;
  double POT_miss1E10;

  TTree * my_event_;
  
  void Initialize_event();

  double EventWeight; // Spine reweight using Steven's tool  
  
  bool MC_beamNeutrino; // MCTruth beam origin
  bool MC_FV; // MCTruth vertex in FV = true, out of FV = false
  bool MC_if_in_active; // MCTruth vertex in active volume = true, out of active volume = false
  int MC_ccnc; // MCTruth cc = 0 or nc = 1
  int MC_nupdg; // MCTruth nupdg; numu = 14, nue = 12
  int MC_int_mode; // https://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
  double MC_nu_E; // MCTruth nu energy
  double MC_nuVtxX; // MCTruth nu vtx X
  double MC_nuVtxY; // MCTruth nu vtx Y
  double MC_nuVtxZ; // MCTruth nu vtx Z
  int MC_nMuon; // Number of muon(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nElectron; // Number of eletron(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nNeutron; // Number of neutron(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_below260; // Number of proton(s) (p<260) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_above260; // Number of proton(s) (p >= 260) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPi0; // Number of pi0(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus_below65; // Number of pi plus(s) (p < 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus_above65; // Number of pi plus(s) (p > 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus_below65; // Number of pi minus(s) (p < 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus_above65; // Number of pi minus(s) (p > 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  std::vector<int> MC_Primary_PDG; // PDG of neutrino daughters
  std::vector<double> MC_Primary_Mom; // Momemtum of neutrino daughters
  std::vector<double> MC_proton_true_Mom_above260; // Momentum of proton above 260 MeV
  std::vector<double> MC_muon_true_Mom; // True Momentum of muon
  std::vector<double> MC_muon_true_cos_theta; // True cos theta of muon
  std::vector<double> MC_muon_true_phi; // True phi of muon
  std::vector<double> MC_proton_true_cos_theta; // True cos theta of proton
  std::vector<double> MC_proton_true_phi; // True phi of proton

  std::vector<int> Ghost_PDG; // pdg code of the pfp which has no track or shower associated; No elements ideally

  double Genie_Q2;
  double Genie_q2;
  double Genie_W;
  int Genie_nNeutron_preFSI;// before FSI 
  int Genie_nProton_preFSI;// before FSI 
  int Genie_nPi0_preFSI;// before FSI 
  int Genie_nPiPlus_preFSI;// before FSI 
  int Genie_nPiMinus_preFSI;// before FSI 

  int TopologyType;// The topology of true neutrino interaction + FSI products after Geant4
  double cos_ang_muon_proton; // cosine of the angle in between muon and proton
  double dist_muon_proton; // distance between muon start and proton start
  double len_Muon_0pi0p; // length of true muon
  double len_Muon_0pi1p; // length of true muon
  double len_Proton_0pi1p; // length of true proton

  int n_pfp_nuDaughters = 0; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks = 0; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers = 0; // number of showers asssociated to pfp neutrino daughters

  int nr_granddau_shw = 0;
  int nr_granddau_trk = 0;
  int nr_granddau = 0;
  std::vector<int> MC_granddau_pdg;
  std::vector<double> granddau_trk_len;
  std::vector<double> granddau_shw_len;

  bool if_1track = false; // If selected based on the reco info

  double flash_matching_chi2 = -999; //Chi2 of flash matching in each neutrino slice
  double flash_YCenter = -999;
  double flash_YWidth = -999;
  double flash_ZCenter = -999;
  double flash_ZWidth = -999;
  double flash_TotalPE = -999;

  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime
  int Nr_trk_asso_crthit = 0; // Number of CRT hits associated to the track
  double trk_crt_time = -999;// The CRT time of the hit which matched to the track
  bool if_trk_CRT_out_Beam = false; // Check if a track matches with out of beam CRT hit(s)
  bool evt_CRTveto = false; // If CRT veto, eliminate the events for contained (70PE threshold)
  bool evt_CRTveto_100 = false; // If CRT veto, eliminate the events for contained (100PE threshold)

  double trk_vtx_x = -999;
  double trk_vtx_y = -999;
  double trk_vtx_z = -999;
  double trk_start_x = -999;
  double trk_start_y = -999;
  double trk_start_z = -999;
  double trk_end_x = -999;
  double trk_end_y = -999;
  double trk_end_z = -999;

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

  double bestMCSLL_NoSCE = -999;
  double fwdMCSLL_NoSCE = -999;
  double bwdMCSLL_NoSCE = -999;

  double missing_PT_range = -999;
  double missing_PT_MCS = -999;

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

  double avg_dEdx_LargeHit_pl0 = -999;
  double avg_dEdx_LargeHit_pl1 = -999;
  double avg_dEdx_LargeHit_pl2 = -999;

  double dEdx_pl0_mid = -999;
  double dEdx_pl1_mid = -999;
  double dEdx_pl2_mid = -999;

  double charge_std_bin0 = -999; // the multiplication of charge in bin 0 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin1 = -999; // the multiplication of charge in bin 1 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin2 = -999; // the multiplication of charge in bin 2 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin3 = -999; // the multiplication of charge in bin 3 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin4 = -999; // the multiplication of charge in bin 4 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin5 = -999; // the multiplication of charge in bin 5 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin6 = -999; // the multiplication of charge in bin 6 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_std_bin7 = -999; // the multiplication of charge in bin 7 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections

  double charge_avg_bin0 = -999; // the multiplication of charge in bin 0 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin1 = -999; // the multiplication of charge in bin 1 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin2 = -999; // the multiplication of charge in bin 2 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin3 = -999; // the multiplication of charge in bin 3 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin4 = -999; // the multiplication of charge in bin 4 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin5 = -999; // the multiplication of charge in bin 5 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin6 = -999; // the multiplication of charge in bin 6 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  double charge_avg_bin7 = -999; // the multiplication of charge in bin 7 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections

  double vtx_hit_distance = -999;// Distance of track vertex and the closest hit spacepoints

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
 
  bool if_fwd_true = true; // If fwd by the reco true vertex distance
  bool if_fwd_MCS = true; // If using forward MCS direction judge
  bool if_fwd_dEdx1020 = true; // If fwd by the reco dEdx 10 hits (should use for contained)
  bool if_fwd_dEdxhalf = true; // If fwd by the reco dEdx half of the hits (should use for contained)

  bool                                IsMC;
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
  std::string                         m_FlashLabel;

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


SingleMuon::SingleMuon(fhicl::ParameterSet const& pset)
  : 
  EDAnalyzer{pset},
  IsMC(pset.get<bool>("IsMC")),
  UsingCRT(pset.get<bool>("UsingCRT")),
  fBeamStart(pset.get<double>("BeamStart")),
  fBeamEnd(pset.get<double>("BeamEnd")),
  fDTOffset(pset.get<double>("DTOffset")),
  fDTOffset_overlay(pset.get<double>("DTOffset_overlay")),
  m_DAQHeaderProducer(pset.get<std::string>("DAQHeaderProducer")),
  m_generatorLabel(pset.get<std::string>("GeneratorLabel")),
  m_geantLabel(pset.get<std::string>("GeantLabel")),
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
  m_FlashLabel(pset.get<std::string>("FlashLabel")),
  _min_track_len{pset.get<double>("MinTrackLength", 0.1)},
  _trk_mom_calculator{_min_track_len}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
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

  //// Get necessary handles
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  std::vector<art::Ptr<simb::GTruth> > GTruthCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
  std::vector<art::Ptr<evwgh::MCEventWeight> > WeightCollection;

  if(IsMC){
    // MC Truth
    art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
    evt.getByLabel(m_generatorLabel, Handle_MCTruth);
    art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

    // Genie Truth
    art::Handle< std::vector<simb::GTruth> > Handle_GTruth;
    evt.getByLabel(m_generatorLabel, Handle_GTruth);
    art::fill_ptr_vector(GTruthCollection, Handle_GTruth);
 
    // MC Particle
    art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;
    evt.getByLabel(m_geantLabel, Handle_MCParticle);
    art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

    // Reweight
    art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight;
    if(evt.getByLabel("eventweight4to4aFix", Handle_Weight)){
      art::fill_ptr_vector(WeightCollection, Handle_Weight);
      std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection.at(0)->fWeight;
      const std::vector<double> &weights = evtwgt_map.at("splines_general_Spline");
      EventWeight = weights.front();
    }
    else{
      EventWeight = 1;
    }
  }

  // Slice
  art::Handle<std::vector<recob::Slice> > Handle_Slice;
  evt.getByLabel(m_pandoraLabel, Handle_Slice);
  //std::vector< art::Ptr<recob::Slice> > AllSliceCollection;
  //art::fill_ptr_vector(AllSliceCollection, Handle_Slice);

  //// SpacePoint
  //art::Handle<std::vector<recob::SpacePoint> > Handle_SpacePoint;
  //evt.getByLabel(m_pandoraLabel, Handle_SpacePoint);
  //std::vector< art::Ptr<recob::SpacePoint> > AllSpacePointCollection;
  //art::fill_ptr_vector(AllSpacePointCollection, Handle_SpacePoint);

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
  //std::vector<int> Ghost_PDG; // pdg code of the pfp which has no track or shower associated; No elements ideally
  std::vector<float> vTrk_len;

  TVector3 true_nuVtx; // useful to calculate the distance to the true vtx
  TVector3 Trk_vtx;
  TVector3 Trk_start;
  TVector3 Trk_end;
  TVector3 Trk_vtx_SCEcorr;
  TVector3 Trk_start_SCEcorr; 
  TVector3 Trk_end_SCEcorr;
 
  //Constants
  const simb::Origin_t Neutrino_Origin = simb::kBeamNeutrino;

  if(IsMC){
    //------- Get part of the generator neutrino info
    //Initiate the variables
    MC_beamNeutrino = false;
    MC_nupdg = -99999;
    MC_ccnc = -99999;
    MC_FV = false;

    MC_nMuon = 0;
    MC_nElectron = 0;
    MC_nNeutron = 0;
    MC_nProton_below260 = 0;
    MC_nProton_above260 = 0;
    MC_nPi0 = 0;
    MC_nPiPlus_below65 = 0;
    MC_nPiPlus_above65 = 0;
    MC_nPiMinus_below65 = 0;
    MC_nPiMinus_above65 = 0;

    Genie_nNeutron_preFSI = 0;
    Genie_nProton_preFSI = 0;
    Genie_nPi0_preFSI = 0;
    Genie_nPiPlus_preFSI = 0;
    Genie_nPiMinus_preFSI = 0;
    
    for(unsigned int i_mc = 0; i_mc < MCTruthCollection.size(); i_mc++){
      if (MCTruthCollection[i_mc]->Origin() == Neutrino_Origin) MC_beamNeutrino = true;
      MC_int_mode = MCTruthCollection[i_mc]->GetNeutrino().Mode();
      MC_nupdg = MCTruthCollection[i_mc]->GetNeutrino().Nu().PdgCode();
      MC_ccnc = MCTruthCollection[i_mc]->GetNeutrino().CCNC();

      MC_nu_E = MCTruthCollection[i_mc]->GetNeutrino().Nu().E();
      MC_nuVtxX = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vx();
      MC_nuVtxY = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vy();
      MC_nuVtxZ = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vz();
      true_nuVtx.SetXYZ(MC_nuVtxX, MC_nuVtxY, MC_nuVtxZ);
      MC_FV = _fiducial_volume.VertexInFV(true_nuVtx);
      MC_if_in_active = _fiducial_volume.VertexInActive(true_nuVtx);
    }

    // Loop all the MCParticles to determine the true topology (all the MCParticles are from the neutrino events in overlay)
    // Not necessary all the Genie particles go through the geant4 stage?
    TVector3 MuonDir;
    TVector3 ProtonDir;
    TVector3 MuonStart;
    TVector3 ProtonStart;
    TVector3 MuonEnd;
    TVector3 ProtonEnd;
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
            MuonDir.SetXYZ(MCParticleCollection[i_mcp]->Px(), MCParticleCollection[i_mcp]->Py(),MCParticleCollection[i_mcp]->Pz());
            MuonStart.SetXYZ(MCParticleCollection[i_mcp]->Vx(), MCParticleCollection[i_mcp]->Vy(),MCParticleCollection[i_mcp]->Vz());
            MuonEnd.SetXYZ(MCParticleCollection[i_mcp]->EndX(), MCParticleCollection[i_mcp]->EndY(),MCParticleCollection[i_mcp]->EndZ());
            MC_muon_true_cos_theta.push_back(cos((MuonEnd - MuonStart).Theta()));
            MC_muon_true_phi.push_back((MuonEnd - MuonStart).Phi());
          } 
          // electron
          if(MCParticleCollection[i_mcp]->PdgCode() == 11) MC_nElectron++;
          // neutron
          if(MCParticleCollection[i_mcp]->PdgCode() == 2112) MC_nNeutron++;
          // proton
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() < 0.260) MC_nProton_below260++;
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() >= 0.260) {
            MC_nProton_above260++; 
            MC_proton_true_Mom_above260.push_back(MCParticleCollection[i_mcp]->P());
            ProtonDir.SetXYZ(MCParticleCollection[i_mcp]->Px(), MCParticleCollection[i_mcp]->Py(),MCParticleCollection[i_mcp]->Pz());
            ProtonStart.SetXYZ(MCParticleCollection[i_mcp]->Vx(), MCParticleCollection[i_mcp]->Vy(),MCParticleCollection[i_mcp]->Vz());
            ProtonEnd.SetXYZ(MCParticleCollection[i_mcp]->EndX(), MCParticleCollection[i_mcp]->EndY(),MCParticleCollection[i_mcp]->EndZ());
          }
          // pion0
          if(MCParticleCollection[i_mcp]->PdgCode() == 111) MC_nPi0++;
          // pion+
          if(MCParticleCollection[i_mcp]->PdgCode() == 211 && MCParticleCollection[i_mcp]->P() < 0.065) MC_nPiPlus_below65++;
          if(MCParticleCollection[i_mcp]->PdgCode() == 211 && MCParticleCollection[i_mcp]->P() >= 0.065) MC_nPiPlus_above65++;
          // pion-
          if(MCParticleCollection[i_mcp]->PdgCode() == -211 && MCParticleCollection[i_mcp]->P() < 0.065) MC_nPiMinus_below65++;
          if(MCParticleCollection[i_mcp]->PdgCode() == -211 && MCParticleCollection[i_mcp]->P() >= 0.065) MC_nPiMinus_above65++;
        }
      }
    }
 
    Topology topology;
    TopologyType = topology.TopologyLabel(MC_nMuon, MC_nElectron, MC_nPiPlus_above65, MC_nPiPlus_below65, MC_nPiMinus_above65, MC_nPiMinus_below65, MC_nPi0, MC_nProton_above260, MC_nProton_below260, MC_nupdg, MC_ccnc, MC_beamNeutrino, MC_FV);
    if(TopologyType == 1){
      len_Muon_0pi0p = (MuonStart - MuonEnd).Mag();
    }
    // If it is cc0pi1p then there should be only information of 1p
    if(TopologyType == 2){
      cos_ang_muon_proton = (MuonDir * ProtonDir) / (MuonDir.Mag() * ProtonDir.Mag());
      dist_muon_proton = (MuonStart - ProtonStart).Mag();
      len_Muon_0pi1p = (MuonStart - MuonEnd).Mag();
      len_Proton_0pi1p = (ProtonStart - ProtonEnd).Mag();
      MC_proton_true_cos_theta.push_back(cos((ProtonEnd - ProtonStart).Theta()));
      MC_proton_true_phi.push_back((ProtonEnd - ProtonStart).Phi());
    }    
    // Get Genie info on how many particles produced
    for(unsigned int i_gn = 0; i_gn < GTruthCollection.size(); i_gn++){
      Genie_Q2 = GTruthCollection[i_gn]->fgQ2;
      Genie_q2 = GTruthCollection[i_gn]->fgq2;
      Genie_W = GTruthCollection[i_gn]->fgW;
      Genie_nNeutron_preFSI = GTruthCollection[i_gn]->fNumNeutron;
      Genie_nProton_preFSI = GTruthCollection[i_gn]->fNumProton;
      Genie_nPi0_preFSI = GTruthCollection[i_gn]->fNumPi0;
      Genie_nPiPlus_preFSI = GTruthCollection[i_gn]->fNumPiPlus;
      Genie_nPiMinus_preFSI = GTruthCollection[i_gn]->fNumPiMinus;
    }
  }

  //-------- Get Reco neutrino (pfparticle)
  for(unsigned int i = 0; i < pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
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
        } // finish looping of pfp
      }
     
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
        }

        //-- OpFlash related information
        if(flash_v.size() == 1){
          flash_YCenter = flash_v[0]->YCenter();         
          flash_YWidth = flash_v[0]->YWidth();         
          flash_ZCenter = flash_v[0]->ZCenter();         
          flash_ZWidth = flash_v[0]->ZWidth();         
          flash_TotalPE = flash_v[0]->TotalPE();         
        }
        auto flash_matching_T0 = pfpToT0Asso.at(pfp.key());

        if(flash_matching_T0.size() == 1){
          flash_matching_chi2 = flash_matching_T0.front()->TriggerConfidence();
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
          
          //- CRT info in beam time     
          for (unsigned int i_crt = 0; i_crt < crthit_v.size(); i_crt ++){
            // figure out what plane this hit comes from
            // 3 -> top, 0 -> bottom, 1 -> anode, 2 -> cathode
            double crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);

            if(crt_time >= fBeamStart && crt_time <= fBeamEnd){
              crthit_PE.push_back(crthit_v[i_crt]->peshit);
              crthit_plane.push_back(crthit_v[i_crt]->plane);
              crthit_time.push_back(crt_time);
              Nr_crthit_inBeam++;
            } // If CRT hit in Beam window
          } // CRT hit loop
          if(Nr_crthit_inBeam != (int) crthit_time.size()) {
            throw cet::exception("[Numu0pi0p]") << "Number of crt hits does not match!" << std::endl;
          }

          //- If a track is associated to a CRT hit which is outside of beam window, exclude them (For both contained and exiting)
          auto Track_CRThit = CRTToTrackAsso.at(daughter_Tracks.front().key()); 
          Nr_trk_asso_crthit = Track_CRThit.size();
          if(Track_CRThit.size() > 0){
            for(unsigned int i_trk_hit = 0; i_trk_hit < Track_CRThit.size(); i_trk_hit++){
              trk_crt_time = ((Track_CRThit[i_trk_hit]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
              if(trk_crt_time < fBeamStart || trk_crt_time > fBeamEnd){
                if_trk_CRT_out_Beam = true;
              } // If matched CRT hit out of Beam window
            }
          }

          //- [For cross-check] For contained (Veto if there is any CRT hit in beam window)
          if(flash_v.size() > 0){
            for(unsigned int i_fl = 0; i_fl < flash_v.size(); i_fl++){
              auto CRT_hit = CRThitFlashAsso.at(flash_v[i_fl].key());
              if(CRT_hit.size() == 1){
                evt_CRTveto = true;
                if(CRT_hit.front()->peshit > 100) evt_CRTveto_100 = true;
              } // if CRT veto
            } // loop over flash(es)
          } // if flash exists

        } // Using CRT

        //-- Fill RECO track info (in the naive version this is selected)

        // Add spatial correction to the track start and end
        Trk_vtx = daughter_Tracks.front()->Vertex<TVector3>();
        auto Trk_vtx_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_vtx.X(), Trk_vtx.Y(), Trk_vtx.Z()));
        Trk_vtx_SCEcorr.SetX(Trk_vtx.X() - Trk_vtx_offset.X());
        Trk_vtx_SCEcorr.SetY(Trk_vtx.Y() + Trk_vtx_offset.Y());
        Trk_vtx_SCEcorr.SetZ(Trk_vtx.Z() + Trk_vtx_offset.Z());

        trk_vtx_x = Trk_vtx_SCEcorr.X();
        trk_vtx_y = Trk_vtx_SCEcorr.Y();
        trk_vtx_z = Trk_vtx_SCEcorr.Z();

        Trk_start = daughter_Tracks.front()->Start<TVector3>();
        auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
        Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
        Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
        Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

        trk_start_x = Trk_start_SCEcorr.X();
        trk_start_y = Trk_start_SCEcorr.Y();
        trk_start_z = Trk_start_SCEcorr.Z();

        Trk_end = daughter_Tracks.front()->End<TVector3>();
        auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
        Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X());
        Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
        Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

        trk_end_x = Trk_end_SCEcorr.X();
        trk_end_y = Trk_end_SCEcorr.Y();
        trk_end_z = Trk_end_SCEcorr.Z();

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

        bestMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->bestLogLikelihood();
        fwdMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->fwdLogLikelihood();
        bwdMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(daughter_Tracks.front().key())->bwdLogLikelihood();

        //-- Missing PT
        missing_PT_range = mom_Range_mu * sin(trk_theta);
        missing_PT_MCS = bestMCS * sin(trk_theta);

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
        dEdx_pl0_start_half = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
        dEdx_pl0_end_half = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
        dEdx_pl1_start_half = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
        dEdx_pl1_end_half = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
        dEdx_pl2_start_half = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
        dEdx_pl2_end_half = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;

        // dEdx_1020
        if (dEdx_pl0.size()<=30) {
          dEdx_pl0_start1020 = dEdx_pl0_start_half;
          dEdx_pl0_end1020 = dEdx_pl0_end_half;
        }
        else{
          dEdx_pl0_start1020 = std::accumulate(dEdx_pl0.end() - 20, dEdx_pl0.end() - 10, 0.) / 10.;
          dEdx_pl0_end1020 = std::accumulate(dEdx_pl0.begin() + 10, dEdx_pl0.begin() + 20, 0.) / 10.;
        }
        if (dEdx_pl1.size()<=30) {
          dEdx_pl1_start1020 = dEdx_pl1_start_half;
          dEdx_pl1_end1020 = dEdx_pl1_end_half;
        }
        else{
          dEdx_pl1_start1020 = std::accumulate(dEdx_pl1.end() - 20, dEdx_pl1.end() - 10, 0.) / 10.;
          dEdx_pl1_end1020 = std::accumulate(dEdx_pl1.begin() + 10, dEdx_pl1.begin() + 20, 0.) / 10.;
        }
        if (dEdx_pl2.size()<=30) {
          dEdx_pl2_start1020 = dEdx_pl2_start_half;
          dEdx_pl2_end1020 = dEdx_pl2_end_half;
        }
        else{
          dEdx_pl2_start1020 = std::accumulate(dEdx_pl2.end() - 20, dEdx_pl2.end() - 10, 0.) / 10.;
          dEdx_pl2_end1020 = std::accumulate(dEdx_pl2.begin() + 10, dEdx_pl2.begin() + 20, 0.) / 10.;
        }

        dEdx_pl0_1020_ratio = dEdx_pl0_end1020 / (dEdx_pl0_end1020 + dEdx_pl0_start1020);
        dEdx_pl0_half_ratio = dEdx_pl0_end_half / (dEdx_pl0_end_half + dEdx_pl0_start_half);

        dEdx_pl1_1020_ratio = dEdx_pl1_end1020 / (dEdx_pl1_end1020 + dEdx_pl1_start1020);
        dEdx_pl1_half_ratio = dEdx_pl1_end_half / (dEdx_pl1_end_half + dEdx_pl1_start_half);

        dEdx_pl2_1020_ratio = dEdx_pl2_end1020 / (dEdx_pl2_end1020 + dEdx_pl2_start1020);
        dEdx_pl2_half_ratio = dEdx_pl2_end_half / (dEdx_pl2_end_half + dEdx_pl2_start_half);

        // - Large hits (the three largest dE/dx)
        auto copy_dEdx_pl0 = dEdx_pl0;
        auto copy_dEdx_pl1 = dEdx_pl1;
        auto copy_dEdx_pl2 = dEdx_pl2;

        std::sort(copy_dEdx_pl0.begin(), copy_dEdx_pl0.end());
        std::sort(copy_dEdx_pl1.begin(), copy_dEdx_pl1.end());
        std::sort(copy_dEdx_pl2.begin(), copy_dEdx_pl2.end());
        // pl 0        
        if(hits_dEdx_size_pl0 < 3){
          for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl0; i_hit++){
            avg_dEdx_LargeHit_pl0 += copy_dEdx_pl0[i_hit];
          }
          avg_dEdx_LargeHit_pl0 = avg_dEdx_LargeHit_pl0 / hits_dEdx_size_pl0;
        }
        else{
          avg_dEdx_LargeHit_pl0 = (copy_dEdx_pl0[hits_dEdx_size_pl0 - 1] + copy_dEdx_pl0[hits_dEdx_size_pl0 - 2] + copy_dEdx_pl0[hits_dEdx_size_pl0 - 3]) / 3;
        }
        // pl 1
        if(hits_dEdx_size_pl1 < 3){
          for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl1; i_hit++){
            avg_dEdx_LargeHit_pl1 += copy_dEdx_pl1[i_hit];
          }
          avg_dEdx_LargeHit_pl1 = avg_dEdx_LargeHit_pl1 / hits_dEdx_size_pl1;
        }
        else{
          avg_dEdx_LargeHit_pl1 = (copy_dEdx_pl1[hits_dEdx_size_pl1 - 1] + copy_dEdx_pl1[hits_dEdx_size_pl1 - 2] + copy_dEdx_pl1[hits_dEdx_size_pl1 - 3]) / 3;
        }
        // pl 2
        if(hits_dEdx_size_pl2 < 3){
          for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl2; i_hit++){
            avg_dEdx_LargeHit_pl2 += copy_dEdx_pl2[i_hit];
          }
          avg_dEdx_LargeHit_pl2 = avg_dEdx_LargeHit_pl2 / hits_dEdx_size_pl2;
        }
        else{
          avg_dEdx_LargeHit_pl2 = (copy_dEdx_pl2[hits_dEdx_size_pl2 - 1] + copy_dEdx_pl2[hits_dEdx_size_pl2 - 2] + copy_dEdx_pl2[hits_dEdx_size_pl2 - 3]) / 3;
        }

        //- dEdx of the middle part of the track
        int nr_pl0_mid = dEdx_pl0.size()/3;
        int nr_pl1_mid = dEdx_pl1.size()/3;
        int nr_pl2_mid = dEdx_pl2.size()/3;
        if(nr_pl0_mid > 0){ //pl0
          dEdx_pl0_mid = std::accumulate(dEdx_pl0.begin() + nr_pl0_mid , dEdx_pl0.begin() + 2*nr_pl0_mid, 0.) / nr_pl0_mid;
        }
        else {
          dEdx_pl0_mid = 0;
        }
        if(nr_pl1_mid > 0){ //pl1
          dEdx_pl1_mid = std::accumulate(dEdx_pl1.begin() + nr_pl1_mid , dEdx_pl1.begin() + 2*nr_pl1_mid, 0.) / nr_pl1_mid;
        }
        else {
          dEdx_pl1_mid = 0;
        }
        if(nr_pl2_mid > 0){ //pl2
          dEdx_pl2_mid = std::accumulate(dEdx_pl2.begin() + nr_pl2_mid , dEdx_pl2.begin() + 2*nr_pl2_mid, 0.) / nr_pl2_mid;
        }
        else {
          dEdx_pl2_mid = 0;
        }

        //-- Charge linearity
        double radius = (Trk_start_SCEcorr - Trk_end_SCEcorr).Mag();
        double sec = radius / 8;
        double min_dist = 9999;
        double charge_sec[8] = {0};
        int nr_charge[8] = {0};
        std::vector<std::vector<double>> vec_charge;
        vec_charge.resize(8);
        // divide the sphere into [8] sections
        auto assoSlice = SliceToPFPAsso.at(pfp.key()); // vector
        auto assoHits = HitToSliceAsso.at(assoSlice.front().key());

        // TODO: If hits is belong to track, use trajectory point; Otherwise use spacepoints
        for (unsigned int i_hit = 0; i_hit < assoHits.size(); i_hit++){
          auto assoSpacePoints = SpacePointToHitAsso.at(assoHits[i_hit].key());
          if(assoSpacePoints.size()==1){
            TVector3 SP = (TVector3)assoSpacePoints.front()->XYZ();
            auto SP_offset = SCE->GetCalPosOffsets(geo::Point_t(SP.X(), SP.Y(), SP.Z()));
            TVector3 SP_SCEcorr;
            SP_SCEcorr.SetX(SP.X() - SP_offset.X());
            SP_SCEcorr.SetY(SP.Y() + SP_offset.Y());
            SP_SCEcorr.SetZ(SP.Z() + SP_offset.Z());

            if ((SP_SCEcorr - Trk_start_SCEcorr).Mag() < min_dist){
              min_dist = (SP_SCEcorr - Trk_start_SCEcorr).Mag();
            }
            for(unsigned int i_sec = 0; i_sec < 8; i_sec++){
              if((SP_SCEcorr - Trk_start_SCEcorr).Mag() < sec * (i_sec + 1) && (SP_SCEcorr - Trk_start_SCEcorr).Mag() >= sec * i_sec){
                vec_charge[i_sec].push_back(assoHits[i_hit]->Integral()/100);
                charge_sec[i_sec] += assoHits[i_hit]->Integral()/100;     
                nr_charge[i_sec]++;
              }
            }
          }
          else if(assoSpacePoints.size()>1){
            std::cout<<"Number of space points: "<< assoSpacePoints.size()<<", which exceeds the normal quota."<<std::endl;
          }
        } // Finish loop the associated hits in the neutrino slice
        vtx_hit_distance = min_dist;
        double avg = 0;
        double std = 0;
        double stddev[8] = {0};
        for(int n = 0; n < 8; n++){
          avg += charge_sec[n];
        }
        avg = avg / 8.;

        for(int n = 0; n < 8; n++){
          std += (charge_sec[n] - avg) * (charge_sec[n] - avg);
        }
        for(int n = 0; n < 8; n++){
          stddev[n] = std - (charge_sec[n] - avg) * (charge_sec[n] - avg);
          stddev[n] = sqrt(1./7 * stddev[n]);
        }
        
        charge_std_bin0 = (charge_sec[0] - avg) / stddev[0];
        charge_std_bin1 = (charge_sec[1] - avg) / stddev[1];
        charge_std_bin2 = (charge_sec[2] - avg) / stddev[2];
        charge_std_bin3 = (charge_sec[3] - avg) / stddev[3];
        charge_std_bin4 = (charge_sec[4] - avg) / stddev[4];
        charge_std_bin5 = (charge_sec[5] - avg) / stddev[5];
        charge_std_bin6 = (charge_sec[6] - avg) / stddev[6];
        charge_std_bin7 = (charge_sec[7] - avg) / stddev[7];

        charge_avg_bin0 = (charge_sec[0] - avg) / avg;
        charge_avg_bin1 = (charge_sec[1] - avg) / avg;
        charge_avg_bin2 = (charge_sec[2] - avg) / avg;
        charge_avg_bin3 = (charge_sec[3] - avg) / avg;
        charge_avg_bin4 = (charge_sec[4] - avg) / avg;
        charge_avg_bin5 = (charge_sec[5] - avg) / avg;
        charge_avg_bin6 = (charge_sec[6] - avg) / avg;
        charge_avg_bin7 = (charge_sec[7] - avg) / avg;

        //-- PID
        PID pid;
        pid.Chi2(PIDTotrackAsso,daughter_Tracks.front(), Trk_start_SCEcorr, Trk_end_SCEcorr,hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2);
        PID_Chi2Mu_3pl = pid.PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
        PID_Chi2P_3pl = pid.PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID
        PID_Chi2Pi_3pl = pid.PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID
        PID_Chi2K_3pl = pid.PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

        /////////////////
        // Fill TRUE info
        //////////////// 
        if(IsMC){
          std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(daughter_Tracks.front().key());
          BackTrackerTruthMatch backtrackertruthmatch;
          backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
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
            }
          } // If MC particle exists
          reco_MC_dist_vtx = (true_nuVtx - Trk_start_SCEcorr).Mag();
          reco_MC_dist_vtx_noSCE = (true_nuVtx - Trk_start).Mag();
        }
       
        //--Directional Info
        // Check the directional info of the track by MCS
        if (bestMCSLL == fwdMCSLL) if_fwd_MCS = true;
        else if_fwd_MCS = false;

        // Check the direction info of the track by dEdx
        if (dEdx_pl2_start_half < dEdx_pl2_end_half) if_fwd_dEdxhalf = true;
        else if_fwd_dEdxhalf = false;
        if (dEdx_pl2_start1020 < dEdx_pl2_end1020) if_fwd_dEdx1020 = true;
        else if_fwd_dEdx1020 = false;

        if(IsMC){
          // Check the directional info of the track by true reco vertex distance
          TVector3 vtx(true_start_x, true_start_y, true_start_z);
          TVector3 D_start = Trk_start_SCEcorr - vtx;
          TVector3 D_end = Trk_end_SCEcorr - vtx;
          if (D_start.Mag() < D_end.Mag()) if_fwd_true = true;
          else if_fwd_true = false;
        }

      }
      
    }
  }

  my_event_->Fill();

  if(IsMC){
    MC_Primary_PDG.clear();
    MC_Primary_Mom.clear();
    MC_proton_true_Mom_above260.clear();
    MC_muon_true_Mom.clear();
    MC_muon_true_cos_theta.clear();
    MC_muon_true_phi.clear();
    MC_proton_true_cos_theta.clear();
    MC_proton_true_phi.clear();
    Ghost_PDG.clear();
 
    cos_ang_muon_proton = -999;
    dist_muon_proton = -999;
    len_Muon_0pi0p = -999;
    len_Muon_0pi1p = -999;
    len_Proton_0pi1p = -999;

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

  crthit_PE.clear();
  crthit_plane.clear();
  crthit_time.clear();
  Nr_crthit_inBeam = 0;
  Nr_trk_asso_crthit = 0;
  trk_crt_time = -999;
  evt_CRTveto = false;
  evt_CRTveto_100 = false;

  trk_vtx_x = -999;
  trk_vtx_y = -999;
  trk_vtx_z = -999;
  trk_start_x = -999;
  trk_start_y = -999;
  trk_start_z = -999;
  trk_end_x = -999;
  trk_end_y = -999;
  trk_end_z = -999;

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

  bestMCSLL_NoSCE = -999;
  fwdMCSLL_NoSCE = -999;
  bwdMCSLL_NoSCE = -999;

  missing_PT_range = -999;
  missing_PT_MCS = -999;

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

  avg_dEdx_LargeHit_pl0 = -999;
  avg_dEdx_LargeHit_pl1 = -999;
  avg_dEdx_LargeHit_pl2 = -999;

  dEdx_pl0_mid = -999;
  dEdx_pl1_mid = -999;
  dEdx_pl2_mid = -999;

  charge_std_bin0 = -999; // the multiplication of charge in bin 0 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin1 = -999; // the multiplication of charge in bin 1 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin2 = -999; // the multiplication of charge in bin 2 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin3 = -999; // the multiplication of charge in bin 3 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin4 = -999; // the multiplication of charge in bin 4 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin5 = -999; // the multiplication of charge in bin 5 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin6 = -999; // the multiplication of charge in bin 6 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_std_bin7 = -999; // the multiplication of charge in bin 7 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections

  charge_avg_bin0 = -999; // the multiplication of charge in bin 0 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin1 = -999; // the multiplication of charge in bin 1 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin2 = -999; // the multiplication of charge in bin 2 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin3 = -999; // the multiplication of charge in bin 3 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin4 = -999; // the multiplication of charge in bin 4 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin5 = -999; // the multiplication of charge in bin 5 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin6 = -999; // the multiplication of charge in bin 6 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections
  charge_avg_bin7 = -999; // the multiplication of charge in bin 7 diff to the avg of 8 sections in the slice to the standard deviation of charge in 8 sections

  vtx_hit_distance = -999;// Distance of track vertex and the closest hit spacepoints

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
  
  if(IsMC){
    my_event_->Branch("EventWeight", &EventWeight);
    my_event_->Branch("TopologyType", &TopologyType);
    my_event_->Branch("cos_ang_muon_proton", &cos_ang_muon_proton);
    my_event_->Branch("dist_muon_proton", &dist_muon_proton);
    my_event_->Branch("len_Muon_0pi0p", &len_Muon_0pi0p);
    my_event_->Branch("len_Muon_0pi1p", &len_Muon_0pi1p);
    my_event_->Branch("len_Proton_0pi1p", &len_Proton_0pi1p);
    my_event_->Branch("MC_beamNeutrino", &MC_beamNeutrino);
    my_event_->Branch("MC_int_mode", &MC_int_mode);
    my_event_->Branch("MC_nupdg", &MC_nupdg);
    my_event_->Branch("MC_ccnc", &MC_ccnc);
    my_event_->Branch("MC_nu_E", &MC_nu_E);
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
    my_event_->Branch("MC_nPiPlus_below65", &MC_nPiPlus_below65);
    my_event_->Branch("MC_nPiPlus_above65", &MC_nPiPlus_above65);
    my_event_->Branch("MC_nPiMinus_below65", &MC_nPiMinus_below65);
    my_event_->Branch("MC_nPiMinus_above65", &MC_nPiMinus_above65);
    my_event_->Branch("MC_Primary_PDG", &MC_Primary_PDG);
    my_event_->Branch("MC_Primary_Mom", &MC_Primary_Mom);
    my_event_->Branch("MC_proton_true_Mom_above260", &MC_proton_true_Mom_above260);
    my_event_->Branch("MC_muon_true_Mom", &MC_muon_true_Mom);
    my_event_->Branch("MC_muon_true_cos_theta", &MC_muon_true_cos_theta);
    my_event_->Branch("MC_muon_true_phi", &MC_muon_true_phi);
    my_event_->Branch("MC_proton_true_cos_theta", &MC_proton_true_cos_theta);
    my_event_->Branch("MC_proton_true_phi", &MC_proton_true_phi);
    my_event_->Branch("Ghost_PDG", &Ghost_PDG);

    my_event_->Branch("Genie_Q2", &Genie_Q2);
    my_event_->Branch("Genie_q2", &Genie_q2);
    my_event_->Branch("Genie_W", &Genie_W);
    my_event_->Branch("Genie_nNeutron_preFSI", &Genie_nNeutron_preFSI);
    my_event_->Branch("Genie_nProton_preFSI", &Genie_nProton_preFSI);
    my_event_->Branch("Genie_nPi0_preFSI", &Genie_nPi0_preFSI);
    my_event_->Branch("Genie_nPiPlus_preFSI", &Genie_nPiPlus_preFSI);
    my_event_->Branch("Genie_nPiMinus_preFSI", &Genie_nPiMinus_preFSI);

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
  }

  my_event_->Branch("n_pfp_nuDaughters", &n_pfp_nuDaughters);
  my_event_->Branch("n_dau_tracks", &n_dau_tracks);
  my_event_->Branch("n_dau_showers", &n_dau_showers);
  
  my_event_->Branch("nr_granddau_shw", &nr_granddau_shw);
  my_event_->Branch("nr_granddau_trk", &nr_granddau_trk);
  my_event_->Branch("nr_granddau", &nr_granddau);
  my_event_->Branch("granddau_trk_len", &granddau_trk_len);
  my_event_->Branch("granddau_shw_len", &granddau_shw_len);
  
  my_event_->Branch("flash_matching_chi2", &flash_matching_chi2);
  my_event_->Branch("flash_YCenter", &flash_YCenter);
  my_event_->Branch("flash_YWidth", &flash_YWidth);
  my_event_->Branch("flash_ZCenter", &flash_ZCenter);
  my_event_->Branch("flash_ZWidth", &flash_ZWidth);
  my_event_->Branch("flash_TotalPE", &flash_TotalPE);

  my_event_->Branch("evt_CRTveto", &evt_CRTveto);
  my_event_->Branch("evt_CRTveto_100", &evt_CRTveto_100);
  my_event_->Branch("crthit_PE", &crthit_PE);
  my_event_->Branch("crthit_plane", &crthit_plane);
  my_event_->Branch("crthit_time", &crthit_time);
  my_event_->Branch("Nr_crthit_inBeam", &Nr_crthit_inBeam);
  my_event_->Branch("Nr_trk_asso_crthit", &Nr_trk_asso_crthit);
  my_event_->Branch("trk_crt_time", &trk_crt_time);
  
  my_event_->Branch("trk_vtx_x", &trk_vtx_x);
  my_event_->Branch("trk_vtx_y", &trk_vtx_y);
  my_event_->Branch("trk_vtx_z", &trk_vtx_z);
  my_event_->Branch("trk_start_x", &trk_start_x);
  my_event_->Branch("trk_start_y", &trk_start_y);
  my_event_->Branch("trk_start_z", &trk_start_z);
  my_event_->Branch("trk_end_x", &trk_end_x);
  my_event_->Branch("trk_end_y", &trk_end_y);
  my_event_->Branch("trk_end_z", &trk_end_z);

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
  my_event_->Branch("bestMCSLL_NoSCE", &bestMCSLL_NoSCE);
  my_event_->Branch("fwdMCSLL_NoSCE", &fwdMCSLL_NoSCE);
  my_event_->Branch("bwdMCSLL_NoSCE", &bwdMCSLL_NoSCE);

  my_event_->Branch("missing_PT_range", &missing_PT_range);
  my_event_->Branch("missing_PT_MCS", &missing_PT_MCS);

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

  my_event_->Branch("avg_dEdx_LargeHit_pl0", &avg_dEdx_LargeHit_pl0);
  my_event_->Branch("avg_dEdx_LargeHit_pl1", &avg_dEdx_LargeHit_pl1);
  my_event_->Branch("avg_dEdx_LargeHit_pl2", &avg_dEdx_LargeHit_pl2);

  my_event_->Branch("dEdx_pl0_mid", &dEdx_pl0_mid);
  my_event_->Branch("dEdx_pl1_mid", &dEdx_pl1_mid);
  my_event_->Branch("dEdx_pl2_mid", &dEdx_pl2_mid);

  my_event_->Branch("charge_std_bin0", &charge_std_bin0);
  my_event_->Branch("charge_std_bin1", &charge_std_bin1);
  my_event_->Branch("charge_std_bin2", &charge_std_bin2);
  my_event_->Branch("charge_std_bin3", &charge_std_bin3);
  my_event_->Branch("charge_std_bin4", &charge_std_bin4);
  my_event_->Branch("charge_std_bin5", &charge_std_bin5);
  my_event_->Branch("charge_std_bin6", &charge_std_bin6);
  my_event_->Branch("charge_std_bin7", &charge_std_bin7);

  my_event_->Branch("charge_avg_bin0", &charge_avg_bin0);
  my_event_->Branch("charge_avg_bin1", &charge_avg_bin1);
  my_event_->Branch("charge_avg_bin2", &charge_avg_bin2);
  my_event_->Branch("charge_avg_bin3", &charge_avg_bin3);
  my_event_->Branch("charge_avg_bin4", &charge_avg_bin4);
  my_event_->Branch("charge_avg_bin5", &charge_avg_bin5);
  my_event_->Branch("charge_avg_bin6", &charge_avg_bin6);
  my_event_->Branch("charge_avg_bin7", &charge_avg_bin7);

  my_event_->Branch("vtx_hit_distance", &vtx_hit_distance);

  my_event_->Branch("if_fwd_true", &if_fwd_true);
  my_event_->Branch("if_fwd_MCS", &if_fwd_MCS);
  my_event_->Branch("if_fwd_dEdx1020", &if_fwd_dEdx1020);
  my_event_->Branch("if_fwd_dEdxhalf", &if_fwd_dEdxhalf);

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

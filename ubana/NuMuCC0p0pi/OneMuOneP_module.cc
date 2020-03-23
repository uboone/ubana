////////////////////////////////////////////////////////////////////////
// Class:       OneMuOneP
// Plugin Type: analyzer (art v3_01_02)
// File:        OneMuOneP_module.cc
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
  int MC_nProton_below255; // Number of proton(s) (p<255) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_above255; // Number of proton(s) (p >= 255) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPi0; // Number of pi0(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus_below65; // Number of pi plus(s) (p < 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus_above65; // Number of pi plus(s) (p > 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus_below65; // Number of pi minus(s) (p < 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus_above65; // Number of pi minus(s) (p > 65MeV) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  std::vector<int> MC_Primary_PDG; // PDG of neutrino daughters
  std::vector<double> MC_Primary_Mom; // Momemtum of neutrino daughters
  std::vector<double> MC_proton_true_Mom_above255; // Momentum of proton above 255 MeV
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

  double flash_matching_chi2 = -999; //Chi2 of flash matching in each neutrino slice
  double flash_YCenter = -999;
  double flash_YWidth = -999;
  double flash_ZCenter = -999;
  double flash_ZWidth = -999;
  double flash_TotalPE = -999;

  int n_pfp_nuDaughters;
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters

  bool if_2tracks = false; // if there are two primary tracks and nothing

  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime
  
  std::vector<int> Nr_trk_asso_crthit(2, -999); // Number of CRT hits associated to the track
  std::vector<double> trk_crt_time(2, -999);// The CRT time of the hit which matched to the track
  std::vector<bool> if_trk_CRT_out_Beam(2, false); // Check if a track matches with out of beam CRT hit(s)

  std::vector<bool> trk_OutOfTime(2, false);

  std::vector<bool> vtx_InFV(2, false);
  std::vector<bool> trk_contained(2, false);

  std::vector<double> trk_phi(2, -999);
  std::vector<double> trk_theta(2, -999);
  std::vector<double> trk_costheta(2, -999);

  std::vector<double> theta_pl2(2, -999);
  std::vector<double> theta_pl1(2, -999);
  std::vector<double> theta_pl0(2, -999);
  std::vector<double> phi_readout(2, -999);

  std::vector<double> sin2_theta_pl2(2, -999);
  std::vector<double> sin2_theta_pl1(2, -999);
  std::vector<double> sin2_theta_pl0(2, -999);
  std::vector<double> sin2_phi_readout(2, -999);

  std::vector<double> bestMCS(2, -999);
  std::vector<double> bestMCSLL(2, -999);
  std::vector<double> fwdMCS(2, -999);
  std::vector<double> fwdMCSLL(2, -999);
  std::vector<double> bwdMCS(2, -999);
  std::vector<double> bwdMCSLL(2, -999);

  std::vector<double> bestMCSLL_NoSCE(2, -999);
  std::vector<double> fwdMCSLL_NoSCE(2, -999);
  std::vector<double> bwdMCSLL_NoSCE(2, -999);

  std::vector<double> Trk_length(2, -999);
  std::vector<double> Mom_Range_mu(2, -999);
  std::vector<double> Mom_Range_p(2, -999);
  std::vector<double> Mom_Range_pi(2, -999);
 
  std::vector<float> dEdx_pl0; // dE/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<float> dEdx_pl1; // dE/dx of the selected (muon) track from plane 1
  std::vector<float> dEdx_pl2; // dE/dx of the selected (muon) track from plane 2 (collection)
  std::vector<float> resRange_pl0; // range from a hit to the end of the selected track end
  std::vector<float> resRange_pl1; // range from a hit to the end of the selected track end
  std::vector<float> resRange_pl2; // range from a hit to the end of the selected track end 

  std::vector<double> avg_dEdx_LargeHit_pl0(2, -999);
  std::vector<double> avg_dEdx_LargeHit_pl1(2, -999);
  std::vector<double> avg_dEdx_LargeHit_pl2(2, -999);

  std::vector<double> dEdx_pl0_start_half(2, -999);
  std::vector<double> dEdx_pl0_end_half(2, -999);
  std::vector<double> dEdx_pl1_start_half(2, -999);
  std::vector<double> dEdx_pl1_end_half(2, -999);
  std::vector<double> dEdx_pl2_start_half(2, -999);
  std::vector<double> dEdx_pl2_end_half(2, -999);

  std::vector<double> dEdx_pl0_start1020(2, -999);
  std::vector<double> dEdx_pl0_end1020(2, -999);
  std::vector<double> dEdx_pl1_start1020(2, -999);
  std::vector<double> dEdx_pl1_end1020(2, -999);
  std::vector<double> dEdx_pl2_start1020(2, -999);
  std::vector<double> dEdx_pl2_end1020(2, -999);

  std::vector<double> dEdx_pl2_1020_ratio(2, -999);
  std::vector<double> dEdx_pl2_half_ratio(2, -999);

  std::vector<double> dEdx_pl0_mid(2, -999);
  std::vector<double> dEdx_pl1_mid(2, -999);
  std::vector<double> dEdx_pl2_mid(2, -999);

  std::vector<double> PID_Chi2Mu_3pl(2, -999);
  std::vector<double> PID_Chi2P_3pl(2, -999);
  std::vector<double> PID_Chi2Pi_3pl(2, -999);
  std::vector<double> PID_Chi2K_3pl(2, -999);

  std::vector<double> trk_cosmic_percent(2, -999);
  std::vector<double> trk_purity(2, -999);
  std::vector<double> trk_completeness(2, -999);

  std::vector<bool> if_cosmic(2, true);
  std::vector<bool> if_matchMu(2, false);
  
  std::vector<double> true_mom(2, -999);
  std::vector<double> true_start_x(2, -999);
  std::vector<double> true_start_y(2, -999);
  std::vector<double> true_start_z(2, -999);
  std::vector<double> true_end_x(2, -999);
  std::vector<double> true_end_y(2, -999);
  std::vector<double> true_end_z(2, -999);
  std::vector<double> true_trk_phi(2, -999);
  std::vector<double> true_trk_theta(2, -999);
  std::vector<double> true_trk_costheta(2, -999);
  std::vector<double> true_trk_theta_yz(2, -999);
  std::vector<double> true_trk_costheta_yz(2, -999);
  std::vector<double> true_trk_theta_xz(2, -999);
  std::vector<double> true_trk_costheta_xz(2, -999);
  std::vector<double> true_trk_length(2, -999);
  std::vector<double> true_trk_PDG(2, -999);
  std::vector<bool> true_trk_ifcontained(2, false);
  std::vector<bool> true_vtxFV(2, false);

  int muon_idx = -999;
  int proton_idx = -999;

  double vtx_x = -999;
  double vtx_y = -999;
  double vtx_z = -999;

  double muon_cand_phi = -999;
  double muon_cand_theta = -999;
  double muon_cand_costheta = -999;
  double muon_cand_mom_MCS = -999;
  double muon_cand_mom_range = -999;

  double proton_cand_phi = -999;
  double proton_cand_theta = -999;
  double proton_cand_costheta = -999;
  double proton_cand_mom_range = -999;

  double missing_PT = -999;

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
  std::vector<TVector3> Trk_start(2);
  std::vector<TVector3> Trk_end(2);
  std::vector<TVector3> Trk_start_SCEcorr(2); 
  std::vector<TVector3> Trk_end_SCEcorr(2);
 
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
    MC_nProton_below255 = 0;
    MC_nProton_above255 = 0;
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
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() < 0.255) MC_nProton_below255++;
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() >= 0.255) {
            MC_nProton_above255++; 
            MC_proton_true_Mom_above255.push_back(MCParticleCollection[i_mcp]->P());
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
    TopologyType = topology.TopologyLabel(MC_nMuon, MC_nElectron, MC_nPiPlus_above65, MC_nPiPlus_below65, MC_nPiMinus_above65, MC_nPiMinus_below65, MC_nPi0, MC_nProton_above255, MC_nProton_below255, MC_nupdg, MC_ccnc, MC_beamNeutrino, MC_FV);
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

      // Get Neutrino slice flash matching score (chi2)
      auto flash_matching_T0 = pfpToT0Asso.at(pfp.key());

      if(flash_matching_T0.size() == 1){
        flash_matching_chi2 = flash_matching_T0.front()->TriggerConfidence();
      }
    
      if(flash_v.size() == 1){
        flash_YCenter = flash_v[0]->YCenter();
        flash_YWidth = flash_v[0]->YWidth();
        flash_ZCenter = flash_v[0]->ZCenter();
        flash_ZWidth = flash_v[0]->ZWidth();
        flash_TotalPE = flash_v[0]->TotalPE();
      }

      // For CC0pi0p, we only consider the case with the number of neutrino daughters less than 6
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
        } // finish looping of daughter pfp
      } // if nr_pfp daughter < 6
     
      //number of tracks and showers
      n_dau_tracks = daughter_Tracks.size();
      n_dau_showers = daughter_Showers.size();

      ////////////////////
      // Selection info
      ////////////////////
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
        if(UsingCRT){
          // overlay is also basically data, using ts0
          for (unsigned int i_crt = 0; i_crt < crthit_v.size(); i_crt ++){
            // figure out what plane this hit comes from
            // 3 -> top, 0 -> bottom, 1 -> anode, 2 -> cathode
            double crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
            if(crt_time >= fBeamStart && crt_time <= fBeamEnd){
              crthit_PE.push_back(crthit_v[i_crt]->peshit);
              crthit_plane.push_back(crthit_v[i_crt]->plane);
              crthit_time.push_back(crt_time);
              Nr_crthit_inBeam++; // The total number of CRT hits in beam should be less than the number of exiting tracks
            } // If CRT hit in Beam window
          } // CRT hit loop
          if(Nr_crthit_inBeam != (int) crthit_time.size()) {
            throw cet::exception("[Numu0pi0p]") << "Number of CRT hits in beam time does not match!" << std::endl;
          }
        }

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
            Nr_trk_asso_crthit[i_trk] = Track_CRThit.size();
            if(Nr_trk_asso_crthit[i_trk] > 0){
              for(unsigned int i_trk_hit = 0; i_trk_hit < Track_CRThit.size(); i_trk_hit++){
                trk_crt_time[i_trk] = ((Track_CRThit[i_trk_hit]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
                if(trk_crt_time[i_trk] < fBeamStart || trk_crt_time[i_trk] > fBeamEnd){
                  // For 2 tracks, as long as one of them has associated CRT hit out of time, reject!
                  if_trk_CRT_out_Beam[i_trk] = true;
                } // If matched CRT hit out of Beam window
              } // Loop over the associated CRT hits
            }  // If there is associated CRT hits
          } // CRT info

          //------ Reco info

          //-- Track start and end (The track start of the muon candidate track will be the vertex)
          Trk_start[i_trk] = daughter_Tracks[i_trk]->Start<TVector3>();
          auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start[i_trk].X(), Trk_start[i_trk].Y(), Trk_start[i_trk].Z()));
          Trk_start_SCEcorr[i_trk].SetX(Trk_start[i_trk].X() - Trk_start_offset.X());
          Trk_start_SCEcorr[i_trk].SetY(Trk_start[i_trk].Y() + Trk_start_offset.Y());
          Trk_start_SCEcorr[i_trk].SetZ(Trk_start[i_trk].Z() + Trk_start_offset.Z());

          Trk_end[i_trk] = daughter_Tracks[i_trk]->End<TVector3>();
          auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end[i_trk].X(), Trk_end[i_trk].Y(), Trk_end[i_trk].Z()));
          Trk_end_SCEcorr[i_trk].SetX(Trk_end[i_trk].X() - Trk_end_offset.X());
          Trk_end_SCEcorr[i_trk].SetY(Trk_end[i_trk].Y() + Trk_end_offset.Y());
          Trk_end_SCEcorr[i_trk].SetZ(Trk_end[i_trk].Z() + Trk_end_offset.Z());

          //-- Track out of time
          if( Trk_start_SCEcorr[i_trk].X() < 0 || Trk_start_SCEcorr[i_trk].X() > 2. * geo->DetHalfWidth() || Trk_end_SCEcorr[i_trk].X() < 0 || Trk_end_SCEcorr[i_trk].X() > 2. * geo->DetHalfWidth()){
            trk_OutOfTime[i_trk] = true;
          }
          else{
            trk_OutOfTime[i_trk] = false;
          }

          //-- If contained
          vtx_InFV[i_trk] = _fiducial_volume.VertexInFV(Trk_start_SCEcorr[i_trk]);
          trk_contained[i_trk] = _fiducial_volume.TrackContain(Trk_start_SCEcorr[i_trk], Trk_end_SCEcorr[i_trk]);

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
        
          //- Get the three largest hits (dE/dx)
          auto copy_dEdx_pl0 = dEdx_pl0;
          auto copy_dEdx_pl1 = dEdx_pl1;
          auto copy_dEdx_pl2 = dEdx_pl2;

          std::sort(copy_dEdx_pl0.begin(), copy_dEdx_pl0.end());
          std::sort(copy_dEdx_pl1.begin(), copy_dEdx_pl1.end());
          std::sort(copy_dEdx_pl2.begin(), copy_dEdx_pl2.end());  
          // pl 0
          if(hits_dEdx_size_pl0 < 3){
            for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl0; i_hit++){
              avg_dEdx_LargeHit_pl0[i_trk] += copy_dEdx_pl0[i_hit];
            }
            avg_dEdx_LargeHit_pl0[i_trk] = avg_dEdx_LargeHit_pl0[i_trk] / hits_dEdx_size_pl0;
          }
          else{
            avg_dEdx_LargeHit_pl0[i_trk] = (copy_dEdx_pl0[hits_dEdx_size_pl0 - 1] + copy_dEdx_pl0[hits_dEdx_size_pl0 - 2] + copy_dEdx_pl0[hits_dEdx_size_pl0 - 3]) / 3;
          }
          // pl 1
          if(hits_dEdx_size_pl1 < 3){
            for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl1; i_hit++){
              avg_dEdx_LargeHit_pl1[i_trk] += copy_dEdx_pl1[i_hit];
            }
            avg_dEdx_LargeHit_pl1[i_trk] = avg_dEdx_LargeHit_pl1[i_trk] / hits_dEdx_size_pl1;
          }
          else{
            avg_dEdx_LargeHit_pl1[i_trk] = (copy_dEdx_pl1[hits_dEdx_size_pl1 - 1] + copy_dEdx_pl1[hits_dEdx_size_pl1 - 2] + copy_dEdx_pl1[hits_dEdx_size_pl1 - 3]) / 3;
          }
          // pl 2
          if(hits_dEdx_size_pl2 < 3){
            for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl2; i_hit++){
              avg_dEdx_LargeHit_pl2[i_trk] += copy_dEdx_pl2[i_hit];
            }
            avg_dEdx_LargeHit_pl2[i_trk] = avg_dEdx_LargeHit_pl2[i_trk] / hits_dEdx_size_pl2;
          }
          else{
            avg_dEdx_LargeHit_pl2[i_trk] = (copy_dEdx_pl2[hits_dEdx_size_pl2 - 1] + copy_dEdx_pl2[hits_dEdx_size_pl2 - 2] + copy_dEdx_pl2[hits_dEdx_size_pl2 - 3]) / 3;
          }

          //- Average dE/dx close to the start and the end of the track
          dEdx_pl0_start_half[i_trk] = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
          dEdx_pl0_end_half[i_trk] = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
          dEdx_pl1_start_half[i_trk] = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
          dEdx_pl1_end_half[i_trk] = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
          dEdx_pl2_start_half[i_trk] = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
          dEdx_pl2_end_half[i_trk] = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;

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

          dEdx_pl2_1020_ratio[i_trk] = dEdx_pl2_end1020[i_trk] / (dEdx_pl2_end1020[i_trk] + dEdx_pl2_start1020[i_trk]);
          dEdx_pl2_half_ratio[i_trk] = dEdx_pl2_end_half[i_trk] / (dEdx_pl2_end_half[i_trk] + dEdx_pl2_start_half[i_trk]);

          //- Get dEdx of the middle part in the track
          int nr_pl0_mid = dEdx_pl0.size()/3;
          int nr_pl1_mid = dEdx_pl1.size()/3;
          int nr_pl2_mid = dEdx_pl2.size()/3;
          if(nr_pl0_mid > 0){ //pl0
            dEdx_pl0_mid[i_trk] = std::accumulate(dEdx_pl0.begin() + nr_pl0_mid , dEdx_pl0.begin() + 2*nr_pl0_mid, 0.) / nr_pl0_mid;
          }
          else {
            dEdx_pl0_mid[i_trk] = 0;
          }
          if(nr_pl1_mid > 0){ //pl1
            dEdx_pl1_mid[i_trk] = std::accumulate(dEdx_pl1.begin() + nr_pl1_mid , dEdx_pl1.begin() + 2*nr_pl1_mid, 0.) / nr_pl1_mid;
          }
          else {
            dEdx_pl1_mid[i_trk] = 0;
          }
          if(nr_pl2_mid > 0){ //pl2
            dEdx_pl2_mid[i_trk] = std::accumulate(dEdx_pl2.begin() + nr_pl2_mid , dEdx_pl2.begin() + 2*nr_pl2_mid, 0.) / nr_pl2_mid;
          }
          else {
            dEdx_pl2_mid[i_trk] = 0;
          }

          //--PID
          PID pid;
          pid.Chi2(PIDTotrackAsso,daughter_Tracks[i_trk], Trk_start_SCEcorr[i_trk], Trk_end_SCEcorr[i_trk],hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2);

          PID_Chi2Mu_3pl[i_trk] = pid.PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
          PID_Chi2P_3pl[i_trk] = pid.PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID
          PID_Chi2Pi_3pl[i_trk] = pid.PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID
          PID_Chi2K_3pl[i_trk] = pid.PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

          /////////////////
          // Fill True info from reco-truth matching
          ////////////////
          // The corresponding MC particle is not necessary from numu, cc and primary process at the moment
          if(IsMC){
            std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(daughter_Tracks[i_trk].key());
            BackTrackerTruthMatch backtrackertruthmatch;
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
              if_cosmic[i_trk] = false;
              if(MCparticle->PdgCode() == 13) if_matchMu[i_trk] = true;
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
            }
          }
          
        } // Finish looping the 2 tracks

        ///////////////////////
        // Decide the vertex position, missing PT, which track is proton and which is muon
        //////////////////////
        //
        //-- Track particle ID (Set the threshold to be PID_Chi2P_3pl 95)
        // If we don't use track length for the moment
        if(PID_Chi2P_3pl[0] <= 95 && PID_Chi2P_3pl[1] > 95){
          muon_idx = 0;
          proton_idx = 1;
        }
        if(PID_Chi2P_3pl[0] > 95 && PID_Chi2P_3pl[1] <= 95){
          muon_idx = 1;
          proton_idx = 0;
        }
        if(PID_Chi2P_3pl[0] <= 95 && PID_Chi2P_3pl[1] <= 95){
          if(dEdx_pl2_mid[0] >= dEdx_pl2_mid[1]){
            if(dEdx_pl2_mid[0] > 3.5){
              muon_idx = 1;
              proton_idx = 0;
            }
            else{
              muon_idx = -1;
              proton_idx = -1;
            }
          }

          if(dEdx_pl2_mid[0] < dEdx_pl2_mid[1]){
            if(dEdx_pl2_mid[1] > 3.5){
              muon_idx = 0;
              proton_idx = 1;
            }
            else{
              muon_idx = -1;
              proton_idx = -1;
            }
          }
        }
        if(PID_Chi2P_3pl[0] > 95 && PID_Chi2P_3pl[1] > 95){
          // Maybe check the muon PID and decide the
          if(dEdx_pl2_mid[0] >= dEdx_pl2_mid[1]){
            if(dEdx_pl2_mid[0] > 3.5){
              muon_idx = 1;
              proton_idx = 0;
            }
            else{
              muon_idx = -1;
              proton_idx = -1;
            }
          }

          if(dEdx_pl2_mid[0] < dEdx_pl2_mid[1]){
            if(dEdx_pl2_mid[1] > 3.5){
              muon_idx = 0;
              proton_idx = 1;
            }
            else{
              muon_idx = -1;
              proton_idx = -1;
            }
          }
        }

        //-- Vertex and Others
        if(muon_idx >= 0 && proton_idx >= 0){
          vtx_x = Trk_start_SCEcorr[muon_idx].X();
          vtx_y = Trk_start_SCEcorr[muon_idx].Y();
          vtx_z = Trk_start_SCEcorr[muon_idx].Z();

          muon_cand_phi = daughter_Tracks[muon_idx]->Phi();
          muon_cand_theta = daughter_Tracks[muon_idx]->Theta();
          muon_cand_costheta = cos(daughter_Tracks[muon_idx]->Theta());
          muon_cand_mom_MCS = bestMCS[muon_idx];
          muon_cand_mom_range = Mom_Range_mu[muon_idx];

          proton_cand_phi = daughter_Tracks[proton_idx]->Phi();
          proton_cand_theta = daughter_Tracks[proton_idx]->Theta();
          proton_cand_costheta = cos(daughter_Tracks[proton_idx]->Theta());
          proton_cand_mom_range = Mom_Range_p[proton_idx];
        }
        
        //-- Missing PT (no need to decide the PID of the tracks)
        // PT = PT(track 1) + PT(track 2) = sqrt((PT1_x + PT2_x)^2 + (PT1_y + PT2_y)^2)
        if(muon_idx >= 0 && proton_idx >= 0){
          float muon_PT_x = bestMCS[muon_idx] * sin(daughter_Tracks[muon_idx]->Theta()) * cos(daughter_Tracks[muon_idx]->Phi());
          float muon_PT_y = bestMCS[muon_idx] * sin(daughter_Tracks[muon_idx]->Theta()) * sin(daughter_Tracks[muon_idx]->Phi());
          float proton_PT_x = Mom_Range_p[proton_idx] * sin(daughter_Tracks[proton_idx]->Theta()) * cos(daughter_Tracks[proton_idx]->Phi());
          float proton_PT_y = Mom_Range_p[proton_idx] * sin(daughter_Tracks[proton_idx]->Theta()) * sin(daughter_Tracks[proton_idx]->Phi());

          missing_PT = sqrt((muon_PT_x + proton_PT_x) * (muon_PT_x + proton_PT_x) + (muon_PT_y + proton_PT_y) * (muon_PT_y + proton_PT_y));
        }

      } // If it is 2-primary-track selection
    } // If it is neutrino slice
  } // End of Loop over all pfparticles

  my_event_->Fill();

  flash_matching_chi2 = -999; 
  flash_YCenter = -999;
  flash_YWidth = -999;
  flash_ZCenter = -999;
  flash_ZWidth = -999;
  flash_TotalPE = -999;

  if_2tracks = false;

  crthit_PE.clear();
  crthit_plane.clear();
  crthit_time.clear();
  Nr_crthit_inBeam = 0;

  std::fill_n(Nr_trk_asso_crthit, 2, -999); 
  std::fill_n(trk_crt_time, 2, -999);
  std::fill_n(if_trk_CRT_out_Beam, 2, false);


  std::fill_n(trk_OutOfTime, 2, false);

  std::fill_n(vtx_InFV, 2, false);
  std::fill_n(trk_contained, 2, false);

  std::fill_n(trk_phi, 2, -999);
  std::fill_n(trk_theta, 2, -999);
  std::fill_n(trk_costheta, 2, -999);

  std::fill_n(theta_pl0, 2, -999);
  std::fill_n(theta_pl1, 2, -999);
  std::fill_n(theta_pl2, 2, -999);
  std::fill_n(phi_readout, 2, -999);

  std::fill_n(sin2_theta_pl0, 2, -999);
  std::fill_n(sin2_theta_pl1, 2, -999);
  std::fill_n(sin2_theta_pl2, 2, -999);
  std::fill_n(sin2_phi_readout, 2, -999);

  std::fill_n(bestMCS, 2, -999);
  std::fill_n(bestMCSLL, 2, -999);
  std::fill_n(fwdMCS, 2, -999);
  std::fill_n(fwdMCSLL, 2, -999);
  std::fill_n(bwdMCS, 2, -999);
  std::fill_n(bwdMCSLL, 2, -999);

  std::fill_n(bestMCSLL_NoSCE, 2, -999);
  std::fill_n(fwdMCSLL_NoSCE, 2, -999);
  std::fill_n(bwdMCSLL_NoSCE, 2, -999);

  
  std::fill_n(Trk_length, 2, -999);
  std::fill_n(Mom_Range_mu, 2, -999);
  std::fill_n(Mom_Range_p, 2, -999);
  std::fill_n(Mom_Range_pi, 2, -999);
 
  dEdx_pl0.clear(); 
  dEdx_pl1.clear(); 
  dEdx_pl2.clear(); 
  resRange_pl0.clear();
  resRange_pl1.clear();
  resRange_pl2.clear();

  std::fill_n(avg_dEdx_LargeHit_pl0, 2, -999);
  std::fill_n(avg_dEdx_LargeHit_pl1, 2, -999);
  std::fill_n(avg_dEdx_LargeHit_pl2, 2, -999);
 
  std::fill_n(dEdx_pl0_start_half, 2, -999);
  std::fill_n(dEdx_pl0_end_half, 2, -999);
  std::fill_n(dEdx_pl1_start_half, 2, -999);
  std::fill_n(dEdx_pl1_end_half, 2, -999);
  std::fill_n(dEdx_pl2_start_half, 2, -999);
  std::fill_n(dEdx_pl2_end_half, 2, -999);

  std::fill_n(dEdx_pl0_start1020, 2, -999);
  std::fill_n(dEdx_pl0_end1020, 2, -999);
  std::fill_n(dEdx_pl1_start1020, 2, -999);
  std::fill_n(dEdx_pl1_end1020, 2, -999);
  std::fill_n(dEdx_pl2_start1020, 2, -999);
  std::fill_n(dEdx_pl2_end1020, 2, -999);

  std::fill_n(dEdx_pl0_mid, 2, -999);
  std::fill_n(dEdx_pl1_mid, 2, -999);
  std::fill_n(dEdx_pl2_mid, 2, -999);

  std::fill_n(PID_Chi2Mu_3pl, 2, -999);
  std::fill_n(PID_Chi2P_3pl, 2, -999);
  std::fill_n(PID_Chi2Pi_3pl, 2, -999);
  std::fill_n(PID_Chi2K_3pl, 2, -999);

  std::fill_n(trk_cosmic_percent, 2, -999);
  std::fill_n(trk_purity, 2, -999);
  std::fill_n(trk_completeness, 2, -999);

  std::fill_n(if_matchMu, 2, false);
  std::fill_n(if_cosmic, 2, true);

  std::fill_n(true_mom, 2, -999);  
  std::fill_n(true_start_x, 2, -999);  
  std::fill_n(true_start_y, 2, -999);  
  std::fill_n(true_start_z, 2, -999);
  std::fill_n(true_end_x, 2, -999);
  std::fill_n(true_end_y, 2, -999);
  std::fill_n(true_end_z, 2, -999);  
  std::fill_n(true_trk_ifcontained, 2, false);  
  std::fill_n(true_vtxFV, 2, false);  
  std::fill_n(true_trk_phi, 2, -999);  
  std::fill_n(true_trk_theta, 2, -999);  
  std::fill_n(true_trk_costheta, 2, -999);  
  std::fill_n(true_trk_theta_yz, 2, -999);  
  std::fill_n(true_trk_costheta_yz, 2, -999);  
  std::fill_n(true_trk_theta_xz, 2, -999);  
  std::fill_n(true_trk_costheta_xz, 2, -999);  
  std::fill_n(true_trk_length, 2, -999);  
  std::fill_n(true_trk_PDG, 2, -999);  

  muon_idx = -999;
  proton_idx = -999;

  vtx_x = -999;
  vtx_y = -999;
  vtx_z = -999;

  muon_cand_phi = -999;
  muon_cand_theta = -999;
  muon_cand_costheta = -999;
  muon_cand_mom_MCS = -999;
  muon_cand_mom_range = -999;

  proton_cand_phi = -999;
  proton_cand_theta = -999;
  proton_cand_costheta = -999;
  proton_cand_mom_range = -999;

  missing_PT = -999;

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
    my_event_->Branch("MC_nProton_below255", &MC_nProton_below255);
    my_event_->Branch("MC_nProton_above255", &MC_nProton_above255);
    my_event_->Branch("MC_nPi0", &MC_nPi0);
    my_event_->Branch("MC_nPiPlus_below65", &MC_nPiPlus_below65);
    my_event_->Branch("MC_nPiPlus_above65", &MC_nPiPlus_above65);
    my_event_->Branch("MC_nPiMinus_below65", &MC_nPiMinus_below65);
    my_event_->Branch("MC_nPiMinus_above65", &MC_nPiMinus_above65);
    my_event_->Branch("MC_Primary_PDG", &MC_Primary_PDG);
    my_event_->Branch("MC_Primary_Mom", &MC_Primary_Mom);
    my_event_->Branch("MC_proton_true_Mom_above255", &MC_proton_true_Mom_above255);
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

  }
  my_event_->Branch("flash_matching_chi2", &flash_matching_chi2);
  my_event_->Branch("flash_YCenter", &flash_YCenter);
  my_event_->Branch("flash_YWidth", &flash_YWidth);
  my_event_->Branch("flash_ZCenter", &flash_ZCenter);
  my_event_->Branch("flash_ZWidth", &flash_ZWidth);
  my_event_->Branch("flash_TotalPE", &flash_TotalPE);

  my_event_->Branch("n_pfp_nuDaughters", &n_pfp_nuDaughters);
  my_event_->Branch("n_dau_tracks", &n_dau_tracks);
  my_event_->Branch("n_dau_showers", &n_dau_showers);
  
  my_event_->Branch("if_2tracks", &if_2tracks);

  my_event_->Branch("crthit_PE", &crthit_PE);
  my_event_->Branch("crthit_plane", &crthit_plane);
  my_event_->Branch("crthit_time", &crthit_time);
  my_event_->Branch("Nr_crthit_inBeam", &Nr_crthit_inBeam);

  my_event_->Branch("Nr_trk_asso_crthit", &Nr_trk_asso_crthit);
  my_event_->Branch("trk_crt_time", &trk_crt_time);
  my_event_->Branch("if_trk_CRT_out_Beam", &if_trk_CRT_out_Beam);

  my_event_->Branch("trk_OutOfTime", &trk_OutOfTime);

  my_event_->Branch("vtx_InFV", &vtx_InFV);
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
  my_event_->Branch("bestMCSLL_NoSCE", &bestMCSLL_NoSCE);
  my_event_->Branch("fwdMCSLL_NoSCE", &fwdMCSLL_NoSCE);
  my_event_->Branch("bwdMCSLL_NoSCE", &bwdMCSLL_NoSCE);

  my_event_->Branch("Trk_length", &Trk_length);
  my_event_->Branch("Mom_Range_mu", &Mom_Range_mu);
  my_event_->Branch("Mom_Range_p", &Mom_Range_p);
  my_event_->Branch("Mom_Range_pi", &Mom_Range_pi);

  my_event_->Branch("avg_dEdx_LargeHit_pl0", &avg_dEdx_LargeHit_pl0);
  my_event_->Branch("avg_dEdx_LargeHit_pl1", &avg_dEdx_LargeHit_pl1);
  my_event_->Branch("avg_dEdx_LargeHit_pl2", &avg_dEdx_LargeHit_pl2);

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

  my_event_->Branch("dEdx_pl0_mid", &dEdx_pl0_mid);
  my_event_->Branch("dEdx_pl1_mid", &dEdx_pl1_mid);
  my_event_->Branch("dEdx_pl2_mid", &dEdx_pl2_mid);

  my_event_->Branch("PID_Chi2Mu_3pl", &PID_Chi2Mu_3pl);
  my_event_->Branch("PID_Chi2P_3pl", &PID_Chi2P_3pl);
  my_event_->Branch("PID_Chi2Pi_3pl", &PID_Chi2Pi_3pl);
  my_event_->Branch("PID_Chi2K_3pl", &PID_Chi2K_3pl);
  
  my_event_->Branch("trk_cosmic_percent", &trk_cosmic_percent);
  my_event_->Branch("trk_purity", &trk_purity);
  my_event_->Branch("trk_completeness", &trk_completeness);

  my_event_->Branch("if_cosmic", &if_cosmic);
  my_event_->Branch("if_matchMu", &if_matchMu);

  my_event_->Branch("muon_idx", &muon_idx);
  my_event_->Branch("proton_idx", &proton_idx);

  my_event_->Branch("vtx_x", &vtx_x);
  my_event_->Branch("vtx_y", &vtx_y);
  my_event_->Branch("vtx_z", &vtx_z);

  my_event_->Branch("muon_cand_phi", &muon_cand_phi);
  my_event_->Branch("muon_cand_theta", &muon_cand_theta);
  my_event_->Branch("muon_cand_costheta", &muon_cand_costheta);
  my_event_->Branch("muon_cand_mom_MCS", &muon_cand_mom_MCS);
  my_event_->Branch("muon_cand_mom_range", &muon_cand_mom_range);

  my_event_->Branch("proton_cand_phi", &proton_cand_phi);
  my_event_->Branch("proton_cand_theta", &proton_cand_theta);
  my_event_->Branch("proton_cand_costheta", &proton_cand_costheta);
  my_event_->Branch("proton_cand_mom_range", &proton_cand_mom_range);

  my_event_->Branch("missing_PT", &missing_PT);
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

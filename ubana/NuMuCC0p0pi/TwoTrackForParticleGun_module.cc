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

#include <tgmath.h>

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
  //detinfo::DetectorProperties const* detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  //detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();

  TTree * POTtree;
  int run, subrun;
  double POT_miss1E10;

  TTree * my_event_;
  
  void Initialize_event();

  std::vector<double> All_true_PDG;//Track pdg
  std::vector<double> All_true_mom;//True momentum of muon track in the every event
  std::vector<double> All_true_mu_mom;//True muon momentum of muon track in the every event
  std::vector<double> All_true_hadron_mom;//True p or pi momentum of muon track in the every event
  std::vector<double> All_true_start_x;//True start of muon track (X)
  std::vector<double> All_true_start_y;//True start of muon track (Y)
  std::vector<double> All_true_start_z;//True start of muon track (Z)
  std::vector<double> All_true_end_x;//True end of muon track (X)
  std::vector<double> All_true_end_y;//True end of muon track (Y)
  std::vector<double> All_true_end_z;//True end of muon track (Z)
  std::vector<bool> All_true_trk_ifcontained; // True track if contained or not
  std::vector<bool> All_true_vtxFV; // True track start in FV or not
  std::vector<double> All_true_trk_phi;
  std::vector<double> All_true_trk_theta;
  std::vector<double> All_true_trk_costheta;
  std::vector<double> All_true_mu_costheta;//True muon momentum of muon track in the every event
  std::vector<double> All_true_hadron_costheta;//True p or pi momentum of muon track in the every event
  std::vector<double> All_true_trk_theta_yz;
  std::vector<double> All_true_trk_costheta_yz;
  std::vector<double> All_true_trk_theta_xz;
  std::vector<double> All_true_trk_costheta_xz;
  std::vector<double> All_true_trk_length;

  int N_MCP; // number of MC particles in this event
  int vtx_dist; // distance in between tracks (default = -999)
  double cos_ang_muon_hadron; // angle in between muon and the hadron (default = -999)

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

  bool if_2tracks = false; // if there are two primary tracks and nothing

  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime
  double only_crthit_time_inBeam = -999;
  bool if_t0_crt_time_inBeam_match = true;
  
  bool evt_CRTveto = false; // If CRT veto, eliminate the events for contained (70PE threshold)
  bool evt_CRTveto_100 = false; // If CRT veto, eliminate the events for contained (100PE threshold)

  std::vector<int> Nr_trk_asso_crthit = {0, 0}; // Number of CRT hits associated to the track
  std::vector<double> trk_crt_time = {-999, -999};// The CRT time of the hit which matched to the track
  std::vector<bool> if_trk_CRT_out_Beam = {false, false}; // Check if a track matches with out of beam CRT hit(s)
  std::vector<bool> if_t0_trk_crt_time_match = {true, true};

  std::vector<double> crthit_time_T0corr; // Time of CRT hits
  int Nr_crthit_inBeam_T0corr = 0; // Number of CRT hits in beamtime
  double only_crthit_time_inBeam_T0corr = -999;
  bool if_t0_crt_time_inBeam_match_T0corr = true;

  std::vector<double> trk_crt_time_T0corr = {-999, -999};// The CRT time of the hit which matched to the track
  std::vector<bool> if_trk_CRT_out_Beam_T0corr = {false, false}; // Check if a track matches with out of beam CRT hit(s)
  std::vector<bool> if_t0_trk_crt_time_match_T0corr = {true, true};

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

  std::vector<double> avg_dEdx_LargeHit_pl0 = {-999, -999};
  std::vector<double> avg_dEdx_LargeHit_pl1 = {-999, -999};
  std::vector<double> avg_dEdx_LargeHit_pl2 = {-999, -999};

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

  std::vector<double> dEdx_pl0_mid = {-999, -999};
  std::vector<double> dEdx_pl1_mid = {-999, -999};
  std::vector<double> dEdx_pl2_mid = {-999, -999};

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

  bool                                IsMC;
  bool                                T0Corr;
  bool                                UsingCRT;
  double                              fBeamStart;
  double                              fBeamEnd;
  double                              fDTOffset;
  double                              fDTOffset_overlay;
  //std::string                         m_DAQHeaderProducer;
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

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


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
  //m_DAQHeaderProducer(pset.get<std::string>("DAQHeaderProducer")),
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
  m_CRTcorrT0Label(pset.get<std::string>("CRTcorrT0Label")),
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
  //// Prepare X offset for position correction (SCE)
  //double xtimeoffset = 0;
  //if(IsMC){
  //  auto const& mct_h = evt.getValidHandle<std::vector<simb::MCTruth> >("generator");
  //  auto gen = mct_h->at(0);
  //  double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
  //  xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);
  //}
  //else{
  //  xtimeoffset = 0;
  //}

  ////// Get necessary handles
  //std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  //std::vector<art::Ptr<simb::GTruth> > GTruthCollection;
  //std::vector<art::Ptr<evwgh::MCEventWeight> > WeightCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;

  if(IsMC){
   // // MC Truth
   // art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
   // evt.getByLabel(m_generatorLabel, Handle_MCTruth);
   // art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

   // // Genie Truth
   // art::Handle< std::vector<simb::GTruth> > Handle_GTruth;
   // evt.getByLabel(m_generatorLabel, Handle_GTruth);
   // art::fill_ptr_vector(GTruthCollection, Handle_GTruth);
 
    // MC Particle
    art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;
    evt.getByLabel(m_geantLabel, Handle_MCParticle);
    art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

   // // Reweight
   // art::Handle< std::vector<evwgh::MCEventWeight> > Handle_Weight;
   // if(evt.getByLabel("eventweight4to4aFix", Handle_Weight)){
   //   art::fill_ptr_vector(WeightCollection, Handle_Weight);
   //   std::map<std::string, std::vector<double>> evtwgt_map = WeightCollection.at(0)->fWeight;
   //   const std::vector<double> &weights = evtwgt_map.at("splines_general_Spline");
   //   EventWeight = weights.front();
   // }
   // else{
   //   EventWeight = 1;
   // }
  }

  // Slice
  art::Handle<std::vector<recob::Slice> > Handle_Slice;
  evt.getByLabel(m_pandoraLabel, Handle_Slice);
  //std::vector< art::Ptr<recob::Slice> > AllSliceCollection;
  //art::fill_ptr_vector(AllSliceCollection, Handle_Slice);

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

  //// DAQ
  //art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
  //evt.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);

  // CRT
  art::Handle<std::vector<crt::CRTHit>> Handle_crthit;
  evt.getByLabel(m_CRTHitLabel, Handle_crthit);
  std::vector<art::Ptr<crt::CRTHit>> crthit_v;
  art::fill_ptr_vector(crthit_v, Handle_crthit);

  // CRT corr T0
  //std::vector<art::Ptr<anab::T0>> crtT0_v;
  //if(T0Corr){
  //  art::Handle<std::vector<anab::T0>> Handle_crtcorrT0;
  //  evt.getByLabel(m_CRTcorrT0Label, Handle_crtcorrT0);
  //  art::fill_ptr_vector(crtT0_v, Handle_crtcorrT0);
  //}

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
  std::vector<TVector3> trk_startdir(2);
 
  //Constants
  //const simb::Origin_t Neutrino_Origin = simb::kBeamNeutrino;

  N_MCP = 0;
  TVector3 mu_start;
  TVector3 hadron_start;
  TVector3 mu_dir;
  TVector3 hadron_dir;
  TVector3 reco_mu_start;
  TVector3 reco_hadron_start;
  TVector3 All_true_start;
  TVector3 All_true_end;

  if(IsMC){
    //------- Get part of the generated particle info
    for(unsigned int i_mcp = 0; i_mcp < MCParticleCollection.size(); i_mcp++){
      if(MCParticleCollection[i_mcp]->Process() == "primary"){

        N_MCP++;
        auto MCP = MCParticleCollection[i_mcp];
        auto AllTrueTrackPos = MCP->EndPosition() - MCP->Position();
        All_true_PDG.push_back(MCP->PdgCode());
        All_true_mom.push_back(MCP->P());
        if(MCP->PdgCode()==13){
          All_true_mu_mom.push_back(MCP->P());
        }
        if(MCP->PdgCode()!=13){
          All_true_hadron_mom.push_back(MCP->P());
        }
        All_true_start_x.push_back(MCP->Position().X());
        All_true_start_y.push_back(MCP->Position().Y());
        All_true_start_z.push_back(MCP->Position().Z());
        All_true_end_x.push_back(MCP->EndPosition().X());
        All_true_end_y.push_back(MCP->EndPosition().Y());
        All_true_end_z.push_back(MCP->EndPosition().Z());

        All_true_start.SetXYZ(All_true_start_x.back(), All_true_start_y.back(), All_true_start_z.back());
        All_true_end.SetXYZ(All_true_end_x.back(), All_true_end_y.back(), All_true_end_z.back());
        All_true_trk_ifcontained.push_back(_fiducial_volume.TrackContain(All_true_start, All_true_end));
        //All_true_vtxFV.push_back(_fiducial_volume.VertexInFV(All_true_start));

        All_true_trk_phi.push_back(AllTrueTrackPos.Phi());
        All_true_trk_theta.push_back(AllTrueTrackPos.Theta());
        All_true_trk_costheta.push_back(cos(AllTrueTrackPos.Theta()));
        if(MCP->PdgCode()==13){
          All_true_mu_costheta.push_back(cos(AllTrueTrackPos.Theta()));
          All_true_vtxFV.push_back(_fiducial_volume.VertexInFV(All_true_start));
        }
        if(MCP->PdgCode()!=13){
          All_true_hadron_costheta.push_back(cos(AllTrueTrackPos.Theta()));
        }
        All_true_trk_theta_yz.push_back(std::atan2(AllTrueTrackPos.Y(), AllTrueTrackPos.Z()));
        All_true_trk_costheta_yz.push_back(cos(All_true_trk_theta_yz.back()));
        All_true_trk_theta_xz.push_back(std::atan2(AllTrueTrackPos.X(), AllTrueTrackPos.Z()));
        All_true_trk_costheta_xz.push_back(cos(All_true_trk_theta_xz.back()));
        All_true_trk_length.push_back(sqrt(AllTrueTrackPos.X()*AllTrueTrackPos.X() + AllTrueTrackPos.Y()*AllTrueTrackPos.Y() + AllTrueTrackPos.Z()*AllTrueTrackPos.Z())); // An estimation of true track length

        if(MCP->PdgCode() == 13){
          mu_start = All_true_start;
          mu_dir.SetXYZ(MCP->Px(), MCP->Py(), MCP->Pz());
        }
        if(MCP->PdgCode() != 13){ // In this edition, it's either a proton or pion
          hadron_start = All_true_start;
          hadron_dir.SetXYZ(MCP->Px(), MCP->Py(), MCP->Pz());
        }
      }
    }
    if(N_MCP == 2){ // assuming if there are two primary MC paritlces, they must be one muon and one othe
      vtx_dist = (mu_start - hadron_start).Mag();
      cos_ang_muon_hadron = (mu_dir * hadron_dir) / (mu_dir.Mag() * hadron_dir.Mag());
    }
    else{
      vtx_dist = -999;
      cos_ang_muon_hadron = -999;
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
        trigger_time = flash_matching_T0.front()->Time();
      }
   
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
          //if(assoTrack.empty() && assoShower.empty()){
          //  Ghost_PDG.push_back(dau_pfp->PdgCode());
          //}
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
        //double evt_timeGPS_nsec = 0.;
        //if(!rawHandle_DAQHeader.isValid()) {
        //   throw cet::exception("[Numu0pi0p]") << "The DAQ header cannot be located." << std::endl;
        // }
        //raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
        //art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
        //evt_timeGPS_nsec = evtTimeGPS.timeLow();

        //double CRTT0corr = 0;
        //if(T0Corr){
        //  if(crtT0_v.size() == 1){
        //    CRTT0corr = crtT0_v.front()->Time();
        //  }
        //  else{
        //    CRTT0corr = 0;
        //  }
        //}

        if(UsingCRT){
          // overlay is also basically data, using ts0
          for (unsigned int i_crt = 0; i_crt < crthit_v.size(); i_crt ++){
            // figure out what plane this hit comes from
            // 3 -> top, 0 -> bottom, 1 -> anode, 2 -> cathode
            //double crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
            double crt_time = crthit_v[i_crt]->ts1_ns / 1000.;

            crthit_PE.push_back(crthit_v[i_crt]->peshit);
            crthit_plane.push_back(crthit_v[i_crt]->plane);
            crthit_time.push_back(crt_time);

            if(crt_time >= fBeamStart && crt_time <= fBeamEnd && crthit_v[i_crt]->peshit > 70){
              Nr_crthit_inBeam++; // The total number of CRT hits in beam should be less than the number of exiting tracks
              only_crthit_time_inBeam = crt_time;
            } // If CRT hit in Beam window

            // with CRT T0 time correction
            if(T0Corr){
              //double crt_time_T0corr = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
              double crt_time_T0corr = crthit_v[i_crt]->ts1_ns / 1000.;
              crthit_time_T0corr.push_back(crt_time_T0corr);

              if(crt_time_T0corr >= fBeamStart && crt_time_T0corr <= fBeamEnd && crthit_v[i_crt]->peshit > 70){
                Nr_crthit_inBeam_T0corr++;
                only_crthit_time_inBeam_T0corr = crt_time_T0corr;
              } // If CRT hit in Beam window
            }

          } // CRT hit loop

          if(Nr_crthit_inBeam == 1){
            if(trigger_time != -999){
              //if(only_crthit_time_inBeam - trigger_time> -0.2){
              if(abs(only_crthit_time_inBeam - trigger_time) > 1){
                if_t0_crt_time_inBeam_match = false;
              }
            }
          }

          if(T0Corr){
            if(Nr_crthit_inBeam_T0corr == 1){
              if(trigger_time != -999){
                if(abs(only_crthit_time_inBeam_T0corr - trigger_time) > 1){
                  if_t0_crt_time_inBeam_match_T0corr = false;
                }
              }
            }
          }

          //- For contained (Veto if there is any CRT hit within 1us of the flash which is in the beam window)
          if(flash_v.size() > 0){
            for(unsigned int i_fl = 0; i_fl < flash_v.size(); i_fl++){
              auto CRT_hit = CRThitFlashAsso.at(flash_v[i_fl].key());
              if(CRT_hit.size() > 0){
                if(CRT_hit.front()->peshit > 70) evt_CRTveto = true;
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
                if(Track_CRThit[i_trk_crt]->peshit > 70){
                  // w/o T0 correction
                  //trk_crt_time[i_trk] = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
                  trk_crt_time[i_trk] = Track_CRThit[i_trk_crt]->ts1_ns / 1000.;
                  // w/ T0 correction
                  if(T0Corr){
                    //trk_crt_time_T0corr[i_trk] = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
                    trk_crt_time_T0corr[i_trk] = Track_CRThit[i_trk_crt]->ts1_ns / 1000.;
                  }

                  Nr_trk_asso_crthit[i_trk]++;
                }

                if(trk_crt_time[i_trk] < fBeamStart || trk_crt_time[i_trk] > fBeamEnd){
                  // For 2 tracks, as long as one of them has associated CRT hit out of time, reject!
                  if_trk_CRT_out_Beam[i_trk] = true;
                } // If matched CRT hit out of Beam window

                if(T0Corr){
                  if(trk_crt_time_T0corr[i_trk] < fBeamStart || trk_crt_time_T0corr[i_trk] > fBeamEnd){
                    if_trk_CRT_out_Beam_T0corr[i_trk] = true;
                  } // If matched CRT hit out of Beam window
                }

              } // Loop over the associated CRT hits
            }  // If there is associated CRT hits

            // - crt time and t0 time match
            if(trigger_time != -999 && trk_crt_time[i_trk] != -999){
              //if((trk_crt_time[i_trk] - trigger_time) == 0 || ((trk_crt_time[i_trk] - trigger_time) > -0.52 && (trk_crt_time[i_trk] - trigger_time) < -0.42)){
              if(abs(trk_crt_time[i_trk] - trigger_time) < 1){
                if_t0_trk_crt_time_match[i_trk] =  true;
              }
              else{
                if_t0_trk_crt_time_match[i_trk] = false;
              }
            }

            if(T0Corr){
              if(trigger_time != -999 && trk_crt_time_T0corr[i_trk] != -999){
                if(abs(trk_crt_time_T0corr[i_trk] - trigger_time) < 1){
                  if_t0_trk_crt_time_match_T0corr[i_trk] =  true;
                }
                else{
                  if_t0_trk_crt_time_match_T0corr[i_trk] = false;
                }
              }
            }

          } // CRT info

          //------ Reco info

          //-- Track start and end (The track start of the muon candidate track will be the vertex)
          Trk_start[i_trk] = daughter_Tracks[i_trk]->Start<TVector3>();
          auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start[i_trk].X(), Trk_start[i_trk].Y(), Trk_start[i_trk].Z()));
          //Trk_start_SCEcorr[i_trk].SetX(Trk_start[i_trk].X() - Trk_start_offset.X() + xtimeoffset + 0.6);
          Trk_start_SCEcorr[i_trk].SetX(Trk_start[i_trk].X() - Trk_start_offset.X());
          Trk_start_SCEcorr[i_trk].SetY(Trk_start[i_trk].Y() + Trk_start_offset.Y());
          Trk_start_SCEcorr[i_trk].SetZ(Trk_start[i_trk].Z() + Trk_start_offset.Z());

          trk_start_x[i_trk] = Trk_start_SCEcorr[i_trk].X();
          trk_start_y[i_trk] = Trk_start_SCEcorr[i_trk].Y();
          trk_start_z[i_trk] = Trk_start_SCEcorr[i_trk].Z();

          Trk_end[i_trk] = daughter_Tracks[i_trk]->End<TVector3>();
          auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end[i_trk].X(), Trk_end[i_trk].Y(), Trk_end[i_trk].Z()));
          //Trk_end_SCEcorr[i_trk].SetX(Trk_end[i_trk].X() - Trk_end_offset.X() + xtimeoffset + 0.6);
          Trk_end_SCEcorr[i_trk].SetX(Trk_end[i_trk].X() - Trk_end_offset.X());
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
            if(hits_dEdx_size_pl0 != 0){
              for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl0; i_hit++){
                avg_dEdx_LargeHit_pl0[i_trk] += copy_dEdx_pl0[i_hit];
              }
              avg_dEdx_LargeHit_pl0[i_trk] = avg_dEdx_LargeHit_pl0[i_trk] / hits_dEdx_size_pl0;
            }
            else{
              avg_dEdx_LargeHit_pl0[i_trk] = 0;
            }
          }
          else{
            avg_dEdx_LargeHit_pl0[i_trk] = (copy_dEdx_pl0[hits_dEdx_size_pl0 - 1] + copy_dEdx_pl0[hits_dEdx_size_pl0 - 2] + copy_dEdx_pl0[hits_dEdx_size_pl0 - 3]) / 3;
          }
          // pl 1
          if(hits_dEdx_size_pl1 < 3){
            if(hits_dEdx_size_pl1 != 0){
              for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl1; i_hit++){
                avg_dEdx_LargeHit_pl1[i_trk] += copy_dEdx_pl1[i_hit];
              }
              avg_dEdx_LargeHit_pl1[i_trk] = avg_dEdx_LargeHit_pl1[i_trk] / hits_dEdx_size_pl1;
            }
            else{
              avg_dEdx_LargeHit_pl1[i_trk] = 0;
            }
          }
          else{
            avg_dEdx_LargeHit_pl1[i_trk] = (copy_dEdx_pl1[hits_dEdx_size_pl1 - 1] + copy_dEdx_pl1[hits_dEdx_size_pl1 - 2] + copy_dEdx_pl1[hits_dEdx_size_pl1 - 3]) / 3;
          }
          // pl 2
          if(hits_dEdx_size_pl2 < 3){
            if(hits_dEdx_size_pl2 != 0){
              for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl2; i_hit++){
                avg_dEdx_LargeHit_pl2[i_trk] += copy_dEdx_pl2[i_hit];
              }
              avg_dEdx_LargeHit_pl2[i_trk] = avg_dEdx_LargeHit_pl2[i_trk] / hits_dEdx_size_pl2;
            }
            else{
              avg_dEdx_LargeHit_pl2[i_trk] = 0;
            }
          }
          else{
            avg_dEdx_LargeHit_pl2[i_trk] = (copy_dEdx_pl2[hits_dEdx_size_pl2 - 1] + copy_dEdx_pl2[hits_dEdx_size_pl2 - 2] + copy_dEdx_pl2[hits_dEdx_size_pl2 - 3]) / 3;
          }

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
              }
            }
          }
          
        } // Finish looping the 2 tracks

        ///////////////////////
        // Decide the vertex position, PT PL, which track is proton and which is muon
        //////////////////////
        //
        //-- Track particle ID (Set the threshold to be PID_Chi2P_3pl 95)
        // If we don't use track length for the moment
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
          // PT = PT(track 1) + PT(track 2) = sqrt((PT1_x + PT2_x)^2 + (PT1_y + PT2_y)^2)
          // PL = PL(track 1) + PL(track 2) = abs(PL1_z + PL2_z)
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
        
      } // If it is 2-primary-track selection
    } // If it is neutrino slice
  } // End of Loop over all pfparticles

  my_event_->Fill();

  if(IsMC){
    All_true_PDG.clear();
    All_true_mom.clear();
    All_true_mu_mom.clear();
    All_true_hadron_mom.clear();
    All_true_start_x.clear();
    All_true_start_y.clear();
    All_true_start_z.clear();
    All_true_end_x.clear();
    All_true_end_y.clear();
    All_true_end_z.clear();
    All_true_trk_ifcontained.clear();
    All_true_vtxFV.clear();
    All_true_trk_phi.clear();
    All_true_trk_theta.clear();
    All_true_trk_costheta.clear();
    All_true_mu_costheta.clear();
    All_true_hadron_costheta.clear();
    All_true_trk_theta_yz.clear();
    All_true_trk_costheta_yz.clear();
    All_true_trk_theta_xz.clear();
    All_true_trk_costheta_xz.clear();
    All_true_trk_length.clear();

    vtx_dist = -999;
    cos_ang_muon_hadron = -999;
  }

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

  evt_CRTveto = false;
  evt_CRTveto_100 = false;

  Nr_trk_asso_crthit = {0, 0}; 
  trk_crt_time = {-999, -999};
  if_trk_CRT_out_Beam = {false, false};
  if_t0_trk_crt_time_match = {true, true};

  if(T0Corr){
    crthit_time_T0corr.clear();
    Nr_crthit_inBeam_T0corr = 0;
    only_crthit_time_inBeam_T0corr = -999;
    if_t0_crt_time_inBeam_match_T0corr = true;
  
    trk_crt_time_T0corr = {-999, -999};
    if_trk_CRT_out_Beam_T0corr = {false, false};
    if_t0_trk_crt_time_match_T0corr = {true, true};
  }

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

  avg_dEdx_LargeHit_pl0 = {-999, -999};
  avg_dEdx_LargeHit_pl1 = {-999, -999};
  avg_dEdx_LargeHit_pl2 = {-999, -999};
 
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

  dEdx_pl0_mid = {-999, -999};
  dEdx_pl1_mid = {-999, -999};
  dEdx_pl2_mid = {-999, -999};

  PID_Chi2Mu_3pl = {-999, -999};
  PID_Chi2P_3pl = {-999, -999};
  PID_Chi2Pi_3pl = {-999, -999};
  PID_Chi2K_3pl = {-999, -999};

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
    my_event_->Branch("All_true_PDG", &All_true_PDG);
    my_event_->Branch("All_true_mom", &All_true_mom);
    my_event_->Branch("All_true_mu_mom", &All_true_mu_mom);
    my_event_->Branch("All_true_hadron_mom", &All_true_hadron_mom);
    my_event_->Branch("All_true_start_x", &All_true_start_x);
    my_event_->Branch("All_true_start_y", &All_true_start_y);
    my_event_->Branch("All_true_start_z", &All_true_start_z);
    my_event_->Branch("All_true_end_x", &All_true_end_x);
    my_event_->Branch("All_true_end_y", &All_true_end_y);
    my_event_->Branch("All_true_end_z", &All_true_end_z);
    my_event_->Branch("All_true_trk_ifcontained", &All_true_trk_ifcontained);
    my_event_->Branch("All_true_vtxFV", &All_true_vtxFV);
    my_event_->Branch("All_true_trk_phi", &All_true_trk_phi);
    my_event_->Branch("All_true_trk_theta", &All_true_trk_theta);
    my_event_->Branch("All_true_trk_costheta", &All_true_trk_costheta);
    my_event_->Branch("All_true_mu_costheta", &All_true_mu_costheta);
    my_event_->Branch("All_true_hadron_costheta", &All_true_hadron_costheta);
    my_event_->Branch("All_true_trk_theta_yz", &All_true_trk_theta_yz);
    my_event_->Branch("All_true_trk_costheta_yz", &All_true_trk_costheta_yz);
    my_event_->Branch("All_true_trk_theta_xz", &All_true_trk_theta_xz);
    my_event_->Branch("All_true_trk_costheta_xz", &All_true_trk_costheta_xz);
    my_event_->Branch("All_true_trk_length", &All_true_trk_length);
  
    my_event_->Branch("vtx_dist", &vtx_dist);
    my_event_->Branch("cos_ang_muon_hadron", &cos_ang_muon_hadron);
    my_event_->Branch("N_MCP", &N_MCP);

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

  my_event_->Branch("evt_CRTveto", &evt_CRTveto);
  my_event_->Branch("evt_CRTveto_100", &evt_CRTveto_100);

  my_event_->Branch("Nr_trk_asso_crthit", &Nr_trk_asso_crthit);
  my_event_->Branch("trk_crt_time", &trk_crt_time);
  my_event_->Branch("if_trk_CRT_out_Beam", &if_trk_CRT_out_Beam);
  my_event_->Branch("if_t0_trk_crt_time_match", &if_t0_trk_crt_time_match);

  if(T0Corr){
    my_event_->Branch("crthit_time_T0corr", &crthit_time_T0corr);
    my_event_->Branch("Nr_crthit_inBeam_T0corr", &Nr_crthit_inBeam_T0corr);
    my_event_->Branch("only_crthit_time_inBeam_T0corr", &only_crthit_time_inBeam_T0corr);
    my_event_->Branch("if_t0_crt_time_inBeam_match_T0corr", &if_t0_crt_time_inBeam_match_T0corr);
  
    my_event_->Branch("trk_crt_time_T0corr", &trk_crt_time_T0corr);
    my_event_->Branch("if_trk_CRT_out_Beam_T0corr", &if_trk_CRT_out_Beam_T0corr);
    my_event_->Branch("if_t0_trk_crt_time_match_T0corr", &if_t0_trk_crt_time_match_T0corr);
  }

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
  my_event_->Branch("if_matchP", &if_matchP);
  my_event_->Branch("if_matchPrimary", &if_matchPrimary);

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

  if(IsMC){
    my_event_->Branch("if_match_mu_p", &if_match_mu_p);
    my_event_->Branch("if_match_mu_p_flipped", &if_match_mu_p_flipped);
  }

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
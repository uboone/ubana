////////////////////////////////////////////////////////////////////////
// Class:       EveryTrack
// Plugin Type: analyzer (art v3_01_02)
// File:        EveryTrack_module.cc
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

class EveryTrack;


class EveryTrack : public art::EDAnalyzer {
public:
  explicit EveryTrack(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EveryTrack(EveryTrack const&) = delete;
  EveryTrack(EveryTrack&&) = delete;
  EveryTrack& operator=(EveryTrack const&) = delete;
  EveryTrack& operator=(EveryTrack&&) = delete;

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

  double EventWeight = 1; // Spine reweight using Steven's tool  

  double flash_matching_chi2 = -999; //Chi2 of flash matching in each neutrino slice
  double trigger_time = -999;
  double flash_YCenter = -999;
  double flash_YWidth = -999;
  double flash_ZCenter = -999;
  double flash_ZWidth = -999;
  double flash_TotalPE = -999;
  double flash_time = -999;

  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime
  bool if_t0_crt_time_inBeam_match = true;

  int Nr_trk_asso_crthit = 0; // Number of CRT hits associated to the track
  double trk_crt_time = -999;// The CRT time of the hit which matched to the track
  bool if_trk_CRT_out_Beam = false; // Check if a track matches with out of beam CRT hit(s)
  bool if_t0_trk_crt_time_match = true;

  std::vector<double> crthit_time_T0corr; // Time of CRT hits
  int Nr_crthit_inBeam_T0corr = 0; // Number of CRT hits in beamtime
  bool if_t0_crt_time_inBeam_match_T0corr = true;

  double trk_crt_time_T0corr = -999;// The CRT time of the hit which matched to the track
  bool if_trk_CRT_out_Beam_T0corr = false; // Check if a track matches with out of beam CRT hit(s)
  bool if_t0_trk_crt_time_match_T0corr = true;
 
  std::vector<bool> if_trk_NeutrinoSlice_primary;
 
  std::vector<double> trk_vtx_x;
  std::vector<double> trk_vtx_y;
  std::vector<double> trk_vtx_z;
  std::vector<double> trk_start_x;
  std::vector<double> trk_start_y;
  std::vector<double> trk_start_z;
  std::vector<double> trk_end_x;
  std::vector<double> trk_end_y;
  std::vector<double> trk_end_z;

  std::vector<double> trk_vtx_x_noSCE;
  std::vector<double> trk_vtx_y_noSCE;
  std::vector<double> trk_vtx_z_noSCE;
  std::vector<double> trk_start_x_noSCE;
  std::vector<double> trk_start_y_noSCE;
  std::vector<double> trk_start_z_noSCE;
  std::vector<double> trk_end_x_noSCE;
  std::vector<double> trk_end_y_noSCE;
  std::vector<double> trk_end_z_noSCE;

  std::vector<bool> trk_contained;
  std::vector<bool> vtx_InFV; 

  std::vector<bool> trk_contained_noSCE;
  std::vector<bool> vtx_InFV_noSCE;

  std::vector<double> sin2_theta_pl0;
  std::vector<double> sin2_theta_pl1;
  std::vector<double> sin2_theta_pl2;
  std::vector<double> sin2_phi_readout;

  std::vector<double> trk_phi;
  std::vector<double> trk_theta;
  std::vector<double> trk_costheta;

  std::vector<double> trk_length;

  std::vector<double> mom_Range_mu;
  std::vector<double> mom_Range_p; 
  std::vector<double> mom_Range_pi;

  std::vector<double> bestMCS;
  std::vector<double> bestMCSLL;
  std::vector<double> fwdMCS;
  std::vector<double> fwdMCSLL;
  std::vector<double> bwdMCS;
  std::vector<double> bwdMCSLL;

  std::vector<double> bestMCSLL_NoSCE;
  std::vector<double> fwdMCSLL_NoSCE;
  std::vector<double> bwdMCSLL_NoSCE;

  std::vector<double> PID_Chi2Mu_3pl;
  std::vector<double> PID_Chi2P_3pl; 
  std::vector<double> PID_Chi2Pi_3pl;
  std::vector<double> PID_Chi2K_3pl;

  std::vector<double> dEdx_pl0_start_half;
  std::vector<double> dEdx_pl1_start_half;
  std::vector<double> dEdx_pl2_start_half;
  std::vector<double> dEdx_pl0_end_half;
  std::vector<double> dEdx_pl1_end_half;
  std::vector<double> dEdx_pl2_end_half;

  std::vector<double> dEdx_pl0_start1020;
  std::vector<double> dEdx_pl1_start1020;
  std::vector<double> dEdx_pl2_start1020;
  std::vector<double> dEdx_pl0_end1020;
  std::vector<double> dEdx_pl1_end1020;
  std::vector<double> dEdx_pl2_end1020;

  std::vector<double> dEdx_pl0_1020_ratio;
  std::vector<double> dEdx_pl0_half_ratio;
  std::vector<double> dEdx_pl1_1020_ratio;
  std::vector<double> dEdx_pl1_half_ratio;
  std::vector<double> dEdx_pl2_1020_ratio;
  std::vector<double> dEdx_pl2_half_ratio;

  std::vector<double> dEdx_pl0_mid;
  std::vector<double> dEdx_pl1_mid;
  std::vector<double> dEdx_pl2_mid;

  std::vector<double> avg_dEdx_LargeHit_pl0;
  std::vector<double> avg_dEdx_LargeHit_pl1;
  std::vector<double> avg_dEdx_LargeHit_pl2;

  std::vector<double> dEdx_median_pl0;
  std::vector<double> dEdx_median_pl1;
  std::vector<double> dEdx_median_pl2;

  std::vector<double> dEdx_pl0_start_half_clean;
  std::vector<double> dEdx_pl1_start_half_clean;
  std::vector<double> dEdx_pl2_start_half_clean;
  std::vector<double> dEdx_pl0_end_half_clean;
  std::vector<double> dEdx_pl1_end_half_clean;
  std::vector<double> dEdx_pl2_end_half_clean;

  std::vector<double> dEdx_pl0_start1020_clean;
  std::vector<double> dEdx_pl1_start1020_clean;
  std::vector<double> dEdx_pl2_start1020_clean;
  std::vector<double> dEdx_pl0_end1020_clean;
  std::vector<double> dEdx_pl1_end1020_clean;
  std::vector<double> dEdx_pl2_end1020_clean;

  std::vector<double> dEdx_pl0_1020_ratio_clean;
  std::vector<double> dEdx_pl0_half_ratio_clean;
  std::vector<double> dEdx_pl1_1020_ratio_clean;
  std::vector<double> dEdx_pl1_half_ratio_clean;
  std::vector<double> dEdx_pl2_1020_ratio_clean;
  std::vector<double> dEdx_pl2_half_ratio_clean;

  std::vector<double> dEdx_pl0_mid_clean;
  std::vector<double> dEdx_pl1_mid_clean;
  std::vector<double> dEdx_pl2_mid_clean;

  std::vector<double> avg_dEdx_LargeHit_pl0_clean;
  std::vector<double> avg_dEdx_LargeHit_pl1_clean;
  std::vector<double> avg_dEdx_LargeHit_pl2_clean;

  std::vector<double> dEdx_median_pl0_clean;
  std::vector<double> dEdx_median_pl1_clean;
  std::vector<double> dEdx_median_pl2_clean;

  std::vector<double> trk_cosmic_percent;
  std::vector<double> trk_purity;
  std::vector<double> trk_completeness;
  std::vector<bool> if_cosmic;
  std::vector<bool> if_matchPrimary;
  std::vector<bool> if_matchMu;

  std::vector<double> true_mom;//True momentum of muon track in the every event
  std::vector<double> true_start_x;//True start of muon track (X)
  std::vector<double> true_start_y;//True start of muon track (Y)
  std::vector<double> true_start_z;//True start of muon track (Z)
  std::vector<double> true_end_x;//True end of muon track (X)
  std::vector<double> true_end_y;//True end of muon track (Y)
  std::vector<double> true_end_z;//True end of muon track (Z)
  std::vector<double> true_trk_phi;//True phi of muon track
  std::vector<double> true_trk_theta;//True theta of muon track
  std::vector<double> true_trk_costheta;//True cos(theta) of muon track
  std::vector<double> true_trk_theta_yz;
  std::vector<double> true_trk_costheta_yz;
  std::vector<double> true_trk_theta_xz;
  std::vector<double> true_trk_costheta_xz;
  std::vector<double> true_trk_length;//True track length (distance from the start to the end point)
  std::vector<double> true_trk_PDG;//Track pdg
  std::vector<bool> true_trk_ifcontained; // True track if contained or not
  std::vector<bool> true_vtxFV; // True track if contained or not

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

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


EveryTrack::EveryTrack(fhicl::ParameterSet const& pset)
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

void EveryTrack::analyze(art::Event const& evt)
{

  //// Get necessary handles
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
  std::vector<art::Ptr<evwgh::MCEventWeight> > WeightCollection;

  if(IsMC){
    // MC Truth
    art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
    evt.getByLabel(m_generatorLabel, Handle_MCTruth);
    art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);
 
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

  // CRT corr T0
  art::Handle<std::vector<anab::T0>> Handle_crtcorrT0;
  evt.getByLabel(m_CRTcorrT0Label, Handle_crtcorrT0);
  std::vector<art::Ptr<anab::T0>> crtT0_v;
  art::fill_ptr_vector(crtT0_v, Handle_crtcorrT0);

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

  std::vector<art::Ptr<recob::Track> > daughter_Tracks;

  TVector3 Trk_vtx;
  TVector3 Trk_start;
  TVector3 Trk_end;
  TVector3 Trk_vtx_SCEcorr;
  TVector3 Trk_start_SCEcorr; 
  TVector3 Trk_end_SCEcorr;

  //-- OpFlash
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

  //-- CRT
  for (unsigned int i_crt = 0; i_crt < crthit_v.size(); i_crt ++){
    double crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);

    crthit_PE.push_back(crthit_v[i_crt]->peshit);
    crthit_plane.push_back(crthit_v[i_crt]->plane);
    crthit_time.push_back(crt_time);

    // only count if the CRT hit has signal above 70 PE
    if(crt_time >= fBeamStart && crt_time <= fBeamEnd && crthit_v[i_crt]->peshit > 70){
      Nr_crthit_inBeam++;
    } // If CRT hit in Beam window

    // with CRT T0 time correction
    if(T0Corr){
      double crt_time_T0corr = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
      crthit_time_T0corr.push_back(crt_time_T0corr);

      if(crt_time_T0corr >= fBeamStart && crt_time_T0corr <= fBeamEnd && crthit_v[i_crt]->peshit > 70){
        Nr_crthit_inBeam_T0corr++;
      } // If CRT hit in Beam window

    }
  }

  nr_pfp = pfParticle_v.size();

  for(unsigned int i_pfp = 0; i_pfp < pfParticle_v.size(); i_pfp++){
    auto pfp = pfParticle_v[i_pfp];
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      nr_pfp_nuDaughters = pfp->NumDaughters();
      for(int j = 0; j < n_pfp_nuDaughters; j++){
        auto Iterator = pfParticleIdMap.find(pfp->Daughters().at(j));
        auto dau_pfp = Iterator->second;

        auto assoTrack = pfpToTrackAsso.at(dau_pfp.key());
        if(assoTrack.size()==1){
          daughter_Tracks.push_back(assoTrack.front());
        }
        if(assoTrack.size()>1){
          throw cet::exception("[Numu0pi0p]") << "PFParticle has >1 track!" << std::endl;
        }

        auto assoShower = pfpToShowerAsso.at(dau_pfp.key()); // vector
        if(assoShower.size()==1){
          daughter_Showers.push_back(assoShower.front());
        }
        if(assoShower.size()>1){
          throw cet::exception("[Numu0pi0p]") << "PFParticle has >1 shower!" << std::endl;
        }

      }
    }
  }
  nr_dau_tracks = daughter_Tracks.size();
  nr_dau_showers = daughter_Showers.size();

  //-------- Get All tracks
  for(unsigned int i_trk = 0; i_trk < AllTrackCollection.size(); i_trk++){

    for(unsigned int i_pri_trk = 0; i_pri_trk < daughter_Tracks.size(); i_pri_trk++){
      if(AllTrackCollection[i_trk] == daughter_Tracks[i_pri_trk]){
        if_trk_NeutrinoSlice_primary = true;
        break;
      }
      else{
        if_trk_NeutrinoSlice_primary = false;
      }
    }

    // Track CRT
    auto Track_CRThit = CRTToTrackAsso.at(AllTrackCollection[i_trk].key());
    if(Track_CRThit.size() > 0){
      for(unsigned int i_trk_crt = 0; i_trk_crt < Track_CRThit.size(); i_trk_crt++){
        if(Track_CRThit[i_trk_crt]->peshit > 70){
          // w/o T0 correction
          trk_crt_time = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
          // w/ T0 correction
          if(T0Corr){
            trk_crt_time_T0corr = ((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
            Nr_trk_asso_crthit++;
          }
        }
      }
    }

    // Add spatial correction to the track start and end
    Trk_vtx = AllTrackCollection[i_trk]->Vertex<TVector3>();
    auto Trk_vtx_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_vtx.X(), Trk_vtx.Y(), Trk_vtx.Z()));
    Trk_vtx_SCEcorr.SetX(Trk_vtx.X() - Trk_vtx_offset.X());
    Trk_vtx_SCEcorr.SetY(Trk_vtx.Y() + Trk_vtx_offset.Y());
    Trk_vtx_SCEcorr.SetZ(Trk_vtx.Z() + Trk_vtx_offset.Z());

    trk_vtx_x = Trk_vtx_SCEcorr.X();
    trk_vtx_y = Trk_vtx_SCEcorr.Y();
    trk_vtx_z = Trk_vtx_SCEcorr.Z();

    trk_vtx_x_noSCE = Trk_vtx.X();
    trk_vtx_y_noSCE = Trk_vtx.Y();
    trk_vtx_z_noSCE = Trk_vtx.Z();

    Trk_start = AllTrackCollection[i_trk]->Start<TVector3>();
    auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
    Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
    Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
    Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

    trk_start_x = Trk_start_SCEcorr.X();
    trk_start_y = Trk_start_SCEcorr.Y();
    trk_start_z = Trk_start_SCEcorr.Z();

    trk_start_x_noSCE = Trk_start.X();
    trk_start_y_noSCE = Trk_start.Y();
    trk_start_z_noSCE = Trk_start.Z();

    Trk_end = AllTrackCollection[i_trk]->End<TVector3>();
    auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
    Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X());
    Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
    Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

    trk_end_x = Trk_end_SCEcorr.X();
    trk_end_y = Trk_end_SCEcorr.Y();
    trk_end_z = Trk_end_SCEcorr.Z();

    trk_end_x_noSCE = Trk_end.X();
    trk_end_y_noSCE = Trk_end.Y();
    trk_end_z_noSCE = Trk_end.Z();

    //-- track containment
    trk_contained = _fiducial_volume.TrackContain(Trk_start_SCEcorr, Trk_end_SCEcorr);
    vtx_InFV = _fiducial_volume.VertexInFV(Trk_start_SCEcorr);

    trk_contained_noSCE = _fiducial_volume.TrackContain(Trk_start, Trk_end);
    vtx_InFV_noSCE = _fiducial_volume.VertexInFV(Trk_start);
        
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
    trk_phi = AllTrackCollection[i_trk]->Phi();
    trk_theta = AllTrackCollection[i_trk]->Theta();
    trk_costheta = cos(AllTrackCollection[i_trk]->Theta());

    //-- track length
    trk_length = AllTrackCollection[i_trk]->Length();

    //-- Momentum
    //- range
    mom_Range_mu = _trk_mom_calculator.GetTrackMomentum(trk_length, 13);
    mom_Range_p = _trk_mom_calculator.GetTrackMomentum(trk_length, 2212);
    mom_Range_pi = _trk_mom_calculator.GetTrackMomentum(trk_length, 211);

    //- MCS
    bestMCS =  mcsfitresult_mu_v.at(AllTrackCollection[i_trk].key())->bestMomentum();
    bestMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[i_trk].key())->bestLogLikelihood();
    fwdMCS =  mcsfitresult_mu_v.at(AllTrackCollection[i_trk].key())->fwdMomentum();
    fwdMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[i_trk].key())->fwdLogLikelihood();
    bwdMCS =  mcsfitresult_mu_v.at(AllTrackCollection[i_trk].key())->bwdMomentum();
    bwdMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[i_trk].key())->bwdLogLikelihood();

    bestMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(AllTrackCollection[i_trk].key())->bestLogLikelihood();
    fwdMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(AllTrackCollection[i_trk].key())->fwdLogLikelihood();
    bwdMCSLL_NoSCE = mcsfitresult_mu_NoSCE_v.at(AllTrackCollection[i_trk].key())->bwdLogLikelihood();

    //-- dE/dx
    auto assoCal = trackToCalAsso.at(AllTrackCollection[i_trk].key());
    if(assoCal.size()!=3){
      throw cet::exception("[Numu0pi0p]") << "Where are the three planes for the calorimetry!" << std::endl;
    }
   
    std::vector<float> dEdx_pl0 = assoCal[0]->dEdx();
    std::vector<float> dQdx_pl0 = assoCal[0]->dQdx();

    std::vector<float> dEdx_pl1 = assoCal[1]->dEdx();
    std::vector<float> dQdx_pl1 = assoCal[1]->dQdx();

    std::vector<float> dEdx_pl2 = assoCal[2]->dEdx();
    std::vector<float> dQdx_pl2 = assoCal[2]->dQdx();

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

    dEdx_pl0_1020_ratio = dEdx_pl0_end1020 / (dEdx_pl0_end1020 + dEdx_pl0_start1020);
    dEdx_pl0_half_ratio = dEdx_pl0_end_half / (dEdx_pl0_end_half + dEdx_pl0_start_half);

    dEdx_pl1_1020_ratio = dEdx_pl1_end1020 / (dEdx_pl1_end1020 + dEdx_pl1_start1020);
    dEdx_pl1_half_ratio = dEdx_pl1_end_half / (dEdx_pl1_end_half + dEdx_pl1_start_half);

    dEdx_pl2_1020_ratio = dEdx_pl2_end1020 / (dEdx_pl2_end1020 + dEdx_pl2_start1020);
    dEdx_pl2_half_ratio = dEdx_pl2_end_half / (dEdx_pl2_end_half + dEdx_pl2_start_half);

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

    // - median dE/dx
    // pl 0
    if(fmod(copy_dEdx_pl0.size() / 2., 1) == 0){
      dEdx_median_pl0 = copy_dEdx_pl0[copy_dEdx_pl0.size() / 2];
    }
    else if (fmod(copy_dEdx_pl0.size() / 2., 1) == 0.5){
      dEdx_median_pl0 = (copy_dEdx_pl0[copy_dEdx_pl0.size() / 2] + copy_dEdx_pl0[copy_dEdx_pl0.size() / 2 + 1]) / 2.;
    }
    // pl 1
    if(fmod(copy_dEdx_pl1.size() / 2., 1) == 0){
      dEdx_median_pl1 = copy_dEdx_pl1[copy_dEdx_pl1.size() / 2];
    }
    else if (fmod(copy_dEdx_pl1.size() / 2., 1) == 0.5){
      dEdx_median_pl1 = (copy_dEdx_pl1[copy_dEdx_pl1.size() / 2] + copy_dEdx_pl1[copy_dEdx_pl1.size() / 2 + 1]) / 2.;
    }
    // pl 2
    if(fmod(copy_dEdx_pl2.size() / 2., 1) == 0){
      dEdx_median_pl2 = copy_dEdx_pl2[copy_dEdx_pl2.size() / 2];
    }
    else if (fmod(copy_dEdx_pl2.size() / 2., 1) == 0.5){
      dEdx_median_pl2 = (copy_dEdx_pl2[copy_dEdx_pl2.size() / 2] + copy_dEdx_pl2[copy_dEdx_pl2.size() / 2 + 1]) / 2.;
    }

    //-- Clean dE/dx, remove dE/dx is below 0.5 and above 70 MeV (artificial)
    std::vector<float> dEdx_pl0_clean;
    std::vector<float> dEdx_pl1_clean;
    std::vector<float> dEdx_pl2_clean;

    for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl0; i_hit++){
      if(dEdx_pl0[i_hit] > 0.5 && dEdx_pl0[i_hit] < 70){
        dEdx_pl0_clean.push_back(dEdx_pl0[i_hit]);
      }
    }

    for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl1; i_hit++){
      if(dEdx_pl1[i_hit] > 0.5 && dEdx_pl1[i_hit] < 70){
        dEdx_pl1_clean.push_back(dEdx_pl1[i_hit]);
      }
    }

    for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl2; i_hit++){
      if(dEdx_pl2[i_hit] > 0.5 && dEdx_pl2[i_hit] < 70){
        dEdx_pl2_clean.push_back(dEdx_pl2[i_hit]);
      }
    }

    auto hits_dEdx_size_pl0_clean = dEdx_pl0_clean.size();
    auto hits_dEdx_size_pl1_clean = dEdx_pl1_clean.size();
    auto hits_dEdx_size_pl2_clean = dEdx_pl2_clean.size();

    auto half_size_pl0_clean = hits_dEdx_size_pl0_clean / 2;
    auto half_size_pl1_clean = hits_dEdx_size_pl1_clean / 2;
    auto half_size_pl2_clean = hits_dEdx_size_pl2_clean / 2;

    //- ratio
    // dEdx half
    dEdx_pl0_start_half_clean = std::accumulate(dEdx_pl0_clean.end() - half_size_pl0_clean, dEdx_pl0_clean.end(), 0.) / half_size_pl0_clean;
    dEdx_pl0_end_half_clean = std::accumulate(dEdx_pl0_clean.begin(), dEdx_pl0_clean.begin() + half_size_pl0_clean, 0. ) / half_size_pl0_clean;
    dEdx_pl1_start_half_clean = std::accumulate(dEdx_pl1_clean.end() - half_size_pl1_clean, dEdx_pl1_clean.end(), 0.) / half_size_pl1_clean;
    dEdx_pl1_end_half_clean = std::accumulate(dEdx_pl1_clean.begin(), dEdx_pl1_clean.begin() + half_size_pl1_clean, 0. ) / half_size_pl1_clean;
    dEdx_pl2_start_half_clean = std::accumulate(dEdx_pl2_clean.end() - half_size_pl2_clean, dEdx_pl2_clean.end(), 0.) / half_size_pl2_clean;
    dEdx_pl2_end_half_clean = std::accumulate(dEdx_pl2_clean.begin(), dEdx_pl2_clean.begin() + half_size_pl2_clean, 0. ) / half_size_pl2_clean;

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

    dEdx_pl0_1020_ratio_clean = dEdx_pl0_end1020_clean / (dEdx_pl0_end1020_clean + dEdx_pl0_start1020_clean);
    dEdx_pl0_half_ratio_clean = dEdx_pl0_end_half_clean / (dEdx_pl0_end_half_clean + dEdx_pl0_start_half_clean);

    dEdx_pl1_1020_ratio_clean = dEdx_pl1_end1020_clean / (dEdx_pl1_end1020_clean + dEdx_pl1_start1020_clean);
    dEdx_pl1_half_ratio_clean = dEdx_pl1_end_half_clean / (dEdx_pl1_end_half_clean + dEdx_pl1_start_half_clean);

    dEdx_pl2_1020_ratio_clean = dEdx_pl2_end1020_clean / (dEdx_pl2_end1020_clean + dEdx_pl2_start1020_clean);
    dEdx_pl2_half_ratio_clean = dEdx_pl2_end_half_clean / (dEdx_pl2_end_half_clean + dEdx_pl2_start_half_clean);

    //- dEdx of the middle part of the track
    int nr_pl0_mid_clean = dEdx_pl0_clean.size()/3;
    int nr_pl1_mid_clean = dEdx_pl1_clean.size()/3;
    int nr_pl2_mid_clean = dEdx_pl2_clean.size()/3;
    if(nr_pl0_mid_clean > 0){ //pl0
      dEdx_pl0_mid_clean = std::accumulate(dEdx_pl0_clean.begin() + nr_pl0_mid_clean , dEdx_pl0_clean.begin() + 2*nr_pl0_mid_clean, 0.) / nr_pl0_mid_clean;
    }
    else {
      dEdx_pl0_mid_clean = 0;
    }
    if(nr_pl1_mid_clean > 0){ //pl1
      dEdx_pl1_mid_clean = std::accumulate(dEdx_pl1_clean.begin() + nr_pl1_mid_clean , dEdx_pl1_clean.begin() + 2*nr_pl1_mid_clean, 0.) / nr_pl1_mid_clean;
    }
    else {
      dEdx_pl1_mid_clean = 0;
    }
    if(nr_pl2_mid_clean > 0){ //pl2
      dEdx_pl2_mid_clean = std::accumulate(dEdx_pl2_clean.begin() + nr_pl2_mid_clean , dEdx_pl2_clean.begin() + 2*nr_pl2_mid_clean, 0.) / nr_pl2_mid_clean;
    }
    else {
      dEdx_pl2_mid_clean = 0;
    }

    // - Large hits (the three largest dE/dx)
     auto copy_dEdx_pl0_clean = dEdx_pl0_clean;
     auto copy_dEdx_pl1_clean = dEdx_pl1_clean;
     auto copy_dEdx_pl2_clean = dEdx_pl2_clean;

     std::sort(copy_dEdx_pl0_clean.begin(), copy_dEdx_pl0_clean.end());
     std::sort(copy_dEdx_pl1_clean.begin(), copy_dEdx_pl1_clean.end());
     std::sort(copy_dEdx_pl2_clean.begin(), copy_dEdx_pl2_clean.end());
     // pl 0
     if(hits_dEdx_size_pl0_clean < 3){
       for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl0_clean; i_hit++){
         avg_dEdx_LargeHit_pl0_clean += copy_dEdx_pl0_clean[i_hit];
       }
       avg_dEdx_LargeHit_pl0_clean = avg_dEdx_LargeHit_pl0_clean / hits_dEdx_size_pl0_clean;
     }
     else{
       avg_dEdx_LargeHit_pl0_clean = (copy_dEdx_pl0_clean[hits_dEdx_size_pl0_clean - 1] + copy_dEdx_pl0_clean[hits_dEdx_size_pl0_clean - 2] + copy_dEdx_pl0_clean[hits_dEdx_size_pl0_clean - 3]) / 3;
     }
     // pl 1
     if(hits_dEdx_size_pl1_clean < 3){
       for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl1_clean; i_hit++){
         avg_dEdx_LargeHit_pl1_clean += copy_dEdx_pl1_clean[i_hit];
       }
       avg_dEdx_LargeHit_pl1_clean = avg_dEdx_LargeHit_pl1_clean / hits_dEdx_size_pl1_clean;
     }
     else{
       avg_dEdx_LargeHit_pl1_clean = (copy_dEdx_pl1_clean[hits_dEdx_size_pl1_clean - 1] + copy_dEdx_pl1_clean[hits_dEdx_size_pl1_clean - 2] + copy_dEdx_pl1_clean[hits_dEdx_size_pl1_clean - 3]) / 3;
     }
     // pl 2
     if(hits_dEdx_size_pl2_clean < 3){
       for(unsigned int i_hit = 0; i_hit < hits_dEdx_size_pl2_clean; i_hit++){
         avg_dEdx_LargeHit_pl2_clean += copy_dEdx_pl2_clean[i_hit];
       }
       avg_dEdx_LargeHit_pl2_clean = avg_dEdx_LargeHit_pl2_clean / hits_dEdx_size_pl2_clean;
     }
     else{
       avg_dEdx_LargeHit_pl2_clean = (copy_dEdx_pl2_clean[hits_dEdx_size_pl2_clean - 1] + copy_dEdx_pl2_clean[hits_dEdx_size_pl2_clean - 2] + copy_dEdx_pl2_clean[hits_dEdx_size_pl2_clean - 3]) / 3;
     }

     // - median dE/dx
     // pl 0
     if(fmod(copy_dEdx_pl0_clean.size() / 2., 1) == 0){
      dEdx_median_pl0_clean = copy_dEdx_pl0_clean[copy_dEdx_pl0_clean.size() / 2];
    }
    else if (fmod(copy_dEdx_pl0_clean.size() / 2., 1) == 0.5){
      dEdx_median_pl0_clean = (copy_dEdx_pl0_clean[copy_dEdx_pl0_clean.size() / 2] + copy_dEdx_pl0_clean[copy_dEdx_pl0_clean.size() / 2 + 1]) / 2.;
    }
    // pl 1
    if(fmod(copy_dEdx_pl1_clean.size() / 2., 1) == 0){
      dEdx_median_pl1_clean = copy_dEdx_pl1_clean[copy_dEdx_pl1_clean.size() / 2];
    }
    else if (fmod(copy_dEdx_pl1_clean.size() / 2., 1) == 0.5){
      dEdx_median_pl1_clean = (copy_dEdx_pl1_clean[copy_dEdx_pl1_clean.size() / 2] + copy_dEdx_pl1_clean[copy_dEdx_pl1_clean.size() / 2 + 1]) / 2.;
    }
    // pl 2
    if(fmod(copy_dEdx_pl2_clean.size() / 2., 1) == 0){
      dEdx_median_pl2_clean = copy_dEdx_pl2_clean[copy_dEdx_pl2_clean.size() / 2];
    }
    else if (fmod(copy_dEdx_pl2_clean.size() / 2., 1) == 0.5){
      dEdx_median_pl2_clean = (copy_dEdx_pl2_clean[copy_dEdx_pl2_clean.size() / 2] + copy_dEdx_pl2_clean[copy_dEdx_pl2_clean.size() / 2 + 1]) / 2.;
    }

    //-- PID
     PID pid;
     pid.Chi2(PIDTotrackAsso,AllTrackCollection[i_trk], Trk_start_SCEcorr, Trk_end_SCEcorr,hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2);
     PID_Chi2Mu_3pl = pid.PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
     PID_Chi2P_3pl = pid.PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID
     PID_Chi2Pi_3pl = pid.PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID
     PID_Chi2K_3pl = pid.PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

     /////////////////
     // Fill TRUE info
     ////////////////
     if(IsMC){
       std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[i_trk].key());
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
     }

     std::cout<<"Tree Fill"<<std::endl;
     my_event_->Fill();

     if_trk_NeutrinoSlice_primary = false;

     if(IsMC){
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

       trk_cosmic_percent = -999;
       trk_purity = -999;
       trk_completeness = -999;
       if_cosmic = false;
       if_matchPrimary = false;
       if_matchMu = false;
     }
       
    
  } // End of track loop 

}

void EveryTrack::Initialize_event()
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

    my_event_->Branch("if_trk_NeutrinoSlice_primary", &if_trk_NeutrinoSlice_primary);

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
    
    my_event_->Branch("if_cosmic", &if_cosmic);
    my_event_->Branch("if_matchPrimary", &if_matchPrimary);
    my_event_->Branch("if_matchMu", &if_matchMu);
    my_event_->Branch("trk_cosmic_percent", &trk_cosmic_percent);
    my_event_->Branch("trk_purity", &trk_purity);
    my_event_->Branch("trk_completeness", &trk_completeness);
  }

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

  my_event_->Branch("trk_contained", &trk_contained);
  my_event_->Branch("vtx_InFV", &vtx_InFV);
  my_event_->Branch("trk_contained_noSCE", &trk_contained_noSCE);
  my_event_->Branch("vtx_InFV_noSCE", &vtx_InFV_noSCE);

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

  my_event_->Branch("PID_Chi2Mu_3pl", &PID_Chi2Mu_3pl);
  my_event_->Branch("PID_Chi2P_3pl", &PID_Chi2P_3pl);
  my_event_->Branch("PID_Chi2Pi_3pl", &PID_Chi2Pi_3pl);
  my_event_->Branch("PID_Chi2K_3pl", &PID_Chi2K_3pl);

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

  my_event_->Branch("dEdx_pl0_mid", &dEdx_pl0_mid);
  my_event_->Branch("dEdx_pl1_mid", &dEdx_pl1_mid);
  my_event_->Branch("dEdx_pl2_mid", &dEdx_pl2_mid);

  my_event_->Branch("avg_dEdx_LargeHit_pl0", &avg_dEdx_LargeHit_pl0);
  my_event_->Branch("avg_dEdx_LargeHit_pl1", &avg_dEdx_LargeHit_pl1);
  my_event_->Branch("avg_dEdx_LargeHit_pl2", &avg_dEdx_LargeHit_pl2);

  my_event_->Branch("dEdx_median_pl0", &dEdx_median_pl0);
  my_event_->Branch("dEdx_median_pl1", &dEdx_median_pl1);
  my_event_->Branch("dEdx_median_pl2", &dEdx_median_pl2);

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

  my_event_->Branch("dEdx_pl0_1020_ratio_clean", &dEdx_pl0_1020_ratio_clean);
  my_event_->Branch("dEdx_pl0_half_ratio_clean", &dEdx_pl0_half_ratio_clean);
  my_event_->Branch("dEdx_pl1_1020_ratio_clean", &dEdx_pl1_1020_ratio_clean);
  my_event_->Branch("dEdx_pl1_half_ratio_clean", &dEdx_pl1_half_ratio_clean);
  my_event_->Branch("dEdx_pl2_1020_ratio_clean", &dEdx_pl2_1020_ratio_clean);
  my_event_->Branch("dEdx_pl2_half_ratio_clean", &dEdx_pl2_half_ratio_clean);

  my_event_->Branch("dEdx_pl0_mid_clean", &dEdx_pl0_mid_clean);
  my_event_->Branch("dEdx_pl1_mid_clean", &dEdx_pl1_mid_clean);
  my_event_->Branch("dEdx_pl2_mid_clean", &dEdx_pl2_mid_clean);

  my_event_->Branch("avg_dEdx_LargeHit_pl0_clean", &avg_dEdx_LargeHit_pl0_clean);
  my_event_->Branch("avg_dEdx_LargeHit_pl1_clean", &avg_dEdx_LargeHit_pl1_clean);
  my_event_->Branch("avg_dEdx_LargeHit_pl2_clean", &avg_dEdx_LargeHit_pl2_clean);

  my_event_->Branch("dEdx_median_pl0_clean", &dEdx_median_pl0_clean);
  my_event_->Branch("dEdx_median_pl1_clean", &dEdx_median_pl1_clean);
  my_event_->Branch("dEdx_median_pl2_clean", &dEdx_median_pl2_clean);


}

void EveryTrack::endSubRun(art::SubRun const &sr){

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


void EveryTrack::beginJob()
{
  // Implementation of optional member function here.
  Initialize_event();
}

void EveryTrack::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(EveryTrack)

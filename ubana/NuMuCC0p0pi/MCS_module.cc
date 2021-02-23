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
  // Get time offset for x space charge correction
  detinfo::DetectorProperties const* detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();

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

  bool evt_CRTveto = false; // If CRT veto, eliminate the events for contained (70PE threshold)
  bool evt_CRTveto_100 = false; // If CRT veto, eliminate the events for contained (100PE threshold)

  double CRTT0corr = 0;
  std::vector<double> crthit_PE; // The photonelectrons of CRT hits which are in beam window
  std::vector<double> crthit_plane; // Plane of CRT hits
  std::vector<double> crthit_time; // Time of CRT hits
  int Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime

  std::vector<int> Nr_trk_asso_crthit; // Number of CRT hits associated to the track
  std::vector<double> trk_crt_time;// The CRT time of the hit which matched to the track
  std::vector<double> trk_crt_time_T0asso;// The CRT time of the hit which matched to the track

 
  int nr_pfp;
  int nr_pfp_nuDaughters;
  int nr_dau_tracks;
  int nr_dau_showers;

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

  std::vector<double> bestMCS_bulk;
  std::vector<double> bestMCSLL_bulk;
  std::vector<double> fwdMCS_bulk;
  std::vector<double> fwdMCSLL_bulk;
  std::vector<double> bwdMCS_bulk;
  std::vector<double> bwdMCSLL_bulk;

  std::vector<double> bestMCS_NoSCE;
  std::vector<double> bestMCSLL_NoSCE;
  std::vector<double> fwdMCS_NoSCE;
  std::vector<double> fwdMCSLL_NoSCE;
  std::vector<double> bwdMCS_NoSCE;
  std::vector<double> bwdMCSLL_NoSCE;

  std::vector<double> PID_Chi2Mu_3pl;
  std::vector<double> PID_Chi2P_3pl; 
  std::vector<double> PID_Chi2Pi_3pl;
  std::vector<double> PID_Chi2K_3pl;

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

  //int event = 0;


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
  std::string                         m_MCSmuBulkProducerLabel;
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
  m_MCSmuBulkProducerLabel(pset.get<std::string>("MCSmuBulkProducerLabel")),
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

  // Initialise X offset for postion correction (SCE)
  double xtimeoffset = 0;

  if(IsMC){
    auto const&  mct_h = evt.getValidHandle<std::vector<simb::MCTruth> >("generator");
    auto gen = mct_h->at(0);
    double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
    xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);
  }
  else{
    xtimeoffset = 0;
  }

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

  //MCS momentum (bulk position correction)
  art::Handle<std::vector<recob::MCSFitResult> > Handle_mcsfitresult_mu_bulk;
  evt.getByLabel(m_MCSmuBulkProducerLabel, Handle_mcsfitresult_mu_bulk);
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_bulk_v;
  art::fill_ptr_vector(mcsfitresult_mu_bulk_v, Handle_mcsfitresult_mu_bulk);

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
  std::vector<art::Ptr<anab::T0>> crtT0_v;
  if(T0Corr){
    art::Handle<std::vector<anab::T0>> Handle_crtcorrT0;
    evt.getByLabel(m_CRTcorrT0Label, Handle_crtcorrT0);
    art::fill_ptr_vector(crtT0_v, Handle_crtcorrT0);
  }

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

  // CRT T0 Track association
  art::FindManyP<anab::T0> CRTT0ToTrackAsso(Handle_TPCtrack,evt,CRT_TrackAssLabel);

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
  std::vector<art::Ptr<recob::Shower> > daughter_Showers;

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

      //CRT veto
      auto CRT_hit = CRThitFlashAsso.at(flash_v[i_flash].key());
      if(CRT_hit.size() == 1){
        if(CRT_hit.front()->peshit > 70) evt_CRTveto = true;
        if(CRT_hit.front()->peshit > 100) evt_CRTveto_100 = true;
      } // if CRT veto
    }
  }


  //-- CRT
  double evt_timeGPS_nsec = 0.;
  if(!rawHandle_DAQHeader.isValid()) {
     std::cout << "Could not locate DAQ header." << std::endl;
   }
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  evt_timeGPS_nsec = evtTimeGPS.timeLow();

  //std::cout<<"evt_timeGPS_nsec: "<< evt_timeGPS_nsec<<std::endl;

  if(T0Corr){
    if(crtT0_v.size() == 1){
      CRTT0corr = crtT0_v.front()->Time();
    }
    else{
      CRTT0corr = 0;
    }
  }

  for (unsigned int i_crt = 0; i_crt < crthit_v.size(); i_crt ++){
    double crt_time;
    if(T0Corr){
      crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.);
    }
    else{
      crt_time = ((crthit_v[i_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
    }

    crthit_PE.push_back(crthit_v[i_crt]->peshit);
    crthit_plane.push_back(crthit_v[i_crt]->plane);
    crthit_time.push_back(crt_time);

    // only count if the CRT hit has signal above 70 PE
    if(crt_time >= fBeamStart && crt_time <= fBeamEnd && crthit_v[i_crt]->peshit > 70){
      Nr_crthit_inBeam++;
    } // If CRT hit in Beam window

  }

  nr_pfp = pfParticle_v.size();

  for(unsigned int i_pfp = 0; i_pfp < pfParticle_v.size(); i_pfp++){
    auto pfp = pfParticle_v[i_pfp];
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      nr_pfp_nuDaughters = pfp->NumDaughters();
      for(int j = 0; j < nr_pfp_nuDaughters; j++){
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

      auto flash_matching_T0 = pfpToT0Asso.at(pfp.key());
      if(flash_matching_T0.size() == 1){
        flash_matching_chi2 = flash_matching_T0.front()->TriggerConfidence();
        trigger_time = flash_matching_T0.front()->Time();
      }
    }
  }
  nr_dau_tracks = daughter_Tracks.size();
  nr_dau_showers = daughter_Showers.size();

  //-------- Get All tracks in neutrino slice
  for(unsigned int i_trk = 0; i_trk < daughter_Tracks.size(); i_trk++){

    // Track CRT
    auto Track_CRThit = CRTToTrackAsso.at(daughter_Tracks[i_trk].key());
    auto Track_CRTT0 = CRTT0ToTrackAsso.at(daughter_Tracks[i_trk].key());
    if(Track_CRTT0.size() > 0){
      for(unsigned int i_trk_crt = 0; i_trk_crt < Track_CRThit.size(); i_trk_crt++){
        trk_crt_time_T0asso.push_back(Track_CRTT0[i_trk_crt]->Time());
      }
    }

    if(Track_CRThit.size() > 0){
      for(unsigned int i_trk_crt = 0; i_trk_crt < Track_CRThit.size(); i_trk_crt++){
        if(Track_CRThit[i_trk_crt]->peshit > 70){
          // w/o T0 correction
          trk_crt_time.push_back(((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.));
          //std::cout<<"Track_CRThit[i_trk_crt]->ts0_ns: "<< Track_CRThit[i_trk_crt]->ts0_ns<<std::endl;
          //std::cout<<"trk_crt_time: "<< trk_crt_time.back()<<std::endl;
          // w/ T0 correction
          if(T0Corr){
            trk_crt_time.push_back(((Track_CRThit[i_trk_crt]->ts0_ns - evt_timeGPS_nsec + CRTT0corr) / 1000.));
          }
        }
      }
    }
    Nr_trk_asso_crthit.push_back(trk_crt_time.size());

    // Add spatial correction to the track start and end
    Trk_start = daughter_Tracks[i_trk]->Start<TVector3>();
    auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
    Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X() + xtimeoffset + 0.6);
    Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
    Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

    trk_start_x.push_back(Trk_start_SCEcorr.X());
    trk_start_y.push_back(Trk_start_SCEcorr.Y());
    trk_start_z.push_back(Trk_start_SCEcorr.Z());

    trk_start_x_noSCE.push_back(Trk_start.X());
    trk_start_y_noSCE.push_back(Trk_start.Y());
    trk_start_z_noSCE.push_back(Trk_start.Z());

    Trk_end = daughter_Tracks[i_trk]->End<TVector3>();
    auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
    Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X() + xtimeoffset + 0.6);
    Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
    Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

    trk_end_x.push_back(Trk_end_SCEcorr.X());
    trk_end_y.push_back(Trk_end_SCEcorr.Y());
    trk_end_z.push_back(Trk_end_SCEcorr.Z());

    trk_end_x_noSCE.push_back(Trk_end.X());
    trk_end_y_noSCE.push_back(Trk_end.Y());
    trk_end_z_noSCE.push_back(Trk_end.Z());

    //-- track containment
    trk_contained.push_back(_fiducial_volume.TrackContain(Trk_start_SCEcorr, Trk_end_SCEcorr));
    vtx_InFV.push_back(_fiducial_volume.VertexInFV(Trk_start_SCEcorr));

    trk_contained_noSCE.push_back(_fiducial_volume.TrackContain(Trk_start, Trk_end));
    vtx_InFV_noSCE.push_back(_fiducial_volume.VertexInFV(Trk_start));
        
    //-- Angle wrt to the readout; theta is to the wire; phi is to x
    auto Trk_pos = Trk_end_SCEcorr - Trk_start_SCEcorr;
    double theta_pl2 = std::atan2(Trk_pos.Z(), Trk_pos.Y()); // atan2(y,x)
    double theta_pl1 = theta_pl2 + M_PI/3; // If plan1 is -60 degree to Y, looking from outside to the TPC
    double theta_pl0 = theta_pl2 - M_PI/3; // If plan0 is +60 degree to Y, looking from outside to the TPC
    double phi_readout = std::atan2(Trk_pos.X(), Trk_pos.Y());

    sin2_theta_pl2.push_back(sin(theta_pl2) * sin(theta_pl2));
    sin2_theta_pl1.push_back(sin(theta_pl1) * sin(theta_pl1));
    sin2_theta_pl0.push_back(sin(theta_pl0) * sin(theta_pl0));
    sin2_phi_readout.push_back(sin(phi_readout) * sin(phi_readout));

    //-- track angle
    trk_phi.push_back(daughter_Tracks[i_trk]->Phi());
    trk_theta.push_back(daughter_Tracks[i_trk]->Theta());
    trk_costheta.push_back(cos(daughter_Tracks[i_trk]->Theta()));

    //-- track length
    trk_length.push_back(daughter_Tracks[i_trk]->Length());


    //-- Momentum
    //- range
    mom_Range_mu.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length.back(), 13));
    mom_Range_p.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length.back(), 2212));
    mom_Range_pi.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length.back(), 211));

    //- MCS
    bestMCS.push_back(mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bestMomentum());
    bestMCSLL.push_back(mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bestLogLikelihood());
    fwdMCS.push_back(mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->fwdMomentum());
    fwdMCSLL.push_back(mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->fwdLogLikelihood());
    bwdMCS.push_back(mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bwdMomentum());
    bwdMCSLL.push_back(mcsfitresult_mu_v.at(daughter_Tracks[i_trk].key())->bwdLogLikelihood());

    bestMCS_bulk.push_back(mcsfitresult_mu_bulk_v.at(daughter_Tracks[i_trk].key())->bestMomentum());
    bestMCSLL_bulk.push_back(mcsfitresult_mu_bulk_v.at(daughter_Tracks[i_trk].key())->bestLogLikelihood());
    fwdMCS_bulk.push_back(mcsfitresult_mu_bulk_v.at(daughter_Tracks[i_trk].key())->fwdMomentum());
    fwdMCSLL_bulk.push_back(mcsfitresult_mu_bulk_v.at(daughter_Tracks[i_trk].key())->fwdLogLikelihood());
    bwdMCS_bulk.push_back(mcsfitresult_mu_bulk_v.at(daughter_Tracks[i_trk].key())->bwdMomentum());
    bwdMCSLL_bulk.push_back(mcsfitresult_mu_bulk_v.at(daughter_Tracks[i_trk].key())->bwdLogLikelihood());

    bestMCS_NoSCE.push_back(mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bestMomentum());
    bestMCSLL_NoSCE.push_back(mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bestLogLikelihood());
    fwdMCS_NoSCE.push_back(mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->fwdMomentum());
    fwdMCSLL_NoSCE.push_back(mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->fwdLogLikelihood());
    bwdMCS_NoSCE.push_back(mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bwdMomentum());
    bwdMCSLL_NoSCE.push_back(mcsfitresult_mu_NoSCE_v.at(daughter_Tracks[i_trk].key())->bwdLogLikelihood());

    std::cout<<"(1) bestMCS: "<< bestMCS.back()<<", bestMCSLL: "<< bestMCSLL.back()<<", fwd: "<< fwdMCS.back() << ", fwdLL: "<< fwdMCSLL.back()<<", bwd: "<<bwdMCS.back()<<", bwdLL: "<<bwdMCSLL.back()<<std::endl;
    std::cout<<"(2) bestMCS: "<< bestMCS_bulk.back()<<", bestMCSLL: "<< bestMCSLL_bulk.back()<<", fwd: "<< fwdMCS_bulk.back() << ", fwdLL: "<< fwdMCSLL_bulk.back()<<", bwd: "<<bwdMCS_bulk.back()<<", bwdLL: "<<bwdMCSLL_bulk.back()<<std::endl;
    std::cout<<"(3) bestMCS: "<< bestMCS_NoSCE.back()<<", bestMCSLL: "<< bestMCSLL_NoSCE.back()<<", fwd: "<< fwdMCS_NoSCE.back() << ", fwdLL: "<< fwdMCSLL_NoSCE.back()<<", bwd: "<<bwdMCS_NoSCE.back()<<", bwdLL: "<<bwdMCSLL_NoSCE.back()<<std::endl;

    //-- PID
     PID pid;
     PID_Chi2Mu_3pl.push_back(pid.PID_Chi2Mu_3pl); // Chi2 of muon assumption of 3 planes in PID
     PID_Chi2P_3pl.push_back(pid.PID_Chi2P_3pl); // Chi2 of proton assumption of 3 planes in PID
     PID_Chi2Pi_3pl.push_back(pid.PID_Chi2Pi_3pl); // Chi2 of pion assumption of 3 planes in PID
     PID_Chi2K_3pl.push_back(pid.PID_Chi2K_3pl); // Chi2 of kaon assumption of 3 planes in PID

     /////////////////
     // Fill TRUE info
     ////////////////
     if(IsMC){
       std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(daughter_Tracks[i_trk].key());
       BackTrackerTruthMatch backtrackertruthmatch;
       backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
       auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
       trk_cosmic_percent.push_back(backtrackertruthmatch.ReturnCosmicPercent());
       trk_purity.push_back(backtrackertruthmatch.ReturnPurity());
       trk_completeness.push_back(backtrackertruthmatch.ReturnCompleteness());

       if(!MCparticle){
         if_cosmic.push_back(true);
         if_matchPrimary.push_back(false);
         if_matchMu.push_back(false);

         //std::cout<<"MC particle does not exist!"<<std::endl;

         true_mom.push_back(-999);
         true_start_x.push_back(-999);
         true_start_y.push_back(-999);
         true_start_z.push_back(-999);
         true_end_x.push_back(-999);
         true_end_y.push_back(-999);
         true_end_z.push_back(-999);
         true_trk_ifcontained.push_back(-999);
         true_vtxFV.push_back(-999);
         true_trk_phi.push_back(-999);
         true_trk_theta.push_back(-999);
         true_trk_costheta.push_back(-999);
         true_trk_theta_yz.push_back(-999);
         true_trk_costheta_yz.push_back(-999);
         true_trk_theta_xz.push_back(-999);
         true_trk_costheta_xz.push_back(-999);
         true_trk_length.push_back(-999); // An estimation of true track length
         true_trk_PDG.push_back(-999);
       }
       else{
         if(trk_cosmic_percent.back() > 0.8){
           if_cosmic.push_back(true);
           if_matchPrimary.push_back(false);
           if_matchMu.push_back(false);
           
           true_mom.push_back(-999);
           true_start_x.push_back(-999);
           true_start_y.push_back(-999);
           true_start_z.push_back(-999);
           true_end_x.push_back(-999);
           true_end_y.push_back(-999);
           true_end_z.push_back(-999);
           true_trk_ifcontained.push_back(-999);
           true_vtxFV.push_back(-999);
           true_trk_phi.push_back(-999);
           true_trk_theta.push_back(-999);
           true_trk_costheta.push_back(-999);
           true_trk_theta_yz.push_back(-999);
           true_trk_costheta_yz.push_back(-999);
           true_trk_theta_xz.push_back(-999);
           true_trk_costheta_xz.push_back(-999);
           true_trk_length.push_back(-999); // An estimation of true track length
           true_trk_PDG.push_back(-999);

         }
         else{
           if_cosmic.push_back(false);
           if(MCparticle->Process() == "primary" && trk_purity.back() > 0.2 && trk_completeness.back() > 0.2){
             if_matchPrimary.push_back(true);
             if(MCparticle->PdgCode() == 13){
               if_matchMu.push_back(true);
             }
             else{
               if_matchMu.push_back(false);
             }
           }
           TVector3 TrueTrackPos(MCparticle->Px(), MCparticle->Py(), MCparticle->Pz());// The initial momentum represent the angle of true track
           true_mom.push_back(MCparticle->P());
           true_start_x.push_back(MCparticle->Position().X());
           true_start_y.push_back(MCparticle->Position().Y());
           true_start_z.push_back(MCparticle->Position().Z());
           true_end_x.push_back(MCparticle->EndPosition().X());
           true_end_y.push_back(MCparticle->EndPosition().Y());
           true_end_z.push_back(MCparticle->EndPosition().Z());
           TVector3 true_start(MCparticle->Position().X(), MCparticle->Position().Y(), MCparticle->Position().Z());
           TVector3 true_end(MCparticle->EndPosition().X(), MCparticle->EndPosition().Y(), MCparticle->EndPosition().Z());
           true_trk_ifcontained.push_back(_fiducial_volume.TrackContain(true_start, true_end));
           true_vtxFV.push_back(_fiducial_volume.VertexInFV(true_start));
           true_trk_phi.push_back(TrueTrackPos.Phi());
           true_trk_theta.push_back(TrueTrackPos.Theta());
           true_trk_costheta.push_back(cos(TrueTrackPos.Theta()));
           true_trk_theta_yz.push_back(std::atan2(TrueTrackPos.Y(), TrueTrackPos.Z()));
           true_trk_costheta_yz.push_back(cos(true_trk_theta_yz.back()));
           true_trk_theta_xz.push_back(std::atan2(TrueTrackPos.X(), TrueTrackPos.Z()));
           true_trk_costheta_xz.push_back(cos(true_trk_theta_xz.back()));
           true_trk_length.push_back((true_start - true_end).Mag()); // An estimation of true track length
           true_trk_PDG.push_back(MCparticle->PdgCode());
         }
       } // If MC particle exists
     } //IsMC

  } // End of track loop 

  my_event_->Fill();

  flash_YCenter = -999;
  flash_YWidth = -999;
  flash_ZCenter = -999;
  flash_ZWidth = -999;
  flash_TotalPE = -999;
  flash_time = -999;

  flash_matching_chi2 = -999;
  trigger_time = -999;

  evt_CRTveto = false;
  evt_CRTveto_100 = false;

  CRTT0corr = 0;
  crthit_PE.clear(); // The photonelectrons of CRT hits which are in beam window
  crthit_plane.clear(); // Plane of CRT hits
  crthit_time.clear(); // Time of CRT hits
  Nr_crthit_inBeam = 0; // Number of CRT hits in beamtime

  Nr_trk_asso_crthit.clear(); // Number of CRT hits associated to the track
  trk_crt_time.clear();// The CRT time of the hit which matched to the track
  trk_crt_time_T0asso.clear();

  daughter_Tracks.clear();
  daughter_Showers.clear();

  if_trk_NeutrinoSlice_primary.clear();

  trk_vtx_x.clear();
  trk_vtx_y.clear();
  trk_vtx_z.clear();
  trk_start_x.clear();
  trk_start_y.clear();
  trk_start_z.clear();
  trk_end_x.clear();
  trk_end_y.clear();
  trk_end_z.clear();

  trk_vtx_x_noSCE.clear();
  trk_vtx_y_noSCE.clear();
  trk_vtx_z_noSCE.clear();
  trk_start_x_noSCE.clear();
  trk_start_y_noSCE.clear();
  trk_start_z_noSCE.clear();
  trk_end_x_noSCE.clear();
  trk_end_y_noSCE.clear();
  trk_end_z_noSCE.clear();

  trk_contained.clear();
  vtx_InFV.clear();

  trk_contained_noSCE.clear();
  vtx_InFV_noSCE.clear();

  sin2_theta_pl0.clear();
  sin2_theta_pl1.clear();
  sin2_theta_pl2.clear();
  sin2_phi_readout.clear();

  trk_phi.clear();
  trk_theta.clear();
  trk_costheta.clear();

  trk_length.clear();

  mom_Range_mu.clear();
  mom_Range_p.clear();
  mom_Range_pi.clear();

  bestMCS.clear();
  bestMCSLL.clear();
  fwdMCS.clear();
  fwdMCSLL.clear();
  bwdMCS.clear();
  bwdMCSLL.clear();

  bestMCS_bulk.clear();
  bestMCSLL_bulk.clear();
  fwdMCS_bulk.clear();
  fwdMCSLL_bulk.clear();
  bwdMCS_bulk.clear();
  bwdMCSLL_bulk.clear();
 
  bestMCS_NoSCE.clear();
  bestMCSLL_NoSCE.clear();
  fwdMCS_NoSCE.clear();
  fwdMCSLL_NoSCE.clear();
  bwdMCS_NoSCE.clear();
  bwdMCSLL_NoSCE.clear();

  PID_Chi2Mu_3pl.clear();
  PID_Chi2P_3pl.clear();
  PID_Chi2Pi_3pl.clear();
  PID_Chi2K_3pl.clear();

  if(IsMC){
    true_mom.clear();
    true_start_x.clear();
    true_start_y.clear();
    true_start_z.clear();
    true_end_x.clear();
    true_end_y.clear();
    true_end_z.clear();
    true_trk_phi.clear();
    true_trk_theta.clear();
    true_trk_costheta.clear();
    true_trk_theta_yz.clear();
    true_trk_costheta_yz.clear();
    true_trk_theta_xz.clear();
    true_trk_costheta_xz.clear();
    true_trk_length.clear();
    true_trk_PDG.clear();
    true_trk_ifcontained.clear();
    true_vtxFV.clear();

    trk_cosmic_percent.clear();
    trk_purity.clear();
    trk_completeness.clear();
    if_cosmic.clear();
    if_matchPrimary.clear();
    if_matchMu.clear();
  }
       
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

  my_event_->Branch("flash_YCenter", &flash_YCenter);
  my_event_->Branch("flash_YWidth", &flash_YWidth);
  my_event_->Branch("flash_ZCenter", &flash_ZCenter);
  my_event_->Branch("flash_ZWidth", &flash_ZWidth);
  my_event_->Branch("flash_TotalPE", &flash_TotalPE);
  my_event_->Branch("flash_time", &flash_time);
  my_event_->Branch("flash_matching_chi2", &flash_matching_chi2);
  my_event_->Branch("trigger_time", &trigger_time);

  my_event_->Branch("evt_CRTveto", &evt_CRTveto);
  my_event_->Branch("evt_CRTveto_100", &evt_CRTveto_100);

  my_event_->Branch("crthit_PE", &crthit_PE);
  my_event_->Branch("crthit_plane", &crthit_plane);
  my_event_->Branch("crthit_time", &crthit_time);
  my_event_->Branch("Nr_crthit_inBeam", &Nr_crthit_inBeam);
  my_event_->Branch("Nr_trk_asso_crthit", &Nr_trk_asso_crthit);
  my_event_->Branch("trk_crt_time", &trk_crt_time);
  my_event_->Branch("trk_crt_time_T0asso", &trk_crt_time_T0asso);

  my_event_->Branch("nr_pfp", &nr_pfp);
  my_event_->Branch("nr_pfp_nuDaughters", &nr_pfp_nuDaughters);
  my_event_->Branch("nr_dau_tracks", &nr_dau_tracks);
  my_event_->Branch("nr_dau_showers", &nr_dau_showers);

  my_event_->Branch("if_trk_NeutrinoSlice_primary", &if_trk_NeutrinoSlice_primary);

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
  my_event_->Branch("bestMCS_bulk", &bestMCS_bulk);
  my_event_->Branch("bestMCSLL_bulk", &bestMCSLL_bulk);
  my_event_->Branch("fwdMCS_bulk", &fwdMCS_bulk);
  my_event_->Branch("fwdMCSLL_bulk", &fwdMCSLL_bulk);
  my_event_->Branch("bwdMCS_bulk", &bwdMCS_bulk);
  my_event_->Branch("bwdMCSLL_bulk", &bwdMCSLL_bulk);
  my_event_->Branch("bestMCS_NoSCE", &bestMCS_NoSCE);
  my_event_->Branch("bestMCSLL_NoSCE", &bestMCSLL_NoSCE);
  my_event_->Branch("fwdMCS_NoSCE", &fwdMCS_NoSCE);
  my_event_->Branch("fwdMCSLL_NoSCE", &fwdMCSLL_NoSCE);
  my_event_->Branch("bwdMCS_NoSCE", &bwdMCS_NoSCE);
  my_event_->Branch("bwdMCSLL_NoSCE", &bwdMCSLL_NoSCE);

  my_event_->Branch("PID_Chi2Mu_3pl", &PID_Chi2Mu_3pl);
  my_event_->Branch("PID_Chi2P_3pl", &PID_Chi2P_3pl);
  my_event_->Branch("PID_Chi2Pi_3pl", &PID_Chi2Pi_3pl);
  my_event_->Branch("PID_Chi2K_3pl", &PID_Chi2K_3pl);

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

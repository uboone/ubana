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
//#include "lardataobj/AnalysisBase/PlaneIDBitsetHelperFunctions.h"

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
  
  bool MC_beamNeutrino; // MCTruth beam origin
  bool MC_FV; // MCTruth vertex in FV = true, out of FV = false
  int MC_ccnc; // MCTruth cc = 0 or nc = 1
  int MC_nupdg; // MCTruth nupdg; numu = 14, nue = 12
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

  int Genie_nNeutron_preFSI;// before FSI 
  int Genie_nProton_preFSI;// before FSI 
  int Genie_nPi0_preFSI;// before FSI 
  int Genie_nPiPlus_preFSI;// before FSI 
  int Genie_nPiMinus_preFSI;// before FSI 

  int TopologyType;// The topology of true neutrino interaction + FSI products after Geant4

  double flash_matching_chi2; //Chi2 of flash matching in each neutrino slice

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

  bool evt_CRTveto = false; // If CRT veto, eliminate the events for contained (70PE threshold)
  bool evt_CRTveto_100 = false; // If CRT veto, eliminate the events for contained (100PE threshold)

  bool if_selected = false; // If selected based on the reco info
  bool if_matchMu = false; // If the selected track matched with true muon from numu cc
  bool if_cosmic = true; // Check if a track is cosmic or not by if it has an associated MCParticle

  bool if_broken = false; // if find broken track
  bool if_newTrkThroughGoing = false; // if the new track is through going
  std::vector<double> trk_broken_len;// length of the merged broken track (distance of two furtherest track ends)
  std::vector<int> trk_broken_nr_merged;// number of broken tracks merged in

  std::vector<double> mom_bestMCS_mu;//MCS best momentum of muon track in the every event
  std::vector<double> mom_bestMCS_ll_mu;//Likelihood of MCS best momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_mu;//MCS forward momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_ll_mu;//Likelihood of MCS forward momentum of muon track in the every event
  std::vector<double> mom_bwdMCS_mu;//MCS backward momentum of muon track in the every event
  std::vector<double> mom_bwdMCS_ll_mu;//Likelihood of MCS backward momentum of muon track in the every event
  std::vector<double> mom_Range_mu;//Range momentum of muon track in the every event
  std::vector<double> mom_Range_p;//Range momentum of proton track in the every event
  std::vector<double> mom_Range_pi;//Range momentum of pion track in the every event
  std::vector<double> mom_Range_mu_noSCE;//Range momentum of muon track in the every event
  std::vector<double> mom_Range_p_noSCE;//Range momentum of proton track in the every event
  std::vector<double> mom_Range_pi_noSCE;//Range momentum of pion track in the every event
  std::vector<double> mom_range_PID_avg_noSCE;//Range momentum of tracks based on their PID particle type using 3 pls

  std::vector<double> vtx_x;//Reconstructed vtx x in the every event
  std::vector<double> vtx_y;//Reconstructed vtx y in the every event
  std::vector<double> vtx_z;//Reconstructed vtx z in the every event
  std::vector<double> vtx_x_MCS;//Reconstructed vtx x in the every event
  std::vector<double> vtx_y_MCS;//Reconstructed vtx y in the every event
  std::vector<double> vtx_z_MCS;//Reconstructed vtx z in the every event
  std::vector<double> start_x;//Reconstructed start x in the every event
  std::vector<double> start_y;//Reconstructed start y in the every event
  std::vector<double> start_z;//Reconstructed start z in the every event
  std::vector<double> end_x;//Reconstructed end x in the every event
  std::vector<double> end_y;//Reconstructed end y in the every event
  std::vector<double> end_z;//Reconstructed end z in the every event
  std::vector<double> trk_phi;//Reconstructed track phi in the every event
  std::vector<double> trk_theta;//Reconstructed track theta in the every event
  std::vector<double> trk_costheta;//Reconstructed track cos(theta) in the every event
  std::vector<double> trk_phi_MCS;//Reconstructed track phi in the every event
  std::vector<double> trk_theta_MCS;//Reconstructed track theta in the every event
  std::vector<double> trk_costheta_MCS;//Reconstructed track cos(theta) in the every event
  std::vector<double> trk_length_pl0;//Range momentum of muon track in the every event
  std::vector<double> trk_length_pl1;//Range momentum of muon track in the every event
  std::vector<double> trk_length_pl2;//Range momentum of muon track in the every event
  std::vector<double> trk_length_avg;//Range momentum of muon track in the every event
  std::vector<double> trk_length_noSCE;//Range momentum of muon track in the every event
  std::vector<bool> trk_ifcontained;//to check if the track is contained or not
  std::vector<bool> trk_OutOfTime;//to check if either of the track end is out of X boundary
  std::vector<bool> vtx_FV;//to check if the vertex is in FV or not
  std::vector<bool> vtx_MCS_FV;//to check if the vertex is in FV or not
  std::vector<double> trk_end_theta_yz;
  std::vector<double> trk_end_costheta_yz;
  std::vector<double> trk_end_theta_xz;
  std::vector<double> trk_end_costheta_xz;
  std::vector<double> trk_theta_yz;
  std::vector<double> trk_costheta_yz;
  std::vector<double> trk_theta_xz;
  std::vector<double> trk_costheta_xz;

  std::vector<bool> old_trk_ifcontained;//to check if the track is contained or not
  std::vector<bool> old_vtx_FV;//to check if the vertex is in FV or not

  int n_pfp_nuDaughters; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters

  int hits_dEdx_size_pl0;
  int hits_dEdx_size_pl1;
  int hits_dEdx_size_pl2;

  std::vector<float> dEdx_pl0; // dE/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<float> dEdx_pl1; // dE/dx of the selected (muon) track from plane 1
  std::vector<float> dEdx_pl2; // dE/dx of the selected (muon) track from plane 2 (collection)
  std::vector<float> dQdx_pl0; // dQ/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<float> dQdx_pl1; // dQ/dx of the selected (muon) track from plane 1
  std::vector<float> dQdx_pl2; // dQ/dx of the selected (muon) track from plane 2 (collection)
  std::vector<float> resRange_pl0; // range from a hit to the end of the selected track end
  std::vector<float> resRange_pl1; // range from a hit to the end of the selected track end
  std::vector<float> resRange_pl2; // range from a hit to the end of the selected track end
  float dEdx_pl0_start_half; // average dEdx of start half hits of pl 0
  float dEdx_pl1_start_half; // average dEdx of start half hits of pl 0
  float dEdx_pl2_start_half; // average dEdx of start half hits of pl 0
  float dEdx_pl0_end_half; // average dEdx of end half hits of pl 0
  float dEdx_pl1_end_half; // average dEdx of end half hits of pl 0
  float dEdx_pl2_end_half; // average dEdx of end half hits of pl 0
  float dEdx_pl0_start10; // average dEdx of first 10 hit of pl 0
  float dEdx_pl1_start10; // average dEdx of first 10 hit of pl 0
  float dEdx_pl2_start10; // average dEdx of first 10 hit of pl 0
  float dEdx_pl0_end10; // average dEdx of end 10 hit of pl 0
  float dEdx_pl1_end10; // average dEdx of end 10 hit of pl 0
  float dEdx_pl2_end10; // average dEdx of end 10 hit of pl 0
  float dEdx_pl0_start1020; // average dEdx of first 10 hit of pl 0
  float dEdx_pl1_start1020; // average dEdx of first 10 hit of pl 0
  float dEdx_pl2_start1020; // average dEdx of first 10 hit of pl 0
  float dEdx_pl0_end1020; // average dEdx of end 10 hit of pl 0
  float dEdx_pl1_end1020; // average dEdx of end 10 hit of pl 0
  float dEdx_pl2_end1020; // average dEdx of end 10 hit of pl 0

  float dEdx_pl2_1020_ratio; // dEdx_pl2_end1020/(dEdx_pl2_end1020 + dEdx_pl2_start1020)

  double PID_Chi2Mu_pl0; // Chi2 of muon assumption of plane 0 in PID
  double PID_Chi2Mu_pl1; // Chi2 of muon assumption of plane 1 in PID
  double PID_Chi2Mu_pl2; // Chi2 of muon assumption of plane 2 in PID
  double PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID 
 
  double PID_Chi2P_pl0; // Chi2 of proton assumption of plane 0 in PID
  double PID_Chi2P_pl1; // Chi2 of proton assumption of plane 1 in PID
  double PID_Chi2P_pl2; // Chi2 of proton assumption of plane 2 in PID
  double PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID

  double PID_Chi2Pi_pl0; // Chi2 of pion assumption of plane 0 in PID
  double PID_Chi2Pi_pl1; // Chi2 of pion assumption of plane 1 in PID
  double PID_Chi2Pi_pl2; // Chi2 of pion assumption of plane 2 in PID
  double PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID
  
  double PID_Chi2K_pl0; // Chi2 of kaon assumption of plane 0 in PID
  double PID_Chi2K_pl1; // Chi2 of kaon assumption of plane 1 in PID
  double PID_Chi2K_pl2; // Chi2 of kaon assumption of plane 2 in PID
  double PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

  int PID_Pdg_3pl; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
  int PID_Pdg_pl2;
  int PID_Pdg_pl1;
  int PID_Pdg_pl0;
  double PID_avg_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions
  double PID_pl2_Chi2;
  double PID_pl1_Chi2;
  double PID_pl0_Chi2;

  int BestPlane_PID;
  bool Pl2_for_PID;
  bool Pl1_for_PID;
  bool Pl0_for_PID;

  bool if_fwd_MCS; // If using forward MCS direction judge
  bool if_fwd_true; // If fwd by the reco true vertex distance
  bool if_fwd_dEdx10; // If fwd by the reco dEdx 10 hits (should use for contained)
  bool if_fwd_dEdx1020; // If fwd by the reco dEdx 10 hits (should use for contained)
  bool if_fwd_dEdxhalf; // If fwd by the reco dEdx half of the hits (should use for contained)

  bool                                IsMC;
  bool                                UsingCRT;
  std::string                         m_generatorLabel;
  std::string                         m_geantLabel;
  std::string                         m_pandoraLabel;
  std::string                         m_T0ProducerLabel;
  std::string                         m_hitProducerLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_showerProducerLabel;
  std::string                         m_MCSmuProducerLabel;
  std::string                         m_calorimetryProducerLabel;
  std::string                         Hits_TrackAssLabel;
  std::string                         PID_TrackAssLabel;
  std::string                         m_CRTVetoLabel;
  std::string                         m_FlashLabel;

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


SingleMuon::SingleMuon(fhicl::ParameterSet const& pset)
  : 
  EDAnalyzer{pset},
  IsMC(pset.get<bool>("IsMC")),
  UsingCRT(pset.get<bool>("UsingCRT")),
  m_generatorLabel(pset.get<std::string>("GeneratorLabel")),
  m_geantLabel(pset.get<std::string>("GeantLabel")),
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_T0ProducerLabel(pset.get<std::string>("T0ProducerLabel")),
  m_hitProducerLabel(pset.get<std::string>("HitProducerLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel")),
  m_MCSmuProducerLabel(pset.get<std::string>("MCSmuProducerLabel")),
  m_calorimetryProducerLabel(pset.get<std::string>("calorimetryProducerLabel")),
  m_CRTVetoLabel(pset.get<std::string>("CRTVetoLabel")),
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
}

void SingleMuon::analyze(art::Event const& evt)
{

  //// Get necessary handles
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  std::vector<art::Ptr<simb::GTruth> > GTruthCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
  //art::FindMany<crt::CRTHit> CRThitFlashAsso;

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
  }

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

  //if(UsingCRT){
  //  // CRT Hit - Flash association
  //  CRThitFlashAsso(Handle_opflash, evt, m_CRTVetoLabel);
  //  //art::FindMany<crt::CRTHit> CRThitFlashAsso(Handle_opflash, evt, m_CRTVetoLabel);
  //}

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
  std::vector<int> Ghost_PDG; // pdg code of the pfp which has no track or shower associated; No elements ideally
  std::vector<float> vTrk_len;

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
    
    for(int i_mc = 0; i_mc < (int) MCTruthCollection.size(); i_mc++){
      if (MCTruthCollection[i_mc]->Origin() == Neutrino_Origin) MC_beamNeutrino = true;
      MC_nupdg = MCTruthCollection[i_mc]->GetNeutrino().Nu().PdgCode();
      MC_ccnc = MCTruthCollection[i_mc]->GetNeutrino().CCNC();

      MC_nuVtxX = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vx();
      MC_nuVtxY = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vy();
      MC_nuVtxZ = MCTruthCollection[i_mc]->GetNeutrino().Nu().Vz();
      TVector3 true_nuVtx(MC_nuVtxX, MC_nuVtxY, MC_nuVtxZ);
      MC_FV = _fiducial_volume.InFV(true_nuVtx);
    }

    // Loop all the MCParticles to determine the true topology (all the MCParticles are from the neutrino events in overlay)
    // Not necessary all the Genie particles go through the geant4 stage?
    if (MC_ccnc == 0 && MC_nupdg == 14 && MC_beamNeutrino == true){
      for(int i_mcp = 0; i_mcp < (int) MCParticleCollection.size(); i_mcp++){
        if(MCParticleCollection[i_mcp]->Process() == "primary"){
          // PDG and momemtum of neutrino daughters
          MC_Primary_PDG.push_back(MCParticleCollection[i_mcp]->PdgCode());
          MC_Primary_Mom.push_back(MCParticleCollection[i_mcp]->P());

          // muon
          if(MCParticleCollection[i_mcp]->PdgCode() == 13) MC_nMuon++;
          // electron
          if(MCParticleCollection[i_mcp]->PdgCode() == 11) MC_nElectron++;
          // neutron
          if(MCParticleCollection[i_mcp]->PdgCode() == 2112) MC_nNeutron++;
          // proton
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() < 0.255) MC_nProton_below255++;
          if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() >= 0.255) MC_nProton_above255++;
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
    
    // Get Genie info on how many particles produced
    for(int i_gn = 0; i_gn < (int) GTruthCollection.size(); i_gn++){
      Genie_nNeutron_preFSI = GTruthCollection[i_gn]->fNumNeutron;
      Genie_nProton_preFSI = GTruthCollection[i_gn]->fNumProton;
      Genie_nPi0_preFSI = GTruthCollection[i_gn]->fNumPi0;
      Genie_nPiPlus_preFSI = GTruthCollection[i_gn]->fNumPiPlus;
      Genie_nPiMinus_preFSI = GTruthCollection[i_gn]->fNumPiMinus;
    }
  }

  //-------- Get Reco neutrino (pfparticle)
  for(int i = 0; i < (int) pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      n_pfp_nuDaughters = pfp->NumDaughters();

      // Get Neutrino slice flash matching score (chi2)
      auto flash_matching_T0 = pfpToT0Asso.at(pfp.key());

      if(flash_matching_T0.size() == 1){
        flash_matching_chi2 = flash_matching_T0.front()->TriggerConfidence();
      }

      // For CC0pi0p, we only consider the case with the number of neutrino daughters less than 4
      if(n_pfp_nuDaughters < 4){
        // Get the pointer for the daughters of the neutrino
        for(int j = 0; j< n_pfp_nuDaughters; j++){
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

      //Todo: temperary version
      // Selection and Fill in Info
      if(n_dau_tracks == 1 && n_dau_showers == 0){

        if(UsingCRT){
          if(flash_v.size() > 0){
            for(int i_fl = 0; i_fl < (int) flash_v.size(); i_fl++){
              auto CRT_hit = CRThitFlashAsso.at(flash_v[i_fl].key());
              std::cout<<"CRT_hit.size(): "<< CRT_hit.size()<<std::endl;
              if(CRT_hit.size() == 1){
                evt_CRTveto = true;
                if(CRT_hit.front()->peshit > 100) evt_CRTveto_100 = true;
              } // if CRT veto
            } // loop over flash(es)
          } // if flash exists
        } // Using CRT

        //-- Fill RECO track info (in the naive version this is selected)
        if_selected = true;

        // Add spatial correction to the track start and end
        TVector3 Trk_start = daughter_Tracks.front()->Vertex<TVector3>();
        auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
        TVector3 Trk_start_SCEcorr;
        Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
        Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
        Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

        TVector3 Trk_end = daughter_Tracks.front()->End<TVector3>();
        auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
        TVector3 Trk_end_SCEcorr;
        Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X());
        Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
        Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

        //-- if either of the track end is out of the time, label them 
        if(Trk_start_SCEcorr.X() < 0 || Trk_start_SCEcorr.X() > 2. * geo->DetHalfWidth() || Trk_end_SCEcorr.X() < 0 || Trk_end_SCEcorr.X() > 2. * geo->DetHalfWidth()) trk_OutOfTime.push_back(true);
        else trk_OutOfTime.push_back(false);

        bool trk_contained = _fiducial_volume.InFV(Trk_start_SCEcorr, Trk_end_SCEcorr);
        old_trk_ifcontained.push_back(trk_contained);       
       
        bool vtx_InFV = _fiducial_volume.InFV(Trk_start_SCEcorr);
        old_vtx_FV.push_back(vtx_InFV);

        //-- Preliminary broken track searching (For the moment, only support 2 track merging)
        BrokenTrack brokentrack;
        brokentrack.MatchTracks(daughter_Tracks.front(), AllTrackCollection);
        if(brokentrack.NewTrk()){
          if_broken = true;
          trk_broken_len.push_back(brokentrack.TrkLen());
          trk_broken_nr_merged.push_back(brokentrack.NumberMergedTracks());

          TVector3 trk_end1 = brokentrack.TrkEnd1();
          TVector3 trk_end2 = brokentrack.TrkEnd2();

          trk_ifcontained.push_back(_fiducial_volume.InFV(trk_end1, trk_end2));
          
          if(!_fiducial_volume.InFV(trk_end1) && !_fiducial_volume.InFV(trk_end2)){ 
            vtx_FV.push_back(false);
            if_newTrkThroughGoing = true;
            std::cout<<"trk_end1_ X: "<< trk_end1.X()<<", Y: "<<trk_end1.Y()<<", Z: "<<trk_end1.Z()<<std::endl;
            std::cout<<"trk_end2_ X: "<< trk_end2.X()<<", Y: "<<trk_end2.Y()<<", Z: "<<trk_end2.Z()<<std::endl;
          }
          else vtx_FV.push_back(vtx_InFV);
        } 
        else{
          trk_ifcontained.push_back(trk_contained);
          vtx_FV.push_back(vtx_InFV);
        }

        double bestMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bestMomentum();
        double bestMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bestLogLikelihood();
        double fwdMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->fwdMomentum();
        double fwdMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->fwdLogLikelihood();
        double bwdMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bwdMomentum();
        double bwdMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bwdLogLikelihood();

        mom_bestMCS_mu.push_back(bestMCS);
        mom_bestMCS_ll_mu.push_back(bestMCSLL);
        mom_fwdMCS_mu.push_back(fwdMCS);
        mom_fwdMCS_ll_mu.push_back(fwdMCSLL);
        mom_bwdMCS_mu.push_back(bwdMCS);
        mom_bwdMCS_ll_mu.push_back(bwdMCSLL);

        // Track length and range momentum
        double Trk_length_noSCE = daughter_Tracks.front()->Length();
        auto assoCal = trackToCalAsso.at(daughter_Tracks.front().key());
        double Trk_length_pl0 = assoCal[0]->Range();  //pandoracali has spatial correction
        double Trk_length_pl1 = assoCal[1]->Range();  //pandoracali has spatial correction
        double Trk_length_pl2 = assoCal[2]->Range();  //pandoracali has spatial correction
        double Trk_length_avg = 0;
        int valid_pl = 0;
        for (int i_pl = 0; i_pl < (int) assoCal.size(); i_pl++){
          if(assoCal[i_pl]->Range() > 0){
            Trk_length_avg += assoCal[i_pl]->Range();
            valid_pl++;
          }
        }
        Trk_length_avg = Trk_length_avg / valid_pl;

        trk_length_pl0.push_back(Trk_length_pl0); // track length with spatial correction
        trk_length_pl1.push_back(Trk_length_pl1); // track length with spatial correction
        trk_length_pl2.push_back(Trk_length_pl2); // track length with spatial correction
        trk_length_avg.push_back(Trk_length_avg); // track length with spatial correction
        trk_length_noSCE.push_back(Trk_length_noSCE);

        // Usual case, the vertex is the single track start
        vtx_x.push_back(Trk_start_SCEcorr.X());
        vtx_y.push_back(Trk_start_SCEcorr.Y());
        vtx_z.push_back(Trk_start_SCEcorr.Z());
        start_x.push_back(Trk_start_SCEcorr.X());
        start_y.push_back(Trk_start_SCEcorr.Y());
        start_z.push_back(Trk_start_SCEcorr.Z());
        end_x.push_back(Trk_end_SCEcorr.X());
        end_y.push_back(Trk_end_SCEcorr.Y());
        end_z.push_back(Trk_end_SCEcorr.Z());
        
        trk_phi.push_back(daughter_Tracks.front()->Phi());
        trk_theta.push_back(daughter_Tracks.front()->Theta());
        trk_costheta.push_back(cos(daughter_Tracks.front()->Theta()));
    
        double Mom_Range_mu = _trk_mom_calculator.GetTrackMomentum(Trk_length_pl2, 13);
        double Mom_Range_mu_noSCE = _trk_mom_calculator.GetTrackMomentum(Trk_length_noSCE, 13);
        mom_Range_mu.push_back(Mom_Range_mu);
        mom_Range_mu_noSCE.push_back(Mom_Range_mu_noSCE);

        double Mom_Range_p = _trk_mom_calculator.GetTrackMomentum(Trk_length_pl2, 2212);
        double Mom_Range_p_noSCE = _trk_mom_calculator.GetTrackMomentum(Trk_length_noSCE, 2212);
        mom_Range_p.push_back(Mom_Range_p);
        mom_Range_p_noSCE.push_back(Mom_Range_p_noSCE);

        double Mom_Range_pi = _trk_mom_calculator.GetTrackMomentum(Trk_length_pl2, 211);
        double Mom_Range_pi_noSCE = _trk_mom_calculator.GetTrackMomentum(Trk_length_noSCE, 211);
        mom_Range_pi.push_back(Mom_Range_pi);
        mom_Range_pi_noSCE.push_back(Mom_Range_pi_noSCE);

        //Calorimetry Info
        // pandoracaliSCE has E-field and spatial correction
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
 
        hits_dEdx_size_pl0 = dEdx_pl0.size();
        hits_dEdx_size_pl1 = dEdx_pl1.size();
        hits_dEdx_size_pl2 = dEdx_pl2.size();
  
        int half_size_pl0 = hits_dEdx_size_pl0 / 2;
        int half_size_pl1 = hits_dEdx_size_pl1 / 2;
        int half_size_pl2 = hits_dEdx_size_pl2 / 2;
 
        dEdx_pl0_start_half = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
        dEdx_pl0_end_half = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
        dEdx_pl1_start_half = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
        dEdx_pl1_end_half = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
        dEdx_pl2_start_half = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
        dEdx_pl2_end_half = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;
        // dEdx_10
        if (dEdx_pl0.size()<=10) {
          dEdx_pl0_start10 = dEdx_pl0_start_half;
          dEdx_pl0_end10 = dEdx_pl0_end_half;
        }
        else{ 
          dEdx_pl0_start10 = std::accumulate(dEdx_pl0.end() - 10, dEdx_pl0.end(), 0.) / 10.;
          dEdx_pl0_end10 = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + 10, 0.) / 10.;
        }
        if (dEdx_pl1.size()<=10) {
          dEdx_pl1_start10 = dEdx_pl1_start_half;
          dEdx_pl1_end10 = dEdx_pl1_end_half;
        }
        else{
          dEdx_pl1_start10 = std::accumulate(dEdx_pl1.end() - 10, dEdx_pl1.end(), 0.) / 10.;
          dEdx_pl1_end10 = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + 10, 0.) / 10.;
        }
        if (dEdx_pl2.size()<=10) {
          dEdx_pl2_start10 = dEdx_pl2_start_half;
          dEdx_pl2_end10 = dEdx_pl2_end_half;
        }
        else{
          dEdx_pl2_start10 = std::accumulate(dEdx_pl2.end() - 10, dEdx_pl2.end(), 0.) / 10.;
          dEdx_pl2_end10 = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + 10, 0.) / 10.;
        }
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

        dEdx_pl2_1020_ratio = dEdx_pl2_end1020 / (dEdx_pl2_end1020 + dEdx_pl2_start1020);

        // Gain PID info of the track
        if(!PIDTotrackAsso.isValid()){
          throw cet::exception("[Numu0pi0p]") << "No matched PID - track information!" << std::endl;
        }
        // Get projected angle wrt to the wires (docdb 23008)
        TVector3 End_Dir = daughter_Tracks.front()->EndDirection<TVector3>();
        trk_end_theta_yz.push_back(std::atan2(End_Dir.Y(), End_Dir.Z()));
        trk_end_costheta_yz.push_back(cos(std::atan2(End_Dir.Y(), End_Dir.Z())));
        trk_end_theta_xz.push_back(std::atan2(End_Dir.X(), End_Dir.Z()));
        trk_end_costheta_xz.push_back(cos(std::atan2(End_Dir.X(), End_Dir.Z())));
  
        auto Trk_pos = Trk_end_SCEcorr - Trk_start_SCEcorr;
        trk_theta_yz.push_back(std::atan2(Trk_pos.Y(), Trk_pos.Z()));
        trk_costheta_yz.push_back(cos(std::atan2(Trk_pos.Y(), Trk_pos.Z())));
        trk_theta_xz.push_back(std::atan2(Trk_pos.X(), Trk_pos.Z()));
        trk_costheta_xz.push_back(cos(std::atan2(Trk_pos.X(), Trk_pos.Z())));

        //double theta_pl2 = std::atan2(End_Dir.Z(), End_Dir.Y()); // atan2(y,x)
        double theta_pl2 = std::atan2(Trk_pos.Z(), Trk_pos.Y()); // atan2(y,x)
        double theta_pl1 = theta_pl2 + M_PI/3; // If plan1 is -60 degree to Y, looking from outside to the TPC
        double theta_pl0 = theta_pl2 - M_PI/3; // If plan0 is +60 degree to Y, looking from outside to the TPC
        int w2 = 0; int w1 = 0; int w0 = 0;
        double sin2_pl2 = sin(theta_pl2) * sin(theta_pl2);
        double sin2_pl1 = sin(theta_pl1) * sin(theta_pl1);
        double sin2_pl0 = sin(theta_pl0) * sin(theta_pl0);
        if (sin2_pl2 >= 0.5) w2 = 1;
        if (sin2_pl1 >= 0.5) w1 = 1;
        if (sin2_pl0 >= 0.5) w0 = 1;

        // PID
        auto trkPID = PIDTotrackAsso.at(daughter_Tracks.front().key());
        if (trkPID.size() == 0){
          std::cout << "No PID information for this selected track!" << std::endl;
        }
        std::vector<anab::sParticleIDAlgScores> vAlg_PID = trkPID.front()->ParticleIDAlgScores();
        double PIDChi2_mu[3] = {-999,-999,-999};
        double PIDChi2_p[3] = {-999,-999,-999};
        double PIDChi2_pi[3] = {-999,-999,-999};
        double PIDChi2_K[3] = {-999,-999,-999};
        for(int i_Alg_PID = 0; i_Alg_PID < (int) vAlg_PID.size(); i_Alg_PID++){
          anab::sParticleIDAlgScores Alg_PID = vAlg_PID.at(i_Alg_PID);
          for(int id_pl = 0; id_pl < 3; id_pl++){
            if (Alg_PID.fPlaneMask.test(id_pl) && Alg_PID.fAlgName == "Chi2"){
               if (abs(Alg_PID.fAssumedPdg) == 13){ // muon
                 PIDChi2_mu[id_pl] = Alg_PID.fValue;
               }
               if (abs(Alg_PID.fAssumedPdg) == 2212){ // proton
                 PIDChi2_p[id_pl] = Alg_PID.fValue;
               }
               if (abs(Alg_PID.fAssumedPdg) == 211){ // pion
                 PIDChi2_pi[id_pl] = Alg_PID.fValue;
               }
               if (abs(Alg_PID.fAssumedPdg) == 321){ // kaon
                 PIDChi2_K[id_pl] = Alg_PID.fValue;
               }
            }
          }
        }
	PID_Chi2Mu_pl0 = PIDChi2_mu[0];
	PID_Chi2Mu_pl1 = PIDChi2_mu[1];
	PID_Chi2Mu_pl2 = PIDChi2_mu[2];
        int ww0 = w0; int ww1 = w1; int ww2 = w2; // copy from the origin
        if (PID_Chi2Mu_pl0 < 0) ww0 = 0;
        if (PID_Chi2Mu_pl1 < 0) ww1 = 0;
        if (PID_Chi2Mu_pl2 < 0) ww2 = 0;
        if (ww0 + ww1 + ww2 == 0) PID_Chi2Mu_3pl = -999;
        else PID_Chi2Mu_3pl = (ww0 * PID_Chi2Mu_pl0 + ww1 * PID_Chi2Mu_pl1 + ww2 * PID_Chi2Mu_pl2) / (ww0 + ww1 + ww2);

        PID_Chi2P_pl0 = PIDChi2_p[0];
        PID_Chi2P_pl1 = PIDChi2_p[1];
        PID_Chi2P_pl2 = PIDChi2_p[2];
        ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
        if (PID_Chi2P_pl0 < 0) ww0 = 0;
        if (PID_Chi2P_pl1 < 0) ww1 = 0;
        if (PID_Chi2P_pl2 < 0) ww2 = 0;
        if (ww0 + ww1 + ww2 == 0) PID_Chi2P_3pl = -999;
        else PID_Chi2P_3pl = (ww0 * PID_Chi2P_pl0 + ww1 * PID_Chi2P_pl1 + ww2 * PID_Chi2P_pl2) / (ww0 + ww1 + ww2);

        PID_Chi2Pi_pl0 = PIDChi2_pi[0];
        PID_Chi2Pi_pl1 = PIDChi2_pi[1];
        PID_Chi2Pi_pl2 = PIDChi2_pi[2];
        ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
        if (PID_Chi2Pi_pl0 < 0) ww0 = 0;
        if (PID_Chi2Pi_pl1 < 0) ww1 = 0;
        if (PID_Chi2Pi_pl2 < 0) ww2 = 0;
        if (ww0 + ww1 + ww2 == 0) PID_Chi2Pi_3pl = -999;
        else PID_Chi2Pi_3pl = (ww0 * PID_Chi2Pi_pl0 + ww1 * PID_Chi2Pi_pl1 + ww2 * PID_Chi2Pi_pl2) / (ww0 + ww1 + ww2);

        PID_Chi2K_pl0 = PIDChi2_K[0];
        PID_Chi2K_pl1 = PIDChi2_K[1];
        PID_Chi2K_pl2 = PIDChi2_K[2];
        ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
        if (PID_Chi2K_pl0 < 0) ww0 = 0;
        if (PID_Chi2K_pl1 < 0) ww1 = 0;
        if (PID_Chi2K_pl2 < 0) ww2 = 0;
        if (ww0 + ww1 + ww2 == 0) PID_Chi2K_3pl = -999;
        else PID_Chi2K_3pl = (ww0 * PID_Chi2K_pl0 + ww1 * PID_Chi2K_pl1 + ww2 * PID_Chi2K_pl2) / (ww0 + ww1 + ww2);

        std::vector<int> Nhit_3pl = {hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2};
        BestPlane_PID = std::max_element(Nhit_3pl.begin(), Nhit_3pl.end()) - Nhit_3pl.begin();
        if (BestPlane_PID == 0){
          Pl0_for_PID = true;
          Pl1_for_PID = false;
          Pl2_for_PID = false;
        }
        if (BestPlane_PID == 1){
          Pl0_for_PID = false;
          Pl1_for_PID = true;
          Pl2_for_PID = false;
        }
        if (BestPlane_PID == 2){
          Pl0_for_PID = false;
          Pl1_for_PID = false;
          Pl2_for_PID = true;
        }

        //-- Get minimum Chi2 and there corresponding particle type
        std::vector<double> PIDChi2_avg;// It follows the order of muon, proton, pion, kaon
        if (PID_Chi2Mu_3pl < 0) PIDChi2_avg.push_back(9999);
        else PIDChi2_avg.push_back(PID_Chi2Mu_3pl);
        if (PID_Chi2P_3pl < 0) PIDChi2_avg.push_back(9999);
        else PIDChi2_avg.push_back(PID_Chi2P_3pl);
        if (PID_Chi2Pi_3pl < 0) PIDChi2_avg.push_back(9999);
        else PIDChi2_avg.push_back(PID_Chi2Pi_3pl);
        if (PID_Chi2K_3pl < 0) PIDChi2_avg.push_back(9999);
        else PIDChi2_avg.push_back(PID_Chi2K_3pl);

        PID_avg_Chi2 = *std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end());
        int ID_PID = std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end()) - PIDChi2_avg.begin();
        if (ID_PID == 0) PID_Pdg_3pl = 13;
        if (ID_PID == 1) PID_Pdg_3pl = 2212;
        if (ID_PID == 2) PID_Pdg_3pl = 211;
        if (ID_PID == 3) PID_Pdg_3pl = 321;

        mom_range_PID_avg_noSCE.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_length_noSCE, PID_Pdg_3pl));
       
        // Use Plane 2
        std::vector<double> PIDChi2_pl2;// It follows the order of muon, proton, pion, kaon
        if (PIDChi2_mu[2] < 0) PIDChi2_pl2.push_back(9999);
        else PIDChi2_pl2.push_back(PIDChi2_mu[2]);
        if (PIDChi2_p[2] < 0) PIDChi2_pl2.push_back(9999);
        else PIDChi2_pl2.push_back(PIDChi2_p[2]);
        if (PIDChi2_pi[2] < 0) PIDChi2_pl2.push_back(9999);
        else PIDChi2_pl2.push_back(PIDChi2_pi[2]);
        if (PIDChi2_K[2] < 0) PIDChi2_pl2.push_back(9999);
        else PIDChi2_pl2.push_back(PIDChi2_K[2]);
  
        PID_pl2_Chi2 = *std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end());
        int ID_PID_pl2 = std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end()) - PIDChi2_pl2.begin();
        if (ID_PID_pl2 == 0) PID_Pdg_pl2 = 13;
        if (ID_PID_pl2 == 1) PID_Pdg_pl2 = 2212;
        if (ID_PID_pl2 == 2) PID_Pdg_pl2 = 211;
        if (ID_PID_pl2 == 3) PID_Pdg_pl2 = 321;

        // Use plane 1
        std::vector<double> PIDChi2_pl1;// It follows the order of muon, proton, pion, kaon
        if (PIDChi2_mu[1] < 0) PIDChi2_pl1.push_back(9999);
        else PIDChi2_pl1.push_back(PIDChi2_mu[1]);
        if (PIDChi2_p[1] < 0) PIDChi2_pl1.push_back(9999);
        else PIDChi2_pl1.push_back(PIDChi2_p[1]);
        if (PIDChi2_pi[1] < 0) PIDChi2_pl1.push_back(9999);
        else PIDChi2_pl1.push_back(PIDChi2_pi[1]);
        if (PIDChi2_K[1] < 0) PIDChi2_pl1.push_back(9999);
        else PIDChi2_pl1.push_back(PIDChi2_K[1]);

        PID_pl1_Chi2 = *std::min_element(PIDChi2_pl1.begin(), PIDChi2_pl1.end());
        int ID_PID_pl1 = std::min_element(PIDChi2_pl1.begin(), PIDChi2_pl1.end()) - PIDChi2_pl1.begin();
        if (ID_PID_pl1 == 0) PID_Pdg_pl1 = 13;
        if (ID_PID_pl1 == 1) PID_Pdg_pl1 = 2212;
        if (ID_PID_pl1 == 2) PID_Pdg_pl1 = 211;
        if (ID_PID_pl1 == 3) PID_Pdg_pl1 = 321;

        // Use plane 0
        std::vector<double> PIDChi2_pl0;// It follows the order of muon, proton, pion, kaon
        if (PIDChi2_mu[0] < 0) PIDChi2_pl0.push_back(9999);
        else PIDChi2_pl0.push_back(PIDChi2_mu[0]);
        if (PIDChi2_p[0] < 0) PIDChi2_pl0.push_back(9999);
        else PIDChi2_pl0.push_back(PIDChi2_p[0]);
        if (PIDChi2_pi[0] < 0) PIDChi2_pl0.push_back(9999);
        else PIDChi2_pl0.push_back(PIDChi2_pi[0]);
        if (PIDChi2_K[0] < 0) PIDChi2_pl0.push_back(9999);
        else PIDChi2_pl0.push_back(PIDChi2_K[0]);
  
        PID_pl0_Chi2 = *std::min_element(PIDChi2_pl0.begin(), PIDChi2_pl0.end());
        int ID_PID_pl0 = std::min_element(PIDChi2_pl0.begin(), PIDChi2_pl0.end()) - PIDChi2_pl0.begin();
        if (ID_PID_pl0 == 0) PID_Pdg_pl0 = 13;
        if (ID_PID_pl0 == 1) PID_Pdg_pl0 = 2212;
        if (ID_PID_pl0 == 2) PID_Pdg_pl0 = 211;
        if (ID_PID_pl0 == 3) PID_Pdg_pl0 = 321;
 
        //-- Fill TRUE info, if the track is from numu cc muon
        if(IsMC){
          std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(daughter_Tracks.front().key());
          BackTrackerTruthMatch backtrackertruthmatch;
          backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
          auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
          if(!MCparticle){
            if_cosmic = true;
            std::cout<<"MC particle does not exist!"<<std::endl;
          }
          else{
            if_cosmic = false;
            if(MCparticle->PdgCode() == 13) if_matchMu = true;
            auto TrueTrackPos = MCparticle->EndPosition() - MCparticle->Position();
            true_mom.push_back(MCparticle->P());
            true_start_x.push_back(MCparticle->Position().X());
            true_start_y.push_back(MCparticle->Position().Y());
            true_start_z.push_back(MCparticle->Position().Z());
            true_end_x.push_back(MCparticle->EndPosition().X());
            true_end_y.push_back(MCparticle->EndPosition().Y());
            true_end_z.push_back(MCparticle->EndPosition().Z());
            TVector3 true_start(MCparticle->Position().X(), MCparticle->Position().Y(), MCparticle->Position().Z());
            TVector3 true_end(MCparticle->EndPosition().X(), MCparticle->EndPosition().Y(), MCparticle->EndPosition().Z());
            true_trk_ifcontained.push_back(_fiducial_volume.InFV(true_start, true_end));
            true_vtxFV.push_back(_fiducial_volume.InFV(true_start));
            true_trk_phi.push_back(TrueTrackPos.Phi());
            true_trk_theta.push_back(TrueTrackPos.Theta());
            true_trk_costheta.push_back(cos(TrueTrackPos.Theta()));
            true_trk_theta_yz.push_back(std::atan2(TrueTrackPos.Y(), TrueTrackPos.Z()));
            true_trk_costheta_yz.push_back(cos(true_trk_theta_yz.back()));
            true_trk_theta_xz.push_back(std::atan2(TrueTrackPos.X(), TrueTrackPos.Z()));
            true_trk_costheta_xz.push_back(cos(true_trk_theta_xz.back()));
            true_trk_length.push_back(sqrt(TrueTrackPos.X()*TrueTrackPos.X() + TrueTrackPos.Y()*TrueTrackPos.Y() + TrueTrackPos.Z()*TrueTrackPos.Z())); // An estimation of true track length
            true_trk_PDG.push_back(MCparticle->PdgCode());
            
          }
        }
       
        //Directional Info
        // Check the directional info of the track by MCS
        if (mom_bestMCS_ll_mu == mom_fwdMCS_ll_mu) if_fwd_MCS = true;
        else if_fwd_MCS = false;

        // Check the direction info of the track by dEdx
        if (dEdx_pl2_start_half < dEdx_pl2_end_half) if_fwd_dEdxhalf = true;
        else if_fwd_dEdxhalf = false;
        if (dEdx_pl2_start10 < dEdx_pl2_end10) if_fwd_dEdx10 = true;
        else if_fwd_dEdx10 = false;
        if (dEdx_pl2_start1020 < dEdx_pl2_end1020) if_fwd_dEdx1020 = true;
        else if_fwd_dEdx1020 = false;

        if(IsMC){
          // Check the directional info of the track by true reco vertex distance
          TVector3 vtx(true_start_x[0], true_start_y[0], true_start_z[0]);
          TVector3 D_start = Trk_start_SCEcorr - vtx;
          TVector3 D_end = Trk_end_SCEcorr - vtx;
          if (D_start.Mag() < D_end.Mag()) if_fwd_true = true;
          else if_fwd_true = false;
        }

        if(if_fwd_MCS){
          vtx_x_MCS.push_back(vtx_x.back());
          vtx_y_MCS.push_back(vtx_y.back());
          vtx_z_MCS.push_back(vtx_z.back());
  
          vtx_MCS_FV.push_back(vtx_FV.back());
       
          trk_phi_MCS.push_back(trk_phi.back());
          trk_theta_MCS.push_back(trk_theta.back());
          trk_costheta_MCS.push_back(trk_costheta.back());
        }
        else{
          vtx_x_MCS.push_back(end_x.back());
          vtx_y_MCS.push_back(end_y.back());
          vtx_z_MCS.push_back(end_z.back());

          TVector3 vtx_MCS(end_x.back(), end_y.back(), end_z.back());          
          vtx_MCS_FV.push_back(_fiducial_volume.InFV(vtx_MCS));

          trk_phi_MCS.push_back(M_PI + trk_phi.back());
          trk_theta_MCS.push_back(M_PI + trk_theta.back());
          trk_costheta_MCS.push_back(cos(trk_theta_MCS.back()));
        }
      }
      ///////////Plan of implementation.......
      //
      //////1 track + 0 shower
      // What is the direction of the track? (dE/dx ramp up for bragg peak, MCS likelihood, vertex activity? If there's any daughter of the track? What are the daughters?)
      // If the vertex is out the TPC, cosmic.
      // If the vertex is in the TPC, can it be broken track? CRT info? nearby track?
      // Is this track actually an integration of different tracks (dE/dx discontinuity)
      // If yes, are the two tracks going same direction? (Is it a daughter of muon)
      // Is it or are they muon-like from PID?
      /////If at least one of them is muon-like, futher investigation!
      //
      //////1 track + 1 shower
      // Is the shower michel electron like (low energy < 60 MeV?)?
      // Is the track muon-like?
      // What is the direction? Similar checks as the previous (1 track + 0 shower) case..
      //// Only process to next step if this topology is from a muon decay to michel electron which makes a shower and the muon has a vertex in the TPC FV
      //
      ///////1 track + 2 shower
      // What is the energy of the showers? assuming large shower can only comes from neutrino interaction..which is not my topology
      // If it's weirdly low energy..check the track following the procedure above..Determine vertix position!
      //
      //
      //////2 tracks + 0 shower
      // Check the colinearity of the 2 tracks? 
      // If colinear, are they the same particle or not? Where's the vertex
      // If not colinear, Is one of them muon-like? The other electron-like (daughter of muon)? Or do they like cosmic coincidence?
      //
      /////2 tracks + 1 shower
      // Check what's the shower energy?
      // If the shower has weirdly low energy, perform the above procedures to check if it worth to investigate further
      // If yes, check the shower position wrt the position of tracks..
      //
      ///// 3 tracks
      // Check the colinearity of all three of them and either two of them.
      // If there's at least a pair of colinear tracks, check the domenancy of the pair tracks and one track.
      // This could be cosmic coincidence. Check if it’s actually two tracks which one of them is cosmic
      //                        
      ///// 4 tracks
      // Is this actually cosmic coincidence…. Actually two tracks cross like four tracks? And one is neutrino induced muon and the other is cosmic induced muon?
      //
      ///////////////

      // If there's only one track, is it muon?
      // What is the direction? 
      /////////// (select if it's a muon) PID needed here
      // Are the rest showers or ghost in the wrong positions of the family tree?
      ////////// If they are electrons? (PID maybe)
      ////////// The vertex activity, the shower vs track length, the shower vs track momentum
    }
  }

  my_event_->Fill();

  evt_CRTveto = false;
  evt_CRTveto_100 = false;

  if_cosmic = true;
  if_matchMu = false;
  if_selected = false;
  if_broken = false;
  if_newTrkThroughGoing = false;

  if(IsMC){
    MC_Primary_PDG.clear();
    MC_Primary_Mom.clear();
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
  }

  daughter_Tracks.clear();
  daughter_Showers.clear();
  Track_PDG.clear();
  Shower_PDG.clear();

  trk_broken_len.clear();
  trk_broken_nr_merged.clear();

  mom_bestMCS_mu.clear();
  mom_bestMCS_ll_mu.clear();
  mom_fwdMCS_mu.clear();
  mom_fwdMCS_ll_mu.clear();
  mom_bwdMCS_mu.clear();
  mom_bwdMCS_ll_mu.clear();
  mom_Range_mu.clear();
  mom_Range_mu_noSCE.clear();
  mom_Range_p.clear();
  mom_Range_p_noSCE.clear();
  mom_Range_pi.clear();
  mom_Range_pi_noSCE.clear();
  mom_range_PID_avg_noSCE.clear();

  vtx_x.clear();
  vtx_y.clear();
  vtx_z.clear();
  vtx_x_MCS.clear();
  vtx_y_MCS.clear();
  vtx_z_MCS.clear();
  start_x.clear();
  start_y.clear();
  start_z.clear();
  end_x.clear();
  end_y.clear();
  end_z.clear();
  trk_phi.clear();
  trk_theta.clear();
  trk_costheta.clear();
  trk_phi_MCS.clear();
  trk_theta_MCS.clear();
  trk_costheta_MCS.clear();
  trk_length_pl0.clear();
  trk_length_pl1.clear();
  trk_length_pl2.clear();
  trk_length_avg.clear();
  trk_length_noSCE.clear();
  trk_ifcontained.clear();
  trk_OutOfTime.clear();
  vtx_FV.clear();
  vtx_MCS_FV.clear();
  trk_end_theta_yz.clear();
  trk_end_costheta_yz.clear();
  trk_end_theta_xz.clear();
  trk_end_costheta_xz.clear();
  trk_theta_yz.clear();
  trk_costheta_yz.clear();
  trk_theta_xz.clear();
  trk_costheta_xz.clear();
  
  old_trk_ifcontained.clear();
  old_vtx_FV.clear();

  dEdx_pl0.clear();
  dEdx_pl1.clear();
  dEdx_pl2.clear();
  dQdx_pl0.clear();
  dQdx_pl1.clear();
  dQdx_pl2.clear();
  resRange_pl0.clear();
  resRange_pl1.clear();
  resRange_pl2.clear();
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
    my_event_->Branch("TopologyType", &TopologyType);
    my_event_->Branch("MC_beamNeutrino", &MC_beamNeutrino);
    my_event_->Branch("MC_nupdg", &MC_nupdg);
    my_event_->Branch("MC_ccnc", &MC_ccnc);
    my_event_->Branch("MC_nuVtxX", &MC_nuVtxX);
    my_event_->Branch("MC_nuVtxY", &MC_nuVtxY);
    my_event_->Branch("MC_nuVtxZ", &MC_nuVtxZ);
    my_event_->Branch("MC_FV", &MC_FV);
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

  my_event_->Branch("n_dau_tracks", &n_dau_tracks);
  my_event_->Branch("n_dau_showers", &n_dau_showers);
  
  my_event_->Branch("flash_matching_chi2", &flash_matching_chi2);

  my_event_->Branch("evt_CRTveto", &evt_CRTveto);
  my_event_->Branch("evt_CRTveto_100", &evt_CRTveto_100);

  my_event_->Branch("if_cosmic", &if_cosmic);
  my_event_->Branch("if_matchMu", &if_matchMu);
  my_event_->Branch("if_selected", &if_selected);
  my_event_->Branch("if_broken", &if_broken);
  my_event_->Branch("if_newTrkThroughGoing", &if_newTrkThroughGoing);

  my_event_->Branch("trk_broken_len", &trk_broken_len);
  my_event_->Branch("trk_broken_nr_merged", &trk_broken_nr_merged);

  my_event_->Branch("mom_bestMCS_mu", &mom_bestMCS_mu);
  my_event_->Branch("mom_bestMCS_ll_mu", &mom_bestMCS_ll_mu);
  my_event_->Branch("mom_fwdMCS_mu", &mom_fwdMCS_mu);
  my_event_->Branch("mom_fwdMCS_ll_mu", &mom_fwdMCS_ll_mu);
  my_event_->Branch("mom_bwdMCS_mu", &mom_bwdMCS_mu);
  my_event_->Branch("mom_bwdMCS_ll_mu", &mom_bwdMCS_ll_mu);
  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
  my_event_->Branch("mom_Range_p", &mom_Range_p);
  my_event_->Branch("mom_Range_pi", &mom_Range_pi);
  my_event_->Branch("mom_Range_mu_noSCE", &mom_Range_mu_noSCE);
  my_event_->Branch("mom_Range_p_noSCE", &mom_Range_p_noSCE);
  my_event_->Branch("mom_Range_pi_noSCE", &mom_Range_pi_noSCE);
  my_event_->Branch("mom_range_PID_avg_noSCE", &mom_range_PID_avg_noSCE);
  my_event_->Branch("vtx_x", &vtx_x);
  my_event_->Branch("vtx_y", &vtx_y);
  my_event_->Branch("vtx_z", &vtx_z);
  my_event_->Branch("vtx_x_MCS", &vtx_x_MCS);
  my_event_->Branch("vtx_y_MCS", &vtx_y_MCS);
  my_event_->Branch("vtx_z_MCS", &vtx_z_MCS);
  my_event_->Branch("start_x", &start_x);
  my_event_->Branch("start_y", &start_y);
  my_event_->Branch("start_z", &start_z);
  my_event_->Branch("end_x", &end_x);
  my_event_->Branch("end_y", &end_y);
  my_event_->Branch("end_z", &end_z);
  my_event_->Branch("trk_phi", &trk_phi);
  my_event_->Branch("trk_theta", &trk_theta);
  my_event_->Branch("trk_costheta", &trk_costheta);
  my_event_->Branch("trk_phi_MCS", &trk_phi_MCS);
  my_event_->Branch("trk_theta_MCS", &trk_theta_MCS);
  my_event_->Branch("trk_costheta_MCS", &trk_costheta_MCS);
  my_event_->Branch("trk_length_pl0", &trk_length_pl0);
  my_event_->Branch("trk_length_pl1", &trk_length_pl1);
  my_event_->Branch("trk_length_pl2", &trk_length_pl2);
  my_event_->Branch("trk_length_avg", &trk_length_avg);
  my_event_->Branch("trk_length_noSCE", &trk_length_noSCE);
  my_event_->Branch("trk_ifcontained", &trk_ifcontained);
  my_event_->Branch("trk_OutOfTime", &trk_OutOfTime);
  my_event_->Branch("vtx_FV", &vtx_FV);
  my_event_->Branch("vtx_MCS_FV", &vtx_MCS_FV);

  my_event_->Branch("old_vtx_FV", &old_vtx_FV);
  my_event_->Branch("old_trk_ifcontained", &old_trk_ifcontained);

  my_event_->Branch("hits_dEdx_size_pl0", &hits_dEdx_size_pl0);
  my_event_->Branch("hits_dEdx_size_pl1", &hits_dEdx_size_pl1);
  my_event_->Branch("hits_dEdx_size_pl2", &hits_dEdx_size_pl2);

  my_event_->Branch("dEdx_pl0", &dEdx_pl0);
  my_event_->Branch("dEdx_pl1", &dEdx_pl1);
  my_event_->Branch("dEdx_pl2", &dEdx_pl2);
  my_event_->Branch("dQdx_pl0", &dQdx_pl0);
  my_event_->Branch("dQdx_pl1", &dQdx_pl1);
  my_event_->Branch("dQdx_pl2", &dQdx_pl2);
  my_event_->Branch("resRange_pl0", &resRange_pl0);
  my_event_->Branch("resRange_pl1", &resRange_pl1);
  my_event_->Branch("resRange_pl2", &resRange_pl2);
  my_event_->Branch("dEdx_pl0_start_half", &dEdx_pl0_start_half);
  my_event_->Branch("dEdx_pl1_start_half", &dEdx_pl1_start_half);
  my_event_->Branch("dEdx_pl2_start_half", &dEdx_pl2_start_half);
  my_event_->Branch("dEdx_pl0_end_half", &dEdx_pl0_end_half);
  my_event_->Branch("dEdx_pl1_end_half", &dEdx_pl1_end_half);
  my_event_->Branch("dEdx_pl2_end_half", &dEdx_pl2_end_half);
  my_event_->Branch("dEdx_pl0_start10", &dEdx_pl0_start10);
  my_event_->Branch("dEdx_pl1_start10", &dEdx_pl1_start10);
  my_event_->Branch("dEdx_pl2_start10", &dEdx_pl2_start10);
  my_event_->Branch("dEdx_pl0_end10", &dEdx_pl0_end10);
  my_event_->Branch("dEdx_pl1_end10", &dEdx_pl1_end10);
  my_event_->Branch("dEdx_pl2_end10", &dEdx_pl2_end10);
  my_event_->Branch("dEdx_pl0_start1020", &dEdx_pl0_start1020);
  my_event_->Branch("dEdx_pl1_start1020", &dEdx_pl1_start1020);
  my_event_->Branch("dEdx_pl2_start1020", &dEdx_pl2_start1020);
  my_event_->Branch("dEdx_pl0_end1020", &dEdx_pl0_end1020);
  my_event_->Branch("dEdx_pl1_end1020", &dEdx_pl1_end1020);
  my_event_->Branch("dEdx_pl2_end1020", &dEdx_pl2_end1020);
  
  my_event_->Branch("dEdx_pl2_1020_ratio", &dEdx_pl2_1020_ratio);
  
  my_event_->Branch("PID_Chi2Mu_pl0", &PID_Chi2Mu_pl0);
  my_event_->Branch("PID_Chi2Mu_pl1", &PID_Chi2Mu_pl1);
  my_event_->Branch("PID_Chi2Mu_pl2", &PID_Chi2Mu_pl2);
  my_event_->Branch("PID_Chi2Mu_3pl", &PID_Chi2Mu_3pl);
  my_event_->Branch("PID_Chi2P_pl0", &PID_Chi2P_pl0);
  my_event_->Branch("PID_Chi2P_pl1", &PID_Chi2P_pl1);
  my_event_->Branch("PID_Chi2P_pl2", &PID_Chi2P_pl2);
  my_event_->Branch("PID_Chi2P_3pl", &PID_Chi2P_3pl);
  my_event_->Branch("PID_Chi2Pi_pl0", &PID_Chi2Pi_pl0);
  my_event_->Branch("PID_Chi2Pi_pl1", &PID_Chi2Pi_pl1);
  my_event_->Branch("PID_Chi2Pi_pl2", &PID_Chi2Pi_pl2);
  my_event_->Branch("PID_Chi2Pi_3pl", &PID_Chi2Pi_3pl);
  my_event_->Branch("PID_Chi2K_pl0", &PID_Chi2K_pl0);
  my_event_->Branch("PID_Chi2K_pl1", &PID_Chi2K_pl1);
  my_event_->Branch("PID_Chi2K_pl2", &PID_Chi2K_pl2);
  my_event_->Branch("PID_Chi2K_3pl", &PID_Chi2K_3pl);

  my_event_->Branch("PID_Pdg_3pl", &PID_Pdg_3pl);
  my_event_->Branch("PID_Pdg_pl2", &PID_Pdg_pl2);
  my_event_->Branch("PID_Pdg_pl1", &PID_Pdg_pl1);
  my_event_->Branch("PID_Pdg_pl0", &PID_Pdg_pl0);
  my_event_->Branch("PID_avg_Chi2", &PID_avg_Chi2);
  my_event_->Branch("PID_pl2_Chi2", &PID_pl2_Chi2);
  my_event_->Branch("PID_pl1_Chi2", &PID_pl1_Chi2);
  my_event_->Branch("PID_pl0_Chi2", &PID_pl0_Chi2);

  my_event_->Branch("Pl0_for_PID", &Pl0_for_PID);
  my_event_->Branch("Pl1_for_PID", &Pl1_for_PID);
  my_event_->Branch("Pl2_for_PID", &Pl2_for_PID);
  my_event_->Branch("BestPlane_PID", &BestPlane_PID);

  my_event_->Branch("if_fwd_MCS", &if_fwd_MCS);
  my_event_->Branch("if_fwd_true", &if_fwd_true);
  my_event_->Branch("if_fwd_dEdx10", &if_fwd_dEdx10);
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

//int SingleMuon::Topology(int Nmuons, int Nelectrons, int Nprotons, int Npiplus, int Npiminus, int Npi0, int Nprotons_abTH, int Nprotons_blTH, int PDG, int CCNC,  bool ifcosmic, bool ifbeam, bool ifFV){
//  // In the current version, Proton Momentum threshold is set to be 300MeV, which should be reviewed soon
//  // 1. NuMuCC0pi0p in FV
//  if (Nmuons > 0 && Npi0 == 0 && Npiplus == 0 && Npiminus == 0 && Nprotons_abTH == 0 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 1;
//  // 2. NuMuCC0pi1p in FV
//  if (Nmuons > 0 && Npi0 == 0 && Npiplus == 0 && Npiminus == 0 && Nprotons_abTH == 1 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 2;
//  // 3. NuMuCC0pi2p in FV
//  if (Nmuons > 0 && Npi0 == 0 && Npiplus == 0 && Npiminus == 0 && Nprotons_abTH == 2 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 3;
//  // 4. NuMuCC0piNp in FV
//  if (Nmuons > 0 && Npi0 == 0 && Npiplus == 0 && Npiminus == 0 && Nprotons_abTH > 2 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 4; 
//  // 5. NuMuCC1pi+Xp in FV
//  if (Nmuons > 0 && Npi0 == 0 && Npiplus == 1 && Npiminus == 0 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 5;
//  // 6. NuMuCC1pi-Xp in FV
//  if (Nmuons > 0 && Npi0 == 0 && Npiplus == 0 && Npiminus == 1 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 6;
//  // 7. NuMuCC1pi0Xp in FV
//  if (Nmuons > 0 && Npi0 == 1 && Npiplus == 0 && Npiminus == 0 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 7;
//  // 8. NuMuCCNpiXp in FV
//  if (Nmuons > 0 && (Npi0 + Npiplus + Npiminus) > 0 && ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == 14) return 8;
//  // 9. Anti NuMu CC in FV
//  if (ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && PDG == -14) return 9;
//  // 10. Nue / Anti-Nue CC in FV
//  if (ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 0 && abs(PDG) == 12) return 10;
//  // 11. NC
//  if (ifcosmic == 0 && ifbeam == 1 && ifFV == 1 && CCNC = 1) return 11;
//  // 12. true neutrino vertex out of FV
//  if (ifcosmic == 0 && ifbeam == 1 && ifFV == 0) return 12;
//  // 13. cosmic
//  if (ifcosmic == 1) return 13;
//  // 14. other
//  else return 14;
//}

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

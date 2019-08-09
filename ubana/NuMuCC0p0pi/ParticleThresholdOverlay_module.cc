////////////////////////////////////////////////////////////////////////
// Class:       ParticleThreshold
// Plugin Type: analyzer (art v3_01_02)
// File:        ParticleThreshold_module.cc
//
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

class ParticleThreshold;


class ParticleThreshold : public art::EDAnalyzer {
public:
  explicit ParticleThreshold(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParticleThreshold(ParticleThreshold const&) = delete;
  ParticleThreshold(ParticleThreshold&&) = delete;
  ParticleThreshold& operator=(ParticleThreshold const&) = delete;
  ParticleThreshold& operator=(ParticleThreshold&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
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
  
  std::vector<double> All_true_PDG;//Track pdg
  std::vector<double> All_true_mom;//True momentum of muon track in the every event
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
  std::vector<double> All_true_trk_theta_yz;
  std::vector<double> All_true_trk_costheta_yz;
  std::vector<double> All_true_trk_theta_xz;
  std::vector<double> All_true_trk_costheta_xz;
  std::vector<double> All_true_trk_length; 

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

  std::vector<bool> if_cosmic; // Check if a track is cosmic or not by if it has an associated MCParticle

  std::vector<double> mom_bestMCS_mu;//MCS best momentum of muon track in the every event
  std::vector<double> mom_bestMCS_ll_mu;//Likelihood of MCS best momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_mu;//MCS forward momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_ll_mu;//Likelihood of MCS forward momentum of muon track in the every event
  std::vector<double> mom_bkwdMCS_mu;//MCS backward momentum of muon track in the every event
  std::vector<double> mom_bkwdMCS_ll_mu;//Likelihood of MCS backward momentum of muon track in the every event
  std::vector<double> mom_Range_mu;//Range momentum of muon track in the every event
  std::vector<double> mom_Range_p;//Range momentum of proton track in the every event
  std::vector<double> mom_Range_pi;//Range momentum of pion track in the every event
  std::vector<double> mom_range_PID_avg;//Range momentum of tracks based on their PID particle type using 3 pls
  std::vector<double> mom_range_PID_pl2;//Range momentum of tracks based on their PID particle type using pl 2
  std::vector<double> mom_range_truePDG;//Range momentum of tracks based on their true particle type
  std::vector<double> mom_Range_mu_noSCE;//Range momentum of muon track in the every event
  std::vector<double> mom_Range_p_noSCE;//Range momentum of proton track in the every event
  std::vector<double> mom_Range_pi_noSCE;//Range momentum of pion track in the every event
  std::vector<double> mom_range_PID_avg_noSCE;//Range momentum of tracks based on their PID particle type using 3 pls
  std::vector<double> mom_range_PID_pl2_noSCE;//Range momentum of tracks based on their PID particle type using pl 2
  std::vector<double> mom_range_truePDG_noSCE;//Range momentum of tracks based on their true particle type

  std::vector<double> vtx_x;//Reconstructed vtx x in the every event
  std::vector<double> vtx_y;//Reconstructed vtx y in the every event
  std::vector<double> vtx_z;//Reconstructed vtx z in the every event
  std::vector<double> start_x;//Reconstructed start x in the every event
  std::vector<double> start_y;//Reconstructed start y in the every event
  std::vector<double> start_z;//Reconstructed start z in the every event
  std::vector<double> end_x;//Reconstructed end x in the every event
  std::vector<double> end_y;//Reconstructed end y in the every event
  std::vector<double> end_z;//Reconstructed end z in the every event
  std::vector<double> trk_phi;//Reconstructed track phi in the every event
  std::vector<double> trk_theta;//Reconstructed track theta in the every event
  std::vector<double> trk_costheta;//Reconstructed track cos(theta) in the every event
  std::vector<double> trk_length_pl0;//Range momentum of muon track in the every event
  std::vector<double> trk_length_pl1;//Range momentum of muon track in the every event
  std::vector<double> trk_length_pl2;//Range momentum of muon track in the every event
  std::vector<double> trk_length_avg;//Range momentum of muon track in the every event
  std::vector<double> trk_length_noSCE;//Range momentum of muon track in the every event no spatial correction
  std::vector<bool> trk_ifcontained;//to check if the track is contained or not
  std::vector<bool> vtx_FV;//to check if the vertex is in FV or not
  std::vector<double> trk_end_theta_yz;
  std::vector<double> trk_end_costheta_yz;
  std::vector<double> trk_end_theta_xz;
  std::vector<double> trk_end_costheta_xz;
  std::vector<double> trk_theta_yz;
  std::vector<double> trk_costheta_yz;
  std::vector<double> trk_theta_xz;
  std::vector<double> trk_costheta_xz;
 
  int Ntrack; // number of tracks in this event

  std::vector<int> hits_dEdx_size_pl0;
  std::vector<int> hits_dEdx_size_pl1;
  std::vector<int> hits_dEdx_size_pl2;

  std::vector<std::vector<float>> dEdx_pl0; // dE/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<std::vector<float>> dEdx_pl1; // dE/dx of the selected (muon) track from plane 1
  std::vector<std::vector<float>> dEdx_pl2; // dE/dx of the selected (muon) track from plane 2 (collection)
  std::vector<std::vector<float>> dQdx_pl0; // dQ/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<std::vector<float>> dQdx_pl1; // dQ/dx of the selected (muon) track from plane 1
  std::vector<std::vector<float>> dQdx_pl2; // dQ/dx of the selected (muon) track from plane 2 (collection)
  std::vector<std::vector<float>> resRange_pl0; // range from a hit to the end of the selected track end
  std::vector<std::vector<float>> resRange_pl1; // range from a hit to the end of the selected track end
  std::vector<std::vector<float>> resRange_pl2; // range from a hit to the end of the selected track end
  std::vector<float> dEdx_pl0_start_half; // average dEdx of start half hits of pl 0
  std::vector<float> dEdx_pl1_start_half; // average dEdx of start half hits of pl 0
  std::vector<float> dEdx_pl2_start_half; // average dEdx of start half hits of pl 0
  std::vector<float> dEdx_pl0_end_half; // average dEdx of end half hits of pl 0
  std::vector<float> dEdx_pl1_end_half; // average dEdx of end half hits of pl 0
  std::vector<float> dEdx_pl2_end_half; // average dEdx of end half hits of pl 0
  std::vector<float> dEdx_pl0_start10; // average dEdx of first 10 hit of pl 0
  std::vector<float> dEdx_pl1_start10; // average dEdx of first 10 hit of pl 0
  std::vector<float> dEdx_pl2_start10; // average dEdx of first 10 hit of pl 0
  std::vector<float> dEdx_pl0_end10; // average dEdx of end 10 hit of pl 0
  std::vector<float> dEdx_pl1_end10; // average dEdx of end 10 hit of pl 0
  std::vector<float> dEdx_pl2_end10; // average dEdx of end 10 hit of pl 0
  std::vector<float> dEdx_pl0_start1020; // average dEdx of first 10 hit of pl 0
  std::vector<float> dEdx_pl1_start1020; // average dEdx of first 10 hit of pl 0
  std::vector<float> dEdx_pl2_start1020; // average dEdx of first 10 hit of pl 0
  std::vector<float> dEdx_pl0_end1020; // average dEdx of end 10 hit of pl 0
  std::vector<float> dEdx_pl1_end1020; // average dEdx of end 10 hit of pl 0
  std::vector<float> dEdx_pl2_end1020; // average dEdx of end 10 hit of pl 0

  std::vector<double> PID_Chi2Mu_pl0; // Chi2 of muon assumption of plane 0 in PID
  std::vector<double> PID_Chi2Mu_pl1; // Chi2 of muon assumption of plane 1 in PID
  std::vector<double> PID_Chi2Mu_pl2; // Chi2 of muon assumption of plane 2 in PID
  std::vector<double> PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID

  std::vector<double> PID_Chi2P_pl0; // Chi2 of proton assumption of plane 0 in PID
  std::vector<double> PID_Chi2P_pl1; // Chi2 of proton assumption of plane 1 in PID
  std::vector<double> PID_Chi2P_pl2; // Chi2 of proton assumption of plane 2 in PID
  std::vector<double> PID_Chi2P_3pl; // Chi2 of muon assumption of 3 planes in PID

  std::vector<double> PID_Chi2Pi_pl0; // Chi2 of pion assumption of plane 0 in PID
  std::vector<double> PID_Chi2Pi_pl1; // Chi2 of pion assumption of plane 1 in PID
  std::vector<double> PID_Chi2Pi_pl2; // Chi2 of pion assumption of plane 2 in PID
  std::vector<double> PID_Chi2Pi_3pl; // Chi2 of muon assumption of 3 planes in PID

  std::vector<double> PID_Chi2K_pl0; // Chi2 of kaon assumption of plane 0 in PID
  std::vector<double> PID_Chi2K_pl1; // Chi2 of kaon assumption of plane 1 in PID
  std::vector<double> PID_Chi2K_pl2; // Chi2 of kaon assumption of plane 2 in PID
  std::vector<double> PID_Chi2K_3pl; // Chi2 of muon assumption of 3 planes in PID

  std::vector<int> PID_Pdg_allPlane; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
  std::vector<int> PID_Pdg_pl2; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
  std::vector<int> PID_Pdg_pl1;
  std::vector<int> PID_Pdg_pl0;
  std::vector<double> PID_avg_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions
  std::vector<double> PID_pl2_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions
  std::vector<double> PID_pl1_Chi2;
  std::vector<double> PID_pl0_Chi2;
  
  std::vector<bool> Pl0_for_PID;
  std::vector<bool> Pl1_for_PID;
  std::vector<bool> Pl2_for_PID;
  std::vector<int> BestPlane_PID;

  std::vector<bool> if_fwd_MCS; // If using forward MCS direction judge
  std::vector<bool> if_fwd_true; // If fwd by the reco true vertex distance
  std::vector<bool> if_fwd_dEdx10; // If fwd by the reco dEdx 10 hits (should use for contained)
  std::vector<bool> if_fwd_dEdx1020; // If fwd by the reco dEdx 10 hits (should use for contained)
  std::vector<bool> if_fwd_dEdxhalf; // If fwd by the reco dEdx half of the hits (should use for contained)

  bool                                IsMC;
  bool                                IfAll;
  std::string                         m_generatorLabel;
  std::string                         m_geantLabel;
  std::string                         m_pandoraLabel;
  std::string                         m_hitProducerLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_showerProducerLabel;
  std::string                         m_MCSmuProducerLabel;
  std::string                         m_calorimetryProducerLabel;
  //std::string                         m_SCEcorr_MCSmuProducerLabel;
  std::string                         Hits_TrackAssLabel;
  std::string                         PID_TrackAssLabel;

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


ParticleThreshold::ParticleThreshold(fhicl::ParameterSet const& pset)
  : 
  EDAnalyzer{pset},
  IsMC(pset.get<bool>("IsMC")),
  IfAll(pset.get<bool>("IfAll")),
  m_generatorLabel(pset.get<std::string>("GeneratorLabel")),
  m_geantLabel(pset.get<std::string>("GeantLabel")),
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_hitProducerLabel(pset.get<std::string>("HitProducerLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel")),
  m_MCSmuProducerLabel(pset.get<std::string>("MCSmuProducerLabel")),
  m_calorimetryProducerLabel(pset.get<std::string>("calorimetryProducerLabel")),
  //m_SCEcorr_MCSmuProducerLabel(pset.get<std::string>("SCEcorr_MCSmuProducerLabel")),
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

void ParticleThreshold::analyze(art::Event const& evt)
{

  //// Get necessary handles
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;

  if(IsMC){
    // MC Truth
    art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
    evt.getByLabel(m_generatorLabel, Handle_MCTruth);
    art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

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

  //pfp track association
  art::FindManyP<recob::Track> pfpToTrackAsso(Handle_pfParticle, evt, m_trackProducerLabel);

  //pfp shower association
  art::FindManyP<recob::Shower> pfpToShowerAsso(Handle_pfParticle, evt, m_showerProducerLabel);
  
  //track calorimetry association 
  art::FindManyP<anab::Calorimetry> trackToCalAsso(Handle_TPCtrack, evt, m_calorimetryProducerLabel);

  //PID
  art::FindManyP<anab::ParticleID> PIDTotrackAsso(Handle_TPCtrack,evt,PID_TrackAssLabel);

  if(IfAll){
    //------Loop over all MCParticles
    for(int i_mcp = 0; i_mcp < (int) MCParticleCollection.size(); i_mcp++){
      auto MCP = MCParticleCollection[i_mcp];
      auto AllTrueTrackPos = MCP->EndPosition() - MCP->Position();
      All_true_PDG.push_back(MCP->PdgCode());
      All_true_mom.push_back(MCP->P());
      All_true_start_x.push_back(MCP->Position().X());
      All_true_start_y.push_back(MCP->Position().Y());
      All_true_start_z.push_back(MCP->Position().Z());
      All_true_end_x.push_back(MCP->EndPosition().X());
      All_true_end_y.push_back(MCP->EndPosition().Y());
      All_true_end_z.push_back(MCP->EndPosition().Z());

      TVector3 All_true_start(All_true_start_x.back(), All_true_start_y.back(), All_true_start_z.back());
      TVector3 All_true_end(All_true_end_x.back(), All_true_end_y.back(), All_true_end_z.back());
      All_true_trk_ifcontained.push_back(_fiducial_volume.InFV(All_true_start, All_true_end));
      All_true_vtxFV.push_back(_fiducial_volume.InFV(All_true_start));

      All_true_trk_phi.push_back(AllTrueTrackPos.Phi());
      All_true_trk_theta.push_back(AllTrueTrackPos.Theta());
      All_true_trk_costheta.push_back(cos(AllTrueTrackPos.Theta()));
      All_true_trk_theta_yz.push_back(std::atan2(AllTrueTrackPos.Y(), AllTrueTrackPos.Z()));
      All_true_trk_costheta_yz.push_back(cos(All_true_trk_theta_yz.back()));
      All_true_trk_theta_xz.push_back(std::atan2(AllTrueTrackPos.X(), AllTrueTrackPos.Z()));
      All_true_trk_costheta_xz.push_back(cos(All_true_trk_theta_xz.back()));
      All_true_trk_length.push_back(sqrt(AllTrueTrackPos.X()*AllTrueTrackPos.X() + AllTrueTrackPos.Y()*AllTrueTrackPos.Y() + AllTrueTrackPos.Z()*AllTrueTrackPos.Z())); // An estimation of true track length
    }
  }
  //------Loop over all tracks
  Ntrack = AllTrackCollection.size();
  for(int trk_id = 0; trk_id < Ntrack; trk_id++){

    //---Fill in True MCParticle information
    std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[trk_id].key());
    BackTrackerTruthMatch backtrackertruthmatch;
    backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
    auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
    if(!MCparticle){
      if_cosmic.push_back(true);
      std::cout<<"MC particle does not exist!"<<std::endl;
    }
    else{
      if_cosmic.push_back(false);
      auto TrueTrackPos = MCparticle->EndPosition() - MCparticle->Position();
      true_mom.push_back(MCparticle->P());
      true_start_x.push_back(MCparticle->Position().X());
      true_start_y.push_back(MCparticle->Position().Y());
      true_start_z.push_back(MCparticle->Position().Z());
      true_end_x.push_back(MCparticle->EndPosition().X());
      true_end_y.push_back(MCparticle->EndPosition().Y());
      true_end_z.push_back(MCparticle->EndPosition().Z());
      TVector3 true_start(true_start_x.back(), true_start_y.back(), true_start_z.back());
      TVector3 true_end(true_end_x.back(), true_end_y.back(), true_end_z.back());
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

      //---Fill in Reco information
      // Correct the position of the track ends
      TVector3 Trk_start = AllTrackCollection[trk_id]->Vertex<TVector3>();
      auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
      TVector3 Trk_start_SCEcorr;
      Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
      Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
      Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

      TVector3 Trk_end = AllTrackCollection[trk_id]->End<TVector3>();
      auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
      TVector3 Trk_end_SCEcorr;
      Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X());
      Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
      Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

      trk_ifcontained.push_back(_fiducial_volume.InFV(Trk_start_SCEcorr, Trk_end_SCEcorr));

      vtx_FV.push_back(_fiducial_volume.InFV(Trk_start_SCEcorr));

      mom_bestMCS_mu.push_back(mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestMomentum());
      mom_bestMCS_ll_mu.push_back(mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestLogLikelihood());

      mom_fwdMCS_mu.push_back(mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdMomentum());
      mom_fwdMCS_ll_mu.push_back(mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdLogLikelihood());

      mom_bkwdMCS_mu.push_back(mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdMomentum());
      mom_bkwdMCS_ll_mu.push_back(mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdLogLikelihood());
    
      // track length with no SCE 
      trk_length_noSCE.push_back(AllTrackCollection[trk_id]->Length());
      auto assoCal = trackToCalAsso.at(AllTrackCollection[trk_id].key());
      trk_length_pl0.push_back(assoCal[0]->Range());  //pandoracali has spatial correction
      trk_length_pl1.push_back(assoCal[1]->Range());  //pandoracali has spatial correction
      trk_length_pl2.push_back(assoCal[2]->Range());  //pandoracali has spatial correction
      double Trk_length_avg = 0;
      int valid_pl = 0;
      for (int i_pl = 0; i_pl < (int) assoCal.size(); i_pl++){
        if(assoCal[i_pl]->Range() > 0){
          Trk_length_avg += assoCal[i_pl]->Range();
          valid_pl++;
        }
      }
      trk_length_avg.push_back(Trk_length_avg / valid_pl);
      
      vtx_x.push_back(Trk_start_SCEcorr.X());
      vtx_y.push_back(Trk_start_SCEcorr.Y());
      vtx_z.push_back(Trk_start_SCEcorr.Z());
      start_x.push_back(Trk_start_SCEcorr.X());
      start_y.push_back(Trk_start_SCEcorr.Y());
      start_z.push_back(Trk_start_SCEcorr.Z());
      end_x.push_back(Trk_end_SCEcorr.X());
      end_y.push_back(Trk_end_SCEcorr.Y());
      end_z.push_back(Trk_end_SCEcorr.Z());

      trk_phi.push_back(AllTrackCollection[trk_id]->Phi());
      trk_theta.push_back(AllTrackCollection[trk_id]->Theta());
      trk_costheta.push_back(cos(AllTrackCollection[trk_id]->Theta()));

      if(trk_length_pl2.back() < 0){
         mom_Range_mu.push_back(-999);
         mom_Range_mu_noSCE.push_back(-999);
         mom_Range_p.push_back(-999);
         mom_Range_p_noSCE.push_back(-999);
         mom_Range_pi.push_back(-999);
         mom_Range_pi_noSCE.push_back(-999);
      }
      else{
        mom_Range_mu.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length_pl2.back(), 13));
        mom_Range_mu_noSCE.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length_noSCE.back(), 13));
        mom_Range_p.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length_pl2.back(), 2212));
        mom_Range_p_noSCE.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length_noSCE.back(), 2212));
        mom_Range_pi.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length_pl2.back(), 211));
        mom_Range_pi_noSCE.push_back(_trk_mom_calculator.GetTrackMomentum(trk_length_noSCE.back(), 211));
      }
 
      //Get Calorimety of the track
      if(assoCal.size()!=3){
        throw cet::exception("[Numu0pi0p]") << "Where are the three planes for the calorimetry!" << std::endl;
      }
      // induction 0 = 0, induction 1 = 1, collection = 2 (Check if the plane ID is correct)
      // assoCal[id_pl]->PlaneID().Plane == 2 (collection)
      // The vector is ordered by residual range from small to big (track end to track start)
      dEdx_pl0.push_back(assoCal[0]->dEdx());
      dQdx_pl0.push_back(assoCal[0]->dQdx());
      resRange_pl0.push_back(assoCal[0]->ResidualRange());

      dEdx_pl1.push_back(assoCal[1]->dEdx());
      dQdx_pl1.push_back(assoCal[1]->dQdx());
      resRange_pl1.push_back(assoCal[1]->ResidualRange());

      dEdx_pl2.push_back(assoCal[2]->dEdx());
      dQdx_pl2.push_back(assoCal[2]->dQdx());
      resRange_pl2.push_back(assoCal[2]->ResidualRange());

      int Nhits_pl0 = dEdx_pl0.back().size();
      int Nhits_pl1 = dEdx_pl1.back().size();
      int Nhits_pl2 = dEdx_pl2.back().size();

      int half_size_pl0 = Nhits_pl0 / 2;
      int half_size_pl1 = Nhits_pl1 / 2;
      int half_size_pl2 = Nhits_pl2 / 2;

      hits_dEdx_size_pl0.push_back(Nhits_pl0);
      hits_dEdx_size_pl1.push_back(Nhits_pl1);
      hits_dEdx_size_pl2.push_back(Nhits_pl2);

      dEdx_pl0_start_half.push_back(std::accumulate(dEdx_pl0.back().end() - half_size_pl0, dEdx_pl0.back().end(), 0.) / half_size_pl0);
      dEdx_pl0_end_half.push_back(std::accumulate(dEdx_pl0.back().begin(), dEdx_pl0.back().begin() + half_size_pl0, 0. ) / half_size_pl0);
      dEdx_pl1_start_half.push_back(std::accumulate(dEdx_pl1.back().end() - half_size_pl1, dEdx_pl1.back().end(), 0.) / half_size_pl1);
      dEdx_pl1_end_half.push_back(std::accumulate(dEdx_pl1.back().begin(), dEdx_pl1.back().begin() + half_size_pl1, 0. ) / half_size_pl1);
      dEdx_pl2_start_half.push_back(std::accumulate(dEdx_pl2.back().end() - half_size_pl2, dEdx_pl2.back().end(), 0.) / half_size_pl2);
      dEdx_pl2_end_half.push_back(std::accumulate(dEdx_pl2.back().begin(), dEdx_pl2.back().begin() + half_size_pl2, 0. ) / half_size_pl2);

      // dEdx_10
      if (dEdx_pl0.back().size()<10) {
        dEdx_pl0_start10.push_back(dEdx_pl0_start_half.back());
        dEdx_pl0_end10.push_back(dEdx_pl0_end_half.back());
      }
      else{
        dEdx_pl0_start10.push_back(std::accumulate(dEdx_pl0.back().end() - 10, dEdx_pl0.back().end(), 0.) / 10.);
        dEdx_pl0_end10.push_back(std::accumulate(dEdx_pl0.back().begin(), dEdx_pl0.back().begin() + 10, 0.) / 10.);
      }
      if (dEdx_pl1.back().size()<10) {
        dEdx_pl1_start10.push_back(dEdx_pl1_start_half.back());
        dEdx_pl1_end10.push_back(dEdx_pl1_end_half.back());
      }
      else{
        dEdx_pl1_start10.push_back(std::accumulate(dEdx_pl1.back().end() - 10, dEdx_pl1.back().end(), 0.) / 10.);
        dEdx_pl1_end10.push_back(std::accumulate(dEdx_pl1.back().begin(), dEdx_pl1.back().begin() + 10, 0.) / 10.);
      }
      if (dEdx_pl2.back().size()<10) {
        dEdx_pl2_start10.push_back(dEdx_pl2_start_half.back());
        dEdx_pl2_end10.push_back(dEdx_pl2_end_half.back());
      }
      else{
        dEdx_pl2_start10.push_back(std::accumulate(dEdx_pl2.back().end() - 10, dEdx_pl2.back().end(), 0.) / 10.);
        dEdx_pl2_end10.push_back(std::accumulate(dEdx_pl2.back().begin(), dEdx_pl2.back().begin() + 10, 0.) / 10.);
      }
      // dEdx_1020
      if (dEdx_pl0.back().size()<30) {
        dEdx_pl0_start1020.push_back(dEdx_pl0_start_half.back());
        dEdx_pl0_end1020.push_back(dEdx_pl0_end_half.back());
      }
      else{
        dEdx_pl0_start1020.push_back(std::accumulate(dEdx_pl0.back().end() - 20, dEdx_pl0.back().end() - 10, 0.) / 10.);
        dEdx_pl0_end1020.push_back(std::accumulate(dEdx_pl0.back().begin() + 10, dEdx_pl0.back().begin() + 20, 0.) / 10.);
      }
      if (dEdx_pl1.back().size()<30) {
        dEdx_pl1_start1020.push_back(dEdx_pl1_start_half.back());
        dEdx_pl1_end1020.push_back(dEdx_pl1_end_half.back());
      }
      else{
        dEdx_pl1_start1020.push_back(std::accumulate(dEdx_pl1.back().end() - 20, dEdx_pl1.back().end() - 10, 0.) / 10.);
        dEdx_pl1_end1020.push_back(std::accumulate(dEdx_pl1.back().begin() + 10, dEdx_pl1.back().begin() + 20, 0.) / 10.);
      }
      if (dEdx_pl2.back().size()<30) {
        dEdx_pl2_start1020.push_back(dEdx_pl2_start_half.back());
        dEdx_pl2_end1020.push_back(dEdx_pl2_end_half.back());
      }
      else{
        dEdx_pl2_start1020.push_back(std::accumulate(dEdx_pl2.back().end() - 20, dEdx_pl2.back().end() - 10, 0.) / 10.);
        dEdx_pl2_end1020.push_back(std::accumulate(dEdx_pl2.back().begin() + 10, dEdx_pl2.back().begin() + 20, 0.) / 10.);
      }


      //--- Gain PID info of the track
      if(!PIDTotrackAsso.isValid()){
        throw cet::exception("[Numu0pi0p]") << "No matched PID - track information!" << std::endl;
      }
      // Get projected angle wrt to the wires (docdb 23008)
      TVector3 End_Dir = AllTrackCollection[trk_id]->EndDirection<TVector3>();
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
      auto trkPID = PIDTotrackAsso.at(AllTrackCollection[trk_id].key());
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
      PID_Chi2Mu_pl0.push_back(PIDChi2_mu[0]);
      PID_Chi2Mu_pl1.push_back(PIDChi2_mu[1]);
      PID_Chi2Mu_pl2.push_back(PIDChi2_mu[2]);
      int ww0 = w0; int ww1 = w1; int ww2 = w2; // copy from the origin 
      if (PIDChi2_mu[0] < 0) ww0 = 0;
      if (PIDChi2_mu[1] < 0) ww1 = 0;
      if (PIDChi2_mu[2] < 0) ww2 = 0;
      if (ww0 + ww1 + ww2 == 0) PID_Chi2Mu_3pl.push_back(-999);
      else PID_Chi2Mu_3pl.push_back((ww0 * PIDChi2_mu[0] + ww1 * PIDChi2_mu[1] + ww2 * PIDChi2_mu[2]) / (ww0 + ww1 + ww2));

      PID_Chi2P_pl0.push_back(PIDChi2_p[0]);
      PID_Chi2P_pl1.push_back(PIDChi2_p[1]);
      PID_Chi2P_pl2.push_back(PIDChi2_p[2]);
      ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
      if (PIDChi2_p[0] < 0) ww0 = 0;
      if (PIDChi2_p[1] < 0) ww1 = 0;
      if (PIDChi2_p[2] < 0) ww2 = 0;
      if (ww0 + ww1 + ww2 == 0) PID_Chi2P_3pl.push_back(-999);
      else PID_Chi2P_3pl.push_back((ww0 * PIDChi2_p[0] + ww1 * PIDChi2_p[1] + ww2 * PIDChi2_p[2]) / (ww0 + ww1 + ww2));

      PID_Chi2Pi_pl0.push_back(PIDChi2_pi[0]);
      PID_Chi2Pi_pl1.push_back(PIDChi2_pi[1]);
      PID_Chi2Pi_pl2.push_back(PIDChi2_pi[2]);
      ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
      if (PIDChi2_pi[0] < 0) ww0 = 0;
      if (PIDChi2_pi[1] < 0) ww1 = 0;
      if (PIDChi2_pi[2] < 0) ww2 = 0;
      if (ww0 + ww1 + ww2 == 0) PID_Chi2Pi_3pl.push_back(-999);
      else PID_Chi2Pi_3pl.push_back((ww0 * PIDChi2_pi[0] + ww1 * PIDChi2_pi[1] + ww2 * PIDChi2_pi[2]) / (ww0 + ww1 + ww2));

      PID_Chi2K_pl0.push_back(PIDChi2_K[0]);
      PID_Chi2K_pl1.push_back(PIDChi2_K[1]);
      PID_Chi2K_pl2.push_back(PIDChi2_K[2]);
      ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
      if (PIDChi2_K[0] < 0) ww0 = 0;
      if (PIDChi2_K[1] < 0) ww1 = 0;
      if (PIDChi2_K[2] < 0) ww2 = 0;
      if (ww0 + ww1 + ww2 == 0) PID_Chi2K_3pl.push_back(-999);
      else PID_Chi2K_3pl.push_back((ww0 * PIDChi2_K[0] + ww1 * PIDChi2_K[1] + ww2 * PIDChi2_K[2]) / (ww0 + ww1 + ww2));

      std::vector<int> Nhit_3pl = {Nhits_pl0, Nhits_pl1, Nhits_pl2};//at similar condition 2> 0 > 1
      int bestpl = std::max_element(Nhit_3pl.begin(), Nhit_3pl.end()) - Nhit_3pl.begin();
      if (bestpl == 0){
        if (Nhits_pl0 == Nhits_pl2){
          bestpl = 2;
        }
        else{
          Pl0_for_PID.push_back(true);
          Pl1_for_PID.push_back(false);
          Pl2_for_PID.push_back(false);
        }
      }
      if (bestpl == 1){
        if (Nhits_pl1 == Nhits_pl2){
          bestpl = 2;
        }
        else{
          Pl0_for_PID.push_back(false);
          Pl1_for_PID.push_back(true);
          Pl2_for_PID.push_back(false);
        }
      }
      if (bestpl == 2){
        Pl0_for_PID.push_back(false);
        Pl1_for_PID.push_back(false);
        Pl2_for_PID.push_back(true);
      }
      BestPlane_PID.push_back(bestpl);

      //-- Get minimum Chi2 and there corresponding particle type
      std::vector<double> PIDChi2_avg;// It follows the order of muon, proton, pion, kaon
      if (PID_Chi2Mu_3pl.back() < 0) PIDChi2_avg.push_back(9999);
      else PIDChi2_avg.push_back(PID_Chi2Mu_3pl.back());
      if (PID_Chi2P_3pl.back() < 0) PID_Chi2Mu_3pl.push_back(9999);
      else PIDChi2_avg.push_back(PID_Chi2P_3pl.back());
      if (PID_Chi2Pi_3pl.back() < 0) PIDChi2_avg.push_back(9999);
      else PIDChi2_avg.push_back(PID_Chi2Pi_3pl.back());
      if (PID_Chi2K_3pl.back() < 0) PIDChi2_avg.push_back(9999);
      else PIDChi2_avg.push_back(PID_Chi2K_3pl.back());
      
      PID_avg_Chi2.push_back(*std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end()));
      int ID_PID = std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end()) - PIDChi2_avg.begin();
      if (ID_PID == 0) {
        PID_Pdg_allPlane.push_back(13);
      }
      if (ID_PID == 1) {
        PID_Pdg_allPlane.push_back(2212);
      }
      if (ID_PID == 2) {
        PID_Pdg_allPlane.push_back(211);
      }
      if (ID_PID == 3) {
        PID_Pdg_allPlane.push_back(321); 
      }
      double mom_range_PID_avg_value = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2.back(), PID_Pdg_allPlane.back());
      mom_range_PID_avg.push_back(mom_range_PID_avg_value);
      double mom_range_PID_avg_noSCE_value = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE.back(), PID_Pdg_allPlane.back());
      mom_range_PID_avg_noSCE.push_back(mom_range_PID_avg_noSCE_value);

      //-- Use plane 2 only
      std::vector<double> PIDChi2_pl2;// It follows the order of muon, proton, pion, kaon
      if (PIDChi2_mu[2] < 0) PIDChi2_pl2.push_back(9999);
      else PIDChi2_pl2.push_back(PIDChi2_mu[2]);
      if (PIDChi2_p[2] < 0) PIDChi2_pl2.push_back(9999);
      else PIDChi2_pl2.push_back(PIDChi2_p[2]);
      if (PIDChi2_pi[2] < 0) PIDChi2_pl2.push_back(9999);
      else PIDChi2_pl2.push_back(PIDChi2_pi[2]);
      if (PIDChi2_K[2] < 0) PIDChi2_pl2.push_back(9999);
      else PIDChi2_pl2.push_back(PIDChi2_K[2]);

      PID_pl2_Chi2.push_back(*std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end()));
      int ID_PID_pl2 = std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end()) - PIDChi2_pl2.begin();
      if (ID_PID_pl2 == 0) {
        PID_Pdg_pl2.push_back(13);
      }
      if (ID_PID_pl2 == 1) {
        PID_Pdg_pl2.push_back(2212);
      }
      if (ID_PID_pl2 == 2) {
        PID_Pdg_pl2.push_back(211);
      }
      if (ID_PID_pl2 == 3) {
        PID_Pdg_pl2.push_back(321);
      }
      //-- Use plane 1 only
      std::vector<double> PIDChi2_pl1;// It follows the order of muon, proton, pion, kaon
      if (PIDChi2_mu[1] < 0) PIDChi2_pl1.push_back(9999);
      else PIDChi2_pl1.push_back(PIDChi2_mu[1]);
      if (PIDChi2_p[1] < 0) PIDChi2_pl1.push_back(9999);
      else PIDChi2_pl1.push_back(PIDChi2_p[1]);
      if (PIDChi2_pi[1] < 0) PIDChi2_pl1.push_back(9999);
      else PIDChi2_pl1.push_back(PIDChi2_pi[1]);
      if (PIDChi2_K[1] < 0) PIDChi2_pl1.push_back(9999);
      else PIDChi2_pl1.push_back(PIDChi2_K[1]);

      PID_pl1_Chi2.push_back(*std::min_element(PIDChi2_pl1.begin(), PIDChi2_pl1.end()));
      int ID_PID_pl1 = std::min_element(PIDChi2_pl1.begin(), PIDChi2_pl1.end()) - PIDChi2_pl1.begin();
      if (ID_PID_pl1 == 0) {
        PID_Pdg_pl1.push_back(13);
      }
      if (ID_PID_pl1 == 1) {
        PID_Pdg_pl1.push_back(2212);
      }
      if (ID_PID_pl1 == 2) {
        PID_Pdg_pl1.push_back(211);
      }
      if (ID_PID_pl1 == 3) {
        PID_Pdg_pl1.push_back(321);
      }
      //-- Use plane 0 only
      std::vector<double> PIDChi2_pl0;// It follows the order of muon, proton, pion, kaon
      if (PIDChi2_mu[0] < 0) PIDChi2_pl0.push_back(9999);
      else PIDChi2_pl0.push_back(PIDChi2_mu[0]);
      if (PIDChi2_p[0] < 0) PIDChi2_pl0.push_back(9999);
      else PIDChi2_pl0.push_back(PIDChi2_p[0]);
      if (PIDChi2_pi[0] < 0) PIDChi2_pl0.push_back(9999);
      else PIDChi2_pl0.push_back(PIDChi2_pi[0]);
      if (PIDChi2_K[0] < 0) PIDChi2_pl0.push_back(9999);
      else PIDChi2_pl0.push_back(PIDChi2_K[0]);

      PID_pl0_Chi2.push_back(*std::min_element(PIDChi2_pl0.begin(), PIDChi2_pl0.end()));
      int ID_PID_pl0 = std::min_element(PIDChi2_pl0.begin(), PIDChi2_pl0.end()) - PIDChi2_pl0.begin();
      if (ID_PID_pl0 == 0) {
        PID_Pdg_pl0.push_back(13);
      }
      if (ID_PID_pl0 == 1) {
        PID_Pdg_pl0.push_back(2212);
      }
      if (ID_PID_pl0 == 2) {
        PID_Pdg_pl0.push_back(211);
      }
      if (ID_PID_pl0 == 3) {
        PID_Pdg_pl0.push_back(321);
      }

      double mom_range_PID_pl2_value = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2.back(), PID_Pdg_pl2.back());
      mom_range_PID_pl2.push_back(mom_range_PID_pl2_value);
      double mom_range_PID_pl2_noSCE_value = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE.back(), PID_Pdg_pl2.back());
      mom_range_PID_pl2_noSCE.push_back(mom_range_PID_pl2_noSCE_value);

      // Momentum by their true PDG
      double mom_range_truePDG_value = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2.back(), true_trk_PDG.back());
      mom_range_truePDG.push_back(mom_range_truePDG_value);
      double mom_range_truePDG_noSCE_value = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE.back(), true_trk_PDG.back());
      mom_range_truePDG_noSCE.push_back(mom_range_truePDG_noSCE_value);

      // Check the directional info of the track by MCS
      if (mom_bestMCS_ll_mu.back() == mom_fwdMCS_ll_mu.back()) if_fwd_MCS.push_back(true);
      else if_fwd_MCS.push_back(false);

      // Check the directional info of the track by true reco vertex distance
      TVector3 vtx(true_start_x.back(), true_start_y.back(), true_start_z.back());
      TVector3 D_start = Trk_start_SCEcorr - vtx;
      TVector3 D_end = Trk_end_SCEcorr - vtx;
      if (D_start.Mag() < D_end.Mag()) if_fwd_true.push_back(true);
      else if_fwd_true.push_back(false);
 
      // Check the direction info of the track by dEdx
      if (dEdx_pl2_start_half.back() < dEdx_pl2_end_half.back()) if_fwd_dEdxhalf.push_back(true);
      else if_fwd_dEdxhalf.push_back(false);
      if (dEdx_pl2_start10.back() < dEdx_pl2_end10.back()) if_fwd_dEdx10.push_back(true);
      else if_fwd_dEdx10.push_back(false);
      if (dEdx_pl2_start1020.back() < dEdx_pl2_end1020.back()) if_fwd_dEdx1020.push_back(true);
      else if_fwd_dEdx1020.push_back(false);
     
    } // MCParticle
  } // The end of track loops

  my_event_->Fill();
  if(IfAll){
    All_true_PDG.clear();
    All_true_mom.clear();
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
    All_true_trk_theta_yz.clear();
    All_true_trk_costheta_yz.clear();
    All_true_trk_theta_xz.clear();
    All_true_trk_costheta_xz.clear();
    All_true_trk_length.clear();
  }
 
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

  if_cosmic.clear();

  mom_bestMCS_mu.clear();
  mom_bestMCS_ll_mu.clear();
  mom_fwdMCS_mu.clear();
  mom_fwdMCS_ll_mu.clear();
  mom_bkwdMCS_mu.clear();
  mom_bkwdMCS_ll_mu.clear();
  mom_Range_mu.clear();
  mom_Range_p.clear();
  mom_Range_pi.clear();
  mom_range_PID_avg.clear();
  mom_range_PID_pl2.clear();
  mom_range_truePDG.clear();
  mom_Range_mu_noSCE.clear();
  mom_Range_p_noSCE.clear();
  mom_Range_pi_noSCE.clear();
  mom_range_PID_avg_noSCE.clear();
  mom_range_PID_pl2_noSCE.clear();
  mom_range_truePDG_noSCE.clear();

  vtx_x.clear();
  vtx_y.clear();
  vtx_z.clear();
  start_x.clear();
  start_y.clear();
  start_z.clear();
  end_x.clear();
  end_y.clear();
  end_z.clear();
  trk_phi.clear();
  trk_theta.clear();
  trk_costheta.clear();
  trk_length_pl0.clear();
  trk_length_pl1.clear();
  trk_length_pl2.clear();
  trk_length_noSCE.clear();
  trk_ifcontained.clear();
  vtx_FV.clear();
  trk_end_theta_yz.clear();
  trk_end_costheta_yz.clear();
  trk_end_theta_xz.clear();
  trk_end_costheta_xz.clear();
  trk_theta_yz.clear();
  trk_costheta_yz.clear();
  trk_theta_xz.clear();
  trk_costheta_xz.clear();

  hits_dEdx_size_pl0.clear();
  hits_dEdx_size_pl1.clear();
  hits_dEdx_size_pl2.clear();

  dEdx_pl0.clear();
  dEdx_pl1.clear();
  dEdx_pl2.clear();
  dQdx_pl0.clear();
  dQdx_pl1.clear();
  dQdx_pl2.clear();
  resRange_pl0.clear();
  resRange_pl1.clear();
  resRange_pl2.clear();

  dEdx_pl0_start_half.clear();
  dEdx_pl1_start_half.clear();
  dEdx_pl2_start_half.clear();
  dEdx_pl0_end_half.clear();
  dEdx_pl1_end_half.clear();
  dEdx_pl2_end_half.clear();
  dEdx_pl0_start10.clear();
  dEdx_pl1_start10.clear();
  dEdx_pl2_start10.clear();
  dEdx_pl0_end10.clear();
  dEdx_pl1_end10.clear();
  dEdx_pl2_end10.clear();
  dEdx_pl0_start1020.clear();
  dEdx_pl1_start1020.clear();
  dEdx_pl2_start1020.clear();
  dEdx_pl0_end1020.clear();
  dEdx_pl1_end1020.clear();
  dEdx_pl2_end1020.clear();

  PID_Chi2Mu_pl0.clear();
  PID_Chi2Mu_pl1.clear();
  PID_Chi2Mu_pl2.clear();
  PID_Chi2Mu_3pl.clear();

  PID_Chi2P_pl0.clear();
  PID_Chi2P_pl1.clear();
  PID_Chi2P_pl2.clear();
  PID_Chi2P_3pl.clear();

  PID_Chi2Pi_pl0.clear();
  PID_Chi2Pi_pl1.clear();
  PID_Chi2Pi_pl2.clear();
  PID_Chi2Pi_3pl.clear();
  
  PID_Chi2K_pl0.clear();
  PID_Chi2K_pl1.clear();
  PID_Chi2K_pl2.clear();
  PID_Chi2K_3pl.clear();

  PID_Pdg_allPlane.clear();
  PID_Pdg_pl2.clear();
  PID_Pdg_pl1.clear();
  PID_Pdg_pl0.clear();
  PID_avg_Chi2.clear();
  PID_pl2_Chi2.clear();
  PID_pl1_Chi2.clear();
  PID_pl0_Chi2.clear();

  Pl0_for_PID.clear();
  Pl1_for_PID.clear();
  Pl2_for_PID.clear();
  BestPlane_PID.clear();

  if_fwd_MCS.clear();
  if_fwd_true.clear();
  if_fwd_dEdx10.clear();
  if_fwd_dEdx1020.clear();
  if_fwd_dEdxhalf.clear();
}

void ParticleThreshold::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;

  // Make a tree to store selection information
  my_event_ = tfs->make<TTree>("tree","tree");

  if(IfAll){
    my_event_->Branch("All_true_PDG", &All_true_PDG);
    my_event_->Branch("All_true_mom", &All_true_mom);
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
    my_event_->Branch("All_true_trk_theta_yz", &All_true_trk_theta_yz);
    my_event_->Branch("All_true_trk_costheta_yz", &All_true_trk_costheta_yz);
    my_event_->Branch("All_true_trk_theta_xz", &All_true_trk_theta_xz);
    my_event_->Branch("All_true_trk_costheta_xz", &All_true_trk_costheta_xz);
    my_event_->Branch("All_true_trk_length", &All_true_trk_length);
  }
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

  my_event_->Branch("mom_bestMCS_mu", &mom_bestMCS_mu);
  my_event_->Branch("mom_bestMCS_ll_mu", &mom_bestMCS_ll_mu);
  my_event_->Branch("mom_fwdMCS_mu", &mom_fwdMCS_mu);
  my_event_->Branch("mom_fwdMCS_ll_mu", &mom_fwdMCS_ll_mu);
  my_event_->Branch("mom_bkwdMCS_mu", &mom_bkwdMCS_mu);
  my_event_->Branch("mom_bkwdMCS_ll_mu", &mom_bkwdMCS_ll_mu);
  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
  my_event_->Branch("mom_Range_p", &mom_Range_p);
  my_event_->Branch("mom_Range_pi", &mom_Range_pi);
  my_event_->Branch("mom_range_PID_avg", &mom_range_PID_avg);
  my_event_->Branch("mom_range_PID_pl2", &mom_range_PID_pl2);
  my_event_->Branch("mom_range_truePDG", &mom_range_truePDG);

  my_event_->Branch("mom_Range_mu_noSCE", &mom_Range_mu_noSCE);
  my_event_->Branch("mom_Range_p_noSCE", &mom_Range_p_noSCE);
  my_event_->Branch("mom_Range_pi_noSCE", &mom_Range_pi_noSCE);
  my_event_->Branch("mom_range_PID_avg_noSCE", &mom_range_PID_avg_noSCE);
  my_event_->Branch("mom_range_PID_pl2_noSCE", &mom_range_PID_pl2_noSCE);
  my_event_->Branch("mom_range_truePDG_noSCE", &mom_range_truePDG_noSCE);

  my_event_->Branch("vtx_x", &vtx_x);
  my_event_->Branch("vtx_y", &vtx_y);
  my_event_->Branch("vtx_z", &vtx_z);
  my_event_->Branch("start_x", &start_x);
  my_event_->Branch("start_y", &start_y);
  my_event_->Branch("start_z", &start_z);
  my_event_->Branch("end_x", &end_x);
  my_event_->Branch("end_y", &end_y);
  my_event_->Branch("end_z", &end_z);
  my_event_->Branch("trk_phi", &trk_phi);
  my_event_->Branch("trk_theta", &trk_theta);
  my_event_->Branch("trk_costheta", &trk_costheta);
  my_event_->Branch("trk_length_pl0", &trk_length_pl0);
  my_event_->Branch("trk_length_pl1", &trk_length_pl1);
  my_event_->Branch("trk_length_pl2", &trk_length_pl2);
  my_event_->Branch("trk_length_avg", &trk_length_avg);
  my_event_->Branch("trk_length_noSCE", &trk_length_noSCE);
  my_event_->Branch("trk_ifcontained", &trk_ifcontained);
  my_event_->Branch("vtx_FV", &vtx_FV);
  my_event_->Branch("trk_end_theta_yz", &trk_end_theta_yz);
  my_event_->Branch("trk_end_costheta_yz", &trk_end_costheta_yz);
  my_event_->Branch("trk_end_theta_xz", &trk_end_theta_xz);
  my_event_->Branch("trk_end_costheta_xz", &trk_end_costheta_xz);
  my_event_->Branch("trk_theta_yz", &trk_theta_yz);
  my_event_->Branch("trk_costheta_yz", &trk_costheta_yz);
  my_event_->Branch("trk_theta_xz", &trk_theta_xz);
  my_event_->Branch("trk_costheta_xz", &trk_costheta_xz);

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

  my_event_->Branch("PID_Pdg_allPlane", &PID_Pdg_allPlane);
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


void ParticleThreshold::beginJob()
{
  // Implementation of optional member function here.
  Initialize_event();
}

void ParticleThreshold::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleThreshold)

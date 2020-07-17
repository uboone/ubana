////////////////////////////////////////////////////////////////////////
// Class:       SingleParticleThreshold
// Plugin Type: analyzer (art v3_01_02)
// File:        SingleParticleThreshold_module.cc
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
#include "lardataobj/RecoBase/SpacePoint.h"
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
#include "PID.h"

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
  std::vector<double> All_true_trk_theta_yz;
  std::vector<double> All_true_trk_costheta_yz;
  std::vector<double> All_true_trk_theta_xz;
  std::vector<double> All_true_trk_costheta_xz;
  std::vector<double> All_true_trk_length;
 
  bool hadron_reconstructed; // if the event has proton or pion reconstructed (1mu1p or 1mu1pi)
  bool mu_reconstructed; // if the event has muon reconstructed (1mu1p or 1mu1pi)
 
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

  std::vector<double> dist_reco_hadron_start_true_vtx; // The distance of reco p(pi) to true vtx distance
  std::vector<double> dist_reco_mu_start_true_vtx; // The distance of reco mu to true vtx distance
  std::vector<double> dist_reco_mu_hadron_start; // The distance of reco p(pi) and mu start distance

  std::vector<double> mom_Range_mu;//Range momentum of muon track in the every event
  std::vector<double> mom_Range_p;//Range momentum of proton track in the every event
  std::vector<double> mom_Range_pi;//Range momentum of pion track in the every event

  std::vector<double> true_mu_mom;//True momentum of muon track in the every event (reconstructed)
  std::vector<double> true_hadron_mom;//True momentum of hadron track in the every event (reconstructed)
  std::vector<double> true_mu_purity; // Reco-True purity in the muon track
  std::vector<double> true_hadron_purity; // Reco-True purity in the hadron track
  std::vector<double> true_mu_completeness; // Reco-True completeness in the muon track
  std::vector<double> true_hadron_completeness; // Reco-True completeness in the hadron track

  std::vector<double> start_mu_x;//Reconstructed start x in the every event
  std::vector<double> start_mu_y;//Reconstructed start y in the every event
  std::vector<double> start_mu_z;//Reconstructed start z in the every event
  std::vector<double> end_mu_x;//Reconstructed end x in the every event
  std::vector<double> end_mu_y;//Reconstructed end y in the every event
  std::vector<double> end_mu_z;//Reconstructed end z in the every event
  std::vector<double> trk_mu_phi;//Reconstructed track phi in the every event
  std::vector<double> trk_mu_theta;//Reconstructed track theta in the every event
  std::vector<double> trk_mu_costheta;//Reconstructed track cos(theta) in the every event
  std::vector<double> trk_mu_length;//Range momentum of muon track in the every event no spatial correction
  std::vector<bool> trk_mu_ifcontained;//to check if the track is contained or not
  std::vector<bool> vtx_mu_FV;//to check if the vertex is in FV or not
  std::vector<double> trk_mu_theta_yz;
  std::vector<double> trk_mu_costheta_yz;
  std::vector<double> trk_mu_theta_xz;
  std::vector<double> trk_mu_costheta_xz;

  std::vector<double> start_hadron_x;//Reconstructed start x in the every event
  std::vector<double> start_hadron_y;//Reconstructed start y in the every event
  std::vector<double> start_hadron_z;//Reconstructed start z in the every event
  std::vector<double> end_hadron_x;//Reconstructed end x in the every event
  std::vector<double> end_hadron_y;//Reconstructed end y in the every event
  std::vector<double> end_hadron_z;//Reconstructed end z in the every event
  std::vector<double> trk_hadron_phi;//Reconstructed track phi in the every event
  std::vector<double> trk_hadron_theta;//Reconstructed track theta in the every event
  std::vector<double> trk_hadron_costheta;//Reconstructed track cos(theta) in the every event
  std::vector<double> trk_hadron_length;//Range momentum of muon track in the every event no spatial correction
  std::vector<bool> trk_hadron_ifcontained;//to check if the track is contained or not
  std::vector<bool> vtx_hadron_FV;//to check if the vertex is in FV or not
  std::vector<double> trk_hadron_theta_yz;
  std::vector<double> trk_hadron_costheta_yz;
  std::vector<double> trk_hadron_theta_xz;
  std::vector<double> trk_hadron_costheta_xz;
 
  int Ntrack; // number of tracks in this event
  int N_primary_tracks; // number of tracks matched with primary MC particles
  int N_MCP; // number of MC particles in this event
  int vtx_dist; // distance in between tracks (default = 99999)
  double cos_ang_muon_hadron; // angle in between muon and the hadron (default = 99999)

  std::vector<std::vector<float>> dEdx_pl0; // dE/dx of the selected (muon) track from plane 0 (closest to drift)
  std::vector<std::vector<float>> dEdx_pl1; // dE/dx of the selected (muon) track from plane 1
  std::vector<std::vector<float>> dEdx_pl2; // dE/dx of the selected (muon) track from plane 2 (collection)
  std::vector<std::vector<float>> resRange_pl0; // range from a hit to the end of the selected track end
  std::vector<std::vector<float>> resRange_pl1; // range from a hit to the end of the selected track end
  std::vector<std::vector<float>> resRange_pl2; // range from a hit to the end of the selected track end

  std::vector<double> PID_mu_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
  std::vector<double> PID_mu_Chi2P_3pl; // Chi2 of muon assumption of 3 planes in PID
  std::vector<double> PID_mu_Chi2Pi_3pl; // Chi2 of muon assumption of 3 planes in PID
  std::vector<double> PID_mu_Chi2K_3pl; // Chi2 of muon assumption of 3 planes in PID

  std::vector<double> PID_hadron_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID
  std::vector<double> PID_hadron_Chi2P_3pl; // Chi2 of muon assumption of 3 planes in PID
  std::vector<double> PID_hadron_Chi2Pi_3pl; // Chi2 of muon assumption of 3 planes in PID
  std::vector<double> PID_hadron_Chi2K_3pl; // Chi2 of muon assumption of 3 planes in PID

  std::vector<bool> if_mu_fwd; // If the reconstructed muon has the same direction as the simulated one
  std::vector<bool> if_hadron_fwd; // If the reconstructed hadron has the same direction as the simulated one

  std::vector<double> hit_SP_x; // The x coordinate of the spacepoint of hits
  std::vector<double> hit_SP_y; // The y coordinate of the spacepoint of hits
  std::vector<double> hit_SP_z; // The z coordinate of the spacepoint of hits
  std::vector<double> hit_charge; // The charge of hits

  std::string                         m_geantLabel;
  std::string                         m_pandoraLabel;
  std::string                         m_hitProducerLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_showerProducerLabel;
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
  m_geantLabel(pset.get<std::string>("GeantLabel")),
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_hitProducerLabel(pset.get<std::string>("HitProducerLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel")),
  m_calorimetryProducerLabel(pset.get<std::string>("calorimetryProducerLabel")),
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
  // Hit
  art::Handle<std::vector<recob::Hit> > Handle_Hit;
  evt.getByLabel(m_hitProducerLabel, Handle_Hit);
  std::vector< art::Ptr<recob::Hit> > AllHitCollection;
  art::fill_ptr_vector(AllHitCollection, Handle_Hit);

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

  // SpacePoint Hit association
  art::FindManyP<recob::SpacePoint> SpacePointToHitAsso(Handle_Hit, evt, m_pandoraLabel);
  
  //PID
  art::FindManyP<anab::ParticleID> PIDTotrackAsso(Handle_TPCtrack,evt,PID_TrackAssLabel);

  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;

  art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;
  evt.getByLabel(m_geantLabel, Handle_MCParticle);
  art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

  //------Loop over all MCParticles
  N_MCP = 0;
  TVector3 mu_start;
  TVector3 hadron_start;
  TVector3 mu_dir;
  TVector3 hadron_dir;
  TVector3 reco_mu_start;
  TVector3 reco_hadron_start;
  TVector3 All_true_start;
  TVector3 All_true_end;

  for(int i_mcp = 0; i_mcp < (int) MCParticleCollection.size(); i_mcp++){
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
      All_true_vtxFV.push_back(_fiducial_volume.VertexInFV(All_true_start));

      All_true_trk_phi.push_back(AllTrueTrackPos.Phi());
      All_true_trk_theta.push_back(AllTrueTrackPos.Theta());
      All_true_trk_costheta.push_back(cos(AllTrueTrackPos.Theta()));
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
    vtx_dist = 99999;
    cos_ang_muon_hadron = 99999;
  }

  //------Loop over all tracks
  Ntrack = AllTrackCollection.size();
  N_primary_tracks = 0;
  int trk_key = -9999; // it should later be filled with the key of the only/last track in the loop 
  for(int trk_id = 0; trk_id < Ntrack; trk_id++){

    //---Fill in True MCParticle information
    std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[trk_id].key());
    BackTrackerTruthMatch backtrackertruthmatch;
    backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
    auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
    if(!MCparticle){
      std::cout<<"MC particle does not exist!"<<std::endl;
    }
    else{
      if(MCparticle->Process() == "primary"){
        N_primary_tracks++;
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
        true_trk_ifcontained.push_back(_fiducial_volume.TrackContain(true_start, true_end));
        true_vtxFV.push_back(_fiducial_volume.VertexInFV(true_start));
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

      //---Fill in Reco information
      // Correct the position of the track ends
      TVector3 Trk_start = AllTrackCollection[trk_id]->Vertex<TVector3>();
      auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
      TVector3 Trk_start_SCEcorr;
      Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X() + 0.6);
      Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
      Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

      TVector3 Trk_end = AllTrackCollection[trk_id]->End<TVector3>();
      auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
      TVector3 Trk_end_SCEcorr;
      Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X() + 0.6);
      Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
      Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

      // induction 0 = 0, induction 1 = 1, collection = 2 (Check if the plane ID is correct)
      // assoCal[id_pl]->PlaneID().Plane == 2 (collection)
      // The vector is ordered by residual range from small to big (track end to track start)
      auto assoCal = trackToCalAsso.at(AllTrackCollection[trk_id].key());
      int hits_dEdx_size_pl0 = assoCal[0]->dEdx().size();
      int hits_dEdx_size_pl1 = assoCal[1]->dEdx().size();
      int hits_dEdx_size_pl2 = assoCal[2]->dEdx().size();
      PID pid;
      pid.Chi2(PIDTotrackAsso,AllTrackCollection[trk_id], Trk_start_SCEcorr, Trk_end_SCEcorr,hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2);

      if(MCparticle->Process() == "primary" && MCparticle->PdgCode() == 13){
        true_mu_mom.push_back(MCparticle->P());
        true_mu_purity.push_back(backtrackertruthmatch.ReturnPurity());
        true_mu_completeness.push_back(backtrackertruthmatch.ReturnCompleteness());

        mu_reconstructed = true;
        reco_mu_start = Trk_start_SCEcorr;
        trk_mu_ifcontained.push_back(_fiducial_volume.TrackContain(Trk_start_SCEcorr, Trk_end_SCEcorr));
        vtx_mu_FV.push_back(_fiducial_volume.VertexInFV(Trk_start_SCEcorr));
        trk_mu_length.push_back(AllTrackCollection[trk_id]->Length());
 
        start_mu_x.push_back(Trk_start_SCEcorr.X());
        start_mu_y.push_back(Trk_start_SCEcorr.Y());
        start_mu_z.push_back(Trk_start_SCEcorr.Z());
        end_mu_x.push_back(Trk_end_SCEcorr.X());
        end_mu_y.push_back(Trk_end_SCEcorr.Y());
        end_mu_z.push_back(Trk_end_SCEcorr.Z());

        trk_mu_phi.push_back(AllTrackCollection[trk_id]->Phi());
        trk_mu_theta.push_back(AllTrackCollection[trk_id]->Theta());
        trk_mu_costheta.push_back(cos(AllTrackCollection[trk_id]->Theta()));
       
        mom_Range_mu.push_back(_trk_mom_calculator.GetTrackMomentum(trk_mu_length.back(), 13));

        PID_mu_Chi2Mu_3pl.push_back(pid.PID_Chi2Mu_3pl); // Chi2 of muon assumption of 3 planes in PID
        PID_mu_Chi2P_3pl.push_back(pid.PID_Chi2P_3pl); // Chi2 of proton assumption of 3 planes in PID
        PID_mu_Chi2Pi_3pl.push_back(pid.PID_Chi2Pi_3pl); // Chi2 of pion assumption of 3 planes in PID
        PID_mu_Chi2K_3pl.push_back(pid.PID_Chi2K_3pl); // Chi2 of kaon assumption of 3 planes in PID

        // Check the directional info of the track by true reco vertex distance
        TVector3 vtx(All_true_start_x.back(), All_true_start_y.back(), All_true_start_z.back());
        TVector3 D_start = Trk_start_SCEcorr - vtx;
        TVector3 D_end = Trk_end_SCEcorr - vtx;
        if (D_start.Mag() < D_end.Mag()) if_mu_fwd.push_back(true);
        else if_mu_fwd.push_back(false);
      }
      if(MCparticle->Process() == "primary" && MCparticle->PdgCode() != 13){
        true_hadron_mom.push_back(MCparticle->P());
        true_hadron_purity.push_back(backtrackertruthmatch.ReturnPurity());
        true_hadron_completeness.push_back(backtrackertruthmatch.ReturnCompleteness());

        hadron_reconstructed = true;
        reco_hadron_start = Trk_start_SCEcorr;
        trk_hadron_ifcontained.push_back(_fiducial_volume.TrackContain(Trk_start_SCEcorr, Trk_end_SCEcorr));
        vtx_hadron_FV.push_back(_fiducial_volume.VertexInFV(Trk_start_SCEcorr));
        trk_hadron_length.push_back(AllTrackCollection[trk_id]->Length());

        start_hadron_x.push_back(Trk_start_SCEcorr.X());
        start_hadron_y.push_back(Trk_start_SCEcorr.Y());
        start_hadron_z.push_back(Trk_start_SCEcorr.Z());
        end_hadron_x.push_back(Trk_end_SCEcorr.X());
        end_hadron_y.push_back(Trk_end_SCEcorr.Y());
        end_hadron_z.push_back(Trk_end_SCEcorr.Z());

        trk_hadron_phi.push_back(AllTrackCollection[trk_id]->Phi());
        trk_hadron_theta.push_back(AllTrackCollection[trk_id]->Theta());
        trk_hadron_costheta.push_back(cos(AllTrackCollection[trk_id]->Theta()));

        if(MCparticle->PdgCode() == 2212){
          mom_Range_p.push_back(_trk_mom_calculator.GetTrackMomentum(trk_hadron_length.back(), 2212));
        }
        if(abs(MCparticle->PdgCode()) == 211){
          mom_Range_pi.push_back(_trk_mom_calculator.GetTrackMomentum(trk_hadron_length.back(), 211));
        }

        PID_hadron_Chi2Mu_3pl.push_back(pid.PID_Chi2Mu_3pl); // Chi2 of muon assumption of 3 planes in PID
        PID_hadron_Chi2P_3pl.push_back(pid.PID_Chi2P_3pl); // Chi2 of proton assumption of 3 planes in PID
        PID_hadron_Chi2Pi_3pl.push_back(pid.PID_Chi2Pi_3pl); // Chi2 of pion assumption of 3 planes in PID
        PID_hadron_Chi2K_3pl.push_back(pid.PID_Chi2K_3pl); // Chi2 of kaon assumption of 3 planes in PID

        // Check the directional info of the track by true reco vertex distance
        TVector3 vtx(All_true_start_x.back(), All_true_start_y.back(), All_true_start_z.back());
        TVector3 D_start = Trk_start_SCEcorr - vtx;
        TVector3 D_end = Trk_end_SCEcorr - vtx;
        if (D_start.Mag() < D_end.Mag()) if_hadron_fwd.push_back(true);
        else if_hadron_fwd.push_back(false);
      }

      trk_key = trk_id;
     
    } // if it's MC particle
  } // loop of reconstructed tracks
 
  // If there is only one primary MC particle reconstructed, read the dE/dx
  if(N_primary_tracks == 1){

    auto assoCal = trackToCalAsso.at(AllTrackCollection[trk_key].key());

    dEdx_pl0.push_back(assoCal[0]->dEdx());
    resRange_pl0.push_back(assoCal[0]->ResidualRange());

    dEdx_pl1.push_back(assoCal[1]->dEdx());
    resRange_pl1.push_back(assoCal[1]->ResidualRange());

    dEdx_pl2.push_back(assoCal[2]->dEdx());
    resRange_pl2.push_back(assoCal[2]->ResidualRange());
  }

  // Store the charge and position (space point) of the hits
  for(unsigned int i_hit = 0; i_hit < AllHitCollection.size(); i_hit++){
    auto assoSpacePoints = SpacePointToHitAsso.at(AllHitCollection[i_hit].key());
    if(assoSpacePoints.size()==1){
      TVector3 SP = (TVector3) assoSpacePoints.front()->XYZ();
      auto SP_offset = SCE->GetCalPosOffsets(geo::Point_t(SP.X(), SP.Y(), SP.Z()));
      TVector3 SP_SCEcorr;
      SP_SCEcorr.SetX(SP.X() - SP_offset.X());
      SP_SCEcorr.SetY(SP.Y() + SP_offset.Y());
      SP_SCEcorr.SetZ(SP.Z() + SP_offset.Z());
      hit_SP_x.push_back(SP_SCEcorr.X());
      hit_SP_y.push_back(SP_SCEcorr.Y());
      hit_SP_z.push_back(SP_SCEcorr.Z());
    }
    hit_charge.push_back(AllHitCollection[i_hit]->Integral());
  }

  // The distance of reco and true "vertex" 
  if(mu_reconstructed){    
    dist_reco_mu_start_true_vtx.push_back((All_true_start - reco_mu_start).Mag());
  }
  if(hadron_reconstructed){
    dist_reco_hadron_start_true_vtx.push_back((All_true_start - reco_hadron_start).Mag());
  }
  if(mu_reconstructed && hadron_reconstructed){
    dist_reco_mu_hadron_start.push_back((reco_mu_start - reco_hadron_start).Mag());
  }

  my_event_->Fill();

  mu_reconstructed = false;
  hadron_reconstructed = false;

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
  All_true_trk_theta_yz.clear();
  All_true_trk_costheta_yz.clear();
  All_true_trk_theta_xz.clear();
  All_true_trk_costheta_xz.clear();
  All_true_trk_length.clear();
 
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

  mom_Range_mu.clear();
  mom_Range_p.clear();
  mom_Range_pi.clear();

  true_mu_mom.clear();
  true_hadron_mom.clear();
  true_mu_purity.clear();
  true_hadron_purity.clear();
  true_mu_completeness.clear();
  true_hadron_completeness.clear();

  start_mu_x.clear();
  start_mu_y.clear();
  start_mu_z.clear();
  end_mu_x.clear();
  end_mu_y.clear();
  end_mu_z.clear();
  trk_mu_phi.clear();
  trk_mu_theta.clear();
  trk_mu_costheta.clear();
  trk_mu_length.clear();
  trk_mu_ifcontained.clear();
  vtx_mu_FV.clear();
  trk_mu_theta_yz.clear();
  trk_mu_costheta_yz.clear();
  trk_mu_theta_xz.clear();
  trk_mu_costheta_xz.clear();

  start_hadron_x.clear();
  start_hadron_y.clear();
  start_hadron_z.clear();
  end_hadron_x.clear();
  end_hadron_y.clear();
  end_hadron_z.clear();
  trk_hadron_phi.clear();
  trk_hadron_theta.clear();
  trk_hadron_costheta.clear();
  trk_hadron_length.clear();
  trk_hadron_ifcontained.clear();
  vtx_hadron_FV.clear();
  trk_hadron_theta_yz.clear();
  trk_hadron_costheta_yz.clear();
  trk_hadron_theta_xz.clear();
  trk_hadron_costheta_xz.clear();

  PID_mu_Chi2Mu_3pl.clear();
  PID_mu_Chi2P_3pl.clear();
  PID_mu_Chi2Pi_3pl.clear();
  PID_mu_Chi2K_3pl.clear();

  PID_hadron_Chi2Mu_3pl.clear();
  PID_hadron_Chi2P_3pl.clear();
  PID_hadron_Chi2Pi_3pl.clear();
  PID_hadron_Chi2K_3pl.clear();

  if_mu_fwd.clear();
  if_hadron_fwd.clear();

  dEdx_pl0.clear();
  dEdx_pl1.clear();
  dEdx_pl2.clear();
  resRange_pl0.clear();
  resRange_pl1.clear();
  resRange_pl2.clear();

  hit_SP_x.clear();
  hit_SP_y.clear();
  hit_SP_z.clear();
  hit_charge.clear();
}

void ParticleThreshold::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;

  // Make a tree to store selection information
  my_event_ = tfs->make<TTree>("tree","tree");

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
  my_event_->Branch("All_true_trk_theta_yz", &All_true_trk_theta_yz);
  my_event_->Branch("All_true_trk_costheta_yz", &All_true_trk_costheta_yz);
  my_event_->Branch("All_true_trk_theta_xz", &All_true_trk_theta_xz);
  my_event_->Branch("All_true_trk_costheta_xz", &All_true_trk_costheta_xz);
  my_event_->Branch("All_true_trk_length", &All_true_trk_length);

  my_event_->Branch("vtx_dist", &vtx_dist);
  my_event_->Branch("cos_ang_muon_hadron", &cos_ang_muon_hadron);
  my_event_->Branch("N_MCP", &N_MCP);
  my_event_->Branch("Ntrack", &Ntrack);
  my_event_->Branch("N_primary_tracks", &N_primary_tracks);

  my_event_->Branch("hadron_reconstructed", &hadron_reconstructed);
  my_event_->Branch("mu_reconstructed", &mu_reconstructed);

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

  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
  my_event_->Branch("mom_Range_p", &mom_Range_p);
  my_event_->Branch("mom_Range_pi", &mom_Range_pi);
  
  my_event_->Branch("dist_reco_hadron_start_true_vtx", &dist_reco_hadron_start_true_vtx);
  my_event_->Branch("dist_reco_mu_start_true_vtx", &dist_reco_mu_start_true_vtx);
  my_event_->Branch("dist_reco_mu_hadron_start", &dist_reco_mu_hadron_start);

  my_event_->Branch("true_mu_mom", &true_mu_mom);
  my_event_->Branch("true_hadron_mom", &true_hadron_mom);
  my_event_->Branch("true_mu_purity", &true_mu_purity);
  my_event_->Branch("true_hadron_purity", &true_hadron_purity);
  my_event_->Branch("true_mu_completeness", &true_mu_completeness);
  my_event_->Branch("true_hadron_completeness", &true_hadron_completeness);

  my_event_->Branch("start_mu_x", &start_mu_x);
  my_event_->Branch("start_mu_y", &start_mu_y);
  my_event_->Branch("start_mu_z", &start_mu_z);
  my_event_->Branch("end_mu_x", &end_mu_x);
  my_event_->Branch("end_mu_y", &end_mu_y);
  my_event_->Branch("end_mu_z", &end_mu_z);
  my_event_->Branch("trk_mu_phi", &trk_mu_phi);
  my_event_->Branch("trk_mu_theta", &trk_mu_theta);
  my_event_->Branch("trk_mu_costheta", &trk_mu_costheta);
  my_event_->Branch("trk_mu_length", &trk_mu_length);
  my_event_->Branch("trk_mu_ifcontained", &trk_mu_ifcontained);
  my_event_->Branch("vtx_mu_FV", &vtx_mu_FV);
  my_event_->Branch("trk_mu_theta_yz", &trk_mu_theta_yz);
  my_event_->Branch("trk_mu_costheta_yz", &trk_mu_costheta_yz);
  my_event_->Branch("trk_mu_theta_xz", &trk_mu_theta_xz);
  my_event_->Branch("trk_mu_costheta_xz", &trk_mu_costheta_xz);

  my_event_->Branch("start_hadron_x", &start_hadron_x);
  my_event_->Branch("start_hadron_y", &start_hadron_y);
  my_event_->Branch("start_hadron_z", &start_hadron_z);
  my_event_->Branch("end_hadron_x", &end_hadron_x);
  my_event_->Branch("end_hadron_y", &end_hadron_y);
  my_event_->Branch("end_hadron_z", &end_hadron_z);
  my_event_->Branch("trk_hadron_phi", &trk_hadron_phi);
  my_event_->Branch("trk_hadron_theta", &trk_hadron_theta);
  my_event_->Branch("trk_hadron_costheta", &trk_hadron_costheta);
  my_event_->Branch("trk_hadron_length", &trk_hadron_length);
  my_event_->Branch("trk_hadron_ifcontained", &trk_hadron_ifcontained);
  my_event_->Branch("vtx_hadron_FV", &vtx_hadron_FV);
  my_event_->Branch("trk_hadron_theta_yz", &trk_hadron_theta_yz);
  my_event_->Branch("trk_hadron_costheta_yz", &trk_hadron_costheta_yz);
  my_event_->Branch("trk_hadron_theta_xz", &trk_hadron_theta_xz);
  my_event_->Branch("trk_hadron_costheta_xz", &trk_hadron_costheta_xz);

  my_event_->Branch("PID_mu_Chi2Mu_3pl", &PID_mu_Chi2Mu_3pl);
  my_event_->Branch("PID_mu_Chi2P_3pl", &PID_mu_Chi2P_3pl);
  my_event_->Branch("PID_mu_Chi2Pi_3pl", &PID_mu_Chi2Pi_3pl);
  my_event_->Branch("PID_mu_Chi2K_3pl", &PID_mu_Chi2K_3pl);

  my_event_->Branch("PID_hadron_Chi2Mu_3pl", &PID_hadron_Chi2Mu_3pl);
  my_event_->Branch("PID_hadron_Chi2P_3pl", &PID_hadron_Chi2P_3pl);
  my_event_->Branch("PID_hadron_Chi2Pi_3pl", &PID_hadron_Chi2Pi_3pl);
  my_event_->Branch("PID_hadron_Chi2K_3pl", &PID_hadron_Chi2K_3pl);

  my_event_->Branch("if_mu_fwd", &if_mu_fwd);
  my_event_->Branch("if_hadron_fwd", &if_hadron_fwd);

  my_event_->Branch("dEdx_pl0", &dEdx_pl0);
  my_event_->Branch("dEdx_pl1", &dEdx_pl1);
  my_event_->Branch("dEdx_pl2", &dEdx_pl2);
  my_event_->Branch("resRange_pl0", &resRange_pl0);
  my_event_->Branch("resRange_pl1", &resRange_pl1);
  my_event_->Branch("resRange_pl2", &resRange_pl2);

  my_event_->Branch("hit_SP_x", &hit_SP_x);
  my_event_->Branch("hit_SP_y", &hit_SP_y);
  my_event_->Branch("hit_SP_z", &hit_SP_z);
  my_event_->Branch("hit_charge", &hit_charge);
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

////////////////////////////////////////////////////////////////////////
// Class:       CCthreshold
// Plugin Type: analyzer (art v3_01_02)
// File:        CCthreshold_module.cc
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

class CCthreshold;


class CCthreshold : public art::EDAnalyzer {
public:
  explicit CCthreshold(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CCthreshold(CCthreshold const&) = delete;
  CCthreshold(CCthreshold&&) = delete;
  CCthreshold& operator=(CCthreshold const&) = delete;
  CCthreshold& operator=(CCthreshold&&) = delete;

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
 
  bool MC_beamNeutrino; // MCTruth beam origin
  bool MC_FV; // MCTruth vertex in FV = true, out of FV = false
  int MC_ccnc; // MCTruth cc = 0 or nc = 1
  int MC_nupdg; // MCTruth nupdg; numu = 14, nue = 12
  double MC_nuVtxX; // MCTruth nu vtx X
  double MC_nuVtxY; // MCTruth nu vtx Y
  double MC_nuVtxZ; // MCTruth nu vtx Z
 
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


CCthreshold::CCthreshold(fhicl::ParameterSet const& pset)
  : 
  EDAnalyzer{pset},
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

void CCthreshold::analyze(art::Event const& evt)
{

  //// Get necessary handles
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;

  // MC Truth
  art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
  evt.getByLabel(m_generatorLabel, Handle_MCTruth);
  art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

  // MC Particle
  art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;
  evt.getByLabel(m_geantLabel, Handle_MCParticle);
  art::fill_ptr_vector(MCParticleCollection, Handle_MCParticle);

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

  //Constants
  const simb::Origin_t Neutrino_Origin = simb::kBeamNeutrino;
  
  //focus on numu cc beam neutrino event
  MC_beamNeutrino = false;
  MC_nupdg = -99999;
  MC_ccnc = -99999;
  MC_FV = false;

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

  // Get all the neutrino daughters
  if(MC_ccnc == 0 && MC_nupdg == 14 && MC_beamNeutrino == true){
    //------Loop over all MCParticles
    for(int i_mcp = 0; i_mcp < (int) MCParticleCollection.size(); i_mcp++){
      if(MCParticleCollection[i_mcp]->Process() == "primary"){

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

        // If a reconstructed track corresponds to a MC particle which matches numu primary daughter, store the truth info of the reco track
        for(int trk_id = 0; trk_id < (int) AllTrackCollection.size(); trk_id++){
          // Backtracker
          std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[trk_id].key());
          BackTrackerTruthMatch backtrackertruthmatch;
          backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
          auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
          if(MCparticle && MCparticle == MCP){
            std::cout<<"Track matches with primary MC particle!"<<std::endl;
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
          } // If MC particles match
        } // Finish Looping tracks
      } // If Primary
    } // Loop of MC particles
  } // If beam numu CC 

  my_event_->Fill();
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

void CCthreshold::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;

  // Make a tree to store selection information
  my_event_ = tfs->make<TTree>("tree","tree");

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


void CCthreshold::beginJob()
{
  // Implementation of optional member function here.
  Initialize_event();
}

void CCthreshold::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CCthreshold)

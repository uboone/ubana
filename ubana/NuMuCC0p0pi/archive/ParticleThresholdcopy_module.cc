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
  std::vector<double> true_trk_length;//True track length (distance from the start to the end point) 
  std::vector<double> trk_pdg;//Track pdg 

  bool if_selected; // If selected based on the reco info
  bool if_matchMu; // If the selected track matched with true muon from numu cc
  bool if_cosmic; // Check if a track is cosmic or not by if it has an associated MCParticle

  std::vector<double> mom_bestMCS_mu;//MCS best momentum of muon track in the every event
  std::vector<double> mom_bestMCS_ll_mu;//Likelihood of MCS best momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_mu;//MCS forward momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_ll_mu;//Likelihood of MCS forward momentum of muon track in the every event
  std::vector<double> mom_bkwdMCS_mu;//MCS backward momentum of muon track in the every event
  std::vector<double> mom_bkwdMCS_ll_mu;//Likelihood of MCS backward momentum of muon track in the every event
  std::vector<double> mom_Range_mu;//Range momentum of muon track in the every event
  std::vector<double> mom_range_PID;//Range momentum of tracks based on their PID particle type
  std::vector<double> mom_range_truePDG;//Range momentum of tracks based on their true particle type
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
  std::vector<double> trk_length;//Range momentum of muon track in the every event
  std::vector<bool> trk_ifcontained;//to check if the track is contained or not
  std::vector<bool> vtx_FV;//to check if the vertex is in FV or not
  int Ntrack; // number of tracks in this event

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

  double PID_Chi2Mu_pl0; // Chi2 of muon assumption of plane 0 in PID
  double PID_Chi2Mu_pl1; // Chi2 of muon assumption of plane 1 in PID
  double PID_Chi2Mu_pl2; // Chi2 of muon assumption of plane 2 in PID
  double PID_Chi2P_pl0; // Chi2 of proton assumption of plane 0 in PID
  double PID_Chi2P_pl1; // Chi2 of proton assumption of plane 1 in PID
  double PID_Chi2P_pl2; // Chi2 of proton assumption of plane 2 in PID
  double PID_Chi2Pi_pl0; // Chi2 of pion assumption of plane 0 in PID
  double PID_Chi2Pi_pl1; // Chi2 of pion assumption of plane 1 in PID
  double PID_Chi2Pi_pl2; // Chi2 of pion assumption of plane 2 in PID
  double PID_Chi2K_pl0; // Chi2 of kaon assumption of plane 0 in PID
  double PID_Chi2K_pl1; // Chi2 of kaon assumption of plane 1 in PID
  double PID_Chi2K_pl2; // Chi2 of kaon assumption of plane 2 in PID
  int PID_Pdg; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
  double PID_avg_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions

  bool                                IsMC;
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

  //------Loop over all tracks
  Ntrack = AllTrackCollection.size();
  for(int trk_id = 0; trk_id < Ntrack; trk_id++){
    auto assoCal = trackToCalAsso.at(AllTrackCollection[trk_id].key());
    double Trk_Length = assoCal.front()->Range();  //It is said track length in pandoracali has spatial correction
    trk_length.push_back(Trk_Length); // track length with spatial correction

    //---Fill in True MCParticle information
    if(IsMC){
      std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[trk_id].key());
      BackTrackerTruthMatch backtrackertruthmatch;
      backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
      auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
      if(!MCparticle){
        if_cosmic = true;
        std::cout<<"MC particle does not exist!"<<std::endl;
      }
      else{
        if_cosmic = false;
        auto TrueTrackPos = MCparticle->EndPosition() - MCparticle->Position();
        true_mom.push_back(MCparticle->P());
        true_start_x.push_back(MCparticle->Position().X());
        true_start_y.push_back(MCparticle->Position().Y());
        true_start_z.push_back(MCparticle->Position().Z());
        true_end_x.push_back(MCparticle->EndPosition().X());
        true_end_y.push_back(MCparticle->EndPosition().Y());
        true_end_z.push_back(MCparticle->EndPosition().Z());
        true_trk_phi.push_back(TrueTrackPos.Phi());
        true_trk_theta.push_back(TrueTrackPos.Theta());
        true_trk_costheta.push_back(cos(TrueTrackPos.Theta()));
        true_trk_length.push_back(sqrt(TrueTrackPos.X()*TrueTrackPos.X() + TrueTrackPos.Y()*TrueTrackPos.Y() + TrueTrackPos.Z()*TrueTrackPos.Z())); // An estimation of true track length
        trk_pdg.push_back(MCparticle->PdgCode());

        //momentum by true PDG
        if (abs(MCparticle->PdgCode()) == 13){
          mom_range_truePDG.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 13));
        }
        if (abs(MCparticle->PdgCode()) == 2212){
          mom_range_truePDG.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 2212));
        }
        if (abs(MCparticle->PdgCode()) == 211){
          mom_range_truePDG.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 211));
        }
        if (abs(MCparticle->PdgCode()) == 321){
          mom_range_truePDG.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 321));
        }
      }
    }
   
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

    bool trk_contained = _fiducial_volume.InFV(Trk_start_SCEcorr, Trk_end_SCEcorr);
    trk_ifcontained.push_back(trk_contained);

    bool vtx_InFV = _fiducial_volume.InFV(Trk_start_SCEcorr);
    vtx_FV.push_back(vtx_InFV);

    double bestMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestMomentum();
    double bestMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestLogLikelihood();
    mom_bestMCS_mu.push_back(bestMCS);
    mom_bestMCS_ll_mu.push_back(bestMCSLL);

    double fwdMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdMomentum();
    double fwdMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdLogLikelihood();
    mom_fwdMCS_mu.push_back(fwdMCS);
    mom_fwdMCS_ll_mu.push_back(fwdMCSLL);

    double bkwdMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdMomentum();
    double bkwdMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdLogLikelihood();
    mom_bkwdMCS_mu.push_back(bkwdMCS);
    mom_bkwdMCS_ll_mu.push_back(bkwdMCSLL);

    //auto assoCal = trackToCalAsso.at(AllTrackCollection[trk_id].key());
    //double Trk_Length = assoCal.front()->Range();  //It is said track length in pandoracali has spatial correction
    //trk_length.push_back(Trk_Length); // track length with spatial correction

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
    trk_costheta.push_back(AllTrackCollection[trk_id]->Theta());

    double RangeMom_Mu_assumption = _trk_mom_calculator.GetTrackMomentum(Trk_Length, 13);
    mom_Range_mu.push_back(RangeMom_Mu_assumption);
  
    //Get Calorimety of the track
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

    int half_size_pl0 = dEdx_pl0.size() / 2;
    int half_size_pl1 = dEdx_pl1.size() / 2;
    int half_size_pl2 = dEdx_pl2.size() / 2;

    dEdx_pl0_start_half = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
    dEdx_pl0_end_half = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
    dEdx_pl1_start_half = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
    dEdx_pl1_end_half = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
    dEdx_pl2_start_half = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
    dEdx_pl2_end_half = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;
    if (dEdx_pl0.size()<10) {
      dEdx_pl0_start10 = dEdx_pl0_start_half;
      dEdx_pl0_end10 = dEdx_pl0_end_half;
    }
    else{
      dEdx_pl0_start10 = std::accumulate(dEdx_pl0.end() - 10, dEdx_pl0.end(), 0.) / 10.;
      dEdx_pl0_end10 = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + 10, 0.) / 10.;
    }
    if (dEdx_pl1.size()<10) {
      dEdx_pl1_start10 = dEdx_pl1_start_half;
      dEdx_pl1_end10 = dEdx_pl1_end_half;
    }
    else{
      dEdx_pl1_start10 = std::accumulate(dEdx_pl1.end() - 10, dEdx_pl1.end(), 0.) / 10.;
      dEdx_pl1_end10 = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + 10, 0.) / 10.;
    }
    if (dEdx_pl2.size()<10) {
      dEdx_pl2_start10 = dEdx_pl2_start_half;
      dEdx_pl2_end10 = dEdx_pl2_end_half;
    }
    else{
      dEdx_pl2_start10 = std::accumulate(dEdx_pl2.end() - 10, dEdx_pl2.end(), 0.) / 10.;
      dEdx_pl2_end10 = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + 10, 0.) / 10.;
    }

    // Gain PID info of the track
    if(!PIDTotrackAsso.isValid()){
      throw cet::exception("[Numu0pi0p]") << "No matched PID - track information!" << std::endl;
    }
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
    PID_Chi2Mu_pl0 = PIDChi2_mu[0];
    PID_Chi2Mu_pl1 = PIDChi2_mu[1];
    PID_Chi2Mu_pl2 = PIDChi2_mu[2];
    PID_Chi2P_pl0 = PIDChi2_p[0];
    PID_Chi2P_pl1 = PIDChi2_p[1];
    PID_Chi2P_pl2 = PIDChi2_p[2];
    PID_Chi2Pi_pl0 = PIDChi2_pi[0];
    PID_Chi2Pi_pl1 = PIDChi2_pi[1];
    PID_Chi2Pi_pl2 = PIDChi2_pi[2];
    PID_Chi2K_pl0 = PIDChi2_K[0];
    PID_Chi2K_pl1 = PIDChi2_K[1];
    PID_Chi2K_pl2 = PIDChi2_K[2];

    std::vector<double> PIDChi2_avg (4, 0.); // It follows the order of muon, proton, pion, kaon
    for(int id_pl = 0; id_pl < 3; id_pl++){
      PIDChi2_avg[0] += PIDChi2_mu[id_pl];
      PIDChi2_avg[1] += PIDChi2_p[id_pl];
      PIDChi2_avg[2] += PIDChi2_pi[id_pl];
      PIDChi2_avg[3] += PIDChi2_K[id_pl];
    }

    for(int index = 0; index < (int) PIDChi2_avg.size(); index++){
      PIDChi2_avg[index] = PIDChi2_avg[index] / 3; // average the Chi2 of each assumption of the 3 planes
    }
    PID_avg_Chi2 = *std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end());
    int ID_PID = std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end()) - PIDChi2_avg.begin();
    if (ID_PID == 0) {
      PID_Pdg = 13;
      mom_range_PID.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 13));
    }
    if (ID_PID == 1) {
      PID_Pdg = 2212;
      mom_range_PID.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 2212));
    }
    if (ID_PID == 2) {
      PID_Pdg = 211;
      mom_range_PID.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 211));
    }
    if (ID_PID == 3) {
      PID_Pdg = 321; 
      mom_range_PID.push_back(_trk_mom_calculator.GetTrackMomentum(Trk_Length, 321));
    }

    my_event_->Fill();

  } // The end of track loops

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
    true_trk_length.clear();
    trk_pdg.clear();
  }

  mom_bestMCS_mu.clear();
  mom_bestMCS_ll_mu.clear();
  mom_Range_mu.clear();
  mom_range_PID.clear();
  mom_range_truePDG.clear();
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
  trk_length.clear();
  trk_ifcontained.clear();
  vtx_FV.clear();
  
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

void ParticleThreshold::Initialize_event()
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
    my_event_->Branch("true_trk_length", &true_trk_length);
    my_event_->Branch("trk_pdg", &trk_pdg);
  }

  my_event_->Branch("if_cosmic", &if_cosmic);
  my_event_->Branch("if_matchMu", &if_matchMu);
  my_event_->Branch("if_selected", &if_selected);

  my_event_->Branch("mom_bestMCS_mu", &mom_bestMCS_mu);
  my_event_->Branch("mom_bestMCS_ll_mu", &mom_bestMCS_ll_mu);
  my_event_->Branch("mom_fwdMCS_mu", &mom_fwdMCS_mu);
  my_event_->Branch("mom_fwdMCS_ll_mu", &mom_fwdMCS_ll_mu);
  my_event_->Branch("mom_bkwdMCS_mu", &mom_bkwdMCS_mu);
  my_event_->Branch("mom_bkwdMCS_ll_mu", &mom_bkwdMCS_ll_mu);
  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
  my_event_->Branch("mom_range_PID", &mom_range_PID);
  my_event_->Branch("mom_range_truePDG", &mom_range_truePDG);
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
  my_event_->Branch("trk_length", &trk_length);
  my_event_->Branch("trk_ifcontained", &trk_ifcontained);
  my_event_->Branch("vtx_FV", &vtx_FV);

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
  
  my_event_->Branch("PID_Chi2Mu_pl0", &PID_Chi2Mu_pl0);
  my_event_->Branch("PID_Chi2Mu_pl1", &PID_Chi2Mu_pl1);
  my_event_->Branch("PID_Chi2Mu_pl2", &PID_Chi2Mu_pl2);
  my_event_->Branch("PID_Chi2P_pl0", &PID_Chi2P_pl0);
  my_event_->Branch("PID_Chi2P_pl1", &PID_Chi2P_pl1);
  my_event_->Branch("PID_Chi2P_pl2", &PID_Chi2P_pl2);
  my_event_->Branch("PID_Chi2Pi_pl0", &PID_Chi2Pi_pl0);
  my_event_->Branch("PID_Chi2Pi_pl1", &PID_Chi2Pi_pl1);
  my_event_->Branch("PID_Chi2Pi_pl2", &PID_Chi2Pi_pl2);
  my_event_->Branch("PID_Chi2K_pl0", &PID_Chi2K_pl0);
  my_event_->Branch("PID_Chi2K_pl1", &PID_Chi2K_pl1);
  my_event_->Branch("PID_Chi2K_pl2", &PID_Chi2K_pl2);
  my_event_->Branch("PID_Pdg", &PID_Pdg);
  my_event_->Branch("PID_avg_Chi2", &PID_avg_Chi2);
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

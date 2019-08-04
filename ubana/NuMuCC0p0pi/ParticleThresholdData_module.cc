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
  
  double mom_bestMCS_mu;//MCS best momentum of muon track in the every event
  double mom_bestMCS_ll_mu;//Likelihood of MCS best momentum of muon track in the every event
  double mom_fwdMCS_mu;//MCS forward momentum of muon track in the every event
  double mom_fwdMCS_ll_mu;//Likelihood of MCS forward momentum of muon track in the every event
  double mom_bkwdMCS_mu;//MCS backward momentum of muon track in the every event
  double mom_bkwdMCS_ll_mu;//Likelihood of MCS backward momentum of muon track in the every event
  double mom_Range_mu;//Range momentum of muon track in the every event
  double mom_Range_p;//Range momentum of proton track in the every event
  double mom_Range_pi;//Range momentum of pion track in the every event
  double mom_range_PID_avg;//Range momentum of tracks based on their PID particle type using 3 pls
  double mom_range_PID_pl2;//Range momentum of tracks based on their PID particle type using pl 2
  
  double mom_bestMCS_mu_noSCE;//MCS best momentum of muon track in the every event
  double mom_bestMCS_ll_mu_noSCE;//Likelihood of MCS best momentum of muon track in the every event
  double mom_fwdMCS_mu_noSCE;//MCS forward momentum of muon track in the every event
  double mom_fwdMCS_ll_mu_noSCE;//Likelihood of MCS forward momentum of muon track in the every event
  double mom_bkwdMCS_mu_noSCE;//MCS backward momentum of muon track in the every event
  double mom_bkwdMCS_ll_mu_noSCE;//Likelihood of MCS backward momentum of muon track in the every event
  double mom_Range_mu_noSCE;//Range momentum of muon track in the every event
  double mom_Range_p_noSCE;//Range momentum of proton track in the every event
  double mom_Range_pi_noSCE;//Range momentum of pion track in the every event
  double mom_range_PID_avg_noSCE;//Range momentum of tracks based on their PID particle type using 3 pls
  double mom_range_PID_pl2_noSCE;//Range momentum of tracks based on their PID particle type using pl 2

  double vtx_x;//Reconstructed vtx x in the every event
  double vtx_y;//Reconstructed vtx y in the every event
  double vtx_z;//Reconstructed vtx z in the every event
  double start_x;//Reconstructed start x in the every event
  double start_y;//Reconstructed start y in the every event
  double start_z;//Reconstructed start z in the every event
  double end_x;//Reconstructed end x in the every event
  double end_y;//Reconstructed end y in the every event
  double end_z;//Reconstructed end z in the every event
  double trk_phi;//Reconstructed track phi in the every event
  double trk_theta;//Reconstructed track theta in the every event
  double trk_costheta;//Reconstructed track cos(theta) in the every event
  double trk_length_pl0;//plane 0 length of muon track in the every event
  double trk_length_pl1;//plane 1 length of muon track in the every event
  double trk_length_pl2;//plane 2 length of muon track in the every event
  double trk_length_avg;//average length of muon track over 3 planes in the every event
  double trk_length_noSCE;//Range momentum of muon track in the every event no spatial correction
  bool trk_ifcontained;//to check if the track is contained or not
  bool vtx_FV;//to check if the vertex is in FV or not
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
  int PID_Pdg_allPlane; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
  int PID_Pdg_pl2; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
  double PID_avg_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions
  double PID_pl2_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions

  bool if_fwd_MCS; // If using forward MCS direction judge
  bool if_fwd_dEdx10; // If fwd by the reco dEdx 10 hits (should use for contained)
  bool if_fwd_dEdxhalf; // If fwd by the reco dEdx half of the hits (should use for contained)

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

    trk_ifcontained = _fiducial_volume.InFV(Trk_start_SCEcorr, Trk_end_SCEcorr);

    vtx_FV = _fiducial_volume.InFV(Trk_start_SCEcorr);

    mom_bestMCS_mu =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestMomentum();
    mom_bestMCS_ll_mu =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestLogLikelihood();

    mom_fwdMCS_mu =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdMomentum();
    mom_fwdMCS_ll_mu =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdLogLikelihood();

    mom_bkwdMCS_mu =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdMomentum();
    mom_bkwdMCS_ll_mu =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdLogLikelihood();
    
    // track length with no SCE 
    trk_length_noSCE = AllTrackCollection[trk_id]->Length();
    auto assoCal = trackToCalAsso.at(AllTrackCollection[trk_id].key());
    trk_length_pl0 = assoCal[0]->Range();  //pandoracali has spatial correction
    trk_length_pl1 = assoCal[1]->Range();  //pandoracali has spatial correction
    trk_length_pl2 = assoCal[2]->Range();  //pandoracali has spatial correction
    trk_length_avg = 0;
    int valid_pl = 0;
    for (int i_pl = 0; i_pl < (int) assoCal.size(); i_pl++){
      if(assoCal[i_pl]->Range() > 0){
        trk_length_avg += assoCal[i_pl]->Range();
        valid_pl++;
      }
    }
    trk_length_avg = trk_length_avg / valid_pl;
    
    vtx_x = Trk_start_SCEcorr.X();
    vtx_y = Trk_start_SCEcorr.Y();
    vtx_z = Trk_start_SCEcorr.Z();
    start_x = Trk_start_SCEcorr.X();
    start_y = Trk_start_SCEcorr.Y();
    start_z = Trk_start_SCEcorr.Z();
    end_x = Trk_end_SCEcorr.X();
    end_y = Trk_end_SCEcorr.Y();
    end_z = Trk_end_SCEcorr.Z();

    trk_phi = AllTrackCollection[trk_id]->Phi();
    trk_theta = AllTrackCollection[trk_id]->Theta();
    trk_costheta = AllTrackCollection[trk_id]->Theta();

    mom_Range_mu = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2, 13);
    mom_Range_mu_noSCE = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE, 13);
 
    mom_Range_p = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2, 2212);
    mom_Range_p_noSCE = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE, 2212);

    mom_Range_pi = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2, 211);
    mom_Range_pi_noSCE = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE, 211);

    // Get PID info
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

    // Naively use all planes
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
      PID_Pdg_allPlane = 13;
    }
    if (ID_PID == 1) {
      PID_Pdg_allPlane = 2212;
    }
    if (ID_PID == 2) {
      PID_Pdg_allPlane = 211;
    }
    if (ID_PID == 3) {
      PID_Pdg_allPlane = 321;
    }
    mom_range_PID_avg = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2, PID_Pdg_allPlane);
    mom_range_PID_avg_noSCE = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE, PID_Pdg_allPlane);

    //-- Use plane 2 only
    std::vector<double> PIDChi2_pl2 {PID_Chi2Mu_pl2, PID_Chi2P_pl2, PID_Chi2Pi_pl2, PID_Chi2K_pl2}; // It follows the order of muon, proton, pion, kaon

    PID_pl2_Chi2 = *std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end());
    int ID_PID_pl2 = std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end()) - PIDChi2_pl2.begin();
    if (ID_PID_pl2 == 0) {
      PID_Pdg_pl2 = 13;
    }
    if (ID_PID_pl2 == 1) {
      PID_Pdg_pl2 = 2212;
    }
    if (ID_PID_pl2 == 2) {
      PID_Pdg_pl2 = 211;
    }
    if (ID_PID_pl2 == 3) {
      PID_Pdg_pl2 = 321;
    }
    mom_range_PID_pl2 = _trk_mom_calculator.GetTrackMomentum(trk_length_pl2, PID_Pdg_pl2);
    mom_range_PID_pl2_noSCE = _trk_mom_calculator.GetTrackMomentum(trk_length_noSCE, PID_Pdg_pl2);

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

    // Check the directional info of the track by MCS
    if (mom_bestMCS_ll_mu == mom_fwdMCS_ll_mu) if_fwd_MCS = true;
    else if_fwd_MCS = false;

    // Check the direction info of the track by dEdx
    if (dEdx_pl2_start_half < dEdx_pl2_end_half) if_fwd_dEdxhalf = true;
    else if_fwd_dEdxhalf = false;
    if (dEdx_pl2_start10 < dEdx_pl2_end10) if_fwd_dEdx10 = true;
    else if_fwd_dEdx10 = false;
    
    my_event_->Fill();

    dEdx_pl0.clear();
    dEdx_pl1.clear();
    dEdx_pl2.clear();
    dQdx_pl0.clear();
    dQdx_pl1.clear();
    dQdx_pl2.clear();
    resRange_pl0.clear();
    resRange_pl1.clear();
    resRange_pl2.clear();

  } // The end of track loops
}

void ParticleThreshold::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;

  // Make a tree to store selection information
  my_event_ = tfs->make<TTree>("tree","tree");

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

  my_event_->Branch("mom_bestMCS_mu_noSCE", &mom_bestMCS_mu_noSCE);
  my_event_->Branch("mom_bestMCS_ll_mu_noSCE", &mom_bestMCS_ll_mu_noSCE);
  my_event_->Branch("mom_fwdMCS_mu_noSCE", &mom_fwdMCS_mu_noSCE);
  my_event_->Branch("mom_fwdMCS_ll_mu_noSCE", &mom_fwdMCS_ll_mu_noSCE);
  my_event_->Branch("mom_bkwdMCS_mu_noSCE", &mom_bkwdMCS_mu_noSCE);
  my_event_->Branch("mom_bkwdMCS_ll_mu_noSCE", &mom_bkwdMCS_ll_mu_noSCE);
  my_event_->Branch("mom_Range_mu_noSCE", &mom_Range_mu_noSCE);
  my_event_->Branch("mom_Range_p_noSCE", &mom_Range_p_noSCE);
  my_event_->Branch("mom_Range_pi_noSCE", &mom_Range_pi_noSCE);
  my_event_->Branch("mom_range_PID_avg_noSCE", &mom_range_PID_avg_noSCE);
  my_event_->Branch("mom_range_PID_pl2_noSCE", &mom_range_PID_pl2_noSCE);

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
  my_event_->Branch("PID_Pdg_allPlane", &PID_Pdg_allPlane);
  my_event_->Branch("PID_Pdg_pl2", &PID_Pdg_pl2);
  my_event_->Branch("PID_avg_Chi2", &PID_avg_Chi2);
  my_event_->Branch("PID_pl2_Chi2", &PID_pl2_Chi2);
  
  my_event_->Branch("if_fwd_MCS", &if_fwd_MCS);
  my_event_->Branch("if_fwd_dEdx10", &if_fwd_dEdx10);
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

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

#include "nusimdata/SimulationBase/simb.h"
//#include "nutools/G4Base/PrimaryParticleInformation.h"
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
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  ::ubana::FiducialVolume _fiducial_volume;
 
  //::art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<geo::Geometry> geo;

  art::ServiceHandle<art::TFileService> tfs;
  TTree * my_event_;
  
  void Initialize_event();

  bool MC_beamNeutrino; // MCTruth beam origin
  bool MC_FV; // MCTruth vertex in FV = true, out of FV = false
  int MC_ccnc; // MCTruth cc = 0 or nc = 1
  int MC_nupdg; // MCTruth nupdg; numu = 14, nue = 12
  double MC_nuVtxX; // MCTruth nu vtx X
  double MC_nuVtxY; // MCTruth nu vtx Y
  double MC_nuVtxZ; // MCTruth nu vtx Z
  int MC_nNeutron; // Number of neutron(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_belowTH; // Number of proton(s) (p<200) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_middle; // Number of proton(s) (200<p<300) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nProton_aboveTH; // Number of proton(s) (p > 300) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPi0; // Number of pi0(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiPlus; // Number of pi plus(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int MC_nPiMinus; // Number of pi minus(s) from MCParticles, neutrino interaction + FSI for cc events (NC: default value)
  int Genie_nNeutron_preFSI;// before FSI 
  int Genie_nProton_preFSI;// before FSI 
  int Genie_nPi0_preFSI;// before FSI 
  int Genie_nPiPlus_preFSI;// before FSI 
  int Genie_nPiMinus_preFSI;// before FSI 

  bool if_selected; // If selected based on the reco info
  bool if_matchMu; // If the selected track matched with true muon from numu cc
  bool if_cosmic; // Check if a track is cosmic or not by if it has an associated MCParticle

  std::vector<double> true_mom_mu;//True momentum of muon track in the every event
  //std::vector<double> true_vtx_x;//True vertex of muon track (X)
  //std::vector<double> true_vtx_y;//True vertex of muon track (Y)
  //std::vector<double> true_vtx_z;//True vertex of muon track (Z)
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

  std::vector<double> mom_bestMCS_mu;//MCS best momentum of muon track in the every event
  std::vector<double> mom_bestMCS_ll_mu;//Likelihood of MCS best momentum of muon track in the every event
  //std::vector<double> mom_fwdMCS_mu;//MCS forward momentum of muon track in the every event
  //std::vector<double> mom_fwdMCS_ll_mu;//Likelihood of MCS forward momentum of muon track in the every event
  //std::vector<double> mom_bwdMCS_mu;//MCS backward momentum of muon track in the every event
  //std::vector<double> mom_bwdMCS_ll_mu;//Likelihood of MCS backward momentum of muon track in the every event
 // std::vector<double> mom_bestMCS_SCEcorr_mu;//MCS best momentum with SCE correction of muon track in the every event
 // std::vector<double> mom_bestMCS_SCEcorr_ll_mu;//Likelihood of MCS best momentum with SCE correction of muon track in the every event
 // std::vector<double> mom_fwdMCS_SCEcorr_mu;//MCS forward momentum with SCE correction of muon track in the every event
 // std::vector<double> mom_fwdMCS_SCEcorr_ll_mu;//Likelihood of MCS forward momentum with SCE correction of muon track in the every event
 // std::vector<double> mom_bwdMCS_SCEcorr_mu;//MCS backward momentum with SCE correction of muon track in the every event
 // std::vector<double> mom_bwdMCS_SCEcorr_ll_mu;//Likelihood of MCS backward momentum with SCE correction of muon track in the every event
  std::vector<double> mom_Range_mu;//Range momentum of muon track in the every event
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
  //int ntrack;//number of tracks in this event
  //int nshower;//number of shower in this event
  int n_pfp_nuDaughters; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters

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

  double _min_track_len;

  ::trkf::TrackMomentumCalculator _trk_mom_calculator;
};


SingleMuon::SingleMuon(fhicl::ParameterSet const& pset)
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
}

void SingleMuon::analyze(art::Event const& evt)
{

  //// Get necessary handles
  // MC Truth
  art::Handle< std::vector<simb::MCTruth> > Handle_MCTruth;
  evt.getByLabel(m_generatorLabel, Handle_MCTruth);
  std::vector<art::Ptr<simb::MCTruth> > MCTruthCollection;
  art::fill_ptr_vector(MCTruthCollection, Handle_MCTruth);

  // Genie Truth
  art::Handle< std::vector<simb::GTruth> > Handle_GTruth;
  evt.getByLabel(m_generatorLabel, Handle_GTruth);
  std::vector<art::Ptr<simb::GTruth> > GTruthCollection;
  art::fill_ptr_vector(GTruthCollection, Handle_GTruth);
 
  // MC Particle
  art::Handle< std::vector<simb::MCParticle> > Handle_MCParticle;
  evt.getByLabel(m_geantLabel, Handle_MCParticle);
  std::vector<art::Ptr<simb::MCParticle> > MCParticleCollection;
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

  //MCS momentum
  art::Handle<std::vector<recob::MCSFitResult> > Handle_mcsfitresult_mu;
  evt.getByLabel(m_MCSmuProducerLabel, Handle_mcsfitresult_mu);
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_v;
  art::fill_ptr_vector(mcsfitresult_mu_v, Handle_mcsfitresult_mu);

//  //MCS corrected momentum
//  art::Handle<std::vector<recob::MCSFitResult> > Handle_mcsfitresult_SCEcorr_mu;
//  evt.getByLabel(m_SCEcorr_MCSmuProducerLabel, Handle_mcsfitresult_SCEcorr_mu);
//  std::vector<art::Ptr<recob::MCSFitResult>> SCEcorr_mcsfitresult_mu_v;
//  art::fill_ptr_vector(SCEcorr_mcsfitresult_mu_v, Handle_mcsfitresult_SCEcorr_mu);

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

  //const std::vector< art::Ptr<recob::PFParticle> > &slc_PFP_v(slc_PFP_assn_v.at(SliceCollection.at(slc_id).key()));
  //art::FindManyP<recob::Slice> slc_PFP_assn_v(rawHandle_PFParticle, evt, data_label_assoSlc_);

  //std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(this_track_ptr.key());
  //BackTrackerTruthMatch backtrackertruthmatch;
  //backtrackertruthmatch.MatchToMCParticle(hit_handle,evt,trk_hits_ptrs);

  //RecoTruth recotruth;
  //recotruth.TrackToMCParticle(evt, Handle_Hit, m_trackProducerLabel);

  //ntrack = AllTrackCollection.size();
  //nshower = AllShowerCollection.size();

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

  //------- Get part of the generator neutrino info
  //Initiate the variables
  MC_beamNeutrino = false;
  MC_nupdg = -99999;
  MC_ccnc = -99999;
  MC_FV = false;

  MC_nNeutron = 0;
  MC_nProton_belowTH = 0;
  MC_nProton_middle = 0;
  MC_nProton_aboveTH = 0;
  MC_nPi0 = 0;
  MC_nPiPlus = 0;
  MC_nPiMinus = 0;

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
  if (MC_ccnc == 0){
    for(int i_mcp = 0; i_mcp < (int) MCParticleCollection.size(); i_mcp++){
      if(MCParticleCollection[i_mcp]->Process() == "primary"){
        // neutron
        if(MCParticleCollection[i_mcp]->PdgCode() == 2112) MC_nNeutron++;
        // proton
        if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() < 0.2) MC_nProton_belowTH++;
        if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() >= 0.2 && MCParticleCollection[i_mcp]->P() <= 0.3) MC_nProton_middle++;
        if(MCParticleCollection[i_mcp]->PdgCode() == 2212 && MCParticleCollection[i_mcp]->P() > 0.3) MC_nProton_aboveTH++;
        // pion0
        if(MCParticleCollection[i_mcp]->PdgCode() == 111) MC_nPi0++;
        // pion+
        if(MCParticleCollection[i_mcp]->PdgCode() == 211) MC_nPiPlus++;
        // pion-
        if(MCParticleCollection[i_mcp]->PdgCode() == -211) MC_nPiMinus++;
      }
    }
  }
  // Get Genie info on how many particles produced
  for(int i_gn = 0; i_gn < (int) GTruthCollection.size(); i_gn++){
    Genie_nNeutron_preFSI = GTruthCollection[i_gn]->fNumNeutron;
    Genie_nProton_preFSI = GTruthCollection[i_gn]->fNumProton;
    Genie_nPi0_preFSI = GTruthCollection[i_gn]->fNumPi0;
    Genie_nPiPlus_preFSI = GTruthCollection[i_gn]->fNumPiPlus;
    Genie_nPiMinus_preFSI = GTruthCollection[i_gn]->fNumPiMinus;
  }

  //-------- Get Reco neutrino (pfparticle)
  for(int i = 0; i < (int) pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      n_pfp_nuDaughters = pfp->NumDaughters();
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
      if_cosmic = true;
      if_matchMu = false;
      if_selected = false;

      //Todo: temperary version
      // Selection and Fill in Info
      if(n_dau_tracks == 1){
        //-- Fill RECO track info (in the naive version this is selected)
        if_selected = true;

        bool trk_contained = _fiducial_volume.InFV(daughter_Tracks.front()->Vertex<TVector3>(), daughter_Tracks.front()->End<TVector3>());
        trk_ifcontained.push_back(trk_contained);       
       
        bool vtx_contained = _fiducial_volume.InFV(daughter_Tracks.front()->Vertex<TVector3>());
        vtx_FV.push_back(vtx_contained);
 
        double bestMCS =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bestMomentum();
        double bestMCSLL =  mcsfitresult_mu_v.at(daughter_Tracks.front().key())->bestLogLikelihood();
        mom_bestMCS_mu.push_back(bestMCS);
        mom_bestMCS_ll_mu.push_back(bestMCSLL);

        auto assoCal = trackToCalAsso.at(daughter_Tracks.front().key());
        double Trk_Length = assoCal.front()->Range();  //It is said track length in pandoracali has spatial correction       
        //double trk_len = daughter_Tracks.front()->Length(); //TODO: The track length may be uncorrected. Get the corrected track length from pandoracali
        trk_length.push_back(Trk_Length);
        
        vtx_x.push_back(daughter_Tracks.front()->Vertex().X());
        vtx_y.push_back(daughter_Tracks.front()->Vertex().Y());
        vtx_z.push_back(daughter_Tracks.front()->Vertex().Z());
        start_x.push_back(daughter_Tracks.front()->Start().X());
        start_y.push_back(daughter_Tracks.front()->Start().Y());
        start_z.push_back(daughter_Tracks.front()->Start().Z());
        end_x.push_back(daughter_Tracks.front()->End().X());
        end_y.push_back(daughter_Tracks.front()->End().Y());
        end_z.push_back(daughter_Tracks.front()->End().Z());
        
        trk_phi.push_back(daughter_Tracks.front()->Phi());
        trk_theta.push_back(daughter_Tracks.front()->Theta());
        trk_costheta.push_back(cos(daughter_Tracks.front()->Theta()));
    
        double RangeMom = _trk_mom_calculator.GetTrackMomentum(Trk_Length, 13);
        mom_Range_mu.push_back(RangeMom);


        //-- Fill TRUE info, if the track is from numu cc muon
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
          if(MCparticle->PdgCode() == 13){
            if_matchMu = true;

            auto TrueTrackPos = MCparticle->EndPosition() - MCparticle->Position();
            true_mom_mu.push_back(MCparticle->P());
            true_start_x.push_back(MCparticle->Position().X());
            true_start_y.push_back(MCparticle->Position().Y());
            true_start_z.push_back(MCparticle->Position().Z());
            true_end_x.push_back(MCparticle->EndPosition().X());
            true_end_y.push_back(MCparticle->EndPosition().Y());
            true_end_z.push_back(MCparticle->EndPosition().Z());
            true_trk_phi.push_back(TrueTrackPos.Phi());
            true_trk_theta.push_back(TrueTrackPos.Theta());
            true_trk_costheta.push_back(cos(TrueTrackPos.Theta()));
            true_trk_length.push_back(sqrt(TrueTrackPos.X()*TrueTrackPos.X() + TrueTrackPos.Y()*TrueTrackPos.Y() + TrueTrackPos.Z()*TrueTrackPos.Z()));
            trk_pdg.push_back(MCparticle->PdgCode());
          }
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
      if(n_dau_tracks==1){
        // Get track length, and do we use it later?
        ////// Add SCE correction?
        vTrk_len.push_back(daughter_Tracks.front()->Length());
        // Get directional info to see if the track start is in the TPC / FV by dE/dx
        // pandoracaliSCE has E-field and spatial correction
        auto assoCal = trackToCalAsso.at(daughter_Tracks.front().key()); // take the only track
        std::cout<<"calo size: "<< assoCal.size()<<std::endl;
        if(assoCal.size()!=3){
          throw cet::exception("[Numu0pi0p]") << "Where are the three planes for the calorimetry!" << std::endl;
        }
        int id_collection = -99;
        for(int id_pl = 0; id_pl < 3; id_pl++){
          if(assoCal[id_pl]->PlaneID().Plane == 2){
            id_collection = id_pl;
          }
        }
        std::cout<<"collection plane id: "<<id_collection<<std::endl;
       // if(assoCal.empty()){
       //   throw cet::exception("[Numu0pi0p]") << "The only track don't have corresponding calorimetry." << std::endl;
       // }
       // if(assoCal.size() > 1){
       //   throw cet::exception("[Numu0pi0p]") << "The only track associates more than a set of calorimetry." << std::endl;
       // }
        std::vector<float> vdEdx = assoCal.front()->dEdx();
        std::vector<float> vdQdx = assoCal.front()->dQdx();
        std::vector<float> vresRange = assoCal.front()->ResidualRange();
        float Trk_Length = assoCal.front()->Range();  //It is said track length in pandoracali has spatial correction

        int size_dEdx = vdEdx.size();
        int unit = size_dEdx / 4;
        // Average 4 hits and compare 4 groups of average dEdx
        std::vector<float> v_dEdx_seg;
        if(size_dEdx > 16){
          for(int i_seg = 0; i_seg < 4; i_seg++){

            v_dEdx_seg.push_back(0.25 * (vdEdx[size_dEdx - unit * i_seg - 1]
                                       + vdEdx[size_dEdx - unit * i_seg - 2]
                                       + vdEdx[size_dEdx - unit * i_seg - 3]
                                       + vdEdx[size_dEdx - unit * i_seg - 4]));
          }
          auto sorted_v_dEdx_seg = v_dEdx_seg;
          sort(sorted_v_dEdx_seg.begin(), sorted_v_dEdx_seg.end());

          //If the direction of the track is correct, dEdx_seg 1 should represent track end, so dEdx_seg1 > dEdx_seg4
          //std::cout<<"seg1: "<<dEdx_seg1<<", seg2: "<<dEdx_seg2<<", seg3: "<<dEdx_seg3<<", seg4: "<<dEdx_seg4<<std::endl;
          //if() 
        }
        std::cout<<"Spatial corrected length: "<< Trk_Length<<", raw track length: "<<daughter_Tracks.front()->Length()<<std::endl;
        //std::cout<<"size of Range: "<< vresRange.size()<<", size of dEdx: "<<vdEdx.size()<<std::endl;
        //for(int kk = 0; kk < size_dEdx; kk++){
        //  std::cout<<"Residual: "<< vresRange[kk]<<", dEdx: "<< vdEdx[kk]<<std::endl;
        //}
      }


    }
  }

/*
  //
  for(int trk_id = 0; trk_id < ntrack; trk_id++){
    // Back track matching
    std::vector<art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(AllTrackCollection[trk_id].key());
    BackTrackerTruthMatch backtrackertruthmatch;
    backtrackertruthmatch.MatchToMCParticle(Handle_Hit,evt,trk_hits_ptrs);
    auto MCparticle = backtrackertruthmatch.ReturnMCParticle();
    if(!MCparticle){
      //std::cout<<"MC particle does not exist!"<<std::endl;
    }
    else{
     
      auto TrueTrackPos = MCparticle->EndPosition() - MCparticle->Position();
      true_mom_mu.push_back(MCparticle->P());
      true_vtx_x.push_back(MCparticle->Vx());
      true_vtx_y.push_back(MCparticle->Vy());
      true_vtx_z.push_back(MCparticle->Vz());
      true_start_x.push_back(MCparticle->Position().X());
      true_start_y.push_back(MCparticle->Position().Y());
      true_start_z.push_back(MCparticle->Position().Z());
      true_end_x.push_back(MCparticle->EndPosition().X());
      true_end_y.push_back(MCparticle->EndPosition().Y());
      true_end_z.push_back(MCparticle->EndPosition().Z());
      true_trk_phi.push_back(TrueTrackPos.Phi());
      true_trk_theta.push_back(TrueTrackPos.Theta());
      true_trk_length.push_back(sqrt(TrueTrackPos.X()*TrueTrackPos.X() + TrueTrackPos.Y()*TrueTrackPos.Y() + TrueTrackPos.Z()*TrueTrackPos.Z()));
      trk_pdg.push_back(MCparticle->PdgCode());
    }
    // Check if the track is contained or not
    bool contained = _fiducial_volume.InFV(AllTrackCollection[trk_id]->Vertex<TVector3>(), AllTrackCollection[trk_id]->End<TVector3>());
    trk_ifcontained.push_back(contained);

    // Fill MCS fitting result
    double bestMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestMomentum();
    double bestMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestLogLikelihood();
    mom_bestMCS_mu.push_back(bestMCS);
    mom_bestMCS_ll_mu.push_back(bestMCSLL);
    double fwdMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdMomentum();
    double fwdMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdLogLikelihood();
    mom_fwdMCS_mu.push_back(fwdMCS);
    mom_fwdMCS_ll_mu.push_back(fwdMCSLL);
    double bwdMCS =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdMomentum();
    double bwdMCSLL =  mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdLogLikelihood();
    mom_bwdMCS_mu.push_back(bwdMCS);
    mom_bwdMCS_ll_mu.push_back(bwdMCSLL);
   
    //double bestMCS_SCEcorr =  SCEcorr_mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestMomentum();
    //double bestMCSLL_SCEcorr =  SCEcorr_mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bestLogLikelihood();
    //mom_bestMCS_SCEcorr_mu.push_back(bestMCS_SCEcorr);
    //mom_bestMCS_SCEcorr_ll_mu.push_back(bestMCSLL_SCEcorr);
    //double fwdMCS_SCEcorr =  SCEcorr_mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdMomentum();
    //double fwdMCSLL_SCEcorr =  SCEcorr_mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->fwdLogLikelihood();
    //mom_fwdMCS_SCEcorr_mu.push_back(fwdMCS_SCEcorr);
    //mom_fwdMCS_SCEcorr_ll_mu.push_back(fwdMCSLL_SCEcorr);
    //double bwdMCS_SCEcorr =  SCEcorr_mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdMomentum();
    //double bwdMCSLL_SCEcorr =  SCEcorr_mcsfitresult_mu_v.at(AllTrackCollection[trk_id].key())->bwdLogLikelihood();
    //mom_bwdMCS_SCEcorr_mu.push_back(bwdMCS_SCEcorr);
    //mom_bwdMCS_SCEcorr_ll_mu.push_back(bwdMCSLL_SCEcorr);

    // Track Length
    double trk_len = AllTrackCollection[trk_id]->Length();
    trk_length.push_back(trk_len);
    // Track vertex
    vtx_x.push_back(AllTrackCollection[trk_id]->Vertex().X());
    vtx_y.push_back(AllTrackCollection[trk_id]->Vertex().Y());
    vtx_z.push_back(AllTrackCollection[trk_id]->Vertex().Z());
    start_x.push_back(AllTrackCollection[trk_id]->Start().X());
    start_y.push_back(AllTrackCollection[trk_id]->Start().Y());
    start_z.push_back(AllTrackCollection[trk_id]->Start().Z());
    end_x.push_back(AllTrackCollection[trk_id]->End().X());
    end_y.push_back(AllTrackCollection[trk_id]->End().Y());
    end_z.push_back(AllTrackCollection[trk_id]->End().Z());

    // Track angles
    trk_phi.push_back(AllTrackCollection[trk_id]->Phi());
    trk_theta.push_back(AllTrackCollection[trk_id]->Theta());

    // Range based momentum (Give muon pdg = 13), track length in [cm] ?
    double RangeMom = _trk_mom_calculator.GetTrackMomentum(trk_len, 13);
    mom_Range_mu.push_back(RangeMom);
  }
*/


  my_event_->Fill();

  true_mom_mu.clear();
  //true_vtx_x.clear();
  //true_vtx_y.clear();
  //true_vtx_z.clear();
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

  mom_bestMCS_mu.clear();
  mom_bestMCS_ll_mu.clear();
  //mom_fwdMCS_mu.clear();
  //mom_fwdMCS_ll_mu.clear();
  //mom_bwdMCS_mu.clear();
  //mom_bwdMCS_ll_mu.clear();
  //mom_bestMCS_SCEcorr_mu.clear();
  //mom_bestMCS_SCEcorr_ll_mu.clear();
  //mom_fwdMCS_SCEcorr_mu.clear();
  //mom_fwdMCS_SCEcorr_ll_mu.clear();
  //mom_bwdMCS_SCEcorr_mu.clear();
  //mom_bwdMCS_SCEcorr_ll_mu.clear();
  mom_Range_mu.clear();
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
}

void SingleMuon::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;
  // Make a tree to store momentum information
  my_event_ = tfs->make<TTree>("tree","tree");

  my_event_->Branch("MC_beamNeutrino", &MC_beamNeutrino);
  my_event_->Branch("MC_nupdg", &MC_nupdg);
  my_event_->Branch("MC_ccnc", &MC_ccnc);
  my_event_->Branch("MC_nuVtxX", &MC_nuVtxX);
  my_event_->Branch("MC_nuVtxY", &MC_nuVtxY);
  my_event_->Branch("MC_nuVtxZ", &MC_nuVtxZ);
  my_event_->Branch("MC_FV", &MC_FV);
  my_event_->Branch("MC_nNeutron", &MC_nNeutron);
  my_event_->Branch("MC_nProton_belowTH", &MC_nProton_belowTH);
  my_event_->Branch("MC_nProton_middle", &MC_nProton_middle);
  my_event_->Branch("MC_nProton_aboveTH", &MC_nProton_aboveTH);
  my_event_->Branch("MC_nPi0", &MC_nPi0);
  my_event_->Branch("MC_nPiPlus", &MC_nPiPlus);
  my_event_->Branch("MC_nPiMinus", &MC_nPiMinus);
  my_event_->Branch("Genie_nNeutron_preFSI", &Genie_nNeutron_preFSI);
  my_event_->Branch("Genie_nProton_preFSI", &Genie_nProton_preFSI);
  my_event_->Branch("Genie_nPi0_preFSI", &Genie_nPi0_preFSI);
  my_event_->Branch("Genie_nPiPlus_preFSI", &Genie_nPiPlus_preFSI);
  my_event_->Branch("Genie_nPiMinus_preFSI", &Genie_nPiMinus_preFSI);

  my_event_->Branch("true_mom_mu", &true_mom_mu);
  //my_event_->Branch("true_vtx_x", &true_vtx_x);
  //my_event_->Branch("true_vtx_y", &true_vtx_y);
  //my_event_->Branch("true_vtx_z", &true_vtx_z);
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

  my_event_->Branch("n_dau_tracks", &n_dau_tracks);
  my_event_->Branch("n_dau_showers", &n_dau_showers);

  my_event_->Branch("if_cosmic", &if_cosmic);
  my_event_->Branch("if_matchMu", &if_matchMu);
  my_event_->Branch("if_selected", &if_selected);

  my_event_->Branch("mom_bestMCS_mu", &mom_bestMCS_mu);
  my_event_->Branch("mom_bestMCS_ll_mu", &mom_bestMCS_ll_mu);
  //my_event_->Branch("mom_fwdMCS_mu", &mom_fwdMCS_mu);
  //my_event_->Branch("mom_fwdMCS_ll_mu", &mom_fwdMCS_ll_mu);
  //my_event_->Branch("mom_bwdMCS_mu", &mom_bwdMCS_mu);
  //my_event_->Branch("mom_bwdMCS_ll_mu", &mom_bwdMCS_ll_mu);
  //my_event_->Branch("mom_bestMCS_SCEcorr_mu", &mom_bestMCS_SCEcorr_mu);
  //my_event_->Branch("mom_bestMCS_SCEcorr_ll_mu", &mom_bestMCS_SCEcorr_ll_mu);
  //my_event_->Branch("mom_fwdMCS_SCEcorr_mu", &mom_fwdMCS_SCEcorr_mu);
  //my_event_->Branch("mom_fwdMCS_SCEcorr_ll_mu", &mom_fwdMCS_SCEcorr_ll_mu);
  //my_event_->Branch("mom_bwdMCS_SCEcorr_mu", &mom_bwdMCS_SCEcorr_mu);
  //my_event_->Branch("mom_bwdMCS_SCEcorr_ll_mu", &mom_bwdMCS_SCEcorr_ll_mu);
  my_event_->Branch("mom_Range_mu", &mom_Range_mu);
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
  //my_event_->Branch("ntrack", &ntrack);
  //my_event_->Branch("nshower", &nshower);
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

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
 
  std::vector<double> true_mom_mu;//True momentum of muon track in the every event
  std::vector<double> true_vtx_x;//True vertex of muon track (X)
  std::vector<double> true_vtx_y;//True vertex of muon track (Y)
  std::vector<double> true_vtx_z;//True vertex of muon track (Z)
  std::vector<double> true_start_x;//True start of muon track (X)
  std::vector<double> true_start_y;//True start of muon track (Y)
  std::vector<double> true_start_z;//True start of muon track (Z)
  std::vector<double> true_end_x;//True end of muon track (X)
  std::vector<double> true_end_y;//True end of muon track (Y)
  std::vector<double> true_end_z;//True end of muon track (Z)
  std::vector<double> true_trk_phi;//True phi of muon track 
  std::vector<double> true_trk_theta;//True theta of muon track 
  std::vector<double> true_trk_length;//True track length (distance from the start to the end point) 
  std::vector<double> trk_pdg;//Track pdg 

  std::vector<double> mom_bestMCS_mu;//MCS best momentum of muon track in the every event
  std::vector<double> mom_bestMCS_ll_mu;//Likelihood of MCS best momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_mu;//MCS forward momentum of muon track in the every event
  std::vector<double> mom_fwdMCS_ll_mu;//Likelihood of MCS forward momentum of muon track in the every event
  std::vector<double> mom_bwdMCS_mu;//MCS backward momentum of muon track in the every event
  std::vector<double> mom_bwdMCS_ll_mu;//Likelihood of MCS backward momentum of muon track in the every event
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
  std::vector<double> trk_length;//Range momentum of muon track in the every event
  std::vector<bool> trk_ifcontained;//to check if the track is contained or not
  int ntrack;//number of tracks in this event
  int nshower;//number of shower in this event
  int n_pfp_nuDaughters; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters

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

  ntrack = AllTrackCollection.size();
  nshower = AllShowerCollection.size();

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

  //////// Get neutrino
  for(int i = 0; i < (int) pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      std::cout<<"number of neutrino daughters: "<<pfp->NumDaughters()<<std::endl;
      n_pfp_nuDaughters = pfp->NumDaughters();
      if(n_pfp_nuDaughters<4){
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

      ///////////Logics.......
      //
      //////1 track + 0 shower
      // What is the direct
      //
      //
      //
      //
      //
      //
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
        std::cout<<"size of Range: "<< vresRange.size()<<", size of dEdx: "<<vdEdx.size()<<std::endl;
        for(int kk = 0; kk < size_dEdx; kk++){
          std::cout<<"Residual: "<< vresRange[kk]<<", dEdx: "<< vdEdx[kk]<<std::endl;
        }
      }


      /////
      //Reco-True Matching
      /////
//////////////////////////////////////
      if(pfp->NumDaughters()==1){
        // Check if the pfp particle is track or not
        // Check if the pfp particle is muon or not
        //recob::PFParticle daughter = pfp->Daughters().front();
        std::cout<<"type: "<<typeid(pfp->Daughters()).name()<<std::endl;
        //std::cout<<"pdg of the only neutrino daughter: "<< daughter->PdgCode();
        //if(daughter->PdgCode() == 13){
        //}
        // Check if this pfp particle has daughter or not
      }
    }
  }
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

  my_event_->Fill();

  true_mom_mu.clear();
  true_vtx_x.clear();
  true_vtx_y.clear();
  true_vtx_z.clear();
  true_start_x.clear();
  true_start_y.clear();
  true_start_z.clear();
  true_end_x.clear();
  true_end_y.clear();
  true_end_z.clear();
  true_trk_phi.clear();
  true_trk_theta.clear();
  true_trk_length.clear();
  trk_pdg.clear();

  mom_bestMCS_mu.clear();
  mom_bestMCS_ll_mu.clear();
  mom_fwdMCS_mu.clear();
  mom_fwdMCS_ll_mu.clear();
  mom_bwdMCS_mu.clear();
  mom_bwdMCS_ll_mu.clear();
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
  trk_length.clear();
  trk_ifcontained.clear();
}

void SingleMuon::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;
  // Make a tree to store momentum information
  my_event_ = tfs->make<TTree>("tree","tree");

  my_event_->Branch("true_mom_mu", &true_mom_mu);
  my_event_->Branch("true_vtx_x", &true_vtx_x);
  my_event_->Branch("true_vtx_y", &true_vtx_y);
  my_event_->Branch("true_vtx_z", &true_vtx_z);
  my_event_->Branch("true_start_x", &true_start_x);
  my_event_->Branch("true_start_y", &true_start_y);
  my_event_->Branch("true_start_z", &true_start_z);
  my_event_->Branch("true_end_x", &true_end_x);
  my_event_->Branch("true_end_y", &true_end_y);
  my_event_->Branch("true_end_z", &true_end_z);
  my_event_->Branch("true_trk_phi", &true_trk_phi);
  my_event_->Branch("true_trk_theta", &true_trk_theta);
  my_event_->Branch("true_trk_length", &true_trk_length);
  my_event_->Branch("trk_pdg", &trk_pdg);

  my_event_->Branch("mom_bestMCS_mu", &mom_bestMCS_mu);
  my_event_->Branch("mom_bestMCS_ll_mu", &mom_bestMCS_ll_mu);
  my_event_->Branch("mom_fwdMCS_mu", &mom_fwdMCS_mu);
  my_event_->Branch("mom_fwdMCS_ll_mu", &mom_fwdMCS_ll_mu);
  my_event_->Branch("mom_bwdMCS_mu", &mom_bwdMCS_mu);
  my_event_->Branch("mom_bwdMCS_ll_mu", &mom_bwdMCS_ll_mu);
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
  my_event_->Branch("trk_length", &trk_length);
  my_event_->Branch("trk_ifcontained", &trk_ifcontained);
  my_event_->Branch("ntrack", &ntrack);
  my_event_->Branch("nshower", &nshower);
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

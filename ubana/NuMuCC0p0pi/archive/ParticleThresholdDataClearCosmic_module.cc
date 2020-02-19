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
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

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
  lar_pandora::LArPandoraHelper larpandora;

  TTree * POTtree;
  int run, subrun;
  double POT_miss1E10;

  TTree * my_event_;
  
  void Initialize_event();
  
  std::vector<double> start_x;//Reconstructed start x in the every event
  std::vector<double> start_y;//Reconstructed start y in the every event
  std::vector<double> start_z;//Reconstructed start z in the every event
  std::vector<double> end_x;//Reconstructed end x in the every event
  std::vector<double> end_y;//Reconstructed end y in the every event
  std::vector<double> end_z;//Reconstructed end z in the every event
  std::vector<double> trk_phi;//Reconstructed track phi in the every event
  std::vector<double> trk_theta;//Reconstructed track theta in the every event
  std::vector<double> trk_costheta;//Reconstructed track cos(theta) in the every event
  std::vector<double> trk_length_noSCE;//Range momentum of muon track in the every event no spatial correction
  std::vector<bool> trk_ifcontained;//to check if the track is contained or not
  std::vector<bool> vtx_FV;//to check if the vertex is in FV or not
  int NCosmic; // Number of clear cosmic in every event

  std::string                         m_pandoraLabel;
  std::string                         m_MetaProducerLabel;
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
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_MetaProducerLabel(pset.get<std::string>("PFPMetaProducerLabel")),
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

  //PFP Meta Data collection
  lar_pandora::PFParticleVector pfParticles;
  lar_pandora::PFParticlesToMetadata pfParticlesToMetadata;
  larpandora.CollectPFParticleMetadata(evt, m_pandoraLabel, pfParticles, pfParticlesToMetadata);

  //pfp track association
  art::FindManyP<recob::Track> pfpToTrackAsso(Handle_pfParticle, evt, m_trackProducerLabel);

  //pfp shower association
  art::FindManyP<recob::Shower> pfpToShowerAsso(Handle_pfParticle, evt, m_showerProducerLabel);

  // pfp pfpMetaData association (to see if a track is a Clear Cosmic)
  art::FindMany<larpandoraobj::PFParticleMetadata> pfpToMetaAsso(Handle_pfParticle, evt, m_MetaProducerLabel); 
 
  //PID
  art::FindManyP<anab::ParticleID> PIDTotrackAsso(Handle_TPCtrack,evt,PID_TrackAssLabel);

  //------Loop over PFParticles
  for(int pfp_id = 0; pfp_id < (int) pfParticle_v.size(); pfp_id++){
    auto pfp = pfParticle_v[pfp_id];
    if(pfp->IsPrimary()){
      // If pfp corresponds to a clear cosmic
      //auto MetaData = pfpToMetaAsso.at(pfp.key());
    
      const auto itParticle = pfParticlesToMetadata.find(pfp);
      if (itParticle == pfParticlesToMetadata.end()){
        throw cet::exception("WorkshopTrackShowerHelper") << "PFParticle has no metadata" << std::endl;
      }

      if (itParticle->second.size() != 1){
        throw cet::exception("WorkshopTrackShowerHelper") << "PFParticle has mutiple metadata" << std::endl;
      }
     
      const auto metadata = itParticle->second.front();
      const auto propertiesMap = metadata->GetPropertiesMap();
      const auto isClearCosmic = propertiesMap.find("IsClearCosmic");
 
      //if (isClearCosmic->second == 1){
       // std::cout<<"found one clear cosmic"<<std::endl;
      //}

      //if(MetaData.size() == 1 && MetaData.front()->IsClearCosmic() == 1){
      if(isClearCosmic->second == 1){
        std::cout<<"Clear Cosmic!"<<std::endl;
        // Track assoiciated to the clear cosmic
        auto assoTrack = pfpToTrackAsso.at(pfp.key());
        
        TVector3 Trk_start = assoTrack.front()->Vertex<TVector3>();
        auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
        TVector3 Trk_start_SCEcorr;
        Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
        Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
        Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

        TVector3 Trk_end = assoTrack.front()->End<TVector3>();
        auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
        TVector3 Trk_end_SCEcorr;
        Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X());
        Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
        Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());

        bool trk_contained = _fiducial_volume.InFV(Trk_start_SCEcorr, Trk_end_SCEcorr);
        trk_ifcontained.push_back(trk_contained);

        bool vtx_InFV = _fiducial_volume.InFV(Trk_start_SCEcorr);
        vtx_FV.push_back(vtx_InFV);

        start_x.push_back(Trk_start_SCEcorr.X());
        start_y.push_back(Trk_start_SCEcorr.Y());
        start_z.push_back(Trk_start_SCEcorr.Z());
        end_x.push_back(Trk_end_SCEcorr.X());
        end_y.push_back(Trk_end_SCEcorr.Y());
        end_z.push_back(Trk_end_SCEcorr.Z());

        trk_phi.push_back(assoTrack.front()->Phi());
        trk_theta.push_back(assoTrack.front()->Theta());
        trk_costheta.push_back(cos(assoTrack.front()->Theta()));
      } 
    }
  }

  NCosmic = trk_theta.size();
  my_event_->Fill();
  start_x.clear();
  start_y.clear();
  start_z.clear();
  end_x.clear();
  end_y.clear();
  end_z.clear();
  trk_phi.clear();
  trk_theta.clear();
  trk_costheta.clear();
  trk_length_noSCE.clear();
  trk_ifcontained.clear();
  vtx_FV.clear(); 
}

void ParticleThreshold::Initialize_event()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;

  // Make a tree to store selection information
  my_event_ = tfs->make<TTree>("tree","tree");

  my_event_->Branch("start_x", &start_x);
  my_event_->Branch("start_y", &start_y);
  my_event_->Branch("start_z", &start_z);
  my_event_->Branch("end_x", &end_x);
  my_event_->Branch("end_y", &end_y);
  my_event_->Branch("end_z", &end_z);
  my_event_->Branch("trk_phi", &trk_phi);
  my_event_->Branch("trk_theta", &trk_theta);
  my_event_->Branch("trk_costheta", &trk_costheta);
  my_event_->Branch("trk_length_noSCE", &trk_length_noSCE);
  my_event_->Branch("trk_ifcontained", &trk_ifcontained);
  my_event_->Branch("vtx_FV", &vtx_FV);
  my_event_->Branch("NCosmic", &NCosmic);

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

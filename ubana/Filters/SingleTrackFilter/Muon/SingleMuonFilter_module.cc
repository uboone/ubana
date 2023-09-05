////////////////////////////////////////////////////////////////////////
// Class:       SingleMuonFilter
// Plugin Type: filter (art v3_01_02)
// File:        SingleMuonFilter_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
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
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TTree.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"

#include "ubana/Utilities/FiducialVolume.h"

class SingleMuonFilter : public art::EDFilter {
public:
  explicit SingleMuonFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SingleMuonFilter(SingleMuonFilter const&) = delete;
  SingleMuonFilter(SingleMuonFilter&&) = delete;
  SingleMuonFilter& operator=(SingleMuonFilter const&) = delete;
  SingleMuonFilter& operator=(SingleMuonFilter&&) = delete;

  // Required functions.
  //bool filter(art::Event const& evt) override;
  bool filter(art::Event & evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<art::TFileService> tfs;

  ::ubana::FiducialVolume _fiducial_volume;

  spacecharge::SpaceCharge const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  TTree * POTtree;
  int run, subrun;
  double POT_miss1E10;

  TTree * my_event_;
  
  void Initialize_event();

  int n_pfp_nuDaughters; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters


  std::string                         m_pandoraLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_showerProducerLabel;

};


SingleMuonFilter::SingleMuonFilter(fhicl::ParameterSet const& pset): 
  EDFilter{pset},
  _fiducial_volume(pset.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                   geo->DetHalfHeight(),
                   2.*geo->DetHalfWidth(),
                   geo->DetLength()),
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  _fiducial_volume.PrintConfig();
}

//bool SingleMuonFilter::filter(art::Event const& evt)
bool SingleMuonFilter::filter(art::Event & evt)
{

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

  //PFParticle collection
  art::Handle< std::vector<recob::PFParticle> > Handle_pfParticle;
  evt.getByLabel(m_pandoraLabel, Handle_pfParticle);
  std::vector<art::Ptr<recob::PFParticle>> pfParticle_v;
  art::fill_ptr_vector(pfParticle_v, Handle_pfParticle);

  //pfp track association
  art::FindManyP<recob::Track> pfpToTrackAsso(Handle_pfParticle, evt, m_trackProducerLabel);

  //pfp shower association
  art::FindManyP<recob::Shower> pfpToShowerAsso(Handle_pfParticle, evt, m_showerProducerLabel);
  
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
  std::vector<float> vTrk_len;

  TVector3 Trk_start;
  TVector3 Trk_end;
  TVector3 Trk_start_SCEcorr; 
  TVector3 Trk_end_SCEcorr;
 

  // Reject an event if
  // - No neutrino slice
  // - Not 1-track 0-shower topology
  // - Vtx not in FV
 
  //bool if_Nu_Slice = false;
  //-------- Get Reco neutrino (pfparticle)
  for(unsigned int i = 0; i < pfParticle_v.size(); i++){
    auto pfp = pfParticle_v[i];
    if(pfp->IsPrimary() && pfp->PdgCode() == 14){
      //if_Nu_Slice = true;

      n_pfp_nuDaughters = pfp->NumDaughters();

      // For CC0pi0p, we only consider the case with the number of neutrino daughters less than 4
      if(n_pfp_nuDaughters < 4){
        // Get the pointer for the daughters of the neutrino
        for(int j = 0; j < n_pfp_nuDaughters; j++){
          auto Iterator = pfParticleIdMap.find(pfp->Daughters().at(j));
          auto dau_pfp = Iterator->second;
          NeutrinoDaughters.push_back(dau_pfp);
          // Collect pfparticle associated track in a vector
          auto assoTrack = pfpToTrackAsso.at(dau_pfp.key()); // vector
          if(assoTrack.size()==1){
            daughter_Tracks.push_back(assoTrack.front());
          }
          if(assoTrack.size()>1){
            throw cet::exception("[Numu0pi0p]") << "PFParticle has >1 track!" << std::endl;
          }
          // Collect pfparticle associated shower in a vector
          auto assoShower = pfpToShowerAsso.at(dau_pfp.key()); // vector
          if(assoShower.size()==1){
            daughter_Showers.push_back(assoShower.front());
          }
          if(assoShower.size()>1){
            throw cet::exception("[Numu0pi0p]") << "PFParticle has >1 shower!" << std::endl;
          }
        } // finish looping of pfp
      }
     
      //number of tracks and showers
      n_dau_tracks = daughter_Tracks.size();
      n_dau_showers = daughter_Showers.size();

      // Selection and Fill in Info
      if(n_dau_tracks == 1 && n_dau_showers == 0){
        // Add spatial correction to the track start and end
        Trk_start = daughter_Tracks.front()->Vertex<TVector3>();
        auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()), 0);
        Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
        Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
        Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

        bool vtx_InFV = _fiducial_volume.InFV(Trk_start_SCEcorr);

        if(vtx_InFV){
          return true;
        } // if in FV
        else{
          return false;
        } //not in FV
      } // if 1trk, 0shw
      else{
        return false;        
      } // not 1trk, 0 shw
      
    } // if pfp < 4
  } // Finish loop all the pfp particles
  
  return false; // If the function hasn't return before this, it means that there's no neutrino slice

}

void SingleMuonFilter::beginJob()
{
  // Implementation of optional member function here.
}

void SingleMuonFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SingleMuonFilter)

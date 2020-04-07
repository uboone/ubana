//////////////////////////////////////////////////////////////////////
// Class:       SingleMuonFilter
// Plugin Type: filter (art v3_01_02)
// File:        ThreeTrackFilter_module.cc
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
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"

#include "FiducialVolume.h"

class ThreeTrackFilter;


class ThreeTrackFilter : public art::EDFilter {
public:
  explicit ThreeTrackFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ThreeTrackFilter(ThreeTrackFilter const&) = delete;
  ThreeTrackFilter(ThreeTrackFilter&&) = delete;
  ThreeTrackFilter& operator=(ThreeTrackFilter const&) = delete;
  ThreeTrackFilter& operator=(ThreeTrackFilter&&) = delete;

  // Required functions.
  //bool filter(art::Event const& evt) override;
  bool filter(art::Event & evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void ClearLocalData();

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

  int n_pfp_nuDaughters; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters
  int n_good_daughters; //number of daughters that are tracks and whose start is within 5cm of the vertex
  int ntracks;
  int nshowers;
  float vtxcut;
  float reco_nu_vtxx = -9999., reco_nu_vtxy=-9999., reco_nu_vtxz=-9999.;
  bool _debug = false;
  TVector3 trk_good_start;
  TVector3 trk_good_start_SCEcorr;
  TVector3 trk_210_SCEcorr;

  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles,pfparticles1;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::VertexVector pfvertices;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints;

  std::string                         m_pandoraLabel;
  std::string                         m_trackProducerLabel;
  std::string                         m_showerProducerLabel;
  std::string                         m_pfp_producer; //pfp producer   


};


ThreeTrackFilter::ThreeTrackFilter(fhicl::ParameterSet const& pset): 
  EDFilter{pset},
  m_pandoraLabel(pset.get<std::string>("PandoraLabel")),
  m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
  m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel")),
  m_pfp_producer(pset.get<std::string>("PandoraLabel"))

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _fiducial_volume.Configure(pset.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  ntracks = pset.get<int>("NTrack", 3); //fhicl parameter that makes ther required number of tracks == 3
  nshowers = pset.get<int>("NShowers", 0); //"" number of showers == 0
  vtxcut =  pset.get<float>("VtxCut", 5.0); //require the start of each track to be within 5cm of the reconstructed vertex
  _debug = pset.get<bool>("Debug", "true");

}

//bool ThreeTrackFilter::filter(art::Event const& evt)
bool ThreeTrackFilter::filter(art::Event & evt)
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

  //PFP track association
  art::FindManyP<recob::Track> pfpToTrackAsso(Handle_pfParticle, evt, m_trackProducerLabel);

  //PFP shower association
  art::FindManyP<recob::Shower> pfpToShowerAsso(Handle_pfParticle, evt, m_showerProducerLabel);

  //adding a boatload of stuff to get reco vertex
  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  larpandora.CollectPFParticles(evt, m_pfp_producer, pfparticles1, particlesToSpacePoints);
  if(_debug) std::cout<<"[Numu0pi2p] Number of PFParticles in the Event: "<<pfparticles.size()<<std::endl;

  //Stuff in order to get the reconstructed vertex:
  larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
  larpandora.CollectVertices(evt, m_pandoraLabel, pfvertices, particlesToVertices);
  if(_debug) std::cout<<"[Numu0pi2p] Number of RecoNeutrinos: "<<pfneutrinos.size()<<std::endl;
  if(_debug) std::cout<<"[Numu0pi2p] Number of RecoVertices: "<<pfvertices.size() <<std::endl;

  // Get mapping from ID to PFParticle
  std::unordered_map<size_t, art::Ptr<recob::PFParticle> > pfParticleIdMap;
  for (unsigned int i = 0; i < Handle_pfParticle->size(); ++i){
    auto pfp = pfParticle_v[i];
    if (!pfParticleIdMap.emplace(pfp->Self(), pfp).second){
      throw cet::exception("[Numu0pi0p]")<< "Repeated PFParticles!" << std::endl;
    }
  }

  //Define the daughters of neutrino
  std::vector<art::Ptr<recob::PFParticle> > NeutrinoDaughters; //total number of neutrino daughters
  std::vector<art::Ptr<recob::Track> > daughter_Tracks; //number of daugher tracks
  std::vector<art::Ptr<recob::Track> > good_daughter_Tracks; //Daughters that happen to be close to the vertex (i.e. within 5cm)
  std::vector<art::Ptr<recob::Shower> > daughter_Showers; //number of daughter showers
  std::vector<float> vTrk_len;

  ///////////////////////////////////////////////////////////
  // Reject an event if 1) Not in the neutrino slice 2) Not 3-track 0-shower topology 3) Vtx not in FV 4) the 3 track starts are not within 5cm of the reconstructed vertex
  //////////////////////////////////////////////////////////
  
  if (pfparticles.size() == 0){
    if(_debug) std::cout<<"[Numu0pi2p] No Reconstructed PFParticles in the event." << std::endl;
  }
  else {//if there are reconstructed pfparticles continue on your merry way       
    if(_debug) std::cout<< "[Numu0pi2p] Number of neutrino candidates in the event: " <<pfneutrinos.size() <<std::endl;
    if(pfneutrinos.size()==1){ //when there is exactly 1 neutrino candidate in the event                                                                               
	  
      //////////////////////////////////////
      //Get Reco neutrino (pfparticle)
      /////////////////////////////////////
      for(unsigned int i = 0; i < pfParticle_v.size(); i++){
	auto pfp = pfParticle_v[i];
	if(pfp->IsPrimary() && pfp->PdgCode() == 14){
	  //if_Nu_Slice = true;

	  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfp);
	  const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
	  reco_nu_vtxx = neutrino_vtx.X();
	  reco_nu_vtxy = neutrino_vtx.Y();
	  reco_nu_vtxz = neutrino_vtx.Z();
	
	  if(_debug) std::cout<<"[Numu0pi2p] Location of the Reco Vertex: "<<reco_nu_vtxx<<" <-Vx "<<reco_nu_vtxy<<" <-Vy "<<reco_nu_vtxz<<" <-Vz "<<std::endl;
	  
	  n_pfp_nuDaughters = pfp->NumDaughters();

	  // For CC0pi2p, we only consider the case with the number of neutrino daughters less than 4
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

		//All da space charge garbage!
		// Get time offset for x space charge correction
		auto const& detProperties_good = lar::providerFrom<detinfo::DetectorPropertiesService>();
		auto const& detClocks_good = lar::providerFrom<detinfo::DetectorClocksService>();
		auto const& mct_h_good = evt.getValidHandle<std::vector<simb::MCTruth> >("generator");
		auto gen_good = mct_h_good->at(0);
		double g4Ticks_good = detClocks_good->TPCG4Time2Tick(gen_good.GetNeutrino().Nu().T()) + detProperties_good->GetXTicksOffset(0,0,0) - detProperties_good->TriggerOffset();
		double xtimeoffset_good = detProperties_good->ConvertTicksToX(g4Ticks_good,0,0,0);
		trk_good_start = assoTrack.front()->Vertex<TVector3>();
		auto trk_good_start_offset = SCE->GetCalPosOffsets(geo::Point_t(trk_good_start.X(), trk_good_start.Y(), trk_good_start.Z()));
	    
		//auto trk_210_offset = SCE->GetCalPosOffsets(geo::Point_t(210, trk_good_start.Y(), trk_good_start.Z()));
		//trk_210_SCEcorr.SetX((210 - trk_210_offset.X() + xtimeoffset_good) + 0.6);
		//if(_debug) std::cout<<"[Numu0pi2p] Correction to X=210: "<<trk_210_SCEcorr.X()<<std::endl;
		
		trk_good_start_SCEcorr.SetX((trk_good_start.X() - trk_good_start_offset.X() + xtimeoffset_good) + 0.6);
		trk_good_start_SCEcorr.SetY(trk_good_start.Y() + trk_good_start_offset.Y());
		trk_good_start_SCEcorr.SetZ(trk_good_start.Z() + trk_good_start_offset.Z());

		if(_debug) std::cout<<"[Numu0pi2p] Location of the True Vertex: "<<trk_good_start.X()<<" <-Vx "<<trk_good_start.Y()<<" <-Vy "<<trk_good_start.Z()<<" <-Vz "<<std::endl;
		if(_debug) std::cout<<"[Numu0pi2p] Location of the True Vertex w/ SCE Correction: "<<trk_good_start_SCEcorr.X()<<" <-Vx "<<trk_good_start_SCEcorr.Y()<<" <-Vy "<<trk_good_start_SCEcorr.Z()<<" <-Vz "<<std::endl;
	    
		//defining the magnitude of the difference between the track start and the reconstructed neutrino start
		float reco_3d_diff = sqrt(pow((reco_nu_vtxx - trk_good_start_SCEcorr.X()),2) + pow((reco_nu_vtxy - trk_good_start_SCEcorr.Y()),2) + pow((reco_nu_vtxz - trk_good_start_SCEcorr.Z()),2));

		if(_debug) std::cout<<"[Numu0pi2p] Reconstructed Vertex Resolution: "<<reco_3d_diff<<std::endl;

		if(reco_3d_diff < vtxcut){ 
		  good_daughter_Tracks.push_back(assoTrack.front());
		}

	      } //finish loop that looks at the pfp assotiation for tracks

	      if(assoTrack.size()>1){ //means that there is more than one track associated to the daughter pfparticle
		throw cet::exception("[Numu0pi2p]") << "PFParticle has >1 track!" << std::endl;
	      }
	      // Collect pfparticle associated shower in a vector
	      auto assoShower = pfpToShowerAsso.at(dau_pfp.key()); // vector
	      if(assoShower.size()==1){
		daughter_Showers.push_back(assoShower.front());
	      }
	      if(assoShower.size()>1){ //means that there is more than one shower associated to the daughter pfparticle
		throw cet::exception("[Numu0pi2p]") << "PFParticle has >1 shower!" << std::endl;
	      }
	    } // finish looping of pfp
	  } //end of n_pfp > 4
     
      //Number of tracks and showers and good daughters (ie.e tracks that pass the reco cut)
      n_dau_tracks = daughter_Tracks.size();
      n_dau_showers = daughter_Showers.size();
      n_good_daughters = good_daughter_Tracks.size();
                  
      if(_debug) std::cout<<"[Numu0pi2p] Number of Daughter Tracks: "<<n_dau_tracks<<std::endl;
      if(_debug) std::cout<<"[Numu0pi2p] Number of Daughter Showers: "<<n_dau_showers<<std::endl;
      if(_debug) std::cout<<"[Numu0pi2p] Number of Good Daughters: "<<n_good_daughters<<std::endl;

      // Selection and Fill in Info
      if(n_dau_tracks == ntracks && n_dau_showers == nshowers){

	if(_debug) std::cout<<"[Numu0pi2p] Yay!!! Let's fill some shit!"<<std::endl;
	bool vtx_InFV = _fiducial_volume.InFV(trk_good_start_SCEcorr);
	
	if(n_good_daughters == ntracks ) {
	  return true;
	}
	else{
	  return false; //checking to make sure that the 
	}
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
    } //if single neutrino
  
    ClearLocalData();
  } //else loop over all the pfparticles

  return false; // If the function hasn't return before this, it means that there's no neutrino slice
  ClearLocalData();

}

void::ThreeTrackFilter::ClearLocalData()
{

  pfparticles.clear();
  pfparticles1.clear();
  pfneutrinos.clear();
  particlesToVertices.clear(); 
  reco_nu_vtxx=-9999.0;
  reco_nu_vtxy=-9999.0;
  reco_nu_vtxz=-9999.0;
  trk_good_start.Clear();
  trk_good_start_SCEcorr.Clear();
  trk_210_SCEcorr.Clear();

}

void ThreeTrackFilter::beginJob()
{
  // Implementation of optional member function here.
}

void ThreeTrackFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ThreeTrackFilter)

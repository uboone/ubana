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
  bool vtxInFid(float x, float y, float z,float xedge_min, float xedge_max, float yedge_min, float yedge_max, float zedge_min, float zedge_max);

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

  bool _debug = false;

  //Initiializing some important variables
  ////////////////////////////////////////
  int n_pfp_nuDaughters; // number of pfp which are the daughters of the neutrino
  int n_dau_tracks; // number of tracks asssociated to pfp neutrino daughters
  int n_dau_showers; // number of showers asssociated to pfp neutrino daughters
  int n_good_daughters; //number of daughters that are tracks and whose start is within 5cm of the vertex
  float reco_nu_vtxx = -9999., reco_nu_vtxy=-9999., reco_nu_vtxz=-9999.; //reconstructed neutrino vertex location
  TVector3 trk_start; 
  TVector3 trk_end;

  //Initializing some fhicl parameters
  //////////////////////////////////////
  int ntracks; //number of tracks for the filter
  int nshowers; //number of showers for the filter
  float vtxcut; //minimum distrance required from vertex to each track start
  float _border_x_low, _border_x_high, _border_y_low, _border_y_high, _border_z_low, _border_z_high;

  std::string m_pandoraLabel; //pandora
  std::string m_trackProducerLabel; //pandora or pandora track
  std::string m_showerProducerLabel; //pandora or pandora shower
  std::string m_pfp_producer; //pfp producer (pandora)   

  //Gathering all of our larsoft products
  /////////////////////////////////////////
  lar_pandora::LArPandoraHelper larpandora; //grabbing the helper class                                                                                 
  lar_pandora::PFParticleVector pfparticles,pfparticles1; //vector of pfparticles.                                                                      
  lar_pandora::PFParticleVector pfneutrinos; //vector of reconstructed neutrino vertices                                                                
  lar_pandora::PFParticlesToVertices particlesToVertices; //the output map from PFParticle to Vertex Objects                                            
  lar_pandora::VertexVector pfvertices; //the output vector of vertex objects                                                                           
  lar_pandora::PFParticlesToMetadata particlesToMetadata; //collection of PFParticleMetadata                                                            
  lar_pandora::PFParticleMap particleMap; //output mapping between reconstructed particles and particleID                                               
  lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints; //collects PFParticles & associated spacepoints from ART event record (evt in this case) 

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

  _debug = pset.get<bool>("Debug", "true"); //debug statements
  ntracks = pset.get<int>("NTrack", 3); //fhicl parameter that makes ther required number of tracks == 3
  nshowers = pset.get<int>("NShowers", 0); //"" number of showers == 0
  vtxcut =  pset.get<float>("VtxCut", 5.0); //require the start of each track to be within 5cm of the reconstructed vertex

  //Defining the padding from each detector edge
  _border_x_low = pset.get<float>("borderX_low", 5.0); 
  _border_x_high = pset.get<float>("borderX_high", 5.0); 
  _border_y_low = pset.get<float>("borderY_low", 5.0); 
  _border_y_high = pset.get<float>("borderY_high ", 5.0); 
  _border_z_low = pset.get<float>("borderZ_low", 5.0); 
  _border_z_high = pset.get<float>("borderZ_high", 5.0); 

}

//bool ThreeTrackFilter::filter(art::Event const& evt)
bool ThreeTrackFilter::filter(art::Event & evt)
{
  ClearLocalData();//clear the local data before starting

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
  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata); //collect reconstructed PFParticle metadata
  larpandora.BuildPFParticleMap(pfparticles, particleMap); //build particle map for reconstructed particles
  larpandora.CollectPFParticles(evt, m_pfp_producer, pfparticles1, particlesToSpacePoints); //collect reconstructed PFParticles
  larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos); //select reco neutrino particles from list of all reco particles
  larpandora.CollectVertices(evt, m_pandoraLabel, pfvertices, particlesToVertices); //collect reco PFParticles and associated Vertices from ART

  if(_debug) std::cout<<"[Numu0pi2p] Number of PFParticles in the Event: "<<pfparticles.size()<<std::endl;
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
	if(pfp->IsPrimary() && pfp->PdgCode() == abs(14)){
	  //if_Nu_Slice = true;

	  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();

	  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu); //vector of neutrino vertices
	  const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position(); //grab position of the neutrino vertex
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
	      if(assoTrack.size()==1){ //1 track association to the PFParticle in question
          
		daughter_Tracks.push_back(assoTrack.front());	    

		trk_start = assoTrack.front()->Vertex<TVector3>();
	    	trk_end = assoTrack.back()->Vertex<TVector3>();

		if(_debug) std::cout<<"[Numu0pi2p] Location of Reco Vertex: "<<" Vx: "<<trk_start.X()<<" Vy: "<<trk_start.Y()<<" Vz:  "<<trk_start.Z()<<std::endl;
		if(_debug) std::cout<<"[Numu0pi2p] Location of Reco Track End: "<<" Vx: "<<trk_end.X()<<" Vy: "<<trk_end.Y()<<" Vz:  "<<trk_end.Z()<<std::endl;

		//defining the magnitude of the difference between the track start and the reco neutrino start
		float reco_3d_diff_start = sqrt(pow((reco_nu_vtxx - trk_start.X()),2) + pow((reco_nu_vtxy - trk_start.Y()),2) + pow((reco_nu_vtxz - trk_start.Z()),2));

		//defining the magnitude of the difference between the track end and the reco neutrino start
		float reco_3d_diff_end = sqrt(pow((reco_nu_vtxx - trk_end.X()),2) + pow((reco_nu_vtxy - trk_end.Y()),2) + pow((reco_nu_vtxz - trk_end.Z()),2));

		if(_debug) std::cout<<"[Numu0pi2p] Reconstructed Vertex Resolution From Start: "<<reco_3d_diff_start<<std::endl;
		if(_debug) std::cout<<"[Numu0pi2p] Reconstructed Vertex Resolution From End: "<<reco_3d_diff_end<<std::endl;

		if(reco_3d_diff_start < vtxcut || reco_3d_diff_end < vtxcut){ 
		  good_daughter_Tracks.push_back(assoTrack.front());
		}

	      } //If associated track size == 1

	      if(assoTrack.size()>1){ //means that there is more than one track associated to the daughter pfparticle
		throw cet::exception("[Numu0pi2p]") << "PFParticle has >1 track!" << std::endl;
	      }
	     
	      auto assoShower = pfpToShowerAsso.at(dau_pfp.key());
	      if(assoShower.size()==1){ //exactly one association to the shower
		daughter_Showers.push_back(assoShower.front());
	      }
	      if(assoShower.size()>1){ //means that there is more than one shower associated to the daughter pfparticle
		throw cet::exception("[Numu0pi2p]") << "PFParticle has >1 shower!" << std::endl;
	      }
	    } // finishing the for loop over the neurino daughters
	  } //end of n_pfp > 4

       //////////////////////
       //Now for the good stuff
       ////////////////////////
      
      //Number of tracks and showers and good daughters (ie.e tracks that pass the reco cut)
      n_dau_tracks = daughter_Tracks.size();
      n_dau_showers = daughter_Showers.size();
      n_good_daughters = good_daughter_Tracks.size();
                  
      if(_debug) std::cout<<"[Numu0pi2p] Number of Daughter Tracks: "<<n_dau_tracks<<std::endl;
      if(_debug) std::cout<<"[Numu0pi2p] Number of Daughter Showers: "<<n_dau_showers<<std::endl;
      if(_debug) std::cout<<"[Numu0pi2p] Number of Good Daughters: "<<n_good_daughters<<std::endl;

      // Selection and Fill in Info
      //////////////////////////////
      if(n_dau_tracks == ntracks && n_dau_showers == nshowers){

	if(_debug) std::cout<<"[Numu0pi2p] Yay!!! Let's fill some shit!"<<std::endl;
	bool vtxInFV = vtxInFid(reco_nu_vtxx, reco_nu_vtxy, reco_nu_vtxz, _border_x_low, _border_x_high, _border_y_low, _border_y_high, _border_z_low, _border_z_high); 

	//if(_debug) std::cout<<"[Numu0pi2p] Value of vtx_InFV: "<<vtx_InFV<<std::endl;
	if(_debug) std::cout<<"[Numu0pi2p] Value of vtxInFV: "<<vtxInFV<<std::endl;

	//Check and make sure there are exactly 3 tracks that are within the 5cm of the vertex
	if(n_good_daughters == ntracks ) {
	  return true;
	}
	else{
	  return false; 
	}
        
	//Checking if the vertex is indeed within the FV
	if(vtxInFV == true){
          return true;
        } else {
          return false; 
	} 

      }else { //if not 3 track 0 shwr topology
        return false;        
      }

	} // End of the neutrino loop
      } // Finish loop all the pfp particles
    }// if the size of pfpneutrinos==1
  
    ClearLocalData();//clear the local data before leaving
  } //else loop over all the pfparticles: end of the event

  return false; // If the function hasn't return before this, it means that there's no neutrino slice

} //End Program

void::ThreeTrackFilter::ClearLocalData()
{

  pfparticles.clear();
  pfparticles1.clear();
  pfneutrinos.clear();
  particlesToVertices.clear(); 
  reco_nu_vtxx=-9999.0;
  reco_nu_vtxy=-9999.0;
  reco_nu_vtxz=-9999.0;
  trk_start.Clear();
  trk_end.Clear();

}

bool::ThreeTrackFilter::vtxInFid(float x, float y, float z,float xedge_min, float xedge_max, float yedge_min, float yedge_max, float zedge_min, float zedge_max)
{

  //Defining the edges of the detector (hard coding cause I don't know why this is broken                                                                 
  float xmin = 0.0 + xedge_min;
  float xmax = 256.35 - xedge_max;
  float ymin = -116.5 + yedge_min;
  float ymax = 116.5 - yedge_max;
  float zmin = 0.0 + zedge_min;
  float zmax = 1036.8 - zedge_max;

    if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax) ){
    return false;
  } else{
    return true;
  }


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

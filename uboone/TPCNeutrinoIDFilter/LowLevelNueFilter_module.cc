////////////////////////////////////////////////////////////////////////
// Class:       LowLevelNueFilter
// Plugin Type: analyzer (art v2_08_03)
// File:        LowLevelNueFilter_module.cc
//
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larcore/Geometry/Geometry.h"

#include "TH1F.h"

namespace ub {
  class LowLevelNueFilter;
}


class ub::LowLevelNueFilter : public art::EDAnalyzer {
public:
  explicit LowLevelNueFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LowLevelNueFilter(LowLevelNueFilter const &) = delete;
  LowLevelNueFilter(LowLevelNueFilter &&) = delete;
  LowLevelNueFilter & operator = (LowLevelNueFilter const &) = delete;
  LowLevelNueFilter & operator = (LowLevelNueFilter &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fOpFlash_tag;
  art::InputTag fMCtruth_tag;
  art::InputTag fCluster_tag;
  art::InputTag fVertex_tag;
  art::InputTag fPFParticle_tag;

  TH1F* h_flash_per_ev;
  TH1F* h_flash_pe;
  TH1F* h_flash_y;
  TH1F* h_flash_z;
  TH1F* h_flash_time;
  TH1F* h_ophits_per_flash;
  TH1F* h_ophits_per_flash_2pe;

};


ub::LowLevelNueFilter::LowLevelNueFilter(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fOpFlash_tag(p.get<art::InputTag>("OpFlash_tag")),
  fMCtruth_tag(p.get<art::InputTag>("MCtruth_tag")),
  fCluster_tag(p.get<art::InputTag>("Cluster_tag")),
  fVertex_tag(p.get<art::InputTag>("Vertex_tag")),
  fPFParticle_tag(p.get<art::InputTag>("PFParticle_tag"))
 // More initializers here.
{}

void ub::LowLevelNueFilter::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  std::cout << "Processing "
            << "Run " << e.run() << ", "
            << "SubRun "<<e.subRun()<<", "
            << "Event " << e.id().event() << std::endl;

  art::ServiceHandle<geo::Geometry> geo;

  art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
  std::vector<art::Ptr<simb::MCTruth>> MCtruth_vec;
  if (e.getByLabel(fMCtruth_tag, MCtruthHandle)){
    art::fill_ptr_vector(MCtruth_vec, MCtruthHandle);
  }
  
  std::cout<<"MCtruth size is "<<MCtruth_vec.size()<<std::endl;

    int CCNC=999;//cc=0, nc=1
    int mode=999;//QE=0, RES=1, DIS=2, Coherent =3, MEC=10
    int NuPDG=999; //nue=12, nuebar=-12,numu=14,numubar=-14
    float enutruth=0; //in GeV
    float nuvtxx=0;
    float nuvtxy=0;
    float nuvtxz=0;
    float lep_mom_truth=0;//in GeV
    if(MCtruth_vec[0]->NeutrinoSet()){
      NuPDG = MCtruth_vec[0]->GetNeutrino().Nu().PdgCode();
      CCNC = MCtruth_vec[0]->GetNeutrino().CCNC();
      mode = MCtruth_vec[0]->GetNeutrino().Mode();
      enutruth = MCtruth_vec[0]->GetNeutrino().Nu().E();
      nuvtxx = MCtruth_vec[0]->GetNeutrino().Nu().Vx();
      nuvtxy = MCtruth_vec[0]->GetNeutrino().Nu().Vy();
      nuvtxz = MCtruth_vec[0]->GetNeutrino().Nu().Vz();
      lep_mom_truth =  MCtruth_vec[0]->GetNeutrino().Lepton().P();
    }
    std::cout<<"NuPDG = "<<NuPDG<<" mode = "<<mode<<" CCNC = "<<CCNC<<" nue energy = "<<enutruth<<" nu Vertex = ["<<nuvtxx<<" , "<<nuvtxy<<" , "<<nuvtxz<<"]"<<" lepton momentum = "<<lep_mom_truth<<std::endl;

    //Flashes
    // Clusters
    art::Handle<std::vector<recob::OpFlash> > opflashHandle;
    std::vector<art::Ptr<recob::OpFlash> > OpFlash_vec;
    if (e.getByLabel(fOpFlash_tag, opflashHandle))
      art::fill_ptr_vector(OpFlash_vec, opflashHandle);

    for (auto & flash : OpFlash_vec){
      
      float wire0 = FLT_MAX;
      float wire1 = FLT_MIN;
      //Find the 4 corners and convert them to wire numbers
      std::vector<TVector3> points;
      points.push_back(TVector3(0, flash->YCenter()-flash->YWidth(), flash->ZCenter()-flash->ZWidth()));
      points.push_back(TVector3(0, flash->YCenter()-flash->YWidth(), flash->ZCenter()+flash->ZWidth()));
      points.push_back(TVector3(0, flash->YCenter()+flash->YWidth(), flash->ZCenter()-flash->ZWidth()));
      points.push_back(TVector3(0, flash->YCenter()+flash->YWidth(), flash->ZCenter()+flash->ZWidth()));
      for (int pl = 0; pl <3; ++pl){
        geo::PlaneID pid(0, 0, pl);
        for (size_t i = 0; i<points.size(); ++i){
          geo::WireID wireID;
          try{
            wireID = geo->NearestWireID(points[i], pid);
          }
          catch(geo::InvalidWireError const& e) {
            wireID = e.suggestedWireID(); // pick the closest valid wire
          }
          if (wireID.Wire < wire0) wire0 = wireID.Wire;
          if (wireID.Wire > wire1) wire1 = wireID.Wire;
        }
        std::cout<<"Flash "<<flash.key()<<" on plane "<<pl<<" from wire "<<wire0<<" to wire "<<wire1<<std::endl;
      }
    }

    // Clusters
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    std::vector<art::Ptr<recob::Cluster> > Cluster_vec;
    if (e.getByLabel(fCluster_tag, clusterHandle))
      art::fill_ptr_vector(Cluster_vec, clusterHandle);
    
    for( auto const& cluster : Cluster_vec){
      auto ID = cluster->ID();
      auto plane = cluster->Plane().Plane;
      //cout<<"cluster ID is "<<ID<<" in plane "<<plane<<endl;
      if(((int)ID)<0)
	{
	  auto nHits = cluster->NHits();
	  auto StartWire = cluster->StartWire();
	  auto StartTick = cluster->StartTick();
	  auto StartCharge = cluster->StartCharge();
	  auto EndCharge = cluster->EndCharge();
	  auto TotalCharge = cluster->Integral();
          std::cout<<"ClusterID="<<ID<<" plane="<<plane<<"nhits = "<<nHits<<" startwire = "<<StartWire<<" starttick = "<<StartTick<<" StartCharge = "<<StartCharge<<" EndCharge = "<<EndCharge<<" TotalCharge = "<<TotalCharge<<std::endl;
	}
    }

    // Vertices
    art::Handle<std::vector<recob::Vertex> > vertexHandle;
    std::vector<art::Ptr<recob::Vertex> > Vertex_vec;
    if (e.getByLabel(fVertex_tag, vertexHandle))
      art::fill_ptr_vector(Vertex_vec, vertexHandle);

    std::cout<<"There are "<<Vertex_vec.size()<<" vertices in this event."<<std::endl;
    for( auto const& vertex : Vertex_vec){
      auto VertexID = vertex->ID();
      double* vertexXYZ = new double;
      vertex->XYZ(vertexXYZ);
      
      std::cout<<"Vertex ID "<<VertexID<<" XYZ = ["<<vertexXYZ[0]<<","<<vertexXYZ[1]<<","<<vertexXYZ[2]<<"]"<<std::endl;
    }
    
    // PFParticle
    art::Handle<std::vector<recob::PFParticle> > pfparticleHandle;
    std::vector<art::Ptr<recob::PFParticle> > PFParticle_vec;
    if (e.getByLabel(fPFParticle_tag, pfparticleHandle))
      art::fill_ptr_vector(PFParticle_vec, pfparticleHandle);

    std::cout<<"There are "<<PFParticle_vec.size()<<" PFParticle in this event."<<std::endl;
    for( auto const& pfparticle : PFParticle_vec){
      auto pdgcode = pfparticle->PdgCode();
      auto isprimary = pfparticle->IsPrimary();  
      auto kprimary = pfparticle->kPFParticlePrimary;
      std::cout<<"pfparticle pdgcode "<<pdgcode<<" is primary?"<<isprimary<<" kprimary "<<kprimary<<std::endl;
    }



}

void ub::LowLevelNueFilter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  h_flash_per_ev = tfs->make<TH1F>("h_flash_per_ev","OpFlashes per event;N_{flashes};Events / bin",20,-0.5,19.5);
  h_flash_pe = tfs->make<TH1F>("h_flash_pe","Flash PEs; PE; Events / 0.1 PE",100,0,50);
  h_flash_y = tfs->make<TH1F>("h_flash_y","Flash y position; y (cm); Events / 0.1 cm",100,-200,200);
  h_flash_z = tfs->make<TH1F>("h_flash_z","Flash z position; z (cm); Events / 0.1 cm",100,-100,1100);
  h_flash_time = tfs->make<TH1F>("h_flash_time","Flash Time; time (#mus); Events / 0.5 #mus",60,-5,25);
  h_ophits_per_flash = tfs->make<TH1F>("h_ophits_per_flash","OpHits per Flash;N_{optical hits};Events / bin",20,-0.5,19.5);
  h_ophits_per_flash_2pe = tfs->make<TH1F>("h_ophits_per_flash_2pe","OpHits (> 2 PE) per Flash;N_{optical hits};Events / bin",20,-0.5,19.5);
}

DEFINE_ART_MODULE(ub::LowLevelNueFilter)

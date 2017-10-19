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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
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
#include "TH1I.h"

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
  TH1F* h_flash_z_diff_nue;
  TH1F* h_flash_z_diff_cosmic;
  TH1F* h_flash_time_diff_nue;
  TH1F* h_flash_time_diff_cosmic;

  TH1F* h_flashWidth0;
  TH1F* h_flashWidth1;
  TH1F* h_flashWidth2;

  TH1F* h_StartWire_diff_nue_0;
  TH1F* h_StartWire_diff_numu_0;
  TH1F* h_StartWire_diff_cosmic_0;
  TH1F* h_EndWire_diff_nue_0;
  TH1F* h_EndWire_diff_numu_0;
  TH1F* h_EndWire_diff_cosmic_0;
  TH1F* h_CenterWire_diff_nue_0;
  TH1F* h_CenterWire_diff_numu_0;
  TH1F* h_CenterWire_diff_cosmic_0;
  TH1F* h_FlashClusterDistance_nue_0;
  TH1F* h_FlashClusterDistance_cosmic_0;
  TH1F* h_nhits_nue_0;
  TH1F* h_nhits_cosmic_0;

  TH1F* h_StartWire_diff_nue_1;
  TH1F* h_StartWire_diff_numu_1;
  TH1F* h_StartWire_diff_cosmic_1;
  TH1F* h_EndWire_diff_nue_1;
  TH1F* h_EndWire_diff_numu_1;
  TH1F* h_EndWire_diff_cosmic_1;
  TH1F* h_CenterWire_diff_nue_1;
  TH1F* h_CenterWire_diff_numu_1;
  TH1F* h_CenterWire_diff_cosmic_1;
  TH1F* h_FlashClusterDistance_nue_1;
  TH1F* h_FlashClusterDistance_cosmic_1;
  TH1F* h_nhits_nue_1;
  TH1F* h_nhits_cosmic_1;

  TH1F* h_StartWire_diff_nue_2;
  TH1F* h_StartWire_diff_numu_2;
  TH1F* h_StartWire_diff_cosmic_2;
  TH1F* h_EndWire_diff_nue_2;
  TH1F* h_EndWire_diff_numu_2;
  TH1F* h_EndWire_diff_cosmic_2;
  TH1F* h_CenterWire_diff_nue_2;
  TH1F* h_CenterWire_diff_numu_2;
  TH1F* h_CenterWire_diff_cosmic_2;
  TH1F* h_FlashClusterDistance_nue_2;
  TH1F* h_FlashClusterDistance_cosmic_2;
  TH1F* h_nhits_nue_2;
  TH1F* h_nhits_cosmic_2;

  TH1F* h_Stat;
  TH1F* h_ophits_per_flash;
  TH1F* h_ophits_per_flash_2pe;
  TH1F* h_ClusterNegPdgcode;

  void GetTruthInfo(std::vector<art::Ptr<recob::Hit>> hits, int& origin, int& pdgcode);
  bool inFV(double x, double y, double z);

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

  //Define a bunch of flag
  bool FlashInBeamWindow = false;
  bool TruthSignal = false;
  bool ClusterFlashMatchNue0 = false;
  bool ClusterFlashMatchCosmic0 = false;
  bool ClusterFlashMatch0 = false;
  bool ClusterFlashMatchNue1 = false;
  bool ClusterFlashMatchCosmic1 = false;
  bool ClusterFlashMatch1=false;
  bool ClusterFlashMatchNue2 = false;
  bool ClusterFlashMatchCosmic2 = false;
  bool ClusterFlashMatch2=false;
  //bool AtLeastTwoNegativeCluster = false;
  //bool ClusterMatchFlash = false;
  //bool HaveEPFParticle = false;
  //bool Have3DShower = false;
  //bool PassFilter = false;
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
    //h_Stat->AddBinContent(1);
    h_Stat->SetBinContent(1,h_Stat->GetBinContent(1)+1);
    if(inFV(nuvtxx,nuvtxy,nuvtxz) && CCNC==0 && NuPDG==12){
      //h_Stat->AddBinContent(2);
      h_Stat->SetBinContent(2,h_Stat->GetBinContent(2)+1);
      TruthSignal = true;
    }
    //Flashes
    art::Handle<std::vector<recob::OpFlash> > opflashHandle;
    std::vector<art::Ptr<recob::OpFlash> > OpFlash_vec;
    if (e.getByLabel(fOpFlash_tag, opflashHandle))
      art::fill_ptr_vector(OpFlash_vec, opflashHandle);
    
    std::vector<TVector3> FlashWire;
    for (auto & flash : OpFlash_vec){
      if(flash->Time()<3 || flash->Time()>5)
	continue;
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
        std::cout<<"There are "<<OpFlash_vec.size()<<", Flash "<<flash.key()<<" on plane "<<pl<<" from wire "<<wire0<<" to wire "<<wire1<<" total PE "<<flash->TotalPE()<<" time "<<flash->Time()<<std::endl;
	FlashWire.push_back(TVector3(pl,wire0,wire1));
	FlashInBeamWindow = true;
	if(pl==0)
	  {h_flashWidth0->Fill(wire1-wire0);}
	if(pl==1)
	  {h_flashWidth1->Fill(wire1-wire0);}
	if(pl==2)
	  {h_flashWidth2->Fill(wire1-wire0);}
      }
    }
    for(size_t i=0;i<FlashWire.size();i++)
      {
	std::cout<<"plane "<<FlashWire[i].X()<<" wire start "<<FlashWire[i].Y()<<" wire end "<<FlashWire[i].Z()<<std::endl;
      }
    if(FlashInBeamWindow==true && TruthSignal==true)
      {
	//h_Stat->AddBinContent(3);
	h_Stat->SetBinContent(3, h_Stat->GetBinContent(3)+1);
      }
    // Clusters
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    std::vector<art::Ptr<recob::Cluster> > Cluster_vec;
    if (e.getByLabel(fCluster_tag, clusterHandle))
      art::fill_ptr_vector(Cluster_vec, clusterHandle);

    art::FindManyP<recob::Hit>  fmhc(clusterHandle, e, fCluster_tag);
    
    for( auto const& cluster : Cluster_vec){
      auto ID = cluster->ID();
      auto plane = cluster->Plane().Plane;
      //cout<<"cluster ID is "<<ID<<" in plane "<<plane<<endl;
      if(((int)ID)<0)
	{
	  auto nHits = cluster->NHits();
	  auto StartWire = cluster->StartWire();
	  auto StartTick = cluster->StartTick();
	  auto EndWire = cluster->EndWire();

	  // auto TotalCharge = cluster->Integral();
          int origin = 0;
          int pdgcode = 0;
          GetTruthInfo(fmhc.at(cluster.key()), origin, pdgcode);
          std::cout<<"ClusterID="<<ID<<" plane="<<plane<<"nhits = "<<nHits<<" startwire = "<<StartWire<<" starttick = "<<StartTick<<" End wire "<<EndWire<< " origin = "<<origin<<" pdgcode = "<<pdgcode<<std::endl;
	  //does the cluster match with the nu flash?
	  for(size_t i=0; i<FlashWire.size();i++)
	    {
	      if(FlashWire[i].X()!=plane)
		continue;
	      //if(FlashWire[i].X()==plane)
	      //{
		  auto ClusterFlashStartDiff = FlashWire[i].Y()-StartWire;
		  auto ClusterFlashEndDiff = FlashWire[i].Z()-EndWire;
		  if(origin==1 && pdgcode == 11)
		    {
		      if(plane==0){
			h_StartWire_diff_nue_0->Fill(ClusterFlashStartDiff);
			h_EndWire_diff_nue_0->Fill(ClusterFlashEndDiff);}
		      if(plane==1){
			h_StartWire_diff_nue_1->Fill(ClusterFlashStartDiff);
			h_EndWire_diff_nue_1->Fill(ClusterFlashEndDiff);}
		      if(plane==2){
			std::cout<<"plane 2 nue cluster start wire is "<<StartWire<<" End wire is "<<EndWire<<" start tick is "<<StartTick<<" nhits"<<nHits<<std::endl;
			h_StartWire_diff_nue_2->Fill(ClusterFlashStartDiff);
			h_EndWire_diff_nue_2->Fill(ClusterFlashEndDiff);}
		      
		    }
		  if(origin==1 && pdgcode == 13)
		     {
		       if(plane==0){
			 h_StartWire_diff_numu_0->Fill(ClusterFlashStartDiff);
			 h_EndWire_diff_numu_0->Fill(ClusterFlashEndDiff);}
		       if(plane==1){
			 h_StartWire_diff_numu_1->Fill(ClusterFlashStartDiff);
			 h_EndWire_diff_numu_1->Fill(ClusterFlashEndDiff);}
		       if(plane==2){
			 h_StartWire_diff_numu_2->Fill(ClusterFlashStartDiff);
			 h_EndWire_diff_numu_2->Fill(ClusterFlashEndDiff);			 
		       }
		       
		     }
		  if(origin == 2)
		    {
		      if(plane==0){
			h_StartWire_diff_cosmic_0->Fill(ClusterFlashStartDiff);
			h_EndWire_diff_cosmic_0->Fill(ClusterFlashEndDiff);}
		      if(plane==1){
			h_StartWire_diff_cosmic_1->Fill(ClusterFlashStartDiff);
			h_EndWire_diff_cosmic_1->Fill(ClusterFlashEndDiff);}
		      if(plane==2){
			h_StartWire_diff_cosmic_2->Fill(ClusterFlashStartDiff);
			h_EndWire_diff_cosmic_2->Fill(ClusterFlashEndDiff);
			std::cout<<"plane 2 nue cluster start wire is "<<StartWire<<" End wire is "<<EndWire<<" start tick is "<<StartTick<<" nhits"<<nHits<<std::endl;
		      }
		    }
		
		  //if(StartWire>FlashWire[i].Y() && EndWire<FlashWire[i].Z())
		  //  {
		  auto ClusterWireCenter = (StartWire + EndWire)/2.;
		  auto FlashWireCenter = (FlashWire[i].Y() + FlashWire[i].Z())/2.;
		  auto ClusterFlashCenterDiff = std::fabs(FlashWireCenter - ClusterWireCenter);
		  double FlashClusterDistance=9999;
		  if(StartWire>FlashWire[i].Y() && StartWire<FlashWire[i].Z() && EndWire>FlashWire[i].Y() && EndWire<FlashWire[i].Z())//flash include cluster both ends
		    {FlashClusterDistance=0;}
		  if(StartWire<=FlashWire[i].Y() || EndWire>=FlashWire[i].Z())
		    {FlashClusterDistance = ClusterFlashCenterDiff;}
		  if(origin==1 && pdgcode == 11)
		    {
		      if(plane==0){
			h_CenterWire_diff_nue_0->Fill(ClusterFlashCenterDiff);
			h_FlashClusterDistance_nue_0->Fill(FlashClusterDistance);
			if(FlashClusterDistance<500)
			  {h_nhits_nue_0->Fill(nHits);
			    std::cout<<"Event Dis nue plane 0 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			  }
		      }
			if(plane==1){
			  h_CenterWire_diff_nue_1->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_nue_1->Fill(FlashClusterDistance);
			  if(FlashClusterDistance<500)
			    {h_nhits_nue_1->Fill(nHits);
			      std::cout<<"Event Dis nue plane 1 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			    }
			}
			if(plane==2){
			  h_CenterWire_diff_nue_2->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_nue_2->Fill(FlashClusterDistance);
			  if(FlashClusterDistance<500)
			    {h_nhits_nue_2->Fill(nHits);
			      std::cout<<"Event Dis nue plane 2 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			    }
			}
		    }
			
		  if(origin==1 && pdgcode == 13)
		    {
		      if(plane==0){
			h_CenterWire_diff_numu_0->Fill(ClusterFlashCenterDiff);}
		      if(plane==1){
			h_CenterWire_diff_numu_1->Fill(ClusterFlashCenterDiff);}
		      if(plane==2){
			h_CenterWire_diff_numu_2->Fill(ClusterFlashCenterDiff);}
		    }
		  if(origin==2)
		    {
		      if(plane==0){
			h_CenterWire_diff_cosmic_0->Fill(ClusterFlashCenterDiff);
			h_FlashClusterDistance_cosmic_0->Fill(FlashClusterDistance);
			if(FlashClusterDistance<500)
			  {h_nhits_cosmic_0->Fill(nHits);
			    std::cout<<"Event Dis cosmic plane 0 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			  }
		      }
		      if(plane==1){
			h_CenterWire_diff_cosmic_1->Fill(ClusterFlashCenterDiff);
			h_FlashClusterDistance_cosmic_1->Fill(FlashClusterDistance);
			if(FlashClusterDistance<500)
			  {h_nhits_cosmic_1->Fill(nHits);
			     std::cout<<"Event Dis cosmic plane 1 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			  }
		      }
		      if(plane==2){
			h_CenterWire_diff_cosmic_2->Fill(ClusterFlashCenterDiff);
			h_FlashClusterDistance_cosmic_0->Fill(FlashClusterDistance);
			if(FlashClusterDistance<500)
			  {h_nhits_cosmic_2->Fill(nHits);
			     std::cout<<"Event Dis cosmic plane 2 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			    
			  }
		      }
		    }
		  if(plane==0 && (ClusterFlashStartDiff>-600 && ClusterFlashEndDiff<0) && ClusterFlashCenterDiff<300)
		    {
		      ClusterFlashMatch0=true;
		      if(origin==1 && pdgcode ==11)
			{ClusterFlashMatchNue0 = true;}
		      if(origin==2)
			{ClusterFlashMatchCosmic0 = true;}
		    }
		  if(plane==1 && (ClusterFlashStartDiff>-700 && ClusterFlashEndDiff<50) && ClusterFlashCenterDiff<350)
		    {
		      ClusterFlashMatch1=true;
		      if(origin==1 && pdgcode ==11)
			{
			  ClusterFlashMatchNue1 = true;
			}
		      if(origin==2)
			{ClusterFlashMatchCosmic1 = true;}
		    }	
		  if(plane==2 && ClusterFlashCenterDiff<800)
		    {
		      ClusterFlashMatch2=true;
		      if(origin==1&&pdgcode==11)
			ClusterFlashMatchNue2 = true;
		      if(origin==2)
			ClusterFlashMatchCosmic2 = true;
		    }
		  //}//flash wire cover the cluster wire
		  //}//flash plane match to the cluster plane
	    }//flash plane loop
	}//ID<0, //origin = 1(beam neutrino), 2 (cosmic),3(supernova),4(single particle) pdgcode = 11(e-), 13(mu-),22(gamma) 
      if(FlashInBeamWindow==true && TruthSignal==true)
	{
	  if(ClusterFlashMatchNue0==true)
	    {//h_Stat->AddBinContent(4);
	      h_Stat->SetBinContent(4, h_Stat->GetBinContent(4)+1); 
	    }
	  if(ClusterFlashMatchCosmic0==true)
	    {//h_Stat->AddBinContent(5);
	      h_Stat->SetBinContent(5,h_Stat->GetBinContent(5)+1);
	    }
	  if(ClusterFlashMatchNue1==true)
	    {//h_Stat->AddBinContent(6);
	      h_Stat->SetBinContent(6,h_Stat->GetBinContent(6)+1);
	    }
	  if(ClusterFlashMatchCosmic1==true)
	    {//h_Stat->AddBinContent(7);
	      h_Stat->SetBinContent(7,h_Stat->GetBinContent(7)+1);
	    }
	  if(ClusterFlashMatchNue2==true)
	    {//h_Stat->AddBinContent(8);
	      h_Stat->SetBinContent(8,h_Stat->GetBinContent(8)+1);
	    }
	  if(ClusterFlashMatchCosmic2==true)
	    {//h_Stat->AddBinContent(9);
	      h_Stat->SetBinContent(9,h_Stat->GetBinContent(9)+1);
	    }
	}
    }//cluster loop
    if(FlashInBeamWindow==true && TruthSignal==true)
      {
	if(ClusterFlashMatch0==true || ClusterFlashMatch1==true || ClusterFlashMatch2==true)
	  {
	    //h_Stat->AddBinContent(10);
	    h_Stat->SetBinContent(10,h_Stat->GetBinContent(10)+1);
	    if(ClusterFlashMatchNue0==true || ClusterFlashMatchNue1==true || ClusterFlashMatchNue2==true)
	      {//h_Stat->AddBinContent(11);
		h_Stat->SetBinContent(11,h_Stat->GetBinContent(11)+1);
	      }
	  }	
	if((ClusterFlashMatch0==true && ClusterFlashMatch1==true) || (ClusterFlashMatch0==true && ClusterFlashMatch2==true) || (ClusterFlashMatch1==true && ClusterFlashMatch2==true))//at least 2 flash matched clusters one from each plane
	  {
	    //h_Stat->AddBinContent(12);
	    h_Stat->SetBinContent(12,h_Stat->GetBinContent(12)+1);
	    if(ClusterFlashMatchNue0==true || ClusterFlashMatchNue1==true || ClusterFlashMatchNue2==true)
	      {
		//h_Stat->AddBinContent(13);
		h_Stat->SetBinContent(13,h_Stat->GetBinContent(13)+1);
	      }
	  }
      }
    

    
    // Vertices, to associate with PFParticle?
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

    art::FindOneP<recob::Vertex> PrimaryVertexFromPFP(pfparticleHandle, e, fPFParticle_tag);

    std::cout<<"There are "<<PFParticle_vec.size()<<" PFParticle in this event."<<std::endl;
    if(TruthSignal == true && FlashInBeamWindow==true){
      for( auto const& pfparticle : PFParticle_vec){
	auto pdgcode = pfparticle->PdgCode();
	std::cout<<"pdgcode is "<<pdgcode<<std::endl;
	if(pdgcode==11 || pdgcode==1111){
	  auto PrimaryVertex = PrimaryVertexFromPFP.at(pfparticle.key());
	  std::cout<<"key = "<<pfparticle.key()<<" pfparticle size "<<PrimaryVertexFromPFP.size()<<std::endl;
	  if(PrimaryVertex)
	    {
	      auto PrimaryVertexID = PrimaryVertex->ID();
	      double* PrimaryVertexXYZ = new double;
	      PrimaryVertex->XYZ(PrimaryVertexXYZ);
	      std::cout<<"pfparticle reco pdgcode "<<pdgcode<<" associated vertex ID is"<<PrimaryVertexID<<" XYZ= ["<<PrimaryVertexXYZ[0]<<","<<PrimaryVertexXYZ[1]<<","<<PrimaryVertexXYZ[2]<<"] !!!!!"<<std::endl;
	      //h_Stat->AddBinContent(14);
	      h_Stat->SetBinContent(14,h_Stat->GetBinContent(14)+1);
	    }
	}
      }
    }
}


void ub::LowLevelNueFilter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  h_flash_per_ev = tfs->make<TH1F>("h_flash_per_ev","OpFlashes per event;N_{flashes};Events / bin",20,-0.5,19.5);
  h_flash_pe = tfs->make<TH1F>("h_flash_pe","Flash PEs; PE; Events / 0.1 PE",100,0,50);
  h_flash_y = tfs->make<TH1F>("h_flash_y","Flash y position; y (cm); Events / 0.1 cm",100,-200,200);
  h_flash_z_diff_nue = tfs->make<TH1F>("h_flash_z_diff_nue","Flash z position; z (cm); Events / 0.1 cm",100,-100,1100);
  h_flash_z_diff_cosmic = tfs->make<TH1F>("h_flash_z_diff_cosmic","Flash z position; z (cm); Events / 0.1 cm",100,-100,1100);
  h_flash_time_diff_nue = tfs->make<TH1F>("h_flash_timdiff_nue","Flash Time; time (#mus); Events / 0.5 #mus",60,-5,25);
  h_flash_time_diff_cosmic = tfs->make<TH1F>("h_flash_time_cosmic","Flash Time; time (#mus); Events / 0.5 #mus",60,-5,25);
  
  h_flashWidth0 = tfs->make<TH1F>("h_flashWidth0","flash width projected to plane 0",200,0,4000);
  h_flashWidth1 = tfs->make<TH1F>("h_flashWidth1","flash width projected to plane 0",200,0,4000);
  h_flashWidth2 = tfs->make<TH1F>("h_flashWidth2","flash width projected to plane 0",200,0,4000);
  
  //plane 0
  h_StartWire_diff_nue_0 = tfs->make<TH1F>("h_StartWire_diff_nue_0","Start wire difference between cluster and beam flash",500,-2000,2000);
  h_EndWire_diff_nue_0 = tfs->make<TH1F>("h_EndWire_diff_nue_0","End wire difference between cluster and beam flash",500,-2000,2000);
  h_StartWire_diff_numu_0 = tfs->make<TH1F>("h_StartWire_diff_numu_0","Start wire difference between cluster and beam flash",100,-2000,2000);
  h_EndWire_diff_numu_0 = tfs->make<TH1F>("h_EndWire_diff_numu_0","End wire difference between cluster and beam flash",100,-2000,2000);
  h_StartWire_diff_cosmic_0 = tfs->make<TH1F>("h_StartWire_diff_cosmic_0","Start wire difference between cluster and beam flash",100,-2000,2000);
  h_EndWire_diff_cosmic_0 = tfs->make<TH1F>("h_EndWire_diff_cosmic_0","End wire difference between cluster and beam flash",100,-2000,2000);
  h_CenterWire_diff_nue_0= tfs->make<TH1F>("h_CenterWire_diff_nue_0","center wire difference between cluster and beam flash",500,0,4000);
  h_CenterWire_diff_numu_0= tfs->make<TH1F>("h_CenterWire_diff_numu_0","center wire difference between cluster and beam flash",100,0,4000);
  h_CenterWire_diff_cosmic_0= tfs->make<TH1F>("h_CenterWire_diff_cosmic_0","center wire difference between cluster and beam flash",100,0,4000);
  h_FlashClusterDistance_nue_0= tfs->make<TH1F>("h_FlashClusterDistance_nue_0","flash cluster distance of nue plane 0",100,0,4000);
  h_FlashClusterDistance_cosmic_0= tfs->make<TH1F>("h_FlashClusterDistance_cosmic_0","flash cluster distance of cosmic plane0",100,0,4000);
  h_nhits_nue_0 = tfs->make<TH1F>("h_nhits_nue_0","number of the hits of the e- cluster with FlashClusterDistance < 500",1000,0,1000);
  h_nhits_cosmic_0 = tfs->make<TH1F>("h_nhits_cosmic_0","number of the hits of the e- cluster with FlashClusterDistance < 500",1000,0,1000);
//plane 1
  h_StartWire_diff_nue_1 = tfs->make<TH1F>("h_StartWire_diff_nue_1","Start wire difference between cluster and beam flash",500,-2000,2000);
  h_EndWire_diff_nue_1 = tfs->make<TH1F>("h_EndWire_diff_nue_1","End wire difference between cluster and beam flash",500,-2000,2000);
  h_StartWire_diff_numu_1 = tfs->make<TH1F>("h_StartWire_diff_numu_1","Start wire difference between cluster and beam flash",100,-2000,2000);
  h_EndWire_diff_numu_1 = tfs->make<TH1F>("h_EndWire_diff_numu_1","End wire difference between cluster and beam flash",100,-2000,2000);
  h_StartWire_diff_cosmic_1 = tfs->make<TH1F>("h_StartWire_diff_cosmic_1","Start wire difference between cluster and beam flash",100,-2000,2000);
  h_EndWire_diff_cosmic_1 = tfs->make<TH1F>("h_EndWire_diff_cosmic_1","End wire difference between cluster and beam flash",100,-2000,2000);
  h_CenterWire_diff_nue_1= tfs->make<TH1F>("h_CenterWire_diff_nue_1","center wire difference between cluster and beam flash",500,0,4000);
  h_CenterWire_diff_numu_1= tfs->make<TH1F>("h_CenterWire_diff_numu_1","center wire difference between cluster and beam flash",100,0,4000);
  h_CenterWire_diff_cosmic_1= tfs->make<TH1F>("h_CenterWire_diff_cosmic_1","center wire difference between cluster and beam flash",100,0,4000);
  h_FlashClusterDistance_nue_1= tfs->make<TH1F>("h_FlashClusterDistance_nue_1","flash cluster distance of nue plane 1",100,0,4000);
  h_FlashClusterDistance_cosmic_1= tfs->make<TH1F>("h_FlashClusterDistance_cosmic_1","flash cluster distance of cosmic plane 1",100,0,4000);
  h_nhits_nue_1 = tfs->make<TH1F>("h_nhits_nue_1","number of the hits of the e- cluster with FlashClusterDistance < 500",1000,0,1000);
  h_nhits_cosmic_1 = tfs->make<TH1F>("h_nhits_cosmic_1","number of the hits of the e- cluster with FlashClusterDistance < 500",1000,0,1000);
  //plane 2
  h_StartWire_diff_nue_2 = tfs->make<TH1F>("h_StartWire_diff_nue_2","Start wire difference between cluster and beam flash",500,-2000,2000);
  h_EndWire_diff_nue_2 = tfs->make<TH1F>("h_EndWire_diff_nue_2","End wire difference between cluster and beam flash",500,-2000,2000);
  h_StartWire_diff_numu_2 = tfs->make<TH1F>("h_StartWire_diff_numu_2","Start wire difference between cluster and beam flash",100,-2000,2000);
  h_EndWire_diff_numu_2 = tfs->make<TH1F>("h_EndWire_diff_numu_2","End wire difference between cluster and beam flash",100,-2000,2000);
  h_StartWire_diff_cosmic_2 = tfs->make<TH1F>("h_StartWire_diff_cosmic_2","Start wire difference between cluster and beam flash",100,-2000,2000);
  h_EndWire_diff_cosmic_2 = tfs->make<TH1F>("h_EndWire_diff_cosmic_2","End wire difference between cluster and beam flash",100,-2000,2000);
  h_CenterWire_diff_nue_2= tfs->make<TH1F>("h_CenterWire_diff_nue_2","center wire difference between cluster and beam flash",500,0,4000);
  h_CenterWire_diff_numu_2= tfs->make<TH1F>("h_CenterWire_diff_numu_2","center wire difference between cluster and beam flash",100,0,4000);
  h_CenterWire_diff_cosmic_2= tfs->make<TH1F>("h_CenterWire_diff_cosmic_2","center wire difference between cluster and beam flash",100,0,4000);
  h_FlashClusterDistance_nue_2= tfs->make<TH1F>("h_FlashClusterDistance_nue_2","flash cluster distance of nue plane 0",100,0,4000);
  h_FlashClusterDistance_cosmic_2= tfs->make<TH1F>("h_FlashClusterDistance_cosmic_2","flash cluster distance of cosmic plane0",100,0,4000);
  h_nhits_nue_2 = tfs->make<TH1F>("h_nhits_nue_2","number of the hits of the e- cluster with FlashClusterDistance < 500",1000,0,1000);
  h_nhits_cosmic_2 = tfs->make<TH1F>("h_nhits_cosmic_2","number of the hits of the e- cluster with FlashClusterDistance < 500",1000,0,1000);


  h_Stat = tfs->make<TH1F>("h_Stat","statistics of each cut",20,0.,20.);
  h_ophits_per_flash = tfs->make<TH1F>("h_ophits_per_flash","OpHits per Flash;N_{optical hits};Events / bin",20,-0.5,19.5);
  h_ophits_per_flash_2pe = tfs->make<TH1F>("h_ophits_per_flash_2pe","OpHits (> 2 PE) per Flash;N_{optical hits};Events / bin",20,-0.5,19.5);
}

void ub::LowLevelNueFilter::GetTruthInfo(std::vector<art::Ptr<recob::Hit>> hits, int& origin, int& pdgcode){
  
  art::ServiceHandle<cheat::BackTracker> bt;

  origin = 0;
  pdgcode = 0;

  std::map<int,double> trkide;
  for(size_t h = 0; h < hits.size(); ++h){
    art::Ptr<recob::Hit> hit = hits[h];
    //std::vector<sim::IDE> ides;
    std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
    
    for(size_t e = 0; e < TrackIDs.size(); ++e){
      trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
    }
  }
  // Work out which IDE despoited the most charge in the hit if there was more than one.
  double maxe = -1;
  double tote = 0;
  int Trackid = 0;
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      //if(pfPartIdx < max_pfparticles) origin=ii->first;
      Trackid=ii->first;
    }
  }
  
  const simb::MCParticle* particle=bt->TrackIDToParticle(Trackid);	    
  if(particle){
    pdgcode = particle->PdgCode();
    origin = bt->TrackIDToMCTruth(Trackid)->Origin();
  }
}
    
bool ub::LowLevelNueFilter::inFV(double x, double y, double z)
{
  geo::GeometryCore const* fGeometry(lar::providerFrom<geo::Geometry>());
  double fDistToEdgeX             = fGeometry->DetHalfWidth()   - 20.;
  double fDistToEdgeY             = fGeometry->DetHalfHeight()  - 20.;
  double fDistToEdgeZ             = fGeometry->DetLength() / 2. - 10.;
  double distInX = x - fGeometry->DetHalfWidth();
  double distInY = y;
  double distInZ = z - 0.5 * fGeometry->DetLength();
  
  if (std::abs(distInX) < fDistToEdgeX && std::abs(distInY) < fDistToEdgeY && std::abs(distInZ) < fDistToEdgeZ) return true;
  
  return false;
}
DEFINE_ART_MODULE(ub::LowLevelNueFilter)

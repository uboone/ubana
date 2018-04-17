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
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"


#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

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

  TH1F* h_NflashInBeamSpill;
  TH1F* h_flashWidth0;
  TH1F* h_flashWidth1;
  TH1F* h_flashWidth2;
  TH1F* h_flashZWidth;

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
  TH1F* h_BNBCosmicStat;
  TH1F* h_DataStat;
  TH1F* h_ophits_per_flash;
  TH1F* h_ophits_per_flash_2pe;
  //TH1F* h_ClusterNegPdgcode;

  TH1F* h_trueNuE_Energy_FidVol;
  TH1F* h_trueLep_Mom_FidVol;
  TH1F* h_trueNuE_Energy_FidVol_10;
  TH1F* h_trueLep_Mom_FidVol_10;
  TH1F* h_trueNuE_Energy_FidVol_11;
  TH1F* h_trueLep_Mom_FidVol_11;
  TH1F* h_trueNuE_Energy_FidVol_12;
  TH1F* h_trueLep_Mom_FidVol_12;
  TH1F* h_trueNuE_Energy_FidVol_13;
  TH1F* h_trueLep_Mom_FidVol_13;
  TH1F* h_trueNuE_Energy_FidVol_14;
  TH1F* h_trueLep_Mom_FidVol_14;
  TH1F* h_trueNuE_Energy_FidVol_15;
  TH1F* h_trueLep_Mom_FidVol_15;

  TH1F* h_pdgcodeNuOthers0;
  TH1F* h_pdgcodeNuOthers1;
  TH1F* h_pdgcodeNuOthers2;
  TH1F* h_pdgcodeCosmic0;
  TH1F* h_pdgcodeCosmic1;
  TH1F* h_pdgcodeCosmic2;

  TH1F* h_StartAngleNue2;
  TH1F* h_StartOpenAngleNue2;
  TH1F* h_EndAngleNue2;
  TH1F* h_EndOpenAngleNue2;
  TH1F* h_StartAngleNuOthers2;
  TH1F* h_StartOpenAngleNuOthers2;
  TH1F* h_EndAngleNuOthers2;
  TH1F* h_EndOpenAngleNuOthers2;
  TH1F* h_StartAngleCosmic2;
  TH1F* h_StartOpenAngleCosmic2;
  TH1F* h_EndAngleCosmic2;
  TH1F* h_EndOpenAngleCosmic2;

  TH2F* h_clusterIntegralEnergy_Nue_0;
  TH2F* h_clusterSummedADCEnergy_Nue_0;
  TH2F* h_clusterIntegralEnergy_NuOthers_0;
  TH2F* h_clusterSummedADCEnergy_NuOthers_0;
  TH2F* h_clusterIntegralEnergy_Cosmic_0;
  TH2F* h_clusterSummedADCEnergy_Cosmic_0;
  TH2F* h_clusterIntegralEnergy_Nue_1;
  TH2F* h_clusterSummedADCEnergy_Nue_1;
  TH2F* h_clusterIntegralEnergy_NuOthers_1;
  TH2F* h_clusterSummedADCEnergy_NuOthers_1;
  TH2F* h_clusterIntegralEnergy_Cosmic_1;
  TH2F* h_clusterSummedADCEnergy_Cosmic_1;
  TH2F* h_clusterIntegralEnergy_Nue_2;
  TH2F* h_clusterSummedADCEnergy_Nue_2;
  TH2F* h_clusterIntegralEnergy_NuOthers_2;
  TH2F* h_clusterSummedADCEnergy_NuOthers_2;
  TH2F* h_clusterIntegralEnergy_Cosmic_2;
  TH2F* h_clusterSummedADCEnergy_Cosmic_2;

  TH1F* h_MatchedClusterFromDifferentPlaneTimeDiff_nue;
  TH1F* h_MatchedClusterFromDifferentPlaneTimeDiff_NuOthers;
  TH1F* h_MatchedClusterFromDifferentPlaneTimeDiff_Cosmic;
  TH1F* h_FlashMatchedClusterSize;
  TH1F* h_FlashMatchedClusterSize_Nue;
  TH1F* h_FlashMatchedClusterSize_NuOthers;
  TH1F* h_FlashMatchedClusterSize_Cosmic;

  TH2F* h_FlashCluster_Angle_nue;
  TH2F* h_FlashCluster_NHits_nue;
  TH2F* h_FlashCluster_Width_nue;
  TH2F* h_FlashCluster_StartTime_nue;
  TH2F* h_FlashCluster_Angle_NuOthers;
  TH2F* h_FlashCluster_NHits_NuOthers;
  TH2F* h_FlashCluster_Width_NuOthers;
  TH2F* h_FlashCluster_StartTime_NuOthers;
  TH2F* h_FlashCluster_Angle_Cosmic;
  TH2F* h_FlashCluster_NHits_Cosmic;
  TH2F* h_FlashCluster_Width_Cosmic;
  TH2F* h_FlashCluster_StartTime_Cosmic;

  TH2F* h_AngleNHits_nue_0;
  TH2F* h_AngleNHits_nue_1;
  TH2F* h_AngleNHits_nue_2;
  TH2F* h_AngleNHits_NuOthers_0;
  TH2F* h_AngleNHits_NuOthers_1;
  TH2F* h_AngleNHits_NuOthers_2;
  TH2F* h_AngleNHits_Cosmic_0;
  TH2F* h_AngleNHits_Cosmic_1;
  TH2F* h_AngleNHits_Cosmic_2;

  //data plots
  TH1F* h_StartWire_diff_data_0;
  TH1F* h_EndWire_diff_data_0;
  TH2F* h_AngleNHits_data_0;
  TH1F* h_StartWire_diff_data_1;
  TH1F* h_EndWire_diff_data_1;
  TH2F* h_AngleNHits_data_1;
  TH1F* h_StartWire_diff_data_2;
  TH1F* h_EndWire_diff_data_2;
  TH2F* h_AngleNHits_data_2;
  TH1F* h_CenterWire_diff_data_0;
  TH1F* h_CenterWire_diff_data_1;
  TH1F* h_CenterWire_diff_data_2;
  TH1F* h_StartAngle_data_goodcluster_0;
  TH1F* h_NHits_data_goodcluster_0;
  TH1F* h_TotalEnergy_data_goodcluster_0;
  TH1F* h_StartAngle_data_goodcluster_1;
  TH1F* h_NHits_data_goodcluster_1;
  TH1F* h_TotalEnergy_data_goodcluster_1;
  TH1F* h_StartAngle_data_goodcluster_2;
  TH1F* h_NHits_data_goodcluster_2;
  TH1F* h_TotalEnergy_data_goodcluster_2;
  
  TH1F* h_TwoCluster_TimeDiff;
  TH1F* h_TwoCluster_TimeDiff_nue;
  TH1F* h_TwoCluster_TimeDiff_Cosmic;
  TH1F* h_NumberOfGoodCluster_data;
  TH1F* h_NumberOfGoodCluster_nue;
  TH1F* h_NumberOfGoodCluster_NuOthers;
  TH1F* h_NumberOfGoodCluster_Cosmic;

  TH2F* h_ClusterCenterVsWidth_UplaneSpecialCut;
  TH2F* h_ClusterCenterVsWidth_VplaneSpecialCut;
  TH2F* h_ClusterCenterVsWidth_YplaneSpecialCut;

  TEfficiency* h_Eff_Ev_10 = 0;
  TEfficiency* h_Eff_Ev_11 = 0;
  TEfficiency* h_Eff_Ev_12 = 0;
  TEfficiency* h_Eff_Ev_13 = 0;
  TEfficiency* h_Eff_lm_10 = 0;
  TEfficiency* h_Eff_lm_11 = 0;
  TEfficiency* h_Eff_lm_12 = 0;
  TEfficiency* h_Eff_lm_13 = 0;
 
  

  void GetTruthInfo(std::vector<art::Ptr<recob::Hit>> hits, int& origin, int& pdgcode);
  bool inFV(double x, double y, double z);
  void doEfficiencies();

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
  bool ClusterFlashMatch0 = false;
  bool ClusterFlashMatch1=false;
  bool ClusterFlashMatch2=false;
  bool NuMuCC = false;
  bool NuMuNC = false;
  bool NuENC = false;
  bool NuECC = false;
  bool NuECCVtxInFV = false;
  bool BNBCosmic=true;
  bool MCSample=true;
  bool DataSample=false;
  bool ClusterFlashMatchNueEvent0=false;
  bool ClusterFlashMatchNueEvent1=false;
  bool ClusterFlashMatchNueEvent2=false;
  bool ClusterFlashMatchNuOthersEvent0=false;
  bool ClusterFlashMatchNuOthersEvent1=false;
  bool ClusterFlashMatchNuOthersEvent2=false;
  bool ClusterFlashMatchCosmicEvent0=false;
  bool ClusterFlashMatchCosmicEvent1=false;
  bool ClusterFlashMatchCosmicEvent2=false;
  bool DropTwoClusterEvent=false;
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
  std::cout<<"MCSample is "<<MCSample<<std::endl;
  if(MCSample == true){
    std::cout<<"MCtruth size is "<<MCtruth_vec.size()<<std::endl;
  }
    int CCNC=999;//cc=0, nc=1
    int mode=999;//QE=0, RES=1, DIS=2, Coherent =3, MEC=10
    int NuPDG=999; //nue=12, nuebar=-12,numu=14,numubar=-14
    float enutruth=0; //in GeV
    float nuvtxx=0;
    float nuvtxy=0;
    float nuvtxz=0;
    float lep_mom_truth=0;//in GeV
    if(MCSample){
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
      
      std::cout<<"NuPDG = "<<NuPDG<<" mode = "<<mode<<" CCNC = "<<CCNC<<" nue energy = "<<enutruth<<" nu Vertex = ["<<nuvtxx<<" , "<<nuvtxy<<" , "<<nuvtxz<<"]"<<" lepton momentum = "<<lep_mom_truth<<std::endl;}
    //h_Stat->AddBinContent(1);
    h_Stat->SetBinContent(1,h_Stat->GetBinContent(1)+1);
    if(DataSample){
      h_DataStat->SetBinContent(1,h_DataStat->GetBinContent(1)+1);
    }
    if(MCSample){
      if(inFV(nuvtxx,nuvtxy,nuvtxz) && CCNC==0 && NuPDG==12){
	//h_Stat->AddBinContent(2);
	h_Stat->SetBinContent(2,h_Stat->GetBinContent(2)+1);
	TruthSignal = true;
	h_trueNuE_Energy_FidVol->Fill(enutruth*1000.);
	h_trueLep_Mom_FidVol->Fill(lep_mom_truth*1000.);
      }
      
      //BNB + Cosmic MC sample
      if(BNBCosmic == true){
	h_BNBCosmicStat->SetBinContent(1,h_BNBCosmicStat->GetBinContent(1)+1);
	if(NuPDG==14 || NuPDG==-14) //mumu events
	  {if(CCNC==0)//charge current
	      {h_BNBCosmicStat->SetBinContent(2,h_BNBCosmicStat->GetBinContent(2)+1);
		NuMuCC=true;
	      }
	    if(CCNC==1)
	      {h_BNBCosmicStat->SetBinContent(3,h_BNBCosmicStat->GetBinContent(3)+1);
		NuMuNC=true;
	      }
	  }
	if(NuPDG==12 || NuPDG==-12)
	  {
	    if(CCNC==1)//NC
	      {h_BNBCosmicStat->SetBinContent(4,h_BNBCosmicStat->GetBinContent(4)+1);
		NuENC=true;
	      }
	    if(CCNC==0)//CC
	      {h_BNBCosmicStat->SetBinContent(5,h_BNBCosmicStat->GetBinContent(5)+1);
		NuECC=true;
		if(inFV(nuvtxx,nuvtxy,nuvtxz))
		  {h_BNBCosmicStat->SetBinContent(6,h_BNBCosmicStat->GetBinContent(6)+1);
		    NuECCVtxInFV=true;
		  }
	      }
	  }
      }//BNBCosmic
    }//MCSample
  

    //Flashes
    art::Handle<std::vector<recob::OpFlash> > opflashHandle;
    std::vector<art::Ptr<recob::OpFlash> > OpFlash_vec;
    if (e.getByLabel(fOpFlash_tag, opflashHandle))
      art::fill_ptr_vector(OpFlash_vec, opflashHandle);
    
    std::vector<TVector3> FlashWire;
    for (auto & flash : OpFlash_vec){
      if(flash->Time()<3 || flash->Time()>5)
	continue;
      //float wire0 = FLT_MAX;
      //float wire1 = FLT_MIN;
      h_flashZWidth->Fill(flash->ZWidth());
      //Find the 4 corners and convert them to wire numbers
      std::vector<TVector3> points;
      points.push_back(TVector3(0, flash->YCenter()-flash->YWidth(), flash->ZCenter()-flash->ZWidth()));
      points.push_back(TVector3(0, flash->YCenter()-flash->YWidth(), flash->ZCenter()+flash->ZWidth()));
      points.push_back(TVector3(0, flash->YCenter()+flash->YWidth(), flash->ZCenter()-flash->ZWidth()));
      points.push_back(TVector3(0, flash->YCenter()+flash->YWidth(), flash->ZCenter()+flash->ZWidth()));
      for (int pl = 0; pl <3; ++pl){
	float wire0 = FLT_MAX;
	float wire1 = FLT_MIN;
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
	auto WireDiff = wire1-wire0;
	FlashWire.push_back(TVector3(pl,wire0+WireDiff/4.,wire1-WireDiff/4.));
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
    if(DataSample && FlashInBeamWindow == true)
      {h_DataStat->SetBinContent(2,h_DataStat->GetBinContent(2)+1);}
    // Clusters
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    std::vector<art::Ptr<recob::Cluster> > Cluster_vec;
    if (e.getByLabel(fCluster_tag, clusterHandle))
      art::fill_ptr_vector(Cluster_vec, clusterHandle);

    art::FindManyP<recob::Hit>  fmhc(clusterHandle, e, fCluster_tag);
 
    std::vector<art::Ptr<recob::Cluster>> FlashMatchedClusterVec;
    for( auto const& cluster : Cluster_vec){
      auto ID = cluster->ID();
      auto plane = cluster->Plane().Plane;
      //cout<<"cluster ID is "<<ID<<" in plane "<<plane<<endl;
      if(DataSample)
	{
	  if(plane==0){
	    h_DataStat->SetBinContent(3,h_DataStat->GetBinContent(3)+1);}
	  if(plane==1){
	    h_DataStat->SetBinContent(4,h_DataStat->GetBinContent(4)+1);}
	  if(plane==2){
	    h_DataStat->SetBinContent(5,h_DataStat->GetBinContent(5)+1);}
	}
      if(((int)ID)>=0)
	continue;
      bool ClusterFlashMatchNue0=false;
      bool ClusterFlashMatchNue1=false;
      bool ClusterFlashMatchNue2=false;
      bool ClusterFlashMatchNuOthers0 = false;
      bool ClusterFlashMatchNuOthers1 = false;
      bool ClusterFlashMatchNuOthers2 = false;
      bool ClusterFlashMatchCosmic0 = false;
      bool ClusterFlashMatchCosmic1 = false;
      bool ClusterFlashMatchCosmic2 = false;
      if(((int)ID)<0)
	{
	  auto nHits = cluster->NHits();
	  auto StartWire = cluster->StartWire();
	  //auto StartTick = cluster->StartTick();
	  auto EndWire = cluster->EndWire();
	 
	  // auto TotalCharge = cluster->Integral();
          int origin = 0;
          int pdgcode = 0;
	  if(DataSample){
	    if(plane==0){
	      h_DataStat->SetBinContent(6,h_DataStat->GetBinContent(6)+1);}
	    if(plane==1){
	      h_DataStat->SetBinContent(7,h_DataStat->GetBinContent(7)+1);}
	    if(plane==2){
	      h_DataStat->SetBinContent(8,h_DataStat->GetBinContent(8)+1);}
	  }

	  if(MCSample){
	    GetTruthInfo(fmhc.at(cluster.key()), origin, pdgcode);}
          //massive clean
	  //std::cout<<"ClusterID="<<ID<<" plane="<<plane<<"nhits = "<<nHits<<" startwire = "<<StartWire<<" starttick = "<<StartTick<<" End wire "<<EndWire<< " origin = "<<origin<<" pdgcode = "<<pdgcode<<std::endl;
	  if(StartWire>=EndWire)
	    {StartWire = cluster->EndWire();
	      EndWire = cluster->StartWire();}
	  //does the cluster match with the nu flash?
	  for(size_t i=0; i<FlashWire.size();i++)
	    {
	      if(FlashWire[i].X()!=plane)
		continue;
	      //if(FlashWire[i].X()==plane)
	      //{
	      
		  auto ClusterFlashStartDiff = FlashWire[i].Y()-StartWire;
		  auto ClusterFlashEndDiff = FlashWire[i].Z()-EndWire;
		  if(plane==0){
		    h_StartWire_diff_data_0->Fill(ClusterFlashStartDiff);
		    h_EndWire_diff_data_0->Fill(ClusterFlashEndDiff);
		    h_AngleNHits_data_0->Fill(cluster->StartAngle(),nHits);
		  }
		  if(plane==1){
		    h_StartWire_diff_data_1->Fill(ClusterFlashStartDiff);
		    h_EndWire_diff_data_1->Fill(ClusterFlashEndDiff);
		    h_AngleNHits_data_1->Fill(cluster->StartAngle(),nHits);
		  }
		  if(plane==2){
		    h_StartWire_diff_data_2->Fill(ClusterFlashStartDiff);
		    h_EndWire_diff_data_2->Fill(ClusterFlashEndDiff);
		    h_AngleNHits_data_2->Fill(cluster->StartAngle(),nHits);
		  }
		  if(MCSample){
		    if(origin==1 && pdgcode == 11)
		      {
			if(plane==0){
			  h_StartWire_diff_nue_0->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_nue_0->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_nue_0->Fill(cluster->StartAngle(),nHits);
			}
			if(plane==1){
			  h_StartWire_diff_nue_1->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_nue_1->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_nue_1->Fill(cluster->StartAngle(),nHits);
			}
			if(plane==2){
			  //massive clean
			  //std::cout<<"plane 2 nue cluster start wire is "<<StartWire<<" End wire is "<<EndWire<<" start tick is "<<StartTick<<" nhits"<<nHits<<std::endl;
			  h_StartWire_diff_nue_2->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_nue_2->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_nue_2->Fill(cluster->StartAngle(),nHits);
			}
			
		      }
		    if(origin==1 && pdgcode != 11)
		      {
			if(plane==0){
			  h_StartWire_diff_numu_0->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_numu_0->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_NuOthers_0->Fill(cluster->StartAngle(),nHits);
			}
			if(plane==1){
			  h_StartWire_diff_numu_1->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_numu_1->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_NuOthers_1->Fill(cluster->StartAngle(),nHits);
			}
			if(plane==2){
			  h_StartWire_diff_numu_2->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_numu_2->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_NuOthers_2->Fill(cluster->StartAngle(),nHits);
			}
			
		      }
		    if(origin == 2)
		      {
			if(plane==0){
			  h_StartWire_diff_cosmic_0->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_cosmic_0->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_Cosmic_0->Fill(cluster->StartAngle(),nHits);
			}
			if(plane==1){
			  h_StartWire_diff_cosmic_1->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_cosmic_1->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_Cosmic_1->Fill(cluster->StartAngle(),nHits);
			}
			if(plane==2){
			  h_StartWire_diff_cosmic_2->Fill(ClusterFlashStartDiff);
			  h_EndWire_diff_cosmic_2->Fill(ClusterFlashEndDiff);
			  h_AngleNHits_Cosmic_2->Fill(cluster->StartAngle(),nHits);
			  //massive clean
			  //std::cout<<"plane 2 nue cluster start wire is "<<StartWire<<" End wire is "<<EndWire<<" start tick is "<<StartTick<<" nhits"<<nHits<<std::endl;
			}
		      }
		  }//MCSample == true
		  
		  auto ClusterWireCenter = (StartWire + EndWire)/2.;
		  auto FlashWireCenter = (FlashWire[i].Y() + FlashWire[i].Z())/2.;
		  auto ClusterFlashCenterDiff = std::fabs(FlashWireCenter - ClusterWireCenter);
		  double FlashClusterDistance=9999;
		  if(StartWire>FlashWire[i].Y() && StartWire<FlashWire[i].Z() && EndWire>FlashWire[i].Y() && EndWire<FlashWire[i].Z())//flash include cluster both ends
		    {FlashClusterDistance=0;}
		  if(StartWire<=FlashWire[i].Y() || EndWire>=FlashWire[i].Z())
		    {FlashClusterDistance = ClusterFlashCenterDiff;}
		  if(plane==0)
		    h_CenterWire_diff_data_0->Fill(ClusterFlashCenterDiff);
		  if(plane==1)
		    h_CenterWire_diff_data_1->Fill(ClusterFlashCenterDiff);
		  if (plane==2)
		    h_CenterWire_diff_data_2->Fill(ClusterFlashCenterDiff);
		  if(MCSample){
		    if(origin==1 && pdgcode == 11)
		      {
			if(plane==0){
			  h_CenterWire_diff_nue_0->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_nue_0->Fill(FlashClusterDistance);
			  //massive clean
			  /*if(FlashClusterDistance<500)
			    {
			      h_nhits_nue_0->Fill(nHits);
			      std::cout<<"Event Dis nue plane 0 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			      }*/
			}
			if(plane==1){
			  h_CenterWire_diff_nue_1->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_nue_1->Fill(FlashClusterDistance);
			  //massive clean
			  /*if(FlashClusterDistance<500)
			    {h_nhits_nue_1->Fill(nHits);
			      std::cout<<"Event Dis nue plane 1 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			      }*/
			}
			if(plane==2){
			  h_CenterWire_diff_nue_2->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_nue_2->Fill(FlashClusterDistance);
			  //massive clean
			  /*if(FlashClusterDistance<500)
			    {h_nhits_nue_2->Fill(nHits);
			      std::cout<<"Event Dis nue plane 2 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			      }*/
			}
		      }
		    
		    if(origin==1 && pdgcode != 11)
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
			  //massive clean
			  /*if(FlashClusterDistance<500)
			    {h_nhits_cosmic_0->Fill(nHits);
			      std::cout<<"Event Dis cosmic plane 0 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			      }*/
			}
			if(plane==1){
			  h_CenterWire_diff_cosmic_1->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_cosmic_1->Fill(FlashClusterDistance);
			  //massive clean
			  /*if(FlashClusterDistance<500)
			    {h_nhits_cosmic_1->Fill(nHits);
			      std::cout<<"Event Dis cosmic plane 1 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;
			      }*/
			}
			if(plane==2){
			  h_CenterWire_diff_cosmic_2->Fill(ClusterFlashCenterDiff);
			  h_FlashClusterDistance_cosmic_2->Fill(FlashClusterDistance);
			  //massive clean
			  /*if(FlashClusterDistance<500)
			  {h_nhits_cosmic_2->Fill(nHits);
			  std::cout<<"Event Dis cosmic plane 2 cluster start wire "<<StartWire<<" end wire "<<EndWire<<"Start tick"<<StartTick<<" nhits"<<nHits<<std::endl;			    
			    }*/
			}
		      }
		  }//if MCSample == true
		  //if(plane==0 && (ClusterFlashStartDiff>-700 && ClusterFlashStartDiff<150) && (ClusterFlashEndDiff>-200 && ClusterFlashEndDiff<700) && ClusterFlashCenterDiff<350 && nHits>=35)//this is using full flash width
		  if(plane==0 && (ClusterFlashStartDiff>-550 && ClusterFlashStartDiff<300) && (ClusterFlashEndDiff>-350 && ClusterFlashEndDiff<520) && ClusterFlashCenterDiff<400 && nHits>=30)
		    {
		      ClusterFlashMatch0=true;
		      h_ClusterCenterVsWidth_UplaneSpecialCut->Fill(ClusterWireCenter,EndWire-StartWire);
		      if(DataSample)
			{h_DataStat->SetBinContent(9,h_DataStat->GetBinContent(9)+1);}
		      FlashMatchedClusterVec.push_back(cluster);
		      h_StartAngle_data_goodcluster_0->Fill(cluster->StartAngle());
		      h_NHits_data_goodcluster_0->Fill(nHits);
		      h_TotalEnergy_data_goodcluster_0->Fill(cluster->Integral()*0.007867);//should use the data calibration constant 
		      if(MCSample){
			if(origin==1 && pdgcode ==11)
			  {ClusterFlashMatchNue0 = true;
			    std::cout<<"plane 0 cluster integral is "<<cluster->Integral()<<" summed ADC is "<<cluster->SummedADC()<<"lepton momentym is "<<lep_mom_truth<<std::endl;
			    h_clusterIntegralEnergy_Nue_0->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.); //*200 e/ADC * 23.6eV/e /10^6(MeV) /0.6 (recombination)
			    h_clusterSummedADCEnergy_Nue_0->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			    
			  }
			if(origin==1 && pdgcode !=11)
			  {ClusterFlashMatchNuOthers0 =true;
			    h_clusterIntegralEnergy_NuOthers_0->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_NuOthers_0->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
			if(origin==2)
			  {ClusterFlashMatchCosmic0 = true;
			    h_clusterIntegralEnergy_Cosmic_0->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_Cosmic_0->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
		      }//MC sample
		      
		    }
		  //if(plane==1 && (ClusterFlashStartDiff>-700 && ClusterFlashStartDiff<150) && (ClusterFlashEndDiff>-200 && ClusterFlashEndDiff<700) && ClusterFlashCenterDiff<350 && nHits>=35)//full flash width
		  
		  bool VPlaneNoise = false;
		  if(plane==1){
		    auto VClusterFullWidth = EndWire-StartWire;
		    if(ClusterWireCenter>=1280 && ClusterWireCenter<1300 && VClusterFullWidth<=30)
		      {VPlaneNoise=true;}
		    if((ClusterWireCenter>=1300 && ClusterWireCenter<1310) && (VClusterFullWidth<=50 && VClusterFullWidth>10))
		      {VPlaneNoise=true;}
		    if((ClusterWireCenter>=1310 && ClusterWireCenter<1320) && (VClusterFullWidth<=70 && VClusterFullWidth>20))
		      {VPlaneNoise=true;}
		    if((ClusterWireCenter>=1320 && ClusterWireCenter<1330) && (VClusterFullWidth<=80 && VClusterFullWidth>60))
		      {VPlaneNoise=true;}
		    if((ClusterWireCenter>=1330 && ClusterWireCenter<1340) && (VClusterFullWidth<=110 && VClusterFullWidth>80))
		      {VPlaneNoise=true;}
		    if((ClusterWireCenter>=1340 && ClusterWireCenter<1350) && (VClusterFullWidth<=120 && VClusterFullWidth>100))
		      {VPlaneNoise=true;}
		    if((ClusterWireCenter>=1350 && ClusterWireCenter<1360) && (VClusterFullWidth<=140 && VClusterFullWidth>120))
		      {VPlaneNoise=true;}		   		       
		  }
		  if(plane==1 && (ClusterFlashStartDiff>-550 && ClusterFlashStartDiff<300) && (ClusterFlashEndDiff>-350 && ClusterFlashEndDiff<520) && ClusterFlashCenterDiff<400 && nHits>=30 && VPlaneNoise==false)
		    {
		      ClusterFlashMatch1=true;
		      h_ClusterCenterVsWidth_VplaneSpecialCut->Fill(ClusterWireCenter,EndWire-StartWire);
		      
		      if(DataSample)
			{h_DataStat->SetBinContent(10,h_DataStat->GetBinContent(10)+1);}
		      FlashMatchedClusterVec.push_back(cluster);
		      h_StartAngle_data_goodcluster_1->Fill(cluster->StartAngle());
		      h_NHits_data_goodcluster_1->Fill(nHits);
		      h_TotalEnergy_data_goodcluster_1->Fill(cluster->Integral()*0.007867);//should use the data calibration constant 
		      if(MCSample){
			if(origin==1 && pdgcode ==11)
			  {
			    ClusterFlashMatchNue1 = true;
			    std::cout<<"plane 1 cluster integral is "<<cluster->Integral()<<" summed ADC is "<<cluster->SummedADC()<<"lepton momentym is "<<lep_mom_truth<<std::endl;
			    h_clusterIntegralEnergy_Nue_1->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_Nue_1->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
			if(origin==1 && pdgcode !=11)
			  {ClusterFlashMatchNuOthers1 =true;
			    h_clusterIntegralEnergy_NuOthers_1->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_NuOthers_1->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
			if(origin==2)
			  {ClusterFlashMatchCosmic1 = true;
			    h_clusterIntegralEnergy_Cosmic_1->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_Cosmic_1->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
		      }//MC sample
		    }	
		  //if(plane==2 && (ClusterFlashStartDiff>-800 && ClusterFlashStartDiff<100)&& (ClusterFlashEndDiff>-150 && ClusterFlashEndDiff<800) && ClusterFlashCenterDiff<400 && (std::fabs(cluster->StartAngle())>1.6 || std::fabs(cluster->StartAngle())<1.3) && nHits>=35)//full width
		  if(plane==2 && (ClusterFlashStartDiff>-580 && ClusterFlashStartDiff<280)&& (ClusterFlashEndDiff>-360 && ClusterFlashEndDiff<500) && ClusterFlashCenterDiff<400 && (std::fabs(cluster->StartAngle())>1.6 || std::fabs(cluster->StartAngle())<1.3) && nHits>=30)
		    {
		      ClusterFlashMatch2=true;
		      h_ClusterCenterVsWidth_YplaneSpecialCut->Fill(ClusterWireCenter,EndWire-StartWire);
		      if(DataSample)
			{h_DataStat->SetBinContent(11,h_DataStat->GetBinContent(11)+1);}
		      FlashMatchedClusterVec.push_back(cluster);
		      h_StartAngle_data_goodcluster_2->Fill(cluster->StartAngle());
		      h_NHits_data_goodcluster_2->Fill(nHits);
		      h_TotalEnergy_data_goodcluster_2->Fill(cluster->Integral()*0.007867);//should use the data calibration constant 
		      if(MCSample){
			if(origin==1&&pdgcode==11)
			  {ClusterFlashMatchNue2 = true;
			    std::cout<<"plane 2 cluster integral is "<<cluster->Integral()<<" summed ADC is "<<cluster->SummedADC()<<"lepton momentym is "<<lep_mom_truth<<std::endl;
			    h_StartAngleNue2->Fill(cluster->StartAngle());
			    h_EndAngleNue2->Fill(cluster->EndAngle());
			    h_StartOpenAngleNue2->Fill(cluster->StartOpeningAngle());
			    h_EndOpenAngleNue2->Fill(cluster->EndOpeningAngle());
			    h_clusterIntegralEnergy_Nue_2->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_Nue_2->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
			if(origin==1 && pdgcode !=11)
			  {ClusterFlashMatchNuOthers2 =true;
			    h_StartAngleNuOthers2->Fill(cluster->StartAngle());
			    h_EndAngleNuOthers2->Fill(cluster->EndAngle());
			    h_StartOpenAngleNuOthers2->Fill(cluster->StartOpeningAngle());
			    h_EndOpenAngleNuOthers2->Fill(cluster->EndOpeningAngle());
			    h_clusterIntegralEnergy_NuOthers_2->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_NuOthers_2->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
			if(origin==2)
			  {ClusterFlashMatchCosmic2 = true;
			    h_StartAngleCosmic2->Fill(cluster->StartAngle());
			    h_EndAngleCosmic2->Fill(cluster->EndAngle());
			    h_StartOpenAngleCosmic2->Fill(cluster->StartOpeningAngle());
			    h_EndOpenAngleCosmic2->Fill(cluster->EndOpeningAngle());
			    h_clusterIntegralEnergy_Cosmic_2->Fill(cluster->Integral()*0.007867,lep_mom_truth*1000.);
			    h_clusterSummedADCEnergy_Cosmic_2->Fill(cluster->SummedADC()*0.007867,lep_mom_truth*1000.);
			  }
		      }
		    }
		  
	    }//flash plane loop
	    
	  if(MCSample){
	    if(BNBCosmic==true && FlashInBeamWindow == true)//BNB Cosmic sample 
	      {
		//origin 1 pdgcode 11 cluster matched with flash
		if(ClusterFlashMatchNue0) //plane 0
		  {
		    if(NuMuCC || NuMuNC || NuENC) //background
		      {h_BNBCosmicStat->SetBinContent(7,h_BNBCosmicStat->GetBinContent(7)+1);}
		    if(NuECC && !NuECCVtxInFV)
		      {h_BNBCosmicStat->SetBinContent(8,h_BNBCosmicStat->GetBinContent(8)+1);}
		    if(NuECCVtxInFV)//signal
		      {h_BNBCosmicStat->SetBinContent(9,h_BNBCosmicStat->GetBinContent(9)+1);
			ClusterFlashMatchNueEvent0=true;
		      }
		  }
		if(ClusterFlashMatchNue1)//plane 1
		  {
		    if(NuMuCC || NuMuNC || NuENC)
		      {h_BNBCosmicStat->SetBinContent(10,h_BNBCosmicStat->GetBinContent(10)+1);}
		    if(NuECC && !NuECCVtxInFV)
		      {h_BNBCosmicStat->SetBinContent(11,h_BNBCosmicStat->GetBinContent(11)+1);}
		    if(NuECCVtxInFV)
		      {h_BNBCosmicStat->SetBinContent(12,h_BNBCosmicStat->GetBinContent(12)+1);
			ClusterFlashMatchNueEvent1=true;
		      }
		  }
		if(ClusterFlashMatchNue2) //plane 2
		  {
		    if(NuMuCC || NuMuNC || NuENC)
		      {h_BNBCosmicStat->SetBinContent(13,h_BNBCosmicStat->GetBinContent(13)+1);}
		    if(NuECC && !NuECCVtxInFV)
		      {h_BNBCosmicStat->SetBinContent(14,h_BNBCosmicStat->GetBinContent(14)+1);}
		    if(NuECCVtxInFV)
		      {h_BNBCosmicStat->SetBinContent(15,h_BNBCosmicStat->GetBinContent(15)+1);
			ClusterFlashMatchNueEvent2=true;
		      }
		  }
		//nu interaction but not from electron cluster matched with flash
		if(ClusterFlashMatchNuOthers0){
		  if(NuMuCC)
		    {h_BNBCosmicStat->SetBinContent(16,h_BNBCosmicStat->GetBinContent(16)+1);		    
		    }
		  if(NuMuNC)
		    {h_BNBCosmicStat->SetBinContent(17,h_BNBCosmicStat->GetBinContent(17)+1);}
		  h_pdgcodeNuOthers0->Fill(pdgcode);	
		  ClusterFlashMatchNuOthersEvent0=true;
		}
		if(ClusterFlashMatchNuOthers1){
		  if(NuMuCC)
		    {h_BNBCosmicStat->SetBinContent(18,h_BNBCosmicStat->GetBinContent(18)+1);		    
		    }
		  if(NuMuNC)
		    {h_BNBCosmicStat->SetBinContent(19,h_BNBCosmicStat->GetBinContent(19)+1);}
		  h_pdgcodeNuOthers1->Fill(pdgcode);
		  ClusterFlashMatchNuOthersEvent1=true;
		}
		if(ClusterFlashMatchNuOthers2){
		  if(NuMuCC)
		    {h_BNBCosmicStat->SetBinContent(20,h_BNBCosmicStat->GetBinContent(20)+1);		    
		    }
		  if(NuMuNC)
		    {h_BNBCosmicStat->SetBinContent(21,h_BNBCosmicStat->GetBinContent(21)+1);}
		  h_pdgcodeNuOthers2->Fill(pdgcode);
		  ClusterFlashMatchNuOthersEvent0=true;
		}
		//cosmic origin 2 cluster matched with flash
		if(ClusterFlashMatchCosmic0) //plane 0
		  {
		    h_BNBCosmicStat->SetBinContent(22,h_BNBCosmicStat->GetBinContent(22)+1);
		    h_pdgcodeCosmic0->Fill(pdgcode);
		    ClusterFlashMatchCosmicEvent0=true;
		  }
		if(ClusterFlashMatchCosmic1)//plane 1
		  {h_BNBCosmicStat->SetBinContent(23,h_BNBCosmicStat->GetBinContent(23)+1);
		    h_pdgcodeCosmic1->Fill(pdgcode);
		    ClusterFlashMatchCosmicEvent1=true;
		  }
		if(ClusterFlashMatchCosmic2) //plane 2
		  {
		    h_BNBCosmicStat->SetBinContent(24,h_BNBCosmicStat->GetBinContent(24)+1);
		    h_pdgcodeCosmic2->Fill(pdgcode);
		    ClusterFlashMatchCosmicEvent2=true;
		  }
		
	      } //end of the BNB+Cosmic sample
	  
	  //this chunk of code was test code for the BNB nue intrinsic + cosmic sample
	  if(FlashInBeamWindow==true /*&& TruthSignal==true*/)
	    {
	      if(ClusterFlashMatchNue0==true) //pass flash cluster matching cut in plane 0 && bt origin is 1 (neutrino) and pdgcode is 11(electron)
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
	    }//end of bnb intrinsic nue + cosmic sample
	  }//end of MCSample==true
    	}//ID<0, //origin = 1(beam neutrino), 2 (cosmic),3(supernova),4(single particle) pdgcode = 11(e-), 13(mu-),22(gamma) 
    }//cluster loop

    //in the event loop test bnb intrinsic nue + cosmic sample
    
    std::cout<<"FlashMatchedClusterVecSize is "<<FlashMatchedClusterVec.size()<<std::endl;   

    //choose the best e- cluster candidate for each plane
    //special case when there are two clusters from two planes, if time don't match, set the drop event flag to true;
    if(FlashMatchedClusterVec.size()==2)
      {
	auto cluster0 = FlashMatchedClusterVec.at(0);
	auto cluster1 = FlashMatchedClusterVec.at(1);
	auto plane0 = cluster0->Plane().Plane;
	auto plane1 = cluster1->Plane().Plane;
	auto StartTime0 = cluster0->StartTick();
	auto StartTime1 = cluster1->StartTick();
	auto TimeDiffTwoClusters = std::fabs(StartTime0-StartTime1);
	if(plane0!=plane1)
	  {
	    h_TwoCluster_TimeDiff->Fill(TimeDiffTwoClusters);
	    if(MCSample)
	      {
		int origin0=0;
		int origin1=0;
		int pdgcode0=0;
		int pdgcode1=0;
		GetTruthInfo(fmhc.at(cluster0.key()),origin0,pdgcode0);
		GetTruthInfo(fmhc.at(cluster1.key()),origin1,pdgcode1);
		if(origin0==1 && origin1==1 && pdgcode0==11 && pdgcode1==11)
		  {h_TwoCluster_TimeDiff_nue->Fill(TimeDiffTwoClusters);}
		if(origin0==2 || origin1==2)
		  {h_TwoCluster_TimeDiff_Cosmic->Fill(TimeDiffTwoClusters);}
	      }
	    if(TimeDiffTwoClusters>500)
	      DropTwoClusterEvent = true;
	  }
      }

    //massive clean maybe??
    if(FlashMatchedClusterVec.size()>1)
      {
	std::cout<<"In the function of dropping clusters!!!"<<std::endl;
	for(unsigned pl=0;pl<3;pl++){
	  for(unsigned int i=0;i<FlashMatchedClusterVec.size();i++)
	    {
	      auto iCluster = FlashMatchedClusterVec.at(i);
	      auto iID = iCluster->ID();
	      auto iPlane = iCluster->Plane().Plane;
	      if(iPlane!=pl)
		continue;
	     
	      int iOrigin=0;
	      int iPdgcode=0;
	      auto iAngle=iCluster->StartAngle();
	      auto inhits = iCluster->NHits();
	      auto iWidth = iCluster->Width();
	      auto iStartTick = iCluster->StartTick();
	      h_NumberOfGoodCluster_data->SetBinContent(pl+1,h_NumberOfGoodCluster_data->GetBinContent(pl+1)+1);
	      if(MCSample)
		{GetTruthInfo(fmhc.at(iCluster.key()), iOrigin, iPdgcode);}
	      std::cout<<"plane="<<pl<<" cluster id="<<iID<<" origin="<<iOrigin<<" pdgcode="<<iPdgcode<<" Start Angle="<<iAngle<<" nHits="<<inhits<<" width="<<iWidth<<" Start time="<<iStartTick<<std::endl;
	      if(MCSample){
		if(iOrigin==1 && iPdgcode ==11)
		  {
		    h_FlashCluster_Angle_nue->Fill(iAngle,pl);
		    h_FlashCluster_NHits_nue->Fill(inhits,pl);
		    h_FlashCluster_Width_nue->Fill(iWidth,pl);
		    h_FlashCluster_StartTime_nue->Fill(iStartTick,pl);
		    h_NumberOfGoodCluster_nue->SetBinContent(pl+1,h_NumberOfGoodCluster_nue->GetBinContent(pl+1)+1);
		  }
		if(iOrigin==1 && iPdgcode!=11)
		  {
		    h_FlashCluster_Angle_NuOthers->Fill(iAngle,pl);
		    h_FlashCluster_NHits_NuOthers->Fill(inhits,pl);
		    h_FlashCluster_Width_NuOthers->Fill(iWidth,pl);
		    h_FlashCluster_StartTime_NuOthers->Fill(iStartTick,pl); 
		    h_NumberOfGoodCluster_NuOthers->SetBinContent(pl+1,h_NumberOfGoodCluster_NuOthers->GetBinContent(pl+1)+1);
		  }
		if(iOrigin==2)
		  {
		    h_FlashCluster_Angle_Cosmic->Fill(iAngle,pl);
		    h_FlashCluster_NHits_Cosmic->Fill(inhits,pl);
		    h_FlashCluster_Width_Cosmic->Fill(iWidth,pl);
		    h_FlashCluster_StartTime_Cosmic->Fill(iStartTick,pl); 
		    h_NumberOfGoodCluster_Cosmic->SetBinContent(pl+1,h_NumberOfGoodCluster_Cosmic->GetBinContent(pl+1)+1);
		  }
	      }
	      
	    }//flash matched cluster vec loop
	}//plane
      }//flash matched cluster vec size > 1
    
    //massive clean maybe??
    if(FlashMatchedClusterVec.size()>1)
      {
	for(unsigned int i =0;i<FlashMatchedClusterVec.size()-1;i++)
	  {
	    auto iCluster = FlashMatchedClusterVec.at(i);
	    auto iPlane = iCluster->Plane().Plane;
	    auto iTime = iCluster->StartTick();
	    int iOrigin=0;
	    int iPdgcode=0;
	    if(MCSample){
	      GetTruthInfo(fmhc.at(iCluster.key()), iOrigin, iPdgcode);
	    }
	    std::cout<<"cluster "<<i<< " is at plane "<<iPlane<<" start time is "<<iTime<<std::endl;
	    if(iOrigin==0 && iPdgcode==0)
	      continue;
	    for(unsigned int j=i+1;j<FlashMatchedClusterVec.size();j++)
	      {
		auto jCluster = FlashMatchedClusterVec.at(j);
		auto jPlane = jCluster->Plane().Plane;
		auto jTime = jCluster->StartTick();
		int jOrigin=0;
		int jPdgcode=0; 
		if(MCSample){
		  GetTruthInfo(fmhc.at(jCluster.key()), jOrigin, jPdgcode);}
		std::cout<<"cluster "<<j<< " is at plane "<<jPlane<<" start time is "<<jTime<<std::endl; 
		if(jPlane==iPlane)
		  continue;
		if(jOrigin==0 && jPdgcode==0)
		  continue;
		auto TimeDiff = std::fabs(iTime-jTime);
		if(MCSample){
		  if(iOrigin==1 && jOrigin==1 && iPdgcode==11 && jPdgcode==11)
		    {h_MatchedClusterFromDifferentPlaneTimeDiff_nue->Fill(TimeDiff);}
		  if(iOrigin==1 && jOrigin==1 && (iPdgcode!=11 || jPdgcode!=11))
		    {h_MatchedClusterFromDifferentPlaneTimeDiff_NuOthers->Fill(TimeDiff);}
		  if(iOrigin==2 || jOrigin==2)
		    {
		      h_MatchedClusterFromDifferentPlaneTimeDiff_Cosmic->Fill(TimeDiff);
		    }
		}
		
	      }//cluster j
	  }//cluster i
      }//FlashMatchedClusterVec loop
    if(FlashInBeamWindow==true /*&& TruthSignal==true*/)
      {
	if(ClusterFlashMatch0==true || ClusterFlashMatch1==true || ClusterFlashMatch2==true)
	  {
	    h_Stat->SetBinContent(10,h_Stat->GetBinContent(10)+1);
	    h_trueNuE_Energy_FidVol_10->Fill(enutruth*1000.);
	    h_trueLep_Mom_FidVol_10->Fill(lep_mom_truth*1000.);
	    if(DataSample)
	      {h_DataStat->SetBinContent(12,h_DataStat->GetBinContent(12)+1);}
	    if(MCSample){
	      if((ClusterFlashMatchNueEvent0==true || ClusterFlashMatchNueEvent1==true || ClusterFlashMatchNueEvent2==true)&& (inFV(nuvtxx,nuvtxy,nuvtxz) && CCNC==0 && NuPDG==12))
		{
		  h_Stat->SetBinContent(11,h_Stat->GetBinContent(11)+1);
		  h_trueNuE_Energy_FidVol_11->Fill(enutruth*1000.);
		  h_trueLep_Mom_FidVol_11->Fill(lep_mom_truth*1000.);
		}
	    }
	  }
	if((enutruth*1000<400) && !((ClusterFlashMatch0==true && ClusterFlashMatch1==true) || (ClusterFlashMatch0==true && ClusterFlashMatch2==true) || (ClusterFlashMatch1==true && ClusterFlashMatch2==true)) && (TruthSignal==true))
	  {
	    std::cout<<"Missing low energy (less than 200MeV) CC nue events !!!???"<<std::endl; 
	  }
	if((ClusterFlashMatch0==true && ClusterFlashMatch1==true) || (ClusterFlashMatch0==true && ClusterFlashMatch2==true) || (ClusterFlashMatch1==true && ClusterFlashMatch2==true))//at least 2 flash matched clusters one from each plane
	  {
	  
	    h_Stat->SetBinContent(12,h_Stat->GetBinContent(12)+1);
	    h_trueNuE_Energy_FidVol_12->Fill(enutruth*1000.);
	    h_trueLep_Mom_FidVol_12->Fill(lep_mom_truth*1000.);
	    if(DataSample)
	      {h_DataStat->SetBinContent(13,h_DataStat->GetBinContent(13)+1);}
	    if(DropTwoClusterEvent==false)
	      {
		h_Stat->SetBinContent(14,h_Stat->GetBinContent(14)+1);
		h_trueNuE_Energy_FidVol_14->Fill(enutruth*1000.);
		h_trueLep_Mom_FidVol_14->Fill(lep_mom_truth*1000.);
		if(DataSample)
		  {h_DataStat->SetBinContent(14,h_DataStat->GetBinContent(14)+1);
		    std::cout<<"????????nue candidate????"<<std::endl;
		  }
	      }
	    if(MCSample){
	      if((ClusterFlashMatchNueEvent0==true || ClusterFlashMatchNueEvent1==true || ClusterFlashMatchNueEvent2==true) && (inFV(nuvtxx,nuvtxy,nuvtxz) && CCNC==0 && NuPDG==12))
		{
		  
		  h_Stat->SetBinContent(13,h_Stat->GetBinContent(13)+1);
		  h_trueNuE_Energy_FidVol_13->Fill(enutruth*1000.);
		  h_trueLep_Mom_FidVol_13->Fill(lep_mom_truth*1000.);
		  if(DropTwoClusterEvent==false)
		    {h_Stat->SetBinContent(15,h_Stat->GetBinContent(15)+1);
		      h_trueNuE_Energy_FidVol_15->Fill(enutruth*1000.);
		      h_trueLep_Mom_FidVol_15->Fill(lep_mom_truth*1000.);
		    }
		}
	    }//if MC sample == true
	  }
      }
    
    //in the event loop test bnb + cosmic sample, only consider the flash cluster matching
    if(FlashInBeamWindow)
      {
	if(ClusterFlashMatch0==true || ClusterFlashMatch1==true || ClusterFlashMatch2==true)//at least one plane cluster pass the cut
	  {
	    h_BNBCosmicStat->SetBinContent(25,h_BNBCosmicStat->GetBinContent(25)+1);
	    if(MCSample){
	      if((ClusterFlashMatchNueEvent0 ==true || ClusterFlashMatchNueEvent1 == true || ClusterFlashMatchNueEvent2==true) && inFV(nuvtxx,nuvtxy,nuvtxz))
		{
		  h_BNBCosmicStat->SetBinContent(26,h_BNBCosmicStat->GetBinContent(26)+1);
		}
	      if(ClusterFlashMatchNuOthersEvent0 ==true || ClusterFlashMatchNuOthersEvent1==true || ClusterFlashMatchNuOthersEvent2==true)
		{	h_BNBCosmicStat->SetBinContent(27,h_BNBCosmicStat->GetBinContent(27)+1);}
	      if(ClusterFlashMatchCosmicEvent0==true || ClusterFlashMatchCosmicEvent1==true || ClusterFlashMatchCosmicEvent2==true)
		{
		  h_BNBCosmicStat->SetBinContent(28,h_BNBCosmicStat->GetBinContent(28)+1);
		}
	    }
	  }
	if((ClusterFlashMatch0==true && ClusterFlashMatch1==true) || (ClusterFlashMatch0==true && ClusterFlashMatch2==true) || (ClusterFlashMatch1==true && ClusterFlashMatch2==true))//at least two plane cluster pass the cut
	  {
	    h_BNBCosmicStat->SetBinContent(29,h_BNBCosmicStat->GetBinContent(29)+1);
	    h_FlashMatchedClusterSize->Fill(FlashMatchedClusterVec.size());
	    if(DropTwoClusterEvent==false)
	      {
		h_BNBCosmicStat->SetBinContent(33,h_BNBCosmicStat->GetBinContent(33)+1);	
	      }
	    if(MCSample){
	      if((ClusterFlashMatchNueEvent0 ==true || ClusterFlashMatchNueEvent1 == true || ClusterFlashMatchNueEvent2==true) && inFV(nuvtxx,nuvtxy,nuvtxz))
		{
		  h_BNBCosmicStat->SetBinContent(30,h_BNBCosmicStat->GetBinContent(30)+1);
		  h_FlashMatchedClusterSize_Nue->Fill(FlashMatchedClusterVec.size());
		  if(DropTwoClusterEvent==false)
		    {h_BNBCosmicStat->SetBinContent(34,h_BNBCosmicStat->GetBinContent(34)+1);}
		}
	      if(ClusterFlashMatchNuOthersEvent0 ==true || ClusterFlashMatchNuOthersEvent1==true || ClusterFlashMatchNuOthersEvent2==true)
		{	
		  h_BNBCosmicStat->SetBinContent(31,h_BNBCosmicStat->GetBinContent(31)+1);
		  h_FlashMatchedClusterSize_NuOthers->Fill(FlashMatchedClusterVec.size());
		  if(DropTwoClusterEvent==false)
		    {
		      h_BNBCosmicStat->SetBinContent(35,h_BNBCosmicStat->GetBinContent(35)+1);	
		    }
		}
	      if(ClusterFlashMatchCosmicEvent0==true || ClusterFlashMatchCosmicEvent1==true || ClusterFlashMatchCosmicEvent2==true)
		{
		  h_BNBCosmicStat->SetBinContent(32,h_BNBCosmicStat->GetBinContent(32)+1);
		  h_FlashMatchedClusterSize_Cosmic->Fill(FlashMatchedClusterVec.size());
		  if(DropTwoClusterEvent==false)
		    {
		      h_BNBCosmicStat->SetBinContent(36,h_BNBCosmicStat->GetBinContent(36)+1);	
		    }
		}
	    }//if MC sample
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
	      h_Stat->SetBinContent(16,h_Stat->GetBinContent(16)+1);
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
  
  h_NflashInBeamSpill = tfs->make<TH1F>("h_NflashInBeamSpill","number of flash in beam spill one for each plane",20,0,20);
  h_flashWidth0 = tfs->make<TH1F>("h_flashWidth0","flash width projected to plane 0",200,0,4000);
  h_flashWidth1 = tfs->make<TH1F>("h_flashWidth1","flash width projected to plane 0",200,0,4000);
  h_flashWidth2 = tfs->make<TH1F>("h_flashWidth2","flash width projected to plane 0",200,0,4000);
  h_flashZWidth = tfs->make<TH1F>("h_flashZWidth","flash width in Z direction",100,0,1000);
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
  h_BNBCosmicStat = tfs->make<TH1F>("h_BNBCosmicStat","statistics of each cut",40,0.,40.);
  h_DataStat = tfs->make<TH1F>("h_DataStat","statistics of each cut",20,0.,20.);
  h_ophits_per_flash = tfs->make<TH1F>("h_ophits_per_flash","OpHits per Flash;N_{optical hits};Events / bin",20,-0.5,19.5);
  h_ophits_per_flash_2pe = tfs->make<TH1F>("h_ophits_per_flash_2pe","OpHits (> 2 PE) per Flash;N_{optical hits};Events / bin",20,-0.5,19.5);

  h_trueNuE_Energy_FidVol = tfs->make<TH1F>("h_trueNuE_Energy_FidVol","true nue energy for all the ccnue vertex inside FidVol",100,0,5000);
  h_trueNuE_Energy_FidVol_10 = tfs->make<TH1F>("h_trueNuE_Energy_FidVol_10","true nue energy for at least one shower cluster",100,0,5000);
  h_trueNuE_Energy_FidVol_11 = tfs->make<TH1F>("h_trueNuE_Energy_FidVol_11","true nue energy for at least one shower cluster from nue e",100,0,5000);
  h_trueNuE_Energy_FidVol_12 = tfs->make<TH1F>("h_trueNuE_Energy_FidVol_12","true nue energy for at least two plane shower clusters",100,0,5000);
  h_trueNuE_Energy_FidVol_13 = tfs->make<TH1F>("h_trueNuE_Energy_FidVol_13","true nue energy for at least two plane shower clusters",100,0,5000);
  h_trueNuE_Energy_FidVol_14 = tfs->make<TH1F>("h_trueNuE_Energy_FidVol_14","true nue energy for at least two plane shower clusters",100,0,5000);
  h_trueNuE_Energy_FidVol_15 = tfs->make<TH1F>("h_trueNuE_Energy_FidVol_15","true nue energy for at least two plane shower clusters",100,0,5000);
  h_trueLep_Mom_FidVol = tfs->make<TH1F>("h_trueLep_Mom_FidVol","true nue energy for all the ccnue vertex inside FidVol",100,0,5000);
  h_trueLep_Mom_FidVol_10 = tfs->make<TH1F>("h_trueLep_Mom_FidVol_10","true e- momentum for at least one shower cluster",100,0,5000);
  h_trueLep_Mom_FidVol_11 = tfs->make<TH1F>("h_trueLep_Mom_FidVol_11","true e- momentum for at least one shower cluster from nue e",100,0,5000);
  h_trueLep_Mom_FidVol_12 = tfs->make<TH1F>("h_trueLep_Mom_FidVol_12","true e- momentum for at least two plane shower clusters",100,0,5000);
  h_trueLep_Mom_FidVol_13 = tfs->make<TH1F>("h_trueLep_Mom_FidVol_13","true e- momentum for at least two plane shower clusters",100,0,5000);
  h_trueLep_Mom_FidVol_14 = tfs->make<TH1F>("h_trueLep_Mom_FidVol_14","true e- momentum for at least two plane shower clusters",100,0,5000);
 h_trueLep_Mom_FidVol_15 = tfs->make<TH1F>("h_trueLep_Mom_FidVol_15","true e- momentum for at least two plane shower clusters",100,0,5000);
  h_pdgcodeNuOthers0 = tfs->make<TH1F>("h_pdgcodeNuOthers0","pdg code origin 1 but not electron plane 0",5000,-2500,2500);
  h_pdgcodeNuOthers1 = tfs->make<TH1F>("h_pdgcodeNuOthers1","pdg code origin 1 but not electron plane 1",5000,-2500,2500);
  h_pdgcodeNuOthers2 = tfs->make<TH1F>("h_pdgcodeNuOthers2","pdg code origin 1 but not electron plane 2",5000,-2500,2500);
  h_pdgcodeCosmic0 = tfs->make<TH1F>("h_pdgcodeCosmic0","pdg code origin 2 plane 0",5000,-2500,2500);
  h_pdgcodeCosmic1 = tfs->make<TH1F>("h_pdgcodeCosmic1","pdg code origin 2 plane 1",5000,-2500,2500);
  h_pdgcodeCosmic2 = tfs->make<TH1F>("h_pdgcodeCosmic2","pdg code origin 2 plane 2",5000,-2500,2500);
  h_StartAngleNue2 = tfs->make<TH1F>("h_StartAngleNue2","flash matched plane 2cluster start angle from nue",36,-3.2,3.2);
  h_StartOpenAngleNue2 = tfs->make<TH1F>("h_StartOpenAngleNue2","flash matched plane 2cluster start opening angle from nue",36,-3.2,3.2);
  h_EndAngleNue2 = tfs->make<TH1F>("h_EndAngleNue2","flash matched plane 2cluster end angle from nue",36,-3.2,3.2);
  h_EndOpenAngleNue2 = tfs->make<TH1F>("h_EndOpenAngleNue2","flash matched plane 2cluster end opening angle from nue",36,-3.2,3.2);
  h_StartAngleNuOthers2 = tfs->make<TH1F>("h_StartAngleNuOthers2","flash matched plane 2cluster start angle from nu products but not e-",36,-3.2,3.2);
  h_StartOpenAngleNuOthers2 = tfs->make<TH1F>("h_StartOpenAngleNuOthers2","flash matched plane 2cluster start opening angle from nu products but not e-",36,-3.2,3.2);
  h_EndAngleNuOthers2 = tfs->make<TH1F>("h_EndAngleNuOthers2","flash matched plane 2cluster end angle from nu products but not e-",36,-3.2,3.2);
  h_EndOpenAngleNuOthers2 = tfs->make<TH1F>("h_EndOpenAngleNuOthers2","flash matched plane 2cluster end opening angle from nu products but not e-",36,-3.2,3.2);
   h_StartAngleCosmic2 = tfs->make<TH1F>("h_StartAngleCosmic2","flash matched plane 2cluster start angle from cosmic",36,-3.2,3.2);
  h_StartOpenAngleCosmic2 = tfs->make<TH1F>("h_StartOpenAngleCosmic2","flash matched plane 2cluster start opening angle from cosmic",36,-3.2,3.2);
  h_EndAngleCosmic2 = tfs->make<TH1F>("h_EndAngleCosmic2","flash matched plane 2cluster end angle from cosmic",36,-3.2,3.2);
  h_EndOpenAngleCosmic2 = tfs->make<TH1F>("h_EndOpenAngleCosmic2","flash matched plane 2cluster end opening angle from cosmic",36,-3.2,3.2);

  h_clusterIntegralEnergy_Nue_0 = tfs->make<TH2F>("h_clusterIntegralEnergy_Nue_0","cluster charge inteigral plane 0 nue",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_Nue_0 = tfs->make<TH2F>("h_clusterSummedADCEnergy_Nue_0","cluster charge SummedADC plane 0 nue",100,0,5000,100,0,5000);
  h_clusterIntegralEnergy_NuOthers_0 = tfs->make<TH2F>("h_clusterIntegralEnergy_NuOthers_0","cluster charge inteigral plane 0 nu products not e-",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_NuOthers_0 = tfs->make<TH2F>("h_clusterSummedADCEnergy_NuOthers_0","cluster charge SummedADC plane 0 nu products not e-",100,0,5000,100,0,5000);
  h_clusterIntegralEnergy_Cosmic_0 = tfs->make<TH2F>("h_clusterIntegralEnergy_Cosmic_0","cluster charge inteigral plane 0 cosmic",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_Cosmic_0 = tfs->make<TH2F>("h_clusterSummedADCEnergy_Cosmic_0","cluster charge SummedADC plane 0 cosmic",100,0,5000,100,0,5000);

  h_clusterIntegralEnergy_Nue_1 = tfs->make<TH2F>("h_clusterIntegralEnergy_Nue_1","cluster charge inteigral plane 1 nue",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_Nue_1 = tfs->make<TH2F>("h_clusterSummedADCEnergy_Nue_1","cluster charge SummedADC plane 1 nue",100,0,5000,100,0,5000);
  h_clusterIntegralEnergy_NuOthers_1 = tfs->make<TH2F>("h_clusterIntegralEnergy_NuOthers_1","cluster charge inteigral plane 1 nu products not e-",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_NuOthers_1 = tfs->make<TH2F>("h_clusterSummedADCEnergy_NuOthers_1","cluster charge SummedADC plane 1 nu products not e-",100,0,5000,100,0,5000);
  h_clusterIntegralEnergy_Cosmic_1 = tfs->make<TH2F>("h_clusterIntegralEnergy_Cosmic_1","cluster charge inteigral plane 1 cosmic",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_Cosmic_1 = tfs->make<TH2F>("h_clusterSummedADCEnergy_Cosmic_1","cluster charge SummedADC plane 1 cosmic",100,0,5000,100,0,5000);

  h_clusterIntegralEnergy_Nue_2 = tfs->make<TH2F>("h_clusterIntegralEnergy_Nue_2","cluster charge inteigral plane 2 nue",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_Nue_2 = tfs->make<TH2F>("h_clusterSummedADCEnergy_Nue_2","cluster charge SummedADC plane 2 nue",100,0,5000,100,0,5000);
  h_clusterIntegralEnergy_NuOthers_2 = tfs->make<TH2F>("h_clusterIntegralEnergy_NuOthers_2","cluster charge inteigral plane 2 nu products not e-",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_NuOthers_2 = tfs->make<TH2F>("h_clusterSummedADCEnergy_NuOthers_2","cluster charge SummedADC plane 2 nu products not e-",100,0,5000,100,0,5000);
  h_clusterIntegralEnergy_Cosmic_2 = tfs->make<TH2F>("h_clusterIntegralEnergy_Cosmic_2","cluster charge inteigral plane 2 cosmic",100,0,5000,100,0,5000);
  h_clusterSummedADCEnergy_Cosmic_2 = tfs->make<TH2F>("h_clusterSummedADCEnergy_Cosmic_2","cluster charge SummedADC plane 2 cosmic",100,0,5000,100,0,5000);

  h_MatchedClusterFromDifferentPlaneTimeDiff_nue = tfs->make<TH1F>("h_MatchedClusterFromDifferentPlaneTimeDiff_nue","time difference between two clusters from different plane but are both electron from nue",500,0,5000);
  h_MatchedClusterFromDifferentPlaneTimeDiff_NuOthers = tfs->make<TH1F>("h_MatchedClusterFromDifferentPlaneTimeDiff_NuOthers","time difference between two clusters from different plane but are both from nu interaction with at least one is not electron",500,0,5000);
  h_MatchedClusterFromDifferentPlaneTimeDiff_Cosmic = tfs->make<TH1F>("h_MatchedClusterFromDifferentPlaneTimeDiff_Cosmic","time difference between two clusters from different plane and at least one is from cosmic",500,0,5000);
  h_FlashMatchedClusterSize = tfs->make<TH1F>("h_FlashMatchedClusterSize","for the case of at least two plane cluster pass, the size of the flashMatched Cluster Vector",10,0,10);
  h_FlashMatchedClusterSize_Nue = tfs->make<TH1F>("h_FlashMatchedClusterSize_Nue","for the case of at least two plane cluster pass, the size of the flashMatched Cluster Vector with at least one cluster is from nue",10,0,10);
  h_FlashMatchedClusterSize_NuOthers = tfs->make<TH1F>("h_FlashMatchedClusterSize_NuOthers","for the case of at least two plane cluster pass, the size of the flashMatched Cluster Vector with at least one cluster is from other nu products",10,0,10);
  h_FlashMatchedClusterSize_Cosmic = tfs->make<TH1F>("h_FlashMatchedClusterSize_Cosmic","for the case of at least two plane cluster pass, the size of the flashMatched Cluster Vector with at least one cluster is from cosmic",10,0,10);
  
  h_FlashCluster_Angle_nue = tfs->make<TH2F>("h_FlashCluster_Angle_nue","flash matched nue e- cluster angle",64,-3.2,3.2,4,-0.5,3.5);
  h_FlashCluster_NHits_nue = tfs->make<TH2F>("h_FlashCluster_NHits_nue","flash matched nue e- cluster nhits",1000,0,1000,4,-0.5,3.5);
  h_FlashCluster_Width_nue = tfs->make<TH2F>("h_FlashCluster_Width_nue","flash matched nue e- cluster width",1000,0,1000,4,-0.5,3.5);
  h_FlashCluster_StartTime_nue = tfs->make<TH2F>("h_FlashCluster_StartTime_nue","flash matched nue e- cluster width",1000,0,10000,4,-0.5,3.5);

  h_FlashCluster_Angle_NuOthers = tfs->make<TH2F>("h_FlashCluster_Angle_NuOthers","flash matched nue e- cluster angle",64,-3.2,3.2,4,-0.5,3.5);
  h_FlashCluster_NHits_NuOthers = tfs->make<TH2F>("h_FlashCluster_NHits_NuOthers","flash matched nue e- cluster nhits",1000,0,1000,4,-0.5,3.5);
  h_FlashCluster_Width_NuOthers = tfs->make<TH2F>("h_FlashCluster_Width_NuOthers","flash matched nue e- cluster width",1000,0,1000,4,-0.5,3.5);
  h_FlashCluster_StartTime_NuOthers = tfs->make<TH2F>("h_FlashCluster_StartTime_NuOthers","flash matched nue e- cluster width",1000,0,10000,4,-0.5,3.5);

  h_FlashCluster_Angle_Cosmic = tfs->make<TH2F>("h_FlashCluster_Angle_Cosmic","flash matched nue e- cluster angle",64,-3.2,3.2,4,-0.5,3.5);
  h_FlashCluster_NHits_Cosmic = tfs->make<TH2F>("h_FlashCluster_NHits_Cosmic","flash matched nue e- cluster nhits",1000,0,1000,4,-0.5,3.5);
  h_FlashCluster_Width_Cosmic = tfs->make<TH2F>("h_FlashCluster_Width_Cosmic","flash matched nue e- cluster width",1000,0,1000,4,-0.5,3.5);
  h_FlashCluster_StartTime_Cosmic = tfs->make<TH2F>("h_FlashCluster_StartTime_Cosmic","flash matched nue e- cluster width",1000,0,10000,4,-0.5,3.5);

  h_AngleNHits_nue_0 = tfs->make<TH2F>("h_AngleNHits_nue_0","before flash matching angle (X) Vs nHits(Y) from nue cluster on plane 0",128,-3.2,3.2,200,0,400);
  h_AngleNHits_nue_1 = tfs->make<TH2F>("h_AngleNHits_nue_1","before flash matching angle (X) Vs nHits(Y) from nue cluster on plane 1",128,-3.2,3.2,200,0,400);
  h_AngleNHits_nue_2 = tfs->make<TH2F>("h_AngleNHits_nue_2","before flash matching angle (X) Vs nHits(Y) from nue cluster on plane 2",128,-3.2,3.2,200,0,400);
  h_AngleNHits_NuOthers_0 = tfs->make<TH2F>("h_AngleNHits_NuOthers_0","before flash matching angle (X) Vs nHits(Y) from other non-e- nu cluster on plane 0",128,-3.2,3.2,200,0,400);
  h_AngleNHits_NuOthers_1 = tfs->make<TH2F>("h_AngleNHits_NuOthers_1","before flash matching angle (X) Vs nHits(Y) from other non-e- nu cluster on plane 1",128,-3.2,3.2,200,0,400);
  h_AngleNHits_NuOthers_2 = tfs->make<TH2F>("h_AngleNHits_NuOthers_2","before flash matching angle (X) Vs nHits(Y) from other non-e- nu cluster on plane 2",128,-3.2,3.2,200,0,400);
   h_AngleNHits_Cosmic_0 = tfs->make<TH2F>("h_AngleNHits_Cosmic_0","before flash matching angle (X) Vs nHits(Y) from cosmic cluster on plane 0",128,-3.2,3.2,200,0,400);
  h_AngleNHits_Cosmic_1 = tfs->make<TH2F>("h_AngleNHits_Cosmic_1","before flash matching angle (X) Vs nHits(Y) from cosmic cluster on plane 1",128,-3.2,3.2,200,0,400);
  h_AngleNHits_Cosmic_2 = tfs->make<TH2F>("h_AngleNHits_Cosmic_2","before flash matching angle (X) Vs nHits(Y) from cosmic cluster on plane 2",128,-3.2,3.2,200,0,400);


  //data
  h_StartWire_diff_data_0 = tfs->make<TH1F>("h_StartWire_diff_data_0","start wire difference between flash and cluster on plane 0 for data",100,-2000,2000);
  h_EndWire_diff_data_0 = tfs->make<TH1F>("h_EndWire_diff_data_0","end wire difference between flash and cluster on plane 0 for data",100,-2000,2000);
  h_AngleNHits_data_0 = tfs->make<TH2F>("h_AngleNHits_data_0","angle Vs nhits for all clusters on plane 0 for data",128,-3.2,3.2,200,0,400);
  h_StartWire_diff_data_1 = tfs->make<TH1F>("h_StartWire_diff_data_1","start wire difference between flash and cluster on plane 1 for data",100,-2000,2000);
  h_EndWire_diff_data_1 = tfs->make<TH1F>("h_EndWire_diff_data_1","end wire difference between flash and cluster on plane 1 for data",100,-2000,2000);
  h_AngleNHits_data_1 = tfs->make<TH2F>("h_AngleNHits_data_1","angle Vs nhits for all clusters on plane 1 for data",128,-3.2,3.2,200,0,400); 
  h_StartWire_diff_data_2 = tfs->make<TH1F>("h_StartWire_diff_data_2","start wire difference between flash and cluster on plane 2 for data",100,-2000,2000);
  h_EndWire_diff_data_2 = tfs->make<TH1F>("h_EndWire_diff_data_2","end wire difference between flash and cluster on plane 2 for data",100,-2000,2000);
  h_AngleNHits_data_2 = tfs->make<TH2F>("h_AngleNHits_data_2","angle Vs nhits for all clusters on plane 2 for data",128,-3.2,3.2,200,0,400);
  h_CenterWire_diff_data_0 = tfs->make<TH1F>("h_CenterWire_diff_data_0","center wire difference between flash and cluster on plane 0 for data",100,0,4000);
  h_CenterWire_diff_data_1 = tfs->make<TH1F>("h_CenterWire_diff_data_1","center wire difference between flash and cluster on plane 1 for data",100,0,4000);
  h_CenterWire_diff_data_2 = tfs->make<TH1F>("h_CenterWire_diff_data_2","center wire difference between flash and cluster on plane 2 for data",100,0,4000);
  h_StartAngle_data_goodcluster_0 = tfs->make<TH1F>("h_StartAngle_data_goodcluster_0","flash matched cluster start angle for data on plane 0",64,-3.2,3.2);
  h_NHits_data_goodcluster_0 = tfs->make<TH1F>("h_NHits_data_goodcluster_0","flash matched cluster nhits for data on plane 0",1000,0,1000);
  h_TotalEnergy_data_goodcluster_0 = tfs->make<TH1F>("h_TotalEnergy_data_goodcluster_0","flash matched cluster total energy for dataon plane 0 ",500,0,5000);
  h_StartAngle_data_goodcluster_1 = tfs->make<TH1F>("h_StartAngle_data_goodcluster_1","flash matched cluster start angle for data on plane 1",64,-3.2,3.2);
  h_NHits_data_goodcluster_1 = tfs->make<TH1F>("h_NHits_data_goodcluster_1","flash matched cluster nhits for data on plane 1",1000,0,1000);
  h_TotalEnergy_data_goodcluster_1 = tfs->make<TH1F>("h_TotalEnergy_data_goodcluster_1","flash matched cluster total energy for dataon plane 1",500,0,5000);
  h_StartAngle_data_goodcluster_2 = tfs->make<TH1F>("h_StartAngle_data_goodcluster_2","flash matched cluster start angle for data on plane 2",64,-3.2,3.2);
  h_NHits_data_goodcluster_2 = tfs->make<TH1F>("h_NHits_data_goodcluster_2","flash matched cluster nhits for data on plane 2",1000,0,1000);
  h_TotalEnergy_data_goodcluster_2 = tfs->make<TH1F>("h_TotalEnergy_data_goodcluster_2","flash matched cluster total energy for dataon plane 2",500,0,5000);
  
  h_TwoCluster_TimeDiff = tfs->make<TH1F>("h_TwoCluster_TimeDiff","two clusters from two planes time difference for all case ",500,0,5000);
  h_TwoCluster_TimeDiff_nue = tfs->make<TH1F>("h_TwoCluster_TimeDiff_nue","two clusters from two planes time difference for both clusters are e- from nue ",500,0,5000);
  h_TwoCluster_TimeDiff_Cosmic= tfs->make<TH1F>("h_TwoCluster_TimeDiff_Cosmic","two clusters from two planes time difference for at least one cluster is from cosmic ",500,0,5000);
  h_NumberOfGoodCluster_data = tfs->make<TH1F>("h_NumberOfGoodCluster_data","number of clusters for plane 0,1,2 (X xbin) for all clusters",4,-0.5,3.5);
  h_NumberOfGoodCluster_nue = tfs->make<TH1F>("h_NumberOfGoodCluster_nue","number of clusters for plane 0,1,2 (X xbin) for e- clusters",4,-0.5,3.5);
  h_NumberOfGoodCluster_NuOthers= tfs->make<TH1F>("h_NumberOfGoodCluster_NuOthers","number of clusters for plane 0,1,2 (X xbin) for other nu clusters",4,-0.5,3.5);
  h_NumberOfGoodCluster_Cosmic = tfs->make<TH1F>("h_NumberOfGoodCluster_Cosmic","number of clusters for plane 0,1,2 (X xbin) for cosmic clusters",4,-0.5,3.5);
  h_ClusterCenterVsWidth_UplaneSpecialCut = tfs->make<TH2F>("h_ClusterCenterVsWidth_UplaneSpecialCut","V plane good cluster center",300,0,3000,100,0,1000);
  h_ClusterCenterVsWidth_VplaneSpecialCut = tfs->make<TH2F>("h_ClusterCenterVsWidth_VplaneSpecialCut","V plane good cluster center",300,0,3000,100,0,1000);
  h_ClusterCenterVsWidth_YplaneSpecialCut = tfs->make<TH2F>("h_ClusterCenterVsWidth_YplaneSpecialCut","V plane good cluster center",300,0,3000,100,0,1000);
}

void ub::LowLevelNueFilter::GetTruthInfo(std::vector<art::Ptr<recob::Hit>> hits, int& origin, int& pdgcode){
  
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  origin = 0;
  pdgcode = 0;

  std::map<int,double> trkide;
  for(size_t h = 0; h < hits.size(); ++h){
    art::Ptr<recob::Hit> hit = hits[h];
    //std::vector<sim::IDE> ides;
    std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
    
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
  
  const simb::MCParticle* particle=pi_serv->TrackIdToParticle_P(Trackid);	    
  if(particle){
    pdgcode = particle->PdgCode();
    origin = pi_serv->TrackIdToMCTruth_P(Trackid)->Origin();
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
void ub::LowLevelNueFilter::doEfficiencies(){
  art::ServiceHandle<art::TFileService> tfs;
  if(TEfficiency::CheckConsistency(*h_trueNuE_Energy_FidVol_10,*h_trueNuE_Energy_FidVol)){
    h_Eff_Ev_10 = tfs->make<TEfficiency>(*h_trueNuE_Energy_FidVol_10,*h_trueNuE_Energy_FidVol);
    TGraphAsymmErrors *grEff_Ev_10 = h_Eff_Ev_10->CreateGraph();
    grEff_Ev_10->Write("grEff_Ev_10");
    h_Eff_Ev_10->Write("h_Eff_Ev_10");
  }
  if(TEfficiency::CheckConsistency(*h_trueNuE_Energy_FidVol_11,*h_trueNuE_Energy_FidVol)){
    h_Eff_Ev_11 = tfs->make<TEfficiency>(*h_trueNuE_Energy_FidVol_11,*h_trueNuE_Energy_FidVol);
    TGraphAsymmErrors *grEff_Ev_11 = h_Eff_Ev_11->CreateGraph();
    grEff_Ev_11->Write("grEff_Ev_11");
    h_Eff_Ev_11->Write("h_Eff_Ev_11");
  }
  if(TEfficiency::CheckConsistency(*h_trueNuE_Energy_FidVol_12,*h_trueNuE_Energy_FidVol)){
    h_Eff_Ev_12 = tfs->make<TEfficiency>(*h_trueNuE_Energy_FidVol_12,*h_trueNuE_Energy_FidVol);
    TGraphAsymmErrors *grEff_Ev_12 = h_Eff_Ev_12->CreateGraph();
    grEff_Ev_12->Write("grEff_Ev_12");
    h_Eff_Ev_12->Write("h_Eff_Ev_12");
  }
  if(TEfficiency::CheckConsistency(*h_trueNuE_Energy_FidVol_13,*h_trueNuE_Energy_FidVol)){
    h_Eff_Ev_13 = tfs->make<TEfficiency>(*h_trueNuE_Energy_FidVol_13,*h_trueNuE_Energy_FidVol);
    TGraphAsymmErrors *grEff_Ev_13 = h_Eff_Ev_13->CreateGraph();
    grEff_Ev_13->Write("grEff_Ev_13");
    h_Eff_Ev_13->Write("h_Eff_Ev_13");
  }

  //lepton momentum
if(TEfficiency::CheckConsistency(*h_trueLep_Mom_FidVol_10,*h_trueLep_Mom_FidVol)){
    h_Eff_lm_10 = tfs->make<TEfficiency>(*h_trueLep_Mom_FidVol_10,*h_trueLep_Mom_FidVol);
    TGraphAsymmErrors *grEff_lm_10 = h_Eff_lm_10->CreateGraph();
    grEff_lm_10->Write("grEff_lm_10");
    h_Eff_lm_10->Write("h_Eff_lm_10");
  }
  if(TEfficiency::CheckConsistency(*h_trueLep_Mom_FidVol_11,*h_trueLep_Mom_FidVol)){
    h_Eff_lm_11 = tfs->make<TEfficiency>(*h_trueLep_Mom_FidVol_11,*h_trueLep_Mom_FidVol);
    TGraphAsymmErrors *grEff_lm_11 = h_Eff_lm_11->CreateGraph();
    grEff_lm_11->Write("grEff_lm_11");
    h_Eff_lm_11->Write("h_Eff_lm_11");
  }
  if(TEfficiency::CheckConsistency(*h_trueLep_Mom_FidVol_12,*h_trueLep_Mom_FidVol)){
    h_Eff_lm_12 = tfs->make<TEfficiency>(*h_trueLep_Mom_FidVol_12,*h_trueLep_Mom_FidVol);
    TGraphAsymmErrors *grEff_lm_12 = h_Eff_lm_12->CreateGraph();
    grEff_lm_12->Write("grEff_lm_12");
    h_Eff_lm_12->Write("h_Eff_lm_12");
  }
  if(TEfficiency::CheckConsistency(*h_trueLep_Mom_FidVol_13,*h_trueLep_Mom_FidVol)){
    h_Eff_lm_13 = tfs->make<TEfficiency>(*h_trueLep_Mom_FidVol_13,*h_trueLep_Mom_FidVol);
    TGraphAsymmErrors *grEff_lm_13 = h_Eff_lm_13->CreateGraph();
    grEff_lm_13->Write("grEff_lm_13");
    h_Eff_lm_13->Write("h_Eff_lm_13");
  }
  
}

DEFINE_ART_MODULE(ub::LowLevelNueFilter)

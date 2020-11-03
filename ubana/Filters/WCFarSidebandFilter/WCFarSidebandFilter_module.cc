////////////////////////////////////////////////////////////////////////
// Class:       WCFarSidebandFilter
// Plugin Type: filter (art v3_02_06)
// File:        WCFarSidebandFilter_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"
#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"
#include "ubobj/WcpPort/NuSelectionKINE.h"

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TFile.h"


class WCFarSidebandFilter;


class WCFarSidebandFilter : public art::EDFilter {
public:
  explicit WCFarSidebandFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCFarSidebandFilter(WCFarSidebandFilter const&) = delete;
  WCFarSidebandFilter(WCFarSidebandFilter&&) = delete;
  WCFarSidebandFilter& operator=(WCFarSidebandFilter const&) = delete;
  WCFarSidebandFilter& operator=(WCFarSidebandFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:
  std::string fPFInputTag;
  bool fBDTvars;
  bool fKINEvars;
  float wcp_numu_score, wcp_nue_score;
  float wcp_reco_enu;

  art::Handle< std::vector<nsm::NuSelectionBDT> > bdthandle;
  std::vector<art::Ptr<nsm::NuSelectionBDT> > bdtvec;
  art::Ptr<nsm::NuSelectionBDT> bdt;
  art::Handle< std::vector<nsm::NuSelectionKINE> > kinehandle;
  std::vector<art::Ptr<nsm::NuSelectionKINE> > kinevec;
  art::Ptr<nsm::NuSelectionKINE> kine;
};


WCFarSidebandFilter::WCFarSidebandFilter(fhicl::ParameterSet const& p) : EDFilter{p} {
  fBDTvars = p.get<bool>("BDTvars");
  fKINEvars = p.get<bool>("KINEvars");
  fPFInputTag = p.get<std::string>("PF_inputtag");
}


bool WCFarSidebandFilter::filter(art::Event& e) {

  bdthandle.clear();
  bdtvec.clear();

  kinehandle.clear();
  kinevec.clear();

  if(fBDTvars){
    //art::Handle< std::vector<nsm::NuSelectionBDT> > bdthandle;
    if (! e.getByLabel(fPFInputTag, bdthandle)) return false;
    //std::vector<art::Ptr<nsm::NuSelectionBDT> > bdtvec;
    art::fill_ptr_vector(bdtvec,bdthandle);
    if(bdtvec.size()!=1) {
      //std::cout<<"WARNING: no set of BDT input variables" << std::endl;
      return false;
    }
    //const nsm::NuSelectionBDT::BDTscores BDTscores = bdtvec[0]->GetBDTscores();
    wcp_numu_score = bdtvec[0]->GetBDTscores().numu_score;
    wcp_nue_score = bdtvec[0]->GetBDTscores().nue_score;
    
    //std::cout<<"BDT vars:  " << wcp_numu_score << ",  " << wcp_nue_score << "   -----------------------------------------------------------" << std::endl;
    //if (wcp_numu_score >= 0.7 && wcp_nue_score <= -1.0 ) { return true; }
    //std::cout<<"nue score: "<<wcp_nue_score<<std::endl;
    // far sideband cut
  }

  if(fKINEvars){
    if (! e.getByLabel(fPFInputTag, kinehandle)) return false;
    art::fill_ptr_vector(kinevec,kinehandle);
    if(kinevec.size()!=1) {
      //std::cout<<"WARNING: >1 set of KINE input variables" << std::endl;
      return false;
    }

    wcp_reco_enu = kinevec[0]->GetKineInfo().kine_reco_Enu;
    //std::cout<<"reco energy: "<<wcp_reco_enu<<std::endl;
  }

  // far sideband cut
  
  if( 
     (wcp_nue_score>=-16 && wcp_nue_score<0 && wcp_reco_enu>=0 && wcp_reco_enu<=2500)
     ||
     (wcp_reco_enu >=800 && wcp_reco_enu <=2500 && wcp_nue_score>=-16 && wcp_nue_score<=16)
      ) {
    std::cout<<"Pass!"<<" -----------------------------------------------------"<<std::endl;
    std::cout<<"nue score: "<<wcp_nue_score<<std::endl;
    std::cout<<"reco energy: "<<wcp_reco_enu<<std::endl;
    return true;
  }

  return false;
}


DEFINE_ART_MODULE(WCFarSidebandFilter)

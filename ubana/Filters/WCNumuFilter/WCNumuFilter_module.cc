////////////////////////////////////////////////////////////////////////
// Class:       WCNumuFilter
// Plugin Type: filter (art v3_02_06)
// File:        WCNumuFilter_module.cc
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

//#include "PID/LLR_PID.h"
//#include "PID/LLRPID_proton_muon_lookup.h"

#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"
#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"

//#include "ubana/WcpPortedAna/Helpers/TruthContainersHelper.h"
//#include "ubana/WcpPortedAna/Helpers/TrackMatchHelper.h"
//#include "ubana/WcpPortedAna/Helpers/ShowerMatchHelper.h"
//#include "ubana/WcpPortedAna/Helpers/ClusterMatchHelper.h"
//#include "ubana/WcpPortedAna/Helpers/VertexMatchHelper.h"

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TFile.h"


class WCNumuFilter;


class WCNumuFilter : public art::EDFilter {
	public:
		explicit WCNumuFilter(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		WCNumuFilter(WCNumuFilter const&) = delete;
		WCNumuFilter(WCNumuFilter&&) = delete;
		WCNumuFilter& operator=(WCNumuFilter const&) = delete;
		WCNumuFilter& operator=(WCNumuFilter&&) = delete;

		// Required functions.
		bool filter(art::Event& e) override;

	private:
		std::string fPFInputTag;
		bool fBDTvars;
		float wcp_numu_score, wcp_nue_score;

		art::Handle< std::vector<nsm::NuSelectionBDT> > bdthandle;
		std::vector<art::Ptr<nsm::NuSelectionBDT> > bdtvec;
		art::Ptr<nsm::NuSelectionBDT> bdt;
};


WCNumuFilter::WCNumuFilter(fhicl::ParameterSet const& p) : EDFilter{p} {
	fBDTvars = p.get<bool>("BDTvars");
	fPFInputTag = p.get<std::string>("PF_inputtag");
}


bool WCNumuFilter::filter(art::Event& e) {

	bdthandle.clear();
	bdtvec.clear();
	if(fBDTvars){
		if (! e.getByLabel(fPFInputTag, bdthandle)) return false;
		art::fill_ptr_vector(bdtvec,bdthandle);
		if(bdtvec.size()!=1) {
			//std::cout<<"WARNING: no set of BDT input variables" << std::endl;
			return false;
		}
		wcp_numu_score = bdtvec[0]->GetBDTscores().numu_score;
		wcp_nue_score = bdtvec[0]->GetBDTscores().nue_score;
		if (wcp_numu_score >= 0.7 && wcp_nue_score <= -1.0 ) { return true; }
		//return true;
	}
	//std::cout<<"WARNING: no BDT vars" << std::endl;
	return false;
}


DEFINE_ART_MODULE(WCNumuFilter)

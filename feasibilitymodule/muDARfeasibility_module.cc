////////////////////////////////////////////////////////////////////////
// Class:       muDARfeasibility
// Plugin Type: analyzer (art v3_01_02)
// File:        muDARfeasibility_module.cc
//
// Generated at Thu Dec 10 16:53:41 2020 by Ohana Rodrigues using cetskelgen
// from cetlib version v3_05_01.
///////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/FindManyInChainP.h"

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// save info associated to common optical filter
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "ubevt/Utilities/PMTRemapService.h"
#include "ubevt/Utilities/PMTRemapProvider.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include <TTree.h>

#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

// For the NuMI trigger objects
#include "lardataobj/RawData/TriggerData.h" 
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h" 

//class muDARfeasibility;
class muDARfeasibility : public art::EDAnalyzer {

	public:
		explicit muDARfeasibility(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		muDARfeasibility(muDARfeasibility const&) = delete;
		muDARfeasibility(muDARfeasibility&&) = delete;
		muDARfeasibility& operator=(muDARfeasibility const&) = delete;
		muDARfeasibility& operator=(muDARfeasibility&&) = delete;

		// Required functions.
		void analyze(art::Event const& e) override;
		void beginJob() override;

		//POT count
		void endSubRun(art::SubRun const &subrun) override;

	private:

		//Truth information variabels
		TTree* _tree_Events;
		double _neutrino_E; //The neutrino energy 
		double _lepton_E; //The son lepton energy
		double _vertex_X; //The X position of the vertex 
		double _vertex_Y; //The Y position of the vertex   
		double _vertex_Z; //The Z position of the vertex 
	
		// POT information variabels  
		TTree *_subrun_tree;
  		int _run_sr; // The run number
		int _sub_sr; // The subRun number
  		float _pot;  // The total amount of POT for the current sub run
	
		// Trigger (common optical filtervariables)                                                        
  		float  _opfilter_pe_beam, _opfilter_pe_beam_tot, _opfilter_pe_veto, _opfilter_pe_veto_tot;

		// NuMI Trigger (common optical filtervariables)
		bool _swtrig_pre_ext, _swtrig_pre, _swtrig_post_ext, _swtrig_post;
};

muDARfeasibility::muDARfeasibility(fhicl::ParameterSet const& p)
	: EDAnalyzer{p}  
{
}

void muDARfeasibility::analyze(art::Event const& e){

	// ************************************************************************
	// *********************** Load MC truth info *****************************
	// ************************************************************************

	art::Handle<art::Assns<simb::MCTruth,simb::MCParticle,
	sim::GeneratedParticleInfo>> MCTruthHandle;
	art::InputTag fMCTruthTag("largeant");
	std::cout << "Is true getbylabel: " << 
	e.getByLabel(fMCTruthTag, MCTruthHandle) << std::endl;
	std::cout << "Size: " << MCTruthHandle->size() << std::endl;
	const art::Ptr<simb::MCTruth> mctruth = MCTruthHandle->at(0).first;
	std::cout << "mctruth nparticles  " << mctruth->NParticles() << std::endl;

	//Retrieving the neutrino information
	auto neutrino = mctruth->GetNeutrino();
    const simb::MCParticle nu = neutrino.Nu();
	_neutrino_E = nu.E();
	
	//Retrieving the lepton information
    const simb::MCParticle lepton = neutrino.Lepton();
	_lepton_E = lepton.E();

	//For the vertex...
	const TLorentzVector vertex = lepton.Position(); 
	_vertex_X = vertex.X(); 
 	_vertex_Y = vertex.Y(); 
	_vertex_Z = vertex.Z(); 
    
  // ***********************************************************************
  // ************************************************************************
  // ******** Trigger result output for NuMI software trigger *************** 
  // ******** Thanks to Owen who gave the code snippet for this! ************
  // ************************************************************************
  
	art::InputTag triggerTag ("swtrigger", "", "DataOverlayOpticalNuMI" );
  	const auto& triggerHandle = e.getValidHandle< raw::ubdaqSoftwareTriggerData >(triggerTag);
  	std::vector<std::string> triggerName = triggerHandle->getListOfAlgorithms();

  for (int j=0; j!=triggerHandle->getNumberOfAlgorithms(); j++){

	if (triggerName[j] == "EXT_NUMIwin_FEMBeamTriggerAlgo"){
    	_swtrig_pre_ext = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_pre_ext: " << _swtrig_pre_ext << std::endl;
    }

    else if (triggerName[j] == "EXT_NUMIwin_2018May_FEMBeamTriggerAlgo"){
    	_swtrig_post_ext = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_post_ext: " << _swtrig_post_ext << std::endl;
    }

    else if (triggerName[j] == "NUMI_FEMBeamTriggerAlgo"){
    	_swtrig_pre = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_pre: " << _swtrig_pre << std::endl;
    }

    else if (triggerName[j] == "NUMI_2018May_FEMBeamTriggerAlgo"){
    	_swtrig_post = triggerHandle->passedAlgo(triggerName[j]);
		std::cout << "_swtrig_post: " << _swtrig_post << std::endl;
    }

    else continue;

    // Print the trigger and the result
    std::cout<<triggerName[j]<<": ";
    std::cout<<triggerHandle->passedAlgo(triggerName[j])<<std::endl;
  }

	// ************************************************************************
	// *************** Load common-optical-filter output **********************
	// ************************************************************************
	art::Handle<uboone::UbooneOpticalFilter> CommonOpticalFilter_h;
	art::InputTag fCommonOpFiltTag("opfiltercommon");

	e.getByLabel(fCommonOpFiltTag, CommonOpticalFilter_h);

	_opfilter_pe_beam     = CommonOpticalFilter_h->PE_Beam();
	_opfilter_pe_beam_tot = CommonOpticalFilter_h->PE_Beam_Total();
	_opfilter_pe_veto     = CommonOpticalFilter_h->PE_Veto();
	_opfilter_pe_veto_tot = CommonOpticalFilter_h->PE_Veto_Total();

	// Implementation of required member function here.
	// Software Trigger
	art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    std::string fSoftwareTriggerModuleLabel = "swtrigger"; 
	
	if (e.getByLabel(fSoftwareTriggerModuleLabel, softwareTriggerHandle)){
		
		std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
		
		for (int i = 0; i < int(algoNames.size()); i++){

			if (algoNames[i] == "BNB_FEMBeamTriggerAlgo"){ 
			
				auto EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[i]); 

				std::cout << "EventPassedSwTrigger: " << EventPassedSwTrigger << std::endl;

				std::cout << "Found BNB_FEMBeamTrigger" << std::endl; 
			}
		}
	}
	else if (e.getByLabel("daq", softwareTriggerHandle)){

		std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();

		for (int i = 0; i < int(algoNames.size()); i++){
			
			if (algoNames[i] == "EXT_BNBwin_FEMBeamTriggerAlgo"||algoNames[i] == "BNB_FEMBeamTriggerAlgo"){
		 		
				auto EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[i]); 

				std::cout << "EventPassedSwTrigger: " << EventPassedSwTrigger << std::endl;

				std::cout << "Found EXT_BNBwin_FEMBeamTrigger and BNB_FEMBeamTrigger" << std::endl; 
			}
		}
	}

	_tree_Events->Fill();

}
//*************************************************************************
//***************** Getting POT normalization info ************************ 
//*************************************************************************
void muDARfeasibility::endSubRun(art::SubRun const &subrun){

 // bool isData=true; 
 // bool isFakeData=false;
  art::InputTag MCTproducer("generator");
  _pot=0; 

 // if((!isData) || (isFakeData)){
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    _pot = subrun.getByLabel(MCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.;
    std::cout << "[muDARfeasibility::endSubRun] Storing POT info!" << std::endl;
 // }

  _run_sr = subrun.run();
  _sub_sr = subrun.subRun();
  _subrun_tree->Fill();

}
//************************* end of block ***************************************

void muDARfeasibility::beginJob(){

	art::ServiceHandle<art::TFileService> tfs;

	_tree_Events = tfs->make<TTree>("Events", "Events");
	_subrun_tree = tfs->make<TTree>("SubRun", "SubRun");

	//trigger branches
	_tree_Events -> Branch("_opfilter_pe_beam", &_opfilter_pe_beam, "_opfilter_pe_beam/F");
	_tree_Events -> Branch("_opfilter_pe_beam_tot", &_opfilter_pe_beam_tot, "_opfilter_pe_beam_tot/F");
	_tree_Events -> Branch("_opfilter_pe_veto", &_opfilter_pe_veto, "_opfilter_pe_veto/F");
	_tree_Events -> Branch("_opfilter_pe_veto_tot", &_opfilter_pe_veto_tot, "_opfilter_pe_veto_tot/F");
	
	//NuMI trigger branches
	_tree_Events->Branch("_swtrig_pre", &_swtrig_pre, "_swtrig_pre/O");
	_tree_Events->Branch("_swtrig_pre_ext", &_swtrig_pre_ext, "_swtrig_pre_ext/O");
	_tree_Events->Branch("_swtrig_post_ext", &_swtrig_post_ext, "_swtrig_post_ext/O");
	_tree_Events->Branch("_swtrig_post", &_swtrig_post, "_swtrig_post/O");

	//truth info particles branches
	_tree_Events -> Branch("neutrino_E", &_neutrino_E, "neutrino_E/D");
	_tree_Events -> Branch("lepton_E", &_lepton_E, "lepton_E/D");
	_tree_Events -> Branch("vertex_X", &_vertex_X, "vertex_X/D");
	_tree_Events -> Branch("vertex_Y", &_vertex_Y, "vertex_Y/D");
	_tree_Events -> Branch("vertex_Z", &_vertex_Z, "vertex_Z/D");

	//POT counting info 
  	_subrun_tree->Branch("run", &_run_sr, "run/I");
  	_subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");
	_subrun_tree->Branch("pot", &_pot, "pot/F");
}

DEFINE_ART_MODULE(muDARfeasibility)

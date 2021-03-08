////////////////////////////////////////////////////////////////////////
// Class:       WireCellEventWeightTree
// Plugin Type: analyzer (art v3_01_02)
// File:        WireCellEventWeightTree_module.cc
//
// Generated at Sun Sept 27 2020 by Wenqiang Gu 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" 
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"
#include "ubobj/WcpPort/NuSelectionKINE.h"

#include "lardataobj/MCBase/MCShower.h"

#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>
#include <functional> // std::placeholders
#include <algorithm>  // std::transform

class WireCellEventWeightTree;


class WireCellEventWeightTree : public art::EDAnalyzer {
public:
  explicit WireCellEventWeightTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireCellEventWeightTree(WireCellEventWeightTree const&) = delete;
  WireCellEventWeightTree(WireCellEventWeightTree&&) = delete;
  WireCellEventWeightTree& operator=(WireCellEventWeightTree const&) = delete;
  WireCellEventWeightTree& operator=(WireCellEventWeightTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // user defined
  void reconfigure(fhicl::ParameterSet const& pset);
  void initOutput();
  void resetOutput();
  void save_weights(art::Event const& e);
  void save_LEEweights(art::Event const& e);
  void FillEventWeights(art::Event const & e, const bool partial = false); 

private:

  // Declare member data here.

  // fcl config
  std::string fSTMLabel;
  std::string fTruthLabel;
  std::string fFileType;
  std::vector<std::string> fMCEventWeightLabels;
  std::string fWeightLabel;
  std::string fWeightLeeLabel;
  bool fSaveWeights;
  bool fSaveLeeWeights;
  bool fSaveFullWeights;
  bool fIsNuMI;

  // output
  TTree* fTreeEval; 
  Int_t           f_run;
  Int_t           f_subRun;
  Int_t           f_event;
 
  Float_t 	  f_weight_spline;
  Float_t	  f_weight_cv;
  Float_t	  f_weight_lee;
  std::map<std::string, std::vector<float>> fmcweight;
  std::unordered_set<std::string> fgenie_knobs;
  Bool_t fmcweight_filled;
 
  Int_t		  f_stm_eventtype;
  Int_t		  f_stm_lowenergy;
  Int_t		  f_stm_LM;
  Int_t		  f_stm_TGM;
  Int_t		  f_stm_STM;
  Int_t		  f_stm_FullDead;
  Float_t	  f_stm_clusterlength;

  Int_t           f_truth_vtxInside;

};


WireCellEventWeightTree::WireCellEventWeightTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fgenie_knobs.insert("All_UBGenie");
  fgenie_knobs.insert("AxFFCCQEshape_UBGenie");
  fgenie_knobs.insert("DecayAngMEC_UBGenie");
  fgenie_knobs.insert("NormCCCOH_UBGenie");
  fgenie_knobs.insert("NormNCCOH_UBGenie");
  fgenie_knobs.insert("RPA_CCQE_UBGenie");
  fgenie_knobs.insert("RootinoFix_UBGenie");
  fgenie_knobs.insert("ThetaDelta2NRad_UBGenie");
  fgenie_knobs.insert("Theta_Delta2Npi_UBGenie");
  fgenie_knobs.insert("TunedCentralValue_UBGenie");
  fgenie_knobs.insert("VecFFCCQEshape_UBGenie");
  fgenie_knobs.insert("XSecShape_CCMEC_UBGenie");
  fgenie_knobs.insert("xsr_scc_Fa3_SCC");
  fgenie_knobs.insert("xsr_scc_Fv3_SCC");
  
  // fcl config
  reconfigure(p);

  // T_eval / event
  initOutput();

}

void WireCellEventWeightTree::reconfigure(fhicl::ParameterSet const& pset)
{
  std::cout<<"------------ WireCellEventWeightTree::reconfigure ----------"<<std::endl;  

  fSaveWeights = pset.get<bool>("SaveWeights", false); // GENIE weights
  fSaveLeeWeights = pset.get<bool>("SaveLeeWeights", false); // LEE weights
  fSaveFullWeights = pset.get<bool>("SaveFullWeights", false); // all weights
  fIsNuMI = pset.get<bool>("IsNuMI", false); // if true, converts NuMI's weights into BNB style
  fSTMLabel = pset.get<std::string>("STMLabel");
  fTruthLabel = pset.get<std::string>("TruthLabel");
  fFileType = pset.get<std::string>("FileType", "empty");
  fMCEventWeightLabels = pset.get<std::vector<std::string>>("MCEventWeightLabels");
  fWeightLabel = pset.get<std::string>("WeightLabel","");
  fWeightLeeLabel = pset.get<std::string>("WeightLeeLabel","");
}

void WireCellEventWeightTree::initOutput()
{
  std::cout<<"------------ WireCellEventWeightTree::initOutput ----------"<<std::endl;  

  art::ServiceHandle<art::TFileService> tfs;
  fTreeEval = tfs->make<TTree>("T_wgt", "T_wgt");
  
  fTreeEval->Branch("run", 			&f_run);
  fTreeEval->Branch("subrun", 			&f_subRun);
  fTreeEval->Branch("event", 			&f_event);
  fTreeEval->Branch("file_type", 		&fFileType);
  
  fTreeEval->Branch("weight_spline", 		&f_weight_spline); //MicroBooNE GENIE tune on top of weight_CV; weight_spline*weight_cv = weight
  fTreeEval->Branch("weight_cv",		&f_weight_cv); //MicroBooNE MCC9 untuned GENIE v3
  fTreeEval->Branch("weight_lee",		&f_weight_lee); //MicroBooNE MCC9 LEE weight 
  fTreeEval->Branch("mcweight", "std::map<std::string, std::vector<float>>" ,&fmcweight);  
  fTreeEval->Branch("mcweight_filled"		,&fmcweight_filled); // pass pattern recognition

}

void WireCellEventWeightTree::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // reset output at the beginning of each event
  resetOutput();

  // read NuMetrics
	std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;
	f_run = e.run();
 	f_subRun = e.subRun();
	f_event = e.id().event();

	art::Handle<std::vector<nsm::NuSelectionSTM> > stm_handle;
	e.getByLabel(fSTMLabel,stm_handle);
	std::vector<art::Ptr<nsm::NuSelectionSTM> > stm_vec;
	art::fill_ptr_vector(stm_vec,stm_handle);
	std::cout<<"--- NuSelectionSTM ---"<<std::endl;
	if(stm_vec.size()>1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	if(stm_vec.size()<1) {
		f_stm_eventtype = -1;
		f_stm_lowenergy = -1;
		f_stm_LM = -1;
		f_stm_TGM = -1;
		f_stm_STM = -1;
		f_stm_FullDead = -1;
		f_stm_clusterlength = -1.0;
	} 
	for(size_t i=0; i<stm_vec.size(); i++){
		art::Ptr<nsm::NuSelectionSTM> s = stm_vec.at(i);
		f_stm_eventtype = s->GetEventType();
		f_stm_lowenergy = s->GetLowEnergy();
		f_stm_LM = s->GetLM();
		f_stm_TGM = s->GetTGM();
		f_stm_STM = s->GetSTM();
		f_stm_FullDead = s->GetFullDead();
		f_stm_clusterlength = s->GetClusterLength();
	}

        f_truth_vtxInside = -1;
	art::Handle<std::vector<nsm::NuSelectionTruth> > truth_handle;
	e.getByLabel(fTruthLabel,truth_handle);
	std::vector<art::Ptr<nsm::NuSelectionTruth> > truth_vec;
	art::fill_ptr_vector(truth_vec,truth_handle);
	std::cout<<"--- NuSelectionTruth  ---"<<std::endl;
	if(truth_vec.size()!=1) {
		std::cout<<"WARNING: >1 in-beam matched TPC activity?!" << std::endl;
		return;
	} 
	for(size_t i=0; i<truth_vec.size(); i++){
		art::Ptr<nsm::NuSelectionTruth> t = truth_vec.at(i);
                f_truth_vtxInside = t->GetIsVtxInside();
	}

	/// save GENIE weights
	if(fSaveWeights) save_weights(e);	
	if(fSaveLeeWeights) save_LEEweights(e);
        if(fSaveFullWeights
           && f_stm_eventtype!=0 // only saves out weights
           && f_stm_lowenergy==0 // after pattern recoginition
           && f_stm_LM==0
           && f_stm_TGM==0
           && f_stm_STM==0
           && f_stm_FullDead==0
           && f_stm_clusterlength>0) { FillEventWeights(e); }
        else if (f_truth_vtxInside) { FillEventWeights(e, true); }
	/// end
	
	
	fTreeEval->Fill();
}

void WireCellEventWeightTree::resetOutput()
{
	// live period within each event
	// maybe redundant here
 	f_weight_spline = -1.0;
	f_weight_cv = -1.0;
	f_weight_lee = -1.0;
        fmcweight.clear();
        fmcweight_filled = false;
	
	f_stm_eventtype = -1;
	f_stm_lowenergy = -1;
	f_stm_LM = -1;
	f_stm_TGM = -1;
	f_stm_STM = -1;
	f_stm_FullDead = -1;
	f_stm_clusterlength = -1;
}

void WireCellEventWeightTree::save_weights(art::Event const& e)
{ 
  double ppfx_cv_UBPPFXCV = 1.0; // for NuMI

  // Use the EventWeight producer label here
  art::Handle<std::vector<evwgh::MCEventWeight> > weightsHandle;
  // e.getByLabel("eventweight", weightsHandle);
  e.getByLabel(fWeightLabel, weightsHandle);
  
  for(size_t i=0; i<weightsHandle->size(); i++){
    const evwgh::MCEventWeight& mc_weights = weightsHandle->at(i);
    // Loop over all of the weights in the MCEventWeight object
    for ( const auto& pair : mc_weights.fWeight ) {
      std::string knob_name = pair.first;
      std::vector<double> weights = pair.second;
      //std::cout<<"Knob name: "<<knob_name<<std::endl; 
      //std::cout<<"Weight size: "<<weights.size()<<std::endl; 

      if( knob_name == "TunedCentralValue_UBGenie"){
          f_weight_cv = weights.at(0);
      }
      if (knob_name == "splines_general_Spline"){
          f_weight_spline = weights.at(0);
      }
      if (knob_name == "ppfx_cv_UBPPFXCV" and fIsNuMI){
          double value = weights.at(0);
          if (not std::isnan(value) and not std::isinf(value)) {
            ppfx_cv_UBPPFXCV = value;
          }
      }
    }
  }

  if (fIsNuMI) {
    f_weight_spline *= ppfx_cv_UBPPFXCV; // absorb NuMI's cv correction into spline
    // std::cout << "weight_spline *= ppfx_cv_UBPPFXCV, where ppfx_cv_UBPPFXCV= " << ppfx_cv_UBPPFXCV << std::endl;
  }
  
  //std::cout<<"cv weight: "<<f_weight_cv<<std::endl;
  //std::cout<<"spline weight: "<<f_weight_spline<<std::endl;

}

void WireCellEventWeightTree::save_LEEweights(art::Event const& e)
{ 
  // Use the EventWeight producer label here
  art::Handle<std::vector<evwgh::MCEventWeight> > weightsHandle;
  // e.getByLabel("eventweightLEE", "", "EventWeightLEE", weightsHandle); // producer, instance, process
  e.getByLabel(fWeightLeeLabel, weightsHandle); // producer, instance, process
  
  // Loop through these objects for each neutrino vertex in the event
  for(size_t i=0; i<weightsHandle->size(); i++){
    const evwgh::MCEventWeight& mc_weights = weightsHandle->at(i);
    // Loop over all of the weights in the MCEventWeight object
    for ( const auto& pair : mc_weights.fWeight ) {
      std::string knob_name = pair.first;
      std::vector<double> weights = pair.second;
      //std::cout<<"Knob name: "<<knob_name<<std::endl; 
      //std::cout<<"Weight size: "<<weights.size()<<std::endl;
      if( knob_name == "eLEE_Combined_Oct2018_LEESignalElectron" ){
          f_weight_lee = weights.at(0);
      }
    }
  }
}

void WireCellEventWeightTree::FillEventWeights(art::Event const & e, const bool partial) {
  // To make the NuMI wights fit the BNB analysis, a few assumptions
  // have been made here:
  // 1. absorb ppfx_cv_UBPPFXCV into spline weight
  // 2. keep ppfx_ms_UBPPFX as it is, set all other ppfx_* vectors to unity
  // 3. re-map the 12 ppfx_* vectors to 12 BNB flux weights
  // 4. create a virtual knob with all values 1 in the vector (BNB has 13 flux knobs)
  double ppfx_cv_UBPPFXCV = 1.0; 
  bool got_ppfx_cv_UBPPFXCV = false;
  if (fIsNuMI) {
    for(auto const& weightLabel: fMCEventWeightLabels) {
      if (got_ppfx_cv_UBPPFXCV) break;
      auto ev_evw = e.getValidHandle<std::vector<evwgh::MCEventWeight>>(weightLabel);
      std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;
      if(ev_evw->size() > 1) {
        std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
    	      << "WARNING: eventweight has more than one entry\n";
      }
    
      for (auto const& x: weight_map) {
        std::string knob = x.first; // key
        std::vector<double> weights = x.second; // value
        if (knob == "ppfx_cv_UBPPFXCV") {
          double value = weights.at(0);
          if (not std::isnan(value) and not std::isinf(value)) {
            ppfx_cv_UBPPFXCV = value;
          }
          got_ppfx_cv_UBPPFXCV = true;
          break;
        }
      }
     
    }
    if (not got_ppfx_cv_UBPPFXCV) {
      std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
  	      << "WARNING: ppfx_cv_UBPPFXCV not found!\n";
    }
  }

  // Loop over weight maps
  for(auto const& weightLabel: fMCEventWeightLabels) {

    // art::ValidHandle<std::vector<evwgh::MCEventWeight>> const & ev_evw =
      // CHANGE HERE DEPENDING ON THE SAMPLES
      // e.getValidHandle<std::vector<evwgh::MCEventWeight>>("eventweight");
    auto ev_evw = e.getValidHandle<std::vector<evwgh::MCEventWeight>>(weightLabel);
    
    std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;
    if(ev_evw->size() > 1) {
      std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
  	      << "WARNING: eventweight has more than one entry\n";
    }
  
    // save wights in float
    for (auto const& x: weight_map) {
      std::string knob = x.first; // key
      std::vector<double> weights = x.second; // value
      std::vector<float> weights_asFloat(weights.begin(), weights.end());

      if (fIsNuMI) {

        if (knob.rfind("ppfx_", 0) == 0 and knob != "ppfx_cv_UBPPFXCV" and knob != "ppfx_ms_UBPPFX") {
          weights_asFloat = std::vector<float>(weights_asFloat.size(), 1.0);
        }

        if (knob == "splines_general_Spline") {
          std::transform(weights_asFloat.begin(), weights_asFloat.end(), weights_asFloat.begin(),
                         std::bind(std::multiplies<float>(), std::placeholders::_1, ppfx_cv_UBPPFXCV));
        }

        // The 12 knobs are re-maped to 13 BNB flux knobs:
        // ppfx_mippk_PPFXMIPPKaon 		=> expskin_FluxUnisim
        // ppfx_mipppi_PPFXMIPPPion 		=> horncurrent_FluxUnisim
        // ppfx_ms_UBPPFX 			=> kminus_PrimaryHadronNormalization
        // ppfx_other_PPFXOther 		=> kplus_PrimaryHadronFeynmanScaling
        // ppfx_targatt_PPFXTargAtten 		=> kzero_PrimaryHadronSanfordWang
        // ppfx_think_PPFXThinKaon 		=> nucleoninexsec_FluxUnisim
        // ppfx_thinmes_PPFXThinMeson 		=> nucleonqexsec_FluxUnisim
        // ppfx_thinn_PPFXThinNuc 		=> nucleontotxsec_FluxUnisim
        // ppfx_thinna_PPFXThinNucA 		=> piminus_PrimaryHadronSWCentralSplineVariation
        // ppfx_thinnpi_PPFXThinNeutronPion 	=> pioninexsec_FluxUnisim
        // ppfx_thinpi_PPFXThinPion 		=> pionqexsec_FluxUnisim
        // ppfx_totabs_PPFXTotAbsorp 		=> piontotxsec_FluxUnisim
        // vector<float> (1, ...) 		=> piplus_PrimaryHadronSWCentralSplineVariation

        if (knob == "ppfx_mippk_PPFXMIPPKaon") knob = "expskin_FluxUnisim";
        else if (knob == "ppfx_mipppi_PPFXMIPPPion") knob = "horncurrent_FluxUnisim";
        else if (knob == "ppfx_ms_UBPPFX") knob = "kminus_PrimaryHadronNormalization";
        else if (knob == "ppfx_other_PPFXOther") knob = "kplus_PrimaryHadronFeynmanScaling";
        else if (knob == "ppfx_targatt_PPFXTargAtten") knob = "kzero_PrimaryHadronSanfordWang";
        else if (knob == "ppfx_think_PPFXThinKaon") knob = "nucleoninexsec_FluxUnisim";
        else if (knob == "ppfx_thinmes_PPFXThinMeson") knob = "nucleonqexsec_FluxUnisim";
        else if (knob == "ppfx_thinn_PPFXThinNuc") knob = "nucleontotxsec_FluxUnisim";
        else if (knob == "ppfx_thinna_PPFXThinNucA") knob = "piminus_PrimaryHadronSWCentralSplineVariation";
        else if (knob == "ppfx_thinnpi_PPFXThinNeutronPion") knob = "pioninexsec_FluxUnisim";
        else if (knob == "ppfx_thinpi_PPFXThinPion") knob = "pionqexsec_FluxUnisim";
        else if (knob == "ppfx_totabs_PPFXTotAbsorp") knob = "piontotxsec_FluxUnisim";
      }

      if (partial && fgenie_knobs.find(knob) == fgenie_knobs.end()) {
          continue;
      }
      if (fmcweight.find(knob) == fmcweight.end()) {
        fmcweight.emplace(knob, weights_asFloat);
      }
      else { // a knob previously registered
        fmcweight[knob].reserve(fmcweight.size() + weights_asFloat.size());
        fmcweight[knob].insert(fmcweight[knob].end(), weights_asFloat.begin(), weights_asFloat.end());
      }

    }
   
    fmcweight_filled = true;
  }

  if (fIsNuMI) {
    int nsize = fmcweight["expskin_FluxUnisim"].size();
    std::vector<float> virtual_knob(nsize, 1.0);
    fmcweight.emplace("piplus_PrimaryHadronSWCentralSplineVariation", virtual_knob);
  }

}


DEFINE_ART_MODULE(WireCellEventWeightTree)

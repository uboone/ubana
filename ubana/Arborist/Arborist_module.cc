////////////////////////////////////////////////////////////////////////
// Class:       Arborist
// Plugin Type: analyzer (art v3_01_02)
// File:        Arborist_module.cc
//
// Generated at Mon Mar 18 11:32:48 2019 by Gray Yarbrough using cetskelgen
// from cetlib version v3_05_01.
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
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
// ART includes
#include "art/Framework/Services/Optional/TFileService.h" // used for ROOT file
#include "canvas/Persistency/Common/FindManyP.h" // used for assns

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT includes
#include "TTree.h"


class Arborist;


class Arborist : public art::EDAnalyzer {
public:
  explicit Arborist(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Arborist(Arborist const&) = delete;
  Arborist(Arborist&&) = delete;
  Arborist& operator=(Arborist const&) = delete;
  Arborist& operator=(Arborist&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void FillTruthInfo(art::Event const & e);
  void FillEventWeights(art::Event const & e);
  void resetVariables();


private:
  art::ServiceHandle< art::TFileService > tfs;

  art::InputTag fEventWeightInputTag;
  std::string fSplineBugFixLabel;
  std::string fRootinoBugFixLabel;
  std::string fTunedCentralValueLabel;
  std::string fLEESignalLabel;

  TTree* eventweight_tree;
  int run;
  int subrun;
  int event;
  int fNeutrinoPDGCode;
  double fTrueNeutrinoEnergy;
  int fInteractionCCNC;
  int fInteractionMode;
  int fInteractionType;
  int fTarget;
  double fTrueNeutrinoBaseline;
  double fSplineBugFixWeight;
  double fRootinoBugFixWeight;
  double fTunedCentralValueWeight;
  double fCombinedCentralValueWeight;
  double fLEESignalWeight;
  std::map<std::string, std::vector<double>> fEventWeightMap;

};


Arborist::Arborist(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fEventWeightInputTag(p.get<art::InputTag>("EventWeightInputTag")),
  fSplineBugFixLabel(p.get<std::string>("SplineBugFixLabel")),
  fRootinoBugFixLabel(p.get<std::string>("RootinoBugFixLabel")),
  fTunedCentralValueLabel(p.get<std::string>("TunedCentralValueLabel")),
  fLEESignalLabel(p.get<std::string>("LEESignalLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
void Arborist::beginJob()
{
 eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
 eventweight_tree->Branch("run", &run);
 eventweight_tree->Branch("subrun", &subrun);
 eventweight_tree->Branch("event", &event);
 eventweight_tree->Branch("nu_pdg", &fNeutrinoPDGCode);
 eventweight_tree->Branch("nu_energy_true", &fTrueNeutrinoEnergy);
 eventweight_tree->Branch("nu_interaction_ccnc", &fInteractionCCNC);
 eventweight_tree->Branch("nu_interaction_mode", &fInteractionMode);
 eventweight_tree->Branch("nu_interaction_type", &fInteractionType);
 eventweight_tree->Branch("nu_target_pdg", &fTarget);
 eventweight_tree->Branch("nu_L_true", &fTrueNeutrinoBaseline);
 eventweight_tree->Branch("spline_weight", &fSplineBugFixWeight);
 eventweight_tree->Branch("rootino_weight", &fRootinoBugFixWeight);
 eventweight_tree->Branch("ub_tune_weight", &fTunedCentralValueWeight);
 eventweight_tree->Branch("xsec_corr_weight", &fCombinedCentralValueWeight);
 eventweight_tree->Branch("lee_weight", &fLEESignalWeight);
 eventweight_tree->Branch("sys_weights", "std::map<std::string, std::vector<double>>", &fEventWeightMap);
}

void Arborist::resetVariables()
{
  run = -1;
  subrun = -1;
  event = -1;
  fNeutrinoPDGCode = 0;
  fTrueNeutrinoEnergy = -1.;
  fInteractionCCNC = -1;
  fInteractionMode = -1;
  fInteractionType = -1;
  fTarget = 0;
  fTrueNeutrinoBaseline = -1.;
  fSplineBugFixWeight = -1.;
  fRootinoBugFixWeight = -1.;
  fTunedCentralValueWeight = -1.;
  fCombinedCentralValueWeight = -1.;
  fLEESignalWeight = -1.;
  fEventWeightMap.clear();
}

void Arborist::FillTruthInfo(art::Event const & e) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth = 
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  fNeutrinoPDGCode = ev_mctruth->front().GetNeutrino().Nu().PdgCode();
  fTrueNeutrinoEnergy = ev_mctruth->front().GetNeutrino().Nu().E() * 1000.;  // convert from GeV to MeV
  fInteractionCCNC = ev_mctruth->front().GetNeutrino().CCNC();
  fInteractionMode = ev_mctruth->front().GetNeutrino().Mode();
  fInteractionType = ev_mctruth->front().GetNeutrino().InteractionType();
  fTarget = ev_mctruth->front().GetNeutrino().Target();

  art::ValidHandle<std::vector<simb::MCFlux>> const & ev_mcflux = 
    e.getValidHandle<std::vector<simb::MCFlux>>("generator");

  fTrueNeutrinoBaseline = ev_mcflux->front().fdk2gen + ev_mcflux->front().fgen2vtx;

}

void Arborist::FillEventWeights(art::Event const & e){
  
  std::vector<art::Handle<std::vector<evwgh::MCEventWeight>>> ev_evw_handles;
  e.getManyByType(ev_evw_handles);
  //std::cout << "INFO: ev_evw_handles.size() = " << ev_evw_handles.size() << std::endl;
  
  std::vector<evwgh::MCEventWeight> ev_evw_vec;
  for ( unsigned int i=0; i<ev_evw_handles.size(); i++ ) {
    if( !ev_evw_handles[i].isValid() ) {
      std::cout << "WARNING: eventweight product " << i << " is not valid!" << std::endl;
      continue;
    }
    ev_evw_vec.push_back(ev_evw_handles[i]->front());
    if(ev_evw_handles[i]->size() > 1) {
      std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
		<< "WARNING: eventweight product " << i << " has more than one entry\n";
    }
  }
  //std::cout << "INFO: ev_evw_vec.size() = " << ev_evw_vec.size() << std::endl;
  
  // Loop over eventweight products...
  for ( unsigned int i=0; i<ev_evw_vec.size(); i++ ) {

    std::map<std::string, std::vector<double>> const & weight_map = ev_evw_vec[i].fWeight;
    
    for ( auto evweight_it = weight_map.begin(); evweight_it != weight_map.end(); evweight_it++ ) {
      //std::cout << "INFO: eventweight product " << i << ", calculator " << evweight_it->first << " has " << evweight_it->second.size() << " universes" << std::endl;
      for ( unsigned int j=0; j<evweight_it->second.size(); j++ )
	fEventWeightMap[evweight_it->first].push_back(evweight_it->second[j]);
    }
    
    // Note: The following assumes...
    //   1) Any eventweight calculator with label fSplineBugFixLabel gives the desired fSplineBugFixWeight
    //   2) Any eventweight calculator with label fRootinoBugFixLabel gives the desired fRootinoBugFixWeight
    //   3) Any eventweight calculator with label fTunedCentralValueLabel gives the desired fTunedCentralValueWeight
    //   4) Any eventweight calculator with label fLEESignalLabel gives the desired fLEESignalWeight
    //   5) If the fCombinedCentralValueWeight is to be calculated, there will be at least one eventweight product with
    //        *all* of fSplineBugFixWeight, fRootinoBugFixWeight, and fTunedCentralValueWeight
    //   6) Any eventweight product with fSplineBugFixWeight, fRootinoBugFixWeight, and fTunedCentralValueWeight
    //        gives the desired fCombinedCentralValueWeight
    // ... this is fine with the standard Sep24 fcl chain, but may not be for other use cases
    
    if ( weight_map.find(fSplineBugFixLabel) != weight_map.end() ) {
      fSplineBugFixWeight = weight_map.find(fSplineBugFixLabel)->second[0];
      fEventWeightMap.erase(fSplineBugFixLabel);
    }
    if ( weight_map.find(fRootinoBugFixLabel) != weight_map.end() ) {
      fRootinoBugFixWeight = weight_map.find(fRootinoBugFixLabel)->second[0];
      fEventWeightMap.erase(fRootinoBugFixLabel);
    }
    if ( weight_map.find(fTunedCentralValueLabel) != weight_map.end() ) {
      fTunedCentralValueWeight = weight_map.find(fTunedCentralValueLabel)->second[0];
      fEventWeightMap.erase(fTunedCentralValueLabel);
    }
    if ( weight_map.find(fLEESignalLabel) != weight_map.end() ) {
      fLEESignalWeight = weight_map.find(fLEESignalLabel)->second[0];
      fEventWeightMap.erase(fLEESignalLabel);
    }
    
    if ( weight_map.find(fSplineBugFixLabel) != weight_map.end()
	 && weight_map.find(fRootinoBugFixLabel) != weight_map.end()
	 && weight_map.find(fTunedCentralValueLabel) != weight_map.end() )
      fCombinedCentralValueWeight = fSplineBugFixWeight * fRootinoBugFixWeight * fTunedCentralValueWeight;
    
  } // End loop over eventweight products
  
  /*
  for ( auto evweight_it = fEventWeightMap.begin(); evweight_it != fEventWeightMap.end(); evweight_it++ )
    std::cout << "INFO: in total, calculator " << evweight_it->first << " has " << evweight_it->second.size() << " universes" << std::endl;
  */

}

void Arborist::analyze(art::Event const& e)
{
  resetVariables();

  run = e.run();
  subrun = e.subRun();
  event = e.event();

  FillTruthInfo(e);
  FillEventWeights(e);
  
  eventweight_tree->Fill();
}

DEFINE_ART_MODULE(Arborist)

#ifndef VALIDATETRIGGERSLIDING_MODULE
#define VALIDATETRIGGERSLIDING_MODULE

#include "ValidateTriggerSliding.h"

ValidateTriggerSliding::ValidateTriggerSliding(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor ValidateTriggerSliding
ValidateTriggerSliding::~ValidateTriggerSliding() {}
void ValidateTriggerSliding::endJob() {}

void ValidateTriggerSliding::beginJob()
{
  _tDataTree = _tfs->make<TTree>("Data","");
  _tDataTree->Branch("run",&_run);
  _tDataTree->Branch("subrun",&_subrun);
  _tDataTree->Branch("event",&_event);
  _tDataTree->Branch("triggerName",&_triggerName);
  _tDataTree->Branch("triggerAlgoPass",&_triggerAlgoPass);
  _tDataTree->Branch("triggerPrescalePass",&_triggerPrescalePass);
  _tDataTree->Branch("triggerPass",&_triggerPass);
}

void ValidateTriggerSliding::ClearVectors()
{
  _triggerName.clear();
  _triggerAlgoPass.clear();
  _triggerPrescalePass.clear();
  _triggerPass.clear();
}

// Example function to get handle and associations
void ValidateTriggerSliding::ExtractTriggerInformation(art::Event const & evt)
{
  // Temporary variables
  bool pass_algo, pass_prescale, pass;

  std::vector<std::string> triggerLabel = {
    std::string("daq"),
    std::string("swtrigger0"),
    std::string("swtrigger1"),
    std::string("swtrigger2"),
    std::string("swtrigger3"),
    std::string("swtrigger4"),
    std::string("swtrigger5"),
    std::string("swtrigger6")
  };

  for (int i=0; i!=int(triggerLabel.size()); i++)
  {
    // Prepare trigger information
    art::InputTag triggerTag {triggerLabel[i]};
    const auto& triggerHandle = evt.getValidHandle< raw::ubdaqSoftwareTriggerData >(triggerTag);
    std::vector<std::string> triggerName = triggerHandle->getListOfAlgorithms();
    for (int j=0; j!=triggerHandle->getNumberOfAlgorithms(); j++)
    {
      // Check trigger
      pass_algo = triggerHandle->passedAlgo(triggerName[j]);
      pass_prescale = triggerHandle->passedPrescaleAlgo(triggerName[j]);
      pass = pass_algo && pass_prescale;
      // Assign values
      _triggerName.push_back(triggerLabel[i]+std::string("_")+triggerName[j]);
      _triggerAlgoPass.push_back(pass_algo);
      _triggerPrescalePass.push_back(pass_prescale);
      _triggerPass.push_back(pass);
    }
  }
}

void ValidateTriggerSliding::analyze(art::Event const & evt)
{
  // Determine event ID
  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();
  ClearVectors();

  ExtractTriggerInformation(evt);
  _tDataTree->Fill();
} // END function analyze

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ValidateTriggerSliding)

#endif // END def ValidateTriggerSliding_module

/******************************************************************************
 * @file TriggerInformationAlg.cxx
 * @brief Retrieve information regarding which triggers have been passed or not, then fill the candidate tree. The module calls only the event, not the actual candidate (e.g. there is redundant information, since if there are multiple candidates in an event they're still going to have identical entries).
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see TriggerInformationAlg.h
 * ****************************************************************************/
#include "TriggerInformationAlg.h"

namespace TriggerInformation
{
  // Constructor/destructor
  TriggerInformationAlg::TriggerInformationAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
  }
  TriggerInformationAlg::~TriggerInformationAlg()
  {}
  void TriggerInformationAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fIsData = pset.get<bool>("IsData");
  }

  // Add trigger information to tree
  void TriggerInformationAlg::AddTriggerInformation(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf)
  {
    // Temporary variables
    bool pass_algo, pass_prescale, pass;

    // Prepare trigger information
    std::string triggerLabel;
    if (fIsData) {triggerLabel = std::string("daq");}
    else {triggerLabel = std::string("swtrigger");}
    art::InputTag triggerTag {triggerLabel};
    const auto& triggerHandle = evt.getValidHandle< raw::ubdaqSoftwareTriggerData >(triggerTag);
    
    std::vector<std::string> triggerName = triggerHandle->getListOfAlgorithms();
    for (int j=0; j!=triggerHandle->getNumberOfAlgorithms(); j++)
    {
      // Check trigger
      pass_algo = triggerHandle->passedAlgo(triggerName[j]);
      pass_prescale = triggerHandle->passedPrescaleAlgo(triggerName[j]);
      pass = pass_algo && pass_prescale;
      // Assign values
      ctf.trigger_triggerName.push_back(triggerLabel+std::string("_")+triggerName[j]);
      ctf.trigger_triggerAlgoPass.push_back(pass_algo);
      ctf.trigger_triggerPrescalePass.push_back(pass_prescale);
      ctf.trigger_triggerPass.push_back(pass);
    }
  } //  END function AddTriggerInformation

} // END namespace TriggerInformation

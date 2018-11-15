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
  {}

  // Add trigger information to tree
  void TriggerInformationAlg::AddTriggerInformation(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf)
  {
    // Temporary variables
    bool pass_algo, pass_prescale, pass;

    // Prepare trigger information
    std::string triggerLabel("swtrigger");
    art::InputTag triggerTag {triggerLabel};
    const auto& triggerHandle = evt.getValidHandle< raw::ubdaqSoftwareTriggerData >(triggerTag);
    
    // Check trigger
    std::string bnbTrigName("BNB_FEMBeamTriggerAlgo");
    pass_algo = triggerHandle->passedAlgo(bnbTrigName);
    pass_prescale = triggerHandle->passedPrescaleAlgo(bnbTrigName);
    pass = pass_algo && pass_prescale;
    ctf.trigger_passedBNB = pass;

    std::string bnb5PETrigName("BNB_2017Dec_SwTrigger5PE_FEMBeamTriggerAlgo");
    pass_algo = triggerHandle->passedAlgo(bnb5PETrigName);
    pass_prescale = triggerHandle->passedPrescaleAlgo(bnb5PETrigName);
    pass = pass_algo && pass_prescale;
    ctf.trigger_passedBNB5PE = pass;

    std::string extTrigName("EXT_BNBwin_FEMBeamTriggerAlgo");
    pass_algo = triggerHandle->passedAlgo(extTrigName);
    pass_prescale = triggerHandle->passedPrescaleAlgo(extTrigName);
    pass = pass_algo && pass_prescale;
    ctf.trigger_passedEXT = pass;

    std::string ext5PETrigName("EXT_BNBwin_2017Dec_SwTrigger5PE_FEMBeamTriggerAlgo");
    pass_algo = triggerHandle->passedAlgo(ext5PETrigName);
    pass_prescale = triggerHandle->passedPrescaleAlgo(ext5PETrigName);
    pass = pass_algo && pass_prescale;
    ctf.trigger_passedEXT5PE = pass;

    std::string hsnTrigName("BNB_HSN_c0_FEMBeamTriggerAlgo");
    pass_algo = triggerHandle->passedAlgo(hsnTrigName);
    pass_prescale = triggerHandle->passedPrescaleAlgo(hsnTrigName);
    pass = pass_algo && pass_prescale;
    ctf.trigger_passedHSN = pass;

    std::string hsnExtTrigName("EXT_HSN_c0_FEMBeamTriggerAlgo");
    pass_algo = triggerHandle->passedAlgo(hsnExtTrigName);
    pass_prescale = triggerHandle->passedPrescaleAlgo(hsnExtTrigName);
    pass = pass_algo && pass_prescale;
    ctf.trigger_passedHSNEXT = pass;
  } //  END function AddTriggerInformation

} // END namespace TriggerInformation

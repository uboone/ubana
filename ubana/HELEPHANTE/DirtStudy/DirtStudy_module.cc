// DEFAULT START BLOCK
#ifndef EXPLORETRUTH_MODULE
#define EXPLORETRUTH_MODULE
#include "DirtStudy.h" 
DirtStudy::DirtStudy(fhicl::ParameterSet const & pset) :
  EDAnalyzer(pset),
  fPartEnergyThreshold(pset.get<float>("PartEnergyThreshold")),
  fPhotEnergyThreshold(pset.get<float>("PhotEnergyThreshold"))
{} // END constructor DirtStudy
DirtStudy::~DirtStudy()
{} // END destructor DirtStudy
void DirtStudy::endJob()
{} // END function endJob



// BEGIN JOB FUNCTION
void DirtStudy::beginJob()
{
  art::ServiceHandle< art::TFileService > tfs;
  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&_run);
  tDataTree->Branch("subrun",&_subrun);
  tDataTree->Branch("event",&_event);
  tDataTree->Branch("trigger_passedBNB",&_trigger_passedBNB);
  tDataTree->Branch("trigger_passedBNB5PE",&_trigger_passedBNB5PE);
  tDataTree->Branch("trigger_passedEXT",&_trigger_passedEXT);
  tDataTree->Branch("trigger_passedEXT5PE",&_trigger_passedEXT5PE);
  tDataTree->Branch("trigger_passedHSN",&_trigger_passedHSN);
  tDataTree->Branch("trigger_passedHSNEXT",&_trigger_passedHSNEXT);
  tDataTree->Branch("neutrino_position",&_neutrino_position);
  tDataTree->Branch("neutrino_time",&_neutrino_time);
  tDataTree->Branch("neutrino_time",&_neutrino_time);
  tDataTree->Branch("neutrino_energy",&_neutrino_energy);
  tDataTree->Branch("particles_time",&_particles_time);
  tDataTree->Branch("particles_energy",&_particles_energy);
  tDataTree->Branch("photons_time",&_photons_time);
  // tDataTree->Branch("photons_energy",&_photons_energy);
} // END function beginJob



// FUNCTION AddG4Information
void DirtStudy::AddG4TimingInformation(art::Event const & evt)
{
  _particles_time.clear();
  _particles_energy.clear();
  _photons_time.clear();
  _photons_energy.clear();
  _neutrino_position = {-999.,-999.,-999.};


  // Get energy and time from neutrino
  art::InputTag truthTag {std::string("generator")};
  auto const & mcTruthHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(truthTag);
  art::Ptr<simb::MCTruth> mct(mcTruthHandle,0);
  const simb::MCParticle & mcNu = mct->GetNeutrino().Nu();
  _neutrino_time = mcNu.T();
  _neutrino_energy = mcNu.E();
  _neutrino_position = {(float) mcNu.Vx(),(float) mcNu.Vy(),(float) mcNu.Vz()};

  // Get energy and time from all particles
  art::InputTag mcParticleTag {std::string("largeant")};
  auto const & mcParticleHandle = evt.getValidHandle<std::vector<simb::MCParticle>>(mcParticleTag);
  for(int i=0; i!=int((*mcParticleHandle).size()); i++)
  {
    art::Ptr<simb::MCParticle> mcp(mcParticleHandle,i);
    if(mcp->E() > fPartEnergyThreshold)
    {
      _particles_time.push_back(mcp->T());
      _particles_energy.push_back(mcp->E());
    }
  }

  // Get energy and time from all photons
  art::InputTag simPhotonsTag {std::string("largeant")};
  auto const & simPhotonsHandle = evt.getValidHandle<std::vector<sim::SimPhotons>>(simPhotonsTag);
  for(int i=0; i!=int((*simPhotonsHandle).size()); i++)
  {
    art::Ptr<sim::SimPhotons> simp(simPhotonsHandle,i);
    for(auto const& photon : (*simp))
    {
      if(photon.Energy > fPhotEnergyThreshold)
      {
        _photons_time.push_back(photon.Time);
        _photons_energy.push_back(photon.Energy);
      } // END if on energy
    } // END loop on photons
  } // END loop on optical channels
} // END function AddG4Information



// FUNCTION AddTriggerInformation
void DirtStudy::AddTriggerInformation(art::Event const & evt)
{
  // Temporary variables
  bool pass_algo, pass_prescale, pass;

  // Prepare trigger information
  art::InputTag triggerTag {std::string("swtrigger")};
  const auto& triggerHandle = evt.getValidHandle< raw::ubdaqSoftwareTriggerData >(triggerTag);

  // Check trigger
  std::string bnbTrigName("BNB_FEMBeamTriggerAlgo");
  pass_algo = triggerHandle->passedAlgo(bnbTrigName);
  pass_prescale = triggerHandle->passedPrescaleAlgo(bnbTrigName);
  pass = pass_algo && pass_prescale;
  _trigger_passedBNB = pass;

  std::string bnb5PETrigName("BNB_2017Dec_SwTrigger5PE_FEMBeamTriggerAlgo");
  pass_algo = triggerHandle->passedAlgo(bnb5PETrigName);
  pass_prescale = triggerHandle->passedPrescaleAlgo(bnb5PETrigName);
  pass = pass_algo && pass_prescale;
  _trigger_passedBNB5PE = pass;

  std::string extTrigName("EXT_BNBwin_FEMBeamTriggerAlgo");
  pass_algo = triggerHandle->passedAlgo(extTrigName);
  pass_prescale = triggerHandle->passedPrescaleAlgo(extTrigName);
  pass = pass_algo && pass_prescale;
  _trigger_passedEXT = pass;

  std::string ext5PETrigName("EXT_BNBwin_2017Dec_SwTrigger5PE_FEMBeamTriggerAlgo");
  pass_algo = triggerHandle->passedAlgo(ext5PETrigName);
  pass_prescale = triggerHandle->passedPrescaleAlgo(ext5PETrigName);
  pass = pass_algo && pass_prescale;
  _trigger_passedEXT5PE = pass;

  std::string hsnTrigName("BNB_HSN_c0_FEMBeamTriggerAlgo");
  pass_algo = triggerHandle->passedAlgo(hsnTrigName);
  pass_prescale = triggerHandle->passedPrescaleAlgo(hsnTrigName);
  pass = pass_algo && pass_prescale;
  _trigger_passedHSN = pass;

  std::string hsnExtTrigName("EXT_HSN_c0_FEMBeamTriggerAlgo");
  pass_algo = triggerHandle->passedAlgo(hsnExtTrigName);
  pass_prescale = triggerHandle->passedPrescaleAlgo(hsnExtTrigName);
  pass = pass_algo && pass_prescale;
  _trigger_passedHSNEXT = pass; 
} // END function AddTriggerInformation



// MAIN FUNCTION ANALYZE
void DirtStudy::analyze(art::Event const & evt)
{
  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();
  AddG4TimingInformation(evt);
  AddTriggerInformation(evt);

  tDataTree->Fill();
} // END function analyze



// DEFAULT END BLOCK
DEFINE_ART_MODULE(DirtStudy)
#endif

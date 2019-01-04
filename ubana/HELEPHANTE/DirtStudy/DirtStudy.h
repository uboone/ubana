#ifndef SINGLEPARTICLEMOMENTUMSTUDY_H
#define SINGLEPARTICLEMOMENTUMSTUDY_H
// c++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <exception>
// root includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
// framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
// art includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
// larsoft object includes
#include "lardataobj/RawData/TriggerData.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimPhotons.h"



// Analyzer class
class DirtStudy : public art::EDAnalyzer
{
public:
  explicit DirtStudy(fhicl::ParameterSet const & pset);
  virtual ~DirtStudy();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void AddTriggerInformation(art::Event const & evt);
  void AddG4TimingInformation(art::Event const & evt);

private:
  // Declare trees and tree variables
  art::ServiceHandle< art::TFileService > tfs;
  TTree *tDataTree;
  float fPartEnergyThreshold, fPhotEnergyThreshold;
  int _run, _subrun, _event;
  bool _trigger_passedBNB, _trigger_passedBNB5PE, _trigger_passedEXT, _trigger_passedEXT5PE, _trigger_passedHSN, _trigger_passedHSNEXT;
  float _neutrino_time, _neutrino_energy;
  std::vector<float> _neutrino_position = {-999,-999,-999};
  std::vector<float> _particles_time, _particles_energy, _photons_time, _photons_energy;

}; // End class DirtStudy
#endif // END def DirtStudy header
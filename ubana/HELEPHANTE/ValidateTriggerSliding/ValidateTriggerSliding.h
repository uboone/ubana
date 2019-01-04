#ifndef QUICKEXPLORER_H
#define QUICKEXPLORER_H

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
#include "TH1.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraph.h"

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


// Analyzer class
class ValidateTriggerSliding : public art::EDAnalyzer
{
public:
  explicit ValidateTriggerSliding(fhicl::ParameterSet const & pset);
  virtual ~ValidateTriggerSliding();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void ClearVectors();
  void ExtractTriggerInformation(art::Event const & evt);
private:
  // Declare trees and tree variables
  art::ServiceHandle< art::TFileService > _tfs;
  TTree *_tDataTree;
  int _run, _subrun, _event;
  std::vector<std::string> _triggerName;
  std::vector<bool> _triggerAlgoPass, _triggerPrescalePass, _triggerPass;
}; // End class ValidateTriggerSliding

#endif // END def ValidateTriggerSliding header
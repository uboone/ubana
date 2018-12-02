#ifndef DEDXALG_H
#define DEDXALG_H

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
#include "canvas/Persistency/Common/Assns.h"
// larsoft object includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
// Auxiliary objects includes
#include "ubana/HELEPHANTE/HelephantFinder/DataObjects/DecayVertex.h"
#include "ubana/HELEPHANTE/HelephantFinder/DataObjects/EventTreeFiller.h"
#include "ubana/HELEPHANTE/HelephantFinder/DataObjects/DrawTreeFiller.h"
#include "ubana/HELEPHANTE/HelephantFinder/DataObjects/CandidateTreeFiller.h"


namespace dEdXInformation
{

  class dEdXInformationAlg
  {
  public:
    dEdXInformationAlg(fhicl::ParameterSet const & pset);
    ~dEdXInformationAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

    // Algorithms
    void AdddEdXInformation(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf,
            AuxVertex::DecayVertex const & dv);
    bool IsTrackFlipped(art::Ptr<anab::Calorimetry> const & calo);
  private:
  // Data product labels settings
  std::string fCalorimetryLabel, fCalibratedLabel;
  bool fSaveAlsoUncalibrateddEdXInformation;
  };

} // END namespace dEdXInformation

#endif

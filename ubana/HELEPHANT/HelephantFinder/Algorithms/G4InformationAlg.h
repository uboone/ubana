#ifndef G4INFORMATIONALG_H
#define G4INFORMATIONALG_H

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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// Auxiliary objects includes
#include "ubana/HELEPHANT/HelephantFinder/DataObjects/DecayVertex.h"
#include "ubana/HELEPHANT/HelephantFinder/DataObjects/EventTreeFiller.h"
#include "ubana/HELEPHANT/HelephantFinder/DataObjects/DrawTreeFiller.h"
#include "ubana/HELEPHANT/HelephantFinder/DataObjects/CandidateTreeFiller.h"


namespace G4Information
{

  class G4InformationAlg
  {
  public:
    G4InformationAlg(fhicl::ParameterSet const & pset);
    ~G4InformationAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

    // Algorithms
    void AddG4Information(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf,
            AuxVertex::DecayVertex const & dv);
    float McpLength(art::Ptr<simb::MCParticle> mcp);
    void BadTruthRecoMatching(
            std::vector< art::Ptr<simb::MCParticle> > const & mcps,
            AuxVertex::DecayVertex const & dv,
            AuxEvent::CandidateTreeFiller & ctf);
    float End2EndDistance(
            int prongID,
            art::Ptr<simb::MCParticle> const & mcp,
            AuxVertex::DecayVertex const & dv);
  private:
  };

} // END namespace G4Information

#endif

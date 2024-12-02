#ifndef RECOTRUTHDISTANCEALG_H
#define RECOTRUTHDISTANCEALG_H

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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
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
// Auxiliary objects includes
#include "ubana/HSNAnalysis/HsnFinder/DataObjects/DecayVertex.h"
#include "ubana/HSNAnalysis/HsnFinder/DataObjects/EventTreeFiller.h"


namespace RecoTruthDistance
{

  class RecoTruthDistanceAlg
  {
  public:
    RecoTruthDistanceAlg(fhicl::ParameterSet const & pset);
    ~RecoTruthDistanceAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

    // Algorithms
    void DetermineRecoTruthDistance(
            art::Event const & evt,
            AuxEvent::EventTreeFiller & etf,
            std::vector<AuxVertex::DecayVertex> & ana_decayVertices);
  private:
    bool fVerbose;
    // microboone services
    const geo::GeometryCore* fGeometry;
  };

} // END namespace RecoTruthDistance

#endif

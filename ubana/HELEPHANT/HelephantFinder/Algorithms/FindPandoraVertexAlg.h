#ifndef FINDPANDORAVERTEXALG_H
#define FINDPANDORAVERTEXALG_H

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



namespace FindPandoraVertex
{

  class FindPandoraVertexAlg
  {
  public:
    FindPandoraVertexAlg(fhicl::ParameterSet const & pset);
    ~FindPandoraVertexAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

    // Algorithms
    void GetPotentialNeutrinoVertices(
            art::Event const & evt,
            AuxEvent::EventTreeFiller & evd,
            std::vector<AuxVertex::DecayVertex> & ana_decayVertices);

  private:
    // fhicl parameters
    std::string fPfpLabel;
    std::string fMcsLabel;
    std::vector<double> fMinTpcBound;
    std::vector<double> fMaxTpcBound;
    bool fVerbose;

    // microboone services
    const geo::GeometryCore* fGeometry;
    const detinfo::DetectorProperties* fDetectorProperties;
  };

} // END namespace FindPandoraVertex

#endif

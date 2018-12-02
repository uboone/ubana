#ifndef EXPLORETRUTH_H
#define EXPLORETRUTH_H

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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/MCBase/MCTrack.h"
// #include "lardataobj/RecoBase/Track.h"
// #include "lardataobj/RecoBase/Shower.h"
// #include "lardataobj/RecoBase/Vertex.h"
// #include "lardataobj/RecoBase/PFParticle.h"
// #include "lardataobj/RecoBase/Wire.h"
// #include "lardataobj/RecoBase/Hit.h"
// #include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubobj/RawData/SparseRawDigit.h"


// Analyzer class
class ExploreTruth : public art::EDAnalyzer
{
public:
  explicit ExploreTruth(fhicl::ParameterSet const & pset);
  virtual ~ExploreTruth();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void PrintInfo(art::Ptr<simb::MCParticle> mcp, std::string type, std::string spaces);
private:
  // Declare trees and tree variables
  TTree *tDataTree;
}; // End class ExploreTruth

#endif // END def ExploreTruth header
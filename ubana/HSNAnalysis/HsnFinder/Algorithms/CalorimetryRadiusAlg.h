#ifndef CALORIMETRYRADIUSALG_H
#define CALORIMETRYRADIUSALG_H

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

// larsoft object includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Auxiliary objects includes
#include "ubana/HSNAnalysis/HsnFinder/DataObjects/DecayVertex.h"
#include "ubana/HSNAnalysis/HsnFinder/DataObjects/EventTreeFiller.h"

namespace CalorimetryRadius
{

  class CalorimetryRadiusAlg
  {
  public:
    CalorimetryRadiusAlg(fhicl::ParameterSet const & pset);
    ~CalorimetryRadiusAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

  // Algorithms
  void PerformCalorimetry(
          art::Event const & evt,
          AuxEvent::EventTreeFiller & evd,
          std::vector<AuxVertex::DecayVertex>& decayVertices);

  // PerformCalorimetry returns
  std::vector<std::vector<art::Ptr<recob::Hit>>> ana_calo_totHitsInMaxRadius; // for each hit, for each dv

  int tree_calo_NumTotHits;
  std::vector<std::vector<float>> tree_calo_totChargeInRadius;
  std::vector<std::vector<float>> tree_calo_prong1ChargeInRadius;
  std::vector<std::vector<float>> tree_calo_prong2ChargeInRadius;
  std::vector<std::vector<float>> tree_calo_caloRatio;

  private:
    // fhicl parameters
    std::string fPfpLabel;
    std::string fHitLabel;
    std::vector<double> fRadiusProfileLimits;
    int fRadiusProfileBins;
    double fChannelNorm;
    double fTickNorm;
    bool fVerbose;
    std::vector<float> profileTicks;

    // microboone services
    const geo::GeometryCore* fGeometry;
  };

} // END namespace CalorimetryRadius

#endif
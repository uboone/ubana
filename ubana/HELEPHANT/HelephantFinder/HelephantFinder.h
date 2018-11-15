#ifndef HELEPHANTFINDER_H
#define HELEPHANTFINDER_H

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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Auxiliary objects includes
#include "Algorithms/FindPandoraVertexAlg.h"
#include "Algorithms/CalorimetryRadiusAlg.h"
#include "Algorithms/ExtractTruthInformationAlg.h"
#include "Algorithms/FlashMatchingAlg.h"
#include "Algorithms/TriggerInformationAlg.h"
#include "Algorithms/G4InformationAlg.h"
#include "DataObjects/DecayVertex.h"
#include "DataObjects/EventTreeFiller.h"
#include "DataObjects/CandidateTreeFiller.h"
#include "DataObjects/DrawTreeFiller.h"



// Analyzer class
class HelephantFinder : public art::EDAnalyzer
{
public:
  explicit HelephantFinder(fhicl::ParameterSet const & pset);
  virtual ~HelephantFinder();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Algorithms
  FindPandoraVertex::FindPandoraVertexAlg fFindPandoraVertexAlg;
  CalorimetryRadius::CalorimetryRadiusAlg fCalorimetryRadiusAlg;
  ExtractTruthInformation::ExtractTruthInformationAlg fExtractTruthInformationAlg;
  FlashMatching::FlashMatchingAlg fFlashMatchingAlg;
  TriggerInformation::TriggerInformationAlg fTriggerInformationAlg;
  G4Information::G4InformationAlg fG4InformationAlg;
  // Meta information
  std::string fInstanceName;
  int fIteration;
  std::vector<double> fRadiusProfileLimits;
  int fRadiusProfileBins;
  double fChannelNorm, fTickNorm;
  bool fVerbose;
  // Geometry settings
  std::vector<double> fMinTpcBound, fMaxTpcBound, fCenterCoordinates;
  // Data product labels settings
  std::string fPfpLabel, fHitLabel, fMcsLabel, fMcTrackLabel, fFlashLabel;
  // Additional tree information
  bool fSaveDrawTree, fSaveTruthDrawTree, fSaveTriggerInformation, fSaveG4Information;
  // Data/MC type dependent settings
  bool fUseTruthDistanceMetric, fIsHSN;

  // Declare services
  geo::GeometryCore const* fGeometry; // Pointer to the Geometry service
  detinfo::DetectorProperties const* fDetectorProperties; // Pointer to the Detector Properties

  // Declare trees
  TTree *metaTree;
  TTree *eventTree;
  TTree *candidateTree;
  TTree *drawTree;

  // Declare tree fillers
  AuxEvent::EventTreeFiller etf;
  AuxEvent::CandidateTreeFiller ctf;
  AuxEvent::DrawTreeFiller dtf;

  // Declare analysis variables
  std::vector<float> profileTicks;

  // Declare pandora analysis variables
  std::vector<AuxVertex::DecayVertex> ana_pandora_decayVertices;
  std::vector<recob::PFParticle const*> ana_pandora_neutrinos, ana_pandora_tracks, ana_pandora_showers;

  // Declare analysis functions
  void ClearData();
}; // End class HelephantFinder

#endif

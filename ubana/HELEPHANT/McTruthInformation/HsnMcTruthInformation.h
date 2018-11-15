#ifndef HSNMCTRUTHINFORMATION_H
#define HSNMCTRUTHINFORMATION_H

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
#include "lardataobj/Simulation/SimChannel.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Analyzer class
class HsnMcTruthInformation : public art::EDAnalyzer
{
public:
  explicit HsnMcTruthInformation(fhicl::ParameterSet const & pset);
  virtual ~HsnMcTruthInformation();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void GetTruthParticles(art::Event const & evt);
private:

  // Declare fhiclcpp variables
  std::string fMcTruthLabel;
  std::string fMcTrackLabel;

  // Declare trees and tree variables
  TTree *tDataTree;
  std::vector<int> pdgCode;
  std::vector<float> Vx, Vy, Vz, T, EndX, EndY, EndZ, EndT, Px, Py, Pz, E, P, Pt, Length, Theta, Phi;
  float Nu_E, Nu_Px, Nu_Py, Nu_Pz, Nu_P, Nu_Theta, Nu_Phi;
  float OpeningAngle, InvariantMass;
  bool Contained;


  // Declare analysis variables
  int run, subrun, event;

  // Declare analysis functions
  void ClearData();
  float TrackLength(std::vector<float> start, std::vector<float> end);
}; // End class HsnMcTruthInformation

#endif // END def HsnMcTruthInformation header
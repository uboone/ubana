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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"


// Analyzer class
class SingleParticleMomentumStudy : public art::EDAnalyzer
{
public:
  explicit SingleParticleMomentumStudy(fhicl::ParameterSet const & pset);
  virtual ~SingleParticleMomentumStudy();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void FillEndProcess(
          art::Ptr<simb::MCParticle> const & mcp,
          std::map< int,art::Ptr<simb::MCParticle> > mcpMap);
  std::map< int,art::Ptr<simb::MCParticle> > BuildMcpMap(
        art::ValidHandle< std::vector<simb::MCParticle> > const & mcHandle);
  std::map< int,art::Ptr<sim::MCTrack> > BuildMctMap(
        art::ValidHandle< std::vector<sim::MCTrack> > const & mctHandle);
  art::Ptr<recob::Track> FindCorrectRecoTrack(
          art::ValidHandle< std::vector<recob::Track> > const & recoHandle,
          art::Ptr<sim::MCTrack> const & primary_mct);
  float RecoTruthScore(
        art::Ptr<sim::MCTrack> const & primary_mct,
        art::Ptr<recob::Track> const & reco_candidate);
  float MCTrackLength(art::Ptr<sim::MCTrack> const & mct);
  float TrackLength(art::Ptr<recob::Track> const & recot);

private:
  // Declare trees and tree variables
  art::ServiceHandle< art::TFileService > tfs;
  TTree *tDataTree;
  int run, subrun, event;
  int truth_nTracks, reco_nTracks;
  std::string endProcess, declaredEndProcess;
  std::vector<float> truth_start, truth_end, reco_start, reco_end;
  float truth_length, reco_length;
  float truth_momentum;
  float rangeTruth_pion_momentum, rangeTruth_muon_momentum;
  float range_pion_momentum, range_muon_momentum;
  float mcs_momentum_best, mcs_momentum_forward;
  float mcs_momentum_mu1, mcs_momentum_mu2, mcs_momentum_mu3, mcs_momentum_mu4;
  float mcs_momentum_pi1, mcs_momentum_pi2, mcs_momentum_pi3, mcs_momentum_pi4;

}; // End class SingleParticleMomentumStudy

#endif // END def SingleParticleMomentumStudy header
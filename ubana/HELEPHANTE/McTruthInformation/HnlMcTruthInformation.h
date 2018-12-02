#ifndef HNLMCTRUTHINFORMATION_H
#define HNLMCTRUTHINFORMATION_H

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
#include "lardataobj/RawData/RawDigit.h"

// Analyzer class
class HnlMcTruthInformation : public art::EDAnalyzer
{
public:
  explicit HnlMcTruthInformation(fhicl::ParameterSet const & pset);
  virtual ~HnlMcTruthInformation();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  std::map< int,art::Ptr<simb::MCParticle> > BuildMcpMap(
        art::ValidHandle< std::vector<simb::MCParticle> > const & mcpHandle);
  std::map< int,art::Ptr<sim::MCTrack> > BuildMctMap(
        art::ValidHandle< std::vector<sim::MCTrack> > const & mctHandle);
  float MCParticleLength(art::Ptr<simb::MCParticle> const & mcp);
  float MCTrackLength(art::Ptr<sim::MCTrack> const & mct);

private:
  // Declare trees and tree variables
  TTree *tDataTree;
  int pdgCode_mu, pdgCode_pi;
  std::string process_mu, endProcess_mu, process_pi, endProcess_pi;
  // MCParticle | Position
  float start_X_mu, start_Y_mu, start_Z_mu, start_T_mu, end_X_mu, end_Y_mu, end_Z_mu, end_T_mu;
  float start_X_pi, start_Y_pi, start_Z_pi, start_T_pi, end_X_pi, end_Y_pi, end_Z_pi, end_T_pi;
  // MCParticle | Momentum
  float start_pX_mu, start_pY_mu, start_pZ_mu, start_p_mu, start_e_mu, end_pX_mu, end_pY_mu, end_pZ_mu, end_p_mu, end_e_mu;
  float start_pX_pi, start_pY_pi, start_pZ_pi, start_p_pi, start_e_pi, end_pX_pi, end_pY_pi, end_pZ_pi, end_p_pi, end_e_pi;
  // MCParticle | Others
  float length_mu, theta_mu, phi_mu;
  float length_pi, theta_pi, phi_pi;
  // MCTrack | Position
  float track_start_X_mu, track_start_Y_mu, track_start_Z_mu, track_start_T_mu, track_end_X_mu, track_end_Y_mu, track_end_Z_mu, track_end_T_mu;
  float track_start_X_pi, track_start_Y_pi, track_start_Z_pi, track_start_T_pi, track_end_X_pi, track_end_Y_pi, track_end_Z_pi, track_end_T_pi;
  // MCTrack | Momentum
  float track_start_pX_mu, track_start_pY_mu, track_start_pZ_mu, track_start_p_mu, track_start_e_mu, track_end_pX_mu, track_end_pY_mu, track_end_pZ_mu, track_end_p_mu, track_end_e_mu;
  float track_start_pX_pi, track_start_pY_pi, track_start_pZ_pi, track_start_p_pi, track_start_e_pi, track_end_pX_pi, track_end_pY_pi, track_end_pZ_pi, track_end_p_pi, track_end_e_pi;
  // MCTrack | Others
  float track_length_mu, track_theta_mu, track_phi_mu;
  float track_length_pi, track_theta_pi, track_phi_pi;
  // Neutrino
  float X_nu, Y_nu, Z_nu, T_nu, pX_nu, pY_nu, pZ_nu, p_nu, e_nu, theta_nu, phi_nu;
  float openingAngle, invariantMass;
  // Others
  bool isContained, hasMCParticle_pi, hasMCParticle_mu, hasMCTrack_pi, hasMCTrack_mu;

  // Declare analysis variables
  int run, subrun, event;

  // Declare analysis functions
  void ClearData();
}; // End class HnlMcTruthInformation

#endif // END def HnlMcTruthInformation header
////////////////////////////////////////////////////////////////////////
// Class:       CCInclusiveTruthFilter
// Plugin Type: filter (art v2_05_01)
// File:        CCInclusiveTruthFilter_module.cc
//
// Generated at Thu Jul 26 09:54:38 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \brief CC-inclusive truth filter
 *
 * \author Adam Lister
 *
 * \email a.lister1@lancaster.ac.uk
 *
 * \notes This module filters true CC-inclusive events within a fiducial volume
 *        defined by the CC-inclusive selection defined in the UBXSec package.
 *
 *        This module is intended for use in filtering events which are 
 *        either CC-inclusive selected or true CC-inclusive inside the fiducial
 *        volume in order to reduce the data analysers need to run over.
 *
 */

// Base includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// UBXSec includes
#include "ubana/UBXSec/Algorithms/FiducialVolume.h"

// cpp includes
#include <memory>

class CCInclusiveTruthFilter;


class CCInclusiveTruthFilter : public art::EDFilter {
public:
  explicit CCInclusiveTruthFilter(fhicl::ParameterSet const & p);

  CCInclusiveTruthFilter(CCInclusiveTruthFilter const &) = delete;
  CCInclusiveTruthFilter(CCInclusiveTruthFilter &&) = delete;
  CCInclusiveTruthFilter & operator = (CCInclusiveTruthFilter const &) = delete;
  CCInclusiveTruthFilter & operator = (CCInclusiveTruthFilter &&) = delete;

  bool filter(art::Event & e) override;

private:

  // initialise services
  art::ServiceHandle< geo::Geometry > geo;

  // initialise classes
  ubana::FiducialVolume fiducialVolume;

  // variables
  int ccnc;
  int pdg;
  bool isInFv;

};


CCInclusiveTruthFilter::CCInclusiveTruthFilter(fhicl::ParameterSet const & p)
  : EDFilter{p}
// Initialize member data here.
{

  fhicl::ParameterSet const p_fv = p.get<fhicl::ParameterSet>("FiducialVolumeSettings");

  fiducialVolume.Configure(p_fv,
      geo->DetHalfHeight(),
      2.*geo->DetHalfWidth(),
      geo->DetLength());

  fiducialVolume.PrintConfig();

}

bool CCInclusiveTruthFilter::filter(art::Event & e)
{

  // get MCTruth information
  art::Handle< std::vector< simb::MCTruth > > mcTruthHandle;
  std::vector< art::Ptr< simb::MCTruth > > mcTruthVec;

  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return false;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);

  // if event has more than one MCTruth, then ignore it. In practice this 
  // throws away some events but it's a very small number for MicroBooNE
  if (mcTruthVec.size() != 1){
    std::cout << "[CCIncTruthFilter] There are " << mcTruthVec.size() << " mcTruth informations" << std::endl;
    return false;
  }

  art::Ptr< simb::MCTruth > mcTruth = mcTruthVec.at(0);

  // get mcNeutrino object
  const simb::MCNeutrino mcNeutrino = mcTruth->GetNeutrino();

  // and get the genie MCParticle associated with the neutrino
  const simb::MCParticle mcParticle_nu = mcNeutrino.Nu();


  // now check it's a true CC inclusive interaction in the fiducial volume

  ccnc = mcNeutrino.CCNC();
  pdg = mcParticle_nu.PdgCode();
  isInFv = fiducialVolume.InFV(mcParticle_nu.Vx(), mcParticle_nu.Vy(), mcParticle_nu.Vz());
 
  std::cout << "[CCIncTruthFilter] nu ccnc : " << ccnc << std::endl;
  std::cout << "[CCIncTruthFilter] nu pdg  : " << pdg << std::endl;
  std::cout << "[CCIncTruthFilter] nu inFV : " << isInFv << std::endl;

  if (ccnc == 0 && isInFv == true && std::abs(pdg) == 14){
    std::cout << "[CCIncTruthFilter] Event Passes." << std::endl;
    return true;
  }
  else{
    std::cout << "[CCIncTruthFilter] Event Fails." << std::endl;
    return false;
  }

}

DEFINE_ART_MODULE(CCInclusiveTruthFilter)

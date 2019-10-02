#ifndef __uboone_TreeInterface__
#define __uboone_TreeInterface__

/**
 * Base interface for importing ROOT trees.
 *
 * \author Z. Pavlovic
 */

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "TLorentzVector.h"

namespace uboone {

/**
 * \class TreeInterface
 * \brief Wrapper to access ROOT trees (pure virtual base class)
 */
class TreeInterface {
public:    
  /**
   * Fill in the MCFlux object from the tree data.
   *
   * \param ientry Tree entry index
   * \param mcflux MCFlux object to fill in (by reference)
   * \returns True if successful
   */
  virtual bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux) { return true; }

  /**
   * Fill in the MCTruth object from the tree data.
   *
   * \param ientry Tree entry index
   * \param mctruth MCTruth object to fill in (by reference)
   * \returns True if successful
   */
  virtual bool FillMCTruth(Long64_t ientry, simb::MCTruth& mctruth) { return true; }

  /**
   * Fill in the GTruth object from the tree data.
   *
   * \param ientry Tree entry index
   * \param gtruth GTruth object to fill in (by reference)
   * \returns True if successful
   */
  virtual bool FillGTruth(Long64_t ientry, simb::GTruth& gtruth) { return true; }

  /** Get the POT from the tree. */
  virtual const float GetPOT() = 0;

  /** Get the number of entries in the tree. */
  virtual const Long64_t GetEntries() = 0;

  /** Get the run number (for the first entry in the tree). */
  virtual const int GetRun() = 0;

  /** Get the neutrino position for the current entry. */
  virtual const TLorentzVector GetNuPosition() = 0;

  /** Get the neutrino momentum for the current entry. */
  virtual const TLorentzVector GetNuMomentum() = 0;
};

}  // namespace uboone

#endif  // __uboone_TreeInterface__


#ifndef __uboone_GSimpleInterface__
#define __uboone_GSimpleInterface__

/**
 * Interface for importing GSimple files.
 *
 * \author Z. Pavlovic
 */

#include "TreeInterface.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"

class TTree;
class TFile;

namespace uboone {

/**
 * \class GSimpleInterface
 * \brief Wrapper to access GSimple flux trees
 */
class GSimpleInterface : public TreeInterface {
public:
  /** Constructor. */
  GSimpleInterface();

  /** Destructor. */
  ~GSimpleInterface();
  
  /** Get the number of entries in the tree. */
  const Long64_t GetEntries() { return fNEntries; }

  /** Get the run number (for the first entry in the tree). */
  const int GetRun() { return fRun; }

  /** Get the POT from the tree. */
  const float GetPOT() { return fPOT; }

  /** Get the neutrino position for the current entry. */
  const TLorentzVector GetNuPosition() { return fNuPos; }

  /** Get the neutrino momentum for the current entry. */
  const TLorentzVector GetNuMomentum() { return fNuMom; }

  /**
   * Set the ROOT file to load from.
   *
   * \param rootFile The input ROOT TFile object
   */
  void SetRootFile(TFile* rootFile);

  /** 
   * Fill in the MCFlux object from the tree data.
   *
   * \param ientry Tree entry index
   * \param mcflux MCFlux object to fill in (by reference)
   * \returns True if successful
   */
  bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

private:
  TTree* fFluxTree;  //!< Tree with flux data
  TTree* fMetaTree;  //!< Tree with flux metadata
  genie::flux::GSimpleNtpEntry* fGSimpleEntry;
  genie::flux::GSimpleNtpNuMI* fGSimpleNuMI;
  genie::flux::GSimpleNtpAux* fGSimpleAux;
  genie::flux::GSimpleNtpMeta* fGSimpleMeta;
  Long64_t fNEntries;  //!< Number of events in the input tree
  int fRun;  //!< Run ID
  float fPOT;  //!< POT counter
  TLorentzVector fNuPos;  //!< Neutrino position vector
  TLorentzVector fNuMom;  //!< Neutrino momentum vector
};

}  // namespace uboone

#endif  // __uboone_GSimpleInterface__


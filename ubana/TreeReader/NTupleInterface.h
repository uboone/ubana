#ifndef __uboone_NTupleInterface__
#define __uboone_NTupleInterface__

/**
 * Interface for importing Ntuples.
 *
 * \author A. Lister, A. Mastbaum
 */

#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "TreeInterface.h"
#include "TFile.h"

class TTree;

namespace uboone {

/**
 * \class NTupleInterface
 * \brief Wrapper to access user TTrees
 */
class NTupleInterface : public TreeInterface {
public:

  /** Constructor. */
  NTupleInterface();

  /** Destructor. */
  ~NTupleInterface();

  /** Get the number of entries in the input tree. */
  const Long64_t GetEntries() { return fNEntries; }

  /** Get the current run number. */
  const int GetRun() { return fRun; }

  /** Get the POT count. */
  const float GetPOT() { return fPOT; }

  /** Get the neutrino position for the current entry. */
  const TLorentzVector GetNuPosition() { return fNuPos; }

  /** Get the neutrino momentum for the current entry. */
  const TLorentzVector GetNuMomentum() { return fNuMom; }

  /** Enable DL tree mode. */
  void SetDLMode(bool mode) { fDLMode = mode; }

  /** Enable verbose output. */
  void SetVerbose(bool verbose) { fVerbose = verbose; }

  /**
   * Set up the input ROOT file.
   *
   * The branch names are FHICL-configurable, with the mapping specified
   * in branchDef. These names are used to set the branch pointers for the
   * data required to populate the truth-level objects.
   *
   * \param rootFile The input ROOT file
   * \param treeName Path to the input TTree object
   * \param branchDef A FHICL config specifying the branch name mapping
   */
  void SetRootFile(TFile* rootFile,
                   TString treeName,
                   fhicl::ParameterSet& branchDef);

  /**
   * Fill in the MCFlux object from the tree data.
   *
   * \param ientry Tree entry index
   * \param mcflux MCFlux object to fill in (by reference)
   * \returns True if successful
   */
  bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

  /**
   * Fill in the MCTruth object from the tree data.
   *
   * \param ientry Tree entry index
   * \param mctruth MCTruth object to fill in (by reference)
   * \returns True if successful
   */
  bool FillMCTruth(Long64_t ientry, simb::MCTruth& mctruth);

  /**
   * Fill in the GTruth object from the tree data.
   *
   * \param ientry Tree entry index
   * \param gtruth GTruth object to fill in (by reference)
   * \returns True if successful
   */
  bool FillGTruth(Long64_t ientry, simb::GTruth& gtruth);

private:

  enum Max { kMaxParticles = 200};  //!< Max. MCParticles in MCTruth

  TTree* fTree;  //!< Input TTree
  int fRun;
  double fPOT;
  bool fVerbose;
  bool fDLMode;  //!< Enable DL tree mode

  // Metadata
  int run;
  int subrun;
  int event;
  Long64_t fNEntries;

  // MCFlux
  double MCFlux_NuPosX;
  double MCFlux_NuPosY;
  double MCFlux_NuPosZ;
  double MCFlux_NuMomX;
  double MCFlux_NuMomY;
  double MCFlux_NuMomZ;
  double MCFlux_NuMomE;
  int MCFlux_ntype;
  int MCFlux_ptype;
  double MCFlux_nimpwt;
  double MCFlux_dk2gen;
  double MCFlux_nenergyn;
  double MCFlux_tpx;
  double MCFlux_tpy;
  double MCFlux_tpz;
  int MCFlux_tptype;
  double MCFlux_vx;
  double MCFlux_vy;
  double MCFlux_vz;

  TLorentzVector fNuPos;  //!< Neutrino vertex vector
  TLorentzVector fNuMom;  //!< Neutrino momentum vector

  // MCTruth
  int MCTruth_neutrino_CCNC;
  int MCTruth_neutrino_mode;
  int MCTruth_neutrino_interactionType;
  int MCTruth_neutrino_target;
  int MCTruth_neutrino_nucleon;
  int MCTruth_neutrino_quark;
  double MCTruth_neutrino_W;
  double MCTruth_neutrino_X;
  double MCTruth_neutrino_Y;
  double MCTruth_neutrino_QSqr;

  int MCTruth_NParticles;
  int MCTruth_particles_TrackId[kMaxParticles];
  int MCTruth_particles_PdgCode[kMaxParticles];
  int MCTruth_particles_Mother[kMaxParticles];
  int MCTruth_particles_StatusCode[kMaxParticles];
  int MCTruth_particles_NumberDaughters[kMaxParticles];
  int MCTruth_particles_Daughters[kMaxParticles][100];
  double MCTruth_particles_Gvx[kMaxParticles];
  double MCTruth_particles_Gvy[kMaxParticles];
  double MCTruth_particles_Gvz[kMaxParticles];
  double MCTruth_particles_Gvt[kMaxParticles];
  double MCTruth_particles_px0[kMaxParticles];
  double MCTruth_particles_py0[kMaxParticles];
  double MCTruth_particles_pz0[kMaxParticles];
  double MCTruth_particles_e0[kMaxParticles];
  int MCTruth_particles_Rescatter[kMaxParticles];
  double MCTruth_particles_polx[kMaxParticles];
  double MCTruth_particles_poly[kMaxParticles];
  double MCTruth_particles_polz[kMaxParticles];

  // GTruth
  int GTruth_ProbePDG;
  bool GTruth_IsSeaQuark;
  int GTruth_IsSeaQuark_int;
  int GTruth_tgtPDG;
  int GTruth_tgtZ;
  int GTruth_tgtA;
  double GTruth_tgtP4x;
  double GTruth_tgtP4y;
  double GTruth_tgtP4z;
  double GTruth_tgtP4E;
  double GTruth_weight;
  double GTruth_probability;
  double GTruth_Xsec;
  double GTruth_DiffXsec;
  int    GTruth_GPhaseSpace;
  double GTruth_vertexX;
  double GTruth_vertexY;
  double GTruth_vertexZ;
  double GTruth_vertexT;
  int GTruth_Gscatter;
  int GTruth_Gint;
  int GTruth_ResNum;
  int GTruth_NumPiPlus;
  int GTruth_NumPi0;
  int GTruth_NumPiMinus;
  int GTruth_NumProton;
  int GTruth_NumNeutron;
  bool GTruth_IsCharm;
  int GTruth_IsCharm_int;
  bool GTruth_IsStrange;
  int GTruth_IsStrange_int;
  int GTruth_StrangeHadronPDG;
  int GTruth_CharmHadronPDG;
  int GTruth_DecayMode;
  double GTruth_gX;
  double GTruth_gY;
  //double GTruth_gZ; //not used, is this an mistake?
  double GTruth_gT;
  double GTruth_gW;
  double GTruth_gQ2;
  double GTruth_gq2;
  double GTruth_ProbeP4x;
  double GTruth_ProbeP4y;
  double GTruth_ProbeP4z;
  double GTruth_ProbeP4E;
  double GTruth_HitNucP4x;
  double GTruth_HitNucP4y;
  double GTruth_HitNucP4z;
  double GTruth_HitNucP4E;
  double GTruth_HitNucPos;
  double GTruth_FShadSystP4x;
  double GTruth_FShadSystP4y;
  double GTruth_FShadSystP4z;
  double GTruth_FShadSystP4E;
};

}  // namespace uboone

#endif // __uboone_NTupleInterface__


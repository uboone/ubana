Tree Reader
===========
This module reads in data from a (semi-)arbitrary ROOT tree, to populate
an `art::Event`.

The primary goal is to load in truth information from post-selection analysis
trees and calculate weights with EventWeight, without need to go back to the
the original art ROOT files.

The initial implementation of this module reads in MCTruth, MCFlux, and GTruth
objects, which are required for flux and cross section reweighting. The output
MCEventWeights can then be attached to the analysis trees, for propagation.

NTuple Input
------------
In Ntuple input mode, the TreeReader source accepts a ROOT file, where the
user specifies the name of the ROOT tree and a list of branch names where the
required data is stored.

An example configuration, where the names of the branches in the input TTree
happen to match the key names, which are fixed.

    source: {
      fileNames: []
      module_type: TreeReader
      skipEvents: 0
      maxEvents: -1
      inputType: "ntuple"
      treeName: "NAME OF TREE"
      branches: {
        MCFlux_NuPosX: "MCFlux_NuPosX"
        MCFlux_NuPosY: "flux_nposx"
        MCFlux_NuPosZ: "some_arbitrary_name"
        MCFlux_NuMomX: "etc"
        # ...
      }
    }

Note that the tree can be inside a `TDirectory`, e.g. `dir/treename`. The
full list of fields follows here:

    #############################################################################
    # KEY                            SOURCE                                TYPE #
    #############################################################################

    # Metadata.
    #
    # Information to identify an event.
    #

    run                              Run ID                                int
    subrun                           Subrun ID                             int
    evt                              Event Number                          int

    # MCFlux.
    #
    # A copy of fields from MCFlux required for reweighting. This includes a
    # few required from MCTruth, which are required in the case of a GSimple
    # file input to "fake" the MCTruth object. For NTuple input, the MCTruth
    # is filled in usng the fields below.
    #
    MCFlux_NuPosX                    MCTruth::fMCNeutrino::fNu.Vx(0)       double
    MCFlux_NuPosY                    MCTruth::fMCNeutrino::fNu.Vy(0)       double
    MCFlux_NuPosZ                    MCTruth::fMCNeutrino::fNu.Vz(0)       double
    MCFlux_NuMomX                    MCTruth::fMCNeutrino::fNu.Px(0)       double
    MCFlux_NuMomY                    MCTruth::fMCNeutrino::fNu.Py(0)       double
    MCFlux_NuMomZ                    MCTruth::fMCNeutrino::fNu.Pz(0)       double
    MCFlux_NuMomE                    MCTruth::fMCNeutrino::fNu.E(0)        double
    MCFlux_ntype                     MCFlux::ntype                         int
    MCFlux_ptype                     MCFlux::ptype                         int
    MCFlux_nimpwt                    MCFlux::nimpwt                        double
    MCFlux_dk2gen                    MCFlux::dk2gen                        double
    MCFlux_nenergyn                  MCFlux::nenergyn                      double
    MCFlux_tpx                       MCFlux::tpx                           double
    MCFlux_tpy                       MCFlux::tpy                           double
    MCFlux_tpz                       MCFlux::tpz                           double
    MCFlux_tptype                    MCFlux::tptype                        int
    MCFlux_vx                        MCFlux::vx                            double
    MCFlux_vy                        MCFlux::vy                            double
    MCFlux_vz                        MCFlux::vz                            double

    # MCTruth.
    #
    # A copy of simb::MCTruth information. Note that this contains a list of
    # MCParticles, which is flattened into a set of arrays. One two-dimensional
    # array, Daughters, contains a list of the track IDs of the daughters of
    # ith MCParticle.
    #
    MCTruth_NParticles                    MCTruth::NParticles                   int
    MCTruth_particles_TrackId             MCTruth::fPartList[i].TrackId         int[]
    MCTruth_particles_PdgCode             MCTruth::fPartList[i].PdgCode         int[]
    MCTruth_particles_Mother              MCTruth::fPartList[i].Mother          int[]
    MCTruth_particles_StatusCode          MCTruth::fPartList[i].StatusCode      int[]
    MCTruth_particles_NumberDaughters     MCTruth::fPartList[i].NumberDaughters int[]
    MCTruth_particles_Daughters           MCTruth::fPartList[i].Daughters       int[][]
    MCTruth_particles_Gvx                 MCTruth::fPartList[i].Gvx             double[]
    MCTruth_particles_Gvy                 MCTruth::fPartList[i].Gvy             double[]
    MCTruth_particles_Gvz                 MCTruth::fPartList[i].Gvz             double[]
    MCTruth_particles_Gvt                 MCTruth::fPartList[i].Gvt             double[]
    MCTruth_particles_px0                 MCTruth::fPartList[i].Px(0)           double[]
    MCTruth_particles_py0                 MCTruth::fPartList[i].Py(0)           double[]
    MCTruth_particles_pz0                 MCTruth::fPartList[i].Pz(0)           double[]
    MCTruth_particles_e0                  MCTruth::fPartList[i].E(0)            double[]
    MCTruth_particles_Rescatter           MCTruth::fPartList[i].Rescatter       int[]
    MCTruth_particles_polx                MCTruth::fPartList[i].Polarization.fX double[]
    MCTruth_particles_poly                MCTruth::fPartList[i].Polarization.fY double[]
    MCTruth_particles_polz                MCTruth::fPartList[i].Polarization.fZ double[]
    MCTruth_neutrino_CCNC                 MCTruth::fMCNeutrino.CCNC             int
    MCTruth_neutrino_mode                 MCTruth::fMCNeutrino.Mode             int
    MCTruth_neutrino_interactionType      MCTruth::fMCNeutrino.InteractionType  int
    MCTruth_neutrino_target               MCTruth::fMCNeutrino.Target           int
    MCTruth_neutrino_nucleon              MCTruth::fMCNeutrino.HitNuc           int
    MCTruth_neutrino_quark                MCTruth::fMCNeutrino.HitQuark         int
    MCTruth_neutrino_W                    MCTruth::fMCNeutrino.W                double
    MCTruth_neutrino_X                    MCTruth::fMCNeutrino.X                double
    MCTruth_neutrino_Y                    MCTruth::fMCNeutrino.Y                double
    MCTruth_neutrino_QSqr                 MCTruth::fMCNeutrino.QSqr             double

    # GTruth.
    #
    # A copy of all fields in the simb::GTruth object.
    #
    GTruth_IsSeaQuark                GTruth::fIsSeaQuark                   bool
    GTruth_tgtPDG                    GTruth::ftgtPDG                       int
    GTruth_weight                    GTruth::fweight                       double
    GTruth_probability               GTruth::fprobability                  double
    GTruth_Xsec                      GTruth::fXsec                         double
    GTruth_DiffXsec                  GTruth::fDiffXsec                    double
    GTruth_vertexX                   GTruth::fVertex.X                     double
    GTruth_vertexY                   GTruth::fVertex.Y                     double
    GTruth_vertexZ                   GTruth::fVertex.Z                     double
    GTruth_vertexT                   GTruth::fVertex.T                     double
    GTruth_Gscatter                  GTruth::fGscatter                     int
    GTruth_Gint                      GTruth::fGint                         int
    GTruth_ResNum                    GTruth::fResNum                       int
    GTruth_NumPiPlus                 GTruth::fNumPiPlus                    int
    GTruth_NumPi0                    GTruth::fNumPi0                       int
    GTruth_NumPiMinus                GTruth::fNumPiMinus                   int
    GTruth_NumProton                 GTruth::fNumProton                    int
    GTruth_NumNeutron                GTruth::fNumNeutron                   int
    GTruth_IsCharm                   GTruth::fIsCharm                      bool
    GTruth_gX                        GTruth::fgX                           double
    GTruth_gY                        GTruth::fgY                           double
    GTruth_gZ                        GTruth::fgZ                           double
    GTruth_gT                        GTruth::fgT                           double
    GTruth_gW                        GTruth::fgW                           double
    GTruth_gQ2                       GTruth::fgQ2                          double
    GTruth_gq2                       GTruth::fgq2                          double
    GTruth_ProbePDG                  GTruth::fProbePDG                     int
    GTruth_ProbeP4x                  GTruth::fProbeP4.x                    double
    GTruth_ProbeP4y                  GTruth::fProbeP4.y                    double
    GTruth_ProbeP4z                  GTruth::fProbeP4.z                    double
    GTruth_ProbeP4E                  GTruth::fProbeP4.E                    double
    GTruth_HitNucP4x                 GTruth::fHitNucP4.x                   double
    GTruth_HitNucP4y                 GTruth::fHitNucP4.y                   double
    GTruth_HitNucP4z                 GTruth::fHitNucP4.z                   double
    GTruth_HitNucP4E                 GTruth::fHitNucP4.E                   double
    GTruth_FShadSystP4x              GTruth::fFShadSystP4.x                double
    GTruth_FShadSystP4y              GTruth::fFShadSystP4.y                double
    GTruth_FShadSystP4z              GTruth::fFShadSystP4.z                double
    GTruth_FShadSystP4E              GTruth::fFShadSystP4.E                double

Authors
-------
This module is based on the FluxReader by Zarko Pavlovic. It was generalized
to handle other analysis tree inputs by Adam Lister and Andy Mastbaum.


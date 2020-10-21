//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 21 07:41:37 2020 by ROOT version 6.12/06
// from TTree OutputTree/Truth Info Tree
// found on file: analysisOutput.root
//////////////////////////////////////////////////////////

#ifndef OutputTree_h
#define OutputTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "TVector3.h"
#include "vector"
//#include "ubana/HyperonProduction/util/SimParticle.h"
#include "vector"
//#include "ubana/HyperonProduction/util/RecoParticle.h"


//set input file here
std::string infilename = "analysisOutput.root";


class OutputTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxLepton = 1;
   static constexpr Int_t kMaxHyperon = 1;
   static constexpr Int_t kMaxPrimaryNucleon = 14;
   static constexpr Int_t kMaxPrimaryPion = 5;
   static constexpr Int_t kMaxDecay = 2;
   static constexpr Int_t kMaxSigmaZeroDecayPhoton = 1;
   static constexpr Int_t kMaxSigmaZeroDecayLambda = 1;
   static constexpr Int_t kMaxTracklikePrimaryDaughters = 7;
   static constexpr Int_t kMaxShowerlikePrimaryDaughters = 3;

   // Declaration of leaf types
   UInt_t          EventID;
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           fileID;
   string          *Mode;
   string          *CCNC;
   Bool_t          InActiveTPC;
   Bool_t          IsHyperon;
   Bool_t          IsSigmaZero;
   Bool_t          IsLambda;
   Bool_t          IsLambdaCharged;
   Bool_t          IsSignal;
   Bool_t          SelectedEvent;
   Double_t        NuEnergy;
   Int_t           NuPDG;
   TVector3        *TruePrimaryVertex;
   Int_t           Lepton_;
   Int_t           Lepton_PDG[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_E[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_Px[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_Py[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_Pz[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_ModMomentum[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_KE[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_StartX[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_StartY[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_StartZ[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_EndX[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_EndY[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_EndZ[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_Travel[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_Theta[kMaxLepton];   //[Lepton_]
   Double_t        Lepton_Phi[kMaxLepton];   //[Lepton_]
   Int_t           Lepton_Origin[kMaxLepton];   //[Lepton_]
   Int_t           Hyperon_;
   Int_t           Hyperon_PDG[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_E[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_Px[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_Py[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_Pz[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_ModMomentum[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_KE[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_StartX[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_StartY[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_StartZ[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_EndX[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_EndY[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_EndZ[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_Travel[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_Theta[kMaxHyperon];   //[Hyperon_]
   Double_t        Hyperon_Phi[kMaxHyperon];   //[Hyperon_]
   Int_t           Hyperon_Origin[kMaxHyperon];   //[Hyperon_]
   Int_t           PrimaryNucleon_;
   Int_t           PrimaryNucleon_PDG[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_E[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_Px[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_Py[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_Pz[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_ModMomentum[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_KE[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_StartX[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_StartY[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_StartZ[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_EndX[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_EndY[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_EndZ[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_Travel[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_Theta[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Double_t        PrimaryNucleon_Phi[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Int_t           PrimaryNucleon_Origin[kMaxPrimaryNucleon];   //[PrimaryNucleon_]
   Int_t           PrimaryPion_;
   Int_t           PrimaryPion_PDG[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_E[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_Px[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_Py[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_Pz[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_ModMomentum[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_KE[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_StartX[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_StartY[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_StartZ[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_EndX[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_EndY[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_EndZ[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_Travel[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_Theta[kMaxPrimaryPion];   //[PrimaryPion_]
   Double_t        PrimaryPion_Phi[kMaxPrimaryPion];   //[PrimaryPion_]
   Int_t           PrimaryPion_Origin[kMaxPrimaryPion];   //[PrimaryPion_]
   TVector3        *DecayVertex;
   Int_t           Decay_;
   Int_t           Decay_PDG[kMaxDecay];   //[Decay_]
   Double_t        Decay_E[kMaxDecay];   //[Decay_]
   Double_t        Decay_Px[kMaxDecay];   //[Decay_]
   Double_t        Decay_Py[kMaxDecay];   //[Decay_]
   Double_t        Decay_Pz[kMaxDecay];   //[Decay_]
   Double_t        Decay_ModMomentum[kMaxDecay];   //[Decay_]
   Double_t        Decay_KE[kMaxDecay];   //[Decay_]
   Double_t        Decay_StartX[kMaxDecay];   //[Decay_]
   Double_t        Decay_StartY[kMaxDecay];   //[Decay_]
   Double_t        Decay_StartZ[kMaxDecay];   //[Decay_]
   Double_t        Decay_EndX[kMaxDecay];   //[Decay_]
   Double_t        Decay_EndY[kMaxDecay];   //[Decay_]
   Double_t        Decay_EndZ[kMaxDecay];   //[Decay_]
   Double_t        Decay_Travel[kMaxDecay];   //[Decay_]
   Double_t        Decay_Theta[kMaxDecay];   //[Decay_]
   Double_t        Decay_Phi[kMaxDecay];   //[Decay_]
   Int_t           Decay_Origin[kMaxDecay];   //[Decay_]
   Int_t           SigmaZeroDecayPhoton_;
   Int_t           SigmaZeroDecayPhoton_PDG[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_E[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_Px[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_Py[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_Pz[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_ModMomentum[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_KE[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_StartX[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_StartY[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_StartZ[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_EndX[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_EndY[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_EndZ[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_Travel[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_Theta[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Double_t        SigmaZeroDecayPhoton_Phi[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Int_t           SigmaZeroDecayPhoton_Origin[kMaxSigmaZeroDecayPhoton];   //[SigmaZeroDecayPhoton_]
   Int_t           SigmaZeroDecayLambda_;
   Int_t           SigmaZeroDecayLambda_PDG[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_E[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_Px[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_Py[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_Pz[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_ModMomentum[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_KE[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_StartX[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_StartY[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_StartZ[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_EndX[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_EndY[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_EndZ[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_Travel[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_Theta[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        SigmaZeroDecayLambda_Phi[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Int_t           SigmaZeroDecayLambda_Origin[kMaxSigmaZeroDecayLambda];   //[SigmaZeroDecayLambda_]
   Double_t        DecayOpeningAngle;
   Double_t        LeptonPionAngle;
   Double_t        LeptonNucleonAngle;
   TVector3        *RecoPrimaryVertex;
   Int_t           NPrimaryTrackDaughters;
   Int_t           NPrimaryShowerDaughters;
   Int_t           TracklikePrimaryDaughters_;
   Int_t           TracklikePrimaryDaughters_PDG[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackShowerScore[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_X[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_Y[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_Z[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_Displacement[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackLength[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackPID[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackProtonChi2[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackMuonChi2[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackKaonChi2[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackPionChi2[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Bool_t          TracklikePrimaryDaughters_HasTruth[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Int_t           TracklikePrimaryDaughters_TrackTruePDG[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTrueE[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTruePx[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTruePy[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTruePz[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTrueModMomentum[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTrueKE[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTrueLength[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Int_t           TracklikePrimaryDaughters_TrackTrueOrigin[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackTrueTotalEdep[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Double_t        TracklikePrimaryDaughters_TrackEdepPurity[kMaxTracklikePrimaryDaughters];   //[TracklikePrimaryDaughters_]
   Int_t           ShowerlikePrimaryDaughters_;
   Int_t           ShowerlikePrimaryDaughters_PDG[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackShowerScore[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_X[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_Y[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_Z[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_Displacement[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackLength[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackPID[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackProtonChi2[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackMuonChi2[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackKaonChi2[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackPionChi2[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Bool_t          ShowerlikePrimaryDaughters_HasTruth[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Int_t           ShowerlikePrimaryDaughters_TrackTruePDG[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTrueE[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTruePx[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTruePy[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTruePz[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTrueModMomentum[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTrueKE[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTrueLength[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Int_t           ShowerlikePrimaryDaughters_TrackTrueOrigin[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackTrueTotalEdep[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Double_t        ShowerlikePrimaryDaughters_TrackEdepPurity[kMaxShowerlikePrimaryDaughters];   //[ShowerlikePrimaryDaughters_]
   Int_t           MuonIndex;

   // List of branches
   TBranch        *b_EventID;   //!
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_fileID;   //!
   TBranch        *b_Mode;   //!
   TBranch        *b_CCNC;   //!
   TBranch        *b_InActiveTPC;   //!
   TBranch        *b_IsHyperon;   //!
   TBranch        *b_IsSigmaZero;   //!
   TBranch        *b_IsLambda;   //!
   TBranch        *b_IsLambdaCharged;   //!
   TBranch        *b_IsSignal;   //!
   TBranch        *b_SelectedEvent;   //!
   TBranch        *b_NuEnergy;   //!
   TBranch        *b_NuPDG;   //!
   TBranch        *b_TruePrimaryVertex;   //!
   TBranch        *b_Lepton_;   //!
   TBranch        *b_Lepton_PDG;   //!
   TBranch        *b_Lepton_E;   //!
   TBranch        *b_Lepton_Px;   //!
   TBranch        *b_Lepton_Py;   //!
   TBranch        *b_Lepton_Pz;   //!
   TBranch        *b_Lepton_ModMomentum;   //!
   TBranch        *b_Lepton_KE;   //!
   TBranch        *b_Lepton_StartX;   //!
   TBranch        *b_Lepton_StartY;   //!
   TBranch        *b_Lepton_StartZ;   //!
   TBranch        *b_Lepton_EndX;   //!
   TBranch        *b_Lepton_EndY;   //!
   TBranch        *b_Lepton_EndZ;   //!
   TBranch        *b_Lepton_Travel;   //!
   TBranch        *b_Lepton_Theta;   //!
   TBranch        *b_Lepton_Phi;   //!
   TBranch        *b_Lepton_Origin;   //!
   TBranch        *b_Hyperon_;   //!
   TBranch        *b_Hyperon_PDG;   //!
   TBranch        *b_Hyperon_E;   //!
   TBranch        *b_Hyperon_Px;   //!
   TBranch        *b_Hyperon_Py;   //!
   TBranch        *b_Hyperon_Pz;   //!
   TBranch        *b_Hyperon_ModMomentum;   //!
   TBranch        *b_Hyperon_KE;   //!
   TBranch        *b_Hyperon_StartX;   //!
   TBranch        *b_Hyperon_StartY;   //!
   TBranch        *b_Hyperon_StartZ;   //!
   TBranch        *b_Hyperon_EndX;   //!
   TBranch        *b_Hyperon_EndY;   //!
   TBranch        *b_Hyperon_EndZ;   //!
   TBranch        *b_Hyperon_Travel;   //!
   TBranch        *b_Hyperon_Theta;   //!
   TBranch        *b_Hyperon_Phi;   //!
   TBranch        *b_Hyperon_Origin;   //!
   TBranch        *b_PrimaryNucleon_;   //!
   TBranch        *b_PrimaryNucleon_PDG;   //!
   TBranch        *b_PrimaryNucleon_E;   //!
   TBranch        *b_PrimaryNucleon_Px;   //!
   TBranch        *b_PrimaryNucleon_Py;   //!
   TBranch        *b_PrimaryNucleon_Pz;   //!
   TBranch        *b_PrimaryNucleon_ModMomentum;   //!
   TBranch        *b_PrimaryNucleon_KE;   //!
   TBranch        *b_PrimaryNucleon_StartX;   //!
   TBranch        *b_PrimaryNucleon_StartY;   //!
   TBranch        *b_PrimaryNucleon_StartZ;   //!
   TBranch        *b_PrimaryNucleon_EndX;   //!
   TBranch        *b_PrimaryNucleon_EndY;   //!
   TBranch        *b_PrimaryNucleon_EndZ;   //!
   TBranch        *b_PrimaryNucleon_Travel;   //!
   TBranch        *b_PrimaryNucleon_Theta;   //!
   TBranch        *b_PrimaryNucleon_Phi;   //!
   TBranch        *b_PrimaryNucleon_Origin;   //!
   TBranch        *b_PrimaryPion_;   //!
   TBranch        *b_PrimaryPion_PDG;   //!
   TBranch        *b_PrimaryPion_E;   //!
   TBranch        *b_PrimaryPion_Px;   //!
   TBranch        *b_PrimaryPion_Py;   //!
   TBranch        *b_PrimaryPion_Pz;   //!
   TBranch        *b_PrimaryPion_ModMomentum;   //!
   TBranch        *b_PrimaryPion_KE;   //!
   TBranch        *b_PrimaryPion_StartX;   //!
   TBranch        *b_PrimaryPion_StartY;   //!
   TBranch        *b_PrimaryPion_StartZ;   //!
   TBranch        *b_PrimaryPion_EndX;   //!
   TBranch        *b_PrimaryPion_EndY;   //!
   TBranch        *b_PrimaryPion_EndZ;   //!
   TBranch        *b_PrimaryPion_Travel;   //!
   TBranch        *b_PrimaryPion_Theta;   //!
   TBranch        *b_PrimaryPion_Phi;   //!
   TBranch        *b_PrimaryPion_Origin;   //!
   TBranch        *b_DecayVertex;   //!
   TBranch        *b_Decay_;   //!
   TBranch        *b_Decay_PDG;   //!
   TBranch        *b_Decay_E;   //!
   TBranch        *b_Decay_Px;   //!
   TBranch        *b_Decay_Py;   //!
   TBranch        *b_Decay_Pz;   //!
   TBranch        *b_Decay_ModMomentum;   //!
   TBranch        *b_Decay_KE;   //!
   TBranch        *b_Decay_StartX;   //!
   TBranch        *b_Decay_StartY;   //!
   TBranch        *b_Decay_StartZ;   //!
   TBranch        *b_Decay_EndX;   //!
   TBranch        *b_Decay_EndY;   //!
   TBranch        *b_Decay_EndZ;   //!
   TBranch        *b_Decay_Travel;   //!
   TBranch        *b_Decay_Theta;   //!
   TBranch        *b_Decay_Phi;   //!
   TBranch        *b_Decay_Origin;   //!
   TBranch        *b_SigmaZeroDecayPhoton_;   //!
   TBranch        *b_SigmaZeroDecayPhoton_PDG;   //!
   TBranch        *b_SigmaZeroDecayPhoton_E;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Px;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Py;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Pz;   //!
   TBranch        *b_SigmaZeroDecayPhoton_ModMomentum;   //!
   TBranch        *b_SigmaZeroDecayPhoton_KE;   //!
   TBranch        *b_SigmaZeroDecayPhoton_StartX;   //!
   TBranch        *b_SigmaZeroDecayPhoton_StartY;   //!
   TBranch        *b_SigmaZeroDecayPhoton_StartZ;   //!
   TBranch        *b_SigmaZeroDecayPhoton_EndX;   //!
   TBranch        *b_SigmaZeroDecayPhoton_EndY;   //!
   TBranch        *b_SigmaZeroDecayPhoton_EndZ;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Travel;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Theta;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Phi;   //!
   TBranch        *b_SigmaZeroDecayPhoton_Origin;   //!
   TBranch        *b_SigmaZeroDecayLambda_;   //!
   TBranch        *b_SigmaZeroDecayLambda_PDG;   //!
   TBranch        *b_SigmaZeroDecayLambda_E;   //!
   TBranch        *b_SigmaZeroDecayLambda_Px;   //!
   TBranch        *b_SigmaZeroDecayLambda_Py;   //!
   TBranch        *b_SigmaZeroDecayLambda_Pz;   //!
   TBranch        *b_SigmaZeroDecayLambda_ModMomentum;   //!
   TBranch        *b_SigmaZeroDecayLambda_KE;   //!
   TBranch        *b_SigmaZeroDecayLambda_StartX;   //!
   TBranch        *b_SigmaZeroDecayLambda_StartY;   //!
   TBranch        *b_SigmaZeroDecayLambda_StartZ;   //!
   TBranch        *b_SigmaZeroDecayLambda_EndX;   //!
   TBranch        *b_SigmaZeroDecayLambda_EndY;   //!
   TBranch        *b_SigmaZeroDecayLambda_EndZ;   //!
   TBranch        *b_SigmaZeroDecayLambda_Travel;   //!
   TBranch        *b_SigmaZeroDecayLambda_Theta;   //!
   TBranch        *b_SigmaZeroDecayLambda_Phi;   //!
   TBranch        *b_SigmaZeroDecayLambda_Origin;   //!
   TBranch        *b_DecayOpeningAngle;   //!
   TBranch        *b_LeptonPionAngle;   //!
   TBranch        *b_LeptonNucleonAngle;   //!
   TBranch        *b_RecoPrimaryVertex;   //!
   TBranch        *b_NPrimaryTrackDaughters;   //!
   TBranch        *b_NPrimaryShowerDaughters;   //!
   TBranch        *b_TracklikePrimaryDaughters_;   //!
   TBranch        *b_TracklikePrimaryDaughters_PDG;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackShowerScore;   //!
   TBranch        *b_TracklikePrimaryDaughters_X;   //!
   TBranch        *b_TracklikePrimaryDaughters_Y;   //!
   TBranch        *b_TracklikePrimaryDaughters_Z;   //!
   TBranch        *b_TracklikePrimaryDaughters_Displacement;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackLength;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackPID;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackProtonChi2;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackMuonChi2;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackKaonChi2;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackPionChi2;   //!
   TBranch        *b_TracklikePrimaryDaughters_HasTruth;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTruePDG;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTrueE;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTruePx;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTruePy;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTruePz;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTrueModMomentum;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTrueKE;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTrueLength;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTrueOrigin;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackTrueTotalEdep;   //!
   TBranch        *b_TracklikePrimaryDaughters_TrackEdepPurity;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_PDG;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackShowerScore;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_X;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_Y;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_Z;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_Displacement;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackLength;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackPID;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackProtonChi2;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackMuonChi2;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackKaonChi2;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackPionChi2;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_HasTruth;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTruePDG;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTrueE;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTruePx;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTruePy;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTruePz;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTrueModMomentum;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTrueKE;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTrueLength;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTrueOrigin;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackTrueTotalEdep;   //!
   TBranch        *b_ShowerlikePrimaryDaughters_TrackEdepPurity;   //!
   TBranch        *b_MuonIndex;   //!

   OutputTree(TTree *tree=0);
   virtual ~OutputTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef OutputTree_cxx
OutputTree::OutputTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infilename.c_str());
      if (!f || !f->IsOpen()) {
         f = new TFile(infilename.c_str());
      }
      std::string dirname = infilename + ":/ana";
      TDirectory * dir = (TDirectory*)f->Get(dirname.c_str());
      dir->GetObject("OutputTree",tree);

   }
   Init(tree);
}

OutputTree::~OutputTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t OutputTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t OutputTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void OutputTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Mode = 0;
   CCNC = 0;
   TruePrimaryVertex = 0;
   DecayVertex = 0;
   RecoPrimaryVertex = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("fileID", &fileID, &b_fileID);
   fChain->SetBranchAddress("Mode", &Mode, &b_Mode);
   fChain->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
   fChain->SetBranchAddress("InActiveTPC", &InActiveTPC, &b_InActiveTPC);
   fChain->SetBranchAddress("IsHyperon", &IsHyperon, &b_IsHyperon);
   fChain->SetBranchAddress("IsSigmaZero", &IsSigmaZero, &b_IsSigmaZero);
   fChain->SetBranchAddress("IsLambda", &IsLambda, &b_IsLambda);
   fChain->SetBranchAddress("IsLambdaCharged", &IsLambdaCharged, &b_IsLambdaCharged);
   fChain->SetBranchAddress("IsSignal", &IsSignal, &b_IsSignal);
   fChain->SetBranchAddress("SelectedEvent", &SelectedEvent, &b_SelectedEvent);
   fChain->SetBranchAddress("NuEnergy", &NuEnergy, &b_NuEnergy);
   fChain->SetBranchAddress("NuPDG", &NuPDG, &b_NuPDG);
   fChain->SetBranchAddress("TruePrimaryVertex", &TruePrimaryVertex, &b_TruePrimaryVertex);
   fChain->SetBranchAddress("Lepton", &Lepton_, &b_Lepton_);
   fChain->SetBranchAddress("Lepton.PDG", Lepton_PDG, &b_Lepton_PDG);
   fChain->SetBranchAddress("Lepton.E", Lepton_E, &b_Lepton_E);
   fChain->SetBranchAddress("Lepton.Px", Lepton_Px, &b_Lepton_Px);
   fChain->SetBranchAddress("Lepton.Py", Lepton_Py, &b_Lepton_Py);
   fChain->SetBranchAddress("Lepton.Pz", Lepton_Pz, &b_Lepton_Pz);
   fChain->SetBranchAddress("Lepton.ModMomentum", Lepton_ModMomentum, &b_Lepton_ModMomentum);
   fChain->SetBranchAddress("Lepton.KE", Lepton_KE, &b_Lepton_KE);
   fChain->SetBranchAddress("Lepton.StartX", Lepton_StartX, &b_Lepton_StartX);
   fChain->SetBranchAddress("Lepton.StartY", Lepton_StartY, &b_Lepton_StartY);
   fChain->SetBranchAddress("Lepton.StartZ", Lepton_StartZ, &b_Lepton_StartZ);
   fChain->SetBranchAddress("Lepton.EndX", Lepton_EndX, &b_Lepton_EndX);
   fChain->SetBranchAddress("Lepton.EndY", Lepton_EndY, &b_Lepton_EndY);
   fChain->SetBranchAddress("Lepton.EndZ", Lepton_EndZ, &b_Lepton_EndZ);
   fChain->SetBranchAddress("Lepton.Travel", Lepton_Travel, &b_Lepton_Travel);
   fChain->SetBranchAddress("Lepton.Theta", Lepton_Theta, &b_Lepton_Theta);
   fChain->SetBranchAddress("Lepton.Phi", Lepton_Phi, &b_Lepton_Phi);
   fChain->SetBranchAddress("Lepton.Origin", Lepton_Origin, &b_Lepton_Origin);
   fChain->SetBranchAddress("Hyperon", &Hyperon_, &b_Hyperon_);
   fChain->SetBranchAddress("Hyperon.PDG", Hyperon_PDG, &b_Hyperon_PDG);
   fChain->SetBranchAddress("Hyperon.E", Hyperon_E, &b_Hyperon_E);
   fChain->SetBranchAddress("Hyperon.Px", Hyperon_Px, &b_Hyperon_Px);
   fChain->SetBranchAddress("Hyperon.Py", Hyperon_Py, &b_Hyperon_Py);
   fChain->SetBranchAddress("Hyperon.Pz", Hyperon_Pz, &b_Hyperon_Pz);
   fChain->SetBranchAddress("Hyperon.ModMomentum", Hyperon_ModMomentum, &b_Hyperon_ModMomentum);
   fChain->SetBranchAddress("Hyperon.KE", Hyperon_KE, &b_Hyperon_KE);
   fChain->SetBranchAddress("Hyperon.StartX", Hyperon_StartX, &b_Hyperon_StartX);
   fChain->SetBranchAddress("Hyperon.StartY", Hyperon_StartY, &b_Hyperon_StartY);
   fChain->SetBranchAddress("Hyperon.StartZ", Hyperon_StartZ, &b_Hyperon_StartZ);
   fChain->SetBranchAddress("Hyperon.EndX", Hyperon_EndX, &b_Hyperon_EndX);
   fChain->SetBranchAddress("Hyperon.EndY", Hyperon_EndY, &b_Hyperon_EndY);
   fChain->SetBranchAddress("Hyperon.EndZ", Hyperon_EndZ, &b_Hyperon_EndZ);
   fChain->SetBranchAddress("Hyperon.Travel", Hyperon_Travel, &b_Hyperon_Travel);
   fChain->SetBranchAddress("Hyperon.Theta", Hyperon_Theta, &b_Hyperon_Theta);
   fChain->SetBranchAddress("Hyperon.Phi", Hyperon_Phi, &b_Hyperon_Phi);
   fChain->SetBranchAddress("Hyperon.Origin", Hyperon_Origin, &b_Hyperon_Origin);
   fChain->SetBranchAddress("PrimaryNucleon", &PrimaryNucleon_, &b_PrimaryNucleon_);
   fChain->SetBranchAddress("PrimaryNucleon.PDG", PrimaryNucleon_PDG, &b_PrimaryNucleon_PDG);
   fChain->SetBranchAddress("PrimaryNucleon.E", PrimaryNucleon_E, &b_PrimaryNucleon_E);
   fChain->SetBranchAddress("PrimaryNucleon.Px", PrimaryNucleon_Px, &b_PrimaryNucleon_Px);
   fChain->SetBranchAddress("PrimaryNucleon.Py", PrimaryNucleon_Py, &b_PrimaryNucleon_Py);
   fChain->SetBranchAddress("PrimaryNucleon.Pz", PrimaryNucleon_Pz, &b_PrimaryNucleon_Pz);
   fChain->SetBranchAddress("PrimaryNucleon.ModMomentum", PrimaryNucleon_ModMomentum, &b_PrimaryNucleon_ModMomentum);
   fChain->SetBranchAddress("PrimaryNucleon.KE", PrimaryNucleon_KE, &b_PrimaryNucleon_KE);
   fChain->SetBranchAddress("PrimaryNucleon.StartX", PrimaryNucleon_StartX, &b_PrimaryNucleon_StartX);
   fChain->SetBranchAddress("PrimaryNucleon.StartY", PrimaryNucleon_StartY, &b_PrimaryNucleon_StartY);
   fChain->SetBranchAddress("PrimaryNucleon.StartZ", PrimaryNucleon_StartZ, &b_PrimaryNucleon_StartZ);
   fChain->SetBranchAddress("PrimaryNucleon.EndX", PrimaryNucleon_EndX, &b_PrimaryNucleon_EndX);
   fChain->SetBranchAddress("PrimaryNucleon.EndY", PrimaryNucleon_EndY, &b_PrimaryNucleon_EndY);
   fChain->SetBranchAddress("PrimaryNucleon.EndZ", PrimaryNucleon_EndZ, &b_PrimaryNucleon_EndZ);
   fChain->SetBranchAddress("PrimaryNucleon.Travel", PrimaryNucleon_Travel, &b_PrimaryNucleon_Travel);
   fChain->SetBranchAddress("PrimaryNucleon.Theta", PrimaryNucleon_Theta, &b_PrimaryNucleon_Theta);
   fChain->SetBranchAddress("PrimaryNucleon.Phi", PrimaryNucleon_Phi, &b_PrimaryNucleon_Phi);
   fChain->SetBranchAddress("PrimaryNucleon.Origin", PrimaryNucleon_Origin, &b_PrimaryNucleon_Origin);
   fChain->SetBranchAddress("PrimaryPion", &PrimaryPion_, &b_PrimaryPion_);
   fChain->SetBranchAddress("PrimaryPion.PDG", PrimaryPion_PDG, &b_PrimaryPion_PDG);
   fChain->SetBranchAddress("PrimaryPion.E", PrimaryPion_E, &b_PrimaryPion_E);
   fChain->SetBranchAddress("PrimaryPion.Px", PrimaryPion_Px, &b_PrimaryPion_Px);
   fChain->SetBranchAddress("PrimaryPion.Py", PrimaryPion_Py, &b_PrimaryPion_Py);
   fChain->SetBranchAddress("PrimaryPion.Pz", PrimaryPion_Pz, &b_PrimaryPion_Pz);
   fChain->SetBranchAddress("PrimaryPion.ModMomentum", PrimaryPion_ModMomentum, &b_PrimaryPion_ModMomentum);
   fChain->SetBranchAddress("PrimaryPion.KE", PrimaryPion_KE, &b_PrimaryPion_KE);
   fChain->SetBranchAddress("PrimaryPion.StartX", PrimaryPion_StartX, &b_PrimaryPion_StartX);
   fChain->SetBranchAddress("PrimaryPion.StartY", PrimaryPion_StartY, &b_PrimaryPion_StartY);
   fChain->SetBranchAddress("PrimaryPion.StartZ", PrimaryPion_StartZ, &b_PrimaryPion_StartZ);
   fChain->SetBranchAddress("PrimaryPion.EndX", PrimaryPion_EndX, &b_PrimaryPion_EndX);
   fChain->SetBranchAddress("PrimaryPion.EndY", PrimaryPion_EndY, &b_PrimaryPion_EndY);
   fChain->SetBranchAddress("PrimaryPion.EndZ", PrimaryPion_EndZ, &b_PrimaryPion_EndZ);
   fChain->SetBranchAddress("PrimaryPion.Travel", PrimaryPion_Travel, &b_PrimaryPion_Travel);
   fChain->SetBranchAddress("PrimaryPion.Theta", PrimaryPion_Theta, &b_PrimaryPion_Theta);
   fChain->SetBranchAddress("PrimaryPion.Phi", PrimaryPion_Phi, &b_PrimaryPion_Phi);
   fChain->SetBranchAddress("PrimaryPion.Origin", PrimaryPion_Origin, &b_PrimaryPion_Origin);
   fChain->SetBranchAddress("DecayVertex", &DecayVertex, &b_DecayVertex);
   fChain->SetBranchAddress("Decay", &Decay_, &b_Decay_);
   fChain->SetBranchAddress("Decay.PDG", Decay_PDG, &b_Decay_PDG);
   fChain->SetBranchAddress("Decay.E", Decay_E, &b_Decay_E);
   fChain->SetBranchAddress("Decay.Px", Decay_Px, &b_Decay_Px);
   fChain->SetBranchAddress("Decay.Py", Decay_Py, &b_Decay_Py);
   fChain->SetBranchAddress("Decay.Pz", Decay_Pz, &b_Decay_Pz);
   fChain->SetBranchAddress("Decay.ModMomentum", Decay_ModMomentum, &b_Decay_ModMomentum);
   fChain->SetBranchAddress("Decay.KE", Decay_KE, &b_Decay_KE);
   fChain->SetBranchAddress("Decay.StartX", Decay_StartX, &b_Decay_StartX);
   fChain->SetBranchAddress("Decay.StartY", Decay_StartY, &b_Decay_StartY);
   fChain->SetBranchAddress("Decay.StartZ", Decay_StartZ, &b_Decay_StartZ);
   fChain->SetBranchAddress("Decay.EndX", Decay_EndX, &b_Decay_EndX);
   fChain->SetBranchAddress("Decay.EndY", Decay_EndY, &b_Decay_EndY);
   fChain->SetBranchAddress("Decay.EndZ", Decay_EndZ, &b_Decay_EndZ);
   fChain->SetBranchAddress("Decay.Travel", Decay_Travel, &b_Decay_Travel);
   fChain->SetBranchAddress("Decay.Theta", Decay_Theta, &b_Decay_Theta);
   fChain->SetBranchAddress("Decay.Phi", Decay_Phi, &b_Decay_Phi);
   fChain->SetBranchAddress("Decay.Origin", Decay_Origin, &b_Decay_Origin);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton", &SigmaZeroDecayPhoton_, &b_SigmaZeroDecayPhoton_);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.PDG", SigmaZeroDecayPhoton_PDG, &b_SigmaZeroDecayPhoton_PDG);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.E", SigmaZeroDecayPhoton_E, &b_SigmaZeroDecayPhoton_E);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Px", SigmaZeroDecayPhoton_Px, &b_SigmaZeroDecayPhoton_Px);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Py", SigmaZeroDecayPhoton_Py, &b_SigmaZeroDecayPhoton_Py);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Pz", SigmaZeroDecayPhoton_Pz, &b_SigmaZeroDecayPhoton_Pz);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.ModMomentum", SigmaZeroDecayPhoton_ModMomentum, &b_SigmaZeroDecayPhoton_ModMomentum);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.KE", SigmaZeroDecayPhoton_KE, &b_SigmaZeroDecayPhoton_KE);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.StartX", SigmaZeroDecayPhoton_StartX, &b_SigmaZeroDecayPhoton_StartX);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.StartY", SigmaZeroDecayPhoton_StartY, &b_SigmaZeroDecayPhoton_StartY);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.StartZ", SigmaZeroDecayPhoton_StartZ, &b_SigmaZeroDecayPhoton_StartZ);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.EndX", SigmaZeroDecayPhoton_EndX, &b_SigmaZeroDecayPhoton_EndX);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.EndY", SigmaZeroDecayPhoton_EndY, &b_SigmaZeroDecayPhoton_EndY);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.EndZ", SigmaZeroDecayPhoton_EndZ, &b_SigmaZeroDecayPhoton_EndZ);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Travel", SigmaZeroDecayPhoton_Travel, &b_SigmaZeroDecayPhoton_Travel);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Theta", SigmaZeroDecayPhoton_Theta, &b_SigmaZeroDecayPhoton_Theta);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Phi", SigmaZeroDecayPhoton_Phi, &b_SigmaZeroDecayPhoton_Phi);
   fChain->SetBranchAddress("SigmaZeroDecayPhoton.Origin", SigmaZeroDecayPhoton_Origin, &b_SigmaZeroDecayPhoton_Origin);
   fChain->SetBranchAddress("SigmaZeroDecayLambda", &SigmaZeroDecayLambda_, &b_SigmaZeroDecayLambda_);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.PDG", SigmaZeroDecayLambda_PDG, &b_SigmaZeroDecayLambda_PDG);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.E", SigmaZeroDecayLambda_E, &b_SigmaZeroDecayLambda_E);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Px", SigmaZeroDecayLambda_Px, &b_SigmaZeroDecayLambda_Px);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Py", SigmaZeroDecayLambda_Py, &b_SigmaZeroDecayLambda_Py);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Pz", SigmaZeroDecayLambda_Pz, &b_SigmaZeroDecayLambda_Pz);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.ModMomentum", SigmaZeroDecayLambda_ModMomentum, &b_SigmaZeroDecayLambda_ModMomentum);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.KE", SigmaZeroDecayLambda_KE, &b_SigmaZeroDecayLambda_KE);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.StartX", SigmaZeroDecayLambda_StartX, &b_SigmaZeroDecayLambda_StartX);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.StartY", SigmaZeroDecayLambda_StartY, &b_SigmaZeroDecayLambda_StartY);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.StartZ", SigmaZeroDecayLambda_StartZ, &b_SigmaZeroDecayLambda_StartZ);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.EndX", SigmaZeroDecayLambda_EndX, &b_SigmaZeroDecayLambda_EndX);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.EndY", SigmaZeroDecayLambda_EndY, &b_SigmaZeroDecayLambda_EndY);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.EndZ", SigmaZeroDecayLambda_EndZ, &b_SigmaZeroDecayLambda_EndZ);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Travel", SigmaZeroDecayLambda_Travel, &b_SigmaZeroDecayLambda_Travel);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Theta", SigmaZeroDecayLambda_Theta, &b_SigmaZeroDecayLambda_Theta);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Phi", SigmaZeroDecayLambda_Phi, &b_SigmaZeroDecayLambda_Phi);
   fChain->SetBranchAddress("SigmaZeroDecayLambda.Origin", SigmaZeroDecayLambda_Origin, &b_SigmaZeroDecayLambda_Origin);
   fChain->SetBranchAddress("DecayOpeningAngle", &DecayOpeningAngle, &b_DecayOpeningAngle);
   fChain->SetBranchAddress("LeptonPionAngle", &LeptonPionAngle, &b_LeptonPionAngle);
   fChain->SetBranchAddress("LeptonNucleonAngle", &LeptonNucleonAngle, &b_LeptonNucleonAngle);
   fChain->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex, &b_RecoPrimaryVertex);
   fChain->SetBranchAddress("NPrimaryTrackDaughters", &NPrimaryTrackDaughters, &b_NPrimaryTrackDaughters);
   fChain->SetBranchAddress("NPrimaryShowerDaughters", &NPrimaryShowerDaughters, &b_NPrimaryShowerDaughters);
   fChain->SetBranchAddress("TracklikePrimaryDaughters", &TracklikePrimaryDaughters_, &b_TracklikePrimaryDaughters_);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.PDG", TracklikePrimaryDaughters_PDG, &b_TracklikePrimaryDaughters_PDG);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackShowerScore", TracklikePrimaryDaughters_TrackShowerScore, &b_TracklikePrimaryDaughters_TrackShowerScore);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.X", TracklikePrimaryDaughters_X, &b_TracklikePrimaryDaughters_X);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.Y", TracklikePrimaryDaughters_Y, &b_TracklikePrimaryDaughters_Y);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.Z", TracklikePrimaryDaughters_Z, &b_TracklikePrimaryDaughters_Z);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.Displacement", TracklikePrimaryDaughters_Displacement, &b_TracklikePrimaryDaughters_Displacement);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackLength", TracklikePrimaryDaughters_TrackLength, &b_TracklikePrimaryDaughters_TrackLength);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackPID", TracklikePrimaryDaughters_TrackPID, &b_TracklikePrimaryDaughters_TrackPID);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackProtonChi2", TracklikePrimaryDaughters_TrackProtonChi2, &b_TracklikePrimaryDaughters_TrackProtonChi2);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackMuonChi2", TracklikePrimaryDaughters_TrackMuonChi2, &b_TracklikePrimaryDaughters_TrackMuonChi2);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackKaonChi2", TracklikePrimaryDaughters_TrackKaonChi2, &b_TracklikePrimaryDaughters_TrackKaonChi2);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackPionChi2", TracklikePrimaryDaughters_TrackPionChi2, &b_TracklikePrimaryDaughters_TrackPionChi2);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.HasTruth", TracklikePrimaryDaughters_HasTruth, &b_TracklikePrimaryDaughters_HasTruth);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTruePDG", TracklikePrimaryDaughters_TrackTruePDG, &b_TracklikePrimaryDaughters_TrackTruePDG);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTrueE", TracklikePrimaryDaughters_TrackTrueE, &b_TracklikePrimaryDaughters_TrackTrueE);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTruePx", TracklikePrimaryDaughters_TrackTruePx, &b_TracklikePrimaryDaughters_TrackTruePx);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTruePy", TracklikePrimaryDaughters_TrackTruePy, &b_TracklikePrimaryDaughters_TrackTruePy);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTruePz", TracklikePrimaryDaughters_TrackTruePz, &b_TracklikePrimaryDaughters_TrackTruePz);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTrueModMomentum", TracklikePrimaryDaughters_TrackTrueModMomentum, &b_TracklikePrimaryDaughters_TrackTrueModMomentum);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTrueKE", TracklikePrimaryDaughters_TrackTrueKE, &b_TracklikePrimaryDaughters_TrackTrueKE);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTrueLength", TracklikePrimaryDaughters_TrackTrueLength, &b_TracklikePrimaryDaughters_TrackTrueLength);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTrueOrigin", TracklikePrimaryDaughters_TrackTrueOrigin, &b_TracklikePrimaryDaughters_TrackTrueOrigin);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackTrueTotalEdep", TracklikePrimaryDaughters_TrackTrueTotalEdep, &b_TracklikePrimaryDaughters_TrackTrueTotalEdep);
   fChain->SetBranchAddress("TracklikePrimaryDaughters.TrackEdepPurity", TracklikePrimaryDaughters_TrackEdepPurity, &b_TracklikePrimaryDaughters_TrackEdepPurity);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters", &ShowerlikePrimaryDaughters_, &b_ShowerlikePrimaryDaughters_);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.PDG", ShowerlikePrimaryDaughters_PDG, &b_ShowerlikePrimaryDaughters_PDG);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackShowerScore", ShowerlikePrimaryDaughters_TrackShowerScore, &b_ShowerlikePrimaryDaughters_TrackShowerScore);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.X", ShowerlikePrimaryDaughters_X, &b_ShowerlikePrimaryDaughters_X);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.Y", ShowerlikePrimaryDaughters_Y, &b_ShowerlikePrimaryDaughters_Y);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.Z", ShowerlikePrimaryDaughters_Z, &b_ShowerlikePrimaryDaughters_Z);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.Displacement", ShowerlikePrimaryDaughters_Displacement, &b_ShowerlikePrimaryDaughters_Displacement);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackLength", ShowerlikePrimaryDaughters_TrackLength, &b_ShowerlikePrimaryDaughters_TrackLength);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackPID", ShowerlikePrimaryDaughters_TrackPID, &b_ShowerlikePrimaryDaughters_TrackPID);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackProtonChi2", ShowerlikePrimaryDaughters_TrackProtonChi2, &b_ShowerlikePrimaryDaughters_TrackProtonChi2);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackMuonChi2", ShowerlikePrimaryDaughters_TrackMuonChi2, &b_ShowerlikePrimaryDaughters_TrackMuonChi2);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackKaonChi2", ShowerlikePrimaryDaughters_TrackKaonChi2, &b_ShowerlikePrimaryDaughters_TrackKaonChi2);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackPionChi2", ShowerlikePrimaryDaughters_TrackPionChi2, &b_ShowerlikePrimaryDaughters_TrackPionChi2);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.HasTruth", ShowerlikePrimaryDaughters_HasTruth, &b_ShowerlikePrimaryDaughters_HasTruth);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTruePDG", ShowerlikePrimaryDaughters_TrackTruePDG, &b_ShowerlikePrimaryDaughters_TrackTruePDG);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTrueE", ShowerlikePrimaryDaughters_TrackTrueE, &b_ShowerlikePrimaryDaughters_TrackTrueE);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTruePx", ShowerlikePrimaryDaughters_TrackTruePx, &b_ShowerlikePrimaryDaughters_TrackTruePx);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTruePy", ShowerlikePrimaryDaughters_TrackTruePy, &b_ShowerlikePrimaryDaughters_TrackTruePy);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTruePz", ShowerlikePrimaryDaughters_TrackTruePz, &b_ShowerlikePrimaryDaughters_TrackTruePz);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTrueModMomentum", ShowerlikePrimaryDaughters_TrackTrueModMomentum, &b_ShowerlikePrimaryDaughters_TrackTrueModMomentum);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTrueKE", ShowerlikePrimaryDaughters_TrackTrueKE, &b_ShowerlikePrimaryDaughters_TrackTrueKE);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTrueLength", ShowerlikePrimaryDaughters_TrackTrueLength, &b_ShowerlikePrimaryDaughters_TrackTrueLength);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTrueOrigin", ShowerlikePrimaryDaughters_TrackTrueOrigin, &b_ShowerlikePrimaryDaughters_TrackTrueOrigin);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackTrueTotalEdep", ShowerlikePrimaryDaughters_TrackTrueTotalEdep, &b_ShowerlikePrimaryDaughters_TrackTrueTotalEdep);
   fChain->SetBranchAddress("ShowerlikePrimaryDaughters.TrackEdepPurity", ShowerlikePrimaryDaughters_TrackEdepPurity, &b_ShowerlikePrimaryDaughters_TrackEdepPurity);
   fChain->SetBranchAddress("MuonIndex", &MuonIndex, &b_MuonIndex);
   Notify();
}

Bool_t OutputTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void OutputTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t OutputTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef OutputTree_cxx

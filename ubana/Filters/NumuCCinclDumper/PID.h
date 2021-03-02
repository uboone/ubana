#ifndef PID_H
#define PID_H

//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "TVector3.h"

class PID {

	public:

	  // Default constructor
	  PID(){}

	  // Default destructor
	  ~PID(){}

          void Chi2(art::FindManyP<anab::ParticleID> PIDTotrackAsso, art::Ptr<recob::Track> track, TVector3 Trk_start_SCEcorr, TVector3 Trk_end_SCEcorr, int hits_dEdx_size_pl0 = -1, int hits_dEdx_size_pl1 = -1, int hits_dEdx_size_pl2 = -1);

          double PID_Chi2Mu_pl0; // Chi2 of muon assumption of plane 0 in PID
          double PID_Chi2Mu_pl1; // Chi2 of muon assumption of plane 1 in PID
          double PID_Chi2Mu_pl2; // Chi2 of muon assumption of plane 2 in PID
          double PID_Chi2Mu_3pl; // Chi2 of muon assumption of 3 planes in PID

          double PID_Chi2P_pl0; // Chi2 of proton assumption of plane 0 in PID
          double PID_Chi2P_pl1; // Chi2 of proton assumption of plane 1 in PID
          double PID_Chi2P_pl2; // Chi2 of proton assumption of plane 2 in PID
          double PID_Chi2P_3pl; // Chi2 of proton assumption of 3 planes in PID

          double PID_Chi2Pi_pl0; // Chi2 of pion assumption of plane 0 in PID
          double PID_Chi2Pi_pl1; // Chi2 of pion assumption of plane 1 in PID
          double PID_Chi2Pi_pl2; // Chi2 of pion assumption of plane 2 in PID
          double PID_Chi2Pi_3pl; // Chi2 of pion assumption of 3 planes in PID

          double PID_Chi2K_pl0; // Chi2 of kaon assumption of plane 0 in PID
          double PID_Chi2K_pl1; // Chi2 of kaon assumption of plane 1 in PID
          double PID_Chi2K_pl2; // Chi2 of kaon assumption of plane 2 in PID
          double PID_Chi2K_3pl; // Chi2 of kaon assumption of 3 planes in PID

          int BestPlane_PID;
          bool Pl2_for_PID;
          bool Pl1_for_PID;
          bool Pl0_for_PID;

          int PID_Pdg_3pl; //[Only fill positive value] The Pdg of the corresponding particle assumption with minimum Chi2
          int PID_Pdg_pl2;
          int PID_Pdg_pl1;
          int PID_Pdg_pl0;
          double PID_avg_Chi2; // Minimum averaged Chi2 of 3 planes among all assumptions
          double PID_pl2_Chi2;
          double PID_pl1_Chi2;
          double PID_pl0_Chi2;



	private:


};

#endif

#ifndef BACKTRACKERTRUTHMATCH_H
#define BACKTRACKERTRUTHMATCH_H

//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

//some ROOT includes
#include "TString.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TStopwatch.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/PtrVectorBase.h"

//"larsoft" object includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

class BackTrackerTruthMatch {

	public:

		// Default constructor
		BackTrackerTruthMatch(){}

		// Default destructor
		~BackTrackerTruthMatch(){}

		void MatchToMCParticle(const art::Handle<std::vector<recob::Hit> >& hit_handle, const art::Event& e, std::vector<art::Ptr<recob::Hit> >& trk_hits_ptrs);
		art::Ptr< simb::MCParticle > ReturnMCParticle();  
		bool ParticleAlreadyMatchedInThisHit(std::vector<int> AlreadyMatched_TrackIDs ,int cTrackID);
		double ReturnPurity();
		double ReturnCompleteness();
		double ReturnTrueAssDepositedEnergy();
		double ReturnTotalTrueDepositedEnergy();
		int ReturnMCParticleID();

	private:

		art::Ptr< simb::MCParticle > fmaxp_me;
		double fpurity = -999;
		double fcompleteness = -999;
		double ftote = 0;
		double fmaxe = 0;
		int fMCParticleID;

};

#endif

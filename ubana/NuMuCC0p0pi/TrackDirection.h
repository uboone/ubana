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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


class TrackDirection {

	public:

		// Default constructor
		TrackDirection(){}

		// Default destructor
		~TrackDirection(){}

                void TrackDir(art::Ptr<recob::Track>& ThisTrack, art::FindManyP<anab::Calorimetry>& trkToCalAsso);

	private:


};

#endif

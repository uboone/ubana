#ifndef BROKENTRACK_H
#define BROKENTRACK_H

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
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <algorithm>

#include "FiducialVolume.h"

class BrokenTrack {

	public:

	  // Default constructor
	  BrokenTrack(){}

	  // Default destructor
	  ~BrokenTrack(){}

          void MatchTracks(art::Ptr<recob::Track>& ThisTrack, std::vector< art::Ptr<recob::Track>>& TrackCollection);

          bool NewTrk();
          int NumberMergedTracks();
 
          double TrkLen();

          TVector3 TrkEnd1();
          TVector3 TrkEnd2();

          //bool NewTrkFV();
          //bool NewTrkContained();

	private:
          spacecharge::SpaceCharge const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

          //::ubana::FiducialVolume _fiducial_volume; 

          TVector3 trk_end1;
          TVector3 trk_end2;

          double trk_length;

          TVector3 trk_temp_end1;
          TVector3 trk_temp_end2;

          double trk_temp_length;
          
          bool newTrk;
          int Nr_mergedTrk;
          //bool newTrk_FV;           
          //bool newTrk_contained;           
};

#endif

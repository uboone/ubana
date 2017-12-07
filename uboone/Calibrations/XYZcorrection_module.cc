#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

const int kMaxTracks  = 100;

using namespace std;

namespace microboone{

class XYZcorrection : public art::EDAnalyzer {
public:

    explicit XYZcorrection(fhicl::ParameterSet const& pset);
    virtual ~XYZcorrection();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
private:
    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    cross_trks;
    Int_t    all_trks;
    Float_t  xprojectedlen[kMaxTracks];
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Int_t    TrkID[kMaxTracks]; 
    Float_t  trkstartcosxyz[kMaxTracks][3];
    Float_t  trkendcosxyz[kMaxTracks][3];
    Int_t    ntrkhits[kMaxTracks][3];
    Float_t  trkdqdx[kMaxTracks][3][3000];
    Float_t  trkdedx[kMaxTracks][3][3000];
    Float_t  trkresrange[kMaxTracks][3][3000];
    Float_t  trkhitx[kMaxTracks][3][3000];
    Float_t  trkhity[kMaxTracks][3][3000];
    Float_t  trkhitz[kMaxTracks][3][3000];
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    std::string fClusterModuleLabel;
    std::string fVertexModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fSaveClusterInfo;
    bool  fSaveClusterHitInfo;
    bool  fSaveGenieInfo;
    bool  fSaveVertexInfo; 
    float fG4minE;
}; 

//========================================================================
XYZcorrection::XYZcorrection(fhicl::ParameterSet const& pset) :
EDAnalyzer(pset),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel","")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel","")     ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","")     ),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel","")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel","")   ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel","")      ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel","")       ),
  fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
  fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),  
  fSaveClusterInfo          (pset.get< bool>("SaveClusterInfo",false)),
  fSaveClusterHitInfo       (pset.get< bool>("SaveClusterHitInfo",false)),
  fSaveGenieInfo            (pset.get< bool>("SaveGenieInfo",false)),
  fSaveVertexInfo           (pset.get< bool>("SaveVertexInfo",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01))  
{
  if (fSaveTrackInfo == false) fSaveCaloInfo = false;
}
 
//========================================================================
XYZcorrection::~XYZcorrection(){
}
//========================================================================

//========================================================================
void XYZcorrection::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  fEventTree->Branch("event", &event,"event/I");
  fEventTree->Branch("evttime",&evttime,"evttime/D");
  fEventTree->Branch("run", &run,"run/I");
  fEventTree->Branch("subrun", &subrun,"surbrun/I");
  fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
  fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
  fEventTree->Branch("cross_trks",&cross_trks,"cross_trks/I");
  fEventTree->Branch("all_trks",&all_trks,"all_trks/I");
  fEventTree->Branch("xprojectedlen",xprojectedlen,"xprojectedlen[all_trks]/F");
  fEventTree->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[cross_trks]/F");
  fEventTree->Branch("trackthetayz",trackthetayz,"trackthetayz[cross_trks]/F");
  fEventTree->Branch("TrkID",TrkID,"TrkID[cross_trks]/I");
  fEventTree->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[cross_trks][3]/F");
  fEventTree->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[cross_trks][3]/F");
  fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[cross_trks][3]/I");
  fEventTree->Branch("trkdqdx",trkdqdx,"trkdqdx[cross_trks][3][3000]/F");
  fEventTree->Branch("trkdedx",trkdedx,"trkdedx[cross_trks][3][3000]/F");
  fEventTree->Branch("trkresrange",trkresrange,"trkresrange[cross_trks][3][3000]/F");
  fEventTree->Branch("trkhitx",trkhitx,"trkhitx[cross_trks][3][3000]/F");
  fEventTree->Branch("trkhity",trkhity,"trkhity[cross_trks][3][3000]/F");
  fEventTree->Branch("trkhitz",trkhitz,"trkhitz[cross_trks][3][3000]/F");
}

//========================================================================
void XYZcorrection::endJob(){     

}

//========================================================================
void XYZcorrection::beginRun(const art::Run&){
  mf::LogInfo("XYZcorrection")<<"begin run..."<<std::endl;
}
//========================================================================

//========================================================================

//========================================================================

void XYZcorrection::analyze( const art::Event& evt){
     reset();  
     
     art::Handle< std::vector<recob::Track> > trackListHandle;
     std::vector<art::Ptr<recob::Track> > tracklist;
  
     if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
  
     art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
     
     run = evt.run();
     subrun = evt.subRun();
     event = evt.id().event();
     art::Timestamp ts = evt.time();
     TTimeStamp tts(ts.timeHigh(), ts.timeLow());
     evttime=tts.AsDouble();
     
     UInt_t year=0;
     UInt_t month=0;
     UInt_t day=0;
     
     year_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);
     
     UInt_t hour=0;
     UInt_t min=0;
     UInt_t sec=0;
     
     hour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);
  
     cross_trks=0;
     all_trks=0;
  
     size_t NTracks = tracklist.size();
     for(size_t i=0; i<NTracks;++i){
         art::Ptr<recob::Track> ptrack(trackListHandle, i);
	 std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);

	 const recob::Track& track = *ptrack;
	 TVector3 pos, dir_start, dir_end, end;
	 pos = track.Vertex();
     	 dir_start = track.VertexDirection();
     	 dir_end   = track.EndDirection();
     	 end = track.End();
	 double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
         double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
     	 float X = std::abs(pos.X()-end.X());
	 all_trks++;
	 xprojectedlen[all_trks-1]=X;
	 if(X>250 && X<270){
	    cross_trks++;
	    trackthetaxz[cross_trks-1]=theta_xz;
	    trackthetayz[cross_trks-1]=theta_yz;
	    TrkID[cross_trks-1]=track.ID();
	    trkstartcosxyz[cross_trks-1][0]=dir_start.X();
	    trkstartcosxyz[cross_trks-1][1]=dir_start.Y();
	    trkstartcosxyz[cross_trks-1][2]=dir_start.Z();
	    trkendcosxyz[cross_trks-1][0]=dir_end.X();
	    trkendcosxyz[cross_trks-1][1]=dir_end.Y();
	    trkendcosxyz[cross_trks-1][2]=dir_end.Z();
	    
	    ////////////// Getting min X //////////////////
				   
	    double minx_p2 = 1e10;
	    double minx_p1 = 1e10;
	    double minx_p0 = 1e10;
	    for(size_t ical = 0; ical<calos.size(); ++ical){
		if(!calos[ical]) continue;
	        if(!calos[ical]->PlaneID().isValid) continue;
	        int planenum = calos[ical]->PlaneID().Plane;
	        if(planenum<0||planenum>2) continue;
	        if(planenum==2){
		   const size_t NHits = calos[ical] -> dEdx().size();
		   for(size_t iHit = 0; iHit < NHits; ++iHit){
		       const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
		       if(TrkPos.X()>-180){
			  if(TrkPos.X()<minx_p2) minx_p2=TrkPos.X();
		       }
	            }
	        }
	        //////////////////// end p2 ////////////////////////////
		if(planenum == 1){
		   const size_t NHits = calos[ical] -> dEdx().size();
		   for(size_t iHit = 0; iHit < NHits; ++iHit){
		       const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
		       if(TrkPos.X()>-180){
			  if(TrkPos.X()<minx_p1) minx_p1=TrkPos.X();
		       }
	            }
	         }
		 ////////////////// end p1 ////////////////////////////
		 if(planenum == 0){
		    const size_t NHits = calos[ical] -> dEdx().size();
		    for(size_t iHit = 0; iHit < NHits; ++iHit){
		        const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
		        if(TrkPos.X()>-180){
			   if(TrkPos.X()<minx_p0) minx_p0=TrkPos.X();
			}
	             }
	         }
		   ///////////////// end p0 ////////////////////////////
	      }
				   
	    ////////////////// End of min X //////////////
				   
	      for(size_t ical = 0; ical<calos.size(); ++ical){
	          if(!calos[ical]) continue;
	          if(!calos[ical]->PlaneID().isValid) continue;
	          int planenum = calos[ical]->PlaneID().Plane;
	          if(planenum<0||planenum>2) continue;
	          const size_t NHits = calos[ical] -> dEdx().size();
		  ntrkhits[cross_trks-1][planenum]=int(NHits);
		  for(size_t iHit = 0; iHit < NHits; ++iHit){
		      const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
		      double minx=0;
		      if(planenum==0) minx=minx_p0;
		      if(planenum==1) minx=minx_p1;
		      if(planenum==2) minx=minx_p2;
	              trkdqdx[cross_trks-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	              trkdedx[cross_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
		      trkresrange[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
		      trkhitx[cross_trks-1][planenum][iHit]=TrkPos.X()-minx;
		      trkhity[cross_trks-1][planenum][iHit]=TrkPos.Y();
		      trkhitz[cross_trks-1][planenum][iHit]=TrkPos.Z();
		 } // loop over iHit..
	     } // loop over ical 2 nd time...
          } // crossing trks...
      } // loop over trks...
   fEventTree->Fill();
} // end of analyze function
	   
 /////////////////// Defintion of reset function ///////////
void XYZcorrection::reset(){
     run = -9999;
     subrun = -9999;
     event = -9999;
     evttime = -9999;
     cross_trks = -9999;
     all_trks = -9999;
     year_month_date=-9999;
     hour_min_sec=-9999;
     for(int i=0; i<kMaxTracks; i++){
         trackthetaxz[i]=-9999;
	 trackthetayz[i]=-9999;
	 TrkID[i]=-9999;
	 xprojectedlen[i]=-9999;
	 trkstartcosxyz[i][0]=-9999;
	 trkstartcosxyz[i][1]=-9999;
	 trkstartcosxyz[i][2]=-9999; 
	 trkendcosxyz[i][0]=-9999;
	 trkendcosxyz[i][1]=-9999;
	 trkendcosxyz[i][2]=-9999;
	 ntrkhits[i][0] = -9999;
	 ntrkhits[i][1] = -9999;
	 ntrkhits[i][2] = -9999;
	 for(int j=0; j<3; j++){
	      for(int k=0; k<3000; k++){
	           trkdqdx[i][j][k]=-9999;
		   trkdedx[i][j][k]=-9999;
		   trkresrange[i][j][k]=-9999;
		   trkhitx[i][j][k]=-9999;
		   trkhity[i][j][k]=-9999;
		   trkhitz[i][j][k]=-9999;
	      }
	 }
     }
}
 //////////////////////// End of definition ///////////////	
	  
DEFINE_ART_MODULE(XYZcorrection)
}



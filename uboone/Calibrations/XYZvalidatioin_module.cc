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
#include "nusimdata/SimulationBase/MCParticle.h"
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
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
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
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

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

const int kMaxTracks=100;

using namespace std;

namespace microboone{

class XYZvalidatioin : public art::EDAnalyzer {
public:

    explicit XYZvalidatioin(fhicl::ParameterSet const& pset);
    virtual ~XYZvalidatioin();

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
    Int_t    stop_mu;
    Float_t  mu_trklen[kMaxTracks]; 
    Float_t  mu_trackthetaxz[kMaxTracks];
    Float_t  mu_trackthetayz[kMaxTracks];
    Int_t    mu_TrkID[kMaxTracks]; 
    Float_t  mu_start_info[kMaxTracks][3];
    Float_t  mu_end_info[kMaxTracks][3];
    Float_t  mu_trkstartcosxyz[kMaxTracks][3];
    Float_t  mu_trkendcosxyz[kMaxTracks][3];
    Int_t    mu_ntrkhits[kMaxTracks][3];
    Float_t  mu_def_pida[kMaxTracks][3];
    Float_t  mu_truth_endE[kMaxTracks];
    Float_t  mu_trkdqdx[kMaxTracks][3][3000];
    Float_t  mu_trkdedx[kMaxTracks][3][3000];
    Float_t  mu_trkresrange[kMaxTracks][3][3000];
    Float_t  mu_trkhitx[kMaxTracks][3][3000];
    Float_t  mu_trkhity[kMaxTracks][3][3000];
    Float_t  mu_trkhitz[kMaxTracks][3][3000];
    Float_t  mu_trkpitch[kMaxTracks][3][3000];
    Float_t  mu_truth_start_info[kMaxTracks][3];
    Float_t  mu_truth_end_info[kMaxTracks][3];
    Int_t    stop_pi;
    Float_t  pi_trklen[kMaxTracks]; 
    Float_t  pi_trackthetaxz[kMaxTracks];
    Float_t  pi_trackthetayz[kMaxTracks];
    Int_t    pi_TrkID[kMaxTracks]; 
    Float_t  pi_start_info[kMaxTracks][3];
    Float_t  pi_end_info[kMaxTracks][3];
    Float_t  pi_trkstartcosxyz[kMaxTracks][3];
    Float_t  pi_trkendcosxyz[kMaxTracks][3];
    Int_t    pi_ntrkhits[kMaxTracks][3];
    Float_t  pi_def_pida[kMaxTracks][3];
    Float_t  pi_truth_endE[kMaxTracks];
    Float_t  pi_trkdqdx[kMaxTracks][3][3000];
    Float_t  pi_trkdedx[kMaxTracks][3][3000];
    Float_t  pi_trkresrange[kMaxTracks][3][3000];
    Float_t  pi_trkhitx[kMaxTracks][3][3000];
    Float_t  pi_trkhity[kMaxTracks][3][3000];
    Float_t  pi_trkhitz[kMaxTracks][3][3000];
    Float_t  pi_trkpitch[kMaxTracks][3][3000];
    Float_t  pi_truth_start_info[kMaxTracks][3];
    Float_t  pi_truth_end_info[kMaxTracks][3];
    Int_t    stop_k;
    Float_t  k_trklen[kMaxTracks];
    Float_t  k_trackthetaxz[kMaxTracks];
    Float_t  k_trackthetayz[kMaxTracks];
    Int_t    k_TrkID[kMaxTracks]; 
    Float_t  k_start_info[kMaxTracks][3];
    Float_t  k_end_info[kMaxTracks][3];
    Float_t  k_trkstartcosxyz[kMaxTracks][3];
    Float_t  k_trkendcosxyz[kMaxTracks][3];
    Int_t    k_ntrkhits[kMaxTracks][3];
    Float_t  k_def_pida[kMaxTracks][3];
    Float_t  k_truth_endE[kMaxTracks];
    Float_t  k_trkdqdx[kMaxTracks][3][3000];
    Float_t  k_trkdedx[kMaxTracks][3][3000];
    Float_t  k_trkresrange[kMaxTracks][3][3000];
    Float_t  k_trkhitx[kMaxTracks][3][3000];
    Float_t  k_trkhity[kMaxTracks][3][3000];
    Float_t  k_trkhitz[kMaxTracks][3][3000];
    Float_t  k_trkpitch[kMaxTracks][3][3000];
    Float_t  k_truth_start_info[kMaxTracks][3];
    Float_t  k_truth_end_info[kMaxTracks][3];
    Int_t    stop_p;
    Float_t  p_trklen[kMaxTracks];
    Float_t  p_trackthetaxz[kMaxTracks];
    Float_t  p_trackthetayz[kMaxTracks];
    Int_t    p_TrkID[kMaxTracks]; 
    Float_t  p_start_info[kMaxTracks][3];
    Float_t  p_end_info[kMaxTracks][3];
    Float_t  p_trkstartcosxyz[kMaxTracks][3];
    Float_t  p_trkendcosxyz[kMaxTracks][3];
    Int_t    p_ntrkhits[kMaxTracks][3];
    Float_t  p_def_pida[kMaxTracks][3];
    Float_t  p_truth_endE[kMaxTracks];
    Float_t  p_trkdqdx[kMaxTracks][3][3000];
    Float_t  p_trkdedx[kMaxTracks][3][3000];
    Float_t  p_trkresrange[kMaxTracks][3][3000];
    Float_t  p_trkhitx[kMaxTracks][3][3000];
    Float_t  p_trkhity[kMaxTracks][3][3000];
    Float_t  p_trkhitz[kMaxTracks][3][3000];
    Float_t  p_trkpitch[kMaxTracks][3][3000];
    Float_t  p_truth_start_info[kMaxTracks][3];
    Float_t  p_truth_end_info[kMaxTracks][3];
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
}; // class XYZvalidatioin

//========================================================================
XYZvalidatioin::XYZvalidatioin(fhicl::ParameterSet const& pset) :
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
XYZvalidatioin::~XYZvalidatioin(){
  //destructor
}
//========================================================================

//========================================================================
void XYZvalidatioin::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  fEventTree->Branch("event", &event,"event/I");
  fEventTree->Branch("run", &run,"run/I");
  fEventTree->Branch("subrun", &subrun,"surbrun/I");
  fEventTree->Branch("stop_mu",&stop_mu,"stop_mu/I");
  fEventTree->Branch("mu_trklen",mu_trklen,"mu_trklen[stop_mu]/F");
  fEventTree->Branch("mu_trackthetaxz",mu_trackthetaxz,"mu_trackthetaxz[stop_mu]/F");
  fEventTree->Branch("mu_trackthetayz",mu_trackthetayz,"mu_trackthetayz[stop_mu]/F");
  fEventTree->Branch("mu_TrkID",mu_TrkID,"mu_TrkID[stop_mu]/I");
  fEventTree->Branch("mu_start_info",mu_start_info,"mu_start_info[stop_mu][3]/F");
  fEventTree->Branch("mu_end_info",mu_end_info,"mu_end_info[stop_mu][3]/F");
  fEventTree->Branch("mu_trkstartcosxyz",mu_trkstartcosxyz,"mu_trkstartcosxyz[stop_mu][3]/F");
  fEventTree->Branch("mu_trkendcosxyz",mu_trkendcosxyz,"mu_trkendcosxyz[stop_mu][3]/F");
  fEventTree->Branch("mu_ntrkhits",mu_ntrkhits,"mu_ntrkhits[stop_mu][3]/I");
  fEventTree->Branch("mu_def_pida",mu_def_pida,"mu_def_pida[stop_mu][3]/F");
  fEventTree->Branch("mu_truth_endE",mu_truth_endE,"mu_truth_endE[stop_mu]/F");
  fEventTree->Branch("mu_trkdqdx",mu_trkdqdx,"mu_trkdqdx[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_trkdedx",mu_trkdedx,"mu_trkdedx[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_trkresrange",mu_trkresrange,"mu_trkresrange[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_trkhitx",mu_trkhitx,"mu_trkhitx[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_trkhity",mu_trkhity,"mu_trkhity[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_trkhitz",mu_trkhitz,"mu_trkhitz[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_trkpitch",mu_trkpitch,"mu_trkpitch[stop_mu][3][3000]/F");
  fEventTree->Branch("mu_truth_start_info",mu_truth_start_info,"mu_truth_start_info[stop_mu][3]/F");
  fEventTree->Branch("mu_truth_end_info",mu_truth_end_info,"mu_truth_end_info[stop_mu][3]/F");
  fEventTree->Branch("stop_pi",&stop_pi,"stop_pi/I");
  fEventTree->Branch("pi_trklen",pi_trklen,"pi_trklen[stop_pi]/F");
  fEventTree->Branch("pi_trackthetaxz",pi_trackthetaxz,"pi_trackthetaxz[stop_pi]/F");
  fEventTree->Branch("pi_trackthetayz",pi_trackthetayz,"pi_trackthetayz[stop_pi]/F");
  fEventTree->Branch("pi_TrkID",pi_TrkID,"pi_TrkID[stop_pi]/I");
  fEventTree->Branch("pi_start_info",pi_start_info,"pi_start_info[stop_pi][3]/F");
  fEventTree->Branch("pi_end_info",pi_end_info,"pi_end_info[stop_pi][3]/F");
  fEventTree->Branch("pi_trkstartcosxyz",pi_trkstartcosxyz,"pi_trkstartcosxyz[stop_pi][3]/F");
  fEventTree->Branch("pi_trkendcosxyz",pi_trkendcosxyz,"pi_trkendcosxyz[stop_pi][3]/F");
  fEventTree->Branch("pi_ntrkhits",pi_ntrkhits,"pi_ntrkhits[stop_pi][3]/I");
  fEventTree->Branch("pi_def_pida",pi_def_pida,"pi_def_pida[stop_pi][3]/F");
  fEventTree->Branch("pi_truth_endE",pi_truth_endE,"pi_truth_endE[stop_pi]/F");
  fEventTree->Branch("pi_trkdqdx",pi_trkdqdx,"pi_trkdqdx[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_trkdedx",pi_trkdedx,"pi_trkdedx[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_trkresrange",pi_trkresrange,"pi_trkresrange[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_trkhitx",pi_trkhitx,"pi_trkhitx[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_trkhity",pi_trkhity,"pi_trkhity[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_trkhitz",pi_trkhitz,"pi_trkhitz[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_trkpitch",pi_trkpitch,"pi_trkpitch[stop_pi][3][3000]/F");
  fEventTree->Branch("pi_truth_start_info",pi_truth_start_info,"pi_truth_start_info[stop_pi][3]/F");
  fEventTree->Branch("pi_truth_end_info",pi_truth_end_info,"pi_truth_end_info[stop_pi][3]/F");
  fEventTree->Branch("stop_k",&stop_k,"stop_k/I");
  fEventTree->Branch("k_trklen",k_trklen,"k_trklen[stop_k]/F");
  fEventTree->Branch("k_trackthetaxz",k_trackthetaxz,"k_trackthetaxz[stop_k]/F");
  fEventTree->Branch("k_trackthetayz",k_trackthetayz,"k_trackthetayz[stop_k]/F");
  fEventTree->Branch("k_TrkID",k_TrkID,"k_TrkID[stop_k]/I");
  fEventTree->Branch("k_start_info",k_start_info,"k_start_info[stop_k][3]/F");
  fEventTree->Branch("k_end_info",k_end_info,"k_end_info[stop_k][3]/F");
  fEventTree->Branch("k_trkstartcosxyz",k_trkstartcosxyz,"k_trkstartcosxyz[stop_k][3]/F");
  fEventTree->Branch("k_trkendcosxyz",k_trkendcosxyz,"k_trkendcosxyz[stop_k][3]/F");
  fEventTree->Branch("k_ntrkhits",k_ntrkhits,"k_ntrkhits[stop_k][3]/I");
  fEventTree->Branch("k_def_pida",k_def_pida,"k_def_pida[stop_k][3]/F");
  fEventTree->Branch("k_truth_endE",k_truth_endE,"k_truth_endE[stop_k]/F");
  fEventTree->Branch("k_trkdqdx",k_trkdqdx,"k_trkdqdx[stop_k][3][3000]/F");
  fEventTree->Branch("k_trkdedx",k_trkdedx,"k_trkdedx[stop_k][3][3000]/F");
  fEventTree->Branch("k_trkresrange",k_trkresrange,"k_trkresrange[stop_k][3][3000]/F");
  fEventTree->Branch("k_trkhitx",k_trkhitx,"k_trkhitx[stop_k][3][3000]/F");
  fEventTree->Branch("k_trkhity",k_trkhity,"k_trkhity[stop_k][3][3000]/F");
  fEventTree->Branch("k_trkhitz",k_trkhitz,"k_trkhitz[stop_k][3][3000]/F");
  fEventTree->Branch("k_trkpitch",k_trkpitch,"k_trkpitch[stop_k][3][3000]/F");
  fEventTree->Branch("k_truth_start_info",k_truth_start_info,"k_truth_start_info[stop_k][3]/F");
  fEventTree->Branch("k_truth_end_info",k_truth_end_info,"k_truth_end_info[stop_k][3]/F");
  fEventTree->Branch("stop_p",&stop_p,"stop_p/I");
  fEventTree->Branch("p_trklen",p_trklen,"p_trklen[stop_p]/F");
  fEventTree->Branch("p_trackthetaxz",p_trackthetaxz,"p_trackthetaxz[stop_p]/F");
  fEventTree->Branch("p_trackthetayz",p_trackthetayz,"p_trackthetayz[stop_p]/F");
  fEventTree->Branch("p_TrkID",p_TrkID,"p_TrkID[stop_p]/I");
  fEventTree->Branch("p_start_info",p_start_info,"p_start_info[stop_p][3]/F");
  fEventTree->Branch("p_end_info",p_end_info,"p_end_info[stop_p][3]/F");
  fEventTree->Branch("p_trkstartcosxyz",p_trkstartcosxyz,"p_trkstartcosxyz[stop_p][3]/F");
  fEventTree->Branch("p_trkendcosxyz",p_trkendcosxyz,"p_trkendcosxyz[stop_p][3]/F");
  fEventTree->Branch("p_ntrkhits",p_ntrkhits,"p_ntrkhits[stop_p][3]/I");
  fEventTree->Branch("p_def_pida",p_def_pida,"p_def_pida[stop_p][3]/F");
  fEventTree->Branch("p_truth_endE",p_truth_endE,"p_truth_endE[stop_p]/F");
  fEventTree->Branch("p_trkdqdx",p_trkdqdx,"p_trkdqdx[stop_p][3][3000]/F");
  fEventTree->Branch("p_trkdedx",p_trkdedx,"p_trkdedx[stop_p][3][3000]/F");
  fEventTree->Branch("p_trkresrange",p_trkresrange,"p_trkresrange[stop_p][3][3000]/F");
  fEventTree->Branch("p_trkhitx",p_trkhitx,"p_trkhitx[stop_p][3][3000]/F");
  fEventTree->Branch("p_trkhity",p_trkhity,"p_trkhity[stop_p][3][3000]/F");
  fEventTree->Branch("p_trkhitz",p_trkhitz,"p_trkhitz[stop_p][3][3000]/F");
  fEventTree->Branch("p_trkpitch",p_trkpitch,"p_trkpitch[stop_p][3][3000]/F");
  fEventTree->Branch("p_truth_start_info",p_truth_start_info,"p_truth_start_info[stop_p][3]/F");
  fEventTree->Branch("p_truth_end_info",p_truth_end_info,"p_truth_end_info[stop_p][3]/F");
}

//========================================================================
void XYZvalidatioin::endJob(){     

}

//========================================================================
void XYZvalidatioin::beginRun(const art::Run&){
  mf::LogInfo("XYZvalidatioin")<<"begin run..."<<std::endl;
}
//========================================================================

//========================================================================

//========================================================================

void XYZvalidatioin::analyze( const art::Event& evt){
     reset();  
     bool isMC = !evt.isRealData();
        
     art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
     std::vector<art::Ptr<simb::MCTruth> > mclist;
     if(isMC){
           if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
           art::fill_ptr_vector(mclist, mctruthListHandle);
	   
     } 
  
     art::Handle< std::vector<recob::Track> > trackListHandle;
     std::vector<art::Ptr<recob::Track> > tracklist;
     if(fSaveTrackInfo){
         if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
         art::fill_ptr_vector(tracklist, trackListHandle);
     }
  
     art::Handle< std::vector<recob::Hit> > hitListHandle;
     std::vector<art::Ptr<recob::Hit> > hitlist;
     if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
        art::fill_ptr_vector(hitlist, hitListHandle);

     //art::ServiceHandle<cheat::BackTracker> bt;
     
     art::ServiceHandle<cheat::BackTrackerService> bt_serv;
     art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
     
     art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
     art::FindManyP<recob::Hit> fmht(trackListHandle,evt,fTrackModuleLabel);
     art::FindManyP<anab::ParticleID> fmpid(trackListHandle, evt, fParticleIDModuleLabel);
     art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> fmhitmc(hitListHandle,evt,"crHitRemovalTruthMatch");
  
     run = evt.run();
     subrun = evt.subRun();
     event = evt.id().event();  
  
     stop_mu=0;
     stop_pi=0;
     stop_k=0;
     stop_p=0;
     
     if(isMC){
        if(fSaveTrackInfo){
           size_t NTracks = tracklist.size();
	   for(size_t i=0; i<NTracks;++i){
	       art::Ptr<recob::Track> ptrack(trackListHandle, i);
	       std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
	       std::vector<art::Ptr<anab::ParticleID>> pids = fmpid.at(i);
	       const recob::Track& track = *ptrack;
	       TVector3 pos, dir_start, dir_end, end;
	       pos = track.Vertex();
     	       dir_start = track.VertexDirection();
     	       dir_end   = track.EndDirection();
     	       end = track.End();
	       double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
               double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
			 
	       if((pos.X()>10)&&(pos.X()<250)&&(end.X()>10)&&(end.X()<250)&&(pos.Y()>-110)&&(pos.Y()<110)&&(end.Y()>-110)&&(end.Y()<110)&&(pos.Z()>10)&&(pos.Z()<1030)&&(end.Z()>10)&&(end.Z()<1030)){
			     
	           float rec_stx=pos.X();
		   float rec_sty=pos.Y();
		   float rec_stz=pos.Z();
		   float rec_enx=end.X();
		   float rec_eny=end.Y();
		   float rec_enz=end.Z();
		   
		   std::vector< art::Ptr<recob::Hit> > allKHits = fmht.at(i);
	           std::map<int,double> trk_mu_ide;
	           for(size_t h = 0; h < allKHits.size(); ++h){
		       art::Ptr<recob::Hit> hit = allKHits[h];
		       auto particles = fmhitmc.at(hit.key());
		       auto hitmatch = fmhitmc.data(hit.key());
		       for(size_t e = 0; e<particles.size(); ++e){
		       	 if (!particles[e]) continue;
			 if (!pi_serv->TrackIdToMotherParticle_P(particles[e]->TrackId())) continue;
		       	 size_t trkid = (pi_serv->TrackIdToMotherParticle_P(particles[e]->TrackId()))->TrackId();
		         trk_mu_ide[trkid] += hitmatch[e]->energy;
		       }
	           }
			 
		   double maxke = -1;
	           double totke = 0;
		   int Track_mu_id = 0;
	           for(std::map<int,double>::iterator ii=trk_mu_ide.begin();ii!=trk_mu_ide.end(); ++ii){
		       totke += ii->second;
		       if((ii->second)>maxke){
		           maxke = ii->second;
		           Track_mu_id=ii->first;
		       }
		    }
		     
		    int origin= -1;
		    const simb::MCParticle* particle=pi_serv->TrackIdToParticle_P(Track_mu_id);
		    int pdg=-1;
		    float end_eng=0;
		    float tr_stx=0;
		    float tr_sty=0;
		    float tr_stz=0;
		    float tr_enx=0;
		    float tr_eny=0;
		    float tr_enz=0;
		    std::string pri("primary");
		    bool isPrimary=0;
		    if(particle){
		       origin=pi_serv->TrackIdToMCTruth_P(Track_mu_id)->Origin();
		       pdg=particle->PdgCode();
		       isPrimary=particle->Process()==pri;
		       end_eng=particle->EndE()*1000;
		       tr_stx=particle->Vx();
		       tr_sty=particle->Vy();
		       tr_stz=particle->Vz();
		       tr_enx=particle->EndX();
		       tr_eny=particle->EndY();
		       tr_enz=particle->EndZ();
		    }
			
	           
		   float true_st_reco_st=TMath::Sqrt((tr_stx-rec_stx)*(tr_stx-rec_stx)+(tr_sty-rec_sty)*(tr_sty-rec_sty)+(tr_stz-rec_stz)*(tr_stz-rec_stz));
		   float true_en_reco_en=TMath::Sqrt((tr_enx-rec_enx)*(tr_enx-rec_enx)+(tr_eny-rec_eny)*(tr_eny-rec_eny)+(tr_enz-rec_enz)*(tr_enz-rec_enz));
		   float true_st_reco_en=TMath::Sqrt((tr_stx-rec_enx)*(tr_stx-rec_enx)+(tr_sty-rec_eny)*(tr_sty-rec_eny)+(tr_stz-rec_enz)*(tr_stz-rec_enz));
		   float true_en_reco_st=TMath::Sqrt((tr_enx-rec_stx)*(tr_enx-rec_stx)+(tr_eny-rec_sty)*(tr_eny-rec_sty)+(tr_enz-rec_stz)*(tr_enz-rec_stz));
		   
		   
		   if(origin==1 && ((true_st_reco_st<5 && true_en_reco_en<5)||(true_st_reco_en<5 && true_en_reco_st<5)) && isPrimary){
		      
		      /////////////////////// Stopping Muons ///////////////////////
		      
		      if(TMath::Abs(pdg)==13 && end_eng<107){
		         stop_mu++;
			 mu_trklen[stop_mu-1]=track.Length();
		         mu_truth_endE[stop_mu-1]=end_eng;
		         mu_start_info[stop_mu-1][0]=pos.X();
		         mu_start_info[stop_mu-1][1]=pos.Y();
		         mu_start_info[stop_mu-1][2]=pos.Z();
		         mu_end_info[stop_mu-1][0]=end.X();
		         mu_end_info[stop_mu-1][1]=end.Y();
		         mu_end_info[stop_mu-1][2]=end.Z();
		         mu_trackthetaxz[stop_mu-1]=theta_xz;
		         mu_trackthetayz[stop_mu-1]=theta_yz;
		         mu_TrkID[stop_mu-1]=track.ID();
		         mu_trkstartcosxyz[stop_mu-1][0]=dir_start.X();
		         mu_trkstartcosxyz[stop_mu-1][1]=dir_start.Y();
		         mu_trkstartcosxyz[stop_mu-1][2]=dir_start.Z();
		         mu_trkendcosxyz[stop_mu-1][0]=dir_end.X();
		         mu_trkendcosxyz[stop_mu-1][1]=dir_end.Y();
		         mu_trkendcosxyz[stop_mu-1][2]=dir_end.Z();
			 mu_truth_start_info[stop_mu-1][0]=tr_stx;
			 mu_truth_start_info[stop_mu-1][1]=tr_sty;
			 mu_truth_start_info[stop_mu-1][2]=tr_stz;
			 mu_truth_end_info[stop_mu-1][0]=tr_enx;
			 mu_truth_end_info[stop_mu-1][1]=tr_eny;
			 mu_truth_end_info[stop_mu-1][2]=tr_enz;
		         for(size_t ical = 0; ical<calos.size(); ++ical){
			     if(!calos[ical]) continue;
	                     if(!calos[ical]->PlaneID().isValid) continue;
	                     int planenum = calos[ical]->PlaneID().Plane;
	                     if(planenum<0||planenum>2) continue;
			     const size_t NHits = calos[ical] -> dEdx().size();
		             mu_ntrkhits[stop_mu-1][planenum]=int(NHits);
			     for(size_t iHit = 0; iHit < NHits; ++iHit){
			         const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
			         mu_trkdqdx[stop_mu-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
			         mu_trkdedx[stop_mu-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
			         mu_trkresrange[stop_mu-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
			         mu_trkhitx[stop_mu-1][planenum][iHit]=TrkPos.X();
		                 mu_trkhity[stop_mu-1][planenum][iHit]=TrkPos.Y();
		                 mu_trkhitz[stop_mu-1][planenum][iHit]=TrkPos.Z();
				 mu_trkpitch[stop_mu-1][planenum][iHit]=(calos[ical] -> TrkPitchVec())[iHit];
			     }
		         }
			   
		         if(fmpid.isValid()){
			    for(size_t ipid = 0; ipid < pids.size(); ++ipid){
			        if(!pids[ipid]->PlaneID().isValid) continue;
     	                        int planenum = pids[ipid]->PlaneID().Plane;
     	                        if(planenum<0||planenum>2) continue;
				mu_def_pida[stop_mu-1][planenum]=float(pids[ipid]->PIDA());
			    }
			 }
		       } // stopping muons...
		     
		       ////////////////////////// Stopping Pions ///////////////////////////
		     
		       if(pdg==211 && end_eng<140){
		          int pion_trk_id=particle->TrackId();
			  int n_daughters=0;
			  const sim::ParticleList& plist = pi_serv->ParticleList();
                          sim::ParticleList::const_iterator itPart = plist.begin(),pend = plist.end();
			  for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
			      const simb::MCParticle* pPart = (itPart++)->second;
			      if(pPart->Mother()==pion_trk_id){
			         n_daughters++;
			      }
			  }
			  if(n_daughters>0){
			  stop_pi++;
			  pi_trklen[stop_pi-1]=track.Length();
		          pi_truth_endE[stop_pi-1]=end_eng;
		          pi_start_info[stop_pi-1][0]=pos.X();
		          pi_start_info[stop_pi-1][1]=pos.Y();
		          pi_start_info[stop_pi-1][2]=pos.Z();
		          pi_end_info[stop_pi-1][0]=end.X();
		          pi_end_info[stop_pi-1][1]=end.Y();
		          pi_end_info[stop_pi-1][2]=end.Z();
		          pi_trackthetaxz[stop_pi-1]=theta_xz;
		          pi_trackthetayz[stop_pi-1]=theta_yz;
		          pi_TrkID[stop_pi-1]=track.ID();
		          pi_trkstartcosxyz[stop_pi-1][0]=dir_start.X();
		          pi_trkstartcosxyz[stop_pi-1][1]=dir_start.Y();
		          pi_trkstartcosxyz[stop_pi-1][2]=dir_start.Z();
		          pi_trkendcosxyz[stop_pi-1][0]=dir_end.X();
		          pi_trkendcosxyz[stop_pi-1][1]=dir_end.Y();
		          pi_trkendcosxyz[stop_pi-1][2]=dir_end.Z();
			  pi_truth_start_info[stop_pi-1][0]=tr_stx;
			  pi_truth_start_info[stop_pi-1][1]=tr_sty;
			  pi_truth_start_info[stop_pi-1][2]=tr_stz;
			  pi_truth_end_info[stop_pi-1][0]=tr_enx;
			  pi_truth_end_info[stop_pi-1][1]=tr_eny;
			  pi_truth_end_info[stop_pi-1][2]=tr_enz;
			  for(size_t ical = 0; ical<calos.size(); ++ical){
			     if(!calos[ical]) continue;
	                     if(!calos[ical]->PlaneID().isValid) continue;
	                     int planenum = calos[ical]->PlaneID().Plane;
	                     if(planenum<0||planenum>2) continue;
			     const size_t NHits = calos[ical] -> dEdx().size();
		             pi_ntrkhits[stop_pi-1][planenum]=int(NHits);
			     for(size_t iHit = 0; iHit < NHits; ++iHit){
			         const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
			         pi_trkdqdx[stop_pi-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
			         pi_trkdedx[stop_pi-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
			         pi_trkresrange[stop_pi-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
			         pi_trkhitx[stop_pi-1][planenum][iHit]=TrkPos.X();
		                 pi_trkhity[stop_pi-1][planenum][iHit]=TrkPos.Y();
		                 pi_trkhitz[stop_pi-1][planenum][iHit]=TrkPos.Z();
				 pi_trkpitch[stop_mu-1][planenum][iHit]=(calos[ical] -> TrkPitchVec())[iHit];
			     }
		         }
			 
			 if(fmpid.isValid()){
			    for(size_t ipid = 0; ipid < pids.size(); ++ipid){
			        if(!pids[ipid]->PlaneID().isValid) continue;
     	                        int planenum = pids[ipid]->PlaneID().Plane;
     	                        if(planenum<0||planenum>2) continue;
				pi_def_pida[stop_pi-1][planenum]=float(pids[ipid]->PIDA());
			    }
			 }
		      
		       } // decay pion
		      } // stopping pions.....
		     
		     //////////////////////////////////// Stopping Kaons ///////////////////////////////////////
		     
		     if(pdg==321 && end_eng<494){
		          stop_k++;
			  k_trklen[stop_k-1]=track.Length();
		          k_truth_endE[stop_k-1]=end_eng;
		          k_start_info[stop_k-1][0]=pos.X();
		          k_start_info[stop_k-1][1]=pos.Y();
		          k_start_info[stop_k-1][2]=pos.Z();
		          k_end_info[stop_k-1][0]=end.X();
		          k_end_info[stop_k-1][1]=end.Y();
		          k_end_info[stop_k-1][2]=end.Z();
		          k_trackthetaxz[stop_k-1]=theta_xz;
		          k_trackthetayz[stop_k-1]=theta_yz;
		          k_TrkID[stop_k-1]=track.ID();
		          k_trkstartcosxyz[stop_k-1][0]=dir_start.X();
		          k_trkstartcosxyz[stop_k-1][1]=dir_start.Y();
		          k_trkstartcosxyz[stop_k-1][2]=dir_start.Z();
		          k_trkendcosxyz[stop_k-1][0]=dir_end.X();
		          k_trkendcosxyz[stop_k-1][1]=dir_end.Y();
		          k_trkendcosxyz[stop_k-1][2]=dir_end.Z();
			  k_truth_start_info[stop_k-1][0]=tr_stx;
			  k_truth_start_info[stop_k-1][1]=tr_sty;
			  k_truth_start_info[stop_k-1][2]=tr_stz;
			  k_truth_end_info[stop_k-1][0]=tr_enx;
			  k_truth_end_info[stop_k-1][1]=tr_eny;
			  k_truth_end_info[stop_k-1][2]=tr_enz;
			  for(size_t ical = 0; ical<calos.size(); ++ical){
			     if(!calos[ical]) continue;
	                     if(!calos[ical]->PlaneID().isValid) continue;
	                     int planenum = calos[ical]->PlaneID().Plane;
	                     if(planenum<0||planenum>2) continue;
			     const size_t NHits = calos[ical] -> dEdx().size();
		             k_ntrkhits[stop_k-1][planenum]=int(NHits);
			     for(size_t iHit = 0; iHit < NHits; ++iHit){
			         const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
			         k_trkdqdx[stop_k-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
			         k_trkdedx[stop_k-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
			         k_trkresrange[stop_k-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
			         k_trkhitx[stop_k-1][planenum][iHit]=TrkPos.X();
		                 k_trkhity[stop_k-1][planenum][iHit]=TrkPos.Y();
		                 k_trkhitz[stop_k-1][planenum][iHit]=TrkPos.Z();
				 k_trkpitch[stop_mu-1][planenum][iHit]=(calos[ical] -> TrkPitchVec())[iHit];
			     }
		         }
			 
			 if(fmpid.isValid()){
			    for(size_t ipid = 0; ipid < pids.size(); ++ipid){
			        if(!pids[ipid]->PlaneID().isValid) continue;
     	                        int planenum = pids[ipid]->PlaneID().Plane;
     	                        if(planenum<0||planenum>2) continue;
				k_def_pida[stop_k-1][planenum]=float(pids[ipid]->PIDA());
			    }
			 }
		      } // stopping kaons.....
		     
		     //////////////////////////////////// Stopping Protons ///////////////////////////////////////
		     
		      if(pdg==2212 && end_eng<940){
		         int proton_trk_id=particle->TrackId();
			 int n_daughters=0;
			 const sim::ParticleList& plist = pi_serv->ParticleList();
                         sim::ParticleList::const_iterator itPart = plist.begin(),pend = plist.end();
			 for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
			     const simb::MCParticle* pPart = (itPart++)->second;
			     if(pPart->Mother()==proton_trk_id){
			        if(pPart->Process()!="conv"&&pPart->Process()!="LowEnConversion"&&pPart->Process()!="Pair"&&pPart->Process()!="compt"&&pPart->Process()!="Compt"&&pPart->Process()!="Brem"&&pPart->Process()!="phot"&&pPart->Process()!="Photo"&&pPart->Process()!="Ion"&&pPart->Process()!="annihil"){
			           n_daughters++;
				}
			     }
			 }
			 
			 if(n_daughters==0){
			    stop_p++;
		            p_trklen[stop_p-1]=track.Length();
			    p_truth_endE[stop_p-1]=end_eng;
		            p_start_info[stop_p-1][0]=pos.X();
		            p_start_info[stop_p-1][1]=pos.Y();
		            p_start_info[stop_p-1][2]=pos.Z();
		            p_end_info[stop_p-1][0]=end.X();
		            p_end_info[stop_p-1][1]=end.Y();
		            p_end_info[stop_p-1][2]=end.Z();
		            p_trackthetaxz[stop_p-1]=theta_xz;
		            p_trackthetayz[stop_p-1]=theta_yz;
		            p_TrkID[stop_p-1]=track.ID();
		            p_trkstartcosxyz[stop_p-1][0]=dir_start.X();
		            p_trkstartcosxyz[stop_p-1][1]=dir_start.Y();
		            p_trkstartcosxyz[stop_p-1][2]=dir_start.Z();
		            p_trkendcosxyz[stop_p-1][0]=dir_end.X();
		            p_trkendcosxyz[stop_p-1][1]=dir_end.Y();
		            p_trkendcosxyz[stop_p-1][2]=dir_end.Z();
			    p_truth_start_info[stop_p-1][0]=tr_stx;
			    p_truth_start_info[stop_p-1][1]=tr_sty;
			    p_truth_start_info[stop_p-1][2]=tr_stz;
			    p_truth_end_info[stop_p-1][0]=tr_enx;
			    p_truth_end_info[stop_p-1][1]=tr_eny;
			    p_truth_end_info[stop_p-1][2]=tr_enz;
			    for(size_t ical = 0; ical<calos.size(); ++ical){
			        if(!calos[ical]) continue;
	                        if(!calos[ical]->PlaneID().isValid) continue;
	                        int planenum = calos[ical]->PlaneID().Plane;
	                        if(planenum<0||planenum>2) continue;
			        const size_t NHits = calos[ical] -> dEdx().size();
		                p_ntrkhits[stop_p-1][planenum]=int(NHits);
			        for(size_t iHit = 0; iHit < NHits; ++iHit){
			            const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
			            p_trkdqdx[stop_p-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
			            p_trkdedx[stop_p-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
			            p_trkresrange[stop_p-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
			            p_trkhitx[stop_p-1][planenum][iHit]=TrkPos.X();
		                    p_trkhity[stop_p-1][planenum][iHit]=TrkPos.Y();
		                    p_trkhitz[stop_p-1][planenum][iHit]=TrkPos.Z();
				    p_trkpitch[stop_mu-1][planenum][iHit]=(calos[ical] -> TrkPitchVec())[iHit];
			        }
		            }
			 
			    if(fmpid.isValid()){
			       for(size_t ipid = 0; ipid < pids.size(); ++ipid){
			           if(!pids[ipid]->PlaneID().isValid) continue;
     	                           int planenum = pids[ipid]->PlaneID().Plane;
     	                           if(planenum<0||planenum>2) continue;
			           p_def_pida[stop_p-1][planenum]=float(pids[ipid]->PIDA());
			       }
			   }
		      
		         } // number of daughters are zero
		      } // stopping protons
		    
		    
		    } // origin 1 and true and end points match....
		  } // reco trk contained....
		} // loop over reco trks...
	     } // reco track info saved...
	  } // is MC
          fEventTree->Fill();
} // end of analyze function
	   
 /////////////////// Defintion of reset function ///////////
void XYZvalidatioin::reset(){
     run=-9999;
     subrun=-9999;
     event=-9999;
     stop_mu=-9999;
     stop_pi=-9999;
     stop_k=-9999;
     stop_p=-9999;
     for(int i=0; i<kMaxTracks; i++){
         mu_trklen[i]=-9999;
	 mu_trackthetaxz[i]=-9999;
	 mu_trackthetayz[i]=-9999;
	 mu_TrkID[i]=-9999;
	 mu_truth_endE[i]=-9999;
	 mu_start_info[i][0]=-9999;
	 mu_start_info[i][1]=-9999;
	 mu_start_info[i][2]=-9999;
	 mu_end_info[i][0]=-9999;
	 mu_end_info[i][1]=-9999;
	 mu_end_info[i][2]=-9999;
	 mu_trkstartcosxyz[i][0]=-9999;
	 mu_trkstartcosxyz[i][1]=-9999;
	 mu_trkstartcosxyz[i][2]=-9999; 
	 mu_trkendcosxyz[i][0]=-9999;
	 mu_trkendcosxyz[i][1]=-9999;
	 mu_trkendcosxyz[i][2]=-9999;
	 mu_ntrkhits[i][0]=-9999;
	 mu_ntrkhits[i][1]=-9999;
	 mu_ntrkhits[i][2]=-9999;
	 mu_def_pida[i][0]=-9999;
	 mu_def_pida[i][1]=-9999;
	 mu_def_pida[i][2]=-9999;
	 mu_truth_start_info[i][0]=-9999;
	 mu_truth_start_info[i][1]=-9999;
	 mu_truth_start_info[i][2]=-9999;
	 mu_truth_end_info[i][0]=-9999;
	 mu_truth_end_info[i][1]=-9999;
	 mu_truth_end_info[i][2]=-9999;
	 pi_trklen[i]=-9999;
	 pi_trackthetaxz[i]=-9999;
	 pi_trackthetayz[i]=-9999;
	 pi_TrkID[i]=-9999;
	 pi_truth_endE[i]=-9999;
	 pi_start_info[i][0]=-9999;
	 pi_start_info[i][1]=-9999;
	 pi_start_info[i][2]=-9999;
	 pi_end_info[i][0]=-9999;
	 pi_end_info[i][1]=-9999;
	 pi_end_info[i][2]=-9999;
	 pi_trkstartcosxyz[i][0]=-9999;
	 pi_trkstartcosxyz[i][1]=-9999;
	 pi_trkstartcosxyz[i][2]=-9999; 
	 pi_trkendcosxyz[i][0]=-9999;
	 pi_trkendcosxyz[i][1]=-9999;
	 pi_trkendcosxyz[i][2]=-9999;
	 pi_ntrkhits[i][0]=-9999;
	 pi_ntrkhits[i][1]=-9999;
	 pi_ntrkhits[i][2]=-9999;
	 pi_def_pida[i][0]=-9999;
	 pi_def_pida[i][1]=-9999;
	 pi_def_pida[i][2]=-9999;
	 pi_truth_start_info[i][0]=-9999;
	 pi_truth_start_info[i][1]=-9999;
	 pi_truth_start_info[i][2]=-9999;
	 pi_truth_end_info[i][0]=-9999;
	 pi_truth_end_info[i][1]=-9999;
	 pi_truth_end_info[i][2]=-9999;
	 k_trklen[i]=-9999;
	 k_trackthetaxz[i]=-9999;
	 k_trackthetayz[i]=-9999;
	 k_TrkID[i]=-9999;
	 k_truth_endE[i]=-9999;
	 k_start_info[i][0]=-9999;
	 k_start_info[i][1]=-9999;
	 k_start_info[i][2]=-9999;
	 k_end_info[i][0]=-9999;
	 k_end_info[i][1]=-9999;
	 k_end_info[i][2]=-9999;
	 k_trkstartcosxyz[i][0]=-9999;
	 k_trkstartcosxyz[i][1]=-9999;
	 k_trkstartcosxyz[i][2]=-9999; 
	 k_trkendcosxyz[i][0]=-9999;
	 k_trkendcosxyz[i][1]=-9999;
	 k_trkendcosxyz[i][2]=-9999;
	 k_ntrkhits[i][0]=-9999;
	 k_ntrkhits[i][1]=-9999;
	 k_ntrkhits[i][2]=-9999;
	 k_def_pida[i][0]=-9999;
	 k_def_pida[i][1]=-9999;
	 k_def_pida[i][2]=-9999;
	 k_truth_start_info[i][0]=-9999;
	 k_truth_start_info[i][1]=-9999;
	 k_truth_start_info[i][2]=-9999;
	 k_truth_end_info[i][0]=-9999;
	 k_truth_end_info[i][1]=-9999;
	 k_truth_end_info[i][2]=-9999;
	 p_trklen[i]=-9999;
	 p_trackthetaxz[i]=-9999;
	 p_trackthetayz[i]=-9999;
	 p_TrkID[i]=-9999;
	 p_truth_endE[i]=-9999;
	 p_start_info[i][0]=-9999;
	 p_start_info[i][1]=-9999;
	 p_start_info[i][2]=-9999;
	 p_end_info[i][0]=-9999;
	 p_end_info[i][1]=-9999;
	 p_end_info[i][2]=-9999;
	 p_trkstartcosxyz[i][0]=-9999;
	 p_trkstartcosxyz[i][1]=-9999;
	 p_trkstartcosxyz[i][2]=-9999; 
	 p_trkendcosxyz[i][0]=-9999;
	 p_trkendcosxyz[i][1]=-9999;
	 p_trkendcosxyz[i][2]=-9999;
	 p_ntrkhits[i][0]=-9999;
	 p_ntrkhits[i][1]=-9999;
	 p_ntrkhits[i][2]=-9999;
	 p_def_pida[i][0]=-9999;
	 p_def_pida[i][1]=-9999;
	 p_def_pida[i][2]=-9999;
	 p_truth_start_info[i][0]=-9999;
	 p_truth_start_info[i][1]=-9999;
	 p_truth_start_info[i][2]=-9999;
	 p_truth_end_info[i][0]=-9999;
	 p_truth_end_info[i][1]=-9999;
	 p_truth_end_info[i][2]=-9999;
	 for(int j=0; j<3; j++){
	     for(int k=0; k<3000; k++){
	         p_trkpitch[i][j][k]=-9999;
		 k_trkpitch[i][j][k]=-9999;
		 pi_trkpitch[i][j][k]=-9999;
		 mu_trkpitch[i][j][k]=-9999;
		 mu_trkdqdx[i][j][k]=-9999;
		 mu_trkdedx[i][j][k]=-9999;
		 mu_trkresrange[i][j][k]=-9999;
		 mu_trkhitx[i][j][k]=-9999;
		 mu_trkhity[i][j][k]=-9999;
		 mu_trkhitz[i][j][k]=-9999;
		 pi_trkdqdx[i][j][k]=-9999;
		 pi_trkdedx[i][j][k]=-9999;
		 pi_trkresrange[i][j][k]=-9999;
		 pi_trkhitx[i][j][k]=-9999;
		 pi_trkhity[i][j][k]=-9999;
		 pi_trkhitz[i][j][k]=-9999;
		 k_trkdqdx[i][j][k]=-9999;
		 k_trkdedx[i][j][k]=-9999;
		 k_trkresrange[i][j][k]=-9999;
		 k_trkhitx[i][j][k]=-9999;
		 k_trkhity[i][j][k]=-9999;
		 k_trkhitz[i][j][k]=-9999;
		 p_trkdqdx[i][j][k]=-9999;
		 p_trkdedx[i][j][k]=-9999;
		 p_trkresrange[i][j][k]=-9999;
		 p_trkhitx[i][j][k]=-9999;
		 p_trkhity[i][j][k]=-9999;
		 p_trkhitz[i][j][k]=-9999;
	     }
	 }
     }
}
 //////////////////////// End of definition ///////////////	
	  
DEFINE_ART_MODULE(XYZvalidatioin)
}



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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

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

const int kMaxTracks  = 1000;
const int nbin = 49; //split the total SCE map into 48 bins -- 49th is overflow.

using namespace std;

namespace microboone{
  
  class XYZSCEcorrectionCrossingTrack : public art::EDAnalyzer {
  public:
    
    explicit XYZSCEcorrectionCrossingTrack(fhicl::ParameterSet const& pset);
    virtual ~XYZSCEcorrectionCrossingTrack();
    
    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    int getBin(float ypos);
    
  private:
    TTree* fEventTree;
    TTree* fACTree;

    TH1D     *Efield_hist;
    TH2D     *Efield_vs_x;

    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    cross_trks;
    Int_t    hasCalo[kMaxTracks][3]={0,0,0};
    Int_t    caloIsValid[kMaxTracks][3]={0,0,0};
    Int_t    planenum_exists[kMaxTracks][3]={0,0,0};
    Int_t    planenumber;
    Int_t    all_trks;
    Float_t  xprojectedlen[kMaxTracks];
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Float_t  trackthetaxy[kMaxTracks];
    Int_t    TrkID[kMaxTracks]; 
    Float_t  trkstartcosxyz[kMaxTracks][3];
    Float_t  trkendcosxyz[kMaxTracks][3];
    Float_t  trackstartpos[kMaxTracks][3];
    Float_t  trackendpos[kMaxTracks][3];
    Int_t    ntrkhits[kMaxTracks][3];
    Int_t    ntrajpoints[kMaxTracks];
    Double_t  trajposx[kMaxTracks][3][3000];
    Double_t  trajposy[kMaxTracks][3][3000];
    Double_t  trajposz[kMaxTracks][3][3000];
    Float_t  trkdqdx[kMaxTracks][3][3000];
    Float_t  trkdedx[kMaxTracks][3][3000];
    Float_t  trkresrange[kMaxTracks][3][3000];
    Float_t  trkhitx[kMaxTracks][3][3000];
    Float_t  trkhity[kMaxTracks][3][3000];
    Float_t  trkhitz[kMaxTracks][3][3000];
    Float_t  trkdqdxSCE[kMaxTracks][3][3000];
    Float_t  trkdedxSCE[kMaxTracks][3][3000];
    Float_t  trkresrangeSCE[kMaxTracks][3][3000];
    Float_t  trkhitxSCE[kMaxTracks][3][3000];
    Float_t  trkhitySCE[kMaxTracks][3][3000];
    Float_t  trkhitzSCE[kMaxTracks][3][3000];
    Float_t  trkhitxoffset[kMaxTracks][3][3000];
    Float_t  trkhityoffset[kMaxTracks][3][3000];
    Float_t  trkhitzoffset[kMaxTracks][3][3000];
    Float_t  trkdqdx_efieldcorr[kMaxTracks][3][3000];
    Float_t  trkdedx_efieldcorr[kMaxTracks][3][3000];
    Float_t  trkresrange_efieldcorr[kMaxTracks][3][3000];
    Float_t  trkRconstant[kMaxTracks][3][3000];
    Float_t  trkRcorr[kMaxTracks][3][3000];
    Float_t  trkhitEx[kMaxTracks][3][3000];
    Float_t  trkhitEy[kMaxTracks][3][3000];
    Float_t  trkhitEz[kMaxTracks][3][3000];
    Float_t  trkhitpitch[kMaxTracks][3][3000];
    Float_t  trkhitpitchSCE[kMaxTracks][3][3000];
    Float_t  E0[kMaxTracks][3][3000];
    Float_t  EfieldCorr[kMaxTracks][3][3000];
    Float_t  trklocaltime[kMaxTracks][3][3000];
    Float_t  trklocaltimeSCE[kMaxTracks][3][3000];
    Float_t  trklocalDriftVel[kMaxTracks][3][3000];
    Float_t  trklocalDriftVelSCE[kMaxTracks][3][3000];
    Float_t  trk_hit_time[kMaxTracks][3][3000];
    Float_t  trk_hit_stime[kMaxTracks][3][3000];
    Float_t  trk_hit_etime[kMaxTracks][3][3000];
    Float_t  trk_hit_rms[kMaxTracks][3][3000];
    Float_t  trk_hit_amplitude[kMaxTracks][3][3000];
    Float_t  trk_hit_summedADC[kMaxTracks][3][3000];
    Float_t  trk_hit_integral[kMaxTracks][3][3000];
    Float_t  trk_hit_gof[kMaxTracks][3][3000];
    Int_t    trk_hit_multiplicity[kMaxTracks][3][3000];
    Int_t    trk_hit_dof[kMaxTracks][3][3000];
    Int_t    trk_hit_localindex[kMaxTracks][3][3000];
    Float_t  trk_hit_sptvx[kMaxTracks][3][3000];
    Float_t  trk_hit_sptvy[kMaxTracks][3][3000];
    Float_t  trk_hit_sptvz[kMaxTracks][3][3000];

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fCalorimetrySCEModuleLabel;
    std::string fParticleIDModuleLabel;
    std::string fSpacePointModuleLabel;
    std::string fClusterModuleLabel;
    std::string fVertexModuleLabel;

    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fSaveClusterInfo;
    bool  fSaveClusterHitInfo;
    bool  fSaveGenieInfo;
    bool  fSaveVertexInfo; 
    bool  ft0_corr; 
    float fG4minE;

    // modified box model parameters for data
    double fModBoxA;
    double fModBoxB;

  }; 
  
  //========================================================================
  XYZSCEcorrectionCrossingTrack::XYZSCEcorrectionCrossingTrack(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel","")        ),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
    fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel","")     ),
    fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","")     ),  
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel","")          ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fCalorimetrySCEModuleLabel   (pset.get< std::string >("CalorimetrySCEModuleLabel","")  ),
    fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel","")   ),
    fSpacePointModuleLabel    (pset.get<std::string>("SpacePointModuleLabel")        ),
    fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel","")      ),
    fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel","")       ),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),  
    fSaveClusterInfo          (pset.get< bool>("SaveClusterInfo",false)),
    fSaveClusterHitInfo       (pset.get< bool>("SaveClusterHitInfo",false)),
    fSaveGenieInfo            (pset.get< bool>("SaveGenieInfo",false)),
    fSaveVertexInfo           (pset.get< bool>("SaveVertexInfo",false)),
    ft0_corr                  (pset.get< bool >("t0_corr","false") ),
    fG4minE                   (pset.get< float>("G4minE",0.01)),
    fModBoxA               (pset.get< double >("ModBoxA")),
    fModBoxB               (pset.get< double >("ModBoxB")) 
  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
  
  //========================================================================
  XYZSCEcorrectionCrossingTrack::~XYZSCEcorrectionCrossingTrack(){
  }
  //========================================================================
  
  //========================================================================
  void XYZSCEcorrectionCrossingTrack::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  
    Efield_hist = tfs->make<TH1D>("Efield_offsets","Efield offsets",20,-0.05,0.05);  
    Efield_vs_x = tfs->make<TH2D>("Efield_vs_x","Efield offsets vs x",26,0.,260,20,-0.05,0.05);  
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fEventTree->Branch("cross_trks",&cross_trks,"cross_trks/I");
    fEventTree->Branch("hasCalo",&hasCalo,"hasCalo[cross_trks][3]/I");
    fEventTree->Branch("caloIsValid",&caloIsValid,"caloIsValid[cross_trks][3]/I");
    fEventTree->Branch("planenum_exists",&planenum_exists,"planenum_exists[cross_trks][3]/I");
    fEventTree->Branch("planenumber",&planenumber,"planenumber/I");
    fEventTree->Branch("all_trks",&all_trks,"all_trks/I");
    fEventTree->Branch("xprojectedlen",xprojectedlen,"xprojectedlen[all_trks]/F");
    fEventTree->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[cross_trks]/F");
    fEventTree->Branch("trackthetayz",trackthetayz,"trackthetayz[cross_trks]/F");
    fEventTree->Branch("trackthetaxy",trackthetaxy,"trackthetaxy[cross_trks]/F");
    fEventTree->Branch("trackstartpos",trackstartpos,"trackstartpos[cross_trks][3]/F");
    fEventTree->Branch("trackendpos",trackendpos,"trackendpos[cross_trks][3]/F");
    fEventTree->Branch("TrkID",TrkID,"TrkID[cross_trks]/I");
    fEventTree->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[cross_trks][3]/F");
    fEventTree->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[cross_trks][3]/F");
    fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[cross_trks][3]/I");
    fEventTree->Branch("ntrajpoints",ntrajpoints,"ntrajpoints[cross_trks]/I");
    fEventTree->Branch("trkdqdx",trkdqdx,"trkdqdx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdedx",trkdedx,"trkdedx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkresrange",trkresrange,"trkresrange[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdqdxSCE",trkdqdxSCE,"trkdqdxSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdedxSCE",trkdedxSCE,"trkdedxSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trkresrangeSCE",trkresrangeSCE,"trkresrangeSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitx",trkhitx,"trkhitx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhity",trkhity,"trkhity[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitz",trkhitz,"trkhitz[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitxSCE",trkhitxSCE,"trkhitxSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitySCE",trkhitySCE,"trkhitySCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitzSCE",trkhitzSCE,"trkhitzSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitEx",trkhitEx,"trkhitEx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitEy",trkhitEy,"trkhitEy[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitEz",trkhitEz,"trkhitEz[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitpitch",trkhitpitch,"trkhitpitch[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitpitchSCE",trkhitpitchSCE,"trkhitpitchSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("E0",E0,"E0[cross_trks][3][3000]/F");
    fEventTree->Branch("EfieldCorr",EfieldCorr,"EfieldCorr[cross_trks][3][3000]/F");
    fEventTree->Branch("trajposx",trajposx,"trajposx[cross_trks][3][3000]/F");
    fEventTree->Branch("trajposy",trajposy,"trajposy[cross_trks][3][3000]/F");
    fEventTree->Branch("trajposz",trajposz,"trajposz[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdqdx_efieldcorr",trkdqdx_efieldcorr,"trkdqdx_efieldcorr[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdedx_efieldcorr",trkdedx_efieldcorr,"trkdedx_efieldcorr[cross_trks][3][3000]/F");
    fEventTree->Branch("trkresrange_efieldcorr",trkresrange_efieldcorr,"trkresrange_efieldcorr[cross_trks][3][3000]/F");
    fEventTree->Branch("trkRconstant",trkRconstant,"trkRconstant[cross_trks][3][3000]/F");
    fEventTree->Branch("trkRcorr",trkRcorr,"trkRcorr[cross_trks][3][3000]/F");
    fEventTree->Branch("trklocaltime",trklocaltime,"trklocaltime[cross_trks][3][3000]/F");
    fEventTree->Branch("trklocaltimeSCE",trklocaltimeSCE,"trklocaltimeSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trklocalDriftVel",trklocalDriftVel,"trklocalDriftVel[cross_trks][3][3000]/F");
    fEventTree->Branch("trklocalDriftVelSCE",trklocalDriftVelSCE,"trklocalDriftVelSCE[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_time",trk_hit_time,"trk_hit_time[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_stime",trk_hit_stime,"trk_hit_stime[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_etime",trk_hit_etime,"trk_hit_etime[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_rms",trk_hit_rms,"trk_hit_rms[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_amplitude",trk_hit_amplitude,"trk_hit_amplitude[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_summedADC",trk_hit_summedADC,"trk_hit_summedADC[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_integral",trk_hit_integral,"trk_hit_integral[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_gof",trk_hit_gof,"trk_hit_gof[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_multiplicity",trk_hit_multiplicity,"trk_hit_multiplicity[cross_trks][3][3000]/I");
    fEventTree->Branch("trk_hit_dof",trk_hit_dof,"trk_hit_dof[cross_trks][3][3000]/I");
    fEventTree->Branch("trk_hit_localindex",trk_hit_localindex,"trk_hit_localindex[cross_trks][3][3000]/I");
    fEventTree->Branch("trk_hit_sptvx",trk_hit_sptvx,"trk_hit_sptvx[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_sptvy",trk_hit_sptvy,"trk_hit_sptvy[cross_trks][3][3000]/F");
    fEventTree->Branch("trk_hit_sptvz",trk_hit_sptvz,"trk_hit_sptvz[cross_trks][3][3000]/F");
  
  //========================================================================
  void XYZSCEcorrectionCrossingTrack::endJob(){     
    
  }
  
  //========================================================================
  void XYZSCEcorrectionCrossingTrack::beginRun(const art::Run&){
    mf::LogInfo("XYZSCEcorrectionCrossingTrack")<<"begin run..."<<std::endl;
  }
  //========================================================================
  
  //========================================================================
  
  //========================================================================

  void XYZSCEcorrectionCrossingTrack::analyze( const art::Event& evt){
    reset();
  
    //space charge service
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
 
    //Efield
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    //Clock
    auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();
    
    art::Handle< std::vector<recob::Track> > trackListHandle;
    art::Handle< std::vector<recob::Track> > trackListHandle_2;
    art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
    
    std::vector<art::Ptr<recob::Track> > tracklist;
    std::vector<art::Ptr<recob::Track> > tracklist_2;
    std::vector<art::Ptr<recob::PFParticle> > pfplist;

    // Get Geometry
    art::ServiceHandle<geo::Geometry const> geom;
    
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    //if(evt.getByLabel("pandoraCosmic",trackListHandle_2)) art::fill_ptr_vector(tracklist_2, trackListHandle_2);
    //if(evt.getByLabel("pandoraCosmic",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    art::FindManyP<anab::Calorimetry> fmcalSCE(trackListHandle, evt, fCalorimetrySCEModuleLabel);
    art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); //this has more information about hit-track association, only available in PMA for now
    //art::FindManyP<anab::T0> trk_t0_assn_v(trackListHandle_2, evt, "pandoraCosmicT0Reco");
    //art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle_2, evt, "pandoraCosmic");
    //art::FindManyP<recob::Track> pft_trk_2_assn(PFPListHandle, evt, "pandoraCosmicKalmanTrack");
    
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
    

    //std::cout << "##### XDrift Velocity = " << detprop->DriftVelocity() << std::endl;
    size_t NTracks = tracklist.size();

    //loop over tracks
    for(size_t i=0; i<NTracks;++i){
       art::Ptr<recob::Track> ptrack(trackListHandle, i);
       std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i), calosSCE=fmcalSCE.at(i);
       const recob::Track& track = *ptrack;
       const auto& pos = track.Vertex();
       const auto& dir_start = track.VertexDirection();
       const auto& dir_end   = track.EndDirection();
       const auto& end = track.End();

       double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
       double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
       double theta_xy = std::atan2(dir_start.X(), dir_start.Y());

       float X = std::abs(pos.X()-end.X());
       all_trks++;
       if(all_trks<=kMaxTracks) xprojectedlen[all_trks-1]=X;
         
       if(X>250. && X<270.){
	 cross_trks++;
         auto traj = track.Trajectory();
         auto NPoints = traj.NumberTrajectoryPoints();
         ntrajpoints[cross_trks-1]=NPoints; 
         for( size_t iPoint = 0; iPoint < NPoints; ++iPoint ){
           auto trajPos = traj.LocationAtPoint(iPoint);
           trajposx[cross_trks-1][0][iPoint]=trajPos.X();
           trajposy[cross_trks-1][1][iPoint]=trajPos.Y();
           trajposz[cross_trks-1][2][iPoint]=trajPos.Z();
         }

	 trackthetaxz[cross_trks-1]=theta_xz;
	 trackthetayz[cross_trks-1]=theta_yz;
	 trackthetaxy[cross_trks-1]=theta_xy;
	 TrkID[cross_trks-1]=track.ID();
	 trkstartcosxyz[cross_trks-1][0]=dir_start.X();
	 trkstartcosxyz[cross_trks-1][1]=dir_start.Y();
	 trkstartcosxyz[cross_trks-1][2]=dir_start.Z();
	 trkendcosxyz[cross_trks-1][0]=dir_end.X();
	 trkendcosxyz[cross_trks-1][1]=dir_end.Y();
	 trkendcosxyz[cross_trks-1][2]=dir_end.Z();
	 
         trackstartpos[cross_trks-1][0]=pos.X();
         trackstartpos[cross_trks-1][1]=pos.Y();
         trackstartpos[cross_trks-1][2]=pos.Z();
         trackendpos[cross_trks-1][0]=end.X();
         trackendpos[cross_trks-1][1]=end.Y();
         trackendpos[cross_trks-1][2]=end.Z();
         
	 ////////////// Getting min X //////////////////
	 
	 double minx_p2 = 1e10;
	 double minx_p1 = 1e10;
	 double minx_p0 = 1e10;
	 double minx_p2SCE = 1e10;
	 double minx_p1SCE = 1e10;
	 double minx_p0SCE = 1e10;
	 for(size_t ical = 0; ical<calos.size(); ++ical){
	   
	   if(!calos[ical]) continue;
	   hasCalo[cross_trks-1][ical]++; 
	   
	   if(!calos[ical]->PlaneID().isValid) continue;
	   caloIsValid[cross_trks-1][ical]++;
           
	   int planenum = calos[ical]->PlaneID().Plane;
	   planenumber=planenum;
	   if(planenum<0||planenum>2) continue;
	   planenum_exists[cross_trks-1][ical]++;
	   if(planenum==2){
	     const size_t NHits = calos[ical] -> dEdx().size();
	     for(size_t iHit = 0; iHit < NHits; ++iHit){
	       const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	       const auto& TrkPosSCE = (calosSCE[ical] -> XYZ())[iHit];
               
	       if(TrkPos.X()>-180){
		 if(TrkPos.X()<minx_p2) minx_p2=TrkPos.X();
	       }
	       if(TrkPosSCE.X()>-180){
		 if(TrkPosSCE.X()<minx_p2SCE) minx_p2SCE=TrkPosSCE.X();
	       }
	     }
	   }
	   //////////////////// end p2 ////////////////////////////
	   if(planenum == 1){
	     const size_t NHits = calos[ical] -> dEdx().size();
	     for(size_t iHit = 0; iHit < NHits; ++iHit){
	       const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	       const auto& TrkPosSCE = (calosSCE[ical] -> XYZ())[iHit];
	       if(TrkPos.X()>-180){
		 if(TrkPos.X()<minx_p1) minx_p1=TrkPos.X();
	       }
	       if(TrkPosSCE.X()>-180){
		 if(TrkPosSCE.X()<minx_p1SCE) minx_p1SCE=TrkPosSCE.X();
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
	       const auto& TrkPosSCE = (calosSCE[ical] -> XYZ())[iHit];
	       if(TrkPosSCE.X()>-180){
		 if(TrkPosSCE.X()<minx_p0SCE) minx_p0SCE=TrkPosSCE.X();
	       }
	     }
	   }
	   ///////////////// end p0 ////////////////////////////
	 }
	 
	 ////////////////// End of min X //////////////
	
         //get hits from tracks
         //std::cout << "getting hits" << std::endl;
         std::vector<std::vector<unsigned int>> hits(calos.size());
         std::vector<art::Ptr<recob::Hit>> allHits = fmht.at(i);
         art::FindManyP<recob::SpacePoint> fmspts(allHits, evt, fSpacePointModuleLabel);
         for (size_t ah = 0; ah < allHits.size(); ++ah) {
           hits[allHits[ah]->WireID().Plane].push_back(ah); 
         }

	 for(size_t ical = 0; ical<calos.size(); ++ical){
           //std::cout << "ical = 0: " << calos[0] -> dEdx().size() << std::endl;
           //std::cout << "ical = 1: " << calos[1] -> dEdx().size() << std::endl;
           //std::cout << "ical = 2: " << calos[2] -> dEdx().size() << std::endl;
	   if(!calos[ical]) continue;
	   if(!calos[ical]->PlaneID().isValid) continue;
	   int planenum = calos[ical]->PlaneID().Plane;
	   if(planenum<0||planenum>2) continue;
	   const size_t NHits = calos[ical] -> dEdx().size();
           //std::cout << "NHits, hits["<<ical<<"] = " << NHits << ", " << hits[planenum].size() << std::endl;
           //std::cout << "allHits size, hits size = " << allHits.size() << ", " << hits.size() << std::endl;
           for(size_t pl=0; pl < hits.size(); pl++ ) std::cout << "hits["<<pl<<"] = " << hits[pl].size() << std::endl;
	   ntrkhits[cross_trks-1][planenum]=int(NHits);
	   for(size_t iHit = 0; iHit < NHits; ++iHit){
	     const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	     const auto& TrkPosSCE = (calosSCE[ical] -> XYZ())[iHit];
	     //const auto& TrkTPindices = (calosSCE[ical] -> TpIndices())[iHit];
             //std::cout << "tpindices iHit =  " << iHit << ", " << TrkTPindices << std::endl;
	     double minx=0;
	     double minxSCE=0;
	     if(planenum==0){
               minx=minx_p0;
               minxSCE=minx_p0SCE;
             }
	     if(planenum==1){ 
               minx=minx_p1;
               minxSCE=minx_p1SCE;
             }
	     if(planenum==2){ 
               minx=minx_p2;
               minxSCE=minx_p2SCE;
             }
             //std::cout << "minx, minxSCE = " << minx << ", " << minxSCE << std::endl;
             if(ft0_corr) minx=0.0;
             if(ft0_corr) minxSCE=0.0;
             // Look at the space charge correction
             // geo::Vector_t sce_poscorr    = SCE->GetCalPosOffsets(geo::Point_t(TrkPos.X()-minx,TrkPos.Y(),TrkPos.Z()));
             // Look at the electric field map correction
             //std::cout << " is t0 corr? no SCE? SCE =  " << ft0_corr << ", " << minx << ", " << minxSCE << std::endl;
             //std::cout << " x, y, z orig =  " << TrkPos.X()-minx << ", " << TrkPos.Y() << ", " << TrkPos.Z() << std::endl;
             //std::cout << " x, y, z SCE =  " << TrkPosSCE.X()-minxSCE << ", " << TrkPosSCE.Y() << ", " << TrkPosSCE.Z() << std::endl;

             geo::Vector_t vEfieldcorr = SCE->GetCalEfieldOffsets(geo::Point_t(TrkPosSCE.X()-minxSCE, TrkPosSCE.Y(), TrkPosSCE.Z()));

             //calculate the recombination factor at default Efield
             double dEdx_mip = 2.1;                     // dE/dx in Mev/cm^2  
             double rho = detprop->Density();           // LAr density in g/cm^3
             double E_field = detprop->Efield();        // Electric Field in the drift region in KV/cm
             //get the Efield offsets in x,y,z
             double Ex_corr = vEfieldcorr.X();
             double Ey_corr = vEfieldcorr.Y();
             double Ez_corr = vEfieldcorr.Z();

             //std::cout << std::setprecision(5) << "Ex_corr, Ey_corr, Ez_corr = " << Ex_corr << ", " << Ey_corr << ", " << Ez_corr << std::endl;
             double E_field_corr = sqrt( ((E_field + E_field*Ex_corr)*(E_field + E_field*Ex_corr))+((E_field*Ey_corr)*(E_field*Ey_corr))+((E_field*Ez_corr)*(E_field*Ez_corr)) ); // Electric Field in the drift region in KV/cm
             double Beta_const = fModBoxB / (rho * E_field);
             double Beta_corr  = fModBoxB / (rho * E_field_corr);
             double Alpha      = fModBoxA;
             double R_const = (log( Alpha + Beta_const*dEdx_mip))/(Beta_const*dEdx_mip);
             //at corrected Efield
             double R_corr = (log( Alpha + Beta_corr*dEdx_mip))/(Beta_corr*dEdx_mip);
             //drift velocity
             double temp   = detprop -> Temperature();
             double DriftVelocity_0 = detprop->DriftVelocity(E_field,temp);
             double DriftVelocity_corr = detprop-> DriftVelocity(E_field_corr,temp);
             double time = detprop->ConvertXToTicks(TrkPos.X()-minx,calos[ical]->PlaneID());
             double timeSCE = detprop->ConvertXToTicks(TrkPosSCE.X()-minxSCE,calos[ical]->PlaneID());

             //get hit time
             //match hit index

             //std::cout << "XYZ calo: (" << TrkPos.X()-minx << ", " << TrkPos.Y() << ", " << TrkPos.Z() << " )" << std::endl;
             //std::cout << "XYZ caloSCE: (" << TrkPosSCE.X()-minxSCE << ", " << TrkPosSCE.Y() << ", " << TrkPosSCE.Z() << " )" << std::endl;
            
             for (size_t i = 0; i < hits[planenum].size(); ++i) {
               //Get space points associated with the hit
               if( allHits[hits[planenum][i]].key() == calos[ical]->TpIndices()[iHit] ){
                 //std::cout << "found matching hit!" << std::endl;
                 //std::cout << "allHits[hits["<<planenum<<"]["<<i<<"]].key(), calos["<<ical<<"] ->TpIndices)["<<iHit<<"] = " << allHits[hits[planenum][i]].key() << ", " << calos[ical]->TpIndices()[iHit] << std::endl;
                 //std::cout << "hits["<<planenum<<"]["<<i<<"] = " << hits[planenum][i] << ", " << std::endl; 
                 trk_hit_time[cross_trks-1][planenum][iHit] = detclock->TPCTick2TrigTime(allHits[hits[planenum][i]]->PeakTime());
                 trk_hit_stime[cross_trks-1][planenum][iHit] = detclock->TPCTick2TrigTime(allHits[hits[planenum][i]]->PeakTimeMinusRMS());
                 trk_hit_etime[cross_trks-1][planenum][iHit] = detclock->TPCTick2TrigTime(allHits[hits[planenum][i]]->PeakTimePlusRMS());
                 trk_hit_rms[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->RMS();
                 trk_hit_amplitude[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->PeakAmplitude();
                 trk_hit_summedADC[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->SummedADC();
                 trk_hit_integral[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->Integral();
                 trk_hit_gof[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->GoodnessOfFit();
                 trk_hit_multiplicity[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->Multiplicity();
                 trk_hit_dof[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->DegreesOfFreedom();
                 trk_hit_localindex[cross_trks-1][planenum][iHit] = allHits[hits[planenum][i]]->LocalIndex();
                 //std::cout << "time, stime, etime = " << trk_hit_time[cross_trks-1][planenum][iHit] << ", " << trk_hit_stime[cross_trks-1][planenum][iHit] << ", " << trk_hit_etime[cross_trks-1][planenum][iHit] << std::endl;
                 std::vector<art::Ptr<recob::SpacePoint>> sptv = fmspts.at(hits[planenum][i]);
                 //std:: cout << "XYZ sptv = (" << sptv[0]->XYZ()[0] << ", " << sptv[0]->XYZ()[1] << ", " << sptv[0]->XYZ()[2] << ", " << std::endl;
                 trk_hit_sptvx[cross_trks-1][planenum][iHit] = sptv[0]->XYZ()[0];
                 trk_hit_sptvy[cross_trks-1][planenum][iHit] = sptv[0]->XYZ()[1];
                 trk_hit_sptvz[cross_trks-1][planenum][iHit] = sptv[0]->XYZ()[2];

               }
             }

             //std::cout << std::setprecision(5) << "Efield, Efield_corr = " << E_field << ", " << E_field_corr << std::endl;
             if(planenum==2)Efield_hist->Fill( -(E_field - (E_field + E_field*Ex_corr)));
             if(planenum==2)Efield_vs_x->Fill( (TrkPos.X()-minx),-(E_field - (E_field + E_field*Ex_corr)));
             //std::cout << std::setprecision(5) << "R_const, R_corr = " << R_const << ", " << R_corr << std::endl;
             //these are all the branches after taking into account the space charge correction
	     trkhitx[cross_trks-1][planenum][iHit]=TrkPos.X()-minx;
	     trkhity[cross_trks-1][planenum][iHit]=TrkPos.Y();
	     trkhitz[cross_trks-1][planenum][iHit]=TrkPos.Z();
	     trkhitxSCE[cross_trks-1][planenum][iHit]=TrkPosSCE.X()-minxSCE;
	     trkhitySCE[cross_trks-1][planenum][iHit]=TrkPosSCE.Y();
	     trkhitzSCE[cross_trks-1][planenum][iHit]=TrkPosSCE.Z();

             //Efield correction
	     trkhitEx[cross_trks-1][planenum][iHit]=E_field + E_field*Ex_corr;
	     trkhitEy[cross_trks-1][planenum][iHit]=E_field*Ey_corr;
	     trkhitEz[cross_trks-1][planenum][iHit]=E_field*Ez_corr;
             E0[cross_trks-1][planenum][iHit]=E_field;
             EfieldCorr[cross_trks-1][planenum][iHit]=E_field_corr;

	     trkdqdx[cross_trks-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	     trkdedx[cross_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	     trkresrange[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	     trkhitpitch[cross_trks-1][planenum][iHit]=(calos[ical]->TrkPitchVec())[iHit];
	     trkdqdxSCE[cross_trks-1][planenum][iHit]=(calosSCE[ical] -> dQdx())[iHit];
	     trkdedxSCE[cross_trks-1][planenum][iHit]=(calosSCE[ical] -> dEdx())[iHit];
	     trkresrangeSCE[cross_trks-1][planenum][iHit]=(calosSCE[ical]->ResidualRange())[iHit];
	     trkhitpitchSCE[cross_trks-1][planenum][iHit]=(calosSCE[ical]->TrkPitchVec())[iHit];

	      //------------------------------- 
	      // Time and Drift Velocity
	      //-------------------------------
	      trklocaltime[cross_trks-1][planenum][iHit]=detclock->TPCTick2TrigTime(time);
	      trklocaltimeSCE[cross_trks-1][planenum][iHit]=detclock->TPCTick2TrigTime(timeSCE);
	      trklocalDriftVel[cross_trks-1][planenum][iHit]=DriftVelocity_0;
	      trklocalDriftVelSCE[cross_trks-1][planenum][iHit]=DriftVelocity_corr;

             //assuming dQ/dx has been space corrected, the dQ/dx from the calorimetry should have the correct X,Y,Z position:
             double dQdx_space_efieldcorr = (calosSCE[ical] -> dQdx())[iHit];
             double dQdx_EfieldCorr = (R_const/R_corr) * dQdx_space_efieldcorr;
	     trkdqdx_efieldcorr[cross_trks-1][planenum][iHit]=dQdx_EfieldCorr;
	     trkdedx_efieldcorr[cross_trks-1][planenum][iHit]=(calosSCE[ical] -> dEdx())[iHit];
	     trkresrange_efieldcorr[cross_trks-1][planenum][iHit]=(calosSCE[ical]->ResidualRange())[iHit];
	     trkRconstant[cross_trks-1][planenum][iHit]=R_const;
	     trkRcorr[cross_trks-1][planenum][iHit]=R_corr;
             std::cout << "dQdx, dQdxSCE, dQdxEfieldCorr = " << trkdqdx[cross_trks-1][planenum][iHit] << ", " << trkdqdxSCE[cross_trks-1][planenum][iHit] << ", " << trkdqdx_efieldcorr[cross_trks-1][planenum][iHit] << std::endl;
	   } // loop over iHit..
	 } // loop over ical2 nd time...
      } // crossing trks...
    } // loop over trks..     
    
    fEventTree->Fill();
  } // end of analyze function
  
  /////////////////// Defintion of reset function ///////////
  void XYZSCEcorrectionCrossingTrack::reset(){
    run = -9999;
    subrun = -9999;
    event = -9999;
    evttime = -9999;
    cross_trks = -9999;
    all_trks = -9999;
    year_month_date=-9999;
    hour_min_sec=-9999;
    //AC_trks=-9999;
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
      ntrkhits[i][0] = -9999;
      ntrkhits[i][1] = -9999;
      ntrkhits[i][2] = -9999;
      ntrajpoints[i] = -9999;
	for(int k=0; k<3000; k++){
	  trkdqdx[i][j][k]=-9999;
	  trkdedx[i][j][k]=-9999;
	  trkresrange[i][j][k]=-9999;
	  trkhitpitch[i][j][k]=-9999;
	  trkhitx[i][j][k]=-9999;
	  trkhity[i][j][k]=-9999;
	  trkhitz[i][j][k]=-9999;
          trkdqdxSCE[i][j][k]=-9999.;
          trkdedxSCE[i][j][k]=-9999.;
          trkresrangeSCE[i][j][k]=-9999.;
          trkhitpitchSCE[i][j][k]=-9999.;
          trkhitxSCE[i][j][k]=-9999.;
          trkhitySCE[i][j][k]=-9999.;
          trkhitzSCE[i][j][k]=-9999.;
          trkhitxoffset[i][j][k]=-9999.;
          trkhityoffset[i][j][k]=-9999.;
          trkhitzoffset[i][j][k]=-9999.;
          trkdqdx_efieldcorr[i][j][k]=-9999.;
          trkdedx_efieldcorr[i][j][k]=-9999.;
          trkresrange_efieldcorr[i][j][k]=-9999.;
          trkRconstant[i][j][k]=-9999.;
          trkRcorr[i][j][k]=-9999.;
          trkhitEx[i][j][k]=-9999.;
          trkhitEy[i][j][k]=-9999.;
          trkhitEz[i][j][k]=-9999.;
          E0[i][j][k]=-9999.;
          EfieldCorr[i][j][k]=-9999.;
	  trklocaltime[i][j][k]=-9999.;
	  trklocaltimeSCE[i][j][k]=-9999.;
	  trklocalDriftVel[i][j][k]=-9999.;
	  trklocalDriftVelSCE[i][j][k]=-9999.;
          trk_hit_time[i][j][k]=-9999.;
          trk_hit_stime[i][j][k]=-9999.;
          trk_hit_etime[i][j][k]=-9999.;
          trk_hit_rms[i][j][k]=-9999.;
          trk_hit_amplitude[i][j][k]=-9999.;
          trk_hit_summedADC[i][j][k]=-9999.;
          trk_hit_integral[i][j][k]=-9999.;
          trk_hit_gof[i][j][k]=-9999.;
          trk_hit_multiplicity[i][j][k]=-9999.;
          trk_hit_dof[i][j][k]=-9999.;
          trk_hit_localindex[i][j][k]=-9999.;
          trk_hit_sptvx[i][j][k]=-9999.;
          trk_hit_sptvy[i][j][k]=-9999.;
          trk_hit_sptvz[i][j][k]=-9999.;
	}
      }
    }
    
  }

  int XYZSCEcorrectionCrossingTrack::getBin(float ypos){

  int bin=49;
  for( int i = 0; i < 48; i++ ){
    float ylow = -120. + 5*i;
    float yhi  = -120. + 5*(i+1);
    bin = i;
    if( ypos > ylow && ypos < yhi )
      break;
  }// close loop

  std::cout << "bin = " << bin << std::endl;
  return bin;

  }//end getBin
  //////////////////////// End of definition ///////////////	
  
  DEFINE_ART_MODULE(XYZSCEcorrectionCrossingTrack)
}

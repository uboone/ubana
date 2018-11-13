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
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
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
    TTree* fACTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Double_t evttime; 
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    cross_trks;
    Int_t    hasCalo[3]={0,0,0};
    Int_t    caloIsValid[3]={0,0,0};
    Int_t    planenum_exists[3]={0,0,0};
    Int_t    planenumber;
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
    //Efield corrected
    Float_t  trkdqdx_efieldcorr[kMaxTracks][3][3000];
    Float_t  trkdedx_efieldcorr[kMaxTracks][3][3000];
    Float_t  trkresrange_efieldcorr[kMaxTracks][3][3000];
    Float_t  trkRconstant[kMaxTracks][3][3000];
    Float_t  trkRcorr[kMaxTracks][3][3000];
    
    Int_t    AC_trks;
    Float_t  AC_trklen[kMaxTracks];
    Int_t    AC_TrkID[kMaxTracks];
    Float_t  AC_xprojectedlen[kMaxTracks];
    Float_t  AC_trackthetaxz[kMaxTracks];
    Float_t  AC_trackthetayz[kMaxTracks];
    Float_t  AC_start[kMaxTracks][3];
    Float_t  AC_end[kMaxTracks][3];
    Float_t  AC_trkstartcosxyz[kMaxTracks][3];
    Float_t  AC_trkendcosxyz[kMaxTracks][3];
    Int_t    AC_ntrkhits[kMaxTracks][3];
    Float_t  AC_trkdqdx[kMaxTracks][3][3000];
    Float_t  AC_trkdedx[kMaxTracks][3][3000];
    Float_t  AC_trkresrange[kMaxTracks][3][3000];
    Float_t  AC_trkhitx[kMaxTracks][3][3000];
    Float_t  AC_trkhity[kMaxTracks][3][3000];
    Float_t  AC_trkhitz[kMaxTracks][3][3000];
    Float_t  AC_t0[kMaxTracks];
    
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

    // modified box model parameters for data
    double fModBoxA;
    double fModBoxB;

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
    fG4minE                   (pset.get< float>("G4minE",0.01)),
    fModBoxA               (pset.get< double >("ModBoxA")),
    fModBoxB               (pset.get< double >("ModBoxB")) 
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
    fEventTree->Branch("hasCalo",&hasCalo,"hasCalo[3]/I");
    fEventTree->Branch("caloIsValid",&caloIsValid,"caloIsValid[3]/I");
    fEventTree->Branch("planenum_exists",&planenum_exists,"planenum_exists[3]/I");
    fEventTree->Branch("planenumber",&planenumber,"planenumber/I");
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
    
    fACTree = tfs->make<TTree>("ACpierce", "Event Tree from AC pierce");
    fACTree->Branch("event", &event,"event/I");
    fACTree->Branch("evttime",&evttime,"evttime/D");
    fACTree->Branch("run", &run,"run/I");
    fACTree->Branch("subrun", &subrun,"surbrun/I");
    fACTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fACTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fACTree->Branch("AC_trks",&AC_trks,"AC_trks/I");
    fACTree->Branch("AC_trklen",&AC_trklen,"AC_trklen[AC_trks]/F");
    fACTree->Branch("AC_TrkID",&AC_TrkID,"AC_TrkID[AC_trks]/I");
    fACTree->Branch("AC_xprojectedlen",AC_xprojectedlen,"AC_xprojectedlen[AC_trks]/F");
    fACTree->Branch("AC_trackthetaxz",AC_trackthetaxz,"AC_trackthetaxz[AC_trks]/F");
    fACTree->Branch("AC_trackthetayz",AC_trackthetayz,"AC_trackthetayz[AC_trks]/F");
    fACTree->Branch("AC_start",AC_start,"AC_start[AC_trks][3]/F");
    fACTree->Branch("AC_end",AC_end,"AC_end[AC_trks][3]/F");
    fACTree->Branch("AC_trkstartcosxyz",AC_trkstartcosxyz,"AC_trkstartcosxyz[AC_trks][3]/F");
    fACTree->Branch("AC_trkendcosxyz",AC_trkendcosxyz,"AC_trkendcosxyz[AC_trks][3]/F");
    fACTree->Branch("AC_ntrkhits",AC_ntrkhits,"AC_ntrkhits[AC_trks][3]/I");
    fACTree->Branch("AC_trkdqdx",AC_trkdqdx,"AC_trkdqdx[AC_trks][3][3000]/F");
    fACTree->Branch("AC_trkdedx",AC_trkdedx,"AC_trkdedx[AC_trks][3][3000]/F");
    fACTree->Branch("AC_trkresrange",AC_trkresrange,"AC_trkresrange[AC_trks][3][3000]/F");
    fACTree->Branch("AC_trkhitx",AC_trkhitx,"AC_trkhitx[AC_trks][3][3000]/F");
    fACTree->Branch("AC_trkhity",AC_trkhity,"AC_trkhity[AC_trks][3][3000]/F");
    fACTree->Branch("AC_trkhitz",AC_trkhitz,"AC_trkhitz[AC_trks][3][3000]/F");
    fACTree->Branch("AC_t0",&AC_t0,"AC_t0[AC_trks]/F");
  }
  
  //========================================================================
  void XYZcorrection::endJob(){     
    
  }
  
  //========================================================================
  void XYZcorrection::beginRun(const art::Run&){
    mf::LogInfo("XYZcorrection")<<"begin run..."<<std::endl;
    mf::LogInfo("XYZcorrection")<<"I am running this locally..."<<std::endl;
  }
  //========================================================================
  
  //========================================================================
  
  //========================================================================

  void XYZcorrection::analyze( const art::Event& evt){
    reset();
  
    //space charge service
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
 
    //Efield
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    art::Handle< std::vector<recob::Track> > trackListHandle;
    art::Handle< std::vector<recob::Track> > trackListHandle_2;
    art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
    
    std::vector<art::Ptr<recob::Track> > tracklist;
    std::vector<art::Ptr<recob::Track> > tracklist_2;
    std::vector<art::Ptr<recob::PFParticle> > pfplist;
    
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    if(evt.getByLabel("pandoraCosmic",trackListHandle_2)) art::fill_ptr_vector(tracklist_2, trackListHandle_2);
    if(evt.getByLabel("pandoraCosmic",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
    
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    art::FindManyP<anab::T0> trk_t0_assn_v(trackListHandle_2, evt, "pandoraCosmicT0Reco");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle_2, evt, "pandoraCosmic");
    art::FindManyP<recob::Track> pft_trk_2_assn(PFPListHandle, evt, "pandoraCosmicKalmanTrack");
    
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
       if(all_trks<=kMaxTracks) xprojectedlen[all_trks-1]=X;
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
	   hasCalo[ical]++; 
	   
	   if(!calos[ical]->PlaneID().isValid) continue;
	   caloIsValid[ical]++;
           
	   int planenum = calos[ical]->PlaneID().Plane;
	   planenumber=planenum;
	   if(planenum<0||planenum>2) continue;
	   planenum_exists[ical]++;
	   if(planenum==2){
	     const size_t NHits = calos[ical] -> dEdx().size();
	     for(size_t iHit = 0; iHit < NHits; ++iHit){
	       const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	       if(TrkPos.X()>-180){
		 if(TrkPos.X()<minx_p2) minx_p2=TrkPos.X();
                 //std::cout << "TrkPos.X(), minx_p2 = " << TrkPos.X() << ", " << minx_p2 << std::endl; 
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
                 //std::cout << "TrkPos.X(), minx_p1 = " << TrkPos.X() << ", " << minx_p1 << std::endl; 
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
                 std::cout << "TrkPos.X(), minx_p0 = " << TrkPos.X() << ", " << minx_p0 << std::endl; 
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
             // Look at the space charge correction
             // geo::Vector_t sce_poscorr    = SCE->GetCalPosOffsets(geo::Point_t(TrkPos.X()-minx,TrkPos.Y(),TrkPos.Z()));
             // Look at the electric field map correction
             geo::Vector_t vEfieldcorr = SCE->GetCalEfieldOffsets(geo::Point_t(TrkPos.X()-minx, TrkPos.Y(), TrkPos.Z()));
             //cout << "voxel points x, y, z = " << TrkPos.X()-minx << ", " << TrkPos.Y() << ", " << TrkPos.Z() << endl;
             //cout << "vEfieldcorr.X(), vEfieldcorr.Y(), vEfieldcorr.Z() = " << vEfieldcorr.X() << ", " << vEfieldcorr.Y() << ", " << vEfieldcorr.Z() << endl; 

             //calculate the recombination factor at default Efield
             double dEdx_mip = 2.1;                     // dE/dx in Mev/cm^2  
             double rho = detprop->Density();           // LAr density in g/cm^3
             double E_field = detprop->Efield();        // Electric Field in the drift region in KV/cm
             //get the Efield offsets in x,y,z
             double Ex_corr = vEfieldcorr.X();
             double Ey_corr = vEfieldcorr.Y();
             double Ez_corr = vEfieldcorr.Z();
             double E_field_corr = sqrt( ((E_field + E_field*Ex_corr)*(E_field + E_field*Ex_corr))+(Ey_corr*Ey_corr)+(Ez_corr*Ez_corr) );        // Electric Field in the drift region in KV/cm
             double Beta_const = fModBoxB / (rho * E_field);
             double Beta_corr  = fModBoxB / (rho * E_field_corr);
             double Alpha      = fModBoxA;
             double R_const = (log( Alpha + Beta_const*dEdx_mip))/(Beta_const*dEdx_mip);
             //at corrected Efield
             double R_corr = (log( Alpha + Beta_corr*dEdx_mip))/(Beta_corr*dEdx_mip);
             std::cout << "Efield, Efield_corr = " << E_field << ", " << E_field_corr << std::endl;
             std::cout << "R_const, R_corr = " << R_const << ", " << R_corr << std::endl;
             //these are all the branches after taking into account the space charge correction
	     trkhitx[cross_trks-1][planenum][iHit]=TrkPos.X()-minx;
	     trkhity[cross_trks-1][planenum][iHit]=TrkPos.Y();
	     trkhitz[cross_trks-1][planenum][iHit]=TrkPos.Z();

	     trkdqdx[cross_trks-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	     trkdedx[cross_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	     trkresrange[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];

             //assuming dQ/dx has been space corrected, the dQ/dx from the calorimetry should have the correct X,Y,Z position:
             double dQdx_space_efieldcorr = (calos[ical] -> dQdx())[iHit];
             double dQdx_EfieldCorr = R_const/R_corr * dQdx_space_efieldcorr;
	     trkdqdx_efieldcorr[cross_trks-1][planenum][iHit]=dQdx_EfieldCorr;
	     trkdedx_efieldcorr[cross_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	     trkresrange_efieldcorr[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	     trkRconstant[cross_trks-1][planenum][iHit]=R_const;
	     trkRcorr[cross_trks-1][planenum][iHit]=R_corr;

	   } // loop over iHit..
	 } // loop over ical2 nd time...
       } // crossing trks...
    } // loop over trks..     
    
    fEventTree->Fill();
    
    AC_trks=0;
    size_t N2Tracks = tracklist_2.size();
    double efield = detprop -> Efield();
    double temp   = detprop -> Temperature();
    double fDriftVelocity = detprop -> DriftVelocity(efield,temp);
    int    nACPT_tracks = 0; 
    
    for(size_t i=0; i<N2Tracks;++i){
      std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(i);
      if(!t0s.size()) continue;
      std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);
      if(!pfps.size()) continue;
      std::vector<art::Ptr<recob::Track>> trks=pft_trk_2_assn.at(pfps[0].key());
      if(!trks.size()) continue;
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(trks[0].key());
      if  (t0s.size() == 1) {	
	nACPT_tracks++;
	auto t0 = t0s.at(0);	
	double t_zero=t0->Time();
	
	std::cout << "The t0 value = " << t_zero << " us." << std::endl;
	
	AC_trks++;
	//art::Ptr<recob::Track> ptrack(trackListHandle,i);
	art::Ptr<recob::Track> ptrack = trks[0];
	const recob::Track& track = *ptrack;
	TVector3 pos, dir_start, dir_end, end;
	pos = track.Vertex();
	dir_start = track.VertexDirection();
	dir_end   = track.EndDirection();
	end = track.End();
	
	double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
	double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
	float X = std::abs(pos.X()-end.X());
	AC_trklen[AC_trks-1]=track.Length();
	AC_t0[AC_trks-1]=t_zero;
	AC_TrkID[AC_trks-1]=track.ID();
	AC_xprojectedlen[AC_trks-1]=X;
	AC_trackthetaxz[AC_trks-1]=theta_xz;
	AC_trackthetayz[AC_trks-1]=theta_yz;
	AC_start[AC_trks-1][0]=pos.X(); AC_start[AC_trks-1][1]=pos.Y(); AC_start[AC_trks-1][2]=pos.Z();
	AC_end[AC_trks-1][0]=end.X(); AC_end[AC_trks-1][1]=end.Y(); AC_end[AC_trks-1][2]=end.Z();
	AC_trkstartcosxyz[AC_trks-1][0]=dir_start.X(); AC_trkstartcosxyz[AC_trks-1][1]=dir_start.Y(); AC_trkstartcosxyz[AC_trks-1][2]=dir_start.Z();
	AC_trkendcosxyz[AC_trks-1][0]=dir_end.X(); AC_trkendcosxyz[AC_trks-1][1]=dir_end.Y(); AC_trkendcosxyz[AC_trks-1][2]=dir_end.Z();
	for(size_t ical = 0; ical<calos.size(); ++ical){
	  if(!calos[ical]) continue;
	  if(!calos[ical]->PlaneID().isValid) continue;
	  int planenum = calos[ical]->PlaneID().Plane;
	  if(planenum<0||planenum>2) continue;
	  const size_t NHits = calos[ical] -> dEdx().size();
	  AC_ntrkhits[AC_trks-1][planenum]=int(NHits);
	  for(size_t iHit = 0; iHit < NHits; ++iHit){
	    const auto& TrkPos=(calos[ical] -> XYZ())[iHit];
	    AC_trkdqdx[AC_trks-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	    AC_trkdedx[AC_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	    AC_trkresrange[AC_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	    AC_trkhitx[AC_trks-1][planenum][iHit]=TrkPos.X()-(t_zero*fDriftVelocity);
	    AC_trkhity[AC_trks-1][planenum][iHit]=TrkPos.Y();
	    AC_trkhitz[AC_trks-1][planenum][iHit]=TrkPos.Z();
	  }
	}
      }
      
    } // loop over pandoracosmic trks
    fACTree->Fill();
    std::cout << "The number of ACPT tracks in this event = " <<  nACPT_tracks << "." << std::endl;
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
    AC_trks=-9999;
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
      AC_trklen[i]=-9999;
      AC_TrkID[i]=-9999;
      AC_xprojectedlen[i]=-9999;
      AC_trackthetaxz[i]=-9999;
      AC_trackthetayz[i]=-9999;
      AC_start[i][0]=-9999;
      AC_start[i][1]=-9999;
      AC_start[i][2]=-9999;
      AC_end[i][0]=-9999;
      AC_end[i][1]=-9999;
      AC_end[i][2]=-9999;
      AC_trkstartcosxyz[i][0]=-9999;
      AC_trkstartcosxyz[i][1]=-9999;
      AC_trkstartcosxyz[i][2]=-9999;
      AC_trkendcosxyz[i][0]=-9999;
      AC_trkendcosxyz[i][1]=-9999;
      AC_trkendcosxyz[i][2]=-9999;
      AC_ntrkhits[i][0]=-9999;
      AC_ntrkhits[i][1]=-9999;
      AC_ntrkhits[i][2]=-9999;
      AC_t0[i]=-9999;
      for(int j=0; j<3; j++){
	for(int k=0; k<3000; k++){
	  trkdqdx[i][j][k]=-9999;
	  trkdedx[i][j][k]=-9999;
	  trkresrange[i][j][k]=-9999;
	  trkhitx[i][j][k]=-9999;
	  trkhity[i][j][k]=-9999;
	  trkhitz[i][j][k]=-9999;
	  AC_trkdqdx[i][j][k]=-9999;
	  AC_trkdedx[i][j][k]=-9999;
	  AC_trkresrange[i][j][k]=-9999;
	  AC_trkhitx[i][j][k]=-9999;
	  AC_trkhity[i][j][k]=-9999;
	  AC_trkhitz[i][j][k]=-9999;
	}
      }
    }
    
    //std::cout << "The number of ACPT tracks in this event = " <<  nACPT_tracks << "." << std::endl;  
  }
  //////////////////////// End of definition ///////////////	
  
  DEFINE_ART_MODULE(XYZcorrection)
}

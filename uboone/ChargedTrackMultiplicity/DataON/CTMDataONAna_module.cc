/**
 * @file    CTMDataONAna_module.cc
 * @brief   Module for charged particle multiplicity analysis
 * @Author: Aleena Rafique (aleena@ksu.edu)
 * 
 ******************************************************************************/

  #ifndef CTMDATAONANA_H
  #define CTMDATAONANA_H
  
  // Framework includes
#include "cetlib_except/exception.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
  
  // LArSoft Includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"  
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
//#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubooneobj/UbooneOpticalFilter.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"

  
  // C++ Includes
  #include <iostream>
  #include <cstring>
  #include <sys/stat.h>
  #include <vector>
  #include "TMatrixD.h"
  #include "TVector3.h"
  
  // root
  #include "TH1D.h"
  #include <TMath.h>
  #include "TVectorD.h"
  #include "TFile.h"
  #include "TGeoManager.h"
  #include "TF1.h"
  #include "TF2.h"
  #include "TGraph.h"
  #include "TMath.h"
  #include "TTree.h"
  #include "TH2F.h"
  #include "TProfile.h"
  #include <TStyle.h>
  #include <TCanvas.h>
  
  //NuMuccinclusive includes
  //#include "uboone/TPCNeutrinoIDFilter/Algorithms/NuMuCCInclusiveAlg.h"
  #include "uboone/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgBase.h"
    //constexpr int kMaxGpar = 50000;
    constexpr int kMaxRecoallpar = 50000;
    constexpr int kMaxRecoselpar = 50000;
    //constexpr int kMaxtruematchpar = 50000;
    //constexpr int kMaxSeconpar = 500000;
    //constexpr int kMaxtruepar = 900000;
    //constexpr int kMaxtrueInAccppar = 50000;
    //constexpr int kMaxtruelongInAccppar = 50000;
  namespace CTMDataONAna {  
   
    class CTMDataONAna : public art::EDAnalyzer
    {  

    public:
 
     explicit CTMDataONAna(fhicl::ParameterSet const &pset);
     virtual ~CTMDataONAna();                        
     
      // The actual method for analyzing the events
     void analyze(art::Event const& evt) ;
     //void beginSubRun(const art::SubRun& sr);
     
     // Allow for fhicl parameters to possibly change during processing...
     void reconfigure(fhicl::ParameterSet const&)  ;
     
     // Called when job begins for definitions of histograms/tuples/etc.   
     void beginJob()  ;
     
     // Called when job completes to deal with output of stuff from beginJob
     void endJob();
    
     // Recover information from the start of a run (if processing across runs)
     void beginRun(const art::Run&);
     
    /// Resize the data strutcure for GENIE particles
    void ResizeGENIE(int nGENIE);

    /// Resize the data strutcure for GEANT particles
    void ResizeGEANT(int nGEANT);

    /// Resize the data strutcure for RECO particles
    void ResizeRECO(int nRECO);

    /// Clear all fields if this object (not the tracker algorithm data)
    void ClearLocalData();
    
    /// Clear all fields
    void Clear();

    /// Fills a container with begin()/end() interface
    template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { std::fill(std::begin(data), std::end(data), value); }

     //backTracker
     void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);

     //Function for calculation KE for diffenert particles and ranges
     float T(int pdg, float range);  

    // Function for calculating flash track distance 
    double  FlashTrackDist(double flash, double start, double end) const;
    
    //Function to see if in FV
    bool inFV(double x, double y, double z) const;
    
    float sx=0,sy=0,sz=0;
    float smin=0;

    TTree* fDataTreeTotal;
    TTree* fDataTreeAllEvts;
    TTree* fDataTreeFlashTag;
    TTree* fDataTreeVrtxinFV;
    TTree* fDataTreeTrkNearVrtx;
    TTree* fDataTreeFlashTrackMatch;
    TTree* fDataTreeCandMuContained;
    TTree* fDataTreeCCSel;         
    TTree* fDataTreePassBadReg;
    TTree* fDataTreePassCosmicReq;
    TTree* fDataTreePassLongTrackVrtxDis;
    TTree* fDataTreePassLongcolhits;     
    TTree* fDataTreeFinalSel;

    int   _frun;
    int   _fsubrun;
    int   _fevent; 
    float _fflashPE;
    float _fflashTime;
    float _fvrtxxPos;
    float _fvrtxyPos;
    float _fvrtxzPos;
    float _fflashtrackmatch;
    float _ftrkVrtxDis;
    float _ftrkcandlength;
    float _fvrtxx;
    float _fvrtxy;
    float _fvrtxz;
    int   _flongtrackID;
    int   _flongtrackflipped;
    float _flongtrackstartx;
    float _flongtrackstarty;	
    float _flongtrackstartz;
    float _flongtrackendx;
    float _flongtrackendy;
    float _flongtrackendz;
    float _flongTrackVrtxDis;
    float _flongtrackTheta;
    float _flongtrackCosTheta;
    float _flongtrackSinTheta;
    float _flongtrackPhi;
    float _flongtrackCosPhi;
    float _flongtrackSinPhi;
    float _flongtrackLength;
    int   _fbackward_track;
    int   _fnlongcolhits;
    float _flongtrackmcsfwdmom;   // track momentum assuming forward direction
    float _flongtrackmcsfwdll;    // MCS likelihood assuming forward direction
    float _flongtrackmcsfwderr;   // MCS uncertainty assuming forward direction
    float _flongtrackmcsbwdmom;   // track momentum assuming backward direction
    float _flongtrackmcsbwdll;    // MCS likelihood assuming backward direction
    float _flongtrackmcsbwderr;   // MCS uncertainty assuming backward direction
    int   _fbadReg_removed;
    int   _fstarty_cut;
    int   _fbrokenTracks;
        
    int _fNRecoallPart;
    float _falltrackVrtxDis[kMaxRecoallpar];
    float _falltrackID[kMaxRecoallpar];
    float _falltrackStartx[kMaxRecoallpar];
    float _falltrackStarty[kMaxRecoallpar];
    float _falltrackStartz[kMaxRecoallpar];
    float _falltrackEndx[kMaxRecoallpar];
    float _falltrackEndy[kMaxRecoallpar];
    float _falltrackEndz[kMaxRecoallpar];
    float _falltrackLength[kMaxRecoallpar];
    float _falltrackTheta[kMaxRecoallpar];
    float _falltrackCosTheta[kMaxRecoallpar];
    float _falltrackSinTheta[kMaxRecoallpar];
    float _falltrackPhi[kMaxRecoallpar];
    float _falltrackCosPhi[kMaxRecoallpar];
    float _falltrackSinPhi[kMaxRecoallpar];    

    int   _fNrecomult;
    float _fseltrackVrtxDis[kMaxRecoselpar];
    float _fselntrackhits[kMaxRecoselpar];
    float _fseltrackID[kMaxRecoselpar];
    float _fseltrackStartx[kMaxRecoselpar];
    float _fseltrackStarty[kMaxRecoselpar];
    float _fseltrackStartz[kMaxRecoselpar];
    float _fseltrackEndx[kMaxRecoselpar];
    float _fseltrackEndy[kMaxRecoselpar];
    float _fseltrackEndz[kMaxRecoselpar];
    float _fseltrackLength[kMaxRecoselpar];
    float _fseltrackTheta[kMaxRecoselpar];
    float _fseltrackCosTheta[kMaxRecoselpar];
    float _fseltrackSinTheta[kMaxRecoselpar];
    float _fseltrackPhi[kMaxRecoselpar];
    float _fseltrackCosPhi[kMaxRecoselpar];
    float _fseltrackSinPhi[kMaxRecoselpar];    

    int   _fPH;
    int   _fMCS;
    int   _fMCS1;
    int   _fMCS2;

    float _fPHratio;    
    float _fMCSratio;
    float _fMCSratio1;
    float _fMCSdiff;
    
    bool   fDoHists;
    
    TH1D*   fNFlashPerEvent;
    TH1D*   fFlashPE;
    TH1D*   fFlashTime;

        	
      private:
  
    const int Minhits = 15; //count all tracks with these minimum number of hits
    const int Disvrtx = 3;  //count tracks within this 3D distance around vertex
    const int MinLongTrackHits = 80; //long track minimum collection plane hits
    const int tau  = 18500000; //lifetime constant (converted 18.5msec to 18500000nsec)
    const float wire_sideband  = 0.9; //dead region wire sideband (0.9cm=3 wires)
           
    std::string              fHitsModuleLabel ;
    std::string              fTrackModuleLabel;
    std::string              fVertexModuleLabel;
    std::string              fGenieGenModuleLabel;
    std::string              fOpFlashModuleLabel;
    std::string              fTrackMCSFitLabel;

    double fDistToEdgeX;             ///< fiducial volume - x
    double fDistToEdgeY;             ///< fiducial volume - y
    double fDistToEdgeZ;             ///< fiducial volume - z

    double fFlashWidth;  
    double fBeamMin;
    double fBeamMax; 
    double fPEThresh; 
    double fMinTrk2VtxDist; 
    double fMinTrackLen; 
    /// @{
    /**
     *  @brief Standard useful properties
     */
    geo::GeometryCore const*            fGeometry;           ///< pointer to the Geometry service
    detinfo::DetectorProperties const*  fDetector;           ///< Pointer to the detector properties
    /// @}
    
   };
  
   //-----------------------------------------------------------------------
    // Constructor
   CTMDataONAna::CTMDataONAna(fhicl::ParameterSet const& pset): 
   EDAnalyzer(pset),
   fGeometry(lar::providerFrom<geo::Geometry>()),
   fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
      this->reconfigure(pset);
  }
  
    //-----------------------------------------------------------------------
    // Destructor
    CTMDataONAna::~CTMDataONAna() 
   {
    }
   //-----------------------------------------------------------------------
   void CTMDataONAna::ClearLocalData()
   {

  std::fill(_falltrackVrtxDis, _falltrackVrtxDis + sizeof(_falltrackVrtxDis)/sizeof(_falltrackVrtxDis[0]), -99999.); 
  std::fill(_falltrackID, _falltrackID + sizeof(_falltrackID)/sizeof(_falltrackID[0]), -99999.); 
  std::fill(_falltrackStartx, _falltrackStartx + sizeof(_falltrackStartx)/sizeof(_falltrackStartx[0]), -99999.); 
  std::fill(_falltrackStarty, _falltrackStarty + sizeof(_falltrackStarty)/sizeof(_falltrackStarty[0]), -99999.); 
  std::fill(_falltrackStartz, _falltrackStartz + sizeof(_falltrackStartz)/sizeof(_falltrackStartz[0]), -99999.); 
  std::fill(_falltrackEndx, _falltrackEndx + sizeof(_falltrackEndx)/sizeof(_falltrackEndx[0]), -99999.);
  std::fill(_falltrackEndy, _falltrackEndy + sizeof(_falltrackEndy)/sizeof(_falltrackEndy[0]), -99999.); 
  std::fill(_falltrackEndz, _falltrackEndz + sizeof(_falltrackEndz)/sizeof(_falltrackEndz[0]), -99999.); 
  std::fill(_falltrackLength, _falltrackLength + sizeof(_falltrackLength)/sizeof(_falltrackLength[0]), -99999.); 
  std::fill(_falltrackTheta, _falltrackTheta + sizeof(_falltrackTheta)/sizeof(_falltrackTheta[0]), -99999.); 
  std::fill(_falltrackCosTheta, _falltrackCosTheta + sizeof(_falltrackCosTheta)/sizeof(_falltrackCosTheta[0]), -99999.); 
  std::fill(_falltrackSinTheta, _falltrackSinTheta + sizeof(_falltrackSinTheta)/sizeof(_falltrackSinTheta[0]), -99999.); 
  std::fill(_falltrackPhi, _falltrackPhi + sizeof(_falltrackPhi)/sizeof(_falltrackPhi[0]), -99999.); 
  std::fill(_falltrackCosPhi, _falltrackCosPhi + sizeof(_falltrackCosPhi)/sizeof(_falltrackCosPhi[0]), -99999.); 
  std::fill(_falltrackSinPhi, _falltrackSinPhi + sizeof(_falltrackSinPhi)/sizeof(_falltrackSinPhi[0]), -99999.); 

  std::fill(_fseltrackVrtxDis, _fseltrackVrtxDis + sizeof(_fseltrackVrtxDis)/sizeof(_fseltrackVrtxDis[0]), -99999.); 
  std::fill(_fselntrackhits,  _fselntrackhits + sizeof(_fselntrackhits)/sizeof(_fselntrackhits[0]), -99999.);  std::fill(_fseltrackStartx, _fseltrackStartx + sizeof(_fseltrackStartx)/sizeof(_fseltrackStartx[0]), -99999.); 
  std::fill(_fseltrackID, _fseltrackID + sizeof(_fseltrackID)/sizeof(_fseltrackID[0]), -99999.); 
  std::fill(_fseltrackStartx, _fseltrackStartx + sizeof(_fseltrackStartx)/sizeof(_fseltrackStartx[0]), -99999.); 
  std::fill(_fseltrackStarty, _fseltrackStarty + sizeof(_fseltrackStarty)/sizeof(_fseltrackStarty[0]), -99999.); 
  std::fill(_fseltrackStartz, _fseltrackStartz + sizeof(_fseltrackStartz)/sizeof(_fseltrackStartz[0]), -99999.); 
  std::fill(_fseltrackEndx, _fseltrackEndx + sizeof(_fseltrackEndx)/sizeof(_fseltrackEndx[0]), -99999.);
  std::fill(_fseltrackEndy, _fseltrackEndy + sizeof(_fseltrackEndy)/sizeof(_fseltrackEndy[0]), -99999.); 
  std::fill(_fseltrackEndz, _fseltrackEndz + sizeof(_fseltrackEndz)/sizeof(_fseltrackEndz[0]), -99999.); 
  std::fill(_fseltrackLength, _fseltrackLength + sizeof(_fseltrackLength)/sizeof(_fseltrackLength[0]), -99999.); 
  std::fill(_fseltrackTheta, _fseltrackTheta + sizeof(_fseltrackTheta)/sizeof(_fseltrackTheta[0]), -99999.); 
  std::fill(_fseltrackCosTheta, _fseltrackCosTheta + sizeof(_fseltrackCosTheta)/sizeof(_fseltrackCosTheta[0]), -99999.); 
  std::fill(_fseltrackSinTheta, _fseltrackSinTheta + sizeof(_fseltrackSinTheta)/sizeof(_fseltrackSinTheta[0]), -99999.); 
  std::fill(_fseltrackPhi, _fseltrackPhi + sizeof(_fseltrackPhi)/sizeof(_fseltrackPhi[0]), -99999.); 
  std::fill(_fseltrackCosPhi, _fseltrackCosPhi + sizeof(_fseltrackCosPhi)/sizeof(_fseltrackCosPhi[0]), -99999.); 
  std::fill(_fseltrackSinPhi, _fseltrackSinPhi + sizeof(_fseltrackSinPhi)/sizeof(_fseltrackSinPhi[0]), -99999.); 

   }
   //-----------------------------------------------------------------------
   void CTMDataONAna::beginJob()
   {
    art::ServiceHandle<art::TFileService> tfs;
   
    if (fDoHists)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fNFlashPerEvent   = tfs->make<TH1D>("NFlashEvent", ";Flash/Event",     200,   0.,  200.);
        fFlashPE          = tfs->make<TH1D>("FlashPE",     ";PE",              100,   0.,  200.);
        fFlashTime        = tfs->make<TH1D>("FlashTime",   ";Flash Time(us)",  100, -10.,   30.);
    }
    fDataTreeTotal = tfs->make<TTree>("fDataTreeTotal","Data Holder");
    fDataTreeAllEvts = tfs->make<TTree>("fDataTreeAllEvts","Data Holder");
    fDataTreeFlashTag = tfs->make<TTree>("fDataTreeFlashTag","Data Holder");
    fDataTreeVrtxinFV = tfs->make<TTree>("fDataTreeVrtxinFV","Data Holder");
    fDataTreeTrkNearVrtx=tfs->make<TTree>("fDataTreeTrkNearVrtx","Data Holder");
    fDataTreeFlashTrackMatch = tfs->make<TTree>("fDataTreeFlashTrackMatch","Data Holder");    
    fDataTreeCandMuContained = tfs->make<TTree>("fDataTreeCandMuContained","Data Holder");    
    fDataTreeCCSel = tfs->make<TTree>("fDataTreeCCSel","Data Holder");
    fDataTreePassBadReg = tfs->make<TTree>("fDataTreePassBadReg","Data Holder");
    fDataTreePassCosmicReq = tfs->make<TTree>("fDataTreePassCosmicReq","Data Holder");
    fDataTreePassLongTrackVrtxDis = tfs->make<TTree>("fDataTreePassLongTrackVrtxDis","Data Holder");
    fDataTreePassLongcolhits = tfs->make<TTree>("fDataTreePassLongcolhits","Data Holder");    
    fDataTreeFinalSel = tfs->make<TTree>("fDataTreeFinalSel","Data Holder");

    fDataTreeTotal->Branch("_frun",&_frun,"_frun/I");
    fDataTreeTotal->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeTotal->Branch("_fevent",&_fevent,"_fevent/I");

    fDataTreeAllEvts->Branch("_frun",&_frun,"_frun/I");
    fDataTreeAllEvts->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeAllEvts->Branch("_fevent",&_fevent,"_fevent/I");
    fDataTreeAllEvts->Branch("_fflashPE",&_fflashPE,"_fflashPE/F");
    fDataTreeAllEvts->Branch("_fflashTime",&_fflashTime,"_fflashTime/F");

    fDataTreeFlashTag->Branch("_frun",&_frun,"_frun/I");
    fDataTreeFlashTag->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeFlashTag->Branch("_fevent",&_fevent,"_fevent/I");
    fDataTreeFlashTag->Branch("_fflashPE",&_fflashPE,"_fflashPE/F");
    fDataTreeFlashTag->Branch("_fflashTime",&_fflashTime,"_fflashTime/F");
    
    fDataTreeVrtxinFV->Branch("_frun",&_frun,"_frun/I");
    fDataTreeVrtxinFV->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeVrtxinFV->Branch("_fevent",&_fevent,"_fevent/I");
    fDataTreeVrtxinFV->Branch("_fvrtxxPos",&_fvrtxxPos,"_fvrtxxPos/F");
    fDataTreeVrtxinFV->Branch("_fvrtxyPos",&_fvrtxyPos,"_fvrtxyPos/F");
    fDataTreeVrtxinFV->Branch("_fvrtxzPos",&_fvrtxzPos,"_fvrtxzPos/F");

    fDataTreeTrkNearVrtx->Branch("_frun",&_frun,"_frun/I");
    fDataTreeTrkNearVrtx->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeTrkNearVrtx->Branch("_fevent",&_fevent,"_fevent/I");
    fDataTreeTrkNearVrtx->Branch("_fvrtxxPos",&_fvrtxxPos,"_fvrtxxPos/F");
    fDataTreeTrkNearVrtx->Branch("_fvrtxyPos",&_fvrtxyPos,"_fvrtxyPos/F");
    fDataTreeTrkNearVrtx->Branch("_fvrtxzPos",&_fvrtxzPos,"_fvrtxzPos/F");
    fDataTreeTrkNearVrtx->Branch("_ftrkVrtxDis",&_ftrkVrtxDis,"_ftrkVrtxDis/F");    

    fDataTreeFlashTrackMatch->Branch("_frun",&_frun,"_frun/I");
    fDataTreeFlashTrackMatch->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeFlashTrackMatch->Branch("_fevent",&_fevent,"_fevent/I");
    fDataTreeFlashTrackMatch->Branch("_fvrtxxPos",&_fvrtxxPos,"_fvrtxxPos/F");
    fDataTreeFlashTrackMatch->Branch("_fvrtxyPos",&_fvrtxyPos,"_fvrtxyPos/F");
    fDataTreeFlashTrackMatch->Branch("_fvrtxzPos",&_fvrtxzPos,"_fvrtxzPos/F");
    fDataTreeFlashTrackMatch->Branch("_ftrkVrtxDis",&_ftrkVrtxDis,"_ftrkVrtxDis/F");    
    fDataTreeFlashTrackMatch->Branch("_ftrkcandlength",&_ftrkcandlength,"_ftrkcandlength/F");
    fDataTreeFlashTrackMatch->Branch("_fflashtrackmatch",&_fflashtrackmatch,"_fflashtrackmatch/F");

    fDataTreeCandMuContained->Branch("_frun",&_frun,"_frun/I");
    fDataTreeCandMuContained->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeCandMuContained->Branch("_fevent",&_fevent,"_fevent/I");
    fDataTreeCandMuContained->Branch("_fvrtxxPos",&_fvrtxxPos,"_fvrtxxPos/F");
    fDataTreeCandMuContained->Branch("_fvrtxyPos",&_fvrtxyPos,"_fvrtxyPos/F");
    fDataTreeCandMuContained->Branch("_fvrtxzPos",&_fvrtxzPos,"_fvrtxzPos/F");
    fDataTreeCandMuContained->Branch("_ftrkVrtxDis",&_ftrkVrtxDis,"_ftrkVrtxDis/F");    
    fDataTreeCandMuContained->Branch("_fflashtrackmatch",&_fflashtrackmatch,"_fflashtrackmatch/F");
    fDataTreeCandMuContained->Branch("_ftrkcandlength",&_ftrkcandlength,"_ftrkcandlength/F");

    fDataTreeCCSel->Branch("_frun",&_frun,"_frun/I");
    fDataTreeCCSel->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeCCSel->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreeCCSel->Branch("_fflashtrackmatch",&_fflashtrackmatch,"_fflashtrackmatch/F");
    fDataTreeCCSel->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreeCCSel->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreeCCSel->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreeCCSel->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreeCCSel->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreeCCSel->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreeCCSel->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");	
    fDataTreeCCSel->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreeCCSel->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreeCCSel->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreeCCSel->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
    fDataTreeCCSel->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreeCCSel->Branch("_flongtrackTheta",&_flongtrackTheta,"_flongtrackTheta/F");
    fDataTreeCCSel->Branch("_flongtrackCosTheta",&_flongtrackCosTheta,"_flongtrackCosTheta/F");
    fDataTreeCCSel->Branch("_flongtrackSinTheta",&_flongtrackSinTheta,"_flongtrackSinTheta/F");
    fDataTreeCCSel->Branch("_flongtrackPhi",&_flongtrackPhi,"_flongtrackPhi/F");
    fDataTreeCCSel->Branch("_flongtrackCosPhi",&_flongtrackCosPhi,"_flongtrackCosPhi/F");
    fDataTreeCCSel->Branch("_flongtrackSinPhi",&_flongtrackSinPhi,"_flongtrackSinPhi/F");
    fDataTreeCCSel->Branch("_flongtrackLength",&_flongtrackLength,"_flongtrackLength/F");
    fDataTreeCCSel->Branch("_flongtrackmcsfwdmom",&_flongtrackmcsfwdmom,"_flongtrackmcsfwdmom/F");
    fDataTreeCCSel->Branch("_flongtrackmcsfwdll",&_flongtrackmcsfwdll,"_flongtrackmcsfwdll/F");
    fDataTreeCCSel->Branch("_flongtrackmcsfwderr",&_flongtrackmcsfwderr,"_flongtrackmcsfwderr/F");
    fDataTreeCCSel->Branch("_flongtrackmcsbwdmom",&_flongtrackmcsbwdmom,"_flongtrackmcsbwdmom/F");
    fDataTreeCCSel->Branch("_flongtrackmcsbwdll",&_flongtrackmcsbwdll,"_flongtrackmcsbwdll/F");
    fDataTreeCCSel->Branch("_flongtrackmcsbwderr",&_flongtrackmcsbwderr,"_flongtrackmcsbwderr/F");
    fDataTreeCCSel->Branch("_fbackward_track",&_fbackward_track,"_fbackward_track/I");
    fDataTreeCCSel->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreeCCSel->Branch("_fbadReg_removed",&_fbadReg_removed,"_fbadReg_removed/I");
    fDataTreeCCSel->Branch("_fstarty_cut",&_fstarty_cut,"_fstarty_cut/I");
    fDataTreeCCSel->Branch("_fbrokenTracks",&_fbrokenTracks,"_fbrokenTracks/I");
    fDataTreeCCSel->Branch("_fNRecoallPart",&_fNRecoallPart,"_fNRecoallPart/I");
/*    fDataTreeCCSel->Branch("_falltrackVrtxDis",&_falltrackVrtxDis,"_falltrackVrtxDis[_fNRecoallPart]/F");	
    fDataTreeCCSel->Branch("_falltrackID",&_falltrackID,"_falltrackID/I"); 
    fDataTreeCCSel->Branch("_falltrackStartx",&_falltrackStartx,"_falltrackStartx[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackStarty",&_falltrackStarty,"_falltrackStarty[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackStartz",&_falltrackStartz,"_falltrackStartz[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackEndx",&_falltrackEndx,"_falltrackEndx[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackEndy",&_falltrackEndy,"_falltrackEndy[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackEndz",&_falltrackEndz,"_falltrackEndz[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackLength",&_falltrackLength,"_falltrackLength[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackTheta",&_falltrackTheta,"_falltrackTheta[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackCosTheta",&_falltrackCosTheta,"_falltrackCosTheta[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackSinTheta",&_falltrackSinTheta,"_falltrackSinTheta[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackPhi",&_falltrackPhi,"_falltrackPhi[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackCosPhi",&_falltrackCosPhi,"_falltrackCosPhi[_fNRecoallPart]/F");
    fDataTreeCCSel->Branch("_falltrackSinPhi",&_falltrackSinPhi,"_falltrackSinPhi[_fNRecoallPart]/F");
*/    fDataTreeCCSel->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreeCCSel->Branch("_fseltrackVrtxDis",&_fseltrackVrtxDis,"_fseltrackVrtxDis[_fNrecomult]/F");	
    fDataTreeCCSel->Branch("_fseltrackID",&_fseltrackID,"_fseltrackID/I"); 
    fDataTreeCCSel->Branch("_fselntrackhits",&_fselntrackhits,"_fselntrackhits[_fNrecomult]/I");
    fDataTreeCCSel->Branch("_fseltrackStartx",&_fseltrackStartx,"_fseltrackStartx[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackStarty",&_fseltrackStarty,"_fseltrackStarty[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackStartz",&_fseltrackStartz,"_fseltrackStartz[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackEndx",&_fseltrackEndx,"_fseltrackEndx[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackEndy",&_fseltrackEndy,"_fseltrackEndy[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackEndz",&_fseltrackEndz,"_fseltrackEndz[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackLength",&_fseltrackLength,"_fseltrackLength[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackTheta",&_fseltrackTheta,"_fseltrackTheta[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackCosTheta",&_fseltrackCosTheta,"_fseltrackCosTheta[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackSinTheta",&_fseltrackSinTheta,"_fseltrackSinTheta[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackPhi",&_fseltrackPhi,"_fseltrackPhi[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackCosPhi",&_fseltrackCosPhi,"_fseltrackCosPhi[_fNrecomult]/F");
    fDataTreeCCSel->Branch("_fseltrackSinPhi",&_fseltrackSinPhi,"_fseltrackSinPhi[_fNrecomult]/F");
     
    fDataTreePassBadReg->Branch("_frun",&_frun,"_frun/I");
    fDataTreePassBadReg->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreePassBadReg->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreePassBadReg->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreePassBadReg->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreePassBadReg->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreePassBadReg->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreePassBadReg->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreePassBadReg->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreePassBadReg->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreePassBadReg->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreePassBadReg->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreePassBadReg->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");	
    fDataTreePassBadReg->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreePassBadReg->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreePassBadReg->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreePassBadReg->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
    
    fDataTreePassCosmicReq->Branch("_frun",&_frun,"_frun/I");
    fDataTreePassCosmicReq->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreePassCosmicReq->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreePassCosmicReq->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreePassCosmicReq->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreePassCosmicReq->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreePassCosmicReq->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreePassCosmicReq->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreePassCosmicReq->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreePassCosmicReq->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreePassCosmicReq->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreePassCosmicReq->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreePassCosmicReq->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");  
    fDataTreePassCosmicReq->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreePassCosmicReq->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreePassCosmicReq->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreePassCosmicReq->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
         
    fDataTreePassLongTrackVrtxDis->Branch("_frun",&_frun,"_frun/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreePassLongTrackVrtxDis->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreePassLongTrackVrtxDis->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreePassLongTrackVrtxDis->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreePassLongTrackVrtxDis->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreePassLongTrackVrtxDis->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");  
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreePassLongTrackVrtxDis->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
    
    fDataTreePassLongcolhits->Branch("_frun",&_frun,"_frun/I");
    fDataTreePassLongcolhits->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreePassLongcolhits->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreePassLongcolhits->Branch("_fflashtrackmatch",&_fflashtrackmatch,"_fflashtrackmatch/F");
    fDataTreePassLongcolhits->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreePassLongcolhits->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreePassLongcolhits->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreePassLongcolhits->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreePassLongcolhits->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreePassLongcolhits->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreePassLongcolhits->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");	
    fDataTreePassLongcolhits->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreePassLongcolhits->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreePassLongcolhits->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreePassLongcolhits->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
    fDataTreePassLongcolhits->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreePassLongcolhits->Branch("_flongtrackTheta",&_flongtrackTheta,"_flongtrackTheta/F");
    fDataTreePassLongcolhits->Branch("_flongtrackCosTheta",&_flongtrackCosTheta,"_flongtrackCosTheta/F");
    fDataTreePassLongcolhits->Branch("_flongtrackSinTheta",&_flongtrackSinTheta,"_flongtrackSinTheta/F");
    fDataTreePassLongcolhits->Branch("_flongtrackPhi",&_flongtrackPhi,"_flongtrackPhi/F");
    fDataTreePassLongcolhits->Branch("_flongtrackCosPhi",&_flongtrackCosPhi,"_flongtrackCosPhi/F");
    fDataTreePassLongcolhits->Branch("_flongtrackSinPhi",&_flongtrackSinPhi,"_flongtrackSinPhi/F");
    fDataTreePassLongcolhits->Branch("_flongtrackLength",&_flongtrackLength,"_flongtrackLength/F");
    fDataTreePassLongcolhits->Branch("_flongtrackmcsfwdmom",&_flongtrackmcsfwdmom,"_flongtrackmcsfwdmom/F");
    fDataTreePassLongcolhits->Branch("_flongtrackmcsfwdll",&_flongtrackmcsfwdll,"_flongtrackmcsfwdll/F");
    fDataTreePassLongcolhits->Branch("_flongtrackmcsfwderr",&_flongtrackmcsfwderr,"_flongtrackmcsfwderr/F");
    fDataTreePassLongcolhits->Branch("_flongtrackmcsbwdmom",&_flongtrackmcsbwdmom,"_flongtrackmcsbwdmom/F");
    fDataTreePassLongcolhits->Branch("_flongtrackmcsbwdll",&_flongtrackmcsbwdll,"_flongtrackmcsbwdll/F");
    fDataTreePassLongcolhits->Branch("_flongtrackmcsbwderr",&_flongtrackmcsbwderr,"_flongtrackmcsbwderr/F");
    fDataTreePassLongcolhits->Branch("_fbackward_track",&_fbackward_track,"_fbackward_track/I");
    fDataTreePassLongcolhits->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreePassLongcolhits->Branch("_fbadReg_removed",&_fbadReg_removed,"_fbadReg_removed/I");
    fDataTreePassLongcolhits->Branch("_fstarty_cut",&_fstarty_cut,"_fstarty_cut/I");
    fDataTreePassLongcolhits->Branch("_fbrokenTracks",&_fbrokenTracks,"_fbrokenTracks/I");
    fDataTreePassLongcolhits->Branch("_fNRecoallPart",&_fNRecoallPart,"_fNRecoallPart/I");
/*    fDataTreePassLongcolhits->Branch("_falltrackVrtxDis",&_falltrackVrtxDis,"_falltrackVrtxDis[_fNRecoallPart]/F");	
    fDataTreePassLongcolhits->Branch("_falltrackID",&_falltrackID,"_falltrackID/I"); 
    fDataTreePassLongcolhits->Branch("_falltrackStartx",&_falltrackStartx,"_falltrackStartx[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackStarty",&_falltrackStarty,"_falltrackStarty[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackStartz",&_falltrackStartz,"_falltrackStartz[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackEndx",&_falltrackEndx,"_falltrackEndx[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackEndy",&_falltrackEndy,"_falltrackEndy[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackEndz",&_falltrackEndz,"_falltrackEndz[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackLength",&_falltrackLength,"_falltrackLength[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackTheta",&_falltrackTheta,"_falltrackTheta[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackCosTheta",&_falltrackCosTheta,"_falltrackCosTheta[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackSinTheta",&_falltrackSinTheta,"_falltrackSinTheta[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackPhi",&_falltrackPhi,"_falltrackPhi[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackCosPhi",&_falltrackCosPhi,"_falltrackCosPhi[_fNRecoallPart]/F");
    fDataTreePassLongcolhits->Branch("_falltrackSinPhi",&_falltrackSinPhi,"_falltrackSinPhi[_fNRecoallPart]/F");
*/    fDataTreePassLongcolhits->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreePassLongcolhits->Branch("_fseltrackVrtxDis",&_fseltrackVrtxDis,"_fseltrackVrtxDis[_fNrecomult]/F");	
    fDataTreePassLongcolhits->Branch("_fseltrackID",&_fseltrackID,"_fseltrackID/I"); 
    fDataTreePassLongcolhits->Branch("_fselntrackhits",&_fselntrackhits,"_fselntrackhits[_fNrecomult]/I");
    fDataTreePassLongcolhits->Branch("_fseltrackStartx",&_fseltrackStartx,"_fseltrackStartx[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackStarty",&_fseltrackStarty,"_fseltrackStarty[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackStartz",&_fseltrackStartz,"_fseltrackStartz[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackEndx",&_fseltrackEndx,"_fseltrackEndx[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackEndy",&_fseltrackEndy,"_fseltrackEndy[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackEndz",&_fseltrackEndz,"_fseltrackEndz[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackLength",&_fseltrackLength,"_fseltrackLength[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackTheta",&_fseltrackTheta,"_fseltrackTheta[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackCosTheta",&_fseltrackCosTheta,"_fseltrackCosTheta[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackSinTheta",&_fseltrackSinTheta,"_fseltrackSinTheta[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackPhi",&_fseltrackPhi,"_fseltrackPhi[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackCosPhi",&_fseltrackCosPhi,"_fseltrackCosPhi[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fseltrackSinPhi",&_fseltrackSinPhi,"_fseltrackSinPhi[_fNrecomult]/F");
    fDataTreePassLongcolhits->Branch("_fPH",&_fPH,"_fPH/I");
    fDataTreePassLongcolhits->Branch("_fMCS",&_fMCS,"_fMCS/I");
    fDataTreePassLongcolhits->Branch("_fMCS1",&_fMCS1,"_fMCS1/I");
    fDataTreePassLongcolhits->Branch("_fMCS2",&_fMCS2,"_fMCS2/I");
    fDataTreePassLongcolhits->Branch("_fPHratio",&_fPHratio,"_fPHratio/F");
    fDataTreePassLongcolhits->Branch("_fMCSratio",&_fMCSratio,"_fMCSratio/F");
    fDataTreePassLongcolhits->Branch("_fMCSratio1",&_fMCSratio1,"_fMCSratio1/F");
    fDataTreePassLongcolhits->Branch("_fMCSdiff",&_fMCSdiff,"_fMCSdiff/F");
    
    fDataTreeFinalSel->Branch("_frun",&_frun,"_frun/I");
    fDataTreeFinalSel->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreeFinalSel->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreeFinalSel->Branch("_fflashtrackmatch",&_fflashtrackmatch,"_fflashtrackmatch/F");
    fDataTreeFinalSel->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreeFinalSel->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreeFinalSel->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreeFinalSel->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreeFinalSel->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreeFinalSel->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreeFinalSel->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");	
    fDataTreeFinalSel->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreeFinalSel->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreeFinalSel->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreeFinalSel->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
    fDataTreeFinalSel->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreeFinalSel->Branch("_flongtrackTheta",&_flongtrackTheta,"_flongtrackTheta/F");
    fDataTreeFinalSel->Branch("_flongtrackCosTheta",&_flongtrackCosTheta,"_flongtrackCosTheta/F");
    fDataTreeFinalSel->Branch("_flongtrackSinTheta",&_flongtrackSinTheta,"_flongtrackSinTheta/F");
    fDataTreeFinalSel->Branch("_flongtrackPhi",&_flongtrackPhi,"_flongtrackPhi/F");
    fDataTreeFinalSel->Branch("_flongtrackCosPhi",&_flongtrackCosPhi,"_flongtrackCosPhi/F");
    fDataTreeFinalSel->Branch("_flongtrackSinPhi",&_flongtrackSinPhi,"_flongtrackSinPhi/F");
    fDataTreeFinalSel->Branch("_flongtrackLength",&_flongtrackLength,"_flongtrackLength/F");
    fDataTreeFinalSel->Branch("_flongtrackmcsfwdmom",&_flongtrackmcsfwdmom,"_flongtrackmcsfwdmom/F");
    fDataTreeFinalSel->Branch("_flongtrackmcsfwdll",&_flongtrackmcsfwdll,"_flongtrackmcsfwdll/F");
    fDataTreeFinalSel->Branch("_flongtrackmcsfwderr",&_flongtrackmcsfwderr,"_flongtrackmcsfwderr/F");
    fDataTreeFinalSel->Branch("_flongtrackmcsbwdmom",&_flongtrackmcsbwdmom,"_flongtrackmcsbwdmom/F");
    fDataTreeFinalSel->Branch("_flongtrackmcsbwdll",&_flongtrackmcsbwdll,"_flongtrackmcsbwdll/F");
    fDataTreeFinalSel->Branch("_flongtrackmcsbwderr",&_flongtrackmcsbwderr,"_flongtrackmcsbwderr/F");
    fDataTreeFinalSel->Branch("_fbackward_track",&_fbackward_track,"_fbackward_track/I");
    fDataTreeFinalSel->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreeFinalSel->Branch("_fbadReg_removed",&_fbadReg_removed,"_fbadReg_removed/I");
    fDataTreeFinalSel->Branch("_fstarty_cut",&_fstarty_cut,"_fstarty_cut/I");
    fDataTreeFinalSel->Branch("_fbrokenTracks",&_fbrokenTracks,"_fbrokenTracks/I");
    fDataTreeFinalSel->Branch("_fNRecoallPart",&_fNRecoallPart,"_fNRecoallPart/I");
/*    fDataTreeFinalSel->Branch("_falltrackVrtxDis",&_falltrackVrtxDis,"_falltrackVrtxDis[_fNRecoallPart]/F");	
    fDataTreeFinalSel->Branch("_falltrackID",&_falltrackID,"_falltrackID/I"); 
    fDataTreeFinalSel->Branch("_falltrackStartx",&_falltrackStartx,"_falltrackStartx[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackStarty",&_falltrackStarty,"_falltrackStarty[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackStartz",&_falltrackStartz,"_falltrackStartz[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackEndx",&_falltrackEndx,"_falltrackEndx[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackEndy",&_falltrackEndy,"_falltrackEndy[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackEndz",&_falltrackEndz,"_falltrackEndz[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackLength",&_falltrackLength,"_falltrackLength[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackTheta",&_falltrackTheta,"_falltrackTheta[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackCosTheta",&_falltrackCosTheta,"_falltrackCosTheta[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackSinTheta",&_falltrackSinTheta,"_falltrackSinTheta[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackPhi",&_falltrackPhi,"_falltrackPhi[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackCosPhi",&_falltrackCosPhi,"_falltrackCosPhi[_fNRecoallPart]/F");
    fDataTreeFinalSel->Branch("_falltrackSinPhi",&_falltrackSinPhi,"_falltrackSinPhi[_fNRecoallPart]/F");
*/    fDataTreeFinalSel->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreeFinalSel->Branch("_fseltrackVrtxDis",&_fseltrackVrtxDis,"_fseltrackVrtxDis[_fNrecomult]/F");	
    fDataTreeFinalSel->Branch("_fseltrackID",&_fseltrackID,"_fseltrackID/I"); 
    fDataTreeFinalSel->Branch("_fselntrackhits",&_fselntrackhits,"_fselntrackhits[_fNrecomult]/I");
    fDataTreeFinalSel->Branch("_fseltrackStartx",&_fseltrackStartx,"_fseltrackStartx[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackStarty",&_fseltrackStarty,"_fseltrackStarty[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackStartz",&_fseltrackStartz,"_fseltrackStartz[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackEndx",&_fseltrackEndx,"_fseltrackEndx[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackEndy",&_fseltrackEndy,"_fseltrackEndy[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackEndz",&_fseltrackEndz,"_fseltrackEndz[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackLength",&_fseltrackLength,"_fseltrackLength[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackTheta",&_fseltrackTheta,"_fseltrackTheta[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackCosTheta",&_fseltrackCosTheta,"_fseltrackCosTheta[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackSinTheta",&_fseltrackSinTheta,"_fseltrackSinTheta[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackPhi",&_fseltrackPhi,"_fseltrackPhi[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackCosPhi",&_fseltrackCosPhi,"_fseltrackCosPhi[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fseltrackSinPhi",&_fseltrackSinPhi,"_fseltrackSinPhi[_fNrecomult]/F");
    fDataTreeFinalSel->Branch("_fPH",&_fPH,"_fPH/I");
    fDataTreeFinalSel->Branch("_fMCS",&_fMCS,"_fMCS/I");
    fDataTreeFinalSel->Branch("_fMCS1",&_fMCS1,"_fMCS1/I");
    fDataTreeFinalSel->Branch("_fMCS2",&_fMCS2,"_fMCS2/I");
    fDataTreeFinalSel->Branch("_fPHratio",&_fPHratio,"_fPHratio/F");
    fDataTreeFinalSel->Branch("_fMCSratio",&_fMCSratio,"_fMCSratio/F");
    fDataTreeFinalSel->Branch("_fMCSratio1",&_fMCSratio1,"_fMCSratio1/F");
    fDataTreeFinalSel->Branch("_fMCSdiff",&_fMCSdiff,"_fMCSdiff/F");
             
   }
   
   //-----------------------------------------------------------------------
   void CTMDataONAna::endJob()
  {

  }
  
   //-----------------------------------------------------------------------
   void CTMDataONAna::beginRun(const art::Run& run)
  {
    
  }
    //-----------------------------------------------------------------------
   void CTMDataONAna::reconfigure(fhicl::ParameterSet const& pset)
   {
    fOpFlashModuleLabel      = pset.get<std::string  > ("OpFlashModuleLabel", "simpleFlashBeam");
    //fOpFlashModuleLabel      = pset.get<std::string  > ("OpFlashModuleLabel", "opflash"); //In MCC7
    fHitsModuleLabel         = pset.get< std::string > ("HitsModuleLabel", "gaushit");
    //fTrackModuleLabel        = pset.get< std::string > ("TrackModuleLabel",  "pandoraNu");
    //fVertexModuleLabel       = pset.get< std::string > ("VertexModuleLabel",  "pandoraNu");
    fTrackModuleLabel        = pset.get< std::string > ("TrackModuleLabel",  "pandoraNu");
    fVertexModuleLabel       = pset.get< std::string > ("VertexModuleLabel",  "pandoraNu");
    fGenieGenModuleLabel     = pset.get< std::string > ("GenieGenModuleLabel", "generator");    
    fTrackMCSFitLabel        = pset.get< std::string > ("TrackMCSFitLabel", "pandoraNuMCSMu");

    fDistToEdgeX             = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",   10.);
    fDistToEdgeY             = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   20.);
    fDistToEdgeZ             = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   10.);
    
    fFlashWidth              = pset.get<double>      ("FlashWidth", 80.);
    fBeamMin                 = pset.get<double>      ("BeamMin", 3.3);
    fBeamMax                 = pset.get<double>      ("BeamMax", 4.9);
    fPEThresh                = pset.get<double>      ("PEThresh", 50.);
    fMinTrk2VtxDist          = pset.get<double>      ("MinTrk2VtxDist", 5.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen", 75.);
   }
 
 //========================================================================	
bool CTMDataONAna::inFV(double x, double y, double z) const
{
    double distInX = x - fGeometry->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * fGeometry->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}
 //=======================================================================
//This function returns the distance between a flash and
//a track (in one dimension, here used only for z direction)
double CTMDataONAna::FlashTrackDist(double flash, double start, double end) const
{
    if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if(flash > end && flash < start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
}

   //-----------------------------------------------------------------------

   void CTMDataONAna::analyze(art::Event const& evt) 
   {
    _frun = evt.run();
    _fsubrun = evt.subRun();
    _fevent = evt.id().event();
    
    fDataTreeTotal->Fill();
    
    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle<std::vector<recob::Vertex>>  vertexVecHandle;
    art::Handle<std::vector<recob::Track>>   trackVecHandle;
    art::Handle<std::vector<recob::OpFlash>> flashListHandle;
    
    evt.getByLabel(fVertexModuleLabel,    vertexVecHandle);
    evt.getByLabel(fTrackModuleLabel,     trackVecHandle);
    
    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    
    if (evt.getByLabel(fOpFlashModuleLabel,flashListHandle))
        art::fill_ptr_vector(flashlist, flashListHandle);
    
    // Require valid handles, otherwise nothing to do
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
    {
        // Recover associations to PFParticles...
        art::FindManyP<recob::PFParticle> trackToPFPartAssns(trackVecHandle,  evt, fTrackModuleLabel);
        
        //----loop over all the flashes and check if there are flashes within the beam
        //window and above the PE threshold
        const recob::OpFlash* flashPtr(0);
        double                flashPE(0);
        double                flashtime(0);
        bool                  flashtag(false);
        
        for(const auto& opFlash : flashlist)
        {
            if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
            {
                flashtag = true;
                
                // Keep track of the largest flash
                if (opFlash->TotalPE() > flashPE)
                {
                    flashPtr = opFlash.get();
                    flashPE = opFlash->TotalPE();
		    flashtime = opFlash->Time();
                }
            }
            
            if (fDoHists)
            {
                fFlashPE->Fill(opFlash->TotalPE(), 1.);
                fFlashTime->Fill(opFlash->Time(), 1.);
            }
         }  //end of loop over all the flashes
        
	_fflashPE = flashPE;
	_fflashTime = flashtime;
        fDataTreeAllEvts->Fill();

        if (fDoHists) fNFlashPerEvent->Fill(flashlist.size(), 1.);
        
        if(flashtag)
        {
	fDataTreeFlashTag->Fill();
            // We need to keep track of the best combination
            // Can we assign art ptrs? I don't think so...
            int    VertexCandidate=-1;
            int    TrackCandidate=-1;
            double TrackCandLength=0;
            double trackstartzcandidate=0;
            double trackstartxcandidate=0;
            double trackstartycandidate=0;
            double trackendzcandidate=0;
            double trackendxcandidate=0;
            double trackendycandidate=0;
	    double trackthetacandidate=0;
            TVector3 trackPos(0,0,0);
            TVector3 trackEnd(0,0,0);
	    TVector3 vertexPos(0,0,0);
	    double trackTheta = -1;
            double trkflipcandidate=0;
	    double vrtxxcandidate=-1;
	    double vrtxycandidate=-1;
	    double vrtxzcandidate=-1;
	    	    
            //-----------------------------------------------------------
            bool vrtxinFV = true;
            for(size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
            {
                // Recover art ptr to vertex
                art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);
                
                // Get the position of the vertex
                // Ultimately we really want the vertex position in a TVector3 object...
                double vertexXYZ[3];
                
                vertex->XYZ(vertexXYZ);
                
                TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);
            	   
                if(inFV(vertexPos.X(),vertexPos.Y(),vertexPos.Z()))
                {
///////////////////////////////Just for event count//////////////////
		    if(vrtxinFV)
		    {
		      _fvrtxxPos= vertexPos.X();
		      _fvrtxyPos= vertexPos.Y();
		      _fvrtxzPos= vertexPos.Z();
		      fDataTreeVrtxinFV->Fill();
		    }  
		    vrtxinFV = false;
///////////////////////////////////////////////////////////////		    
                    // For each vertex we loop over all tracks looking for matching pairs
                    // The outer loop here, then is over one less than all tracks
                    for(size_t trackIdx = 0; trackIdx < trackVecHandle->size(); trackIdx++)
                    {
                        // Work with an art Ptr here
                        art::Ptr<recob::Track> track(trackVecHandle,trackIdx);
                        
                        // so we need to get the track direction sorted out.
                        trackPos = track->Vertex();
                        trackEnd = track->End();
			trackTheta = track->Theta();
                        
                        // Take the closer end---------------------------------
                        double trackToVertexDist = (trackPos - vertexPos).Mag();
			double trkflip = 0;
                        
                        if ((trackEnd - vertexPos).Mag() < trackToVertexDist)
                        {
                            trackPos          = track->End();
                            trackEnd          = track->Vertex();
			    trackTheta        = -(track->Theta());

                            trackToVertexDist = (trackPos - vertexPos).Mag();
			    trkflip= 1;

                        }
                        
                        //--------------------------------------------------------------------------
                        if(trackToVertexDist<fMinTrk2VtxDist)
                        {
                            if((trackEnd-trackPos).Mag()>TrackCandLength)
                            {
                                TrackCandLength = (trackEnd-trackPos).Mag();
                                TrackCandidate=trackIdx;
                                VertexCandidate=vertexIdx;
                                trackstartzcandidate=trackPos.z();
                                trackstartxcandidate=trackPos.x();
                                trackstartycandidate=trackPos.y();
                                trackendzcandidate=trackEnd.z();
                                trackendxcandidate=trackEnd.x();
                                trackendycandidate=trackEnd.y();
				trackthetacandidate=trackTheta;
				vrtxxcandidate=vertexPos.X();
				vrtxycandidate=vertexPos.Y();
				vrtxzcandidate=vertexPos.Z();
				trkflipcandidate=trkflip;
				
                            }
                        } //end of if track distance is within 5cm
                    }  //end of loop over the tracks
                }  //end of if the vertex is contained
            } //end of loop over all the vertex
            
            if(TrackCandidate > -1)
            {      
                _ftrkcandlength = TrackCandLength;
		TVector3 trackPos(trackstartxcandidate,trackstartycandidate,trackstartzcandidate);
		TVector3 trackEnd(trackendxcandidate,trackendycandidate,trackendzcandidate);
		TVector3 vertexPos(vrtxxcandidate,vrtxycandidate,vrtxzcandidate);
	        _ftrkVrtxDis= (trackPos - vertexPos).Mag();
		_flongtrackflipped= trkflipcandidate;
		trackTheta=trackthetacandidate;

	        
                bool flashTrackMatchFlag = FlashTrackDist(flashPtr->ZCenter(), trackstartzcandidate, trackendzcandidate) < fFlashWidth;
                
		_fflashtrackmatch = FlashTrackDist(flashPtr->ZCenter(), trackstartzcandidate, trackendzcandidate);
		
		fDataTreeTrkNearVrtx->Fill();
               
	        // Check to see if we think we have a candidate
		if(flashTrackMatchFlag)
		{
                  fDataTreeFlashTrackMatch->Fill();		
		if(inFV(trackstartxcandidate, trackstartycandidate, trackstartzcandidate) && inFV(trackendxcandidate, trackendycandidate, trackendzcandidate))
		{
                  fDataTreeCandMuContained->Fill();		
                if(TrackCandLength>fMinTrackLen )
                {
	            //fDataTreeCCSel->Fill();
                    // Make an association between the best vertex and the matching tracks
                    art::Ptr<recob::Vertex> vertex(vertexVecHandle,VertexCandidate);
                    art::Ptr<recob::Track>  longtrack(trackVecHandle,TrackCandidate);

              double vrtxxyz[3];
              vertex->XYZ(vrtxxyz);
              //TVector3 vertexPos(vrtxxyz[0],vrtxxyz[1],vrtxxyz[2]);
              _fvrtxx = vrtxxyz[0];
              _fvrtxy = vrtxxyz[1];
              _fvrtxz = vrtxxyz[2];    

              //getting reconstructed track position information
              _flongtrackstartx = trackPos.X();	
              _flongtrackstarty = trackPos.Y();	
              _flongtrackstartz = trackPos.Z();
              _flongtrackendx   = trackEnd.X();
              _flongtrackendy   = trackEnd.Y();	
              _flongtrackendz   = trackEnd.Z();
              TVector3 longTrackStartPos(_flongtrackstartx,_flongtrackstarty,_flongtrackstartz);

      //------------backward-going tracks---------------------------//
      int backward_track = 0;
      if(_flongtrackstartz >= _flongtrackendz)
      {
        backward_track = 1;
      }   
      _fbackward_track=backward_track;
      
      // Calculating long track to vertex 3D distance 
      _flongTrackVrtxDis = (longTrackStartPos-vertexPos).Mag();
      
      // Calculating other kinematics of the long track
      _flongtrackID     = longtrack->ID();
      _flongtrackTheta = trackTheta;
      _flongtrackCosTheta = cos(trackTheta);
      _flongtrackSinTheta = sin(trackTheta);
      _flongtrackPhi = longtrack->Phi();
      _flongtrackCosPhi = cos(longtrack->Phi());
      _flongtrackSinPhi = sin(longtrack->Phi());
      _flongtrackLength = longtrack->Length();
      
      //--------Calculating long track collection plane hits--------//
      const int kMaxCh              = 8256; //maximum number of channels;
      double hits_charge[kMaxCh]    ={-999};
      int hits_wire[kMaxCh]         ={-999};
      double hits_peakT[kMaxCh]     ={-999};
      int  longcolhits =0;
      
      art::Handle<std::vector<recob::Track> > trackCol;	
      evt.getByLabel(fTrackModuleLabel, trackCol);
      std::vector<recob::Track> const& TrackVector(*trackCol);   

      art::Handle< std::vector<recob::Hit> > hitListHandle;
      std::vector<art::Ptr<recob::Hit> > hitlist;
      if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
      art::FindMany<recob::Hit> fmh(trackCol, evt, fTrackModuleLabel);
      std::vector<const recob::Hit* > hits = fmh.at(longtrack.key());	
      
     for(std::vector<const recob::Hit* >::iterator itr = hits.begin(); itr < hits.end(); itr++)
     {  
     // looping over the hits of the selected long track
       if((*itr)->WireID().Plane == 2)
       {
	 hits_peakT[longcolhits]     =   (*itr)->PeakTime();
	 hits_charge[longcolhits]    =   ((*itr)->Integral())/(TMath::Exp(-(hits_peakT[longcolhits]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
         hits_wire[longcolhits]      =   (*itr)->WireID().Wire;
	 longcolhits++;
       }
     }
     _fnlongcolhits=longcolhits;      

      //------------Creating dead region flag---------------------//
      int badReg_removed = 0;
      if((_flongtrackstartz <(5.4-wire_sideband)    ||  _flongtrackstartz >(6.0+wire_sideband)   )
      && (_flongtrackstartz <(52.8-wire_sideband)   ||  _flongtrackstartz >(57.3+wire_sideband)  )
      && (_flongtrackstartz <(91.2-wire_sideband)   ||  _flongtrackstartz >(95.7+wire_sideband)  )
      && (_flongtrackstartz <(120-wire_sideband)    ||  _flongtrackstartz >(124.5+wire_sideband) )
      && (_flongtrackstartz <(244.8-wire_sideband)  ||  _flongtrackstartz >(249.3+wire_sideband) )
      && (_flongtrackstartz <(288-wire_sideband)    ||  _flongtrackstartz >(291.6+wire_sideband) )
      && (_flongtrackstartz <(398.4-wire_sideband)  ||  _flongtrackstartz >(402.9+wire_sideband) )
      && (_flongtrackstartz <(412.8-wire_sideband)  ||  _flongtrackstartz >(417.3+wire_sideband) )
      && (_flongtrackstartz <(662.4-wire_sideband)  ||  _flongtrackstartz >(668.1+wire_sideband) )
      && (_flongtrackstartz <(700.8-wire_sideband)  ||  _flongtrackstartz >(738.9+wire_sideband) )
      && (_flongtrackstartz <(806.4-wire_sideband)  ||  _flongtrackstartz >(810.9+wire_sideband) )
      && (_flongtrackstartz <(820.8-wire_sideband)  ||  _flongtrackstartz >(825.3+wire_sideband) )
      && (_flongtrackstartz <(873.6-wire_sideband)  ||  _flongtrackstartz >(878.1+wire_sideband) )      
      && (_flongtrackendz <(5.4-wire_sideband)      ||  _flongtrackendz >(6.0+wire_sideband)   )
      && (_flongtrackendz <(52.8-wire_sideband)     ||  _flongtrackendz >(57.3+wire_sideband)  )
      && (_flongtrackendz <(91.2-wire_sideband)     ||  _flongtrackendz >(95.7+wire_sideband)  )
      && (_flongtrackendz <(120-wire_sideband)      ||  _flongtrackendz >(124.5+wire_sideband) )
      && (_flongtrackendz <(244.8-wire_sideband)    ||  _flongtrackendz >(249.3+wire_sideband) )
      && (_flongtrackendz <(288-wire_sideband)      ||  _flongtrackendz >(291.6+wire_sideband) )
      && (_flongtrackendz <(398.4-wire_sideband)    ||  _flongtrackendz >(402.9+wire_sideband) )
      && (_flongtrackendz <(412.8-wire_sideband)    ||  _flongtrackendz >(417.3+wire_sideband) )
      && (_flongtrackendz <(662.4-wire_sideband)    ||  _flongtrackendz >(668.1+wire_sideband) )
      && (_flongtrackendz <(700.8-wire_sideband)    ||  _flongtrackendz >(738.9+wire_sideband) )
      && (_flongtrackendz <(806.4-wire_sideband)    ||  _flongtrackendz >(810.9+wire_sideband) )
      && (_flongtrackendz <(820.8-wire_sideband)    ||  _flongtrackendz >(825.3+wire_sideband) )
      && (_flongtrackendz <(873.6-wire_sideband)    ||  _flongtrackendz >(878.1+wire_sideband)    )) //removing cosmics and bad regions         
      { 
        badReg_removed = 1;
      }
      _fbadReg_removed=badReg_removed;
      
      //----------Extra cosmic removal cut-----------------//
      int starty_cut = 0;
      if (_flongtrackstarty<70)  starty_cut = 1;    
      _fstarty_cut=starty_cut;           

     /////-------- MCS fit for long track -----------/////// 

     art::Handle< std::vector<recob::MCSFitResult>  > mcsfitListHandle;
     std::vector<art::Ptr<recob::MCSFitResult> > mcsfitlist;
     if(evt.getByLabel(fTrackMCSFitLabel,mcsfitListHandle))
     art::fill_ptr_vector(mcsfitlist, mcsfitListHandle);

     const recob::MCSFitResult& mcsMu = mcsfitListHandle->at(longtrack.key());
     if (_flongtrackflipped == 0)
     {
     _flongtrackmcsfwdmom = mcsMu.fwdMomentum();  
     _flongtrackmcsfwdll = mcsMu.fwdLogLikelihood();  
     _flongtrackmcsfwderr = mcsMu.fwdMomUncertainty();  
     _flongtrackmcsbwdmom = mcsMu.bwdMomentum();  
     _flongtrackmcsbwdll = mcsMu.bwdLogLikelihood();  
     _flongtrackmcsbwderr = mcsMu.bwdMomUncertainty();
     }
     if(_flongtrackflipped == 1)  
     {
     _flongtrackmcsfwdmom = mcsMu.bwdMomentum();  
     _flongtrackmcsfwdll = mcsMu.bwdLogLikelihood();  
     _flongtrackmcsfwderr = mcsMu.bwdMomUncertainty();  
     _flongtrackmcsbwdmom = mcsMu.fwdMomentum();  
     _flongtrackmcsbwdll = mcsMu.fwdLogLikelihood();  
     _flongtrackmcsbwderr = mcsMu.fwdMomUncertainty();
     }        

      //-----------------Reconstructed, truth macthing, secondary multiplicty---------------------------//
      int recomult=0;
      _fNRecoallPart = TrackVector.size();

      _fNrecomult = 0;//starting value

      size_t nrecopart = 0;
      size_t nselrecopart = 0;
          
      for (size_t it=0; it < TrackVector.size(); it++)
      {
        auto track = TrackVector.at(it);
	TVector3 trackStart = track.Vertex();
        TVector3 trackEnd = track.End();
        _falltrackID[nrecopart] = track.ID();
	_falltrackStartx[nrecopart] = track.Vertex()(0);
	_falltrackStarty[nrecopart] = track.Vertex()(1);
	_falltrackStartz[nrecopart] = track.Vertex()(2);
	_falltrackEndx[nrecopart] = track.End()(0);
	_falltrackEndy[nrecopart] = track.End()(1);
	_falltrackEndz[nrecopart] = track.End()(2);
	_falltrackLength[nrecopart] = track.Length();
	_falltrackTheta[nrecopart] = track.Theta();
	_falltrackCosTheta[nrecopart] = cos(track.Theta());
	_falltrackSinTheta[nrecopart] = sin(track.Theta());
	_falltrackPhi[nrecopart] = track.Phi();
	_falltrackCosPhi[nrecopart] = cos(track.Phi());
	_falltrackSinPhi[nrecopart] = sin(track.Phi());
	
	if((trackStart-vertexPos).Mag()<(trackEnd-vertexPos).Mag())
	{
	_falltrackVrtxDis[nrecopart]=(trackStart-vertexPos).Mag();
	}
	else if((trackEnd-vertexPos).Mag()<=(trackStart-vertexPos).Mag())
	{
	_falltrackVrtxDis[nrecopart]=(trackEnd-vertexPos).Mag();
	}
	int ntrackhits=0;
//        if(min((trackStart-vertexPos).Mag(),(trackEnd-vertexPos).Mag())<=Disvrtx)	
	if((trackStart-vertexPos).Mag()<=Disvrtx || (trackEnd-vertexPos).Mag()<=Disvrtx)
	{
	  art::FindMany<recob::Hit> fmhit(trackCol, evt, fTrackModuleLabel);
          std::vector<const recob::Hit* > trackhits = fmhit.at(track.ID());	
               
          art::FindManyP<recob::Hit> track_hitsall(trackCol, evt, fTrackModuleLabel);
          std::vector<art::Ptr<recob::Hit>> alltrackhits = track_hitsall.at(track.ID());	  

          for(std::vector<const recob::Hit* >::iterator itr = trackhits.begin(); itr < trackhits.end(); itr++)
          {
	    if((*itr)->WireID().Plane == 2)
            {
	      ntrackhits++; 
	    } 
     	  }
	  if(ntrackhits >= Minhits)
	  {
	    recomult++;
	  _fselntrackhits[nselrecopart] = ntrackhits;
	  _fseltrackVrtxDis[nselrecopart]=_falltrackVrtxDis[nrecopart];
          _fseltrackID[nselrecopart] = track.ID();
	  _fseltrackStartx[nselrecopart] = track.Vertex()(0);
	  _fseltrackStarty[nselrecopart] = track.Vertex()(1);
	  _fseltrackStartz[nselrecopart] = track.Vertex()(2);
	  _fseltrackEndx[nselrecopart] = track.End()(0);
	  _fseltrackEndy[nselrecopart] = track.End()(1);
	  _fseltrackEndz[nselrecopart] = track.End()(2);
	  _fseltrackLength[nselrecopart] = track.Length();
	  _fseltrackTheta[nselrecopart] = track.Theta();
	  _fseltrackCosTheta[nselrecopart] = cos(track.Theta());
	  _fseltrackSinTheta[nselrecopart] = sin(track.Theta());
	  _fseltrackPhi[nselrecopart] = track.Phi();
	  _fseltrackCosPhi[nselrecopart] = cos(track.Phi());
	  _fseltrackSinPhi[nselrecopart] = sin(track.Phi());

	  ++nselrecopart;   
	  }
	}
	  ++nrecopart;
      }	
      _fNrecomult=recomult;
       
     //------------------------Wire ordering----------------------//
     int a; double b,c;
   // for forward-going tracks, arrange hits in accending order  
   if(_fbackward_track == 0)  
   {     
     for(int i=0; i<longcolhits; i++)
     {
       for (int j = i + 1; j < longcolhits; j++)
       { 
         if (hits_wire[i] > hits_wire[j])
         {	
	   a = hits_wire[i];
           hits_wire[i] = hits_wire[j];
           hits_wire[j] = a;
	   //get corresponding charge
	   b = hits_charge[i];
           hits_charge[i] = hits_charge[j];
           hits_charge[j] = b;	
	   //get corresponding peak time 
	   c = hits_peakT[i];
           hits_peakT[i] = hits_peakT[j];
           hits_peakT[j] = c;
	 }
       }
     }
   } 
   // for backward-going tracks, arrange hits in decending order
   if(_fbackward_track == 1)  
   {     
     for(int i=0; i<longcolhits; i++)
     {
       for (int j = i + 1; j < longcolhits; j++)
       { 
         if (hits_wire[i] < hits_wire[j])
         {	
	   a = hits_wire[i];
           hits_wire[i] = hits_wire[j];
           hits_wire[j] = a;
	   //get corresponding charge
	   b = hits_charge[i];
           hits_charge[i] = hits_charge[j];
           hits_charge[j] = b;	
	   //get corresponding peak time 
	   c = hits_peakT[i];
           hits_peakT[i] = hits_peakT[j];
           hits_peakT[j] = c;
	 }
       }
     }
   } 
     //------------Filling CC selected tree with all events---------//

     fDataTreeCCSel->Fill();

    //--------------Event Selection Begins----------------//
    
    if(_fbadReg_removed == 1)
    {
      fDataTreePassBadReg->Fill();
      if(_fstarty_cut == 1)
      {
        fDataTreePassCosmicReq->Fill(); 
	if(_flongTrackVrtxDis <= Disvrtx)
	{
	  fDataTreePassLongTrackVrtxDis->Fill();
          if(_fnlongcolhits >= MinLongTrackHits)
	  {
          // Tree declared towards the end	    

    //--------------------------------PH Test----------------------------------------//		
    int PH =-999;			
    double uphitcharge=0, downhitcharge=0;
    //int upwire[20]={-999}, downwire[20]={-999};
    double upchrg[20]={-999}, downchrg[20]={-999};
    double avgupchrg, avgdownchrg;
    int k=0, l=0, prevk=-999, prevl=-999;
    double prevavgupchrg=-999, prevavgdownchrg=-999;
    	
       //for hits towards vertex			 								
      for(int i=0; i < longcolhits;i++)
      {
	if(i>=10 && i<30)    //Add up corrected charge on 20 hits, 10 hits away from upstream end
	{		    
	  uphitcharge += hits_charge[i];     //adding charge on vertex side 20 hits of the track	
	  //std::cout<<"hit_wire: "<<hits_wire[i]<<" "<<"hit_charge: "<<hits_charge[i]<<"\n"; 
	  upchrg[k]    = hits_charge[i];
	  //upwire[k]    = hits_wire[i];
	  //std::cout<<k<<" hit wire: "<<upwire[k]<<" hit_charge: "<<upchrg[k]<<"\n";
	  k++; 
	}
      }
      //-----------------get Truncated mean charge--------------------//
      //for hits towards vertex
      avgupchrg  =  uphitcharge/k;
      do
      {   
        prevk= k;
	//std::cout<<"avg.up charge: "<<avgupchrg<<"prev k: "<<prevk<<"\n";
	k=0; 
	uphitcharge = 0;
        for(int i=0; i<prevk; i++)
	{
          if(upchrg[i]>0.2*avgupchrg && upchrg[i]<2*avgupchrg)
	  {
            uphitcharge += upchrg[i];
            upchrg[k]    = upchrg[i];
            k++;
       	  }
        }
	prevavgupchrg =  avgupchrg;
	avgupchrg  =  uphitcharge/k;
	//std::cout<<"prev avg up charge: "<<prevavgupchrg<<" avg up charge: "<<avgupchrg<<"\n";	
      }			
      while(avgupchrg!=prevavgupchrg);
      
      //for hits away from vertex	
      for(int i=0; i < longcolhits;i++)
      {	 								
	if(i<longcolhits-10 && i>=longcolhits-30)
	{
	  downhitcharge += hits_charge[i];   //adding charge on last 20 hits of the track	
	  //std::cout<<i<<" hit_wire: "<<hits_wire[i]<<" "<<"hit_charge: "<<hits_charge[i]<<"\n";
	  downchrg[l]    = hits_charge[i];
	  //downwire[l]    = hits_wire[i];
	  //std::cout<<l<<" hit wire: "<<downwire[l]<<" hit_charge: "<<downchrg[l]<<"\n";	 
	  l++; 
        }    
      }
      //get Truncated mean
      avgdownchrg  =  downhitcharge/l;
      do
      {   
        prevl= l;
	//std::cout<<"avg down charge: "<<avgdownchrg<<"prev l: "<<prevl<<"\n";
	l=0; 
	downhitcharge =0;
        for(int i=0; i<prevl; i++)
	{
          if(downchrg[i]>0.2*avgdownchrg && downchrg[i]<2*avgdownchrg)
	  {
            downhitcharge += downchrg[i];
            downchrg[l]    = downchrg[i];
            l++;
       	  }
        }
	prevavgdownchrg =  avgdownchrg;
	avgdownchrg  =  downhitcharge/l;
	//std::cout<<"prev avg down charge: "<<prevavgdownchrg<<" avg down charge: "<<avgdownchrg<<"\n";		
      }			
      while(avgdownchrg!=prevavgdownchrg);
   			  
      //std::cout<<"upstream avg charge: "<<avgupchrg<<" downstream avg charge: "<<avgdownchrg<<"\n";

      _fPHratio =  avgdownchrg/avgupchrg;
//      fDataTreePH->Fill();
      ///////////////////////////////
      if(avgdownchrg > avgupchrg)	
      {
        PH=1;
      }	
      else 
      {
        PH=-1;
      }
      _fPH = PH;
    std::cout<<"PH = "<<PH<<"\n";  
      	
   //------------------------------MCS Test-------------------------------------------//
    int MCS=-999;
    int counterup =0;
    int wireup=0, wiresqup=0;
    double wirePtimeup=0, timeup=0;
    double slopeup=-999, yinterceptup=-999;
    int countermid =0;
    int wiremid=0, wiresqmid=0;
    double wirePtimemid=0, timemid=0;
    double slopemid=-999, yinterceptmid=-999;
    int counterdown =0;
    int wiredown=0, wiresqdown=0;
    double wirePtimedown=0, timedown=0;
    int midwire = -999;
    double slopedown=-999, yinterceptdown=-999;
    double fittimeup=-999, fittimemid=-999, fittimedown=-999;
      
    for(int i=10; i <30; i++)
    {
      counterup++;
      wirePtimeup += hits_wire[i]*hits_peakT[i];
      timeup += hits_peakT[i];
      wireup += hits_wire[i];
      wiresqup += hits_wire[i]*hits_wire[i];  
    }   
    slopeup = (wirePtimeup*counterup - timeup*wireup)/(wiresqup*counterup - wireup*wireup);
    yinterceptup = (wiresqup*timeup - wirePtimeup*wireup)/(wiresqup*counterup - wireup*wireup);
    
    midwire = (hits_wire[longcolhits-1]+hits_wire[0])/2;
    
    //making sure midwire-10 to midwire+10 and its fits would exist
    for(int i=0; i <longcolhits; i++)
    {
      for(int j=midwire-10; j<midwire+10;j++)
      {
        if(hits_wire[i]==j)
	{ 
	  countermid++;
          wirePtimemid += hits_wire[i]*hits_peakT[i];
          timemid += hits_peakT[i];
          wiremid += hits_wire[i];
          wiresqmid += hits_wire[i]*hits_wire[i]; 
	}  
      }
    }

    slopemid = (wirePtimemid*countermid - timemid*wiremid)/(wiresqmid*countermid - wiremid*wiremid);
    yinterceptmid = (wiresqmid*timemid - wirePtimemid*wiremid)/(wiresqmid*countermid - wiremid*wiremid);     

    for(int i=longcolhits-30; i <longcolhits-10; i++)
    {
      counterdown++;
      wirePtimedown += hits_wire[i]*hits_peakT[i];
      timedown += hits_peakT[i];
      wiredown += hits_wire[i];
      wiresqdown += hits_wire[i]*hits_wire[i];   
    }
    
    slopedown = (wirePtimedown*counterdown - timedown*wiredown)/(wiresqdown*counterdown - wiredown*wiredown);
    yinterceptdown = (wiresqdown*timedown - wirePtimedown*wiredown)/(wiresqdown*counterdown - wiredown*wiredown);

    if(countermid>=5 && TMath::Abs(hits_wire[30]-hits_wire[10])<=40 && TMath::Abs(hits_wire[longcolhits-10]-hits_wire[longcolhits-30])<=40)
    {    
    //Evaluating fits in the middle of the track
    fittimeup = slopeup*midwire + yinterceptup;
    fittimemid = slopemid*midwire + yinterceptmid;
    fittimedown = slopedown*midwire + yinterceptdown;  
    
    //std::cout<<" midwire "<<midwire<<" counterup "<<counterup<<" countermid "<<countermid<<" counterdown "<<counterdown<<"\n";
    //std::cout<<"slopeup "<<slopeup<<" slopemid "<<slopemid<<" slopedown "<<slopedown<<"\n";
    //std::cout<<"fittimeup "<<fittimeup<<" fittimemid "<<fittimemid<<" fittimedown "<<fittimedown<<"\n";
    //std::cout<<"scatter from upstream end: "<<TMath::Abs(fittimeup-fittimemid)<<" scatter from downstream end: "<<TMath::Abs(fittimedown-fittimemid)<<"\n";

    //------------Testing MCS-------------------//
      /////////////plot//////////////////////
      
      _fMCSratio = TMath::Abs(fittimedown-fittimemid)/TMath::Abs(fittimeup-fittimemid);


      //////////////////////////////////////
      if(TMath::Abs(fittimedown-fittimemid) > TMath::Abs(fittimeup-fittimemid))
      {
        MCS=1;
      }
      else 
      {
        MCS=-1;
      }
       
    }//if(countermid=>5 && TMath::Abs(hits_wire[30]-hits_wire[10])<40 && TMath::Abs(hits_wire[longcolhits-10]-hits_wire[longcolhits-30])<40)
    else
    {
      MCS=0;
    }  


    int MCS1=0, MCS2=0;
//    if(_flongtrackmcsfwdll < _flongtrackmcsbwdll) MCS=1;
//    if(_flongtrackmcsfwdll >= _flongtrackmcsbwdll) MCS=-1;
//    std::cout<<"MCS log Likelihood ratio "<<_flongtrackmcsbwdll/_flongtrackmcsfwdll<<"\n";

    _fMCSdiff = (_flongtrackmcsbwdll)-(_flongtrackmcsfwdll);
    _fMCSratio1 = _flongtrackmcsbwdll/_flongtrackmcsfwdll;
  
    if(_fMCSdiff>0) MCS1=1;
    if(_fMCSdiff<=0) MCS1=-1;

    if(_fMCSratio1>1) MCS2=1;
    if(_fMCSratio1<=1) MCS2=-1;

    _fMCS = MCS;
    _fMCS1 = MCS1;
    _fMCS2 = MCS2;
   
  std::cout<<"MCS = "<<MCS<<"\n"; 

  fDataTreePassLongcolhits->Fill();
  
///////////////////////////////////////////////////////////////////////////////

	//-------------Full Final Selected Sample---------------//
	
      if((PH==1 && MCS==1) || (PH==-1 && MCS==1) || (PH==1 && MCS==-1) || (PH==-1 && MCS==-1))
      {

	  fDataTreeFinalSel->Fill();	  
           }//if((PH==1 && MCS==1) || (PH==-1 && MCS==1) || (PH==1 && MCS==-1) || (PH==-1 && MCS==-1))  
          }//if(_fnlongcolhits >= MinColHits)
         }//if(_flongTrackVrtxDis <= Disvrtx)
        }//if(_fstarty_cut == 1)
       }//if(_fbadReg_removed == 1)  
      }//if(TrackCandLength>fMinTrackLen)
     }//if(inFV(trackstartxcandidate, trackstartycandidate, trackstartzcandidate) && inFV(trackendxcandidate, trackendycandidate, trackendzcandidate))
    }//if(flashTrackMatchFlag)
   }//if(TrackCandidate > -1) 
  }////if(flashtag)
 }//if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)   			
   ClearLocalData();
   return;
  } // end CTMDataONAna::analyze()
   	
   DEFINE_ART_MODULE(CTMDataONAna)
 
 } // namespace CTMDataONAna_module
 
  #endif //CTMDATAONANA_H
  

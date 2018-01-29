/**
 * @file    CTMMCAna_module.cc
 * @brief   Module for charged particle multiplicity analysis
 * @Author: Aleena Rafique (aleena@ksu.edu)
 * 
 ******************************************************************************/

  #ifndef CTMMCANA_H
  #define CTMMCANA_H
  
  // Framework includes
#include "cetlib/exception.h"
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
  
  // LArSoft Includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/WireGeo.h"
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
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
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
#include "uboone/EventWeight/MCEventWeight.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubooneobj/UbooneOpticalFilter.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/MCCheater/BackTracker.h"
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

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
  //#include "uboone/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgBase.h"
    constexpr int kMaxGpar = 50000;
    constexpr int kMaxRecoallpar = 50000;
    constexpr int kMaxRecoselpar = 50000;
    constexpr int kMaxtruematchpar = 50000;
    constexpr int kMaxSeconpar = 500000;
    constexpr int kMaxtruepar = 900000;
    constexpr int kMaxtrueInAccppar = 50000;
    constexpr int kMaxtruelongInAccppar = 50000;
  namespace CTMMCAna {  
   
    class CTMMCAna : public art::EDAnalyzer
    {  

    public:
 
     explicit CTMMCAna(fhicl::ParameterSet const &pset);
     virtual ~CTMMCAna();                        
     
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
    int   _fmctrue_origin;    
    int   _fTrueccnc;
    int   _fTruemode;
    int   _fTrueinttype;
    int   _fTruenupdg;
    float _fTrueenu;
    float _fTrueq2truth;
    float _fTruenuvrtxx;
    float _fTruenuvrtxy;
    float _fTruenuvrtxz;
    float _fTrueSPcorrnuvtxx;
    float _fTrueSPcorrnuvtxy;
    float _fTrueSPcorrnuvtxz;
    int   _ftrueVrtxOutFV;
    int   _fbadReg_removed;
    int   _fstarty_cut;
    int   _fbrokenTracks;
    
    int   _fNGTrueallmult;
    int   _fGtrueparID[kMaxGpar];
    int   _fGtrueparpdg[kMaxGpar];
    int   _fGtrueparStatusCode[kMaxGpar];
    float _fGtrueparTheta[kMaxGpar];
    float _fGtrueparCosTheta[kMaxGpar];	
    float _fGtrueparSinTheta[kMaxGpar];
    float _fGtrueparPhi[kMaxGpar];
    float _fGtrueparCosPhi[kMaxGpar];
    float _fGtrueparSinPhi[kMaxGpar];    
    float _fGtrueparE[kMaxGpar];
    float _fGtrueparMass[kMaxGpar];
    float _fGtrueparKE[kMaxGpar];
    float _fGtrueparEndE[kMaxGpar];
    float _fGtrueparPx[kMaxGpar];
    float _fGtrueparPy[kMaxGpar];
    float _fGtrueparPz[kMaxGpar];
    float _fGtrueparP[kMaxGpar];
    float _fGtrueparStartx[kMaxGpar];
    float _fGtrueparStarty[kMaxGpar];
    float _fGtrueparStartz[kMaxGpar];
    float _fGtrueparEndx[kMaxGpar];
    float _fGtrueparEndy[kMaxGpar];
    float _fGtrueparEndz[kMaxGpar];
    float _fGtrueparSPcorrStartx[kMaxGpar];
    float _fGtrueparSPcorrStarty[kMaxGpar];
    float _fGtrueparSPcorrStartz[kMaxGpar];
    float _fGtrueparSPcorrEndx[kMaxGpar];
    float _fGtrueparSPcorrEndy[kMaxGpar];
    float _fGtrueparSPcorrEndz[kMaxGpar];
    
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

    int _fNTruematchPart;
    int _ftruematchparID[kMaxtruematchpar];
    int _ftruematchparpdg[kMaxtruematchpar];
    int _ftruematchparStatusCode[kMaxtruematchpar];
    float _ftruematchparTheta[kMaxtruematchpar];
    float _ftruematchparCosTheta[kMaxtruematchpar];	
    float _ftruematchparSinTheta[kMaxtruematchpar];
    float _ftruematchparPhi[kMaxtruematchpar];
    float _ftruematchparCosPhi[kMaxtruematchpar];
    float _ftruematchparSinPhi[kMaxtruematchpar];    
    float _ftruematchparE[kMaxtruematchpar];
    float _ftruematchparMass[kMaxtruematchpar];
    float _ftruematchparKE[kMaxtruematchpar];
    float _ftruematchparEndE[kMaxtruematchpar];
    float _ftruematchparPx[kMaxtruematchpar];
    float _ftruematchparPy[kMaxtruematchpar];
    float _ftruematchparPz[kMaxtruematchpar];
    float _ftruematchparP[kMaxtruematchpar];
    float _ftruematchparStartx[kMaxtruematchpar];
    float _ftruematchparStarty[kMaxtruematchpar];
    float _ftruematchparStartz[kMaxtruematchpar];
    float _ftruematchparEndx[kMaxtruematchpar];
    float _ftruematchparEndy[kMaxtruematchpar];
    float _ftruematchparEndz[kMaxtruematchpar];
    float _ftruematchparSPcorrStartx[kMaxtruematchpar];
    float _ftruematchparSPcorrStarty[kMaxtruematchpar];
    float _ftruematchparSPcorrStartz[kMaxtruematchpar];
    float _ftruematchparSPcorrEndx[kMaxtruematchpar];
    float _ftruematchparSPcorrEndy[kMaxtruematchpar];
    float _ftruematchparSPcorrEndz[kMaxtruematchpar];

    int _fNSecondaryselmult;
    int _fSeconpargeantID[kMaxSeconpar];
    int _fSeconparpdg[kMaxSeconpar];
    int _fSeconparStatusCode[kMaxSeconpar];
    float _fSeconparTheta[kMaxSeconpar];
    float _fSeconparCosTheta[kMaxSeconpar];	
    float _fSeconparSinTheta[kMaxSeconpar];
    float _fSeconparPhi[kMaxSeconpar];
    float _fSeconparCosPhi[kMaxSeconpar];
    float _fSeconparSinPhi[kMaxSeconpar];    
    float _fSeconparE[kMaxSeconpar];
    float _fSeconparMass[kMaxSeconpar];
    float _fSeconparKE[kMaxSeconpar];
    float _fSeconparEndE[kMaxSeconpar];
    float _fSeconparPx[kMaxSeconpar];
    float _fSeconparPy[kMaxSeconpar];
    float _fSeconparPz[kMaxSeconpar];
    float _fSeconparP[kMaxSeconpar];
    float _fSeconparStartx[kMaxSeconpar];
    float _fSeconparStarty[kMaxSeconpar];
    float _fSeconparStartz[kMaxSeconpar];
    float _fSeconparEndx[kMaxSeconpar];
    float _fSeconparEndy[kMaxSeconpar];
    float _fSeconparEndz[kMaxSeconpar];
    float _fSeconparSPcorrStartx[kMaxSeconpar];
    float _fSeconparSPcorrStarty[kMaxSeconpar];
    float _fSeconparSPcorrStartz[kMaxSeconpar];
    float _fSeconparSPcorrEndx[kMaxSeconpar];
    float _fSeconparSPcorrEndy[kMaxSeconpar];
    float _fSeconparSPcorrEndz[kMaxSeconpar];

    int _fNTrueallPart;
    int _ftrueparID[kMaxtruepar];
    int _ftrueparpdg[kMaxtruepar];
    int _ftrueparStatusCode[kMaxtruepar];
    float _ftrueparTheta[kMaxtruepar];
    float _ftrueparCosTheta[kMaxtruepar];	
    float _ftrueparSinTheta[kMaxtruepar];
    float _ftrueparPhi[kMaxtruepar];
    float _ftrueparCosPhi[kMaxtruepar];
    float _ftrueparSinPhi[kMaxtruepar];    
    float _ftrueparE[kMaxtruepar];
    float _ftrueparMass[kMaxtruepar];
    float _ftrueparKE[kMaxtruepar];
    float _ftrueparEndE[kMaxtruepar];
    float _ftrueparPx[kMaxtruepar];
    float _ftrueparPy[kMaxtruepar];
    float _ftrueparPz[kMaxtruepar];
    float _ftrueparP[kMaxtruepar];
    float _ftrueparStartx[kMaxtruepar];
    float _ftrueparStarty[kMaxtruepar];
    float _ftrueparStartz[kMaxtruepar];
    float _ftrueparEndx[kMaxtruepar];
    float _ftrueparEndy[kMaxtruepar];
    float _ftrueparEndz[kMaxtruepar];
    float _ftrueparSPcorrStartx[kMaxtruepar];
    float _ftrueparSPcorrStarty[kMaxtruepar];
    float _ftrueparSPcorrStartz[kMaxtruepar];
    float _ftrueparSPcorrEndx[kMaxtruepar];
    float _ftrueparSPcorrEndy[kMaxtruepar];
    float _ftrueparSPcorrEndz[kMaxtruepar];

    int _fNtrueInAccpmult;
    int _ftrueInAccppargeantID[kMaxtrueInAccppar];
    int _ftrueInAccpparpdg[kMaxtrueInAccppar];
    int _ftrueInAccpparStatusCode[kMaxtrueInAccppar];
    float _ftrueInAccpparTheta[kMaxtrueInAccppar];
    float _ftrueInAccpparCosTheta[kMaxtrueInAccppar];	
    float _ftrueInAccpparSinTheta[kMaxtrueInAccppar];
    float _ftrueInAccpparPhi[kMaxtrueInAccppar];
    float _ftrueInAccpparCosPhi[kMaxtrueInAccppar];
    float _ftrueInAccpparSinPhi[kMaxtrueInAccppar];    
    float _ftrueInAccpparE[kMaxtrueInAccppar];
    float _ftrueInAccpparMass[kMaxtrueInAccppar];
    float _ftrueInAccpparKE[kMaxtrueInAccppar];
    float _ftrueInAccpparEndE[kMaxtrueInAccppar];
    float _ftrueInAccpparPx[kMaxtrueInAccppar];
    float _ftrueInAccpparPy[kMaxtrueInAccppar];
    float _ftrueInAccpparPz[kMaxtrueInAccppar];
    float _ftrueInAccpparP[kMaxtrueInAccppar];
    float _ftrueInAccpparStartx[kMaxtrueInAccppar];
    float _ftrueInAccpparStarty[kMaxtrueInAccppar];
    float _ftrueInAccpparStartz[kMaxtrueInAccppar];
    float _ftrueInAccpparEndx[kMaxtrueInAccppar];
    float _ftrueInAccpparEndy[kMaxtrueInAccppar];
    float _ftrueInAccpparEndz[kMaxtrueInAccppar];
    float _ftrueInAccpparSPcorrStartx[kMaxtrueInAccppar];
    float _ftrueInAccpparSPcorrStarty[kMaxtrueInAccppar];
    float _ftrueInAccpparSPcorrStartz[kMaxtrueInAccppar];
    float _ftrueInAccpparSPcorrEndx[kMaxtrueInAccppar];
    float _ftrueInAccpparSPcorrEndy[kMaxtrueInAccppar];
    float _ftrueInAccpparSPcorrEndz[kMaxtrueInAccppar];

    int   _fNtrueInAccpmultmu;
    int   _fNtrueInAccpmultpi;
    int   _fNtrueInAccpmultp;
    int   _fNtrueInAccpmultk;

    int _fNtruelongInAccpmult;
    int _ftruelongInAccppargeantID[kMaxtruelongInAccppar];
    int _ftruelongInAccpparpdg[kMaxtruelongInAccppar];
    int _ftruelongInAccpparStatusCode[kMaxtruelongInAccppar];
    float _ftruelongInAccpparTheta[kMaxtruelongInAccppar];
    float _ftruelongInAccpparCosTheta[kMaxtruelongInAccppar];	
    float _ftruelongInAccpparSinTheta[kMaxtruelongInAccppar];
    float _ftruelongInAccpparPhi[kMaxtruelongInAccppar];
    float _ftruelongInAccpparCosPhi[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSinPhi[kMaxtruelongInAccppar];    
    float _ftruelongInAccpparE[kMaxtruelongInAccppar];
    float _ftruelongInAccpparMass[kMaxtruelongInAccppar];
    float _ftruelongInAccpparKE[kMaxtruelongInAccppar];
    float _ftruelongInAccpparEndE[kMaxtruelongInAccppar];
    float _ftruelongInAccpparPx[kMaxtruelongInAccppar];
    float _ftruelongInAccpparPy[kMaxtruelongInAccppar];
    float _ftruelongInAccpparPz[kMaxtruelongInAccppar];
    float _ftruelongInAccpparP[kMaxtruelongInAccppar];
    float _ftruelongInAccpparStartx[kMaxtruelongInAccppar];
    float _ftruelongInAccpparStarty[kMaxtruelongInAccppar];
    float _ftruelongInAccpparStartz[kMaxtruelongInAccppar];
    float _ftruelongInAccpparEndx[kMaxtruelongInAccppar];
    float _ftruelongInAccpparEndy[kMaxtruelongInAccppar];
    float _ftruelongInAccpparEndz[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSPcorrStartx[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSPcorrStarty[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSPcorrStartz[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSPcorrEndx[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSPcorrEndy[kMaxtruelongInAccppar];
    float _ftruelongInAccpparSPcorrEndz[kMaxtruelongInAccppar];

    int   _fPH;
    int   _fMCS;

    float _fPHratio;    
    float _fMCSratio;
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

    const int Xnegbound = 0;
    const int Xposbound = 256;
    const int Ynegbound = -116;
    const int Yposbound = 116;
    const int Znegbound = 0;
    const int Zposbound = 1036;
    
    const int vXnegbound = 5;
    const int vXposbound = 251;
    const int vYnegbound = -106;
    const int vYposbound = 70;
    const int vZnegbound = 5;
    const int vZposbound = 1031;
    
    const int longXnegbound = 5;
    const int longXposbound = 251;
    const int longYnegbound = -106;
    const int longYposbound = 106;
    const int longZnegbound = 5;
    const int longZposbound = 1031;    
        
    const float lengthnumu = 75.0;
           
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

    // For keeping track of the replacement backtracker
    std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;

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
   CTMMCAna::CTMMCAna(fhicl::ParameterSet const& pset): 
   EDAnalyzer(pset),
   fGeometry(lar::providerFrom<geo::Geometry>()),
   fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
      this->reconfigure(pset);
  }
  
    //-----------------------------------------------------------------------
    // Destructor
    CTMMCAna::~CTMMCAna() 
   {
    }
   //-----------------------------------------------------------------------
   void CTMMCAna::ClearLocalData()
   {

  std::fill(_fGtrueparID, _fGtrueparID + sizeof(_fGtrueparID)/sizeof(_fGtrueparID[0]), -99999.); 
  std::fill(_fGtrueparpdg, _fGtrueparpdg + sizeof(_fGtrueparpdg)/sizeof(_fGtrueparpdg[0]), -99999.); 
  std::fill(_fGtrueparStatusCode, _fGtrueparStatusCode + sizeof(_fGtrueparStatusCode)/sizeof(_fGtrueparStatusCode[0]), -99999.);
  std::fill(_fGtrueparTheta, _fGtrueparTheta + sizeof(_fGtrueparTheta)/sizeof(_fGtrueparTheta[0]), -99999.); 
  std::fill(_fGtrueparCosTheta, _fGtrueparCosTheta + sizeof(_fGtrueparCosTheta)/sizeof(_fGtrueparCosTheta[0]), -99999.);     
  std::fill(_fGtrueparSinTheta, _fGtrueparSinTheta + sizeof(_fGtrueparSinTheta)/sizeof(_fGtrueparSinTheta[0]), -99999.); 
  std::fill(_fGtrueparPhi, _fGtrueparPhi + sizeof(_fGtrueparPhi)/sizeof(_fGtrueparPhi[0]), -99999.); 
  std::fill(_fGtrueparCosPhi, _fGtrueparCosPhi + sizeof(_fGtrueparCosPhi)/sizeof(_fGtrueparCosPhi[0]), -99999.);  
  std::fill(_fGtrueparSinPhi, _fGtrueparSinPhi + sizeof(_fGtrueparSinPhi)/sizeof(_fGtrueparSinPhi[0]), -99999.);  
  std::fill(_fGtrueparE, _fGtrueparE+ sizeof(_fGtrueparE)/sizeof(_fGtrueparE[0]), -99999.); 
  std::fill(_fGtrueparMass, _fGtrueparMass + sizeof(_fGtrueparMass)/sizeof(_fGtrueparMass[0]), -99999.); 
  std::fill(_fGtrueparKE, _fGtrueparKE + sizeof(_fGtrueparKE)/sizeof(_fGtrueparKE[0]), -99999.); 
  std::fill(_fGtrueparEndE, _fGtrueparEndE + sizeof(_fGtrueparEndE)/sizeof(_fGtrueparEndE[0]), -99999.); 
  std::fill(_fGtrueparPx, _fGtrueparPx + sizeof(_fGtrueparPx)/sizeof(_fGtrueparPx[0]), -99999.); 
  std::fill(_fGtrueparPy, _fGtrueparPy + sizeof(_fGtrueparPy)/sizeof(_fGtrueparPy[0]), -99999.); 
  std::fill(_fGtrueparPz, _fGtrueparPz + sizeof(_fGtrueparPz)/sizeof(_fGtrueparPz[0]), -99999.); 
  std::fill(_fGtrueparP, _fGtrueparP + sizeof(_fGtrueparP)/sizeof(_fGtrueparP[0]), -99999.); 
  std::fill(_fGtrueparStartx, _fGtrueparStartx + sizeof(_fGtrueparStartx)/sizeof(_fGtrueparStartx[0]), -99999.); 
  std::fill(_fGtrueparStarty, _fGtrueparStarty + sizeof(_fGtrueparStarty)/sizeof(_fGtrueparStarty[0]), -99999.); 
  std::fill(_fGtrueparStartz, _fGtrueparStartz + sizeof(_fGtrueparStartz)/sizeof(_fGtrueparStartz[0]), -99999.); 
  std::fill(_fGtrueparEndx, _fGtrueparEndx + sizeof(_fGtrueparEndx)/sizeof(_fGtrueparEndx[0]), -99999.); 
  std::fill(_fGtrueparEndy, _fGtrueparEndy + sizeof(_fGtrueparEndy)/sizeof(_fGtrueparEndy[0]), -99999.); 
  std::fill(_fGtrueparEndz, _fGtrueparEndz + sizeof(_fGtrueparEndz)/sizeof(_fGtrueparEndz[0]), -99999.); 
  std::fill(_fGtrueparSPcorrStartx, _fGtrueparSPcorrStartx + sizeof(_fGtrueparSPcorrStartx)/sizeof(_fGtrueparSPcorrStartx[0]), -99999.); 
  std::fill(_fGtrueparSPcorrStarty, _fGtrueparSPcorrStarty + sizeof(_fGtrueparSPcorrStarty)/sizeof(_fGtrueparSPcorrStarty[0]), -99999.); 
  std::fill(_fGtrueparSPcorrStartz, _fGtrueparSPcorrStartz + sizeof(_fGtrueparSPcorrStartz)/sizeof(_fGtrueparSPcorrStartz[0]), -99999.); 
  std::fill(_fGtrueparSPcorrEndx, _fGtrueparSPcorrEndx + sizeof(_fGtrueparSPcorrEndx)/sizeof(_fGtrueparSPcorrEndx[0]), -99999.); 
  std::fill(_fGtrueparSPcorrEndy, _fGtrueparSPcorrEndy + sizeof(_fGtrueparSPcorrEndy)/sizeof(_fGtrueparSPcorrEndy[0]), -99999.); 
  std::fill(_fGtrueparSPcorrEndz, _fGtrueparSPcorrEndz + sizeof(_fGtrueparSPcorrEndz)/sizeof(_fGtrueparSPcorrEndz[0]), -99999.); 

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

  std::fill(_ftruematchparID, _ftruematchparID + sizeof(_ftruematchparID)/sizeof(_ftruematchparID[0]), -99999.); 
  std::fill(_ftruematchparpdg, _ftruematchparpdg + sizeof(_ftruematchparpdg)/sizeof(_ftruematchparpdg[0]), -99999.); 
  std::fill(_ftruematchparStatusCode, _ftruematchparStatusCode + sizeof(_ftruematchparStatusCode)/sizeof(_ftruematchparStatusCode[0]), -99999.); 
  std::fill(_ftruematchparTheta, _ftruematchparTheta + sizeof(_ftruematchparTheta)/sizeof(_ftruematchparTheta[0]), -99999.); 
  std::fill(_ftruematchparCosTheta, _ftruematchparCosTheta + sizeof(_ftruematchparCosTheta)/sizeof(_ftruematchparCosTheta[0]), -99999.); 
  std::fill(_ftruematchparSinTheta, _ftruematchparSinTheta + sizeof(_ftruematchparSinTheta)/sizeof(_ftruematchparSinTheta[0]), -99999.); 
  std::fill(_ftruematchparPhi, _ftruematchparPhi + sizeof(_ftruematchparPhi)/sizeof(_ftruematchparPhi[0]), -99999.); 
  std::fill(_ftruematchparCosPhi, _ftruematchparCosPhi + sizeof(_ftruematchparCosPhi)/sizeof(_ftruematchparCosPhi[0]), -99999.); 
  std::fill(_ftruematchparSinPhi, _ftruematchparSinPhi + sizeof(_ftruematchparSinPhi)/sizeof(_ftruematchparSinPhi[0]), -99999.);  
  std::fill(_ftruematchparE, _ftruematchparE + sizeof(_ftruematchparE)/sizeof(_ftruematchparE[0]), -99999.); 
  std::fill(_ftruematchparMass, _ftruematchparMass + sizeof(_ftruematchparMass)/sizeof(_ftruematchparMass[0]), -99999.); 
  std::fill(_ftruematchparKE, _ftruematchparKE + sizeof(_ftruematchparKE)/sizeof(_ftruematchparKE[0]), -99999.); 
  std::fill(_ftruematchparEndE, _ftruematchparEndE + sizeof(_ftruematchparEndE)/sizeof(_ftruematchparEndE[0]), -99999.); 
  std::fill(_ftruematchparP, _ftruematchparP + sizeof(_ftruematchparP)/sizeof(_ftruematchparP[0]), -99999.); 
  std::fill(_ftruematchparPy, _ftruematchparPy + sizeof(_ftruematchparPy)/sizeof(_ftruematchparPy[0]), -99999.); 
  std::fill(_ftruematchparPz, _ftruematchparPz + sizeof(_ftruematchparPz)/sizeof(_ftruematchparPz[0]), -99999.); 
  std::fill(_ftruematchparP, _ftruematchparP + sizeof(_ftruematchparP)/sizeof(_ftruematchparP[0]), -99999.); 
  std::fill(_ftruematchparStartx, _ftruematchparStartx + sizeof(_ftruematchparStartx)/sizeof(_ftruematchparStartx[0]), -99999.); 
  std::fill(_ftruematchparStarty, _ftruematchparStarty + sizeof(_ftruematchparStarty)/sizeof(_ftruematchparStarty[0]), -99999.); 
  std::fill(_ftruematchparStartz, _ftruematchparStartz + sizeof(_ftruematchparStartz)/sizeof(_ftruematchparStartz[0]), -99999.); 
  std::fill(_ftruematchparEndx, _ftruematchparEndx + sizeof(_ftruematchparEndx)/sizeof(_ftruematchparEndx[0]), -99999.); 
  std::fill(_ftruematchparEndy, _ftruematchparEndy + sizeof(_ftruematchparEndy)/sizeof(_ftruematchparEndy[0]), -99999.); 
  std::fill(_ftruematchparEndz, _ftruematchparEndz + sizeof(_ftruematchparEndz)/sizeof(_ftruematchparEndz[0]), -99999.); 
  std::fill(_ftruematchparSPcorrStartx, _ftruematchparSPcorrStartx + sizeof(_ftruematchparSPcorrStartx)/sizeof(_ftruematchparSPcorrStartx[0]), -99999.); 
  std::fill(_ftruematchparSPcorrStarty, _ftruematchparSPcorrStarty + sizeof(_ftruematchparSPcorrStarty)/sizeof(_ftruematchparSPcorrStarty[0]), -99999.); 
  std::fill(_ftruematchparSPcorrStartz, _ftruematchparSPcorrStartz + sizeof(_ftruematchparSPcorrStartz)/sizeof(_ftruematchparSPcorrStartz[0]), -99999.); 
  std::fill(_ftruematchparSPcorrEndx, _ftruematchparSPcorrEndx + sizeof(_ftruematchparSPcorrEndx)/sizeof(_ftruematchparSPcorrEndx[0]), -99999.); 
  std::fill(_ftruematchparSPcorrEndy, _ftruematchparSPcorrEndy + sizeof(_ftruematchparSPcorrEndy)/sizeof(_ftruematchparSPcorrEndy[0]), -99999.); 
  std::fill(_ftruematchparSPcorrEndz, _ftruematchparSPcorrEndz + sizeof(_ftruematchparSPcorrEndz)/sizeof(_ftruematchparSPcorrEndz[0]), -99999.); 

  std::fill(_fSeconpargeantID, _fSeconpargeantID + sizeof(_fSeconpargeantID)/sizeof(_fSeconpargeantID[0]), -99999.); 
  std::fill(_fSeconparpdg, _fSeconparpdg + sizeof(_fSeconparpdg)/sizeof(_fSeconparpdg[0]), -99999.); 
  std::fill(_fSeconparStatusCode, _fSeconparStatusCode + sizeof(_fSeconparStatusCode)/sizeof(_fSeconparStatusCode[0]), -99999.); 
  std::fill(_fSeconparTheta, _fSeconparTheta + sizeof(_fSeconparTheta)/sizeof(_fSeconparTheta[0]), -99999.); 
  std::fill(_fSeconparCosTheta, _fSeconparCosTheta + sizeof(_fSeconparCosTheta)/sizeof(_fSeconparCosTheta[0]), -99999.); 
  std::fill(_fSeconparSinTheta, _fSeconparSinTheta + sizeof(_fSeconparSinTheta)/sizeof(_fSeconparSinTheta[0]), -99999.); 
  std::fill(_fSeconparPhi, _fSeconparPhi + sizeof(_fSeconparPhi)/sizeof(_fSeconparPhi[0]), -99999.); 
  std::fill(_fSeconparCosPhi, _fSeconparCosPhi + sizeof(_fSeconparCosPhi)/sizeof(_fSeconparCosPhi[0]), -99999.); 
  std::fill(_fSeconparSinPhi, _fSeconparSinPhi + sizeof(_fSeconparSinPhi)/sizeof(_fSeconparSinPhi[0]), -99999.);  
  std::fill(_fSeconparE, _fSeconparE + sizeof(_fSeconparE)/sizeof(_fSeconparE[0]), -99999.); 
  std::fill(_fSeconparMass, _fSeconparMass + sizeof(_fSeconparMass)/sizeof(_fSeconparMass[0]), -99999.); 
  std::fill(_fSeconparKE, _fSeconparKE + sizeof(_fSeconparKE)/sizeof(_fSeconparKE[0]), -99999.); 
  std::fill(_fSeconparEndE, _fSeconparEndE + sizeof(_fSeconparEndE)/sizeof(_fSeconparEndE[0]), -99999.); 
  std::fill(_fSeconparP, _fSeconparP + sizeof(_fSeconparP)/sizeof(_fSeconparP[0]), -99999.); 
  std::fill(_fSeconparPy, _fSeconparPy + sizeof(_fSeconparPy)/sizeof(_fSeconparPy[0]), -99999.); 
  std::fill(_fSeconparPz, _fSeconparPz + sizeof(_fSeconparPz)/sizeof(_fSeconparPz[0]), -99999.); 
  std::fill(_fSeconparP, _fSeconparP + sizeof(_fSeconparP)/sizeof(_fSeconparP[0]), -99999.); 
  std::fill(_fSeconparStartx, _fSeconparStartx + sizeof(_fSeconparStartx)/sizeof(_fSeconparStartx[0]), -99999.); 
  std::fill(_fSeconparStarty, _fSeconparStarty + sizeof(_fSeconparStarty)/sizeof(_fSeconparStarty[0]), -99999.); 
  std::fill(_fSeconparStartz, _fSeconparStartz + sizeof(_fSeconparStartz)/sizeof(_fSeconparStartz[0]), -99999.); 
  std::fill(_fSeconparEndx, _fSeconparEndx + sizeof(_fSeconparEndx)/sizeof(_fSeconparEndx[0]), -99999.); 
  std::fill(_fSeconparEndy, _fSeconparEndy + sizeof(_fSeconparEndy)/sizeof(_fSeconparEndy[0]), -99999.); 
  std::fill(_fSeconparEndz, _fSeconparEndz + sizeof(_fSeconparEndz)/sizeof(_fSeconparEndz[0]), -99999.); 
  std::fill(_fSeconparSPcorrStartx, _fSeconparSPcorrStartx + sizeof(_fSeconparSPcorrStartx)/sizeof(_fSeconparSPcorrStartx[0]), -99999.); 
  std::fill(_fSeconparSPcorrStarty, _fSeconparSPcorrStarty + sizeof(_fSeconparSPcorrStarty)/sizeof(_fSeconparSPcorrStarty[0]), -99999.); 
  std::fill(_fSeconparSPcorrStartz, _fSeconparSPcorrStartz + sizeof(_fSeconparSPcorrStartz)/sizeof(_fSeconparSPcorrStartz[0]), -99999.); 
  std::fill(_fSeconparSPcorrEndx, _fSeconparSPcorrEndx + sizeof(_fSeconparSPcorrEndx)/sizeof(_fSeconparSPcorrEndx[0]), -99999.); 
  std::fill(_fSeconparSPcorrEndy, _fSeconparSPcorrEndy + sizeof(_fSeconparSPcorrEndy)/sizeof(_fSeconparSPcorrEndy[0]), -99999.); 
  std::fill(_fSeconparSPcorrEndz, _fSeconparSPcorrEndz + sizeof(_fSeconparSPcorrEndz)/sizeof(_fSeconparSPcorrEndz[0]), -99999.); 

  std::fill(_ftrueparID, _ftrueparID + sizeof(_ftrueparID)/sizeof(_ftrueparID[0]), -99999.); 
  std::fill(_ftrueparpdg, _ftrueparpdg + sizeof(_ftrueparpdg)/sizeof(_ftrueparpdg[0]), -99999.); 
  std::fill(_ftrueparStatusCode, _ftrueparStatusCode + sizeof(_ftrueparStatusCode)/sizeof(_ftrueparStatusCode[0]), -99999.); 
  std::fill(_ftrueparTheta, _ftrueparTheta + sizeof(_ftrueparTheta)/sizeof(_ftrueparTheta[0]), -99999.); 
  std::fill(_ftrueparCosTheta, _ftrueparCosTheta + sizeof(_ftrueparCosTheta)/sizeof(_ftrueparCosTheta[0]), -99999.); 
  std::fill(_ftrueparSinTheta, _ftrueparSinTheta + sizeof(_ftrueparSinTheta)/sizeof(_ftrueparSinTheta[0]), -99999.); 
  std::fill(_ftrueparPhi, _ftrueparPhi + sizeof(_ftrueparPhi)/sizeof(_ftrueparPhi[0]), -99999.); 
  std::fill(_ftrueparCosPhi, _ftrueparCosPhi + sizeof(_ftrueparCosPhi)/sizeof(_ftrueparCosPhi[0]), -99999.); 
  std::fill(_ftrueparSinPhi, _ftrueparSinPhi + sizeof(_ftrueparSinPhi)/sizeof(_ftrueparSinPhi[0]), -99999.);  
  std::fill(_ftrueparE, _ftrueparE + sizeof(_ftrueparE)/sizeof(_ftrueparE[0]), -99999.); 
  std::fill(_ftrueparMass, _ftrueparMass + sizeof(_ftrueparMass)/sizeof(_ftrueparMass[0]), -99999.); 
  std::fill(_ftrueparKE, _ftrueparKE + sizeof(_ftrueparKE)/sizeof(_ftrueparKE[0]), -99999.); 
  std::fill(_ftrueparEndE, _ftrueparEndE + sizeof(_ftrueparEndE)/sizeof(_ftrueparEndE[0]), -99999.); 
  std::fill(_ftrueparP, _ftrueparP + sizeof(_ftrueparP)/sizeof(_ftrueparP[0]), -99999.); 
  std::fill(_ftrueparPy, _ftrueparPy + sizeof(_ftrueparPy)/sizeof(_ftrueparPy[0]), -99999.); 
  std::fill(_ftrueparPz, _ftrueparPz + sizeof(_ftrueparPz)/sizeof(_ftrueparPz[0]), -99999.); 
  std::fill(_ftrueparP, _ftrueparP + sizeof(_ftrueparP)/sizeof(_ftrueparP[0]), -99999.); 
  std::fill(_ftrueparStartx, _ftrueparStartx + sizeof(_ftrueparStartx)/sizeof(_ftrueparStartx[0]), -99999.); 
  std::fill(_ftrueparStarty, _ftrueparStarty + sizeof(_ftrueparStarty)/sizeof(_ftrueparStarty[0]), -99999.); 
  std::fill(_ftrueparStartz, _ftrueparStartz + sizeof(_ftrueparStartz)/sizeof(_ftrueparStartz[0]), -99999.); 
  std::fill(_ftrueparEndx, _ftrueparEndx + sizeof(_ftrueparEndx)/sizeof(_ftrueparEndx[0]), -99999.); 
  std::fill(_ftrueparEndy, _ftrueparEndy + sizeof(_ftrueparEndy)/sizeof(_ftrueparEndy[0]), -99999.); 
  std::fill(_ftrueparEndz, _ftrueparEndz + sizeof(_ftrueparEndz)/sizeof(_ftrueparEndz[0]), -99999.); 
  std::fill(_ftrueparSPcorrStartx, _ftrueparSPcorrStartx + sizeof(_ftrueparSPcorrStartx)/sizeof(_ftrueparSPcorrStartx[0]), -99999.); 
  std::fill(_ftrueparSPcorrStarty, _ftrueparSPcorrStarty + sizeof(_ftrueparSPcorrStarty)/sizeof(_ftrueparSPcorrStarty[0]), -99999.); 
  std::fill(_ftrueparSPcorrStartz, _ftrueparSPcorrStartz + sizeof(_ftrueparSPcorrStartz)/sizeof(_ftrueparSPcorrStartz[0]), -99999.); 
  std::fill(_ftrueparSPcorrEndx, _ftrueparSPcorrEndx + sizeof(_ftrueparSPcorrEndx)/sizeof(_ftrueparSPcorrEndx[0]), -99999.); 
  std::fill(_ftrueparSPcorrEndy, _ftrueparSPcorrEndy + sizeof(_ftrueparSPcorrEndy)/sizeof(_ftrueparSPcorrEndy[0]), -99999.); 
  std::fill(_ftrueparSPcorrEndz, _ftrueparSPcorrEndz + sizeof(_ftrueparSPcorrEndz)/sizeof(_ftrueparSPcorrEndz[0]), -99999.); 

  std::fill(_ftrueInAccppargeantID, _ftrueInAccppargeantID + sizeof(_ftrueInAccppargeantID)/sizeof(_ftrueInAccppargeantID[0]), -99999.); 
  std::fill(_ftrueInAccpparpdg, _ftrueInAccpparpdg + sizeof(_ftrueInAccpparpdg)/sizeof(_ftrueInAccpparpdg[0]), -99999.); 
  std::fill(_ftrueInAccpparStatusCode, _ftrueInAccpparStatusCode + sizeof(_ftrueInAccpparStatusCode)/sizeof(_ftrueInAccpparStatusCode[0]), -99999.); 
  std::fill(_ftrueInAccpparTheta, _ftrueInAccpparTheta + sizeof(_ftrueInAccpparTheta)/sizeof(_ftrueInAccpparTheta[0]), -99999.); 
  std::fill(_ftrueInAccpparCosTheta, _ftrueInAccpparCosTheta + sizeof(_ftrueInAccpparCosTheta)/sizeof(_ftrueInAccpparCosTheta[0]), -99999.); 
  std::fill(_ftrueInAccpparSinTheta, _ftrueInAccpparSinTheta + sizeof(_ftrueInAccpparSinTheta)/sizeof(_ftrueInAccpparSinTheta[0]), -99999.); 
  std::fill(_ftrueInAccpparPhi, _ftrueInAccpparPhi + sizeof(_ftrueInAccpparPhi)/sizeof(_ftrueInAccpparPhi[0]), -99999.); 
  std::fill(_ftrueInAccpparCosPhi, _ftrueInAccpparCosPhi + sizeof(_ftrueInAccpparCosPhi)/sizeof(_ftrueInAccpparCosPhi[0]), -99999.); 
  std::fill(_ftrueInAccpparSinPhi, _ftrueInAccpparSinPhi + sizeof(_ftrueInAccpparSinPhi)/sizeof(_ftrueInAccpparSinPhi[0]), -99999.);  
  std::fill(_ftrueInAccpparE, _ftrueInAccpparE + sizeof(_ftrueInAccpparE)/sizeof(_ftrueInAccpparE[0]), -99999.); 
  std::fill(_ftrueInAccpparMass, _ftrueInAccpparMass + sizeof(_ftrueInAccpparMass)/sizeof(_ftrueInAccpparMass[0]), -99999.); 
  std::fill(_ftrueInAccpparKE, _ftrueInAccpparKE + sizeof(_ftrueInAccpparKE)/sizeof(_ftrueInAccpparKE[0]), -99999.); 
  std::fill(_ftrueInAccpparEndE, _ftrueInAccpparEndE + sizeof(_ftrueInAccpparEndE)/sizeof(_ftrueInAccpparEndE[0]), -99999.); 
  std::fill(_ftrueInAccpparP, _ftrueInAccpparP + sizeof(_ftrueInAccpparP)/sizeof(_ftrueInAccpparP[0]), -99999.); 
  std::fill(_ftrueInAccpparPy, _ftrueInAccpparPy + sizeof(_ftrueInAccpparPy)/sizeof(_ftrueInAccpparPy[0]), -99999.); 
  std::fill(_ftrueInAccpparPz, _ftrueInAccpparPz + sizeof(_ftrueInAccpparPz)/sizeof(_ftrueInAccpparPz[0]), -99999.); 
  std::fill(_ftrueInAccpparP, _ftrueInAccpparP + sizeof(_ftrueInAccpparP)/sizeof(_ftrueInAccpparP[0]), -99999.); 
  std::fill(_ftrueInAccpparStartx, _ftrueInAccpparStartx + sizeof(_ftrueInAccpparStartx)/sizeof(_ftrueInAccpparStartx[0]), -99999.); 
  std::fill(_ftrueInAccpparStarty, _ftrueInAccpparStarty + sizeof(_ftrueInAccpparStarty)/sizeof(_ftrueInAccpparStarty[0]), -99999.); 
  std::fill(_ftrueInAccpparStartz, _ftrueInAccpparStartz + sizeof(_ftrueInAccpparStartz)/sizeof(_ftrueInAccpparStartz[0]), -99999.); 
  std::fill(_ftrueInAccpparEndx, _ftrueInAccpparEndx + sizeof(_ftrueInAccpparEndx)/sizeof(_ftrueInAccpparEndx[0]), -99999.); 
  std::fill(_ftrueInAccpparEndy, _ftrueInAccpparEndy + sizeof(_ftrueInAccpparEndy)/sizeof(_ftrueInAccpparEndy[0]), -99999.); 
  std::fill(_ftrueInAccpparEndz, _ftrueInAccpparEndz + sizeof(_ftrueInAccpparEndz)/sizeof(_ftrueInAccpparEndz[0]), -99999.); 
  std::fill(_ftrueInAccpparSPcorrStartx, _ftrueInAccpparSPcorrStartx + sizeof(_ftrueInAccpparSPcorrStartx)/sizeof(_ftrueInAccpparSPcorrStartx[0]), -99999.); 
  std::fill(_ftrueInAccpparSPcorrStarty, _ftrueInAccpparSPcorrStarty + sizeof(_ftrueInAccpparSPcorrStarty)/sizeof(_ftrueInAccpparSPcorrStarty[0]), -99999.); 
  std::fill(_ftrueInAccpparSPcorrStartz, _ftrueInAccpparSPcorrStartz + sizeof(_ftrueInAccpparSPcorrStartz)/sizeof(_ftrueInAccpparSPcorrStartz[0]), -99999.); 
  std::fill(_ftrueInAccpparSPcorrEndx, _ftrueInAccpparSPcorrEndx + sizeof(_ftrueInAccpparSPcorrEndx)/sizeof(_ftrueInAccpparSPcorrEndx[0]), -99999.); 
  std::fill(_ftrueInAccpparSPcorrEndy, _ftrueInAccpparSPcorrEndy + sizeof(_ftrueInAccpparSPcorrEndy)/sizeof(_ftrueInAccpparSPcorrEndy[0]), -99999.); 
  std::fill(_ftrueInAccpparSPcorrEndz, _ftrueInAccpparSPcorrEndz + sizeof(_ftrueInAccpparSPcorrEndz)/sizeof(_ftrueInAccpparSPcorrEndz[0]), -99999.); 

  std::fill(_ftruelongInAccppargeantID, _ftruelongInAccppargeantID + sizeof(_ftruelongInAccppargeantID)/sizeof(_ftruelongInAccppargeantID[0]), -99999.); 
  std::fill(_ftruelongInAccpparpdg, _ftruelongInAccpparpdg + sizeof(_ftruelongInAccpparpdg)/sizeof(_ftruelongInAccpparpdg[0]), -99999.); 
  std::fill(_ftruelongInAccpparStatusCode, _ftruelongInAccpparStatusCode + sizeof(_ftruelongInAccpparStatusCode)/sizeof(_ftruelongInAccpparStatusCode[0]), -99999.); 
  std::fill(_ftruelongInAccpparTheta, _ftruelongInAccpparTheta + sizeof(_ftruelongInAccpparTheta)/sizeof(_ftruelongInAccpparTheta[0]), -99999.); 
  std::fill(_ftruelongInAccpparCosTheta, _ftruelongInAccpparCosTheta + sizeof(_ftruelongInAccpparCosTheta)/sizeof(_ftruelongInAccpparCosTheta[0]), -99999.); 
  std::fill(_ftruelongInAccpparSinTheta, _ftruelongInAccpparSinTheta + sizeof(_ftruelongInAccpparSinTheta)/sizeof(_ftruelongInAccpparSinTheta[0]), -99999.); 
  std::fill(_ftruelongInAccpparPhi, _ftruelongInAccpparPhi + sizeof(_ftruelongInAccpparPhi)/sizeof(_ftruelongInAccpparPhi[0]), -99999.); 
  std::fill(_ftruelongInAccpparCosPhi, _ftruelongInAccpparCosPhi + sizeof(_ftruelongInAccpparCosPhi)/sizeof(_ftruelongInAccpparCosPhi[0]), -99999.); 
  std::fill(_ftruelongInAccpparSinPhi, _ftruelongInAccpparSinPhi + sizeof(_ftruelongInAccpparSinPhi)/sizeof(_ftruelongInAccpparSinPhi[0]), -99999.);  
  std::fill(_ftruelongInAccpparE, _ftruelongInAccpparE + sizeof(_ftruelongInAccpparE)/sizeof(_ftruelongInAccpparE[0]), -99999.); 
  std::fill(_ftruelongInAccpparMass, _ftruelongInAccpparMass + sizeof(_ftruelongInAccpparMass)/sizeof(_ftruelongInAccpparMass[0]), -99999.); 
  std::fill(_ftruelongInAccpparKE, _ftruelongInAccpparKE + sizeof(_ftruelongInAccpparKE)/sizeof(_ftruelongInAccpparKE[0]), -99999.); 
  std::fill(_ftruelongInAccpparEndE, _ftruelongInAccpparEndE + sizeof(_ftruelongInAccpparEndE)/sizeof(_ftruelongInAccpparEndE[0]), -99999.); 
  std::fill(_ftruelongInAccpparP, _ftruelongInAccpparP + sizeof(_ftruelongInAccpparP)/sizeof(_ftruelongInAccpparP[0]), -99999.); 
  std::fill(_ftruelongInAccpparPy, _ftruelongInAccpparPy + sizeof(_ftruelongInAccpparPy)/sizeof(_ftruelongInAccpparPy[0]), -99999.); 
  std::fill(_ftruelongInAccpparPz, _ftruelongInAccpparPz + sizeof(_ftruelongInAccpparPz)/sizeof(_ftruelongInAccpparPz[0]), -99999.); 
  std::fill(_ftruelongInAccpparP, _ftruelongInAccpparP + sizeof(_ftruelongInAccpparP)/sizeof(_ftruelongInAccpparP[0]), -99999.); 
  std::fill(_ftruelongInAccpparStartx, _ftruelongInAccpparStartx + sizeof(_ftruelongInAccpparStartx)/sizeof(_ftruelongInAccpparStartx[0]), -99999.); 
  std::fill(_ftruelongInAccpparStarty, _ftruelongInAccpparStarty + sizeof(_ftruelongInAccpparStarty)/sizeof(_ftruelongInAccpparStarty[0]), -99999.); 
  std::fill(_ftruelongInAccpparStartz, _ftruelongInAccpparStartz + sizeof(_ftruelongInAccpparStartz)/sizeof(_ftruelongInAccpparStartz[0]), -99999.); 
  std::fill(_ftruelongInAccpparEndx, _ftruelongInAccpparEndx + sizeof(_ftruelongInAccpparEndx)/sizeof(_ftruelongInAccpparEndx[0]), -99999.); 
  std::fill(_ftruelongInAccpparEndy, _ftruelongInAccpparEndy + sizeof(_ftruelongInAccpparEndy)/sizeof(_ftruelongInAccpparEndy[0]), -99999.); 
  std::fill(_ftruelongInAccpparEndz, _ftruelongInAccpparEndz + sizeof(_ftruelongInAccpparEndz)/sizeof(_ftruelongInAccpparEndz[0]), -99999.); 
  std::fill(_ftruelongInAccpparSPcorrStartx, _ftruelongInAccpparSPcorrStartx + sizeof(_ftruelongInAccpparSPcorrStartx)/sizeof(_ftruelongInAccpparSPcorrStartx[0]), -99999.); 
  std::fill(_ftruelongInAccpparSPcorrStarty, _ftruelongInAccpparSPcorrStarty + sizeof(_ftruelongInAccpparSPcorrStarty)/sizeof(_ftruelongInAccpparSPcorrStarty[0]), -99999.); 
  std::fill(_ftruelongInAccpparSPcorrStartz, _ftruelongInAccpparSPcorrStartz + sizeof(_ftruelongInAccpparSPcorrStartz)/sizeof(_ftruelongInAccpparSPcorrStartz[0]), -99999.); 
  std::fill(_ftruelongInAccpparSPcorrEndx, _ftruelongInAccpparSPcorrEndx + sizeof(_ftruelongInAccpparSPcorrEndx)/sizeof(_ftruelongInAccpparSPcorrEndx[0]), -99999.); 
  std::fill(_ftruelongInAccpparSPcorrEndy, _ftruelongInAccpparSPcorrEndy + sizeof(_ftruelongInAccpparSPcorrEndy)/sizeof(_ftruelongInAccpparSPcorrEndy[0]), -99999.); 
  std::fill(_ftruelongInAccpparSPcorrEndz, _ftruelongInAccpparSPcorrEndz + sizeof(_ftruelongInAccpparSPcorrEndz)/sizeof(_ftruelongInAccpparSPcorrEndz[0]), -99999.); 

   }
   //-----------------------------------------------------------------------
   void CTMMCAna::beginJob()
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
    fDataTreeCCSel->Branch("_fTrueccnc",&_fTrueccnc,"_fTrueccnc/I");
    fDataTreeCCSel->Branch("_fTruemode",&_fTruemode,"_fTruemode/I");
    fDataTreeCCSel->Branch("_fTrueinttype",&_fTrueinttype,"_fTrueinttype/I");
    fDataTreeCCSel->Branch("_fTruenupdg",&_fTruenupdg,"_fTruenupdg/I");
    fDataTreeCCSel->Branch("_fTrueenu",&_fTrueenu,"_fTrueenu/F");
    fDataTreeCCSel->Branch("_fTrueq2truth",&_fTrueq2truth,"_fTrueq2truth/F");
    fDataTreeCCSel->Branch("_fTruenuvrtxx",&_fTruenuvrtxx,"_fTruenuvrtxx/F");
    fDataTreeCCSel->Branch("_fTruenuvrtxy",&_fTruenuvrtxy,"_fTruenuvrtxy/F");
    fDataTreeCCSel->Branch("_fTruenuvrtxz",&_fTruenuvrtxz,"_fTruenuvrtxz/F");
    fDataTreeCCSel->Branch("_fTrueSPcorrnuvtxx",&_fTrueSPcorrnuvtxx,"_fTrueSPcorrnuvtxx/F");
    fDataTreeCCSel->Branch("_fTrueSPcorrnuvtxy",&_fTrueSPcorrnuvtxy,"_fTrueSPcorrnuvtxy/F");
    fDataTreeCCSel->Branch("_fTrueSPcorrnuvtxz",&_fTrueSPcorrnuvtxz,"_fTrueSPcorrnuvtxz/F");
    fDataTreeCCSel->Branch("_fmctrue_origin",&_fmctrue_origin,"_fmctrue_origin/I");
    fDataTreeCCSel->Branch("_ftrueVrtxOutFV",&_ftrueVrtxOutFV,"_ftrueVrtxOutFV/I");
    fDataTreeCCSel->Branch("_fbadReg_removed",&_fbadReg_removed,"_fbadReg_removed/I");
    fDataTreeCCSel->Branch("_fstarty_cut",&_fstarty_cut,"_fstarty_cut/I");
    fDataTreeCCSel->Branch("_fbrokenTracks",&_fbrokenTracks,"_fbrokenTracks/I");
    fDataTreeCCSel->Branch("_fNGTrueallmult",&_fNGTrueallmult,"_fNGTrueallmult/I");   
    fDataTreeCCSel->Branch("_fGtrueparID",&_fGtrueparID,"_fGtrueparID[_fNGTrueallmult]/I");
    fDataTreeCCSel->Branch("_fGtrueparpdg",&_fGtrueparpdg,"_fGtrueparpdg[_fNGTrueallmult]/I");
    fDataTreeCCSel->Branch("_fGtrueparStatusCode",&_fGtrueparStatusCode,"_fGtrueparStatusCode[_fNGTrueallmult]/I");
    fDataTreeCCSel->Branch("_fGtrueparTheta",&_fGtrueparTheta,"_fGtrueparTheta[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparCosTheta",&_fGtrueparCosTheta,"_fGtrueparCosTheta[_fNGTrueallmult]/F");	
    fDataTreeCCSel->Branch("_fGtrueparSinTheta",&_fGtrueparSinTheta,"_fGtrueparSinTheta[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparPhi",&_fGtrueparPhi,"_fGtrueparPhi[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparCosPhi",&_fGtrueparCosPhi,"_fGtrueparCosPhi[_fNGTrueallmult]/F");    
    fDataTreeCCSel->Branch("_fGtrueparSinPhi",&_fGtrueparSinPhi,"_fGtrueparSinPhi[_fNGTrueallmult]/F");    
    fDataTreeCCSel->Branch("_fGtrueparE",&_fGtrueparE,"_fGtrueparE[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparMass",&_fGtrueparMass,"_fGtrueparMass[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparKE",&_fGtrueparKE,"_fGtrueparKE[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparEndE",&_fGtrueparEndE,"_fGtrueparEndE[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparPx",&_fGtrueparPx,"_fGtrueparPx[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparPy",&_fGtrueparPy,"_fGtrueparPy[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparPz",&_fGtrueparPz,"_fGtrueparPz[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparP",&_fGtrueparP,"_fGtrueparP[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparStartx",&_fGtrueparStartx,"_fGtrueparStartx[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparStarty",&_fGtrueparStarty,"_fGtrueparStarty[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparStartz",&_fGtrueparStartz,"_fGtrueparStartz[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparEndx",&_fGtrueparEndx,"_fGtrueparEndx[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparEndy",&_fGtrueparEndy,"_fGtrueparEndy[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparEndz",&_fGtrueparEndz,"_fGtrueparEndz[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparSPcorrStartx",&_fGtrueparSPcorrStartx,"_fGtrueparSPcorrStartx[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparSPcorrStarty",&_fGtrueparSPcorrStarty,"_fGtrueparSPcorrStarty[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparSPcorrStartz",&_fGtrueparSPcorrStartz,"_fGtrueparSPcorrStartz[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparSPcorrEndx",&_fGtrueparSPcorrEndx,"_fGtrueparSPcorrEndx[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparSPcorrEndy",&_fGtrueparSPcorrEndy,"_fGtrueparSPcorrEndy[_fNGTrueallmult]/F");
    fDataTreeCCSel->Branch("_fGtrueparSPcorrEndz",&_fGtrueparSPcorrEndz,"_fGtrueparSPcorrEndz[_fNGTrueallmult]/F");
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
    fDataTreeCCSel->Branch("_fNTruematchPart",&_fNTruematchPart,"_fNTruematchPart/I");   
    fDataTreeCCSel->Branch("_ftruematchparID",&_ftruematchparID,"_ftruematchparID[_fNTruematchPart]/I");
    fDataTreeCCSel->Branch("_ftruematchparpdg",&_ftruematchparpdg,"_ftruematchparpdg[_fNTruematchPart]/I");
    fDataTreeCCSel->Branch("_ftruematchparStatusCode",&_ftruematchparStatusCode,"_ftruematchparStatusCode[_fNTruematchPart]/I");
    fDataTreeCCSel->Branch("_ftruematchparTheta",&_ftruematchparTheta,"_ftruematchparTheta[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparCosTheta",&_ftruematchparCosTheta,"_ftruematchparCosTheta[_fNTruematchPart]/F");	
    fDataTreeCCSel->Branch("_ftruematchparSinTheta",&_ftruematchparSinTheta,"_ftruematchparSinTheta[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparPhi",&_ftruematchparPhi,"_ftruematchparPhi[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparCosPhi",&_ftruematchparCosPhi,"_ftruematchparCosPhi[_fNTruematchPart]/F");    
    fDataTreeCCSel->Branch("_ftruematchparSinPhi",&_ftruematchparSinPhi,"_ftruematchparSinPhi[_fNTruematchPart]/F");	 
    fDataTreeCCSel->Branch("_ftruematchparE",&_ftruematchparE,"_ftruematchparE[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparMass",&_ftruematchparMass,"_ftruematchparMass[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparKE",&_ftruematchparKE,"_ftruematchparKE[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparEndE",&_ftruematchparEndE,"_ftruematchparEndE[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparPx",&_ftruematchparPx,"_ftruematchparPx[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparPy",&_ftruematchparPy,"_ftruematchparPy[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparPz",&_ftruematchparPz,"_ftruematchparPz[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparP",&_ftruematchparP,"_fTruematchparP[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparStartx",&_ftruematchparStartx,"_ftruematchparStartx[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparStarty",&_ftruematchparStarty,"_ftruematchparStarty[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparStartz",&_ftruematchparStartz,"_ftruematchparStartz[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparEndx",&_ftruematchparEndx,"_ftruematchparEndx[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparEndy",&_ftruematchparEndy,"_ftruematchparEndy[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparEndz",&_ftruematchparEndz,"_ftruematchparEndz[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparSPcorrStartx",&_ftruematchparSPcorrStartx,"_ftruematchparSPcorrStartx[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparSPcorrStarty",&_ftruematchparSPcorrStarty,"_ftruematchparSPcorrStarty[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparSPcorrStartz",&_ftruematchparSPcorrStartz,"_ftruematchparSPcorrStartz[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparSPcorrEndx",&_ftruematchparSPcorrEndx,"_ftruematchparSPcorrEndx[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparSPcorrEndy",&_ftruematchparSPcorrEndy,"_ftruematchparSPcorrEndy[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_ftruematchparSPcorrEndz",&_ftruematchparSPcorrEndz,"_ftruematchparSPcorrEndz[_fNTruematchPart]/F");
    fDataTreeCCSel->Branch("_fNSecondaryselmult",&_fNSecondaryselmult,"_fNSecondaryselmult/I");
    fDataTreeCCSel->Branch("_fSeconpargeantID",&_fSeconpargeantID,"_fSeconpargeantID[_fNSecondaryselmult]/I");
    fDataTreeCCSel->Branch("_fSeconparpdg",&_fSeconparpdg,"_fSeconparpdg[_fNSecondaryselmult]/I");
    fDataTreeCCSel->Branch("_fSeconparStatusCode",&_fSeconparStatusCode,"_fSeconparStatusCode[_fNSecondaryselmult]/I");
    fDataTreeCCSel->Branch("_fSeconparTheta",&_fSeconparTheta,"_fSeconparTheta[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparCosTheta",&_fSeconparCosTheta,"_fSeconparCosTheta[_fNSecondaryselmult]/F");    
    fDataTreeCCSel->Branch("_fSeconparSinTheta",&_fSeconparSinTheta,"_fSeconparSinTheta[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparPhi",&_fSeconparPhi,"_fSeconparPhi[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparCosPhi",&_fSeconparCosPhi,"_fSeconparCosPhi[_fNSecondaryselmult]/F");	 
    fDataTreeCCSel->Branch("_fSeconparSinPhi",&_fSeconparSinPhi,"_fSeconparSinPhi[_fNSecondaryselmult]/F");    
    fDataTreeCCSel->Branch("_fSeconparE",&_fSeconparE,"_fSeconparE[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparMass",&_fSeconparMass,"_fSeconparMass[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparKE",&_fSeconparKE,"_fSeconparKE[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparEndE",&_fSeconparEndE,"_fSeconparEndE[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparPx",&_fSeconparPx,"_fSeconparPx[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparPy",&_fSeconparPy,"_fSeconparPy[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparPz",&_fSeconparPz,"_fSeconparPz[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparP",&_fSeconparP,"_fSeconparP[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparStartx",&_fSeconparStartx,"_fSeconparStartx[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparStarty",&_fSeconparStarty,"_fSeconparStarty[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparStartz",&_fSeconparStartz,"_fSeconparStartz[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparEndx",&_fSeconparEndx,"_fSeconparEndx[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparEndy",&_fSeconparEndy,"_fSeconparEndy[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparEndz",&_fSeconparEndz,"_fSeconparEndz[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparSPcorrStartx",&_fSeconparSPcorrStartx,"_fSeconparSPcorrStartx[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparSPcorrStarty",&_fSeconparSPcorrStarty,"_fSeconparSPcorrStarty[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparSPcorrStartz",&_fSeconparSPcorrStartz,"_fSeconparSPcorrStartz[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparSPcorrEndx",&_fSeconparSPcorrEndx,"_fSeconparSPcorrEndx[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparSPcorrEndy",&_fSeconparSPcorrEndy,"_fSeconparSPcorrEndy[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fSeconparSPcorrEndz",&_fSeconparSPcorrEndz,"_fSeconparSPcorrEndz[_fNSecondaryselmult]/F");
    fDataTreeCCSel->Branch("_fNTrueallPart",&_fNTrueallPart,"_fNTrueallPart/I");   
/*  fDataTreeFinalSel->Branch("_ftrueparID",&_ftrueparID,"_ftrueparID[_fNTrueallPart]/I");
    fDataTreeCCSel->Branch("_ftrueparpdg",&_ftrueparpdg,"_ftrueparpdg[_fNTrueallPart]/I");
    fDataTreeCCSel->Branch("_ftrueparStatusCode",&_ftrueparStatusCode,"_ftrueparStatusCode[_fNTrueallPart]/I");
    fDataTreeCCSel->Branch("_ftrueparTheta",&_ftrueparTheta,"_ftrueparTheta[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparCosTheta",&_ftrueparCosTheta,"_ftrueparCosTheta[_fNTrueallPart]/F");	
    fDataTreeCCSel->Branch("_ftrueparSinTheta",&_ftrueparSinTheta,"_ftrueparSinTheta[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparPhi",&_ftrueparPhi,"_ftrueparPhi[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparCosPhi",&_ftrueparCosPhi,"_ftrueparCosPhi[_fNTrueallPart]/F");    
    fDataTreeCCSel->Branch("_ftrueparSinPhi",&_ftrueparSinPhi,"_ftrueparSinPhi[_fNTrueallPart]/F");    
    fDataTreeCCSel->Branch("_ftrueparE",&_ftrueparE,"_ftrueparE[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparMass",&_ftrueparMass,"_ftrueparMass[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparKE",&_ftrueparKE,"_ftrueparKE[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparEndE",&_ftrueparEndE,"_ftrueparEndE[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparPx",&_ftrueparPx,"_ftrueparPx[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparPy",&_ftrueparPy,"_ftrueparPy[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparPz",&_ftrueparPz,"_ftrueparPz[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparP",&_ftrueparP,"_fTrueparP[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparStartx",&_ftrueparStartx,"_ftrueparStartx[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparStarty",&_ftrueparStarty,"_ftrueparStarty[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparStartz",&_ftrueparStartz,"_ftrueparStartz[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparEndx",&_ftrueparEndx,"_ftrueparEndx[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparEndy",&_ftrueparEndy,"_ftrueparEndy[_fNTrueallPart]/F");
    fDataTreeCCSel->Branch("_ftrueparEndz",&_ftrueparEndz,"_ftrueparEndz[_fNTrueallPart]/F");
*/    fDataTreeCCSel->Branch("_fNtrueInAccpmult",&_fNtrueInAccpmult,"_fNtrueInAccpmult/I");
    fDataTreeCCSel->Branch("_ftrueInAccppargeantID",&_ftrueInAccppargeantID,"_ftrueInAccppargeantID[_fNtrueInAccpmult]/I");
    fDataTreeCCSel->Branch("_ftrueInAccpparpdg",&_ftrueInAccpparpdg,"_ftrueInAccpparpdg[_fNtrueInAccpmult]/I");
    fDataTreeCCSel->Branch("_ftrueInAccpparStatusCode",&_ftrueInAccpparStatusCode,"_ftrueInAccpparStatusCode[_fNtrueInAccpmult]/I");
    fDataTreeCCSel->Branch("_ftrueInAccpparTheta",&_ftrueInAccpparTheta,"_ftrueInAccpparTheta[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparCosTheta",&_ftrueInAccpparCosTheta,"_ftrueInAccpparCosTheta[_fNtrueInAccpmult]/F");	
    fDataTreeCCSel->Branch("_ftrueInAccpparSinTheta",&_ftrueInAccpparSinTheta,"_ftrueInAccpparSinTheta[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparPhi",&_ftrueInAccpparPhi,"_ftrueInAccpparPhi[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparCosPhi",&_ftrueInAccpparCosPhi,"_ftrueInAccpparCosPhi[_fNtrueInAccpmult]/F");    
    fDataTreeCCSel->Branch("_ftrueInAccpparSinPhi",&_ftrueInAccpparSinPhi,"_ftrueInAccpparSinPhi[_fNtrueInAccpmult]/F");    
    fDataTreeCCSel->Branch("_ftrueInAccpparE",&_ftrueInAccpparE,"_ftrueInAccpparE[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparMass",&_ftrueInAccpparMass,"_ftrueInAccpparMass[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparKE",&_ftrueInAccpparKE,"_ftrueInAccpparKE[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparEndE",&_ftrueInAccpparEndE,"_ftrueInAccpparEndE[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparPx",&_ftrueInAccpparPx,"_ftrueInAccpparPx[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparPy",&_ftrueInAccpparPy,"_ftrueInAccpparPy[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparPz",&_ftrueInAccpparPz,"_ftrueInAccpparPz[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparP",&_ftrueInAccpparP,"_fTrueInAccpparP[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparStartx",&_ftrueInAccpparStartx,"_ftrueInAccpparStartx[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparStarty",&_ftrueInAccpparStarty,"_ftrueInAccpparStarty[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparStartz",&_ftrueInAccpparStartz,"_ftrueInAccpparStartz[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparEndx",&_ftrueInAccpparEndx,"_ftrueInAccpparEndx[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparEndy",&_ftrueInAccpparEndy,"_ftrueInAccpparEndy[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparEndz",&_ftrueInAccpparEndz,"_ftrueInAccpparEndz[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparSPcorrStartx",&_ftrueInAccpparSPcorrStartx,"_ftrueInAccpparSPcorrStartx[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparSPcorrStarty",&_ftrueInAccpparSPcorrStarty,"_ftrueInAccpparSPcorrStarty[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparSPcorrStartz",&_ftrueInAccpparSPcorrStartz,"_ftrueInAccpparSPcorrStartz[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparSPcorrEndx",&_ftrueInAccpparSPcorrEndx,"_ftrueInAccpparSPcorrEndx[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparSPcorrEndy",&_ftrueInAccpparSPcorrEndy,"_ftrueInAccpparSPcorrEndy[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftrueInAccpparSPcorrEndz",&_ftrueInAccpparSPcorrEndz,"_ftrueInAccpparSPcorrEndz[_fNtrueInAccpmult]/F");
    fDataTreeCCSel->Branch("_fNtrueInAccpmultmu",&_fNtrueInAccpmultmu,"_fNtrueInAccpmultmu/I");
    fDataTreeCCSel->Branch("_fNtrueInAccpmultpi",&_fNtrueInAccpmultpi,"_fNtrueInAccpmultpi/I");
    fDataTreeCCSel->Branch("_fNtrueInAccpmultp",&_fNtrueInAccpmultp,"_fNtrueInAccpmultp/I");
    fDataTreeCCSel->Branch("_fNtrueInAccpmultk",&_fNtrueInAccpmultk,"_fNtrueInAccpmultk/I");
    fDataTreeCCSel->Branch("_fNtruelongInAccpmult",&_fNtruelongInAccpmult,"_fNtruelongInAccpmult/I");
    fDataTreeCCSel->Branch("_ftruelongInAccppargeantID",&_ftruelongInAccppargeantID,"_ftruelongInAccppargeantID[_fNtruelongInAccpmult]/I");
    fDataTreeCCSel->Branch("_ftruelongInAccpparpdg",&_ftruelongInAccpparpdg,"_ftruelongInAccpparpdg[_fNtruelongInAccpmult]/I");
    fDataTreeCCSel->Branch("_ftruelongInAccpparStatusCode",&_ftruelongInAccpparStatusCode,"_ftruelongInAccpparStatusCode[_fNtruelongInAccpmult]/I");
    fDataTreeCCSel->Branch("_ftruelongInAccpparTheta",&_ftruelongInAccpparTheta,"_ftruelongInAccpparTheta[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparCosTheta",&_ftruelongInAccpparCosTheta,"_ftruelongInAccpparCosTheta[_fNtruelongInAccpmult]/F");    
    fDataTreeCCSel->Branch("_ftruelongInAccpparSinTheta",&_ftruelongInAccpparSinTheta,"_ftruelongInAccpparSinTheta[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparPhi",&_ftruelongInAccpparPhi,"_ftruelongInAccpparPhi[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparCosPhi",&_ftruelongInAccpparCosPhi,"_ftruelongInAccpparCosPhi[_fNtruelongInAccpmult]/F");       
    fDataTreeCCSel->Branch("_ftruelongInAccpparSinPhi",&_ftruelongInAccpparSinPhi,"_ftruelongInAccpparSinPhi[_fNtruelongInAccpmult]/F");    
    fDataTreeCCSel->Branch("_ftruelongInAccpparE",&_ftruelongInAccpparE,"_ftruelongInAccpparE[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparMass",&_ftruelongInAccpparMass,"_ftruelongInAccpparMass[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparKE",&_ftruelongInAccpparKE,"_ftruelongInAccpparKE[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparEndE",&_ftruelongInAccpparEndE,"_ftruelongInAccpparEndE[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparPx",&_ftruelongInAccpparPx,"_ftruelongInAccpparPx[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparPy",&_ftruelongInAccpparPy,"_ftruelongInAccpparPy[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparPz",&_ftruelongInAccpparPz,"_ftruelongInAccpparPz[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparP",&_ftruelongInAccpparP,"_ftruelongInAccpparP[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparStartx",&_ftruelongInAccpparStartx,"_ftruelongInAccpparStartx[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparStarty",&_ftruelongInAccpparStarty,"_ftruelongInAccpparStarty[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparStartz",&_ftruelongInAccpparStartz,"_ftruelongInAccpparStartz[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparEndx",&_ftruelongInAccpparEndx,"_ftruelongInAccpparEndx[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparEndy",&_ftruelongInAccpparEndy,"_ftruelongInAccpparEndy[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparEndz",&_ftruelongInAccpparEndz,"_ftruelongInAccpparEndz[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparSPcorrStartx",&_ftruelongInAccpparSPcorrStartx,"_ftruelongInAccpparSPcorrStartx[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparSPcorrStarty",&_ftruelongInAccpparSPcorrStarty,"_ftruelongInAccpparSPcorrStarty[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparSPcorrStartz",&_ftruelongInAccpparSPcorrStartz,"_ftruelongInAccpparSPcorrStartz[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparSPcorrEndx",&_ftruelongInAccpparSPcorrEndx,"_ftruelongInAccpparSPcorrEndx[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparSPcorrEndy",&_ftruelongInAccpparSPcorrEndy,"_ftruelongInAccpparSPcorrEndy[_fNtruelongInAccpmult]/F");
    fDataTreeCCSel->Branch("_ftruelongInAccpparSPcorrEndz",&_ftruelongInAccpparSPcorrEndz,"_ftruelongInAccpparSPcorrEndz[_fNtruelongInAccpmult]/F");
     
    fDataTreePassBadReg->Branch("_frun",&_frun,"_frun/I");
    fDataTreePassBadReg->Branch("_fsubrun",&_fsubrun,"_fsubrun/I");
    fDataTreePassBadReg->Branch("_fevent",&_fevent,"_fevent/I"); 
    fDataTreePassBadReg->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreePassBadReg->Branch("_fNSecondaryselmult",&_fNSecondaryselmult,"_fNSecondaryselmult/I");
    fDataTreePassBadReg->Branch("_fNtrueInAccpmult",&_fNtrueInAccpmult,"_fNtrueInAccpmult/I");
    fDataTreePassBadReg->Branch("_fNtruelongInAccpmult",&_fNtruelongInAccpmult,"_fNtruelongInAccpmult/I");
    fDataTreePassBadReg->Branch("_fNtrueInAccpmultmu",&_fNtrueInAccpmultmu,"_fNtrueInAccpmultmu/I");
    fDataTreePassBadReg->Branch("_fNtrueInAccpmultpi",&_fNtrueInAccpmultpi,"_fNtrueInAccpmultpi/I");
    fDataTreePassBadReg->Branch("_fNtrueInAccpmultp",&_fNtrueInAccpmultp,"_fNtrueInAccpmultp/I");
    fDataTreePassBadReg->Branch("_fNtrueInAccpmultk",&_fNtrueInAccpmultk,"_fNtrueInAccpmultk/I");
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
    fDataTreePassCosmicReq->Branch("_fNSecondaryselmult",&_fNSecondaryselmult,"_fNSecondaryselmult/I");
    fDataTreePassCosmicReq->Branch("_fNtrueInAccpmult",&_fNtrueInAccpmult,"_fNtrueInAccpmult/I");
    fDataTreePassCosmicReq->Branch("_fNtruelongInAccpmult",&_fNtruelongInAccpmult,"_fNtruelongInAccpmult/I");
    fDataTreePassCosmicReq->Branch("_fNtrueInAccpmultmu",&_fNtrueInAccpmultmu,"_fNtrueInAccpmultmu/I");
    fDataTreePassCosmicReq->Branch("_fNtrueInAccpmultpi",&_fNtrueInAccpmultpi,"_fNtrueInAccpmultpi/I");
    fDataTreePassCosmicReq->Branch("_fNtrueInAccpmultp",&_fNtrueInAccpmultp,"_fNtrueInAccpmultp/I");
    fDataTreePassCosmicReq->Branch("_fNtrueInAccpmultk",&_fNtrueInAccpmultk,"_fNtrueInAccpmultk/I");
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
    fDataTreePassLongTrackVrtxDis->Branch("_fNSecondaryselmult",&_fNSecondaryselmult,"_fNSecondaryselmult/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fNtrueInAccpmult",&_fNtrueInAccpmult,"_fNtrueInAccpmult/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fNtruelongInAccpmult",&_fNtruelongInAccpmult,"_fNtruelongInAccpmult/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fNtrueInAccpmultmu",&_fNtrueInAccpmultmu,"_fNtrueInAccpmultmu/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fNtrueInAccpmultpi",&_fNtrueInAccpmultpi,"_fNtrueInAccpmultpi/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fNtrueInAccpmultp",&_fNtrueInAccpmultp,"_fNtrueInAccpmultp/I");
    fDataTreePassLongTrackVrtxDis->Branch("_fNtrueInAccpmultk",&_fNtrueInAccpmultk,"_fNtrueInAccpmultk/I");
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
    fDataTreePassLongcolhits->Branch("_fNrecomult",&_fNrecomult,"_fNrecomult/I");
    fDataTreePassLongcolhits->Branch("_fNSecondaryselmult",&_fNSecondaryselmult,"_fNSecondaryselmult/I");
    fDataTreePassLongcolhits->Branch("_fNtrueInAccpmult",&_fNtrueInAccpmult,"_fNtrueInAccpmult/I");
    fDataTreePassLongcolhits->Branch("_fNtruelongInAccpmult",&_fNtruelongInAccpmult,"_fNtruelongInAccpmult/I");
    fDataTreePassLongcolhits->Branch("_fNtrueInAccpmultmu",&_fNtrueInAccpmultmu,"_fNtrueInAccpmultmu/I");
    fDataTreePassLongcolhits->Branch("_fNtrueInAccpmultpi",&_fNtrueInAccpmultpi,"_fNtrueInAccpmultpi/I");
    fDataTreePassLongcolhits->Branch("_fNtrueInAccpmultp",&_fNtrueInAccpmultp,"_fNtrueInAccpmultp/I");
    fDataTreePassLongcolhits->Branch("_fNtrueInAccpmultk",&_fNtrueInAccpmultk,"_fNtrueInAccpmultk/I");
    fDataTreePassLongcolhits->Branch("_fvrtxx",&_fvrtxx,"_fvrtxx/F");
    fDataTreePassLongcolhits->Branch("_fvrtxy",&_fvrtxy,"_fvrtxy/F");
    fDataTreePassLongcolhits->Branch("_fvrtxz",&_fvrtxz,"_fvrtxz/F");
    fDataTreePassLongcolhits->Branch("_flongTrackVrtxDis",&_flongTrackVrtxDis,"_flongTrackVrtxDis/F");
    fDataTreePassLongcolhits->Branch("_fnlongcolhits",&_fnlongcolhits,"_fnlongcolhits/I");
    fDataTreePassLongcolhits->Branch("_flongtrackID",&_flongtrackID,"_flongtrackID/I"); 
    fDataTreePassLongcolhits->Branch("_flongtrackflipped",&_flongtrackflipped,"_flongtrackflipped/I");
    fDataTreePassLongcolhits->Branch("_flongtrackstartx",&_flongtrackstartx,"_flongtrackstartx/F");
    fDataTreePassLongcolhits->Branch("_flongtrackstarty",&_flongtrackstarty,"_flongtrackstarty/F");  
    fDataTreePassLongcolhits->Branch("_flongtrackstartz",&_flongtrackstartz,"_flongtrackstartz/F");
    fDataTreePassLongcolhits->Branch("_flongtrackendx",&_flongtrackendx,"_flongtrackendx/F");
    fDataTreePassLongcolhits->Branch("_flongtrackendy",&_flongtrackendy,"_flongtrackendy/F");
    fDataTreePassLongcolhits->Branch("_flongtrackendz",&_flongtrackendz,"_flongtrackendz/F");
    
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
    fDataTreeFinalSel->Branch("_fTrueccnc",&_fTrueccnc,"_fTrueccnc/I");
    fDataTreeFinalSel->Branch("_fTruemode",&_fTruemode,"_fTruemode/I");
    fDataTreeFinalSel->Branch("_fTrueinttype",&_fTrueinttype,"_fTrueinttype/I");
    fDataTreeFinalSel->Branch("_fTruenupdg",&_fTruenupdg,"_fTruenupdg/I");
    fDataTreeFinalSel->Branch("_fTrueenu",&_fTrueenu,"_fTrueenu/F");
    fDataTreeFinalSel->Branch("_fTrueq2truth",&_fTrueq2truth,"_fTrueq2truth/F");
    fDataTreeFinalSel->Branch("_fTruenuvrtxx",&_fTruenuvrtxx,"_fTruenuvrtxx/F");
    fDataTreeFinalSel->Branch("_fTruenuvrtxy",&_fTruenuvrtxy,"_fTruenuvrtxy/F");
    fDataTreeFinalSel->Branch("_fTruenuvrtxz",&_fTruenuvrtxz,"_fTruenuvrtxz/F");
    fDataTreeFinalSel->Branch("_fTrueSPcorrnuvtxx",&_fTrueSPcorrnuvtxx,"_fTrueSPcorrnuvtxx/F");
    fDataTreeFinalSel->Branch("_fTrueSPcorrnuvtxy",&_fTrueSPcorrnuvtxy,"_fTrueSPcorrnuvtxy/F");
    fDataTreeFinalSel->Branch("_fTrueSPcorrnuvtxz",&_fTrueSPcorrnuvtxz,"_fTrueSPcorrnuvtxz/F");
    fDataTreeFinalSel->Branch("_fmctrue_origin",&_fmctrue_origin,"_fmctrue_origin/I");
    fDataTreeFinalSel->Branch("_ftrueVrtxOutFV",&_ftrueVrtxOutFV,"_ftrueVrtxOutFV/I");
    fDataTreeFinalSel->Branch("_fbadReg_removed",&_fbadReg_removed,"_fbadReg_removed/I");
    fDataTreeFinalSel->Branch("_fstarty_cut",&_fstarty_cut,"_fstarty_cut/I");
    fDataTreeFinalSel->Branch("_fbrokenTracks",&_fbrokenTracks,"_fbrokenTracks/I");
    fDataTreeFinalSel->Branch("_fNGTrueallmult",&_fNGTrueallmult,"_fNGTrueallmult/I");   
    fDataTreeFinalSel->Branch("_fGtrueparID",&_fGtrueparID,"_fGtrueparID[_fNGTrueallmult]/I");
    fDataTreeFinalSel->Branch("_fGtrueparpdg",&_fGtrueparpdg,"_fGtrueparpdg[_fNGTrueallmult]/I");
    fDataTreeFinalSel->Branch("_fGtrueparStatusCode",&_fGtrueparStatusCode,"_fGtrueparStatusCode[_fNGTrueallmult]/I");
    fDataTreeFinalSel->Branch("_fGtrueparTheta",&_fGtrueparTheta,"_fGtrueparTheta[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparCosTheta",&_fGtrueparCosTheta,"_fGtrueparCosTheta[_fNGTrueallmult]/F");	
    fDataTreeFinalSel->Branch("_fGtrueparSinTheta",&_fGtrueparSinTheta,"_fGtrueparSinTheta[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparPhi",&_fGtrueparPhi,"_fGtrueparPhi[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparCosPhi",&_fGtrueparCosPhi,"_fGtrueparCosPhi[_fNGTrueallmult]/F");    
    fDataTreeFinalSel->Branch("_fGtrueparSinPhi",&_fGtrueparSinPhi,"_fGtrueparSinPhi[_fNGTrueallmult]/F");    
    fDataTreeFinalSel->Branch("_fGtrueparE",&_fGtrueparE,"_fGtrueparE[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparMass",&_fGtrueparMass,"_fGtrueparMass[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparKE",&_fGtrueparKE,"_fGtrueparKE[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparEndE",&_fGtrueparEndE,"_fGtrueparEndE[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparPx",&_fGtrueparPx,"_fGtrueparPx[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparPy",&_fGtrueparPy,"_fGtrueparPy[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparPz",&_fGtrueparPz,"_fGtrueparPz[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparP",&_fGtrueparP,"_fGtrueparP[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparStartx",&_fGtrueparStartx,"_fGtrueparStartx[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparStarty",&_fGtrueparStarty,"_fGtrueparStarty[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparStartz",&_fGtrueparStartz,"_fGtrueparStartz[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparEndx",&_fGtrueparEndx,"_fGtrueparEndx[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparEndy",&_fGtrueparEndy,"_fGtrueparEndy[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparEndz",&_fGtrueparEndz,"_fGtrueparEndz[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparSPcorrStartx",&_fGtrueparSPcorrStartx,"_fGtrueparSPcorrStartx[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparSPcorrStarty",&_fGtrueparSPcorrStarty,"_fGtrueparSPcorrStarty[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparSPcorrStartz",&_fGtrueparSPcorrStartz,"_fGtrueparSPcorrStartz[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparSPcorrEndx",&_fGtrueparSPcorrEndx,"_fGtrueparSPcorrEndx[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparSPcorrEndy",&_fGtrueparSPcorrEndy,"_fGtrueparSPcorrEndy[_fNGTrueallmult]/F");
    fDataTreeFinalSel->Branch("_fGtrueparSPcorrEndz",&_fGtrueparSPcorrEndz,"_fGtrueparSPcorrEndz[_fNGTrueallmult]/F");
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
    fDataTreeFinalSel->Branch("_fNTruematchPart",&_fNTruematchPart,"_fNTruematchPart/I");   
    fDataTreeFinalSel->Branch("_ftruematchparID",&_ftruematchparID,"_ftruematchparID[_fNTruematchPart]/I");
    fDataTreeFinalSel->Branch("_ftruematchparpdg",&_ftruematchparpdg,"_ftruematchparpdg[_fNTruematchPart]/I");
    fDataTreeFinalSel->Branch("_ftruematchparStatusCode",&_ftruematchparStatusCode,"_ftruematchparStatusCode[_fNTruematchPart]/I");
    fDataTreeFinalSel->Branch("_ftruematchparTheta",&_ftruematchparTheta,"_ftruematchparTheta[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparCosTheta",&_ftruematchparCosTheta,"_ftruematchparCosTheta[_fNTruematchPart]/F");	
    fDataTreeFinalSel->Branch("_ftruematchparSinTheta",&_ftruematchparSinTheta,"_ftruematchparSinTheta[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparPhi",&_ftruematchparPhi,"_ftruematchparPhi[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparCosPhi",&_ftruematchparCosPhi,"_ftruematchparCosPhi[_fNTruematchPart]/F");    
    fDataTreeFinalSel->Branch("_ftruematchparSinPhi",&_ftruematchparSinPhi,"_ftruematchparSinPhi[_fNTruematchPart]/F");	 
    fDataTreeFinalSel->Branch("_ftruematchparE",&_ftruematchparE,"_ftruematchparE[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparMass",&_ftruematchparMass,"_ftruematchparMass[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparKE",&_ftruematchparKE,"_ftruematchparKE[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparEndE",&_ftruematchparEndE,"_ftruematchparEndE[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparPx",&_ftruematchparPx,"_ftruematchparPx[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparPy",&_ftruematchparPy,"_ftruematchparPy[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparPz",&_ftruematchparPz,"_ftruematchparPz[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparP",&_ftruematchparP,"_fTruematchparP[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparStartx",&_ftruematchparStartx,"_ftruematchparStartx[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparStarty",&_ftruematchparStarty,"_ftruematchparStarty[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparStartz",&_ftruematchparStartz,"_ftruematchparStartz[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparEndx",&_ftruematchparEndx,"_ftruematchparEndx[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparEndy",&_ftruematchparEndy,"_ftruematchparEndy[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparEndz",&_ftruematchparEndz,"_ftruematchparEndz[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparSPcorrStartx",&_ftruematchparSPcorrStartx,"_ftruematchparSPcorrStartx[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparSPcorrStarty",&_ftruematchparSPcorrStarty,"_ftruematchparSPcorrStarty[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparSPcorrStartz",&_ftruematchparSPcorrStartz,"_ftruematchparSPcorrStartz[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparSPcorrEndx",&_ftruematchparSPcorrEndx,"_ftruematchparSPcorrEndx[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparSPcorrEndy",&_ftruematchparSPcorrEndy,"_ftruematchparSPcorrEndy[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_ftruematchparSPcorrEndz",&_ftruematchparSPcorrEndz,"_ftruematchparSPcorrEndz[_fNTruematchPart]/F");
    fDataTreeFinalSel->Branch("_fNSecondaryselmult",&_fNSecondaryselmult,"_fNSecondaryselmult/I");
    fDataTreeFinalSel->Branch("_fSeconpargeantID",&_fSeconpargeantID,"_fSeconpargeantID[_fNSecondaryselmult]/I");
    fDataTreeFinalSel->Branch("_fSeconparpdg",&_fSeconparpdg,"_fSeconparpdg[_fNSecondaryselmult]/I");
    fDataTreeFinalSel->Branch("_fSeconparStatusCode",&_fSeconparStatusCode,"_fSeconparStatusCode[_fNSecondaryselmult]/I");
    fDataTreeFinalSel->Branch("_fSeconparTheta",&_fSeconparTheta,"_fSeconparTheta[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparCosTheta",&_fSeconparCosTheta,"_fSeconparCosTheta[_fNSecondaryselmult]/F");    
    fDataTreeFinalSel->Branch("_fSeconparSinTheta",&_fSeconparSinTheta,"_fSeconparSinTheta[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparPhi",&_fSeconparPhi,"_fSeconparPhi[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparCosPhi",&_fSeconparCosPhi,"_fSeconparCosPhi[_fNSecondaryselmult]/F");	 
    fDataTreeFinalSel->Branch("_fSeconparSinPhi",&_fSeconparSinPhi,"_fSeconparSinPhi[_fNSecondaryselmult]/F");    
    fDataTreeFinalSel->Branch("_fSeconparE",&_fSeconparE,"_fSeconparE[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparMass",&_fSeconparMass,"_fSeconparMass[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparKE",&_fSeconparKE,"_fSeconparKE[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparEndE",&_fSeconparEndE,"_fSeconparEndE[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparPx",&_fSeconparPx,"_fSeconparPx[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparPy",&_fSeconparPy,"_fSeconparPy[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparPz",&_fSeconparPz,"_fSeconparPz[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparP",&_fSeconparP,"_fSeconparP[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparStartx",&_fSeconparStartx,"_fSeconparStartx[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparStarty",&_fSeconparStarty,"_fSeconparStarty[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparStartz",&_fSeconparStartz,"_fSeconparStartz[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparEndx",&_fSeconparEndx,"_fSeconparEndx[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparEndy",&_fSeconparEndy,"_fSeconparEndy[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparEndz",&_fSeconparEndz,"_fSeconparEndz[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparSPcorrStartx",&_fSeconparSPcorrStartx,"_fSeconparSPcorrStartx[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparSPcorrStarty",&_fSeconparSPcorrStarty,"_fSeconparSPcorrStarty[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparSPcorrStartz",&_fSeconparSPcorrStartz,"_fSeconparSPcorrStartz[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparSPcorrEndx",&_fSeconparSPcorrEndx,"_fSeconparSPcorrEndx[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparSPcorrEndy",&_fSeconparSPcorrEndy,"_fSeconparSPcorrEndy[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fSeconparSPcorrEndz",&_fSeconparSPcorrEndz,"_fSeconparSPcorrEndz[_fNSecondaryselmult]/F");
    fDataTreeFinalSel->Branch("_fNTrueallPart",&_fNTrueallPart,"_fNTrueallPart/I");   
/*  fDataTreeFinalSel->Branch("_ftrueparID",&_ftrueparID,"_ftrueparID[_fNTrueallPart]/I");
    fDataTreeFinalSel->Branch("_ftrueparpdg",&_ftrueparpdg,"_ftrueparpdg[_fNTrueallPart]/I");
    fDataTreeFinalSel->Branch("_ftrueparStatusCode",&_ftrueparStatusCode,"_ftrueparStatusCode[_fNTrueallPart]/I");
    fDataTreeFinalSel->Branch("_ftrueparTheta",&_ftrueparTheta,"_ftrueparTheta[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparCosTheta",&_ftrueparCosTheta,"_ftrueparCosTheta[_fNTrueallPart]/F");	
    fDataTreeFinalSel->Branch("_ftrueparSinTheta",&_ftrueparSinTheta,"_ftrueparSinTheta[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparPhi",&_ftrueparPhi,"_ftrueparPhi[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparCosPhi",&_ftrueparCosPhi,"_ftrueparCosPhi[_fNTrueallPart]/F");    
    fDataTreeFinalSel->Branch("_ftrueparSinPhi",&_ftrueparSinPhi,"_ftrueparSinPhi[_fNTrueallPart]/F");    
    fDataTreeFinalSel->Branch("_ftrueparE",&_ftrueparE,"_ftrueparE[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparMass",&_ftrueparMass,"_ftrueparMass[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparKE",&_ftrueparKE,"_ftrueparKE[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparEndE",&_ftrueparEndE,"_ftrueparEndE[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparPx",&_ftrueparPx,"_ftrueparPx[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparPy",&_ftrueparPy,"_ftrueparPy[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparPz",&_ftrueparPz,"_ftrueparPz[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparP",&_ftrueparP,"_fTrueparP[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparStartx",&_ftrueparStartx,"_ftrueparStartx[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparStarty",&_ftrueparStarty,"_ftrueparStarty[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparStartz",&_ftrueparStartz,"_ftrueparStartz[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparEndx",&_ftrueparEndx,"_ftrueparEndx[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparEndy",&_ftrueparEndy,"_ftrueparEndy[_fNTrueallPart]/F");
    fDataTreeFinalSel->Branch("_ftrueparEndz",&_ftrueparEndz,"_ftrueparEndz[_fNTrueallPart]/F");
*/    fDataTreeFinalSel->Branch("_fNtrueInAccpmult",&_fNtrueInAccpmult,"_fNtrueInAccpmult/I");
    fDataTreeFinalSel->Branch("_ftrueInAccppargeantID",&_ftrueInAccppargeantID,"_ftrueInAccppargeantID[_fNtrueInAccpmult]/I");
    fDataTreeFinalSel->Branch("_ftrueInAccpparpdg",&_ftrueInAccpparpdg,"_ftrueInAccpparpdg[_fNtrueInAccpmult]/I");
    fDataTreeFinalSel->Branch("_ftrueInAccpparStatusCode",&_ftrueInAccpparStatusCode,"_ftrueInAccpparStatusCode[_fNtrueInAccpmult]/I");
    fDataTreeFinalSel->Branch("_ftrueInAccpparTheta",&_ftrueInAccpparTheta,"_ftrueInAccpparTheta[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparCosTheta",&_ftrueInAccpparCosTheta,"_ftrueInAccpparCosTheta[_fNtrueInAccpmult]/F");	
    fDataTreeFinalSel->Branch("_ftrueInAccpparSinTheta",&_ftrueInAccpparSinTheta,"_ftrueInAccpparSinTheta[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparPhi",&_ftrueInAccpparPhi,"_ftrueInAccpparPhi[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparCosPhi",&_ftrueInAccpparCosPhi,"_ftrueInAccpparCosPhi[_fNtrueInAccpmult]/F");    
    fDataTreeFinalSel->Branch("_ftrueInAccpparSinPhi",&_ftrueInAccpparSinPhi,"_ftrueInAccpparSinPhi[_fNtrueInAccpmult]/F");    
    fDataTreeFinalSel->Branch("_ftrueInAccpparE",&_ftrueInAccpparE,"_ftrueInAccpparE[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparMass",&_ftrueInAccpparMass,"_ftrueInAccpparMass[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparKE",&_ftrueInAccpparKE,"_ftrueInAccpparKE[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparEndE",&_ftrueInAccpparEndE,"_ftrueInAccpparEndE[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparPx",&_ftrueInAccpparPx,"_ftrueInAccpparPx[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparPy",&_ftrueInAccpparPy,"_ftrueInAccpparPy[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparPz",&_ftrueInAccpparPz,"_ftrueInAccpparPz[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparP",&_ftrueInAccpparP,"_fTrueInAccpparP[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparStartx",&_ftrueInAccpparStartx,"_ftrueInAccpparStartx[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparStarty",&_ftrueInAccpparStarty,"_ftrueInAccpparStarty[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparStartz",&_ftrueInAccpparStartz,"_ftrueInAccpparStartz[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparEndx",&_ftrueInAccpparEndx,"_ftrueInAccpparEndx[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparEndy",&_ftrueInAccpparEndy,"_ftrueInAccpparEndy[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparEndz",&_ftrueInAccpparEndz,"_ftrueInAccpparEndz[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparSPcorrStartx",&_ftrueInAccpparSPcorrStartx,"_ftrueInAccpparSPcorrStartx[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparSPcorrStarty",&_ftrueInAccpparSPcorrStarty,"_ftrueInAccpparSPcorrStarty[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparSPcorrStartz",&_ftrueInAccpparSPcorrStartz,"_ftrueInAccpparSPcorrStartz[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparSPcorrEndx",&_ftrueInAccpparSPcorrEndx,"_ftrueInAccpparSPcorrEndx[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparSPcorrEndy",&_ftrueInAccpparSPcorrEndy,"_ftrueInAccpparSPcorrEndy[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftrueInAccpparSPcorrEndz",&_ftrueInAccpparSPcorrEndz,"_ftrueInAccpparSPcorrEndz[_fNtrueInAccpmult]/F");
    fDataTreeFinalSel->Branch("_fNtrueInAccpmultmu",&_fNtrueInAccpmultmu,"_fNtrueInAccpmultmu/I");
    fDataTreeFinalSel->Branch("_fNtrueInAccpmultpi",&_fNtrueInAccpmultpi,"_fNtrueInAccpmultpi/I");
    fDataTreeFinalSel->Branch("_fNtrueInAccpmultp",&_fNtrueInAccpmultp,"_fNtrueInAccpmultp/I");
    fDataTreeFinalSel->Branch("_fNtrueInAccpmultk",&_fNtrueInAccpmultk,"_fNtrueInAccpmultk/I");
    fDataTreeFinalSel->Branch("_fNtruelongInAccpmult",&_fNtruelongInAccpmult,"_fNtruelongInAccpmult/I");
    fDataTreeFinalSel->Branch("_ftruelongInAccppargeantID",&_ftruelongInAccppargeantID,"_ftruelongInAccppargeantID[_fNtruelongInAccpmult]/I");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparpdg",&_ftruelongInAccpparpdg,"_ftruelongInAccpparpdg[_fNtruelongInAccpmult]/I");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparStatusCode",&_ftruelongInAccpparStatusCode,"_ftruelongInAccpparStatusCode[_fNtruelongInAccpmult]/I");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparTheta",&_ftruelongInAccpparTheta,"_ftruelongInAccpparTheta[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparCosTheta",&_ftruelongInAccpparCosTheta,"_ftruelongInAccpparCosTheta[_fNtruelongInAccpmult]/F");    
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSinTheta",&_ftruelongInAccpparSinTheta,"_ftruelongInAccpparSinTheta[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparPhi",&_ftruelongInAccpparPhi,"_ftruelongInAccpparPhi[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparCosPhi",&_ftruelongInAccpparCosPhi,"_ftruelongInAccpparCosPhi[_fNtruelongInAccpmult]/F");       
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSinPhi",&_ftruelongInAccpparSinPhi,"_ftruelongInAccpparSinPhi[_fNtruelongInAccpmult]/F");    
    fDataTreeFinalSel->Branch("_ftruelongInAccpparE",&_ftruelongInAccpparE,"_ftruelongInAccpparE[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparMass",&_ftruelongInAccpparMass,"_ftruelongInAccpparMass[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparKE",&_ftruelongInAccpparKE,"_ftruelongInAccpparKE[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparEndE",&_ftruelongInAccpparEndE,"_ftruelongInAccpparEndE[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparPx",&_ftruelongInAccpparPx,"_ftruelongInAccpparPx[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparPy",&_ftruelongInAccpparPy,"_ftruelongInAccpparPy[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparPz",&_ftruelongInAccpparPz,"_ftruelongInAccpparPz[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparP",&_ftruelongInAccpparP,"_ftruelongInAccpparP[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparStartx",&_ftruelongInAccpparStartx,"_ftruelongInAccpparStartx[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparStarty",&_ftruelongInAccpparStarty,"_ftruelongInAccpparStarty[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparStartz",&_ftruelongInAccpparStartz,"_ftruelongInAccpparStartz[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparEndx",&_ftruelongInAccpparEndx,"_ftruelongInAccpparEndx[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparEndy",&_ftruelongInAccpparEndy,"_ftruelongInAccpparEndy[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparEndz",&_ftruelongInAccpparEndz,"_ftruelongInAccpparEndz[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSPcorrStartx",&_ftruelongInAccpparSPcorrStartx,"_ftruelongInAccpparSPcorrStartx[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSPcorrStarty",&_ftruelongInAccpparSPcorrStarty,"_ftruelongInAccpparSPcorrStarty[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSPcorrStartz",&_ftruelongInAccpparSPcorrStartz,"_ftruelongInAccpparSPcorrStartz[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSPcorrEndx",&_ftruelongInAccpparSPcorrEndx,"_ftruelongInAccpparSPcorrEndx[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSPcorrEndy",&_ftruelongInAccpparSPcorrEndy,"_ftruelongInAccpparSPcorrEndy[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_ftruelongInAccpparSPcorrEndz",&_ftruelongInAccpparSPcorrEndz,"_ftruelongInAccpparSPcorrEndz[_fNtruelongInAccpmult]/F");
    fDataTreeFinalSel->Branch("_fPH",&_fPH,"_fPH/I");
    fDataTreeFinalSel->Branch("_fMCS",&_fMCS,"_fMCS/I");
    fDataTreeFinalSel->Branch("_fPHratio",&_fPHratio,"_fPHratio/F");
    fDataTreeFinalSel->Branch("_fMCSratio",&_fMCSratio,"_fMCSratio/F");
    fDataTreeFinalSel->Branch("_fMCSdiff",&_fMCSdiff,"_fMCSdiff/F");
             
   }
   
   //-----------------------------------------------------------------------
   void CTMMCAna::endJob()
  {

  }
  
   //-----------------------------------------------------------------------
   void CTMMCAna::beginRun(const art::Run& run)
  {
    
  }
    //-----------------------------------------------------------------------
   void CTMMCAna::reconfigure(fhicl::ParameterSet const& pset)
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
    fBeamMin                 = pset.get<double>      ("BeamMin", 3.2);
    fBeamMax                 = pset.get<double>      ("BeamMax", 4.8);
    fPEThresh                = pset.get<double>      ("PEThresh", 50.);
    fMinTrk2VtxDist          = pset.get<double>      ("MinTrk2VtxDist", 5.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen", 75.);

  // Get the tool for MC Truth matching
  const fhicl::ParameterSet& truthParams = pset.get<fhicl::ParameterSet>("MCTruthMatching");
  
    if (truthParams.get<std::string>("tool_type") == "AssociationsTruth")
    {
        fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(truthParams));
    }
    else
    {
        fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::BackTrackerTruth(truthParams));
    }

   }
 
 //========================================================================	
bool CTMMCAna::inFV(double x, double y, double z) const
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
double CTMMCAna::FlashTrackDist(double flash, double start, double end) const
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
 //=========================================================================
  void CTMMCAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet)
  {		
  //art::ServiceHandle<cheat::BackTracker> bt;	
  std::map<int,double> trkID_E;	
  for(size_t j = 0; j < track_hits.size(); ++j)
  {	
    art::Ptr<recob::Hit> hit = track_hits[j];
    //const auto& hit = *track_hits[j];
    //std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
    for(size_t k = 0; k < TrackIDs.size(); k++)
    {	
      trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
    }
  }

  double E_em =0.0;	
  double max_E = -999.0;	
  double total_E = 0.0;	
  int TrackID = -999;	
  double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla	
  //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 	
  //!since we are looking for muons/pions/protons this should be enough 	
  if( !trkID_E.size() ) 
  {
    MCparticle = 0;	
    return; //Ghost track???	
  }	
  for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii)
  {	
    total_E += ii->second;	
    if((ii->second)>max_E)
    {	
      partial_E = ii->second;
      max_E = ii->second;
      TrackID = ii->first;
      if( TrackID < 0 ) E_em += ii->second;
    }	
  }
   	
  //MCparticle = bt->TrackIDToParticle(TrackID);		
  MCparticle = fMCTruthMatching->TrackIDToParticle(TrackID);		
  //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of 	
  //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle	
  //we don't want to track gammas or any other EM activity 	
  if( TrackID < 0 ) return;		
  //Efrac = (partial_E+E_em)/total_E;	
  Efrac = (partial_E)/total_E;		
  //completeness	
  double totenergy =0;	
  for(size_t k = 0; k < all_hits.size(); ++k)
  {	
    art::Ptr<recob::Hit> hit = all_hits[k];	
    //std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);	
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);	
    for(size_t l = 0; l < TrackIDs.size(); ++l)
    {	
      if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;	
    }	
  } 	
  Ecomplet = partial_E/totenergy;
}

   //----------------------------------------------------------------------
   float CTMMCAna::T(int pdg, float range)
   {
     float ans;
     if(TMath::Abs(pdg)==13)       {ans=(11.518*(TMath::Power(range,0.555))+0.252*(TMath::Power(range,1.289)))/1000;}//converting to GeV 
     else if(TMath::Abs(pdg)==211) {ans=(13.069*(TMath::Power(range,0.555))+0.233*(TMath::Power(range,1.289)))/1000;} 
     else if(pdg==2212)            {ans=(30.573*(TMath::Power(range,0.555))+0.133*(TMath::Power(range,1.289)))/1000;} 
     else if(TMath::Abs(pdg)==321) {ans=(22.966*(TMath::Power(range,0.555))+0.161*(TMath::Power(range,1.289)))/1000;}	 
     else                          {ans=-999;}
    // std::cout<<"KE limiting values: "<<ans<<"\n";
     return ans;
   }	
   //-----------------------------------------------------------------------

   void CTMMCAna::analyze(art::Event const& evt) 
   {
    _frun = evt.run();
    _fsubrun = evt.subRun();
    _fevent = evt.id().event();
    
    fDataTreeTotal->Fill();

    fMCTruthMatching->Rebuild(evt);
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
      ////////////////----------Back Tracking--------//////////////////////////////
      
      art::FindManyP<recob::Hit> track_hits(trackCol, evt, fTrackModuleLabel);
      std::vector<art::Ptr<recob::Hit>> longtrackhits = track_hits.at(longtrack->ID());
      double tmpEfrac = 0;	
      double tmpEcomplet =0;	
      const simb::MCParticle *particle;	
      truthMatcher( hitlist, longtrackhits, particle, tmpEfrac, tmpEcomplet );	

      //if (!particle) continue;	
      
      //art::ServiceHandle<cheat::BackTracker> bt;
      //const art::Ptr<simb::MCTruth> MCtruth = bt->ParticleToMCTruth(particle);

      const art::Ptr<simb::MCTruth> MCtruth = fMCTruthMatching->ParticleToMCTruth(particle);
      std::string pri("primary");

      _fmctrue_origin=MCtruth->Origin();
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 

      if (MCtruth->NeutrinoSet() && MCtruth->Origin()==1)
      {
        _fTrueccnc=MCtruth->GetNeutrino().CCNC();
	_fTruemode=MCtruth->GetNeutrino().Mode();
	_fTrueinttype=MCtruth->GetNeutrino().InteractionType();
	_fTruenupdg=MCtruth->GetNeutrino().Nu().PdgCode();
	_fTrueenu=MCtruth->GetNeutrino().Nu().E();
	_fTrueq2truth=MCtruth->GetNeutrino().QSqr();
	_fTruenuvrtxx=MCtruth->GetNeutrino().Nu().Vx();
	_fTruenuvrtxy=MCtruth->GetNeutrino().Nu().Vy();
	_fTruenuvrtxz=MCtruth->GetNeutrino().Nu().Vz();
	_fTrueSPcorrnuvtxx = MCtruth->GetNeutrino().Nu().Vx() - SCE->GetPosOffsets(MCtruth->GetNeutrino().Nu().Vx(),MCtruth->GetNeutrino().Nu().Vy(),MCtruth->GetNeutrino().Nu().Vz())[0];
	_fTrueSPcorrnuvtxy = MCtruth->GetNeutrino().Nu().Vy() + SCE->GetPosOffsets(MCtruth->GetNeutrino().Nu().Vx(),MCtruth->GetNeutrino().Nu().Vy(),MCtruth->GetNeutrino().Nu().Vz())[1];
	_fTrueSPcorrnuvtxz = MCtruth->GetNeutrino().Nu().Vz() + SCE->GetPosOffsets(MCtruth->GetNeutrino().Nu().Vx(),MCtruth->GetNeutrino().Nu().Vy(),MCtruth->GetNeutrino().Nu().Vz())[2];
      }  
      else 
      {
        _fTrueccnc=-999;
	_fTruemode=-999;
	_fTrueinttype=-999;
	_fTruenupdg=-999;
	_fTrueenu=-999;
	_fTrueq2truth=-999;
	_fTruenuvrtxx=-999;
	_fTruenuvrtxy=-999;
	_fTruenuvrtxz=-999;
        _fTrueSPcorrnuvtxx=-999;
        _fTrueSPcorrnuvtxx=-999;
        _fTrueSPcorrnuvtxx=-999;
      }  
       
      /////------Checking whether true true vertex is inside its limits or not--------/////////
      int trueVrtxOutFV = 0;
      if(_fTruenuvrtxx < vXnegbound || _fTruenuvrtxx > vXposbound
       ||_fTruenuvrtxy < vYnegbound || _fTruenuvrtxy > vYposbound
       ||_fTruenuvrtxz < vZnegbound || _fTruenuvrtxz > vZposbound)
       {
         trueVrtxOutFV = 1;
       }
       _ftrueVrtxOutFV = trueVrtxOutFV;       
                  
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

      ////////////////Loop over all GENIE particles only////////////////////////////	
      vector <int> genieparID;
      vector <int> genieparpdg;
      vector <int> genieparStatusCode;
      vector <float> genieparCosTheta;
      vector <float> genieparStartz;
      _fNGTrueallmult= 0 ; // in case we have cosmic event
      if (MCtruth->NeutrinoSet())
      {  	  
      _fNGTrueallmult = MCtruth->NParticles();
      size_t genie_particle=0;
        for(Int_t iPartc = 0; iPartc < MCtruth->NParticles(); ++iPartc)
        {       
           const simb::MCParticle& part(MCtruth->GetParticle(iPartc));        
	    _fGtrueparID[genie_particle]=part.TrackId();
	    _fGtrueparpdg[genie_particle]=part.PdgCode();
	    _fGtrueparStatusCode[genie_particle]=part.StatusCode();
	    _fGtrueparTheta[genie_particle]=part.Momentum().Theta();
	    _fGtrueparCosTheta[genie_particle]=TMath::Cos(part.Momentum().Theta());	
	    _fGtrueparSinTheta[genie_particle]=TMath::Sin(part.Momentum().Theta());	
	    _fGtrueparE[genie_particle]=part.E();
	    _fGtrueparMass[genie_particle]=part.Mass();
	    _fGtrueparKE[genie_particle]=part.E()-part.Mass();
	    _fGtrueparEndE[genie_particle]=part.EndE();
            _fGtrueparPx[genie_particle]=part.Px();
	    _fGtrueparPy[genie_particle]=part.Py();
	    _fGtrueparPz[genie_particle]=part.Pz();
	    _fGtrueparP[genie_particle]=part.Momentum().Vect().Mag();
	    _fGtrueparStartx[genie_particle]=part.Vx();
	    _fGtrueparStarty[genie_particle]=part.Vy();
	    _fGtrueparStartz[genie_particle]=part.Vz();
	    _fGtrueparEndx[genie_particle]=part.EndPosition()[0];
	    _fGtrueparEndy[genie_particle]=part.EndPosition()[1];
	    _fGtrueparEndz[genie_particle]=part.EndPosition()[2];
	    _fGtrueparSPcorrStartx[genie_particle]=part.Vx() - SCE->GetPosOffsets(part.Vx(),part.Vy(),part.Vz())[0]; 
	    _fGtrueparSPcorrStarty[genie_particle]=part.Vy() + SCE->GetPosOffsets(part.Vx(),part.Vy(),part.Vz())[1];
	    _fGtrueparSPcorrStartz[genie_particle]=part.Vz() + SCE->GetPosOffsets(part.Vx(),part.Vy(),part.Vz())[2]; 
	    _fGtrueparSPcorrEndx[genie_particle]=part.EndPosition()[0] - SCE->GetPosOffsets(part.EndPosition()[0],part.EndPosition()[1],part.EndPosition()[2])[0]; 
	    _fGtrueparSPcorrEndy[genie_particle]=part.EndPosition()[1] + SCE->GetPosOffsets(part.EndPosition()[0],part.EndPosition()[1],part.EndPosition()[2])[1];
	    _fGtrueparSPcorrEndz[genie_particle]=part.EndPosition()[2] + SCE->GetPosOffsets(part.EndPosition()[0],part.EndPosition()[1],part.EndPosition()[2])[2]; 
	    _fGtrueparPhi[genie_particle]=part.Momentum().Phi();	
	    _fGtrueparCosPhi[genie_particle]=TMath::Cos(part.Momentum().Phi());
	    _fGtrueparSinPhi[genie_particle]=TMath::Sin(part.Momentum().Phi());
	    
            genieparID.push_back(_fGtrueparID[genie_particle]);
            genieparpdg.push_back(_fGtrueparpdg[genie_particle]);
            genieparStatusCode.push_back(_fGtrueparStatusCode[genie_particle]);
            genieparCosTheta.push_back(_fGtrueparCosTheta[genie_particle]);
            genieparStartz.push_back(_fGtrueparStartz[genie_particle]);	 
	    ++genie_particle;
		
	 }//particle loop
       }//MCtruth->NeutrinoSet()       	  
      //-----------------Reconstructed, truth macthing, secondary multiplicty---------------------------//
      vector<int> parID;
      int brkntrkcount=0;
      int recomult=0;
      _fNRecoallPart = TrackVector.size();

      _fNrecomult = 0;//starting value
      _fNTruematchPart = 0;//starting value
      _fNSecondaryselmult = 0;//starting value

      size_t nrecopart = 0;
      size_t nselrecopart = 0;
      size_t ntruematchpart = 0;
      size_t nsecondaryPar = 0;
          
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

	  double tmpEfracm = 0;	
          double tmpEcompletm =0;	
          const simb::MCParticle *mparticle;	

          truthMatcher( hitlist, alltrackhits, mparticle, tmpEfracm, tmpEcompletm );	
          if (!mparticle) continue;
	
	  parID.push_back(mparticle->TrackId());
	  _ftruematchparID[ntruematchpart]=mparticle->TrackId();
	  _ftruematchparpdg[ntruematchpart]=mparticle->PdgCode();
	  _ftruematchparStatusCode[ntruematchpart]=mparticle->StatusCode();
	  _ftruematchparTheta[ntruematchpart]=mparticle->Momentum().Theta();
	  _ftruematchparCosTheta[ntruematchpart]=TMath::Cos(mparticle->Momentum().Theta());	
	  _ftruematchparSinTheta[ntruematchpart]=TMath::Sin(mparticle->Momentum().Theta());	
	  _ftruematchparE[ntruematchpart]=mparticle->E();
	  _ftruematchparMass[ntruematchpart]=mparticle->Mass();
	  _ftruematchparKE[ntruematchpart]=(mparticle->E())-(mparticle->Mass());
	  _ftruematchparEndE[ntruematchpart]=mparticle->EndE();
          _ftruematchparPx[ntruematchpart]=mparticle->Px();
	  _ftruematchparPy[ntruematchpart]=mparticle->Py();
	  _ftruematchparPz[ntruematchpart]=mparticle->Pz();
	  _ftruematchparP[ntruematchpart]=mparticle->Momentum().Vect().Mag();
	  _ftruematchparStartx[ntruematchpart]=mparticle->Vx();
	  _ftruematchparStarty[ntruematchpart]=mparticle->Vy();
	  _ftruematchparStartz[ntruematchpart]=mparticle->Vz();
	  _ftruematchparEndx[ntruematchpart]=mparticle->EndPosition()[0];
	  _ftruematchparEndy[ntruematchpart]=mparticle->EndPosition()[1];
	  _ftruematchparEndz[ntruematchpart]=mparticle->EndPosition()[2];
	  _ftruematchparSPcorrStartx[ntruematchpart]=mparticle->Vx() - SCE->GetPosOffsets(mparticle->Vx(),mparticle->Vy(),mparticle->Vz())[0]; 
	  _ftruematchparSPcorrStarty[ntruematchpart]=mparticle->Vy() + SCE->GetPosOffsets(mparticle->Vx(),mparticle->Vy(),mparticle->Vz())[1];
	  _ftruematchparSPcorrStartz[ntruematchpart]=mparticle->Vz() + SCE->GetPosOffsets(mparticle->Vx(),mparticle->Vy(),mparticle->Vz())[2]; 
	  _ftruematchparSPcorrEndx[ntruematchpart]=mparticle->EndPosition()[0] - SCE->GetPosOffsets(mparticle->EndPosition()[0],mparticle->EndPosition()[1],mparticle->EndPosition()[2])[0]; 
	  _ftruematchparSPcorrEndy[ntruematchpart]=mparticle->EndPosition()[1] + SCE->GetPosOffsets(mparticle->EndPosition()[0],mparticle->EndPosition()[1],mparticle->EndPosition()[2])[1];
	  _ftruematchparSPcorrEndz[ntruematchpart]=mparticle->EndPosition()[2] + SCE->GetPosOffsets(mparticle->EndPosition()[0],mparticle->EndPosition()[1],mparticle->EndPosition()[2])[2]; 
	  _ftruematchparPhi[ntruematchpart]=mparticle->Momentum().Phi();
	  _ftruematchparCosPhi[ntruematchpart]=TMath::Cos(mparticle->Momentum().Phi());
	  _ftruematchparSinPhi[ntruematchpart]=TMath::Sin(mparticle->Momentum().Phi());
	 	  
	  bool selparyes = false;
	  for(size_t accpp=0; accpp < genieparpdg.size(); accpp++)
	  {
	    //std::cout<<_fmparpdg<<" "<<genieparpdg.at(accpp)<<" "<<_fmparStartz<<" "<<genieparStartz.at(accpp)<<"\n";
//	    if(genieparpdg.at(accpp)==_ftruematchparpdg[ntruematchpart] && genieparStatusCode.at(accpp)==_ftruematchparStatusCode[ntruematchpart] &&
//	    genieparCosTheta.at(accpp)==_ftruematchparCosTheta[ntruematchpart] && genieparStartz.at(accpp)==_ftruematchparStartz[ntruematchpart])

	    if(genieparID.at(accpp)==_ftruematchparID[ntruematchpart])
            {
              selparyes = true;
	      //std::cout<<"accpted!"<<"\n";
            }
	  }
	  //std::cout<<selparyes<<"\n";
	  if(selparyes==false)
	  {
	  _fSeconpargeantID[nsecondaryPar]=mparticle->TrackId();
	  _fSeconparpdg[nsecondaryPar]=mparticle->PdgCode();
	  _fSeconparStatusCode[nsecondaryPar]=mparticle->StatusCode();
	  _fSeconparTheta[nsecondaryPar]=mparticle->Momentum().Theta();
	  _fSeconparCosTheta[nsecondaryPar]=TMath::Cos(mparticle->Momentum().Theta());	
	  _fSeconparSinTheta[nsecondaryPar]=TMath::Sin(mparticle->Momentum().Theta());	
	  _fSeconparE[nsecondaryPar]=mparticle->E();
	  _fSeconparMass[nsecondaryPar]=mparticle->Mass();
	  _fSeconparKE[nsecondaryPar]=(mparticle->E())-(mparticle->Mass());
	  _fSeconparEndE[nsecondaryPar]=mparticle->EndE();
          _fSeconparPx[nsecondaryPar]=mparticle->Px();
	  _fSeconparPy[nsecondaryPar]=mparticle->Py();
	  _fSeconparPz[nsecondaryPar]=mparticle->Pz();
	  _fSeconparP[nsecondaryPar]=mparticle->Momentum().Vect().Mag();
	  _fSeconparStartx[nsecondaryPar]=mparticle->Vx();
	  _fSeconparStarty[nsecondaryPar]=mparticle->Vy();
	  _fSeconparStartz[nsecondaryPar]=mparticle->Vz();
	  _fSeconparEndx[nsecondaryPar]=mparticle->EndPosition()[0];
	  _fSeconparEndy[nsecondaryPar]=mparticle->EndPosition()[1];
	  _fSeconparEndz[nsecondaryPar]=mparticle->EndPosition()[2];
	  _fSeconparSPcorrStartx[nsecondaryPar]=mparticle->Vx() - SCE->GetPosOffsets(mparticle->Vx(),mparticle->Vy(),mparticle->Vz())[0]; 
	  _fSeconparSPcorrStarty[nsecondaryPar]=mparticle->Vy() + SCE->GetPosOffsets(mparticle->Vx(),mparticle->Vy(),mparticle->Vz())[1];
	  _fSeconparSPcorrStartz[nsecondaryPar]=mparticle->Vz() + SCE->GetPosOffsets(mparticle->Vx(),mparticle->Vy(),mparticle->Vz())[2]; 
	  _fSeconparSPcorrEndx[nsecondaryPar]=mparticle->EndPosition()[0] - SCE->GetPosOffsets(mparticle->EndPosition()[0],mparticle->EndPosition()[1],mparticle->EndPosition()[2])[0]; 
	  _fSeconparSPcorrEndy[nsecondaryPar]=mparticle->EndPosition()[1] + SCE->GetPosOffsets(mparticle->EndPosition()[0],mparticle->EndPosition()[1],mparticle->EndPosition()[2])[1];
	  _fSeconparSPcorrEndz[nsecondaryPar]=mparticle->EndPosition()[2] + SCE->GetPosOffsets(mparticle->EndPosition()[0],mparticle->EndPosition()[1],mparticle->EndPosition()[2])[2]; 
	  _fSeconparPhi[nsecondaryPar]=mparticle->Momentum().Phi();
	  _fSeconparCosPhi[nsecondaryPar]=TMath::Cos(mparticle->Momentum().Phi());
	  _fSeconparSinPhi[nsecondaryPar]=TMath::Sin(mparticle->Momentum().Phi());
	  ++nsecondaryPar;
	  
	  }
	  ++nselrecopart;
	  ++ntruematchpart;   
	  }
	}
	  ++nrecopart;
      }	
      _fNrecomult=recomult;
      _fNTruematchPart = ntruematchpart;
      _fNSecondaryselmult = nsecondaryPar;
       
       //////////broken tracks counter////////////////	
      _fbrokenTracks = 0; //starting value
      for (unsigned bc1=0; bc1<parID.size(); bc1++)
      {
        for (unsigned bc2=bc1+1; bc2<parID.size(); bc2++)
        {
	  //std::cout<<"track ids "<<parID.at(bc1)<<" "<<parID.at(bc2)<<endl;
	  if(parID.at(bc1)==parID.at(bc2)) 
	  {
	    std::cout<<"broken track loop "<<parID.at(bc1)<<" "<<parID.at(bc2)<<endl;
	    brkntrkcount++;
	  }
	}
      }	    
	  _fbrokenTracks=brkntrkcount;
	  
      
      //-------Calculating true multiplicity within our analysis acceptance-----//
      
      int TrueInAccplongcontained=0, trueInAccpmult=0, trueInAccpmultmu=0, trueInAccpmultpi=0, trueInAccpmultp=0,trueInAccpmultk=0;
      _fNTrueallPart = 0; //in case of cosmic event    
      _fNtrueInAccpmult = 0; //in case of cosmic event  
      _fNtrueInAccpmultmu = 0; //in case of cosmic event  
      _fNtrueInAccpmultpi = 0; //in case of cosmic event  
      _fNtrueInAccpmultp = 0; //in case of cosmic event  
      _fNtrueInAccpmultk = 0; //in case of cosmic event  
      _fNtruelongInAccpmult = 0; //in case of cosmic event  
  
      if (MCtruth->NeutrinoSet())
      {   
         //const sim::ParticleList& plist = bt->ParticleList();
         const sim::ParticleList& plist = fMCTruthMatching->ParticleList();
         _fNTrueallPart = plist.size();
         size_t ntruepart = 0;
	 size_t trueInAccep = 0;
	 size_t truelongInAccep = 0;
	 sim::ParticleList::const_iterator itPart = plist.begin(), pend = plist.end(); // iterator to pairs (track id, particle)
	
         for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart)
	 {
            const simb::MCParticle* part = (itPart++)->second;
            if (!part)
	    {
            throw art::Exception(art::errors::LogicError)<< "GEANT particle #" << iPart << " returned a null pointer";
            }       
	  _ftrueparID[ntruepart]=part->TrackId();
	  _ftrueparpdg[ntruepart]=part->PdgCode();
	  _ftrueparStatusCode[ntruepart]=part->StatusCode();
	  _ftrueparTheta[ntruepart]=part->Momentum().Theta();
	  _ftrueparCosTheta[ntruepart]=TMath::Cos(part->Momentum().Theta());	
	  _ftrueparSinTheta[ntruepart]=TMath::Sin(part->Momentum().Theta());	
	  _ftrueparE[ntruepart]=part->E();
	  _ftrueparMass[ntruepart]=part->Mass();
	  _ftrueparKE[ntruepart]=(part->E())-(part->Mass());
	  _ftrueparEndE[ntruepart]=part->EndE();
          _ftrueparPx[ntruepart]=part->Px();
	  _ftrueparPy[ntruepart]=part->Py();
	  _ftrueparPz[ntruepart]=part->Pz();
	  _ftrueparP[ntruepart]=part->Momentum().Vect().Mag();
	  _ftrueparStartx[ntruepart]=part->Vx();
	  _ftrueparStarty[ntruepart]=part->Vy();
	  _ftrueparStartz[ntruepart]=part->Vz();
	  _ftrueparEndx[ntruepart]=part->EndPosition()[0];
	  _ftrueparEndy[ntruepart]=part->EndPosition()[1];
	  _ftrueparEndz[ntruepart]=part->EndPosition()[2];
	  _ftrueparSPcorrStartx[ntruepart]=part->Vx() - SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[0]; 
	  _ftrueparSPcorrStarty[ntruepart]=part->Vy() + SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[1];
	  _ftrueparSPcorrStartz[ntruepart]=part->Vz() + SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[2]; 
	  _ftrueparSPcorrEndx[ntruepart]=part->EndPosition()[0] - SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[0]; 
	  _ftrueparSPcorrEndy[ntruepart]=part->EndPosition()[1] + SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[1];
	  _ftrueparSPcorrEndz[ntruepart]=part->EndPosition()[2] + SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[2]; 
	  _ftrueparPhi[ntruepart]=part->Momentum().Phi();
	  _ftrueparCosPhi[ntruepart]=TMath::Cos(part->Momentum().Phi());
	  _ftrueparSinPhi[ntruepart]=TMath::Sin(part->Momentum().Phi());
            //---- track acceptance------------// 
//	    if(TMath::Abs(TMath::Cos(part.Momentum().Theta()))>=0.052 && part.StatusCode()==1) 
	    if(part->StatusCode()==1) 
	    {
	    if(Xnegbound<(_fTruenuvrtxx+(4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))*(_ftrueparSinTheta[ntruepart])*(_ftrueparCosPhi[ntruepart])) && Xposbound>(_fTruenuvrtxx+(4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))*(_ftrueparSinTheta[ntruepart])*(_ftrueparCosPhi[ntruepart]))
	    && Ynegbound<(_fTruenuvrtxy+(4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))*(_ftrueparSinTheta[ntruepart])*(_ftrueparSinPhi[ntruepart])) && Yposbound>(_fTruenuvrtxy+(4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))*(_ftrueparSinTheta[ntruepart])*(_ftrueparSinPhi[ntruepart]))
	    && Znegbound<(_fTruenuvrtxz+(4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))*(_ftrueparCosTheta[ntruepart])) && Zposbound>(_fTruenuvrtxz+(4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))*(_ftrueparCosTheta[ntruepart])))
	    {
	      if((TMath::Abs(_ftrueparpdg[ntruepart])==13 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))) 
	      || (TMath::Abs(_ftrueparpdg[ntruepart])==211 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))) 
	      || (_ftrueparpdg[ntruepart]==2212 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))) 
	      || (TMath::Abs(_ftrueparpdg[ntruepart])==321 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],4.5/TMath::Abs(_ftrueparCosTheta[ntruepart]))))
	      {
		trueInAccpmult++;
		if(TMath::Abs(part->PdgCode())==13)trueInAccpmultmu++; 
		if(TMath::Abs(part->PdgCode())==211) trueInAccpmultpi++;
		if(part->PdgCode()==2212) trueInAccpmultp++;
		if(TMath::Abs(part->PdgCode())==321) trueInAccpmultk++;

	        _ftrueInAccppargeantID[trueInAccep]=part->TrackId();
	        _ftrueInAccpparpdg[trueInAccep]=part->PdgCode();
	        _ftrueInAccpparStatusCode[trueInAccep]=part->StatusCode();
	        _ftrueInAccpparTheta[trueInAccep]=part->Momentum().Theta();
	        _ftrueInAccpparCosTheta[trueInAccep]=TMath::Cos(part->Momentum().Theta());	
	        _ftrueInAccpparSinTheta[trueInAccep]=TMath::Sin(part->Momentum().Theta());	
	        _ftrueInAccpparE[trueInAccep]=part->E();
	        _ftrueInAccpparMass[trueInAccep]=part->Mass();
	        _ftrueInAccpparKE[trueInAccep]=(part->E())-(part->Mass());
	        _ftrueInAccpparEndE[trueInAccep]=part->EndE();
                _ftrueInAccpparPx[trueInAccep]=part->Px();
	        _ftrueInAccpparPy[trueInAccep]=part->Py();
	        _ftrueInAccpparPz[trueInAccep]=part->Pz();
	        _ftrueInAccpparP[trueInAccep]=part->Momentum().Vect().Mag();
	        _ftrueInAccpparStartx[trueInAccep]=part->Vx();
	        _ftrueInAccpparStarty[trueInAccep]=part->Vy();
	        _ftrueInAccpparStartz[trueInAccep]=part->Vz();
	        _ftrueInAccpparEndx[trueInAccep]=part->EndPosition()[0];
	        _ftrueInAccpparEndy[trueInAccep]=part->EndPosition()[1];
	        _ftrueInAccpparEndz[trueInAccep]=part->EndPosition()[2];
	        _ftrueInAccpparSPcorrStartx[trueInAccep]=part->Vx() - SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[0]; 
	        _ftrueInAccpparSPcorrStarty[trueInAccep]=part->Vy() + SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[1];
	        _ftrueInAccpparSPcorrStartz[trueInAccep]=part->Vz() + SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[2]; 
	        _ftrueInAccpparSPcorrEndx[trueInAccep]=part->EndPosition()[0] - SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[0]; 
	        _ftrueInAccpparSPcorrEndy[trueInAccep]=part->EndPosition()[1] + SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[1];
	        _ftrueInAccpparSPcorrEndz[trueInAccep]=part->EndPosition()[2] + SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[2]; 
	        _ftrueInAccpparPhi[trueInAccep]=part->Momentum().Phi();
	        _ftrueInAccpparCosPhi[trueInAccep]=TMath::Cos(part->Momentum().Phi());
	        _ftrueInAccpparSinPhi[trueInAccep]=TMath::Sin(part->Momentum().Phi());
		trueInAccep++;
	      }
	    }
	    
	    if(_ftrueparCosPhi[ntruepart]>=0)    sx = (Xposbound-_fTruenuvrtxx)/(_ftrueparSinTheta[ntruepart]*_ftrueparCosPhi[ntruepart]);
	    if(_ftrueparCosPhi[ntruepart]<0)     sx = (Xnegbound-_fTruenuvrtxx)/(_ftrueparSinTheta[ntruepart]*_ftrueparCosPhi[ntruepart]);
	    if(_ftrueparSinPhi[ntruepart]>=0)    sy = (Yposbound-_fTruenuvrtxy)/(_ftrueparSinTheta[ntruepart]*_ftrueparSinPhi[ntruepart]);
	    if(_ftrueparSinPhi[ntruepart]<0)	 sy = (Ynegbound-_fTruenuvrtxy)/(_ftrueparSinTheta[ntruepart]*_ftrueparSinPhi[ntruepart]);
	    if(_ftrueparCosTheta[ntruepart]>=0)  sz = (Zposbound-_fTruenuvrtxz)/(_ftrueparCosTheta[ntruepart]);
	    if(_ftrueparCosTheta[ntruepart]<0)   sz = (Znegbound-_fTruenuvrtxz)/(_ftrueparCosTheta[ntruepart]);
	    
	    smin = TMath::Min(TMath::Min(sx,sy),sz);
	    //std::cout<<"sx: "<<sx<<" sy: "<<sy<<" sz: "<<sz<<" smin: "<<smin<<endl;
	    	             
	      if((TMath::Abs(_ftrueparpdg[ntruepart])==13 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],TMath::Max(lengthnumu,24/TMath::Abs(_ftrueparCosTheta[ntruepart]))) && _ftrueparKE[ntruepart]<T(_ftrueparpdg[ntruepart],smin)) 
	        || (TMath::Abs(_ftrueparpdg[ntruepart])==211 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],TMath::Max(lengthnumu,24/TMath::Abs(_ftrueparCosTheta[ntruepart]))) && _ftrueparKE[ntruepart]<=T(_ftrueparpdg[ntruepart],smin)) 
	        || (_ftrueparpdg[ntruepart]==2212 &&_ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],TMath::Max(lengthnumu,24/TMath::Abs(_ftrueparCosTheta[ntruepart])))&& _ftrueparKE[ntruepart]<=T(_ftrueparpdg[ntruepart],smin)) 
	        || (TMath::Abs(_ftrueparpdg[ntruepart])==321 && _ftrueparKE[ntruepart]>=T(_ftrueparpdg[ntruepart],TMath::Max(lengthnumu,24/TMath::Abs(_ftrueparCosTheta[ntruepart]))) && _ftrueparKE[ntruepart]<=T(_ftrueparpdg[ntruepart],smin)))
	 	{
		if(_ftrueparStartx[ntruepart]>longXnegbound && _ftrueparEndx[ntruepart]>longXnegbound && _ftrueparStartx[ntruepart]<longXposbound && _ftrueparEndx[ntruepart]<longXposbound
		&& _ftrueparStarty[ntruepart]>longYnegbound && _ftrueparEndy[ntruepart]>longYnegbound && _ftrueparStarty[ntruepart]<longYposbound && _ftrueparEndy[ntruepart]<longYposbound
		&& _ftrueparStartz[ntruepart]>longZnegbound && _ftrueparEndz[ntruepart]>longZnegbound && _ftrueparStartz[ntruepart]<longZposbound && _ftrueparEndz[ntruepart]<longZposbound)
		{
		  TrueInAccplongcontained++;

	        _ftruelongInAccppargeantID[truelongInAccep]=part->TrackId();
	        _ftruelongInAccpparpdg[truelongInAccep]=part->PdgCode();
	        _ftruelongInAccpparStatusCode[truelongInAccep]=part->StatusCode();
	        _ftruelongInAccpparTheta[truelongInAccep]=part->Momentum().Theta();
	        _ftruelongInAccpparCosTheta[truelongInAccep]=TMath::Cos(part->Momentum().Theta());	
	        _ftruelongInAccpparSinTheta[truelongInAccep]=TMath::Sin(part->Momentum().Theta()); 
	        _ftruelongInAccpparE[truelongInAccep]=part->E();
	        _ftruelongInAccpparMass[truelongInAccep]=part->Mass();
	        _ftruelongInAccpparKE[truelongInAccep]=(part->E())-(part->Mass());
	        _ftruelongInAccpparEndE[truelongInAccep]=part->EndE();
                _ftruelongInAccpparPx[truelongInAccep]=part->Px();
	        _ftruelongInAccpparPy[truelongInAccep]=part->Py();
	        _ftruelongInAccpparPz[truelongInAccep]=part->Pz();
	        _ftruelongInAccpparP[truelongInAccep]=part->Momentum().Vect().Mag();
	        _ftruelongInAccpparStartx[truelongInAccep]=part->Vx();
	        _ftruelongInAccpparStarty[truelongInAccep]=part->Vy();
	        _ftruelongInAccpparStartz[truelongInAccep]=part->Vz();
	        _ftruelongInAccpparEndx[truelongInAccep]=part->EndPosition()[0];
	        _ftruelongInAccpparEndy[truelongInAccep]=part->EndPosition()[1];
	        _ftruelongInAccpparEndz[truelongInAccep]=part->EndPosition()[2];
	        _ftruelongInAccpparSPcorrStartx[truelongInAccep]=part->Vx() - SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[0]; 
	        _ftruelongInAccpparSPcorrStarty[truelongInAccep]=part->Vy() + SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[1];
	        _ftruelongInAccpparSPcorrStartz[truelongInAccep]=part->Vz() + SCE->GetPosOffsets(part->Vx(),part->Vy(),part->Vz())[2]; 
	        _ftruelongInAccpparSPcorrEndx[truelongInAccep]=part->EndPosition()[0] - SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[0]; 
	        _ftruelongInAccpparSPcorrEndy[truelongInAccep]=part->EndPosition()[1] + SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[1];
	        _ftruelongInAccpparSPcorrEndz[truelongInAccep]=part->EndPosition()[2] + SCE->GetPosOffsets(part->EndPosition()[0],part->EndPosition()[1],part->EndPosition()[2])[2]; 
	        _ftruelongInAccpparPhi[truelongInAccep]=part->Momentum().Phi();
	        _ftruelongInAccpparCosPhi[truelongInAccep]=TMath::Cos(part->Momentum().Phi());
	        _ftruelongInAccpparSinPhi[truelongInAccep]=TMath::Sin(part->Momentum().Phi());
		truelongInAccep++;
		} 
		} 
	      }//if(part->StatusCode()==1) 
	     ++ntruepart;
	  } // particle loop
        } //if neutrino set
        _fNtrueInAccpmult=trueInAccpmult;
	_fNtrueInAccpmultmu=trueInAccpmultmu;
	_fNtrueInAccpmultpi=trueInAccpmultpi;
	_fNtrueInAccpmultp=trueInAccpmultp;
	_fNtrueInAccpmultk=trueInAccpmultk;
	_fNtruelongInAccpmult=TrueInAccplongcontained;

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
	    fDataTreePassLongcolhits->Fill();
	    

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
	  //std::cout<<i<<" upward hit_wire: "<<hits_wire[i]<<" "<<"hit_charge: "<<hits_charge[i]<<"\n"; 
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
	 //std::cout<<i<<" downward hit_wire: "<<hits_wire[i]<<" "<<"hit_charge: "<<hits_charge[i]<<"\n";
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
/*    int MCS=-999;
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
*/

    int MCS=-999;
    if(_flongtrackmcsfwdll < _flongtrackmcsbwdll) MCS=1;
    if(_flongtrackmcsfwdll >= _flongtrackmcsbwdll) MCS=-1;
    //std::cout<<"MCS log Likelihood ratio "<<_flongtrackmcsbwdll/_flongtrackmcsfwdll<<"\n";

    _fMCSdiff = (_flongtrackmcsbwdll)-(_flongtrackmcsfwdll);
    _fMCSratio = _flongtrackmcsbwdll/_flongtrackmcsfwdll;
  
    _fMCS = MCS;
   
  std::cout<<"MCS = "<<MCS<<"\n"; 
  
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
  } // end CTMMCAna::analyze()
   	
   DEFINE_ART_MODULE(CTMMCAna)
 
 } // namespace CTMMCAna_module
 
  #endif //CTMMCANA_H
  

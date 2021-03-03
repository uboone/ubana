////////////////////////////////////////////////////////////////////////
// Class:       ParticleIdAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        ParticleIdAnalyzer_module.cc
//
// Generated at Fri Feb  2 10:56:38 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

// data products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// UBXSec includes
//#include "uboone/UBXSec/DataTypes/SelectionResult.h"
//#include "uboone/UBXSec/DataTypes/TPCObject.h"

// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/fiducialVolume.h"
#include "uboone/ParticleID/Algorithms/PIDA.h"
#include "uboone/ParticleID/Algorithms/VarunaGetdEdx.h"
#include "uboone/ParticleID/Algorithms/Theory_dEdx_resrange.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>

class ParticleIdAnalyzer;


class ParticleIdAnalyzer : public art::EDAnalyzer {
  public:
    explicit ParticleIdAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ParticleIdAnalyzer(ParticleIdAnalyzer const &) = delete;
    ParticleIdAnalyzer(ParticleIdAnalyzer &&) = delete;
    ParticleIdAnalyzer & operator = (ParticleIdAnalyzer const &) = delete;
    ParticleIdAnalyzer & operator = (ParticleIdAnalyzer &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    std::vector<double> fv;

  private:

    art::ServiceHandle< art::TFileService > tfs; 

    std::string fTrackingAlgo;
    std::string fCaloLabel;

    double fCutDistance;
    double fCutFraction;

    bool fUseLibosSelection;
    std::string fLibosSelectionFile;

    // fidvol related
    fidvol::fiducialVolume fid;
    particleid::PIDA pida;

    bool fUseNewCalib;

    // ubxsec related
    //ubana::SelectionResult selRes;

    // cc1munp related
    int cc1munpRun;
    int cc1munpSubRun;
    int cc1munpEvent;
    int cc1munpMuonID;
    int cc1munpProtonID;
    TTree* ttree;

    int nProtons = 0;

    TH1D* hMuonPreCutMean;
    TH1D* hMuonPreCutMedian;
    TH1D* hMuonPreCutKde;
    TH1D* hProtonPreCutMean;
    TH1D* hProtonPreCutMedian;
    TH1D* hProtonPreCutKde;

    TH1D* hPreCutMean;
    TH1D* hPreCutMedian;
    TH1D* hPreCutKde;

    TH1D* hMuonPostCutMean;
    TH1D* hMuonPostCutMedian;
    TH1D* hMuonPostCutKde;
    TH1D* hProtonPostCutMean;
    TH1D* hProtonPostCutMedian;
    TH1D* hProtonPostCutKde;

    TH1D* hPostCutMean;
    TH1D* hPostCutMedian;
    TH1D* hPostCutKde;

    TH2D* hProtonStartYZ;
    TH2D* hMuonStartYZ;

    TH1D* hProtonAlldQdx_uncalib;
    TH1D* hMuonAlldQdx_uncalib;
    TH1D* hProtonAlldQdx_oldcalib;
    TH1D* hMuonAlldQdx_oldcalib;
    TH1D* hProtonAlldQdx_newcalib;
    TH1D* hMuonAlldQdx_newcalib;
  
    TH1D* hAlldQdx_uncalib;
    TH1D* hAlldQdx_oldcalib;
    TH1D* hAlldQdx_newcalib;
  
    TH2D* hProtondQdx_resrange_uncalib;
    TH2D* hMuondQdx_resrange_uncalib;
    TH2D* hProtondQdx_resrange_oldcalib;
    TH2D* hMuondQdx_resrange_oldcalib;
    TH2D* hProtondQdx_resrange_newcalib;
    TH2D* hMuondQdx_resrange_newcalib;
  
    TH2D* hdQdx_resrange_uncalib;
    TH2D* hdQdx_resrange_oldcalib;
    TH2D* hdQdx_resrange_newcalib;

    TH1D* hProtonAlldEdx_oldcalib;
    TH1D* hMuonAlldEdx_oldcalib;
    TH1D* hProtonAlldEdx_newcalib;
    TH1D* hMuonAlldEdx_newcalib;
  
    TH1D* hAlldEdx_oldcalib;
    TH1D* hAlldEdx_newcalib;
  
    TH2D* hProtondEdx_resrange_oldcalib;
    TH2D* hMuondEdx_resrange_oldcalib;
    TH2D* hProtondEdx_resrange_newcalib;
    TH2D* hMuondEdx_resrange_newcalib;
  
    TH2D* hdEdx_resrange_oldcalib;
    TH2D* hdEdx_resrange_newcalib;

    TGraph *g_ThdEdxRR_Proton;
    TGraph *g_ThdEdxRR_Kaon;
    TGraph *g_ThdEdxRR_Pion;
    TGraph *g_ThdEdxRR_Muon;

    std::vector<TH1D*> protonPIDAVals;

};


ParticleIdAnalyzer::ParticleIdAnalyzer(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");
  fCaloLabel = p.get< std::string > ("CalorimetryModule");
  fCutDistance = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction = p.get< double > ("DaughterFinderCutFraction");
  fUseLibosSelection = p.get< bool >("UseLibosSelection",true);
  fLibosSelectionFile = p.get< std::string >("LibosSelectionFile","");
  fUseNewCalib = p.get< bool >("UseNewCalib",true);

  std::cout << "fUseLibosSelection = " << fUseLibosSelection << std::endl;

  fv = fid.setFiducialVolume(fv, p);
  fid.printFiducialVolume(fv);

}

void ParticleIdAnalyzer::beginJob()
{
  if (fUseLibosSelection){
    std::cout << "Opening file..." << std::endl;
    TFile *file = new TFile(fLibosSelectionFile.c_str(), "READ");
    std::cout << "...Done." << std::endl;
    
    if (file->IsZombie())
      std::cout << "File is zombie!" << std::endl;
 
    ttree = (TTree*)file->Get("cc1unpselana/fMC_TrunMean");
    
    ttree->SetBranchAddress("fRun", &cc1munpRun);
    ttree->SetBranchAddress("fSubRun", &cc1munpSubRun);
    ttree->SetBranchAddress("fEvent", &cc1munpEvent);
    ttree->SetBranchAddress("trkmuoncandid", &cc1munpMuonID);
    ttree->SetBranchAddress("trkprotoncandid", &cc1munpProtonID);
  }
    
  hMuonPreCutMean = tfs->make<TH1D>("hMuonPreCutMean", ";;", 100, 0, 30);
  hMuonPreCutMedian = tfs->make<TH1D>("hMuonPreCutMedian", ";;", 100, 0, 30);
  hMuonPreCutKde = tfs->make<TH1D>("hMuonPreCutKde", ";;", 100, 0, 30);
  hProtonPreCutMean = tfs->make<TH1D>("hProtonPreCutMean", ";;", 100, 0, 30);
  hProtonPreCutMedian = tfs->make<TH1D>("hProtonPreCutMedian", ";;", 100, 0, 30);
  hProtonPreCutKde = tfs->make<TH1D>("hProtonPreCutKde", ";;", 100, 0, 30);

  hPreCutMean = tfs->make<TH1D>("hPreCutMean", ";;", 100, 0, 30);
  hPreCutMedian = tfs->make<TH1D>("hPreCutMedian", ";;", 100, 0, 30);
  hPreCutKde = tfs->make<TH1D>("hPreCutKde", ";;", 100, 0, 30);

  hMuonPostCutMean = tfs->make<TH1D>("hMuonPostCutMean", ";;", 100, 0, 30);
  hMuonPostCutMedian = tfs->make<TH1D>("hMuonPostCutMedian", ";;", 100, 0, 30);
  hMuonPostCutKde = tfs->make<TH1D>("hMuonPostCutKde", ";;", 100, 0, 30);
  hProtonPostCutMean = tfs->make<TH1D>("hProtonPostCutMean", ";;", 100, 0, 30);
  hProtonPostCutMedian = tfs->make<TH1D>("hProtonPostCutMedian", ";;", 100, 0, 30);
  hProtonPostCutKde = tfs->make<TH1D>("hProtonPostCutKde", ";;", 100, 0, 30);

  hPostCutMean = tfs->make<TH1D>("hPostCutMean", ";;", 100, 0, 30);
  hPostCutMedian = tfs->make<TH1D>("hPostCutMedian", ";;", 100, 0, 30);
  hPostCutKde = tfs->make<TH1D>("hPostCutKde", ";;", 100, 0, 30);

  hProtonStartYZ = tfs->make<TH2D>("hProtonStartYZ", ";;", 50, 0, 1036, 50, -116.5, 116.5);
  hMuonStartYZ = tfs->make<TH2D>("hMuonStartYZ", ";;", 50, 0, 1036, 50, -116.5, 116.5);

  // dQdx 1D
  hProtonAlldQdx_uncalib = tfs->make<TH1D>("hProtonAlldQdx_uncalib","Uncalibrated (Proton candidates);Total dQ/dx (ADC/cm);No. hits",150,0,1500);
  hMuonAlldQdx_uncalib = tfs->make<TH1D>("hMuonAlldQdx_uncalib","Uncalibrated (Muon candidates);Total dQ/dx (ADC/cm);No. hits",150,0,1500);
  hProtonAlldQdx_oldcalib = tfs->make<TH1D>("hProtonAlldQdx_oldcalib","Old calibration (Proton candidates);Total dQ/dx (e^{-}/cm);No. hits",500,0,300e3);
  hMuonAlldQdx_oldcalib = tfs->make<TH1D>("hMuonAlldQdx_oldcalib","Old calibration (Muon candidates);Total dQ/dx (e^{-}/cm);No. hits",500,0,300e3);
  hProtonAlldQdx_newcalib = tfs->make<TH1D>("hProtonAlldQdx_newcalib","New MCC 8.7 calibration (Proton candidates);Total dQ/dx (e^{-}/cm);No. hits",150,0,1500);
  hMuonAlldQdx_newcalib = tfs->make<TH1D>("hMuonAlldQdx_newcalib","New MCC 8.7 calibration (Muon candidates);Total dQ/dx (e^{-}/cm);No. hits",150,0,1500);
  
  hAlldQdx_uncalib = tfs->make<TH1D>("hAlldQdx_uncalib","Uncalibrated (All tracks);Total dQ/dx (ADC/cm);No. hits",150,0,1500);
  hAlldQdx_oldcalib = tfs->make<TH1D>("hAlldQdx_oldcalib","Old calibration (All tracks);Total dQ/dx (e^{-}/cm);No. hits",500,0,300e3);
  hAlldQdx_newcalib = tfs->make<TH1D>("hAlldQdx_newcalib","New MCC 8.7 calibration (All tracks);Total dQ/dx (e^{-}/cm);No. hits",150,0,1500);


  // dQdx vs residual range
  hProtondQdx_resrange_uncalib = tfs->make<TH2D>("hProtondQdx_resrange_uncalib","Uncalibrated (Proton candidates);Residual range (cm);dQ/dx (ADC/cm)",500,0,50,500,0,1500);
  hMuondQdx_resrange_uncalib = tfs->make<TH2D>("hMuondQdx_resrange_uncalib","Uncalibrated (Muon candidates);Residual range (cm);dQ/dx (ADC/cm)",500,0,50,500,0,1500);
  hProtondQdx_resrange_oldcalib = tfs->make<TH2D>("hProtondQdx_resrange_oldcalib","Old calibration (Proton candidates);Residual range (cm);dQ/dx (e^{-}/cm)",500,0,50,500,0,300e3);
  hMuondQdx_resrange_oldcalib = tfs->make<TH2D>("hMuondQdx_resrange_oldcalib","Old calibration (Muon candidates);Residual range (cm);dQ/dx (e^{-}/cm)",500,0,50,500,0,300e3);
  hProtondQdx_resrange_newcalib = tfs->make<TH2D>("hProtondQdx_resrange_newcalib","New MCC 8.7 calibration (Proton candidates);Residual range (cm);dQ/dx (e^{-}/cm)",500,0,50,500,0,1500);
  hMuondQdx_resrange_newcalib = tfs->make<TH2D>("hMuondQdx_resrange_newcalib","New MCC 8.7 calibration (Muon candidates);Residual range (cm);dQ/dx (e^{-}/cm)",500,0,50,500,0,1500);
  
  hdQdx_resrange_uncalib = tfs->make<TH2D>("hdQdx_resrange_uncalib","Uncalibrated (All tracks);Residual range (cm);dQ/dx (ADC/cm)",500,0,50,500,0,1500);
  hdQdx_resrange_oldcalib = tfs->make<TH2D>("hdQdx_resrange_oldcalib","Old calibration (All tracks);Residual range (cm);dQ/dx (e^{-}/cm)",500,0,50,500,0,300e3);
  hdQdx_resrange_newcalib = tfs->make<TH2D>("hdQdx_resrange_newcalib","New MCC 8.7 calibration (All tracks);Residual range (cm);dQ/dx (e^{-}/cm)",500,0,50,500,0,1500);


  // dEdx 1D
  hProtonAlldEdx_oldcalib = tfs->make<TH1D>("hProtonAlldEdx_oldcalib","Old calibration (Proton candidates);Total dE/dx (units?) calculated from dQ/dx;No. hits",100,0,50);
  hMuonAlldEdx_oldcalib = tfs->make<TH1D>("hMuonAlldEdx_oldcalib","Old calibration (Muon candidates);Total dE/dx (units?) calculated from dQ/dx;No. hits",100,0,50);
  hProtonAlldEdx_newcalib = tfs->make<TH1D>("hProtonAlldEdx_newcalib","New MCC 8.7 calibration (Proton candidates);Total dE/dx (units?) from cali data product;No. hits",100,0,50);
  hMuonAlldEdx_newcalib = tfs->make<TH1D>("hMuonAlldEdx_newcalib","New MCC 8.7 calibration (Muon candidates);Total dE/dx (units?) from cali data product;No. hits",100,0,50);
  
  hAlldEdx_oldcalib = tfs->make<TH1D>("hAlldEdx_oldcalib","Old calibration (All tracks);Total dE/dx (units?) calculated from dQ/dx;No. hits",100,0,50);
  hAlldEdx_newcalib = tfs->make<TH1D>("hAlldEdx_newcalib","New MCC 8.7 calibration (All tracks);Total dE/dx (units?) from cali data product;No. hits",10,0,50);


  // dEdx vs residual range
  hProtondEdx_resrange_oldcalib = tfs->make<TH2D>("hProtondEdx_resrange_oldcalib","Old calibration (Proton candidates);Residual range (cm);dE/dx (units?) calculated from dQ/dx",500,0,50,100,0,50);
  hMuondEdx_resrange_oldcalib = tfs->make<TH2D>("hMuondEdx_resrange_oldcalib","Old calibration (Muon candidates);Residual range (cm);dE/dx (units?) calculated from dQ/dx",500,0,50,100,0,50);
  hProtondEdx_resrange_newcalib = tfs->make<TH2D>("hProtondEdx_resrange_newcalib","New MCC 8.7 calibration (Proton candidates);Residual range (cm);dE/dx (units?) from cali data product",500,0,50,100,0,50);
  hMuondEdx_resrange_newcalib = tfs->make<TH2D>("hMuondEdx_resrange_newcalib","New MCC 8.7 calibration (Muon candidates);Residual range (cm);dE/dx (units?) from cali data product",500,0,50,100,0,50);
  
  hdEdx_resrange_oldcalib = tfs->make<TH2D>("hdEdx_resrange_oldcalib","Old calibration (All tracks);Residual range (cm);dE/dx (units?) calculated from dQ/dx",500,0,50,100,0,50);
  hdEdx_resrange_newcalib = tfs->make<TH2D>("hdEdx_resrange_newcalib","New MCC 8.7 calibration (All tracks);Residual range (cm);dE/dx (units?) from cali data product",500,0,50,100,0,50);


  // Theoretical predictions
  particleid::Theory_dEdx_resrange theorypreds;
  //int npoints = theorypreds.GetNPoints();
  // Why is this not saving?
  g_ThdEdxRR_Proton = tfs->make<TGraph>();
  g_ThdEdxRR_Proton = (TGraph*)theorypreds.g_ThdEdxRR_Proton->Clone("Theory_proton");
  g_ThdEdxRR_Kaon = tfs->make<TGraph>();
  g_ThdEdxRR_Kaon = (TGraph*)theorypreds.g_ThdEdxRR_Kaon->Clone("Theory_kaon");
  g_ThdEdxRR_Pion = tfs->make<TGraph>();
  g_ThdEdxRR_Pion = (TGraph*)theorypreds.g_ThdEdxRR_Pion->Clone("Theory_pion");
  g_ThdEdxRR_Muon = tfs->make<TGraph>();
  g_ThdEdxRR_Muon = (TGraph*)theorypreds.g_ThdEdxRR_Muon->Clone("Theory_muon");

  for (int i = 0; i < 50; i ++){

    TString th1name = Form("protonPIDAVals_%i", i);
    protonPIDAVals.push_back(tfs->make<TH1D>(th1name, ";PIDA vals;", 1000, 0, 50));

  }

  std::cout << "Histograms created." << std::endl;
}

void ParticleIdAnalyzer::analyze(art::Event const & e)
{

  bool isData = e.isRealData();
  bool isSelected = false;

  int fRun = e.run();
  int fSubRun = e.subRun();
  int fEvent = e.event();

  std::cout << "----- " << fRun << "." << fSubRun << "." << fEvent << std::endl;

  std::vector<int> muonIds;
  std::vector<int> protonIds;

  // If using Libo's selection, look for a selected event in this file and get the muon and proton
  // track IDs
  if (fUseLibosSelection){
    
    for (int i = 0; i < ttree->GetEntries(); i++){
      
      ttree->GetEntry(i);
      
      
      if (cc1munpRun == fRun && cc1munpSubRun == fSubRun && cc1munpEvent == fEvent){
	
	std::cout << "Found an event!" << std::endl;
	isSelected = true;
	muonIds.push_back(cc1munpMuonID);
	protonIds.push_back(cc1munpProtonID);
	continue;
      }
      
    }

  if (isSelected == false) return;
  
  }
  
  //
  // get handles to needed information
  //

  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackingAlgo, trackHandle);
  std::vector< art::Ptr< recob::Track > > trackPtrs;
  art::fill_ptr_vector(trackPtrs, trackHandle);

  art::FindManyP< anab::Calorimetry > caloFromTracks(trackHandle, e, "pandoraNucalo");
  std::cout << "trackPtrs.size() = " << trackPtrs.size() << std::endl;
  std::cout << "caloFromTracks.size() = " << caloFromTracks.size() << std::endl;
  for (unsigned int i=0; i< caloFromTracks.size(); i++){
    std::cout << "caloFromTracks.at(" << i << ").size() = " << caloFromTracks.at(i).size() << std::endl;
  }
  art::FindManyP< anab::Calorimetry > *caliFromTracks=nullptr;
  art::Handle< std::vector< anab::Calorimetry > > calihandle;
  std::vector< art::Ptr< anab::Calorimetry > > caliPtrs;
  if (fUseNewCalib){
    caliFromTracks = new art::FindManyP< anab::Calorimetry >(trackHandle, e, "pandoraNucali");
    std::cout << "caliFromTracks->size() = " << caliFromTracks->size() << std::endl;
    for (unsigned int i=0; i< caliFromTracks->size(); i++){
      std::cout << "caliFromTracks->at(" << i << ").size() = " << caliFromTracks->at(i).size() << std::endl;
    }
    
    e.getByLabel("pandoraNucali", calihandle);
    art::fill_ptr_vector(caliPtrs, calihandle);
    std::cout << "caliPtrs.size() = " << caliPtrs.size() << std::endl;
  }

  for (size_t j = 0; j < trackPtrs.size(); j++){

    art::Ptr< recob::Track > track     = trackPtrs.at(j);
    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());
    art::Ptr< anab::Calorimetry > calo;
    for (auto c : caloFromTrack){
      if (!c) continue; // avoid art errors if c doesn't exist
      if (!c->PlaneID().isValid) continue; // check has valid plane
      int planenum = c->PlaneID().Plane;
      if (planenum != 2) continue; // only use calorimetry from collection plane (plane 2)
      calo = c;
    }
    std::cout << "CALO USING PLANE: " << calo->PlaneID() << std::endl;

    std::vector< double > dEdx         = calo->dEdx();
    std::vector< double > dQdx         = calo->dQdx();
    std::vector< double > resRange     = calo->ResidualRange();

    double oldcalibfactor;
    if (isData){ oldcalibfactor = 243; } // Data: multiply by 243
    else { oldcalibfactor = 198; } // MC: multiply by 198

    double VarunaCalibfactor = 0.00507669;
    // Used in Varuna's conversion from dQdx to dEdx, for plane 2 only

    std::vector< art::Ptr<anab::Calorimetry> > caliFromTrack;
    art::Ptr< anab::Calorimetry > cali;
    std::vector< double > dEdx_cali;
    std::vector< double > dQdx_cali;
    std::vector< double > resRange_cali;
    if (fUseNewCalib){
      caliFromTrack = caliFromTracks->at(track->ID());
      std::cout << track->ID() << std::endl;
      for (auto c : caliFromTrack){
	if (!c) continue; // avoid art errors if c doesn't exist
	if (!c->PlaneID().isValid) continue; // check has valid plane
	int planenum = c->PlaneID().Plane;
	if (planenum != 2) continue; // only use calorimetry from collection plane (plane 2)
	cali = c;
      }
      std::cout << "CALI USING PLANE: " << cali->PlaneID() << std::endl;

      dEdx_cali    = cali->dEdx();
      dQdx_cali    = cali->dQdx();
      resRange_cali= cali->ResidualRange();

      // Check: residual range should be the same in both data products
      if (resRange.size() != resRange_cali.size()){
	std::cout << "[ERROR] resRange from calo object size = " << resRange.size() << ", but resRange from cali object size = " << resRange_cali.size() << std::endl;
      }
    } // if caliFromTracks
    
    int trackID          = track->ID();
    int nDaughters       = GetNDaughterTracks((*trackHandle), trackID, fCutDistance, fCutFraction);

    TVector3 trackStart  = track->Vertex();
    TVector3 trackEnd    = track->End();

    double pidaValMean   = pida.getPida(dEdx, resRange, "mean");
    double pidaValMedian = pida.getPida(dEdx, resRange, "median");
    double pidaValKde    = pida.getPida(dEdx, resRange, "kde");

    // --- Fill plots for all particles ---
    // This code will execute whether you are using Libo's selection or not

    hPreCutMean->Fill(pidaValMean);
    hPreCutMedian->Fill(pidaValMean);
    hPreCutKde->Fill(pidaValMean);
    
    for (size_t j = 0; j < resRange.size(); j++){
      hdQdx_resrange_uncalib->Fill(resRange.at(j), dQdx.at(j));
      hdQdx_resrange_oldcalib->Fill(resRange.at(j), dQdx.at(j)*oldcalibfactor);
	
      hAlldQdx_uncalib->Fill(dQdx.at(j));
      hAlldQdx_oldcalib->Fill(dQdx.at(j)*oldcalibfactor);

      hdEdx_resrange_oldcalib->Fill(resRange.at(j), VarunaGetdEdx(dQdx.at(j)/VarunaCalibfactor));

      hAlldEdx_oldcalib->Fill(VarunaGetdEdx(dQdx.at(j)/VarunaCalibfactor));
      
      if (fUseNewCalib){
      hdQdx_resrange_newcalib->Fill(resRange_cali.at(j), dQdx_cali.at(j));
      hAlldQdx_newcalib->Fill(dQdx_cali.at(j));
      hdEdx_resrange_newcalib->Fill(resRange_cali.at(j), dEdx_cali.at(j));
      hAlldEdx_newcalib->Fill(dEdx_cali.at(j));
      }
    }
    
    if (nDaughters == 0 && fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){
      hPostCutMean->Fill(pidaValMean);
      hPostCutMedian->Fill(pidaValMean);
      hPostCutKde->Fill(pidaValMean);
    }

    // --- Fill plots for muon candidates ---
    // If not using Libo's selection, muonIds.size()==0 and this section is skipped
    for (size_t i = 0; i < muonIds.size(); i++){
      
      if (trackID == muonIds.at(i)){
	
	std::cout << ">> Found Candidate Muon!" << std::endl;
	
	hMuonPreCutMean->Fill(pidaValMean);
	hMuonPreCutMedian->Fill(pidaValMedian);
	hMuonPreCutKde->Fill(pidaValKde);
	
	hMuonStartYZ->Fill(trackStart.Z(), trackStart.Y());

        for (size_t j = 0; j < resRange.size(); j++){
          hMuondQdx_resrange_uncalib->Fill(resRange.at(j), dQdx.at(j));
          hMuondQdx_resrange_oldcalib->Fill(resRange.at(j), dQdx.at(j)*oldcalibfactor);
	
	  hMuonAlldQdx_uncalib->Fill(dQdx.at(j));
	  hMuonAlldQdx_oldcalib->Fill(dQdx.at(j)*oldcalibfactor);

	  hMuondEdx_resrange_oldcalib->Fill(resRange.at(j), VarunaGetdEdx(dQdx.at(j)/VarunaCalibfactor));
	  
	  hMuonAlldEdx_oldcalib->Fill(VarunaGetdEdx(dQdx.at(j)/VarunaCalibfactor));

	  if (fUseNewCalib){
	  hMuondQdx_resrange_newcalib->Fill(resRange_cali.at(j), dQdx_cali.at(j));
	  hMuonAlldQdx_newcalib->Fill(dQdx_cali.at(j));
	  hMuondEdx_resrange_newcalib->Fill(resRange_cali.at(j), dEdx_cali.at(j));
	  hMuonAlldEdx_newcalib->Fill(dEdx_cali.at(j));
	  }
        }
	
	if (nDaughters == 0 && fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){
	    
	  hMuonPostCutMean->Fill(pidaValMean);
	  hMuonPostCutMedian->Fill(pidaValMedian);
	  hMuonPostCutKde->Fill(pidaValKde);
	  
	} // end if (contained, no daughters)
	
      } // end if (trackID == muonIds.at(i))
      
    } // end loop over muonIds

    // --- Fill plots for proton candidates ---
    // If not using Libo's selection, protonIds.size()==0 and this section is skipped
    for (size_t i = 0; i < protonIds.size(); i++){
      
      if (trackID == protonIds.at(i)){
	
	std::cout << ">> Found Candidate Proton!" << std::endl;
	
	hProtonPreCutMean->Fill(pidaValMean);
	hProtonPreCutMedian->Fill(pidaValMedian);
	hProtonPreCutKde->Fill(pidaValKde);
	
	hProtonStartYZ->Fill(trackStart.Z(), trackStart.Y());

        for (size_t j = 0; j < resRange.size(); j++){
          hProtondQdx_resrange_uncalib->Fill(resRange.at(j), dQdx.at(j));
          hProtondQdx_resrange_oldcalib->Fill(resRange.at(j), dQdx.at(j)*oldcalibfactor);
	
	  hProtonAlldQdx_uncalib->Fill(dQdx.at(j));
	  hProtonAlldQdx_oldcalib->Fill(dQdx.at(j)*oldcalibfactor);

	  hProtondEdx_resrange_oldcalib->Fill(resRange.at(j), VarunaGetdEdx(dQdx.at(j)/VarunaCalibfactor));

	  hProtonAlldEdx_oldcalib->Fill(VarunaGetdEdx(dQdx.at(j)/VarunaCalibfactor));

	  if (fUseNewCalib){
	  hProtondQdx_resrange_newcalib->Fill(resRange_cali.at(j), dQdx_cali.at(j));
	  hProtonAlldQdx_newcalib->Fill(dQdx_cali.at(j));
	  hProtondEdx_resrange_newcalib->Fill(resRange_cali.at(j), dEdx_cali.at(j));
	  hProtonAlldEdx_newcalib->Fill(dEdx_cali.at(j));
	  }
        }
	
	if (nDaughters == 0 && fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){
	  
	  hProtonPostCutMean->Fill(pidaValMean);
	  hProtonPostCutMedian->Fill(pidaValMedian);
	  hProtonPostCutKde->Fill(pidaValKde);
	  
	  if (nProtons < 50){
	    TString th1name = Form("run%i_event%i_trkid%i", fRun, fEvent, trackID);
	    protonPIDAVals.at(nProtons)->SetName(th1name);
	    for (size_t j = 0; j < dEdx.size(); j++){
	      
	      std::cout << "dEdx "  << j << ": " << dEdx.at(j) << " ResRg: " << resRange.at(j) <<  " particleId: " << dEdx.at(j)*std::pow(resRange.at(j), 0.42) << std::endl;
	      protonPIDAVals.at(nProtons)->Fill(dEdx.at(j)*std::pow(resRange.at(j),0.42));
	    } // end loop over dEdx values
	    

            nProtons++;
          } // end if (nProtons < 50)
        } // end if (contained, no daughters)
	
      } // end if (trackID == protonIds.at(i))
      
    } // end loop over protonIds
    
  } // end loop over tracks
}

void ParticleIdAnalyzer::endJob()
{
}

DEFINE_ART_MODULE(ParticleIdAnalyzer)

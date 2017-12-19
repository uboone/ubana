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
//#include "lardata/RecoObjects/BezierTrack.h"
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

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

using namespace std;

namespace microboone{

//========================================================================
// Length of reconstructed track, trajectory by trajectory.
double length(const recob::Track& track)
{
  double result = 0.;
  TVector3 disp = track.LocationAtPoint(0);
  int n = track.NumberTrajectoryPoints();

  for(int i = 1; i < n; ++i) {
    const TVector3& pos = track.LocationAtPoint(i);
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}

//========================================================================

//========================================================================
double length(const simb::MCParticle& part, TVector3& start, TVector3& end)
{
  // Get geometry.
  auto const* geom = lar::providerFrom<geo::Geometry>();
  
  // Get active volume boundary.
  double xmin = 0.;
  double xmax = 2.*geom->DetHalfWidth();
  double ymin = -geom->DetHalfHeight();
  double ymax = geom->DetHalfHeight();
  double zmin = 0.;
  double zmax = geom->DetLength();
  double vDrift = 160*pow(10,-6);

  double result = 0.;
  TVector3 disp;
  int n = part.NumberTrajectoryPoints();
  bool first = true;

  for(int i = 0; i < n; ++i) {
    // check if the particle is inside a TPC
    double mypos[3] = {part.Vx(i), part.Vy(i), part.Vz(i)};
    if (mypos[0] >= xmin && mypos[0] <= xmax && mypos[1] >= ymin && mypos[1] <= ymax && mypos[2] >= zmin && mypos[2] <= zmax){
      double xGen   = part.Vx(i);
      double tGen   = part.T(i);
      // Doing some manual shifting to account for an interaction not occuring with the beam dump
      // we will reconstruct an x distance different from where the particle actually passed to to the time
      // being different from in-spill interactions
      double newX = xGen+(tGen*vDrift);
      if (newX < -xmax || newX > (2*xmax)) continue;
     
      TVector3 pos(newX,part.Vy(i),part.Vz(i));
      if(first){
        start = pos;
      }
      else {
        disp -= pos;
        result += disp.Mag();
      }
      first = false;
      disp = pos;
      end = pos;
    }
  }
  return result;
}

class Kplane2 : public art::EDAnalyzer {
public:

    explicit Kplane2(fhicl::ParameterSet const& pset);
    virtual ~Kplane2();

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
    Int_t    reco_events;
    Int_t    true_events;
    Int_t    candidate_k; 
    Int_t    pri_mu_id;
    Int_t    pri_k_id;
    Int_t    dec_mu_id;
    Int_t    break_value;
    
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
    bool  fSaveVertexInfo;
    float fG4minE;
    
    TH2F *reco_kaon_cont_plane_two_res_range_dEdX_hist;
    TH2F *true_kaon_cont_plane_best_res_range_dEdX_hist;
    TH1F *reco_de_dx_gradient;
    TH1F *true_de_dx_gradient;
    TH1F *reco_k_mu_open_angle_cosine;
    TH1F *true_k_mu_open_angle_cosine;
    TH1F *reco_kaon_trkln_hist;
    TH1F *true_kaon_trkln_hist;
    TH1F *reco_muon_trkln_hist;
    TH1F *true_muon_trkln_hist;
    TH1F *reco_mu_mom_range;
    TH1F *true_mu_mom_range;
    TH2F *reco_muon_cont_plane_two_res_range_dEdX_hist;
    TH2F *true_muon_cont_plane_best_res_range_dEdX_hist;
    TH1F *reco_k_pida;
    TH1F *true_k_pida;
    TH1F *reco_mu_pida;
    TH1F *true_mu_pida;
    TH1F *reco_k_dedx_ratio;
    TH1F *true_k_dedx_ratio;
    TH1F *reco_mu_dedx_ratio;
    TH1F *true_mu_dedx_ratio;
    TH1F *reco_de_dx_median_gradient;
    TH1F *true_de_dx_median_gradient;
    TH1F *reco_k_median_dedx_ratio;
    TH1F *true_k_median_dedx_ratio;
    TH1F *reco_mu_median_dedx_ratio;
    TH1F *true_mu_median_dedx_ratio;
    TH1F *reco_k_mean_first_dedx;
    TH1F *true_k_mean_first_dedx;
    TH1F *reco_k_mean_last_dedx;
    TH1F *true_k_mean_last_dedx;
    TH1F *reco_mu_mean_first_dedx;
    TH1F *true_mu_mean_first_dedx;
    TH1F *reco_mu_mean_last_dedx;
    TH1F *true_mu_mean_last_dedx;
    TH1F *reco_k_median_first_dedx;
    TH1F *true_k_median_first_dedx;
    TH1F *reco_k_median_last_dedx;
    TH1F *true_k_median_last_dedx;
    TH1F *reco_mu_median_first_dedx;
    TH1F *true_mu_median_first_dedx;
    TH1F *reco_mu_median_last_dedx;
    TH1F *true_mu_median_last_dedx;
    TH1F *reco_k_mean_len_ratio;
    TH1F *true_k_mean_len_ratio;
    TH1F *reco_mu_mean_len_ratio;
    TH1F *true_mu_mean_len_ratio;
    TH1F *reco_k_median_len_ratio;
    TH1F *true_k_median_len_ratio;
    TH1F *reco_mu_median_len_ratio;
    TH1F *true_mu_median_len_ratio;
    TH1F *reco_all_trk_median_dedx;
    TH1F *true_all_trk_median_dedx;
    TH1F *reco_all_trmu_median_dedx;
    TH1F *true_all_trmu_median_dedx;
    TH1F *k_best_plane_hist;
    TH1F *mu_best_plane_hist;
    
    ////////////// Test ///////////
    TH1F *event_k_pida;
    TH1F *event_mu_pida;
    TH1F *event_k_trklen;
    TH1F *event_mu_trklen;
    TH1F *event_de_dx_median_gradient;
    TH1F *event_de_dx_mean_gradient;
    TH1F *event_k_mu_open_angle_cosine;
    TH1F *event_mu_mom_range;
    TH2F *event_kaon_cont_plane_best_res_range_dEdX_hist;
    //////////// End test ////////
    
 }; 

//========================================================================
Kplane2::Kplane2(fhicl::ParameterSet const& pset) :
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
  fSaveVertexInfo           (pset.get< bool>("SaveVertexInfo",false)),  
  fG4minE                   (pset.get< float>("G4minE",0.01))  
{
  if (fSaveTrackInfo == false) fSaveCaloInfo = false;
}
 
//========================================================================
Kplane2::~Kplane2(){
  //destructor
}
//========================================================================

//========================================================================
void Kplane2::beginJob(){
  std::cout<<"job begin..."<<std::endl;

  art::ServiceHandle<art::TFileService> tfs;
  
  reco_kaon_cont_plane_two_res_range_dEdX_hist = tfs->make<TH2F>("reco_kaon_cont_plane_two_res_range_dEdX_hist","",200,0,100,200,0,50);
  reco_kaon_cont_plane_two_res_range_dEdX_hist->GetXaxis()->SetTitle("Residual range R(cm)");
  reco_kaon_cont_plane_two_res_range_dEdX_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  true_kaon_cont_plane_best_res_range_dEdX_hist = tfs->make<TH2F>("true_kaon_cont_plane_best_res_range_dEdX_hist","",200,0,100,200,0,50);
  true_kaon_cont_plane_best_res_range_dEdX_hist->GetXaxis()->SetTitle("Residual range R(cm)");
  true_kaon_cont_plane_best_res_range_dEdX_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  reco_de_dx_gradient = tfs->make<TH1F>("reco_de_dx_gradient","",500,-100,1000);
  reco_de_dx_gradient->GetXaxis()->SetTitle("Energy gradient (MeV)");
  
  true_de_dx_gradient = tfs->make<TH1F>("true_de_dx_gradient","",500,-100,1000);
  true_de_dx_gradient->GetXaxis()->SetTitle("Energy gradient (MeV)");
  
  reco_k_mu_open_angle_cosine = tfs->make<TH1F>("reco_k_mu_open_angle_cosine","",50,-1,1);
  reco_k_mu_open_angle_cosine->GetXaxis()->SetTitle("Cosine of open angle");
  
  true_k_mu_open_angle_cosine = tfs->make<TH1F>("true_k_mu_open_angle_cosine","",50,-1,1);
  true_k_mu_open_angle_cosine->GetXaxis()->SetTitle("Cosine of open angle");
  
  reco_kaon_trkln_hist = tfs->make<TH1F>("reco_kaon_trkln_hist","",100,0,1000);
  reco_kaon_trkln_hist->GetXaxis()->SetTitle("Track length (cm)");
  
  true_kaon_trkln_hist = tfs->make<TH1F>("true_kaon_trkln_hist","",100,0,1000);
  true_kaon_trkln_hist->GetXaxis()->SetTitle("Track length (cm)");
  
  reco_muon_trkln_hist = tfs->make<TH1F>("reco_muon_trkln_hist","",100,0,1000);
  reco_muon_trkln_hist->GetXaxis()->SetTitle("Track length (cm)");
  
  true_muon_trkln_hist = tfs->make<TH1F>("true_muon_trkln_hist","",100,0,1000);
  true_muon_trkln_hist->GetXaxis()->SetTitle("Track length (cm)");
  
  reco_mu_mom_range = tfs->make<TH1F>("reco_mu_mom_range","",100,0,1000);
  reco_mu_mom_range->GetXaxis()->SetTitle("Transverse momentum (GeV/c)");
  
  true_mu_mom_range = tfs->make<TH1F>("true_mu_mom_range","",100,0,1000);
  true_mu_mom_range->GetXaxis()->SetTitle("Transverse momentum (GeV/c)");
  
  reco_muon_cont_plane_two_res_range_dEdX_hist = tfs->make<TH2F>("reco_muon_cont_plane_two_res_range_dEdX_hist","",200,0,100,200,0,50);
  reco_muon_cont_plane_two_res_range_dEdX_hist->GetXaxis()->SetTitle("Residual range R(cm)");
  reco_muon_cont_plane_two_res_range_dEdX_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  true_muon_cont_plane_best_res_range_dEdX_hist = tfs->make<TH2F>("true_muon_cont_plane_best_res_range_dEdX_hist","",200,0,100,200,0,50);
  true_muon_cont_plane_best_res_range_dEdX_hist->GetXaxis()->SetTitle("Residual range R(cm)");
  true_muon_cont_plane_best_res_range_dEdX_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  reco_k_pida = tfs->make<TH1F>("reco_k_pida","",50,0,100);
  reco_k_pida->GetXaxis()->SetTitle("PIDA");
  
  true_k_pida = tfs->make<TH1F>("true_k_pida","",50,0,100);
  true_k_pida->GetXaxis()->SetTitle("PIDA");
  
  reco_mu_pida = tfs->make<TH1F>("reco_mu_pida","",50,0,100);
  reco_mu_pida->GetXaxis()->SetTitle("PIDA");
  
  true_mu_pida = tfs->make<TH1F>("true_mu_pida","",50,0,100);
  true_mu_pida->GetXaxis()->SetTitle("PIDA");
  
  reco_k_dedx_ratio = tfs->make<TH1F>("reco_k_dedx_ratio","",50,0,50);
  reco_k_dedx_ratio->GetXaxis()->SetTitle("dE/dx ratio");
  
  true_k_dedx_ratio = tfs->make<TH1F>("true_k_dedx_ratio","",50,0,50);
  true_k_dedx_ratio->GetXaxis()->SetTitle("dE/dx ratio");
  
  reco_mu_dedx_ratio = tfs->make<TH1F>("reco_mu_dedx_ratio","",50,0,50);
  reco_mu_dedx_ratio->GetXaxis()->SetTitle("dE/dx ratio");
  
  true_mu_dedx_ratio = tfs->make<TH1F>("true_mu_dedx_ratio","",50,0,50);
  true_mu_dedx_ratio->GetXaxis()->SetTitle("dE/dx ratio");
  
  reco_de_dx_median_gradient = tfs->make<TH1F>("reco_de_dx_median_gradient","",500,-100,1000);
  reco_de_dx_median_gradient->GetXaxis()->SetTitle("Median energy gradient (MeV)");
  
  true_de_dx_median_gradient = tfs->make<TH1F>("true_de_dx_median_gradient","",500,-100,1000);
  true_de_dx_median_gradient->GetXaxis()->SetTitle("Median energy gradient (MeV)");
  
  reco_k_median_dedx_ratio = tfs->make<TH1F>("reco_k_median_dedx_ratio","",50,0,50);
  reco_k_median_dedx_ratio->GetXaxis()->SetTitle("Median dE/dx ratio");
  
  true_k_median_dedx_ratio = tfs->make<TH1F>("true_k_median_dedx_ratio","",50,0,50);
  true_k_median_dedx_ratio->GetXaxis()->SetTitle("Median dE/dx ratio");
  
  reco_mu_median_dedx_ratio = tfs->make<TH1F>("reco_mu_median_dedx_ratio","",50,0,50);
  reco_mu_median_dedx_ratio->GetXaxis()->SetTitle("Median dE/dx ratio");
  
  true_mu_median_dedx_ratio = tfs->make<TH1F>("true_mu_median_dedx_ratio","",50,0,50);
  true_mu_median_dedx_ratio->GetXaxis()->SetTitle("Median dE/dx ratio");
  
  reco_k_mean_first_dedx = tfs->make<TH1F>("reco_k_mean_first_dedx","",100,0,200);
  reco_k_mean_first_dedx->GetXaxis()->SetTitle("First dEdx average(MeV/cm)");
  
  true_k_mean_first_dedx = tfs->make<TH1F>("true_k_mean_first_dedx","",100,0,200);
  true_k_mean_first_dedx->GetXaxis()->SetTitle("First dE/dx average(MeV/cm)");
  
  reco_k_mean_last_dedx = tfs->make<TH1F>("reco_k_mean_last_dedx","",100,0,200);
  reco_k_mean_last_dedx->GetXaxis()->SetTitle("Last dE/dx average(MeV/cm)");
  
  true_k_mean_last_dedx = tfs->make<TH1F>("true_k_mean_last_dedx","",100,0,200);
  true_k_mean_last_dedx->GetXaxis()->SetTitle("Last dE/dx average(MeV/cm)");
  
  reco_mu_mean_first_dedx = tfs->make<TH1F>("reco_mu_mean_first_dedx","",100,0,200);
  reco_mu_mean_first_dedx->GetXaxis()->SetTitle("First dE/dx average(MeV/cm)");
  
  true_mu_mean_first_dedx = tfs->make<TH1F>("true_mu_mean_first_dedx","",100,0,200);
  true_mu_mean_first_dedx->GetXaxis()->SetTitle("First dEdx average(cm)");
  
  reco_mu_mean_last_dedx = tfs->make<TH1F>("reco_mu_mean_last_dedx","",100,0,200);
  reco_mu_mean_last_dedx->GetXaxis()->SetTitle("Last dE/dx average(MeV/cm)");
  
  true_mu_mean_last_dedx = tfs->make<TH1F>("true_mu_mean_last_dedx","",100,0,200);
  true_mu_mean_last_dedx->GetXaxis()->SetTitle("Last dE/dx average(MeV/cm)");
  
  reco_k_median_first_dedx = tfs->make<TH1F>("reco_k_median_first_dedx","",100,0,200);
  reco_k_median_first_dedx->GetXaxis()->SetTitle("First dE/dx median(MeV/cm)");
  
  true_k_median_first_dedx = tfs->make<TH1F>("true_k_median_first_dedx","",100,0,200);
  true_k_median_first_dedx->GetXaxis()->SetTitle("First dE/dx median(MeV/cm)");
  
  reco_k_median_last_dedx = tfs->make<TH1F>("reco_k_median_last_dedx","",100,0,200);
  reco_k_median_last_dedx->GetXaxis()->SetTitle("Last dE/dx median(MeV/cm)");
  
  true_k_median_last_dedx = tfs->make<TH1F>("true_k_median_last_dedx","",100,0,200);
  true_k_median_last_dedx->GetXaxis()->SetTitle("Last dE/dx median(MeV/cm)");
  
  reco_mu_median_first_dedx = tfs->make<TH1F>("reco_mu_median_first_dedx","",100,0,200);
  reco_mu_median_first_dedx->GetXaxis()->SetTitle("First dE/dx median(MeV/cm)");
  
  true_mu_median_first_dedx = tfs->make<TH1F>("true_mu_median_first_dedx","",100,0,200);
  true_mu_median_first_dedx->GetXaxis()->SetTitle("First dE/dx median(MeV/cm)");
  
  reco_mu_median_last_dedx = tfs->make<TH1F>("reco_mu_median_last_dedx","",100,0,200);
  reco_mu_median_last_dedx->GetXaxis()->SetTitle("Last dE/dx median(MeV/cm)");
  
  true_mu_median_last_dedx = tfs->make<TH1F>("true_mu_median_last_dedx","",100,0,200);
  true_mu_median_last_dedx->GetXaxis()->SetTitle("Last dE/dx median(MeV/cm)");
  
  reco_k_mean_len_ratio = tfs->make<TH1F>("reco_k_mean_len_ratio","",100,-20,100);
  reco_k_mean_len_ratio->GetXaxis()->SetTitle("Mean dE/dx/traklength");
  
  true_k_mean_len_ratio = tfs->make<TH1F>("true_k_mean_len_ratio","",100,-20,100);
  true_k_mean_len_ratio->GetXaxis()->SetTitle("Mean dE/dx/traklength");
  
  reco_mu_mean_len_ratio = tfs->make<TH1F>("reco_mu_mean_len_ratio","",100,-20,100);
  reco_mu_mean_len_ratio->GetXaxis()->SetTitle("Mean dE/dx/traklength");
  
  true_mu_mean_len_ratio = tfs->make<TH1F>("true_mu_mean_len_ratio","",100,-20,100);
  true_mu_mean_len_ratio->GetXaxis()->SetTitle("Mean dE/dx/traklength");
  
  reco_k_median_len_ratio = tfs->make<TH1F>("reco_k_median_len_ratio","",100,-20,100);
  reco_k_median_len_ratio->GetXaxis()->SetTitle("Median dE/dx/traklength");
  
  true_k_median_len_ratio = tfs->make<TH1F>("true_k_median_len_ratio","",100,-20,100);
  true_k_median_len_ratio->GetXaxis()->SetTitle("Median dE/dx/traklength");
  
  reco_mu_median_len_ratio = tfs->make<TH1F>("reco_mu_median_len_ratio","",100,-20,100);
  reco_mu_median_len_ratio->GetXaxis()->SetTitle("Median dE/dx/traklength");
  
  true_mu_median_len_ratio = tfs->make<TH1F>("true_mu_median_len_ratio","",100,-20,100);
  true_mu_median_len_ratio->GetXaxis()->SetTitle("Median dE/dx/traklength");
  
  reco_all_trk_median_dedx = tfs->make<TH1F>("reco_all_trk_median_dedx","",100,0,200);
  reco_all_trk_median_dedx->GetXaxis()->SetTitle("Median dE/dx");
  
  true_all_trk_median_dedx = tfs->make<TH1F>("true_all_trk_median_dedx","",100,0,200);
  true_all_trk_median_dedx->GetXaxis()->SetTitle("Median dE/dx");
  
  reco_all_trmu_median_dedx = tfs->make<TH1F>("reco_all_trmu_median_dedx","",100,0,200);
  reco_all_trmu_median_dedx->GetXaxis()->SetTitle("Median dE/dx");
  
  true_all_trmu_median_dedx = tfs->make<TH1F>("true_all_trmu_median_dedx","",100,0,200);
  true_all_trmu_median_dedx->GetXaxis()->SetTitle("Median dE/dx");
  
  k_best_plane_hist = tfs->make<TH1F>("k_best_plane_hist","",4,0,4);
  k_best_plane_hist->GetXaxis()->SetTitle("Kaon best plane");
  
  mu_best_plane_hist = tfs->make<TH1F>("mu_best_plane_hist","",4,0,4);
  mu_best_plane_hist->GetXaxis()->SetTitle("Muon best plane");
  
  /////////////// Test ////////////////
  
  event_k_pida = tfs->make<TH1F>("event_k_pida","",50,0,100);
  event_k_pida->GetXaxis()->SetTitle("PIDA");
  
  event_mu_pida = tfs->make<TH1F>("event_mu_pida","",50,0,100);
  event_mu_pida->GetXaxis()->SetTitle("PIDA");
  
  event_k_trklen = tfs->make<TH1F>("event_k_trklen","",100,0,1000);
  event_k_trklen->GetXaxis()->SetTitle("Track length (cm)");
  
  event_mu_trklen = tfs->make<TH1F>("event_mu_trklen","",100,0,1000);
  event_mu_trklen->GetXaxis()->SetTitle("Track length (cm)");
  
  event_de_dx_median_gradient = tfs->make<TH1F>("event_de_dx_median_gradient","",500,-100,1000);
  event_de_dx_median_gradient->GetXaxis()->SetTitle("Median energy gradient (MeV)");
  
  event_de_dx_mean_gradient = tfs->make<TH1F>("event_de_dx_mean_gradient","",500,-100,1000);
  event_de_dx_mean_gradient->GetXaxis()->SetTitle("Mean energy gradient (MeV)");
  
  event_k_mu_open_angle_cosine = tfs->make<TH1F>("event_k_mu_open_angle_cosine","",50,-1,1);
  event_k_mu_open_angle_cosine->GetXaxis()->SetTitle("Cosine of open angle");
  
  event_mu_mom_range = tfs->make<TH1F>("event_mu_mom_range","",100,0,1000);
  event_mu_mom_range->GetXaxis()->SetTitle("Transverse momentum (GeV/c)");
  
  event_kaon_cont_plane_best_res_range_dEdX_hist = tfs->make<TH2F>("event_kaon_cont_plane_best_res_range_dEdX_hist","",200,0,100,200,0,50);
  event_kaon_cont_plane_best_res_range_dEdX_hist->GetXaxis()->SetTitle("Residual range R(cm)");
  event_kaon_cont_plane_best_res_range_dEdX_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  ////////////// End of test /////////
  
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  fEventTree->Branch("event", &event);
  fEventTree->Branch("run", &run);
  fEventTree->Branch("subrun", &subrun);
  fEventTree->Branch("pri_mu_id",&pri_mu_id);
  fEventTree->Branch("pri_k_id",&pri_k_id);
  fEventTree->Branch("dec_mu_id",&dec_mu_id);
  fEventTree->Branch("reco_events",&reco_events,"reco_events/I");
  fEventTree->Branch("true_events",&true_events,"true_events/I");
  fEventTree->Branch("candidate_k",&candidate_k,"candidate_k/I");
  fEventTree->Branch("break_value",&break_value,"break_value/I");
}

//========================================================================
void Kplane2::endJob(){     

}
//========================================================================
void Kplane2::beginRun(const art::Run&){
  mf::LogInfo("Kplane2")<<"begin run..."<<std::endl;
}
//========================================================================

void Kplane2::analyze( const art::Event& evt){
  reset(); 
  
  bool isMC = !evt.isRealData();
        
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (isMC){
     if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
  } 
  
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (fSaveTrackInfo){
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
       art::fill_ptr_vector(tracklist, trackListHandle);
  }
  
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (fSaveClusterInfo){
      if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist,clusterListHandle);
  }
  
  art::ServiceHandle<cheat::BackTracker> bt;
  
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  std::vector<art::Ptr<recob::Vertex> > vertexlist;
  if (fSaveVertexInfo){
      if (evt.getByLabel(fVertexModuleLabel,vertexListHandle))
      art::fill_ptr_vector(vertexlist,vertexListHandle);
  }
  
  art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  art::FindManyP<recob::Hit> fmht(trackListHandle,evt,fTrackModuleLabel);
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();  
  
  bool isdata;
  
  if (evt.isRealData()){
    isdata = 1;
  }
  
  else isdata = 0;
  
  size_t n_reco_evt=0;
  size_t n_true_evt=0;
  size_t n_candidates=0;
  int break_indicator=0;
  
  if (!isdata){
     if (fSaveTrackInfo){
         size_t NTracks = tracklist.size(); 
	 trkf::TrackMomentumCalculator trkm;
	 if(NTracks >=3){
	    for(size_t i=0; i<NTracks;++i){
	         art::Ptr<recob::Track> ptrack(trackListHandle, i);
	         const recob::Track& track = *ptrack;
	         TVector3 pos, dir_start, dir_end, end;
		 double mup_tlen = 0;
		 pos = track.Vertex();
     	         dir_start = track.VertexDirection();
     	         dir_end   = track.EndDirection();
     	         end = track.End();
	         mup_tlen = track.Length();
		 int plane_2_primu=0;
	         int plane_1_primu=0;
	         int plane_0_primu=0;
		 if (mup_tlen >= 60){
		     if (fmcal.isValid()){
		         std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
			 for (size_t ical = 0; ical<calos.size(); ++ical){
			      if (!calos[ical]) continue;
			      if (!calos[ical]->PlaneID().isValid) continue;
			      int planenum = calos[ical]->PlaneID().Plane;
			      if (planenum<0||planenum>2) continue;
			      if (planenum == 2) plane_2_primu = int(calos[ical] -> dEdx().size());
			      if (planenum == 1) plane_1_primu = int(calos[ical] -> dEdx().size());
			      if (planenum == 0) plane_0_primu = int(calos[ical] -> dEdx().size());
			 }
		      }	      
			      if(plane_2_primu>=10 || plane_1_primu>=10 || plane_0_primu>=10){
			         float p_mu_st_x=pos.X();
			         float p_mu_st_y=pos.Y();
			         float p_mu_st_z=pos.Z();
			         float p_mu_en_x=end.X();
			         float p_mu_en_y=end.Y();
			         float p_mu_en_z=end.Z();
				 int prim_mu_index=int(i);
				 for(size_t i_1=0; i_1<NTracks;++i_1){
				     if(int(i_1)==prim_mu_index) continue;
				     art::Ptr<recob::Track> ptrack(trackListHandle, i_1);
				     const recob::Track& track = *ptrack;
				     TVector3 pos, dir_start, dir_end, end;
				     double k_tlen = 0;
				     pos = track.Vertex();
     	                             dir_start = track.VertexDirection();
     	                             dir_end   = track.EndDirection();
     	                             end = track.End();
	                             k_tlen = track.Length();
				     if (k_tlen >=10){
				         std::vector<const anab::Calorimetry*> calos = fmcal.at(i_1);
					 int k_plane_2=0;
					 int k_plane_1=0;
					 int k_plane_0=0;
					 for (size_t ical = 0; ical<calos.size(); ++ical){
					      if (!calos[ical]) continue;
			                      if (!calos[ical]->PlaneID().isValid) continue;
			                      int planenum = calos[ical]->PlaneID().Plane;
			                      if (planenum<0||planenum>2) continue;
			                      if (planenum == 2) k_plane_2 = int(calos[ical] -> dEdx().size());
					      if (planenum == 1) k_plane_1 = int(calos[ical] -> dEdx().size());
					      if (planenum == 0) k_plane_0 = int(calos[ical] -> dEdx().size());
					 }
			                      if(k_plane_2>=10 || k_plane_1>=10 || k_plane_0>=10){
					         if ((pos.X() > 5) && (pos.X() < 251.3) && (end.X() > 5) && (end.X() < 251.3) && (pos.Y() > -111.5) && (pos.Y() < 111.5) && (end.Y() > -111.5) && (end.Y() <111.5) && (pos.Z() > 5) && (pos.Z() < 1031.8) && (end.Z() > 5) && (end.Z() < 1031.8)){
						      float p_k_st_x=pos.X();
						      float p_k_st_y=pos.Y();
						      float p_k_st_z=pos.Z();
						      float p_k_en_x=end.X();
						      float p_k_en_y=end.Y();
						      float p_k_en_z=end.Z();
						      int prim_kaon_index=int(i_1);
						      double st_k_st_p_mu = TMath::Sqrt((p_k_st_x-p_mu_st_x)*(p_k_st_x-p_mu_st_x) + (p_k_st_y-p_mu_st_y)*(p_k_st_y-p_mu_st_y) + (p_k_st_z-p_mu_st_z)*(p_k_st_z-p_mu_st_z));
				                      double st_k_en_p_mu = TMath::Sqrt((p_k_st_x-p_mu_en_x)*(p_k_st_x-p_mu_en_x) + (p_k_st_y-p_mu_en_y)*(p_k_st_y-p_mu_en_y) + (p_k_st_z-p_mu_en_z)*(p_k_st_z-p_mu_en_z));
				                      double en_k_st_p_mu = TMath::Sqrt((p_k_en_x-p_mu_st_x)*(p_k_en_x-p_mu_st_x) + (p_k_en_y-p_mu_st_y)*(p_k_en_y-p_mu_st_y) + (p_k_en_z-p_mu_st_z)*(p_k_en_z-p_mu_st_z));
				                      double en_k_en_p_mu = TMath::Sqrt((p_k_en_x-p_mu_en_x)*(p_k_en_x-p_mu_en_x) + (p_k_en_y-p_mu_en_y)*(p_k_en_y-p_mu_en_y) + (p_k_en_z-p_mu_en_z)*(p_k_en_z-p_mu_en_z));
						      if(st_k_st_p_mu < 10 || st_k_en_p_mu < 10 || en_k_st_p_mu < 10 || en_k_en_p_mu < 10){
						         vector<double>dis_pri_mu_k_vec;
						         dis_pri_mu_k_vec.push_back(st_k_st_p_mu);
				                         dis_pri_mu_k_vec.push_back(st_k_en_p_mu);
				                         dis_pri_mu_k_vec.push_back(en_k_st_p_mu);
				                         dis_pri_mu_k_vec.push_back(en_k_en_p_mu);
							 double min_pri_mu_k_dis=1e10;
						         int pri_mu_pri_k_index=-1;
							 
							 for(int i_3=0; i_3<4; i_3++){
				                             if(dis_pri_mu_k_vec[i_3] < min_pri_mu_k_dis) {
						                min_pri_mu_k_dis=dis_pri_mu_k_vec[i_3];
						                pri_mu_pri_k_index=i_3;
						             }
				                         }
							 
							 for(size_t i_2=0; i_2<NTracks;++i_2){
							     if(int(i_2)==prim_mu_index || int(i_2)==prim_kaon_index) continue;
							        art::Ptr<recob::Track> ptrack(trackListHandle, i_2);
				                                const recob::Track& track = *ptrack;
				                                TVector3 pos, dir_start, dir_end, end;
				                                double mud_tlen = 0;
				                                pos = track.Vertex();
     	                                                        dir_start = track.VertexDirection();
     	                                                        dir_end   = track.EndDirection();
     	                                                        end = track.End();
	                                                        mud_tlen = track.Length();
								if (mud_tlen>50 && mud_tlen<60){ // mud_tlen>50 && mud_tlen<60
								    calos = fmcal.at(i_2);
								    int dmu_plane_2=0;
								    int dmu_plane_1=0;
								    int dmu_plane_0=0;
								    std::vector<const anab::Calorimetry*> calos = fmcal.at(i_2);
								    for (size_t ical = 0; ical<calos.size(); ++ical){
								         if (!calos[ical]) continue;
			                                                 if (!calos[ical]->PlaneID().isValid) continue;
			                                                 int planenum = calos[ical]->PlaneID().Plane;
			                                                 if (planenum<0||planenum>2) continue;
			                                                 if (planenum == 2) dmu_plane_2 = int(calos[ical] -> dEdx().size());
									 if (planenum == 1) dmu_plane_1 = int(calos[ical] -> dEdx().size());
									 if (planenum == 0) dmu_plane_0 = int(calos[ical] -> dEdx().size());
								    }	 
									 if(dmu_plane_2>= 10 || dmu_plane_1>=10 || dmu_plane_0>= 10){
									    if ((pos.X() > 5) && (pos.X() < 251.3) && (end.X() > 5) && (end.X() < 251.3) && (pos.Y() > -111.5) && (pos.Y() < 111.5) && (end.Y() > -111.5) && (end.Y() <111.5) && (pos.Z() > 5) && (pos.Z() < 1031.8) && (end.Z() > 5) && (end.Z() <1031.8)){
									         float d_mu_st_x=pos.X();
						                                 float d_mu_st_y=pos.Y();
						                                 float d_mu_st_z=pos.Z();
						                                 float d_mu_en_x=end.X();
						                                 float d_mu_en_y=end.Y();
						                                 float d_mu_en_z=end.Z();
										 int dec_mu_index=int(i_2);
										 float st_k_st_d_mu = TMath::Sqrt((p_k_st_x-d_mu_st_x)*(p_k_st_x-d_mu_st_x) + (p_k_st_y-d_mu_st_y)*(p_k_st_y-d_mu_st_y) + (p_k_st_z-d_mu_st_z)*(p_k_st_z-d_mu_st_z));  
					                                         float st_k_en_d_mu = TMath::Sqrt((p_k_st_x-d_mu_en_x)*(p_k_st_x-d_mu_en_x) + (p_k_st_y-d_mu_en_y)*(p_k_st_y-d_mu_en_y) + (p_k_st_z-d_mu_en_z)*(p_k_st_z-d_mu_en_z));
					                                         float en_k_st_d_mu = TMath::Sqrt((p_k_en_x-d_mu_st_x)*(p_k_en_x-d_mu_st_x) + (p_k_en_y-d_mu_st_y)*(p_k_en_y-d_mu_st_y) + (p_k_en_z-d_mu_st_z)*(p_k_en_z-d_mu_st_z));
					                                         float en_k_en_d_mu = TMath::Sqrt((p_k_en_x-d_mu_en_x)*(p_k_en_x-d_mu_en_x) + (p_k_en_y-d_mu_en_y)*(p_k_en_y-d_mu_en_y) + (p_k_en_z-d_mu_en_z)*(p_k_en_z-d_mu_en_z));
									         
										////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										 
										 if((pri_mu_pri_k_index==0) || (pri_mu_pri_k_index==1)){
										     if(en_k_st_d_mu < 10 || en_k_en_d_mu < 10){
										     
										        //////////////////////////////////////////// Acccessing calorimetry inormation for Kaons //////////////
											
											int k_best_planenum=-1;
											/*std::vector<const anab::Calorimetry*>*/ calos = fmcal.at(prim_kaon_index); // here
			                                                                int k_plane_0_hits=-1;int k_plane_1_hits=-1;int k_plane_2_hits=-1;
							                                for (size_t ical = 0; ical<calos.size(); ++ical){
							                                     if(!calos[ical]) continue;
			                                                                     if(!calos[ical]->PlaneID().isValid) continue;
			                                                                     int planenum = calos[ical]->PlaneID().Plane;
			                                                                     if(planenum<0||planenum>2) continue;
							                                     if(planenum == 0) k_plane_0_hits = int(calos[ical] -> dEdx().size());
			                                                                     if(planenum == 1) k_plane_1_hits = int(calos[ical] -> dEdx().size());
		                                                                             if(planenum == 2) k_plane_2_hits = int(calos[ical] -> dEdx().size());
							                                } 
						      
						                                        if(k_plane_0_hits != -1 || k_plane_1_hits != -1 || k_plane_2_hits != -1){
							                                if(k_plane_0_hits > k_plane_1_hits && k_plane_0_hits > k_plane_2_hits) k_best_planenum=0;
                                                                                        if(k_plane_1_hits > k_plane_0_hits && k_plane_1_hits > k_plane_2_hits) k_best_planenum=1;
                                                                                        if(k_plane_2_hits > k_plane_0_hits && k_plane_2_hits > k_plane_1_hits) k_best_planenum=2;
                                                                                        if(k_plane_2_hits==k_plane_0_hits && k_plane_2_hits > k_plane_1_hits) k_best_planenum=2;
                                                                                        if(k_plane_2_hits==k_plane_1_hits && k_plane_2_hits > k_plane_0_hits) k_best_planenum=2;
                                                                                        if(k_plane_1_hits==k_plane_0_hits && k_plane_1_hits > k_plane_2_hits) k_best_planenum=0;
                                                                                        if(k_plane_1_hits==k_plane_0_hits && k_plane_1_hits==k_plane_2_hits) k_best_planenum=2;
							                                }
				                                                         
											float dedxk_first_sum=0;
											float dedxk_last_sum=0;
											int firstk_nhits=0;
											int lastk_nhits=0;
											vector<float>dedx_firstk_vec;
											vector<float>dedx_lastk_vec; 
											vector<float>dedx_allk_vec;
											for (size_t ical = 0; ical<calos.size(); ++ical){
											     if(!calos[ical]) continue;
			                                                                     if(!calos[ical]->PlaneID().isValid) continue;
			                                                                     int planenum = calos[ical]->PlaneID().Plane;
			                                                                     if(planenum<0||planenum>2) continue;
											     if(planenum==k_best_planenum){
											         const size_t NHits = calos[ical] -> dEdx().size();
												 for(size_t iHit = 0; iHit < NHits; ++iHit){
												     if((calos[ical]->dEdx())[iHit]>100) continue;
												     dedx_allk_vec.push_back((calos[ical]->dEdx())[iHit]);
												     if((calos[ical]->ResidualRange())[iHit] < 3){
												         dedxk_first_sum=dedxk_first_sum+(calos[ical]->dEdx())[iHit];
													 firstk_nhits++;
													 dedx_firstk_vec.push_back((calos[ical]->dEdx())[iHit]);
												     }
												     if(calos[ical]->Range()-(calos[ical]->ResidualRange())[iHit] < 3){
												        dedxk_last_sum=dedxk_last_sum+(calos[ical]->dEdx())[iHit];
													lastk_nhits++;
													dedx_lastk_vec.push_back((calos[ical]->dEdx())[iHit]);
												     }
												  }
											      }
											 }
											 
											 float firstk_dedx_av=0;
											 float lastk_dedx_av=0;
											 if(firstk_nhits!=0) firstk_dedx_av=float(dedxk_first_sum)/firstk_nhits;
											 if(lastk_nhits!=0) lastk_dedx_av=float(dedxk_last_sum)/lastk_nhits;
											 std::sort(dedx_firstk_vec.begin(),dedx_firstk_vec.end());
	                                                                                 std::sort(dedx_lastk_vec.begin(),dedx_lastk_vec.end());
											 std::sort(dedx_allk_vec.begin(),dedx_allk_vec.end());
											 float med_firstk=0;
											 if(dedx_firstk_vec.size()){
	                                                                                    if(dedx_firstk_vec.size()%2==1){
	                                                                                       med_firstk=dedx_firstk_vec[(dedx_firstk_vec.size()-1)/2];
	                                                                                    }
	                                                                                    if(dedx_firstk_vec.size()%2==0){
	                                                                                       med_firstk=float(dedx_firstk_vec[dedx_firstk_vec.size()/2] + dedx_firstk_vec[(dedx_firstk_vec.size()/2)-1])/2;
	                                                                                    }
	                                                                                  }
											  float med_lastk=0;
											  if(dedx_lastk_vec.size()){
	                                                                                    if(dedx_lastk_vec.size()%2==1){
	                                                                                       med_lastk=dedx_lastk_vec[(dedx_lastk_vec.size()-1)/2];
	                                                                                    }
	                                                                                    if(dedx_lastk_vec.size()%2==0){
	                                                                                       med_lastk=float(dedx_lastk_vec[dedx_lastk_vec.size()/2] + dedx_lastk_vec[(dedx_lastk_vec.size()/2)-1])/2;
	                                                                                    }
	                                                                                  }
											  float med_allk=0;
											  if(dedx_allk_vec.size()){
											     if(dedx_allk_vec.size()%2==1){
											        med_allk=dedx_allk_vec[(dedx_allk_vec.size()-1)/2];
											     }
											     if(dedx_allk_vec.size()%2==0){
											        med_allk=float(dedx_allk_vec[dedx_allk_vec.size()/2] + dedx_allk_vec[(dedx_allk_vec.size()/2)-1])/2;
											     }
											  }
											 /////////////// If you want min max inforamtions insert here /////////////
											 float k_pida=0;
											 float k_pida_sum=0;
											 int k_used_points=0;
											 for (size_t ical = 0; ical<calos.size(); ++ical){
											      if(!calos[ical]) continue;
			                                                                      if(!calos[ical]->PlaneID().isValid) continue;
			                                                                      int planenum = calos[ical]->PlaneID().Plane;
			                                                                      if(planenum<0||planenum>2) continue;
											      if(planenum==k_best_planenum){
											         const size_t NHits = calos[ical] -> dEdx().size();
												 for(size_t iHit = 0; iHit < NHits; ++iHit){
												     if((calos[ical]->dEdx())[iHit]>100) continue;
												     if((calos[ical]->ResidualRange())[iHit] < 30){
												        float res=(calos[ical]->ResidualRange())[iHit];
													k_pida_sum=k_pida_sum+((calos[ical]->dEdx())[iHit])*(TMath::Power(res,0.42));
													k_used_points++;
												     }
												  }
											      }
											 }
											 
											 if(k_used_points!=0)k_pida=float(k_pida_sum)/k_used_points;
											 
											///////////////////////////////////////// End of calorimetry information for kaons /////////////////////////////
											
											//////////////////////////////////////// Accessing calorimetry info of Muons ////////////////////////////////
											
											int mu_best_planenum=-1;
											calos = fmcal.at(dec_mu_index);
											int mu_plane_0_hits=-1;int mu_plane_1_hits=-1;int mu_plane_2_hits=-1;
											for(size_t ical = 0; ical<calos.size(); ++ical){
											    if(!calos[ical]) continue;
			                                                                    if(!calos[ical]->PlaneID().isValid) continue;
			                                                                    int planenum = calos[ical]->PlaneID().Plane;
			                                                                    if(planenum<0||planenum>2) continue;
							                                    if(planenum == 0) mu_plane_0_hits = int(calos[ical] -> dEdx().size());
			                                                                    if(planenum == 1) mu_plane_1_hits = int(calos[ical] -> dEdx().size());
		                                                                            if(planenum == 2) mu_plane_2_hits = int(calos[ical] -> dEdx().size());
											}
											
											if(mu_plane_0_hits != -1 || mu_plane_1_hits != -1 || mu_plane_2_hits != -1){
							                                if(mu_plane_0_hits > mu_plane_1_hits && mu_plane_0_hits > mu_plane_2_hits) mu_best_planenum=0;
                                                                                        if(mu_plane_1_hits > mu_plane_0_hits && mu_plane_1_hits > mu_plane_2_hits) mu_best_planenum=1;
                                                                                        if(mu_plane_2_hits > mu_plane_0_hits && mu_plane_2_hits > mu_plane_1_hits) mu_best_planenum=2;
                                                                                        if(mu_plane_2_hits==mu_plane_0_hits && mu_plane_2_hits > mu_plane_1_hits) mu_best_planenum=2;
                                                                                        if(mu_plane_2_hits==mu_plane_1_hits && mu_plane_2_hits > mu_plane_0_hits) mu_best_planenum=2;
                                                                                        if(mu_plane_1_hits==mu_plane_0_hits && mu_plane_1_hits > mu_plane_2_hits) mu_best_planenum=0;
                                                                                        if(mu_plane_1_hits==mu_plane_0_hits && mu_plane_1_hits==mu_plane_2_hits) mu_best_planenum=2;
							                                }
											
											float dedxmu_first_sum=0;
											float dedxmu_last_sum=0;
											int firstmu_nhits=0;
											int lastmu_nhits=0;
											vector<float>dedx_firstmu_vec;
											vector<float>dedx_lastmu_vec; 
											vector<float>dedx_allmu_vec;
											for (size_t ical = 0; ical<calos.size(); ++ical){
											     if(!calos[ical]) continue;
			                                                                     if(!calos[ical]->PlaneID().isValid) continue;
			                                                                     int planenum = calos[ical]->PlaneID().Plane;
			                                                                     if(planenum<0||planenum>2) continue;
											     if(planenum==mu_best_planenum){
											         const size_t NHits = calos[ical] -> dEdx().size();
												 for(size_t iHit = 0; iHit < NHits; ++iHit){
												     if((calos[ical]->dEdx())[iHit]>100) continue;
												     dedx_allmu_vec.push_back((calos[ical]->dEdx())[iHit]);
												     if((calos[ical]->ResidualRange())[iHit] < 3){
												         dedxmu_first_sum=dedxmu_first_sum+(calos[ical]->dEdx())[iHit];
													 firstmu_nhits++;
													 dedx_firstmu_vec.push_back((calos[ical]->dEdx())[iHit]);
												     }
												     if(calos[ical]->Range()-(calos[ical]->ResidualRange())[iHit] < 3){
												        dedxmu_last_sum=dedxmu_last_sum+(calos[ical]->dEdx())[iHit];
													lastmu_nhits++;
													dedx_lastmu_vec.push_back((calos[ical]->dEdx())[iHit]);
												     }
												  }
											      }
											 }
											
											 float firstmu_dedx_av=0;
											 float lastmu_dedx_av=0;
											 if(firstmu_nhits!=0) firstmu_dedx_av=float(dedxmu_first_sum)/firstmu_nhits;
											 if(lastmu_nhits!=0) lastmu_dedx_av=float(dedxmu_last_sum)/lastmu_nhits;
											 std::sort(dedx_firstmu_vec.begin(),dedx_firstmu_vec.end());
	                                                                                 std::sort(dedx_lastmu_vec.begin(),dedx_lastmu_vec.end());
											 std::sort(dedx_allmu_vec.begin(),dedx_allmu_vec.end());
											 float med_firstmu=0;
											 if(dedx_firstmu_vec.size()){
	                                                                                    if(dedx_firstmu_vec.size()%2==1){
	                                                                                       med_firstmu=dedx_firstmu_vec[(dedx_firstmu_vec.size()-1)/2];
	                                                                                    }
	                                                                                    if(dedx_firstmu_vec.size()%2==0){
	                                                                                       med_firstmu=float(dedx_firstmu_vec[dedx_firstmu_vec.size()/2] + dedx_firstmu_vec[(dedx_firstmu_vec.size()/2)-1])/2;
	                                                                                    }
	                                                                                  }
											  float med_lastmu=0;
											  if(dedx_lastmu_vec.size()){
	                                                                                    if(dedx_lastmu_vec.size()%2==1){
	                                                                                       med_lastmu=dedx_lastmu_vec[(dedx_lastmu_vec.size()-1)/2];
	                                                                                    }
	                                                                                    if(dedx_lastmu_vec.size()%2==0){
	                                                                                       med_lastmu=float(dedx_lastmu_vec[dedx_lastmu_vec.size()/2] + dedx_lastmu_vec[(dedx_lastmu_vec.size()/2)-1])/2;
	                                                                                    }
	                                                                                  }
											  float med_allmu=0;
											  if(dedx_allmu_vec.size()){
											     if(dedx_allmu_vec.size()%2==1){
											        med_allmu=dedx_allmu_vec[(dedx_allmu_vec.size()-1)/2];
											     }
											     if(dedx_allmu_vec.size()%2==0){
											        med_allmu=float(dedx_allmu_vec[dedx_allmu_vec.size()/2] + dedx_allmu_vec[(dedx_allmu_vec.size()/2)-1])/2;
											     }
											  }
											 /////////////// If you want min max inforamtions insert here /////////////
											 float mu_pida=0;
											 float mu_pida_sum=0;
											 int mu_used_points=0;
											 for (size_t ical = 0; ical<calos.size(); ++ical){
											      if(!calos[ical]) continue;
			                                                                      if(!calos[ical]->PlaneID().isValid) continue;
			                                                                      int planenum = calos[ical]->PlaneID().Plane;
			                                                                      if(planenum<0||planenum>2) continue;
											      if(planenum==mu_best_planenum){
											         const size_t NHits = calos[ical] -> dEdx().size();
												 for(size_t iHit = 0; iHit < NHits; ++iHit){
												     if((calos[ical]->dEdx())[iHit]>100) continue;
												     if((calos[ical]->ResidualRange())[iHit] < 30){
												        float res=(calos[ical]->ResidualRange())[iHit];
													mu_pida_sum=mu_pida_sum+((calos[ical]->dEdx())[iHit])*(TMath::Power(res,0.42));
													mu_used_points++;
												     }
												  }
											      }
											 }
											 
											 if(mu_used_points!=0)mu_pida=float(mu_pida_sum)/mu_used_points;
											
											///////////////////////////////////// End of calorimetry info for muons ////////////////////////////////////
											
											//////////////////////////////////// defining k mu angle ///////////////////////////
											
											float average_commen_point_x=0;float average_commen_point_y=0;float average_commen_point_z=0;
									                float av_common_k_cos_alpha=0;float av_common_k_cos_beta=0;float av_common_k_cos_gamma=0;
									                float av_common_mu_cos_alpha=0;float av_common_mu_cos_beta=0;float av_common_mu_cos_gamma=0;
											
											if(en_k_st_d_mu < en_k_en_d_mu){
											   average_commen_point_x = float(p_k_en_x+d_mu_st_x)/2;
										           average_commen_point_y = float(p_k_en_y+d_mu_st_y)/2;
										           average_commen_point_z = float(p_k_en_z+d_mu_st_z)/2;
											   
											   float k_av_dis=TMath::Sqrt((average_commen_point_x-p_k_st_x)*(average_commen_point_x-p_k_st_x) + (average_commen_point_y-p_k_st_y)*(average_commen_point_y-p_k_st_y) +(average_commen_point_z-p_k_st_z)*(average_commen_point_z-p_k_st_z));
										           float mu_av_dis=TMath::Sqrt((d_mu_en_x-average_commen_point_x)*(d_mu_en_x-average_commen_point_x) + (d_mu_en_y-average_commen_point_y)*(d_mu_en_y-average_commen_point_y) + (d_mu_en_z-average_commen_point_z)*(d_mu_en_z-average_commen_point_z));
										
										           av_common_k_cos_alpha  = float(average_commen_point_x-p_k_st_x)/k_av_dis; 
										           av_common_k_cos_beta   = float(average_commen_point_y-p_k_st_y)/k_av_dis; 
										           av_common_k_cos_gamma  = float(average_commen_point_z-p_k_st_z)/k_av_dis; 
										           av_common_mu_cos_alpha = float(d_mu_en_x-average_commen_point_x)/mu_av_dis; 
										           av_common_mu_cos_beta  = float(d_mu_en_y-average_commen_point_y)/mu_av_dis; 
										           av_common_mu_cos_gamma = float(d_mu_en_z-average_commen_point_z)/mu_av_dis; 
											}
											
											if(en_k_st_d_mu > en_k_en_d_mu){
											   average_commen_point_x = float(p_k_en_x+d_mu_en_x)/2;
										           average_commen_point_y = float(p_k_en_y+d_mu_en_y)/2;
										           average_commen_point_z = float(p_k_en_z+d_mu_en_z)/2;
											   
											   float k_av_dis=TMath::Sqrt((average_commen_point_x-p_k_st_x)*(average_commen_point_x-p_k_st_x) + (average_commen_point_y-p_k_st_y)*(average_commen_point_y-p_k_st_y) +(average_commen_point_z-p_k_st_z)*(average_commen_point_z-p_k_st_z));
										           float mu_av_dis=TMath::Sqrt((d_mu_st_x-average_commen_point_x)*(d_mu_st_x-average_commen_point_x) + (d_mu_st_y-average_commen_point_y)*(d_mu_st_y-average_commen_point_y) + (d_mu_st_z-average_commen_point_z)*(d_mu_st_z-average_commen_point_z));
										
										           av_common_k_cos_alpha  = float(average_commen_point_x-p_k_st_x)/k_av_dis; 
										           av_common_k_cos_beta   = float(average_commen_point_y-p_k_st_y)/k_av_dis; 
										           av_common_k_cos_gamma  = float(average_commen_point_z-p_k_st_z)/k_av_dis; 
										           av_common_mu_cos_alpha = float(d_mu_st_x-average_commen_point_x)/mu_av_dis; 
										           av_common_mu_cos_beta  = float(d_mu_st_y-average_commen_point_y)/mu_av_dis; 
										           av_common_mu_cos_gamma = float(d_mu_st_z-average_commen_point_z)/mu_av_dis; 
											}
											
											float av_common_open_angle = TMath::ACos(av_common_k_cos_alpha*av_common_mu_cos_alpha + av_common_k_cos_beta*av_common_mu_cos_beta +av_common_k_cos_gamma*av_common_mu_cos_gamma);
											n_reco_evt++;
											
											if(med_firstk-med_lastmu>=5){ ///////////////////// K/Mu Energy gradient cut here *****************
											
											  if(TMath::Cos(av_common_open_angle) <= 0.9 && TMath::Cos(av_common_open_angle) >=-0.9){ ///////////// K/Mu open angle cut *******************
											     float mu_mom= (trkm.GetTrackMomentum(mud_tlen,13))*1000;
											     if(k_pida > 10 && k_pida < 20){ // K PIDA cut *********** // k_pida > 10 && k_pida < 20
												//if(med_allk >=2 && med_allk <=10){ // K all med dE/dx
												if(mu_pida < 10){ // Mu PIDA cut *********
												//if(med_allmu >=2 && med_allmu <=4){ // Mu all med dE/dx
												  //if(mu_mom*TMath::Sin(av_common_open_angle) < 200 && mu_mom*TMath::Sin(av_common_open_angle)> 150){ //// Mu mom cut ********
											              if(break_indicator==0){
											                  pri_mu_id = prim_mu_index;
													  pri_k_id = prim_kaon_index;
													  dec_mu_id = dec_mu_index;
													  break_indicator++;
											                  /*std::vector<const anab::Calorimetry*>*/ calos = fmcal.at(prim_kaon_index); /// *** here
											                  for (size_t ical = 0; ical<calos.size(); ++ical){
											                       if (!calos[ical]) continue;
			                                                                                       if (!calos[ical]->PlaneID().isValid) continue;
			                                                                                       int planenum = calos[ical]->PlaneID().Plane;
			                                                                                       if (planenum<0||planenum>2) continue;
											                       if (planenum==k_best_planenum){
											                       const size_t NHits = calos[ical] -> dEdx().size();
											                      for(size_t iHit = 0; iHit < NHits; ++iHit){
											                          reco_kaon_cont_plane_two_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												              }
											                    }
											                 }
											                 reco_de_dx_gradient->Fill(firstk_dedx_av-lastmu_dedx_av);
											                 reco_k_mu_open_angle_cosine->Fill(TMath::Cos(av_common_open_angle));
											                 reco_kaon_trkln_hist->Fill(k_tlen);
											                 reco_muon_trkln_hist->Fill(mud_tlen);
											                 reco_mu_mom_range->Fill(mu_mom*TMath::Sin(av_common_open_angle));
													 if(lastk_dedx_av!=0)reco_k_dedx_ratio->Fill(float(firstk_dedx_av)/lastk_dedx_av);
													 if(lastmu_dedx_av!=0)reco_mu_dedx_ratio->Fill(float(firstmu_dedx_av)/lastmu_dedx_av);
													 reco_de_dx_median_gradient->Fill(med_firstk-med_lastmu);
													 if(med_lastk!=0)reco_k_median_dedx_ratio->Fill(float(med_firstk)/med_lastk);
													 if(med_lastmu!=0)reco_mu_median_dedx_ratio->Fill(float(med_firstmu)/med_lastmu);
											                 calos = fmcal.at(dec_mu_index);
											                 for (size_t ical = 0; ical<calos.size(); ++ical){
											                      if (!calos[ical]) continue;
			                                                                                      if (!calos[ical]->PlaneID().isValid) continue;
			                                                                                      int planenum = calos[ical]->PlaneID().Plane;
			                                                                                      if (planenum<0||planenum>2) continue;
											                      if (planenum==mu_best_planenum){
											                          const size_t NHits = calos[ical] -> dEdx().size();
											                          for(size_t iHit = 0; iHit < NHits; ++iHit){
											                              reco_muon_cont_plane_two_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                  }
											                       }
											                   }
											                   reco_k_pida->Fill(k_pida);
											                   reco_mu_pida->Fill(mu_pida);
													   reco_k_mean_first_dedx->Fill(firstk_dedx_av);
													   reco_k_mean_last_dedx->Fill(lastk_dedx_av);
													   reco_mu_mean_first_dedx->Fill(firstmu_dedx_av);
													   reco_mu_mean_last_dedx->Fill(lastmu_dedx_av);
													   reco_k_median_first_dedx->Fill(med_firstk);
													   reco_k_median_last_dedx->Fill(med_lastk);
													   reco_mu_median_first_dedx->Fill(med_firstmu);
													   reco_mu_median_last_dedx->Fill(med_lastmu);
													   reco_k_mean_len_ratio->Fill((float(firstk_dedx_av-lastk_dedx_av))/k_tlen);
													   reco_mu_mean_len_ratio->Fill((float(firstmu_dedx_av-lastmu_dedx_av))/mud_tlen);
													   reco_k_median_len_ratio->Fill((float(med_firstk-med_lastk))/k_tlen);
													   reco_mu_median_len_ratio->Fill((float(med_firstmu-med_lastmu))/mud_tlen);
													   reco_all_trk_median_dedx->Fill(med_allk);
													   reco_all_trmu_median_dedx->Fill(med_allmu);
													   k_best_plane_hist->Fill(k_best_planenum);
													   mu_best_plane_hist->Fill(mu_best_planenum);
											////////////////////////////////// End of k mu angle //////////////////////////////
											
											const sim::ParticleList& plist = bt->ParticleList();
                                                                                        sim::ParticleList::const_iterator itPart = plist.begin(),pend = plist.end();
                                                                                        std::string pri("primary");
										        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
											    const simb::MCParticle* pPart = (itPart++)->second;
	                                                                                    bool isPrimary = pPart->Process()== pri;
											    if (isPrimary && pPart->PdgCode()==321){
											        if (pPart->EndE()*1000 <= 494.){
												    int trk_id  = pPart->TrackId();;
				                                                                    float end_x = pPart->EndPosition()[0];
				                                                                    float end_y = pPart->EndPosition()[1];
				                                                                    float end_z = pPart->EndPosition()[2];
												    sim::ParticleList::const_iterator jtPart = plist.begin(),send = plist.end();
											            for(size_t jPart = 0; (jPart < plist.size()) && (jtPart != send); ++jPart){
												    const simb::MCParticle* sPart = (jtPart++)->second;
												    if ((sPart->Process()=="Decay") && (sPart->Mother()==trk_id) &&(std::abs(sPart->Vx()-end_x)<0.01) && (std::abs(sPart->Vy()-end_y)<0.01) && (std::abs(sPart->Vz()-end_z)<0.01)){
												         if (sPart->PdgCode() == -13 || sPart->PdgCode() == 211){
													      
													      ///////////////////// Test ///////////////
													      
													      event_k_pida->Fill(k_pida);
													      event_mu_pida->Fill(mu_pida);
													      event_k_trklen->Fill(k_tlen);
													      event_mu_trklen->Fill(mud_tlen);
													      event_de_dx_median_gradient->Fill(med_firstk-med_lastmu);
													      event_de_dx_mean_gradient->Fill(firstk_dedx_av-lastmu_dedx_av);
													      event_k_mu_open_angle_cosine->Fill(TMath::Cos(av_common_open_angle));
													      event_mu_mom_range->Fill(mu_mom*TMath::Sin(av_common_open_angle));
													      std::vector<const anab::Calorimetry*> calos = fmcal.at(prim_kaon_index);
											                      for(size_t ical = 0; ical<calos.size(); ++ical){
											                          if(!calos[ical]) continue;
			                                                                                          if(!calos[ical]->PlaneID().isValid) continue;
			                                                                                          int planenum = calos[ical]->PlaneID().Plane;
			                                                                                          if(planenum<0||planenum>2) continue;
											                          if(planenum==k_best_planenum){
											                              const size_t NHits = calos[ical] -> dEdx().size();
											                              for(size_t iHit = 0; iHit < NHits; ++iHit){
													                  event_kaon_cont_plane_best_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                      }
												                   }
											                       }
													      
													      //////////////////// End test ////////////
													      
													      n_candidates++;
													      std::vector< art::Ptr<recob::Hit> > allKHits = fmht.at(prim_kaon_index);
													      std::vector< art::Ptr<recob::Hit> > allMuHits = fmht.at(dec_mu_index);
													      std::map<int,double> trk_k_ide;
													      for(size_t h = 0; h < allKHits.size(); ++h){
													          art::Ptr<recob::Hit> hit = allKHits[h];
													          std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
														  for(size_t e = 0; e < TrackIDs.size(); ++e){
														      trk_k_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
													          }
													       }
													       double maxke = -1;
													       double totke = 0;
													       int Track_k_id = 0;
													       for(std::map<int,double>::iterator ii=trk_k_ide.begin();ii!=trk_k_ide.end(); ++ii){
													           totke += ii->second;
														   if((ii->second)>maxke){
														        maxke = ii->second;
													                Track_k_id=ii->first;
														   }
														}
													        std::map<int,double> trk_mu_ide;
														for(size_t h = 0; h < allMuHits.size(); ++h){
														    art::Ptr<recob::Hit> hit = allMuHits[h];
														    std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
														    for(size_t e = 0; e < TrackIDs.size(); ++e){
															trk_mu_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
														    }
													         }
													         double maxmue = -1;
														 double totmue = 0;
														 int Track_mu_id = 0;
														 for(std::map<int,double>::iterator ii=trk_mu_ide.begin();ii!=trk_mu_ide.end(); ++ii){
														     totmue += ii->second;
														     if((ii->second)>maxmue){
															 maxmue = ii->second;
														         Track_mu_id=ii->first;
														     }
														  }
													          if(Track_k_id==int(pPart->TrackId()) && Track_mu_id==int(sPart->TrackId())){
														     std::cout << "K id : " << Track_k_id << "  " << int(pPart->TrackId()) << std::endl;
														     std::cout << "Mu id : " << Track_mu_id << "  " << int(sPart->TrackId()) << std::endl;
														     
														     
														     n_true_evt++;
														     std::vector<const anab::Calorimetry*> calos = fmcal.at(prim_kaon_index);
											                             for(size_t ical = 0; ical<calos.size(); ++ical){
											                                 if(!calos[ical]) continue;
			                                                                                                 if(!calos[ical]->PlaneID().isValid) continue;
			                                                                                                 int planenum = calos[ical]->PlaneID().Plane;
			                                                                                                 if(planenum<0||planenum>2) continue;
											                                 if(planenum==k_best_planenum){
											                                    const size_t NHits = calos[ical] -> dEdx().size();
											                                    for(size_t iHit = 0; iHit < NHits; ++iHit){
													                        true_kaon_cont_plane_best_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                            }
												                         }
											                             }
														     true_de_dx_gradient->Fill(firstk_dedx_av-lastmu_dedx_av);
														     true_k_mu_open_angle_cosine->Fill(TMath::Cos(av_common_open_angle));
												                     true_kaon_trkln_hist->Fill(k_tlen);
														     true_muon_trkln_hist->Fill(mud_tlen);
														     float mu_mom= (trkm.GetTrackMomentum(mud_tlen,13))*1000;
											                             true_mu_mom_range->Fill(mu_mom*TMath::Sin(av_common_open_angle));
														     if(lastk_dedx_av!=0)true_k_dedx_ratio->Fill(float(firstk_dedx_av)/lastk_dedx_av);
														     if(lastmu_dedx_av!=0)true_mu_dedx_ratio->Fill(float(firstmu_dedx_av)/lastmu_dedx_av);
														     true_de_dx_median_gradient->Fill(med_firstk-med_lastmu);
														     if(med_lastk!=0)true_k_median_dedx_ratio->Fill(float(med_firstk)/med_lastk);
														     if(med_lastmu!=0)true_mu_median_dedx_ratio->Fill(float(med_firstmu)/med_lastmu);
														     calos = fmcal.at(dec_mu_index);
											                             for(size_t ical = 0; ical<calos.size(); ++ical){
											                                 if(!calos[ical]) continue;
			                                                                                                 if(!calos[ical]->PlaneID().isValid) continue;
			                                                                                                 int planenum = calos[ical]->PlaneID().Plane;
			                                                                                                 if(planenum<0||planenum>2) continue;
											                                 if(planenum==mu_best_planenum){
											                                    const size_t NHits = calos[ical] -> dEdx().size();
											                                    for(size_t iHit = 0; iHit < NHits; ++iHit){
													                        true_muon_cont_plane_best_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                            }
												                         }
											                             }
														     true_k_pida->Fill(k_pida);
														     true_mu_pida->Fill(mu_pida);
														     true_k_mean_first_dedx->Fill(firstk_dedx_av);
														     true_k_mean_last_dedx->Fill(lastk_dedx_av);
														     true_mu_mean_first_dedx->Fill(firstmu_dedx_av);
													             true_mu_mean_last_dedx->Fill(lastmu_dedx_av);
														     true_k_median_first_dedx->Fill(med_firstk);
													             true_k_median_last_dedx->Fill(med_lastk);
														     true_mu_median_first_dedx->Fill(med_firstmu);
													             true_mu_median_last_dedx->Fill(med_lastmu);
														     true_k_mean_len_ratio->Fill((float(firstk_dedx_av-lastk_dedx_av))/k_tlen);
														     true_mu_mean_len_ratio->Fill((float(firstmu_dedx_av-lastmu_dedx_av))/mud_tlen);
														     true_k_median_len_ratio->Fill((float(med_firstk-med_lastk))/k_tlen);
														     true_mu_median_len_ratio->Fill((float(med_firstmu-med_lastmu))/mud_tlen);
														     true_all_trk_median_dedx->Fill(med_allk);
														     true_all_trmu_median_dedx->Fill(med_allmu);
													          } // trk g4 match
													       } // getting true mu
													    } // getting scttering particles
													 } // 2 nd loop over geant list for P
												      } // getting stopping K
											           } // getting primary K...
											         } // 1 st loop over geant list for P...
										               } // break indicator....
											     //} // Mu mom cut .....
											    //} // mu all dE/dx cut
											   } // mu PIDA cut ...
											  //} // K all dE/dx cut...
											 } // K PIDA cut....
										     //} // K and Mu min max hit cut....
										   //} // Muon track length cut .....
									          } // K/Mu open angle cut....
									       } // K/Mu energy gradient cut ......
									    } // decay mu connected to end of K....
									} // primary mu is connected to start of K....
									    
                                                                       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
										 
								      if((pri_mu_pri_k_index==2) || (pri_mu_pri_k_index==3)){
								          if(st_k_st_d_mu < 10 || st_k_en_d_mu < 10){
								             //////////////////////////////////////////// Acccessing calorimetry inormation for Kaons //////////////
									     int k_best_planenum=-1;
									     /*std::vector<const anab::Calorimetry*>*/ calos = fmcal.at(prim_kaon_index);
			                                                     int k_plane_0_hits=-1;int k_plane_1_hits=-1;int k_plane_2_hits=-1;
							                     for (size_t ical = 0; ical<calos.size(); ++ical){
							                          if(!calos[ical]) continue;
			                                                          if(!calos[ical]->PlaneID().isValid) continue;
			                                                          int planenum = calos[ical]->PlaneID().Plane;
			                                                          if(planenum<0||planenum>2) continue;
							                          if(planenum == 0) k_plane_0_hits = int(calos[ical] -> dEdx().size());
			                                                          if(planenum == 1) k_plane_1_hits = int(calos[ical] -> dEdx().size());
		                                                                  if(planenum == 2) k_plane_2_hits = int(calos[ical] -> dEdx().size());
							                      } 
						      
						                              if(k_plane_0_hits != -1 || k_plane_1_hits != -1 || k_plane_2_hits != -1){
							                      if(k_plane_0_hits > k_plane_1_hits && k_plane_0_hits > k_plane_2_hits) k_best_planenum=0;
                                                                              if(k_plane_1_hits > k_plane_0_hits && k_plane_1_hits > k_plane_2_hits) k_best_planenum=1;
                                                                              if(k_plane_2_hits > k_plane_0_hits && k_plane_2_hits > k_plane_1_hits) k_best_planenum=2;
                                                                              if(k_plane_2_hits==k_plane_0_hits && k_plane_2_hits > k_plane_1_hits) k_best_planenum=2;
                                                                              if(k_plane_2_hits==k_plane_1_hits && k_plane_2_hits > k_plane_0_hits) k_best_planenum=2;
                                                                              if(k_plane_1_hits==k_plane_0_hits && k_plane_1_hits > k_plane_2_hits) k_best_planenum=0;
                                                                              if(k_plane_1_hits==k_plane_0_hits && k_plane_1_hits==k_plane_2_hits) k_best_planenum=2;
							                      }
				                                                         
									      float dedxk_first_sum=0;
								              float dedxk_last_sum=0;
									      int firstk_nhits=0;
									      int lastk_nhits=0;
									      vector<float>dedx_firstk_vec;
									      vector<float>dedx_lastk_vec;
									      vector<float>dedx_allk_vec; 
									      for (size_t ical = 0; ical<calos.size(); ++ical){
									           if(!calos[ical]) continue;
			                                                           if(!calos[ical]->PlaneID().isValid) continue;
			                                                           int planenum = calos[ical]->PlaneID().Plane;
			                                                           if(planenum<0||planenum>2) continue;
									           if(planenum==k_best_planenum){
										      const size_t NHits = calos[ical] -> dEdx().size();
									              for(size_t iHit = 0; iHit < NHits; ++iHit){
										          if((calos[ical]->dEdx())[iHit]>100) continue;
											  dedx_allk_vec.push_back((calos[ical]->dEdx())[iHit]);    
											  if((calos[ical]->ResidualRange())[iHit] < 3){
											      dedxk_first_sum=dedxk_first_sum+(calos[ical]->dEdx())[iHit];
											      firstk_nhits++;
											      dedx_firstk_vec.push_back((calos[ical]->dEdx())[iHit]);
											   }
											   if(calos[ical]->Range()-(calos[ical]->ResidualRange())[iHit] < 3){
											      dedxk_last_sum=dedxk_last_sum+(calos[ical]->dEdx())[iHit];
											      lastk_nhits++;
											      dedx_lastk_vec.push_back((calos[ical]->dEdx())[iHit]);
											   }
											}
										     }
									          }
											 
										       float firstk_dedx_av=0;
										       float lastk_dedx_av=0;
										       if(firstk_nhits!=0) firstk_dedx_av=float(dedxk_first_sum)/firstk_nhits;
										       if(lastk_nhits!=0) lastk_dedx_av=float(dedxk_last_sum)/lastk_nhits;
										       std::sort(dedx_firstk_vec.begin(),dedx_firstk_vec.end());
	                                                                               std::sort(dedx_lastk_vec.begin(),dedx_lastk_vec.end());
										       std::sort(dedx_allk_vec.begin(),dedx_allk_vec.end());
										       float med_firstk=0;
										       if(dedx_firstk_vec.size()){
	                                                                                  if(dedx_firstk_vec.size()%2==1){
	                                                                                      med_firstk=dedx_firstk_vec[(dedx_firstk_vec.size()-1)/2];
	                                                                                  }
	                                                                                  if(dedx_firstk_vec.size()%2==0){
	                                                                                     med_firstk=float(dedx_firstk_vec[dedx_firstk_vec.size()/2] + dedx_firstk_vec[(dedx_firstk_vec.size()/2)-1])/2;
	                                                                                  }
	                                                                               }
										       float med_lastk=0;
										       if(dedx_lastk_vec.size()){
	                                                                                  if(dedx_lastk_vec.size()%2==1){
	                                                                                     med_lastk=dedx_lastk_vec[(dedx_lastk_vec.size()-1)/2];
	                                                                                  }
	                                                                                  if(dedx_lastk_vec.size()%2==0){
	                                                                                     med_lastk=float(dedx_lastk_vec[dedx_lastk_vec.size()/2] + dedx_lastk_vec[(dedx_lastk_vec.size()/2)-1])/2;
	                                                                                  }
	                                                                               }
										       float med_allk=0;
										       if(dedx_allk_vec.size()){
											  if(dedx_allk_vec.size()%2==1){
											     med_allk=dedx_allk_vec[(dedx_allk_vec.size()-1)/2];
											  }
											  if(dedx_allk_vec.size()%2==0){
											     med_allk=float(dedx_allk_vec[dedx_allk_vec.size()/2] + dedx_allk_vec[(dedx_allk_vec.size()/2)-1])/2;
											  }
										       }
										        /////////////// If you want min max inforamtions insert here /////////////
											float k_pida=0;
											float k_pida_sum=0;
											int k_used_points=0;
											for (size_t ical = 0; ical<calos.size(); ++ical){
											     if(!calos[ical]) continue;
			                                                                     if(!calos[ical]->PlaneID().isValid) continue;
			                                                                     int planenum = calos[ical]->PlaneID().Plane;
			                                                                     if(planenum<0||planenum>2) continue;
											     if(planenum==k_best_planenum){
											        const size_t NHits = calos[ical] -> dEdx().size();
												for(size_t iHit = 0; iHit < NHits; ++iHit){
												    if((calos[ical]->dEdx())[iHit]>100) continue;
												    if((calos[ical]->ResidualRange())[iHit] < 30){
												        float res=(calos[ical]->ResidualRange())[iHit];
													k_pida_sum=k_pida_sum+((calos[ical]->dEdx())[iHit])*(TMath::Power(res,0.42));
													k_used_points++;
												     }
												  }
											      }
											 }
											 
											 if(k_used_points!=0)k_pida=float(k_pida_sum)/k_used_points;		
											
								                         ///////////////////////////////////////// End of calorimetry information for kaons /////////////////////////////
											
											//////////////////////////////////////// Accessing calorimetry info of Muons ////////////////////////////////
											int mu_best_planenum=-1;
											calos = fmcal.at(dec_mu_index);
											int mu_plane_0_hits=-1;int mu_plane_1_hits=-1;int mu_plane_2_hits=-1;
											for(size_t ical = 0; ical<calos.size(); ++ical){
											    if(!calos[ical]) continue;
			                                                                    if(!calos[ical]->PlaneID().isValid) continue;
			                                                                    int planenum = calos[ical]->PlaneID().Plane;
			                                                                    if(planenum<0||planenum>2) continue;
							                                    if(planenum == 0) mu_plane_0_hits = int(calos[ical] -> dEdx().size());
			                                                                    if(planenum == 1) mu_plane_1_hits = int(calos[ical] -> dEdx().size());
		                                                                            if(planenum == 2) mu_plane_2_hits = int(calos[ical] -> dEdx().size());
											}
											
											if(mu_plane_0_hits != -1 || mu_plane_1_hits != -1 || mu_plane_2_hits != -1){
							                                if(mu_plane_0_hits > mu_plane_1_hits && mu_plane_0_hits > mu_plane_2_hits) mu_best_planenum=0;
                                                                                        if(mu_plane_1_hits > mu_plane_0_hits && mu_plane_1_hits > mu_plane_2_hits) mu_best_planenum=1;
                                                                                        if(mu_plane_2_hits > mu_plane_0_hits && mu_plane_2_hits > mu_plane_1_hits) mu_best_planenum=2;
                                                                                        if(mu_plane_2_hits==mu_plane_0_hits && mu_plane_2_hits > mu_plane_1_hits) mu_best_planenum=2;
                                                                                        if(mu_plane_2_hits==mu_plane_1_hits && mu_plane_2_hits > mu_plane_0_hits) mu_best_planenum=2;
                                                                                        if(mu_plane_1_hits==mu_plane_0_hits && mu_plane_1_hits > mu_plane_2_hits) mu_best_planenum=0;
                                                                                        if(mu_plane_1_hits==mu_plane_0_hits && mu_plane_1_hits==mu_plane_2_hits) mu_best_planenum=2;
							                                }
											
											float dedxmu_first_sum=0;
											float dedxmu_last_sum=0;
											int firstmu_nhits=0;
											int lastmu_nhits=0;
											vector<float>dedx_firstmu_vec;
									                vector<float>dedx_lastmu_vec; 
											vector<float>dedx_allmu_vec;
											for (size_t ical = 0; ical<calos.size(); ++ical){
											     if(!calos[ical]) continue;
			                                                                     if(!calos[ical]->PlaneID().isValid) continue;
			                                                                     int planenum = calos[ical]->PlaneID().Plane;
			                                                                     if(planenum<0||planenum>2) continue;
											     if(planenum==mu_best_planenum){
											         const size_t NHits = calos[ical] -> dEdx().size();
												 for(size_t iHit = 0; iHit < NHits; ++iHit){
												     if((calos[ical]->dEdx())[iHit]>100) continue;
												     dedx_allmu_vec.push_back((calos[ical]->dEdx())[iHit]);
												     if((calos[ical]->ResidualRange())[iHit] < 3){
												         dedxmu_first_sum=dedxmu_first_sum+(calos[ical]->dEdx())[iHit];
													 firstmu_nhits++;
													 dedx_firstmu_vec.push_back((calos[ical]->dEdx())[iHit]);
												     }
												     if(calos[ical]->Range()-(calos[ical]->ResidualRange())[iHit] < 3){
												        dedxmu_last_sum=dedxmu_last_sum+(calos[ical]->dEdx())[iHit];
													lastmu_nhits++;
													dedx_lastmu_vec.push_back((calos[ical]->dEdx())[iHit]);
												     }
												  }
											      }
											 }
											
											 float firstmu_dedx_av=0;
											 float lastmu_dedx_av=0;
											 if(firstmu_nhits!=0) firstmu_dedx_av=float(dedxmu_first_sum)/firstmu_nhits;
											 if(lastmu_nhits!=0) lastmu_dedx_av=float(dedxmu_last_sum)/lastmu_nhits;
											 std::sort(dedx_firstmu_vec.begin(),dedx_firstmu_vec.end());
	                                                                                 std::sort(dedx_lastmu_vec.begin(),dedx_lastmu_vec.end());
											 std::sort(dedx_allmu_vec.begin(),dedx_allmu_vec.end());
										         float med_firstmu=0;
										         if(dedx_firstmu_vec.size()){
	                                                                                    if(dedx_firstmu_vec.size()%2==1){
	                                                                                       med_firstmu=dedx_firstmu_vec[(dedx_firstmu_vec.size()-1)/2];
	                                                                                    }
	                                                                                    if(dedx_firstmu_vec.size()%2==0){
	                                                                                       med_firstmu=float(dedx_firstmu_vec[dedx_firstmu_vec.size()/2] + dedx_firstmu_vec[(dedx_firstmu_vec.size()/2)-1])/2;
	                                                                                  }
	                                                                               }
										       float med_lastmu=0;
										       if(dedx_lastmu_vec.size()){
	                                                                                  if(dedx_lastmu_vec.size()%2==1){
	                                                                                     med_lastmu=dedx_lastmu_vec[(dedx_lastmu_vec.size()-1)/2];
	                                                                                  }
	                                                                                  if(dedx_lastmu_vec.size()%2==0){
	                                                                                     med_lastmu=float(dedx_lastmu_vec[dedx_lastmu_vec.size()/2] + dedx_lastmu_vec[(dedx_lastmu_vec.size()/2)-1])/2;
	                                                                                  }
	                                                                               }
										       float med_allmu=0;
										       if(dedx_allmu_vec.size()){
										          if(dedx_allmu_vec.size()%2==1){
											     med_allmu=dedx_allmu_vec[(dedx_allmu_vec.size()-1)/2];
											  }
											  if(dedx_allmu_vec.size()%2==0){
											     med_allmu=float(dedx_allmu_vec[dedx_allmu_vec.size()/2] + dedx_allmu_vec[(dedx_allmu_vec.size()/2)-1])/2;
											  }
										       }
										       /////////////// If you want min max inforamtions insert here /////////////
											 float mu_pida=0;
											 float mu_pida_sum=0;
											 int mu_used_points=0;
											 for (size_t ical = 0; ical<calos.size(); ++ical){
											      if(!calos[ical]) continue;
			                                                                      if(!calos[ical]->PlaneID().isValid) continue;
			                                                                      int planenum = calos[ical]->PlaneID().Plane;
			                                                                      if(planenum<0||planenum>2) continue;
											      if(planenum==mu_best_planenum){
											         const size_t NHits = calos[ical] -> dEdx().size();
												 for(size_t iHit = 0; iHit < NHits; ++iHit){
												     if((calos[ical]->dEdx())[iHit]>100) continue;
												     if((calos[ical]->ResidualRange())[iHit] < 30){
												        float res=(calos[ical]->ResidualRange())[iHit];
													mu_pida_sum=mu_pida_sum+((calos[ical]->dEdx())[iHit])*(TMath::Power(res,0.42));
													mu_used_points++;
												     }
												  }
											      }
											 }
											 
											 if(mu_used_points!=0)mu_pida=float(mu_pida_sum)/mu_used_points;
											
											///////////////////////////////////// End of calorimetry info for muons ///////////////////////////////////
											
											//////////////////////////////////// defining k mu angle ///////////////////////////
											
											float average_commen_point_x=0;float average_commen_point_y=0;float average_commen_point_z=0;
									                float av_common_k_cos_alpha=0;float av_common_k_cos_beta=0;float av_common_k_cos_gamma=0;
									                float av_common_mu_cos_alpha=0;float av_common_mu_cos_beta=0;float av_common_mu_cos_gamma=0;
											
											if(st_k_st_d_mu < st_k_en_d_mu){
											   average_commen_point_x = float(p_k_st_x+d_mu_st_x)/2;
										           average_commen_point_y = float(p_k_st_y+d_mu_st_y)/2;
										           average_commen_point_z = float(p_k_st_z+d_mu_st_z)/2;
											   
											   float k_av_dis=TMath::Sqrt((average_commen_point_x-p_k_en_x)*(average_commen_point_x-p_k_en_x) + (average_commen_point_y-p_k_en_y)*(average_commen_point_y-p_k_en_y)+ (average_commen_point_z-p_k_en_z)*(average_commen_point_z-p_k_en_z));
										           float mu_av_dis=TMath::Sqrt((d_mu_en_x-average_commen_point_x)*(d_mu_en_x-average_commen_point_x) + (d_mu_en_y-average_commen_point_y)*(d_mu_en_y-average_commen_point_y) + (d_mu_en_z-average_commen_point_z)*(d_mu_en_z-average_commen_point_z));
											
											   av_common_k_cos_alpha  = float(average_commen_point_x-p_k_en_x)/k_av_dis; 
										           av_common_k_cos_beta   = float(average_commen_point_y-p_k_en_y)/k_av_dis; 
										           av_common_k_cos_gamma  = float(average_commen_point_z-p_k_en_z)/k_av_dis; 
										           av_common_mu_cos_alpha = float(d_mu_en_x-average_commen_point_x)/mu_av_dis; 
										           av_common_mu_cos_beta  = float(d_mu_en_y-average_commen_point_y)/mu_av_dis; 
										           av_common_mu_cos_gamma = float(d_mu_en_z-average_commen_point_z)/mu_av_dis; 
											}
											
											if(st_k_st_d_mu > st_k_en_d_mu){
											   average_commen_point_x = float(p_k_st_x+d_mu_en_x)/2;
										           average_commen_point_y = float(p_k_st_y+d_mu_en_y)/2;
										           average_commen_point_z = float(p_k_st_z+d_mu_en_z)/2;
										
										           float k_av_dis=TMath::Sqrt((average_commen_point_x-p_k_en_x)*(average_commen_point_x-p_k_en_x) + (average_commen_point_y-p_k_en_y)*(average_commen_point_y-p_k_en_y)+ (average_commen_point_z-p_k_en_z)*(average_commen_point_z-p_k_en_z));
										           float mu_av_dis=TMath::Sqrt((d_mu_st_x-average_commen_point_x)*(d_mu_st_x-average_commen_point_x) + (d_mu_st_y-average_commen_point_y)*(d_mu_st_y-average_commen_point_y) + (d_mu_st_z-average_commen_point_z)*(d_mu_st_z-average_commen_point_z));
										
										           av_common_k_cos_alpha  = float(average_commen_point_x-p_k_en_x)/k_av_dis; 
										           av_common_k_cos_beta   = float(average_commen_point_y-p_k_en_y)/k_av_dis; 
										           av_common_k_cos_gamma  = float(average_commen_point_z-p_k_en_z)/k_av_dis; 
										           av_common_mu_cos_alpha = float(d_mu_st_x-average_commen_point_x)/mu_av_dis; 
										           av_common_mu_cos_beta  = float(d_mu_st_y-average_commen_point_y)/mu_av_dis; 
										           av_common_mu_cos_gamma = float(d_mu_st_z-average_commen_point_z)/mu_av_dis; 
											}
											
											float av_common_open_angle = TMath::ACos(av_common_k_cos_alpha*av_common_mu_cos_alpha + av_common_k_cos_beta*av_common_mu_cos_beta +av_common_k_cos_gamma*av_common_mu_cos_gamma);
											n_reco_evt++;
											
											if(med_firstk-med_lastmu>=5){ ////// K/Mu Energy gradient cut **********************************
											   if(TMath::Cos(av_common_open_angle) <= 0.9 && TMath::Cos(av_common_open_angle) >=-0.9){ /////// K/Mu open angle cut ***********************
											      float mu_mom= (trkm.GetTrackMomentum(mud_tlen,13))*1000;
											       if(k_pida > 10 && k_pida < 20){ //// K PIDA cut ********* // k_pida > 10 && k_pida < 20
											          //if(med_allk >=2 && med_allk <=10){ // all K mean dE/dx
												  if(mu_pida < 10){ //// Mu PIDA cut *********
												  //if(med_allmu >=2 && med_allmu <=4){ // Mu all med dE/dx
												    //if(mu_mom*TMath::Sin(av_common_open_angle) < 200 && mu_mom*TMath::Sin(av_common_open_angle)> 150){ //// Mu mom cut ***********
											                if(break_indicator==0){
													    pri_mu_id = prim_mu_index;
													    pri_k_id = prim_kaon_index;
													    dec_mu_id = dec_mu_index;
													    break_indicator++;
											                    /*std::vector<const anab::Calorimetry*>*/ calos = fmcal.at(prim_kaon_index);
											                    for (size_t ical = 0; ical<calos.size(); ++ical){
											                         if (!calos[ical]) continue;
			                                                                                         if (!calos[ical]->PlaneID().isValid) continue;
			                                                                                         int planenum = calos[ical]->PlaneID().Plane;
			                                                                                         if (planenum<0||planenum>2) continue;
											                         if (planenum==k_best_planenum){
											                         const size_t NHits = calos[ical] -> dEdx().size();
											                         for(size_t iHit = 0; iHit < NHits; ++iHit){
											                             reco_kaon_cont_plane_two_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                 }
											                      }
											                   }
											                   reco_de_dx_gradient->Fill(firstk_dedx_av-lastmu_dedx_av);
											                   reco_k_mu_open_angle_cosine->Fill(TMath::Cos(av_common_open_angle));
											                   reco_kaon_trkln_hist->Fill(k_tlen);
											                   reco_muon_trkln_hist->Fill(mud_tlen);
											                   reco_mu_mom_range->Fill(mu_mom*TMath::Sin(av_common_open_angle));
													   if(lastk_dedx_av!=0)reco_k_dedx_ratio->Fill(float(firstk_dedx_av)/lastk_dedx_av);
													   if(lastmu_dedx_av!=0)reco_mu_dedx_ratio->Fill(float(firstmu_dedx_av)/lastmu_dedx_av);
													   reco_de_dx_median_gradient->Fill(med_firstk-med_lastmu);
													   if(med_lastk!=0)reco_k_median_dedx_ratio->Fill(float(med_firstk)/med_lastk);
													   if(med_lastmu!=0)reco_mu_median_dedx_ratio->Fill(float(med_firstmu)/med_lastmu);
											                   calos = fmcal.at(dec_mu_index);
											                   for (size_t ical = 0; ical<calos.size(); ++ical){
											                        if (!calos[ical]) continue;
			                                                                                        if (!calos[ical]->PlaneID().isValid) continue;
			                                                                                        int planenum = calos[ical]->PlaneID().Plane;
			                                                                                        if (planenum<0||planenum>2) continue;
											                        if (planenum==mu_best_planenum){
											                        const size_t NHits = calos[ical] -> dEdx().size();
											                        for(size_t iHit = 0; iHit < NHits; ++iHit){
											                            reco_muon_cont_plane_two_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                }
											                     }
											                  }
											                  reco_k_pida->Fill(k_pida);
											                  reco_mu_pida->Fill(mu_pida);
										                          reco_k_mean_first_dedx->Fill(firstk_dedx_av);
													  reco_k_mean_last_dedx->Fill(lastk_dedx_av);
													  reco_mu_mean_first_dedx->Fill(firstmu_dedx_av);
													  reco_mu_mean_last_dedx->Fill(lastmu_dedx_av);
													  reco_k_median_first_dedx->Fill(med_firstk);
													  reco_k_median_last_dedx->Fill(med_lastk);
													  reco_mu_median_first_dedx->Fill(med_firstmu);
													  reco_mu_median_last_dedx->Fill(med_lastmu);
													  reco_k_mean_len_ratio->Fill((float(firstk_dedx_av-lastk_dedx_av))/k_tlen);
													  reco_mu_mean_len_ratio->Fill((float(firstmu_dedx_av-lastmu_dedx_av))/mud_tlen);
													  reco_k_median_len_ratio->Fill((float(med_firstk-med_lastk))/k_tlen);
													  reco_mu_median_len_ratio->Fill((float(med_firstmu-med_lastmu))/mud_tlen);
													  reco_all_trk_median_dedx->Fill(med_allk);
													  reco_all_trmu_median_dedx->Fill(med_allmu);
													  k_best_plane_hist->Fill(k_best_planenum);
													  mu_best_plane_hist->Fill(mu_best_planenum);
											                  ////////////////////////////////// End of k mu angle //////////////////////////////
											                  
													  const sim::ParticleList& plist = bt->ParticleList();
                                                                                                          sim::ParticleList::const_iterator itPart = plist.begin(),pend = plist.end();
                                                                                                          std::string pri("primary");
										                          for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
											                      const simb::MCParticle* pPart = (itPart++)->second;
	                                                                                                      bool isPrimary = pPart->Process()== pri;
											                      if (isPrimary && pPart->PdgCode()==321){
											                      if (pPart->EndE()*1000 <= 494.){
												                  int trk_id  = pPart->TrackId();;
				                                                                                  float end_x = pPart->EndPosition()[0];
				                                                                                  float end_y = pPart->EndPosition()[1];
				                                                                                  float end_z = pPart->EndPosition()[2];
												                  sim::ParticleList::const_iterator jtPart = plist.begin(),send = plist.end();
											                          for(size_t jPart = 0; (jPart < plist.size()) && (jtPart != send); ++jPart){
												                      const simb::MCParticle* sPart = (jtPart++)->second;
												                      if ((sPart->Process()=="Decay") && (sPart->Mother()==trk_id) &&(std::abs(sPart->Vx()-end_x)<0.01) && (std::abs(sPart->Vy()-end_y)<0.01) && (std::abs(sPart->Vz()-end_z)<0.01)){
												                           if (sPart->PdgCode() == -13 || sPart->PdgCode() == 211){
													                       
															       ///////////////////////// Test /////////////////////
															       
													                       event_k_pida->Fill(k_pida);
													                       event_mu_pida->Fill(mu_pida);
													                       event_k_trklen->Fill(k_tlen);
													                       event_mu_trklen->Fill(mud_tlen);
													                       event_de_dx_median_gradient->Fill(med_firstk-med_lastmu);
													                       event_de_dx_mean_gradient->Fill(firstk_dedx_av-lastmu_dedx_av);
															       event_k_mu_open_angle_cosine->Fill(TMath::Cos(av_common_open_angle));
													                       event_mu_mom_range->Fill(mu_mom*TMath::Sin(av_common_open_angle));
															       std::vector<const anab::Calorimetry*> calos = fmcal.at(prim_kaon_index);
															       for(size_t ical = 0; ical<calos.size(); ++ical){
											                                           if(!calos[ical]) continue;
			                                                                                                           if(!calos[ical]->PlaneID().isValid) continue;
			                                                                                                           int planenum = calos[ical]->PlaneID().Plane;
			                                                                                                           if(planenum<0||planenum>2) continue;
											                                           if(planenum==k_best_planenum){
											                                              const size_t NHits = calos[ical] -> dEdx().size();
											                                              for(size_t iHit = 0; iHit < NHits; ++iHit){
													                                  event_kaon_cont_plane_best_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                                      }
												                                   }
											                                        }
													      
													                       //////////////////////// End of test ///////////////
															       
															       n_candidates++;
													                       std::vector< art::Ptr<recob::Hit> > allKHits = fmht.at(prim_kaon_index);
													                       std::vector< art::Ptr<recob::Hit> > allMuHits = fmht.at(dec_mu_index);
													                       std::map<int,double> trk_k_ide;
													                       for(size_t h = 0; h < allKHits.size(); ++h){
													                           art::Ptr<recob::Hit> hit = allKHits[h];
													                           std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
														                   for(size_t e = 0; e < TrackIDs.size(); ++e){
														                       trk_k_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
													                           }
													                         }
													                         double maxke = -1;
													                         double totke = 0;
													                         int Track_k_id = 0;
													                         for(std::map<int,double>::iterator ii=trk_k_ide.begin();ii!=trk_k_ide.end(); ++ii){
													                             totke += ii->second;
														                     if((ii->second)>maxke){
														                         maxke = ii->second;
													                                 Track_k_id=ii->first;
														                     }
														                  }
													                          std::map<int,double> trk_mu_ide;
														                  for(size_t h = 0; h < allMuHits.size(); ++h){
														                       art::Ptr<recob::Hit> hit = allMuHits[h];
														                       std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
														                       for(size_t e = 0; e < TrackIDs.size(); ++e){
															                   trk_mu_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
														                       }
													                           }
													                           double maxmue = -1;
														                   double totmue = 0;
														                   int Track_mu_id = 0;
														                   for(std::map<int,double>::iterator ii=trk_mu_ide.begin();ii!=trk_mu_ide.end(); ++ii){
														                       totmue += ii->second;
														                       if((ii->second)>maxmue){
															                   maxmue = ii->second;
														                           Track_mu_id=ii->first;
														                       }
														                    }
													                            if(Track_k_id==int(pPart->TrackId()) && Track_mu_id==int(sPart->TrackId())){
														                       std::cout << "K id : " << Track_k_id << "  " << int(pPart->TrackId()) << std::endl;
														                       std::cout << "Mu id : " << Track_mu_id << "  " << int(sPart->TrackId()) << std::endl;
																       
																       n_true_evt++;
														                       std::vector<const anab::Calorimetry*> calos = fmcal.at(prim_kaon_index);
											                                               for(size_t ical = 0; ical<calos.size(); ++ical){
											                                                   if(!calos[ical]) continue;
			                                                                                                                   if(!calos[ical]->PlaneID().isValid) continue;
			                                                                                                                   int planenum = calos[ical]->PlaneID().Plane;
			                                                                                                                   if(planenum<0||planenum>2) continue;
											                                                      if(planenum==k_best_planenum){
											                                                         const size_t NHits = calos[ical] -> dEdx().size();
											                                                         for(size_t iHit = 0; iHit < NHits; ++iHit){
													                                             true_kaon_cont_plane_best_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                                                 }
												                                              }
											                                                   }
														                           true_de_dx_gradient->Fill(firstk_dedx_av-lastmu_dedx_av);
														                           true_k_mu_open_angle_cosine->Fill(TMath::Cos(av_common_open_angle));
												                                           true_kaon_trkln_hist->Fill(k_tlen);
														                           true_muon_trkln_hist->Fill(mud_tlen);
														                           float mu_mom= (trkm.GetTrackMomentum(mud_tlen,13))*1000;
											                                                   true_mu_mom_range->Fill(mu_mom*TMath::Sin(av_common_open_angle));
																	   if(lastk_dedx_av!=0)true_k_dedx_ratio->Fill(float(firstk_dedx_av)/lastk_dedx_av);
																	   if(lastmu_dedx_av!=0)true_mu_dedx_ratio->Fill(float(firstmu_dedx_av)/lastmu_dedx_av);
																	   true_de_dx_median_gradient->Fill(med_firstk-med_lastmu);
																	   if(med_lastk!=0)true_k_median_dedx_ratio->Fill(float(med_firstk)/med_lastk);
																	   if(med_lastmu!=0)true_mu_median_dedx_ratio->Fill(float(med_firstmu)/med_lastmu);
														                           calos = fmcal.at(dec_mu_index);
											                                                   for(size_t ical = 0; ical<calos.size(); ++ical){
											                                                       if(!calos[ical]) continue;
			                                                                                                                       if(!calos[ical]->PlaneID().isValid) continue;
			                                                                                                                       int planenum = calos[ical]->PlaneID().Plane;
			                                                                                                                       if(planenum<0||planenum>2) continue;
											                                                       if(planenum==mu_best_planenum){
											                                                          const size_t NHits = calos[ical] -> dEdx().size();
											                                                          for(size_t iHit = 0; iHit < NHits; ++iHit){
													                                              true_muon_cont_plane_best_res_range_dEdX_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical]->dEdx())[iHit]);
												                                                  }
												                                               }
											                                                   }
														                           true_k_pida->Fill(k_pida);
														                           true_mu_pida->Fill(mu_pida);
																	   true_k_mean_first_dedx->Fill(firstk_dedx_av);
													                                   true_k_mean_last_dedx->Fill(lastk_dedx_av);
																	   reco_mu_mean_first_dedx->Fill(firstmu_dedx_av);
													                                   reco_mu_mean_last_dedx->Fill(lastmu_dedx_av);
																	   true_k_median_first_dedx->Fill(med_firstk);
													                                   true_k_median_last_dedx->Fill(med_lastk);
																	   true_mu_median_first_dedx->Fill(med_firstmu);
													                                   true_mu_median_last_dedx->Fill(med_lastmu);
																	   true_k_mean_len_ratio->Fill((float(firstk_dedx_av-lastk_dedx_av))/k_tlen);
																	   true_mu_mean_len_ratio->Fill((float(firstmu_dedx_av-lastmu_dedx_av))/mud_tlen);
																	   true_k_median_len_ratio->Fill((float(med_firstk-med_lastk))/k_tlen);
																	   true_mu_median_len_ratio->Fill((float(med_firstmu-med_lastmu))/mud_tlen);
																	   true_all_trk_median_dedx->Fill(med_allk);
																	   true_all_trmu_median_dedx->Fill(med_allmu);
													                                   } // trk g4 match
													                               } // getting true mu
													                           } // getting scttering particles
													                      } // 2 nd loop over geant list for P
											                                 } // getting stopping K
											                             } // getting primary K...
											                         } // 1 st loop over geant list for P...
										                             } // break indicator...
											                 //} // Mu mom cut ....
											                //} // Med all mu dE/dx...
												       } // Mu PIDA .....
											              //} // Med all K dE/dx...
												    } // K PIDA cut ......
											        //} // K/Mu min-max hit cut .....
											     //} // Muon tracklength cut ......
											    } // K/Mu open angle cut .....
										          } // K/Mu energy gradient cut.....
										      } // decay mu connected to start of K
										 } // primary mu is connected to end of K
									         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									    
									    } // decay mu is contained...
									 } // more than 10 hits for dec mu....
								    //} // 3rd loop over ical for i_2 track....
								} // trklen cut for decaying mu
							 } // 3 rd loop over tracks  
						      } // primay K and primay Mu distance matching....
						  } // pri K contained...
					      } // more than 10 hits for Pri K...
					 //} // 2 nd loop over ical for i_1 track....
				     } // pri K >= 10 cm
				 } // 2 nd loop over tracks...
			      } // more thans 10 hits for pri mu...
			 //} // 1 st loop over cal of i track...
		     //} // calorimetry info saved...
		 } // primary mu >= 60 mc ...
	    } // 1 st loop over tracks....
	 } // tracks > 3 in the event....
     } // track infomation is saved.....
  }// is not data (is MC) if statement.........
   
   reco_events = int(n_reco_evt);
   true_events=int(n_true_evt);
   candidate_k=int(n_candidates);
   break_value=int(break_indicator);
   fEventTree->Fill();

} // end of analyze function
	    
/////////////////// Defintion of reset function ///////////

void Kplane2::reset(){
     run = -99999;
     subrun = -99999;
     event = -99999;
     pri_mu_id = -99999;
     pri_k_id = -99999;
     dec_mu_id = -99999;
     reco_events=0;
     true_events=0;
     candidate_k=0;
     break_value=0;
}
/////////////////////// End of definition ///////////////	
	   
	   
DEFINE_ART_MODULE(Kplane2)
  
}



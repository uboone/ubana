////////////////////////////////////////////////////////////////////////
// Class:       PFPProfile
// Plugin Type: analyzer (art v3_03_01)
// File:        PFPProfile_module.cc
//
// Generated at Tue Nov 19 19:03:26 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_08_00.
// 
// Purpose of this module is to study the longitudinal and transverse profiles
// for different particles and use the profile information for particle ID
//
////////////////////////////////////////////////////////////////////////

// art include
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT include
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

constexpr int kNplanes = 3;      //number of wire planes
constexpr int kMaxPFPs = 1000;   //maximum number of PFParticles

class PFPProfile;

class PFPProfile : public art::EDAnalyzer {
public:
  explicit PFPProfile(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PFPProfile(PFPProfile const&) = delete;
  PFPProfile(PFPProfile&&) = delete;
  PFPProfile& operator=(PFPProfile const&) = delete;
  PFPProfile& operator=(PFPProfile&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // user's defined functions
  void project(const double& x0, 
               const double& y0, 
               const double& x1, 
               const double& y1, 
               const double& x, 
               const double& y, 
               double& xp, 
               double& yp);

  void reset();

  // configuration from fcl
  bool fUseMCInfo;

  // TTree
  TTree *fEventTree;
  
  // event information
  int event;
  int run;
  int subrun;

  // MCTruth
  int nuPDG_truth;
  int leptonPDG_truth;
  double leptonEnergy_truth; // [GeV]
  int ccnc_truth;
  double nuVertex_truth[4];
  double nuMomentum_truth[4]; // neutrino incoming momentum [GeV]
  double nuEnergy_truth; // neutrino incoming energy [GeV]


  int npfps;
  float LongProf[kMaxPFPs][3][100];
  float TranProf[kMaxPFPs][3][16];
  float TotalCharge[kMaxPFPs][3];
  int pfpid[kMaxPFPs];
  int trkid[kMaxPFPs];
  int shwid[kMaxPFPs];

};


PFPProfile::PFPProfile(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fUseMCInfo = p.get<bool>("UseMCInfo");
}

void PFPProfile::analyze(art::Event const& e)
{

  reset();

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
 
  // Get mctruth
  art::Handle< std::vector< simb::MCTruth > > mctruthListHandle;
  std::vector< art::Ptr< simb::MCTruth > > mctruthList;
  if (e.getByLabel("generator", mctruthListHandle)) {
    art::fill_ptr_vector(mctruthList, mctruthListHandle);
  }

  // Get all pfparticles
  art::Handle < std::vector < recob::PFParticle > > pfpListHandle;
  std::vector < art::Ptr < recob::PFParticle > > pfpList;
  if (e.getByLabel("pandora", pfpListHandle)) {
    art::fill_ptr_vector(pfpList, pfpListHandle);
  }

  // Get all clusters
  art::Handle < std::vector < recob::Cluster > > cluListHandle;
  std::vector < art::Ptr < recob::Cluster > > cluList;
  if (e.getByLabel("pandora", cluListHandle)) {
    art::fill_ptr_vector(cluList, cluListHandle);
  }


  // Get Track-PFParticle association
  art::FindManyP<recob::Track> fmtpfp(pfpListHandle, e, "pandora");

  // Get Shower-PFParticle association
  art::FindManyP<recob::Shower> fmspfp(pfpListHandle, e, "pandora");

  // Get cluster-PFParticle association
  art::FindManyP<recob::Cluster> fmcpfp(pfpListHandle, e, "pandora");

  // Get hit-cluster association
  art::FindManyP<recob::Hit> fmhc(cluListHandle, e, "pandora");

  // Get Geometry service
  art::ServiceHandle<geo::Geometry const> geom;

  // Get DetectorPropertiesService
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // MCInfo
  if (fUseMCInfo) {
    // mctruth
    for (size_t i = 0; i < mctruthList.size(); ++i) {
      art::Ptr<simb::MCTruth> mctruth = mctruthList[i];
      if (mctruth->NeutrinoSet()) {
        simb::MCNeutrino nu = mctruth->GetNeutrino();
        simb::MCParticle neutrino = nu.Nu(); // incoming neutrino
        simb::MCParticle lepton = nu.Lepton(); // outgoing lepton

        nuPDG_truth = neutrino.PdgCode();
        leptonPDG_truth = lepton.PdgCode();
        leptonEnergy_truth = lepton.E(); // [GeV]
        ccnc_truth = nu.CCNC(); // 0: CC; 1: NC

        const TLorentzVector & vertex = neutrino.Position(0);
        vertex.GetXYZT(nuVertex_truth); 
        
        const TLorentzVector & nu_momentum = neutrino.Momentum(0);
        nu_momentum.GetXYZT(nuMomentum_truth);
        nuEnergy_truth = nuMomentum_truth[3];

        for (int i = 0; i < 4; i++) {
          cout << "nuVertex_truth[" << i << "]: " << nuVertex_truth[i] << endl;
          cout << "nuMomentum_truth[" << i << "]: " << nuMomentum_truth[i] << endl;
        }
      }
    }
  } // end of if (fUseMCInfo)
  


  // Conversion factor from tick to distance
  double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns  

  // Loop over all pfparticles
  npfps = 0;
  for (const auto& pfp : pfpList){

    if (npfps >= kMaxPFPs) continue;
    pfpid[npfps] = pfp.key();

    // Find the vertex and direction of the pfparticle
    double vtx[3] = {0,0,0};
    double dir[3] = {0,0,0};
    bool foundvtxdir = false;

    // First check tracks
    if (fmtpfp.isValid()){
      auto const& tracks = fmtpfp.at(pfp.key());
      if (!tracks.empty()){
        vtx[0] = tracks[0]->Start().X();
        vtx[1] = tracks[0]->Start().Y();
        vtx[2] = tracks[0]->Start().Z();

        dir[0] = tracks[0]->StartDirection().X();
        dir[1] = tracks[0]->StartDirection().Y();
        dir[2] = tracks[0]->StartDirection().Z();
        
        trkid[npfps] = tracks[0].key();
        foundvtxdir = true;
      }
    }
    // If no track is found, check showers
    if (!foundvtxdir && fmspfp.isValid()){
      auto const& showers = fmspfp.at(pfp.key());
      if (!showers.empty()){
        vtx[0] = showers[0]->ShowerStart().X();
        vtx[1] = showers[0]->ShowerStart().Y();
        vtx[2] = showers[0]->ShowerStart().Z();

        dir[0] = showers[0]->Direction().X();
        dir[1] = showers[0]->Direction().Y();
        dir[2] = showers[0]->Direction().Z();
        
        shwid[npfps] = showers[0].key();
        foundvtxdir = true;
      }
    }

    if (foundvtxdir){// foundvtxdir
      // Get hits on each plane
      std::vector<std::vector<art::Ptr<recob::Hit>>> allhits(3);
      if (fmcpfp.isValid()){
        // Get clusters associated with pfparticle
        auto const& clusters = fmcpfp.at(pfp.key());
        for (auto const & cluster : clusters){
          if (fmhc.isValid()){
            // Get hits associated with cluster
            auto const& hits = fmhc.at(cluster.key());
            for (auto const & hit: hits){
              if (bool(hit->WireID())){// wire id valid
                allhits[hit->WireID().Plane].push_back(hit);
              }
            }
          }
        }
      }

      // Loop over all planes
      for (unsigned short pl = 0; pl <geom->Nplanes(); ++pl){
        if (allhits[pl].empty()) continue; //no hits on this plane
        double wirePitch = geom->WirePitch(pl);
        // Project vertex and diection on the wire plane
        double w0 = geom->WireCoordinate(vtx[1], vtx[2], pl, 0, 0);
        double t0 = detprop->ConvertXToTicks(vtx[0], pl, 0, 0);
        double w1 = geom->WireCoordinate(vtx[1]+dir[1], vtx[2]+dir[2], pl, 0, 0);
        double t1 = detprop->ConvertXToTicks(vtx[0]+dir[0], pl, 0, 0);
        // Convert wire number and tick to cm
        double w0_cm = w0*wirePitch;
        double t0_cm = t0*tickToDist;
        double w1_cm = w1*wirePitch;
        double t1_cm = t1*tickToDist;
        // Save 3D information of all hits
        std::vector<double> hitx;
        std::vector<double> hity;
        std::vector<double> hitz;
        std::vector<double> hitdist;
        std::vector<double> hitcharge;
        // Find first point of pfparticle
        double x0 = -DBL_MAX;
        double y0 = -DBL_MAX;
        double z0 = -DBL_MAX;
        double rmin = DBL_MAX;
        // Loop over all hits on the current plane
        for (auto const & hit : allhits[pl]){
          double w_cm = hit->WireID().Wire*wirePitch;
          double t_cm = hit->PeakTime()*tickToDist;
          // Project the hit onto pfparticle projection
          double wp_cm = -DBL_MAX;
          double tp_cm = -DBL_MAX;
          project(w0_cm, t0_cm, w1_cm, t1_cm, w_cm, t_cm, wp_cm, tp_cm);
          // Distance between hit and its projection
          // This will be used for transverse profile
          double dist = sqrt((wp_cm-w_cm)*(wp_cm-w_cm)+(tp_cm-t_cm)*(tp_cm-t_cm));
          // Ratio
          double r = sqrt((wp_cm-w0_cm)*(wp_cm-w0_cm)+(tp_cm-t0_cm)*(tp_cm-t0_cm))/
            sqrt((w1_cm-w0_cm)*(w1_cm-w0_cm)+(t1_cm-t0_cm)*(t1_cm-t0_cm));
          // Determine if the hit is before and after the vertex
          double sign = 1.;
          if ((wp_cm-w0_cm)*(w1_cm-w0_cm)+(tp_cm-t0_cm)*(t1_cm-t0_cm)<0) sign = -1;
          // x,y,z are the 3D coordinates of the hit projection
          double x = vtx[0]+dir[0]*r*sign;
          double y = vtx[1]+dir[1]*r*sign;
          double z = vtx[2]+dir[2]*r*sign;
          if (r*sign<rmin){
            rmin = r*sign;
            x0 = x;
            y0 = y;
            z0 = z;
          }
          hitx.push_back(x);
          hity.push_back(y);
          hitz.push_back(z);
          hitdist.push_back(dist);
          hitcharge.push_back(hit->Integral());
        }
        // Now calculate longitudinal and transverse profiles
        for (size_t i = 0; i<hitx.size(); ++i){
          double Ldist = sqrt((hitx[i]-x0)*(hitx[i]-x0)+
                              (hity[i]-y0)*(hity[i]-y0)+
                              (hitz[i]-z0)*(hitz[i]-z0));
          double Tdist = hitdist[i];
          //if (pfp.key()==0) std::cout<<pfp.key()<<" "<<pl<<" "<<Ldist<<" "<<Tdist<<" "<<hitcharge[i]<<std::endl;
          int iL = int(Ldist/(14./4.)); //0.25 radiation length
          int iT = int(Tdist/0.5);     //0.5 cm
          //std::cout<<npfps<<" "<<pl<<" "<<Ldist<<" "<<iL<<" "<<hitcharge[i]<<std::endl;
          //std::cout<<npfps<<" "<<pl<<" "<<Tdist<<" "<<iT<<" "<<hitcharge[i]<<std::endl;
          if (iL>=0 && iL<100) LongProf[npfps][pl][iL] += hitcharge[i];
          if (iT>=0 && iT<16)  TranProf[npfps][pl][iT] += hitcharge[i];
          TotalCharge[npfps][pl] += hitcharge[i];
        }
      }
      ++ npfps;
    } // if foundvtxdir
    else{
      std::cout<<"Could not find vertex and direction."<<std::endl;
    }
  }// loop over pfps
  fEventTree->Fill();
}

void PFPProfile::project(const double& x0, 
                         const double& y0, 
                         const double& x1, 
                         const double& y1, 
                         const double& x, 
                         const double& y, 
                         double& xp, 
                         double& yp){

  //Project a point on a line, based on
  //http://www.vcskicks.com/code-snippet/point-projection.php

  if (std::abs(x0-x1)<1e-6){
    xp = x0;
    yp = y;
  }
  else{
    double m = (y1 - y0) / (x1 - x0);
    double b = y0 - (m * x0);
    
    xp = (m * y + x - m * b) / (m * m + 1);
    yp = (m * m * y + m * x + b) / (m * m + 1);
  }
}

void PFPProfile::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event");
  fEventTree->Branch("event", &event, "event/I");
  fEventTree->Branch("run", &run, "run/I");
  fEventTree->Branch("subrun", &subrun, "subrun/I");

  fEventTree->Branch("nuPDG_truth", &nuPDG_truth, "nuPDG_truth/I");
  fEventTree->Branch("leptonPDG_truth", &leptonPDG_truth, "leptonPDG_truth/I");
  fEventTree->Branch("leptonEnergy_truth", &leptonEnergy_truth, "leptonEnergy_truth/D");
  fEventTree->Branch("ccnc_truth", &ccnc_truth, "ccnc_truth/I");
  fEventTree->Branch("nuVertex_truth", &nuVertex_truth, "nuVertex_truth[4]/D");
  fEventTree->Branch("nuMomentum_truth", &nuMomentum_truth, "nuMomentum_truth[4]/D");
  fEventTree->Branch("nuEnergy_truth", &nuEnergy_truth, "nuEnergy_truth/D");

  fEventTree->Branch("npfps", &npfps,"npfps/I");
  fEventTree->Branch("pfpid", pfpid, "pfpid[npfps]/I");
  fEventTree->Branch("trkid", trkid, "trkid[npfps]/I");
  fEventTree->Branch("shwid", shwid, "shwid[npfps]/I");
  fEventTree->Branch("TotalCharge", TotalCharge, "TotalCharge[npfps][3]/F");
  fEventTree->Branch("LongProf", LongProf, "LongProf[npfps][3][100]/F");
  fEventTree->Branch("TranProf", TranProf, "TranProf[npfps][3][16]/F");

}

void PFPProfile::reset() {

  run = -99999;
  subrun = -99999;
  event = -99999;
  
  npfps = -99999;
  for (size_t i = 0; i<kMaxPFPs; ++i){
    pfpid[i] = -1;
    trkid[i] = -1;
    shwid[i] = -1;
    for (size_t j = 0; j<3; ++j){
      TotalCharge[i][j] = 0;
      for (size_t k = 0; k<100; ++k) LongProf[i][j][k] = 0;
      for (size_t k = 0; k<16; ++k) TranProf[i][j][k] = 0;
    }
  }
}

DEFINE_ART_MODULE(PFPProfile)

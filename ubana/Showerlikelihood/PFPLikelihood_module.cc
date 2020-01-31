////////////////////////////////////////////////////////////////////////
// Class:       PFPLikelihood
// Plugin Type: analyzer (art v3_03_01)
// File:        PFPLikelihood_module.cc
//
// Generated at Fri Jan  3 08:54:50 2020 by Wanwei Wu using cetskelgen
// from cetlib version v3_08_00.
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
#include "canvas/Persistency/Common/FindMany.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"  // to calculate the time offset for x
#include "larevt/SpaceChargeServices/SpaceChargeService.h" // for SCE correction
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

constexpr int kNenergies = 10;   // selected energies for profile


class PFPLikelihood;

class PFPLikelihood : public art::EDAnalyzer {
public:
  explicit PFPLikelihood(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PFPLikelihood(PFPLikelihood const&) = delete;
  PFPLikelihood(PFPLikelihood&&) = delete;
  PFPLikelihood& operator=(PFPLikelihood const&) = delete;
  PFPLikelihood& operator=(PFPLikelihood&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  bool insideFV(double vertex[3]);

private:

  // Declare member data here.
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
  bool fFVCutOnRecon; // true: use recon; false: use MC
  float fFidVolCutX;
  float fFidVolCutY;
  float fFidVolCutZUP;
  float fFidVolCutZDOWN;
  std::vector<double> fADCtoE; // calibration constants
  double fRecombination;
  bool fUseEnergyCorrection;
  std::vector<double> fEnergyCorrectionSlope;
  std::vector<double> fEnergyCorrectionIntercept;
  float fCutPurity;
  float fCutCompleteness;
  std::string fProfileTemplateFile;

  // Fiducial volume
  float fFidVolXmin;
  float fFidVolXmax;
  float fFidVolYmin;
  float fFidVolYmax;
  float fFidVolZmin;
  float fFidVolZmax;

  // template
  TProfile *profile_long_electron[kNenergies]; // [energies]
  TProfile *profile_long_photon[kNenergies]; // [energies]

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
  int MC_lepton_ID; // primary lepton TrackId

  // PFP
  int npfps;
  float LongProf[kMaxPFPs][3][100];
  float TranProf[kMaxPFPs][3][16];
  float TranProf_1[kMaxPFPs][3][16];
  float TranProf_2[kMaxPFPs][3][16];
  float TranProf_3[kMaxPFPs][3][16];
  float TranProf_4[kMaxPFPs][3][16];
  float TranProf_5[kMaxPFPs][3][16];
  float TotalCharge[kMaxPFPs][3];
  int pfpid[kMaxPFPs];
  int pfpself[kMaxPFPs];
  int trkkey[kMaxPFPs];
  int trkid[kMaxPFPs];
  int shwkey[kMaxPFPs];
  int shwid[kMaxPFPs];
  double pfpvertex_recon[kMaxPFPs][3];
  double pfpvertex_truth[kMaxPFPs][3]; // from the pfp particle start
  double pfpend_truth[kMaxPFPs][3]; // from the pfp particle end
  double pfpdir_recon[kMaxPFPs][3];
  double pfpdir_truth[kMaxPFPs][3]; // from the pfp particle
  double pfp_datahits[kMaxPFPs][3]; // fraction of hits from data in a pfp
  int pfp_primary_e[kMaxPFPs]; // whether a pfp is from primary electron. 0: false; 1:true
  int pfppdg[kMaxPFPs]; // from PdgCode() funciton of PFParticle
  int pfppdg_truth[kMaxPFPs]; // from matching hits and MCParticle
  double pfpcharge[kMaxPFPs][3];
  double pfpenergy_mc[kMaxPFPs]; // pfp particle's kinetic energy [MeV]
  double pfpenergy_recon[kMaxPFPs][3]; // for three planes
  double pfppurity[kMaxPFPs];
  double pfp_purity[kMaxPFPs][3]; // for each planes
  double pfpcompleteness[kMaxPFPs];
  double pfp_completeness[kMaxPFPs][3]; // for each planes
  // likelihood
  double pfp_longlikelihood_electron[kMaxPFPs]; // use plane 2 and electron profile template only
  double pfp_longlikelihood_photon[kMaxPFPs]; // use plane 2 and photon profile template only

};


PFPLikelihood::PFPLikelihood(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fUseMCInfo = p.get<bool>("UseMCInfo");
  fFVCutOnRecon       = p.get<bool>("FVCutOnRecon");
  fFidVolCutX         = p.get<float>("FidVolCutX");
  fFidVolCutY         = p.get<float>("FidVolCutY");
  fFidVolCutZUP         = p.get<float>("FidVolCutZUP");
  fFidVolCutZDOWN         = p.get<float>("FidVolCutZDOWN");
  fADCtoE = p.get<std::vector<double>>("ADCtoE");
  fRecombination = p.get<double>("Recombination");
  fUseEnergyCorrection = p.get<bool>("UseEnergyCorrection");
  fEnergyCorrectionSlope = p.get<std::vector<double>>("EnergyCorrectionSlope");
  fEnergyCorrectionIntercept = p.get<std::vector<double>>("EnergyCorrectionIntercept");
  fCutPurity    = p.get<float>("CutPurity");
  fCutCompleteness = p.get<float>("CutCompleteness");
  fProfileTemplateFile = p.get<std::string>("ProfileTemplateFile");
}

void PFPLikelihood::analyze(art::Event const& e)
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

  // Get all hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector< art::Ptr<recob::Hit> > hitList;
  if(e.getByLabel("gaushit",hitListHandle)){
    art::fill_ptr_vector(hitList, hitListHandle);
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

  // Get MCParticle-hit association
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, e, "gaushitTruthMatch");

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

  // Get ParticleInventoryService
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;

  // hits-MCParticle
  std::vector<std::unordered_map<int, double>> hit_trkide(3);

  // MCInfo
  if (fUseMCInfo) {
    if (mctruthList.size() > 1) {
      std::cout << "Warning:  mctruthList.size() = " << mctruthList.size() << "; more than one mctruth in the event" << std::endl;
    }
    for (size_t i = 0; i < mctruthList.size(); ++i) {
      // mctruth: generator
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

        if (!fFVCutOnRecon) {
          double check_vertex[3];
          check_vertex[0] = nuVertex_truth[0];
          check_vertex[1] = nuVertex_truth[1];
          check_vertex[2] = nuVertex_truth[2];

          if (!insideFV(check_vertex)) {
            cout << "\n **** Interaction is NOT inside the Fiducial Volume. RETURN ****" << endl;
            return;
          }
        }
      }

      // mcparticle: Geant4
      const sim::ParticleList& plist = piserv->ParticleList();
      for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar) {
        simb::MCParticle *particle = 0;
        particle = ipar->second;
        auto & truth = piserv->ParticleToMCTruth_P(particle);
        // primary lepton
        if (truth->Origin() == simb::kBeamNeutrino
            && particle->Mother() == 0
            && std::abs(particle->PdgCode()) >= 11
            && std::abs(particle->PdgCode()) <= 16) { 
          MC_lepton_ID = particle->TrackId();
          break;
        }
      }

      // hit information
      std::vector<simb::MCParticle const *> hit_particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const *> hit_match_vec;
      //std::vector<std::unordered_map<int, double>> hit_trkide(3);
      hit_trkide.clear();
      for (size_t i = 0; i < hitList.size(); ++i) {
        art::Ptr<recob::Hit> hit = hitList[i];
        hit_particle_vec.clear();
        hit_match_vec.clear();
        particles_per_hit.get(hit.key(), hit_particle_vec, hit_match_vec);
        for (size_t i_p = 0; i_p < hit_particle_vec.size(); ++i_p) {
          if (piserv->TrackIdToEveTrackId(hit_particle_vec.at(i_p)->TrackId())) {
            hit_trkide[hit->WireID().Plane][piserv->TrackIdToEveTrackId(hit_particle_vec.at(i_p)->TrackId())] += hit_match_vec[i_p]->energy; // energy deposited by ionization by this track ID [MeV] 
          }
        }
      }

    } // end of for mctruthList.size() i; in principle only one mctruch each simulated event
  } // end of if (fUseMCInfo)
  
  // Conversion factor from tick to distance
  double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns  

  // Loop over all pfparticles
  npfps = 0;
  for (const auto& pfp : pfpList){

    if (npfps >= kMaxPFPs) continue;
    pfpid[npfps] = pfp.key();
    pfpself[npfps] = pfp->Self();
    pfppdg[npfps] = pfp->PdgCode();
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
        
        trkkey[npfps] = tracks[0].key();
        trkid[npfps] = tracks[0]->ID();
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
        
        shwkey[npfps] = showers[0].key();
        shwid[npfps] = showers[0]->ID();
        foundvtxdir = true;
      }
    }

    if (foundvtxdir){// foundvtxdir
      if (fFVCutOnRecon) {
        if (!insideFV(vtx)) {
          cout << "\n **** Interaction is NOT inside the Fiducial Volume. RETURN ****" << endl;
          return;
        }
      }
      // Save reconstructed vertex and direction
      for (int xyz = 0; xyz < 3; xyz++) {
        pfpvertex_recon[npfps][xyz] = vtx[xyz];
        pfpdir_recon[npfps][xyz] = dir[xyz];
      }
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

      // Associate a particle to pfp by considering hits on all planes
      if (fUseMCInfo) {
        std::vector<simb::MCParticle const *> particle_vec;
        std::vector<anab::BackTrackerHitMatchingData const *> match_vec;
        std::unordered_map<int, double> pfphit_trkide;
        std::vector<std::unordered_map<int, double>> plane_pfphit_trkide(3);
        pfphit_trkide.clear();
        plane_pfphit_trkide.clear();
        
        // Do not save data pfp: apply cuts on count_datahits / totalhits 
        for (size_t p = 0; p < 3; p++) {
          int count_datahits = 0;
          for (size_t ihit = 0; ihit < allhits[p].size(); ihit++) {

            pfpcharge[npfps][p] += allhits[p][ihit]->Integral();
            particle_vec.clear();
            match_vec.clear();
            particles_per_hit.get(allhits[p][ihit].key(), particle_vec, match_vec);
            if (particle_vec.empty()) count_datahits++; 
            for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
              if (piserv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId())) {
                pfphit_trkide[piserv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId())] += match_vec[i_p]->energy;
                plane_pfphit_trkide[p][piserv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId())] += match_vec[i_p]->energy;
              }
            }
          }
          if (!allhits[p].empty()) {
            pfp_datahits[npfps][p] = count_datahits * 1.0 / allhits[p].size();
          }
        }

        int pfp_trackID = 0;
        double pfp_particle_E = 0; // pfp energy on all planes from the "particle"
        double pfp_total_E = 0; // pfp energy on all planes
        for (auto const & p_trkide : pfphit_trkide) {
          pfp_total_E += p_trkide.second;
          if (p_trkide.second > pfp_particle_E) {
            pfp_particle_E = p_trkide.second;
            pfp_trackID = p_trkide.first;
          }
        }
     
        double plane_pfp_particle_E[3] = {0, 0, 0};
        double plane_pfp_total_E[3] = {0, 0, 0};
        for (size_t p = 0; p < 3; p++) {
          pfpenergy_recon[npfps][p] = pfpcharge[npfps][p]*fADCtoE[p]/fRecombination*23.6e-6; // [MeV] 
          for (auto const & plane_trkide : plane_pfphit_trkide[p]) {
            plane_pfp_total_E[p] += plane_trkide.second;
            if (plane_trkide.second > plane_pfp_particle_E[p]) {
              plane_pfp_particle_E[p] = plane_trkide.second;
            }
          }
        }

        if (pfp_trackID) {
          if (piserv->TrackIdToParticle_P(pfp_trackID)) {
            const simb::MCParticle* pfp_particle = piserv->TrackIdToParticle_P(pfp_trackID);
            pfppdg_truth[npfps] = pfp_particle->PdgCode();
            if ( (abs(pfp_particle->PdgCode()) == 11) && (pfp_trackID == MC_lepton_ID) ) {
              pfp_primary_e[npfps] = 1;
            }

            pfpvertex_truth[npfps][0] = pfp_particle->Vx();
            pfpvertex_truth[npfps][1] = pfp_particle->Vy();
            pfpvertex_truth[npfps][2] = pfp_particle->Vz();
            pfpend_truth[npfps][0] = pfp_particle->EndX();
            pfpend_truth[npfps][1] = pfp_particle->EndY();
            pfpend_truth[npfps][2] = pfp_particle->EndZ();
            
            double momentum = std::sqrt(pfp_particle->Px()*pfp_particle->Px() + pfp_particle->Py()*pfp_particle->Py() + pfp_particle->Pz()*pfp_particle->Pz());
            pfpdir_truth[npfps][0] = pfp_particle->Px() / momentum;
            pfpdir_truth[npfps][1] = pfp_particle->Py() / momentum;
            pfpdir_truth[npfps][2] = pfp_particle->Pz() / momentum;
            pfpenergy_mc[npfps] = (pfp_particle->E() - pfp_particle->Mass())*1000.0; // kinetic energy [MeV]
            double particle_total_E = hit_trkide[0][pfp_trackID] + hit_trkide[1][pfp_trackID] + hit_trkide[2][pfp_trackID];
            if (pfp_particle_E) {
              pfppurity[npfps] = pfp_particle_E / pfp_total_E;       
              if (particle_total_E) {
                pfpcompleteness[npfps] = pfp_particle_E / particle_total_E;
              }
            }
            for (size_t p=0; p<3; p++) {
              if (plane_pfp_particle_E[p]) {
                pfp_purity[npfps][p] = plane_pfp_particle_E[p] / plane_pfp_total_E[p];
                if (hit_trkide[p][pfp_trackID]) {
                  pfp_completeness[npfps][p] = plane_pfp_particle_E[p] / hit_trkide[p][pfp_trackID];
                }
              }
            }
          }
        }
      } // end of fUseMCInfo

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
          if (iT>=0 && iT<16)  {
            TranProf[npfps][pl][iT] += hitcharge[i];
            if (Ldist/14. < 1) TranProf_1[npfps][pl][iT] += hitcharge[i];
            else if (Ldist/14. < 2) TranProf_2[npfps][pl][iT] += hitcharge[i];
            else if (Ldist/14. < 3) TranProf_3[npfps][pl][iT] += hitcharge[i];
            else if (Ldist/14. < 4) TranProf_4[npfps][pl][iT] += hitcharge[i];
            else if (Ldist/14. < 5) TranProf_5[npfps][pl][iT] += hitcharge[i];
          }
          
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

void PFPLikelihood::project(const double& x0,
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

void PFPLikelihood::beginJob(){

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
  fEventTree->Branch("pfpself", pfpself, "pfpself[npfps]/I");
  fEventTree->Branch("trkkey", trkkey, "trkkey[npfps]/I");
  fEventTree->Branch("trkid", trkid, "trkid[npfps]/I");
  fEventTree->Branch("shwkey", shwkey, "shwkey[npfps]/I");
  fEventTree->Branch("shwid", shwid, "shwid[npfps]/I");
  fEventTree->Branch("pfp_primary_e", pfp_primary_e, "pfp_primary_e[npfps]/I");
  fEventTree->Branch("pfppdg", pfppdg, "pfppdg[npfps]/I");
  fEventTree->Branch("pfppdg_truth", pfppdg_truth, "pfppdg_truth[npfps]/I");
  fEventTree->Branch("pfpvertex_recon", pfpvertex_recon, "pfpvertex_recon[npfps][3]/D");
  fEventTree->Branch("pfpvertex_truth", pfpvertex_truth, "pfpvertex_truth[npfps][3]/D");
  fEventTree->Branch("pfpend_truth", pfpend_truth, "pfpend_truth[npfps][3]/D");
  fEventTree->Branch("pfpdir_recon", pfpdir_recon, "pfpdir_recon[npfps][3]/D");
  fEventTree->Branch("pfpdir_truth", pfpdir_truth, "pfpdir_truth[npfps][3]/D");
  fEventTree->Branch("pfp_datahits", pfp_datahits, "pfp_datahits[npfps][3]/D");
  fEventTree->Branch("pfpenergy_mc", pfpenergy_mc, "pfpenergy_mc[npfps]/D");
  fEventTree->Branch("pfpenergy_recon", pfpenergy_recon, "pfpenergy_recon[npfps][3]/D");
  fEventTree->Branch("pfppurity", pfppurity, "pfppurity[npfps]/D");
  fEventTree->Branch("pfp_purity", pfp_purity, "pfp_purity[npfps][3]/D");
  fEventTree->Branch("pfpcompleteness", pfpcompleteness, "pfpcompleteness[npfps]/D");
  fEventTree->Branch("pfp_completeness", pfp_completeness, "pfp_completeness[npfps][3]/D");
  fEventTree->Branch("TotalCharge", TotalCharge, "TotalCharge[npfps][3]/F");
  fEventTree->Branch("LongProf", LongProf, "LongProf[npfps][3][100]/F");
  fEventTree->Branch("TranProf", TranProf, "TranProf[npfps][3][16]/F");
  fEventTree->Branch("TranProf_1", TranProf_1, "TranProf_1[npfps][3][16]/F");
  fEventTree->Branch("TranProf_2", TranProf_2, "TranProf_2[npfps][3][16]/F");
  fEventTree->Branch("TranProf_3", TranProf_3, "TranProf_3[npfps][3][16]/F");
  fEventTree->Branch("TranProf_4", TranProf_4, "TranProf_4[npfps][3][16]/F");
  fEventTree->Branch("TranProf_5", TranProf_5, "TranProf_5[npfps][3][16]/F");

  fEventTree->Branch("pfp_longlikelihood_electron", pfp_longlikelihood_electron, "pfp_longlikelihood_electron[npfps]/D");
  fEventTree->Branch("pfp_longlikelihood_photon", pfp_longlikelihood_photon, "pfp_longlikelihood_photon[npfps]/D");

  // Get the FV cut information
  // uboone geometry
  double minX = 0.;
  double maxX = 256.;
  double minY = -116.;
  double maxY = 116.;
  double minZ = 0.;
  double maxZ = 1036.;

  fFidVolXmin = minX + fFidVolCutX;
  fFidVolXmax = maxX - fFidVolCutX;
  fFidVolYmin = minY + fFidVolCutY;
  fFidVolYmax = maxY - fFidVolCutY;
  fFidVolZmin = minZ + fFidVolCutZUP;
  fFidVolZmax = maxZ - fFidVolCutZDOWN;
  cout << "Fiducial Volume cuts on vertex:" << endl;
  cout << fFidVolXmin << " =< x =< " << fFidVolXmax << endl;
  cout << fFidVolYmin << " =< y =< " << fFidVolYmax << endl;
  cout << fFidVolZmin << " =< z =< " << fFidVolZmax << endl;

  // Profile Template
  const char * energy_options[10] = {"0_100MeV", "100_200MeV", "200_300MeV","300_400MeV", "400_500MeV","500_600MeV", "600_700MeV","700_800MeV", "800_900MeV","900_1000MeV"};

  TFile *pfpProfileTemplate = TFile::Open(fProfileTemplateFile.c_str());
  for (int j=0; j<kNenergies;j++) {
    profile_long_electron[j] = (TProfile*)pfpProfileTemplate->Get(TString::Format("hprofile_electron_%s_long", energy_options[j]));
    profile_long_photon[j] = (TProfile*)pfpProfileTemplate->Get(TString::Format("hprofile_photon_%s_long", energy_options[j]));
  }

}

bool PFPLikelihood::insideFV( double vertex[3] ) {
  if ( vertex[0] >= fFidVolXmin && vertex[0] <= fFidVolXmax &&
      vertex[1] >= fFidVolYmin && vertex[1] <= fFidVolYmax &&
      vertex[2] >= fFidVolZmin && vertex[2] <= fFidVolZmax ) {
    return true;
  }
  else {
    return false;
  }
}

void PFPLikelihood::reset(){

  run = -99999;
  subrun = -99999;
  event = -99999;

  MC_lepton_ID = 0;

  npfps = -99999;
  for (size_t i = 0; i<kMaxPFPs; ++i){
    pfpid[i] = -1;
    pfpself[i] = -1;
    trkkey[i] = -1;
    trkid[i] = -1;
    shwkey[i] = -1;
    shwid[i] = -1;
    pfp_primary_e[i] = 0;
    pfppdg[i] = -1;
    pfppdg_truth[i] = -1;
    pfpenergy_mc[i] = 0;
    pfppurity[i] = -1.0;
    pfpcompleteness[i] = -1.0;
    pfp_longlikelihood_electron[i] = -99.0;
    pfp_longlikelihood_photon[i] = -99.0;
    for (size_t j = 0; j<3; ++j){
      pfpvertex_recon[i][j] = -99999.0;
      pfpvertex_truth[i][j] = -99999.0;
      pfpend_truth[i][j] = -99999.0;
      pfpdir_recon[i][j] = -99999.0;
      pfpdir_truth[i][j] = -99999.0;
      TotalCharge[i][j] = 0;
      pfpcharge[i][j] = 0;
      pfpenergy_recon[i][j] = 0;
      pfp_datahits[i][j] = 1;
      pfp_purity[i][j] = -1.0;
      pfp_completeness[i][j] = -1.0;
      for (size_t k = 0; k<100; ++k) LongProf[i][j][k] = 0;
      for (size_t k = 0; k<16; ++k) {
        TranProf[i][j][k] = 0;
        TranProf_1[i][j][k] = 0;
        TranProf_2[i][j][k] = 0;
        TranProf_3[i][j][k] = 0;
        TranProf_4[i][j][k] = 0;
        TranProf_5[i][j][k] = 0;
      }
    }
  }
}

DEFINE_ART_MODULE(PFPLikelihood)

////////////////////////////////////////////////////////////////////////
// Class:       ShowerSelection
// Plugin Type: analyzer (art v3_02_06)
// File:        ShowerSelection_module.cc
//
// Generated at Thu Sep 19 14:37:23 2019 by Wanwei Wu using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"  // to calculate the time offset for x
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" // for SCE correction
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include <TROOT.h>
#include <TStyle.h>

using namespace std;


class ShowerSelection;


class ShowerSelection : public art::EDAnalyzer {
public:
  explicit ShowerSelection(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerSelection(ShowerSelection const&) = delete;
  ShowerSelection(ShowerSelection&&) = delete;
  ShowerSelection& operator=(ShowerSelection const&) = delete;
  ShowerSelection& operator=(ShowerSelection&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // user's defined functions
  void beginJob() override;
  
  void makeDataProducts();
  bool insideFV(double vertex[3]);
  void saveShowerInformation(art::Ptr<recob::Shower> myShower, art::FindManyP<recob::Hit> myfmshhits, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> myparticles_per_hit, std::vector<std::unordered_map<int, double>> myallhit_trkide, art::Event const& e);
  void reset();

private:

  // Declare member data here.

  // read from fcl: Input tag
  std::string fMCTruthModuleLabel;
  std::string fHitModuleLabel;
  std::string fSliceModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fPFParticleSliceAssnsModuleLabel;
  std::string fPFParticleShowerAssnsModuleLabel;
  std::string fShowerModuleLabel;
  std::string fHitMatchingDataModuleLabel;

  // read from fcl: user's defined parameters
  float fFidVolCutX; // Fiducial Volume cut [cm]
  float fFidVolCutY;
  float fFidVolCutZ;
  bool fSaveOnlyLargestShowerPerSlice; // true: only save the largest shower in a slice; false: save all showers in a slice 


  // Fiducial Volume parameters
  float fFidVolXmin;
  float fFidVolXmax;
  float fFidVolYmin;
  float fFidVolYmax;
  float fFidVolZmin;
  float fFidVolZmax;
  
  // art service handle
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // TTree
  TTree *fEventTree;
  
  // event information
  int event;
  int run;
  int subrun;

  // MCTruth
  int nuPDG_truth;
  int leptonPDG_truth;
  double leptonEnergy_truth;
  int ccnc_truth;
  double nuVertex_truth[4];
  double nuMomentum_truth[4]; // neutrino incoming momentum
  double nuEnergy_truth; // neutrino incoming energy
  int MC_lepton_ID; // here, electron
  std::vector<double> ehit_totalEnergy; // [MeV]

  // slice information
  std::vector<int> slice_id;

  // shower information
  int numOfShowers;
  std::vector<int> shower_id;
  std::vector<double> sh_start_X;
  std::vector<double> sh_start_Y;
  std::vector<double> sh_start_Z;
  std::vector<double> sh_direction_X;
  std::vector<double> sh_direction_Y;
  std::vector<double> sh_direction_Z;
  std::vector<std::vector<int>> sh_hasPrimary_e; // whether a shower is from primary electron
  std::vector<std::vector<double>> sh_ehit_energy; // for 3 planes
  std::vector<std::vector<double>> sh_allhit_energy;
  std::vector<std::vector<int>> showerParticlePDG;
  std::vector<std::vector<double>> showerEnergyMC_Rel; // relativistic energy
  std::vector<std::vector<double>> showerEnergyMC;
  std::vector<std::vector<double>> showerEnergyReco; // use reco energy from plane 2
  std::vector<std::vector<double>> sh_start_X_MC;
  std::vector<std::vector<double>> sh_start_Y_MC;
  std::vector<std::vector<double>> sh_start_Z_MC;
  std::vector<std::vector<double>> sh_direction_X_MC;
  std::vector<std::vector<double>> sh_direction_Y_MC;
  std::vector<std::vector<double>> sh_direction_Z_MC;
  std::vector<std::vector<double>> sh_purity;
  std::vector<std::vector<double>> sh_completeness;

  // test information
  int test_count_slice;
  int test_count_pfp;
  int test_count_shower;
  //int test_count_shower_not_insideFV;
  int test_count_hit;
};


ShowerSelection::ShowerSelection(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fMCTruthModuleLabel = p.get<std::string>("MCTruthModuleLabel");
  fHitModuleLabel     = p.get<std::string>("HitModuleLabel");
  fSliceModuleLabel     = p.get<std::string>("SliceModuleLabel");
  fPFParticleModuleLabel     = p.get<std::string>("PFParticleModuleLabel");
  fPFParticleSliceAssnsModuleLabel = p.get<std::string>("PFParticleSliceAssnsModuleLabel");
  fPFParticleShowerAssnsModuleLabel = p.get<std::string>("PFParticleShowerAssnsModuleLabel");
  fShowerModuleLabel  = p.get<std::string>("ShowerModuleLabel");
  fHitMatchingDataModuleLabel = p.get<std::string>("HitMatchingDataModuleLabel");
  fFidVolCutX         = p.get<float>("FidVolCutX");
  fFidVolCutY         = p.get<float>("FidVolCutY");
  fFidVolCutZ         = p.get<float>("FidVolCutZ");
  fSaveOnlyLargestShowerPerSlice  = p.get<bool>("SaveOnlyLargestShowerPerSlice");

  this->makeDataProducts();
}

void ShowerSelection::makeDataProducts() {
  // fEventTree
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
  fEventTree->Branch("ehit_totalEnergy", &ehit_totalEnergy);

  fEventTree->Branch("slice_id", &slice_id);

  fEventTree->Branch("numOfShowers", &numOfShowers, "numOfShowers/I");
  fEventTree->Branch("shower_id", &shower_id);
  fEventTree->Branch("sh_hasPrimary_e", &sh_hasPrimary_e);
  fEventTree->Branch("sh_ehit_energy", &sh_ehit_energy);
  fEventTree->Branch("sh_allhit_energy", &sh_allhit_energy);
  fEventTree->Branch("showerParticlePDG", &showerParticlePDG);
  fEventTree->Branch("showerEnergyMC_Rel", &showerEnergyMC_Rel);
  fEventTree->Branch("showerEnergyMC", &showerEnergyMC);
  fEventTree->Branch("showerEnergyReco", &showerEnergyReco);
  fEventTree->Branch("sh_start_X", &sh_start_X);
  fEventTree->Branch("sh_start_Y", &sh_start_Y);
  fEventTree->Branch("sh_start_Z", &sh_start_Z);
  fEventTree->Branch("sh_direction_X", &sh_direction_X);
  fEventTree->Branch("sh_direction_Y", &sh_direction_Y);
  fEventTree->Branch("sh_direction_Z", &sh_direction_Z);
  fEventTree->Branch("sh_start_X_MC", &sh_start_X_MC);
  fEventTree->Branch("sh_start_Y_MC", &sh_start_Y_MC);
  fEventTree->Branch("sh_start_Z_MC", &sh_start_Z_MC);
  fEventTree->Branch("sh_direction_X_MC", &sh_direction_X_MC);
  fEventTree->Branch("sh_direction_Y_MC", &sh_direction_Y_MC);
  fEventTree->Branch("sh_direction_Z_MC", &sh_direction_Z_MC);
  fEventTree->Branch("sh_purity", &sh_purity);
  fEventTree->Branch("sh_completeness", &sh_completeness);

}

void ShowerSelection::beginJob() {
  mf::LogInfo("ShowerSelection") << ".... begin of job ...." << endl;
  // Get geometry: Fiducial Volume
  auto const* geo = lar::providerFrom<geo::Geometry>();
  double minX = 1e9; // [cm]
  double maxX = -1e9;
  double minY = 1e9;
  double maxY = -1e9;
  double minZ = 1e9;
  double maxZ = -1e9;

  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);

    if (minX > world[0] - geo->DetHalfWidth(i))  minX = world[0]-geo->DetHalfWidth(i);
    if (maxX < world[0] + geo->DetHalfWidth(i))  maxX = world[0]+geo->DetHalfWidth(i);
    if (minY > world[1] - geo->DetHalfHeight(i)) minY = world[1]-geo->DetHalfHeight(i);
    if (maxY < world[1] + geo->DetHalfHeight(i)) maxY = world[1]+geo->DetHalfHeight(i);
    if (minZ > world[2] - geo->DetLength(i)/2.)  minZ = world[2]-geo->DetLength(i)/2.;
    if (maxZ < world[2] + geo->DetLength(i)/2.)  maxZ = world[2]+geo->DetLength(i)/2.;
  }

  fFidVolXmin = minX + fFidVolCutX;
  fFidVolXmax = maxX - fFidVolCutX;
  fFidVolYmin = minY + fFidVolCutY;
  fFidVolYmax = maxY - fFidVolCutY;
  fFidVolZmin = minZ + fFidVolCutZ;
  fFidVolZmax = maxZ - 2 * fFidVolCutZ; // for downstream z, we use a cut of 60 cm away from the edge 

  cout << "\n....Fiducial Volume (length unit: cm):\n"
    << "\t" << fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
    << "\t" << fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
    << "\t" << fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";


}

bool ShowerSelection::insideFV( double vertex[3] ) {
  double x = vertex[0];
  double y = vertex[1];
  double z = vertex[2];

  if (x >= fFidVolXmin && x <= fFidVolXmax &&
      y >= fFidVolYmin && y <= fFidVolYmax &&
      z >= fFidVolZmin && z <= fFidVolZmax) {
    return true;
  }
  else {
    return false;
  }
}

void ShowerSelection::saveShowerInformation(art::Ptr<recob::Shower> myShower, art::FindManyP<recob::Hit> myfmshhits, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> myparticles_per_hit, std::vector<std::unordered_map<int, double>> myallhit_trkide, art::Event const& e) {

  // SCE correction for the position in order to compare true and reco
  // ideally we would like to correct the reco and then compare it to truth
  // However, the SCE correction implemented here is to correct the true position and then compare to reco.
  
  // Get space charge correction
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  // Get time offset for x space charge correction
  auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const& mct_h = e.getValidHandle< std::vector<simb::MCTruth> >("generator");
  auto gen = mct_h->at(0);
  double g4Ticks = detClocks->TPCG4Time2Tick(gen.GetNeutrino().Nu().T()) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
  double xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);


  cout << "myShower->ID(): " << myShower->ID() << endl;
  cout << "myShower->Energy().at(2): " << myShower->Energy().at(2) << endl;
  cout << "myShower.key(): " << myShower.key() << endl;

  shower_id.push_back(myShower.key());
  sh_start_X.push_back(myShower->ShowerStart().X());
  sh_start_Y.push_back(myShower->ShowerStart().Y());
  sh_start_Z.push_back(myShower->ShowerStart().Z());
  sh_direction_X.push_back(myShower->Direction().X()); // Direction(): direction cosines at the start of the shower
  sh_direction_Y.push_back(myShower->Direction().Y());
  sh_direction_Z.push_back(myShower->Direction().Z());

  // shower hits
  std::vector<art::Ptr<recob::Hit>> sh_hits;
  sh_hits = myfmshhits.at(myShower.key());

  std::vector<simb::MCParticle const *> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const *> match_vec;
  std::vector<std::unordered_map<int, double>> trkide(3);
  std::unordered_map<int, double> shhit_trkide;

  for (size_t ihit = 0; ihit < sh_hits.size(); ++ihit) {
    test_count_hit += 1;

    art::Ptr<recob::Hit> hit = sh_hits[ihit];
    particle_vec.clear();
    match_vec.clear();
    myparticles_per_hit.get(hit.key(), particle_vec, match_vec);
    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
      if (!(pi_serv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId()))) {
        cout << "........www...................." << endl;
        cout << match_vec[i_p]->energy;
      }
      // here, from the trackID of the MCParticle, we can find the EveTrackId (the main mother at very beginning. We should use this to classify a shower)
      trkide[hit->WireID().Plane][pi_serv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId())] += match_vec[i_p]->energy; // match_vec[i_p]->numElectrons;
      shhit_trkide[pi_serv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId())] += match_vec[i_p]->energy; // match_vec[i_p]->numElectrons;
    }
  } // end of loop ihit

  int shhit_trackID = -999;
  double shhit_maxenergy = -99.0;
  //double shhit_partial_E = 0.0; // shower energy from the "shower particle"
  //double shhit_total_E = 0.0; // total energy of the shower
  for (auto const & p_id : shhit_trkide) {
    //shhit_total_E += p_id.second;
    if (p_id.second > shhit_maxenergy) {
      //shhit_partial_E = p_id.second;
      shhit_trackID = p_id.first;
      shhit_maxenergy = p_id.second;
    }
  }

  // to find the shower particle, we usually use hits on all planes, i.e., using the folloing:
  if (pi_serv->TrackIdToMotherParticle_P(shhit_trackID)) {
    cout << "pi_serv->TrackIdToMotherParticle_P(shhit_trackID)->PdgCode(): " << pi_serv->TrackIdToMotherParticle_P(shhit_trackID)->PdgCode() << endl;
  }

  std::vector<double> sh_vertexMC_x(3);
  std::vector<double> sh_vertexMC_y(3);
  std::vector<double> sh_vertexMC_z(3);
  std::vector<double> sh_momentumMC_x(3); // start momentum
  std::vector<double> sh_momentumMC_y(3); // start momentum
  std::vector<double> sh_momentumMC_z(3); // start momentum
  std::vector<int> sh_particlePDG(3);
  std::vector<int> sh_withPrimary_e(3);
  std::vector<double> sh_energyMC_Rel(3);
  std::vector<double> sh_energyMC(3);
  std::vector<double> sh_energyReco(3);
  std::vector<double> shower_purity(3);
  std::vector<double> shower_completeness(3);
  sh_vertexMC_x.clear();
  sh_vertexMC_y.clear();
  sh_vertexMC_z.clear();
  sh_momentumMC_x.clear();
  sh_momentumMC_y.clear();
  sh_momentumMC_z.clear();
  sh_particlePDG.clear();
  sh_withPrimary_e.clear();
  sh_energyMC_Rel.clear();
  sh_energyMC.clear();
  sh_energyReco.clear();
  shower_purity.clear();
  shower_completeness.clear();

  for (int p = 0; p < 3; p++) {
    int sh_trackID = -999;
    double sh_partial_E = 0.0; // shower energy from the "shower particle"
    double sh_total_E = 0.0; // total energy of the shower
    double maxenergy = -99.0;
    // find shower particle, here, for validation we use hits on each plane
    for (auto const & p_id : trkide[p]) {
      sh_total_E += p_id.second;
      if (p_id.second > maxenergy) {
        sh_partial_E = p_id.second;
        sh_trackID = p_id.first;
        maxenergy = p_id.second;
      }
    }

    // if (allhit_trkide[p][sh_trackID] <= 0.0 || sh_total_E <= 0.0) continue;

    if (pi_serv->TrackIdToMotherParticle_P(sh_trackID)) {
      cout << "pi_serv->TrackIdToMotherParticle_P(sh_trackID)->PdgCode(): " << pi_serv->TrackIdToMotherParticle_P(sh_trackID)->PdgCode() << endl;
      cout << "pi_serv->TrackIdToMotherParticle_P(sh_trackID)->Mother(): " << pi_serv->TrackIdToMotherParticle_P(sh_trackID)->Mother() << endl;

      const simb::MCParticle* sh_particle;
      sh_particle = pi_serv->TrackIdToParticle_P(sh_trackID);
      cout << "sh_particle->PdgCode(): " << sh_particle->PdgCode() << endl;
      cout << "sh_particle->Mother(): " << sh_particle->Mother() << endl;

      // apply SCE corrections to MC vertex
      auto sce_offset = SCE->GetPosOffsets(geo::Point_t(sh_particle->Vx(),sh_particle->Vy(),sh_particle->Vx()));

      sh_vertexMC_x.push_back(sh_particle->Vx() - sce_offset.X() + xtimeoffset + 0.6);
      // The 0.6 cm term accounts for the different coordinate center between WC (Y wires) and LArSoft (U wires)
      sh_vertexMC_y.push_back(sh_particle->Vy() + sce_offset.Y());
      sh_vertexMC_z.push_back(sh_particle->Vz() + sce_offset.Z());
      
      sh_momentumMC_x.push_back(sh_particle->Px() * 1000.0); // [MeV]
      sh_momentumMC_y.push_back(sh_particle->Py() * 1000.0);
      sh_momentumMC_z.push_back(sh_particle->Pz() * 1000.0);

      sh_particlePDG.push_back(sh_particle->PdgCode());
      sh_energyMC_Rel.push_back(sh_particle->E() * 1000.0); // [MeV]
      sh_energyMC.push_back((sh_particle->E() - sh_particle->Mass())*1000.0); // kinetic energy [MeV]

      double particle_total_E = myallhit_trkide[p][sh_trackID]; // all energy deposits of "shower particle" in the event

      if (sh_total_E > 0) {
        shower_purity.push_back(sh_partial_E / sh_total_E);
      }
      else {
        shower_purity.push_back(0.0);
      }
      if (particle_total_E > 0) {
        shower_completeness.push_back(sh_partial_E / particle_total_E);
      }
      else {
        shower_completeness.push_back(0.0);
      }
    }
    else {
      cout << ".................wwww..............." << endl;
      sh_vertexMC_x.push_back(-99999.0);
      sh_vertexMC_y.push_back(-99999.0);
      sh_vertexMC_z.push_back(-99999.0);
      sh_momentumMC_x.push_back(-99999.0);
      sh_momentumMC_y.push_back(-99999.0);
      sh_momentumMC_z.push_back(-99999.0);
      sh_particlePDG.push_back(-99999);
      sh_energyMC_Rel.push_back(-99999.0);
      sh_energyMC.push_back(-99999.0);
      shower_purity.push_back(-99999.0);
      shower_completeness.push_back(-99999.0);
    }

    sh_energyReco.push_back(myShower->Energy().at(p));

    // primary electron shower
    if (sh_trackID == MC_lepton_ID) {
      sh_withPrimary_e.push_back(1);
    }
    else {
      sh_withPrimary_e.push_back(0);
    }

  } // end of loop p

  sh_hasPrimary_e.push_back(sh_withPrimary_e);
  sh_start_X_MC.push_back(sh_vertexMC_x);
  sh_start_Y_MC.push_back(sh_vertexMC_y);
  sh_start_Z_MC.push_back(sh_vertexMC_z);
  sh_direction_X_MC.push_back(sh_momentumMC_x);
  sh_direction_Y_MC.push_back(sh_momentumMC_y);
  sh_direction_Z_MC.push_back(sh_momentumMC_z);
  showerParticlePDG.push_back(sh_particlePDG);
  showerEnergyMC_Rel.push_back(sh_energyMC_Rel);
  showerEnergyMC.push_back(sh_energyMC);
  showerEnergyReco.push_back(sh_energyReco);
  sh_purity.push_back(shower_purity);
  sh_completeness.push_back(shower_completeness);
}

void ShowerSelection::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  
  // reset some variables
  reset();

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  
  cout << "e.isRealData(): " << e.isRealData() << endl;

  // art handle: mctruth
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector< art::Ptr<simb::MCTruth> > mctruthlist;
  if (e.getByLabel(fMCTruthModuleLabel, mctruthListHandle)) {
    art::fill_ptr_vector(mctruthlist, mctruthListHandle);
  }

  // art handle: hit list
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector< art::Ptr<recob::Hit> > hitlist;
  if(e.getByLabel(fHitModuleLabel,hitListHandle)){
    art::fill_ptr_vector(hitlist, hitListHandle);
  }

  // art handle: slice list
  art::Handle< std::vector<recob::Slice> > sliceListHandle;
  std::vector<art::Ptr<recob::Slice>> slicelist;
  if ( e.getByLabel(fSliceModuleLabel, sliceListHandle) ) {
    art::fill_ptr_vector(slicelist, sliceListHandle);
  }

  // art handle: pfparticle list
  art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfplist;
  if ( e.getByLabel(fPFParticleModuleLabel, pfpListHandle) ) {
    art::fill_ptr_vector(pfplist, pfpListHandle);
  }

  // art handle: shower list
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower>> showerlist;
  if ( e.getByLabel(fShowerModuleLabel, showerListHandle) ) {
    art::fill_ptr_vector(showerlist, showerListHandle);
  }

  art::FindManyP<recob::PFParticle> fm_pfp_slice_assns(sliceListHandle, e, fPFParticleSliceAssnsModuleLabel); // PFParticle associated with Slice

  art::FindManyP<recob::Shower> fm_pfp_shower_assns(pfpListHandle, e, fPFParticleShowerAssnsModuleLabel);
  
  art::FindManyP<recob::Hit> fmshhits (showerListHandle, e, fShowerModuleLabel);

  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitListHandle, e, fHitMatchingDataModuleLabel);

  // ------- MCTruch: Generator -------------
  if ( mctruthlist.size() > 1 ) cout << ">>  mctruthlist.size(): " << mctruthlist.size() << endl;
  // ideally one event is corresponding to one mctruth
  for (size_t i = 0; i < mctruthlist.size(); ++i) {
    art::Ptr<simb::MCTruth> mctruth = mctruthlist[i];
    if ( mctruth->NeutrinoSet() ) { // NeutrinoSet(): whether the neutrino information has been set
      simb::MCNeutrino nu = mctruth->GetNeutrino();// GetNeutrino(): reference to neutrino info
      simb::MCParticle neutrino = nu.Nu(); // Nu(): the incoming neutrino
      simb::MCParticle lepton = nu.Lepton();// Lepton(): the outgoing lepton

      nuPDG_truth = neutrino.PdgCode();
      leptonPDG_truth = lepton.PdgCode();
      leptonEnergy_truth = lepton.E();
      ccnc_truth = nu.CCNC(); // CCNC = 0 stands for CC; CCNC = 1 strands for NC

      const TLorentzVector& vertex =neutrino.Position(0);
      vertex.GetXYZT(nuVertex_truth); // for vertex, nuVertex_truth[3] is set to 0

      const TLorentzVector& nu_momentum = nu.Nu().Momentum(0);
      nu_momentum.GetXYZT(nuMomentum_truth);
      nuEnergy_truth = nuMomentum_truth[3];

      // check if inside Fiducial Volume
      double check_vertex[3];
      check_vertex[0] = nuVertex_truth[0];
      check_vertex[1] = nuVertex_truth[1];
      check_vertex[2] = nuVertex_truth[2];
      if (!insideFV(check_vertex)) {
      // check if inside Fiducial Volume, now switch to use reconstructed vertex
      //if (!insideFV(nuVertex_truth)) {
        cout << "\n **** Interaction is NOT inside the Fiducial Volume. RETURN ****" << endl;
        return;
      }
    }
  } // end of loop mctruthlist i

  //------------ MCTruth: Geant4 -------------------
  //  // Note: generator level MCPartilceList is different from Geant4 level MCParticleList.  MCParticleList(Geant4) contains all particles in MCParticleList(generator) but their the indexes (TrackIds) are different.
  simb::MCParticle *MC_lepton = NULL; // Geant4 level
  const sim::ParticleList& plist = pi_serv->ParticleList();
  // note here plist is empty for the overlay sample.
  cout << "plist.size():" << plist.size() << endl;
  for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar) {
    simb::MCParticle *particle = 0;
    particle = ipar->second; // first is index(TrackId), second is value (point address)
    auto & truth = pi_serv->ParticleToMCTruth_P(particle);
    // beam neutrino only, only look at electron
    if ( truth->Origin()==simb::kBeamNeutrino && std::abs(particle->PdgCode()) == 11 &&  particle->Mother()==0){ // primary lepton; Mother() = 0 means e^{-} for v+n=e^{-}+p
      MC_lepton_ID = particle->TrackId();
      MC_lepton = particle;
    }
  } // end of loop plist ipar

  cout << "MC_lepton_ID: " << MC_lepton_ID << endl;
  cout << "MC_lepton->TrackId(): " << MC_lepton->TrackId() << endl;
  cout << "leptonEnergy_truth: " << leptonEnergy_truth << endl;

  // -------- hit information --------------
  //   //find all hits from the primary electron; build a map for total hit charges corresponding to MC_lepton_ID (Geant4)

  std::vector<simb::MCParticle const *> allhit_particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const *> allhit_match_vec;
  std::vector<std::unordered_map<int, double>> allhit_trkide(3);
  //allhit_trkide.clear();

  for (size_t i = 0; i < hitlist.size(); ++i) {
    art::Ptr<recob::Hit> hit = hitlist[i];
    allhit_particle_vec.clear();
    allhit_match_vec.clear();
    particles_per_hit.get(hit.key(), allhit_particle_vec, allhit_match_vec);
    for (size_t i_p = 0; i_p < allhit_particle_vec.size(); ++i_p) {
      allhit_trkide[hit->WireID().Plane][pi_serv->TrackIdToEveTrackId(allhit_particle_vec.at(i_p)->TrackId())] += allhit_match_vec[i_p]->energy; // match_vec[i_p]->numElectrons;
    }
  } // end of loop hitlist i

  for (int p = 0; p < 3; p++) {
    ehit_totalEnergy.push_back(allhit_trkide[p][MC_lepton_ID]);
    cout << "allhit_trkide[p][MC_lepton_ID]: " << allhit_trkide[p][MC_lepton_ID] << endl;
  }

  // --------- shower information -----------
  cout << "showerlist.size(): " << showerlist.size() << endl;
  numOfShowers = showerlist.size();
  // clear the vectors
  slice_id.clear();
  shower_id.clear();
  sh_hasPrimary_e.clear();
  sh_ehit_energy.clear();
  sh_allhit_energy.clear();
  showerParticlePDG.clear();
  showerEnergyMC_Rel.clear();
  showerEnergyMC.clear();
  showerEnergyReco.clear();
  sh_start_X.clear();
  sh_start_Y.clear();
  sh_start_Z.clear();
  sh_direction_X.clear();
  sh_direction_Y.clear();
  sh_direction_Z.clear();
  sh_start_X_MC.clear();
  sh_start_Y_MC.clear();
  sh_start_Z_MC.clear();
  sh_direction_X_MC.clear();
  sh_direction_Y_MC.clear();
  sh_direction_Z_MC.clear();
  sh_purity.clear();
  sh_completeness.clear();

  // ------- shower associated with slice ----------
  cout << "slicelist.size(): " << slicelist.size() << endl;

/*
  std::vector<art::Ptr<recob::Shower>> showerVec;
  int test_count_2 = 0;
  cout << "pfplist.size(): " << pfplist.size() << endl;
  for (size_t j = 0; j < pfplist.size(); ++j) {
    if (fm_pfp_shower_assns.isValid()) {
      std::vector<art::Ptr<recob::Shower>> showers = fm_pfp_shower_assns.at(j);
      for (size_t k = 0; k < showers.size(); ++k) {
        test_count_2 += 1;
        cout << showers[k]->ID() << endl;
        showerVec.push_back(showers[k]);
      }
    }

  } // end of pfs.size() j

  cout << "test_count_2: " << test_count_2 << endl;

  for (size_t i = 0; i < showerVec.size(); ++i){
    cout << showerVec[i]->ID() << endl;
  }
*/

  
  //std::vector<double> slice_maxenergyShower_energy;
  //std::vector<int> slice_maxenergyShower_ID;
  //slice_maxenergyShower_energy.clear();
  //slice_maxenergyShower_ID.clear();
  
  for (size_t i = 0; i < slicelist.size(); ++i) {
    cout << "slicelist[i].key(): " << slicelist[i].key() << endl;

    test_count_slice += 1;
    if (fm_pfp_slice_assns.isValid()) {
      std::vector<art::Ptr<recob::PFParticle>> pfs = fm_pfp_slice_assns.at(i);
      cout << "slice " << i << ": " << endl;
      
      double maxenergy_sh_energy = -999.0;
      int maxenergy_sh_ID = -999;
      art::Ptr<recob::Shower> largestShower;

      for (size_t j = 0; j < pfs.size(); ++j) {
        test_count_pfp += 1;
        if (fm_pfp_shower_assns.isValid()) {
          std::vector<art::Ptr<recob::Shower>> showers = fm_pfp_shower_assns.at(pfs[j].key()); // here, should use pfs[j].key() to get the right showers.
          for (size_t k = 0; k < showers.size(); ++k) {
            test_count_shower += 1;
            //cout << "shower ID: " <<  showers[k]->ID() << endl;
            //cout << showers[k].key() << endl;

          //  // check if inside Fiducial Volume, now switch to use reconstructed vertex
          //  double check_shower_vertex[3];
          //  check_shower_vertex[0] = showers[k]->ShowerStart().X();
          //  check_shower_vertex[1] = showers[k]->ShowerStart().Y();
          //  check_shower_vertex[2] = showers[k]->ShowerStart().Z();
          //  if (!insideFV(check_shower_vertex)) {
          //    test_count_shower_not_insideFV += 1;
          //    cout << "\n **** shower is NOT inside the Fiducial Volume. RETURN ****" << endl;
          //    continue; // here, should use continue instead of return; we want to check all the showers
          //  }

            double sh_energy = showers[k]->Energy().at(2); 
            if (sh_energy < 0) {
              cout << ".....www: shower energy less than 0 ---->sh_energy: " << sh_energy << endl;
            }
            
            if (sh_energy > maxenergy_sh_energy) {
              maxenergy_sh_energy = sh_energy;
              maxenergy_sh_ID = showers[k]->ID();
              largestShower = showers[k];
            }
          
            // save all showers in a slice
            if (!fSaveOnlyLargestShowerPerSlice) {
              slice_id.push_back(slicelist[i].key());
              saveShowerInformation(showers[k], fmshhits, particles_per_hit, allhit_trkide, e);
            }

          } // end of for shower.size() k
        } // end of if fm_pfp_shower_assns.isValid()

      } // end of pfs.size() j
      
      //slice_maxenergyShower_energy.push_back(maxenergy_sh_energy);
      //slice_maxenergyShower_ID.push_back(maxenergy_sh_ID);
      if (fSaveOnlyLargestShowerPerSlice && largestShower) {
        cout << "maxenergy_sh_energy: " << maxenergy_sh_energy << endl;
        cout << "maxenergy_sh_ID: " << maxenergy_sh_ID << endl;
        slice_id.push_back(slicelist[i].key());
        saveShowerInformation(largestShower, fmshhits, particles_per_hit, allhit_trkide, e);
      }
    }  // end of if fm_pfp_slice_assns.isValid()
  } // end of loop slicelist.size() i

  cout << "test_count_slice: " << test_count_slice << endl;
  cout << "test_count_pfp: " << test_count_pfp << endl;
  cout << "test_count_shower: " << test_count_shower << endl;
  //cout << "test_count_shower_not_insideFV: " << test_count_shower_not_insideFV << endl;
  cout << "test_count_hit: " << test_count_hit << endl;
 
 //if (test_count_shower_not_insideFV == test_count_shower) {
 //   return;
 //}

  fEventTree->Fill();

}

void ShowerSelection::reset() {
  run = -999;
  subrun = -999;
  event = -999;

  test_count_slice = 0;
  test_count_pfp = 0;
  test_count_shower = 0;
  //test_count_shower_not_insideFV = 0;
  test_count_hit = 0;
}


DEFINE_ART_MODULE(ShowerSelection)

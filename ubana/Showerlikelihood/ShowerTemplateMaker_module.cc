////////////////////////////////////////////////////////////////////////
// Class:       ShowerTemplateMaker
// Plugin Type: analyzer (art v3_02_06)
// File:        ShowerTemplateMaker_module.cc
//
// Generated at Mon Oct  7 14:52:37 2019 by Wanwei Wu using cetskelgen
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
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TROOT.h>
#include <TStyle.h>

using namespace std;


class ShowerTemplateMaker;


class ShowerTemplateMaker : public art::EDAnalyzer {
  public:
    explicit ShowerTemplateMaker(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ShowerTemplateMaker(ShowerTemplateMaker const&) = delete;
    ShowerTemplateMaker(ShowerTemplateMaker&&) = delete;
    ShowerTemplateMaker& operator=(ShowerTemplateMaker const&) = delete;
    ShowerTemplateMaker& operator=(ShowerTemplateMaker&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // user's defined functions
    void beginJob() override;

    void makeDataProducts();
    bool insideFV(double vertex[3]);
    void reset();
    void saveShowerInformation(art::Ptr<recob::Shower> myShower, art::FindManyP<recob::Hit> myfmshhits, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> myparticles_per_hit, std::vector<std::unordered_map<int, double>> myallhit_trkide, art::Event const& e);
    void showerProfile(std::vector< art::Ptr<recob::Hit> >shhits, TVector3 shvertex, TVector3 shdirection, double shlength, double shenergy, int shpdg); 
    void convertXYZtoGlobalWireTDC(TVector3 pointXYZ,  double& pointGlobalWire, double& pointGlobalTDC);
    void fillElectronHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5);
    void fillPhotonHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5);
    void fillProtonHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5);
    void fillPionHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5);
    void fillOthersHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5);

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
    
    // for MC, we may skip the lifetime correction; for data, we should use lifetime from database 
    calo::CalorimetryAlg fCalorimetryAlg;

    // read from fcl: user's defined parameters
    float fFidVolCutX; // Fiducial Volume cut [cm]
    float fFidVolCutY;
    float fFidVolCutZ;
    bool fFVCutOnRecon; // true: use FV cut on reconstructed vertex?; false: FV cut on true interaction vertex 
    double fModBoxA;
    double fModBoxB;
    std::vector<double> fADCtoE; // calibration constants
    double fRecombination;
    bool fSaveOnlyLargestShowerPerSlice; // true: only save the largest shower in a slice; false: save all showers in a slice 
    float fShowerCutPurity;
    float fShowerCutCompleteness;
    bool fShowerEnergyCorrection;
    std::vector<double> fShowerEnergyCorrectionSlope; // for 3 planes
    std::vector<double> fShowerEnergyCorrectionIntercept; // for 3 planes


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

    // slice information
    std::vector<int> slice_id;

    // shower information
    std::vector<int> shower_id;
    std::vector<double> sh_start_X; // [cm]
    std::vector<double> sh_start_Y; // [cm]
    std::vector<double> sh_start_Z; // [cm]
    std::vector<double> sh_direction_X; // Direction(): direction cosines at start of the shower
    std::vector<double> sh_direction_Y;
    std::vector<double> sh_direction_Z;
    std::vector<double> sh_length; // [cm]
    // shower particle: particle associated with the shower
    std::vector<int> showerParticlePDG;
    std::vector<int> sh_hasPrimary_e; // whether a shower is from primary electron
    std::vector<double> showerEnergyMC; // shower particle's kinetic energy [MeV]
    std::vector<double> sh_purity;
    std::vector<double> sh_completeness;
    std::vector<std::vector<double>> showerEnergyReco;

    std::vector<std::vector<double>> showerHitLongDist;
    std::vector<std::vector<double>> showerHitTranDist;
    std::vector<std::vector<double>> showerHitChargeQ;
    // convert length to radiation length of the material, for LAr X0 = 14.0 cm
    const double X0 = 14.0;

    // histogram information
    const int LBINS = 100;
    const double LMIN = 0.0;
    const double LMAX = 25.0;

    const int TBINS = 20;
    const double TMIN = -8.0;
    const double TMAX = 8.0;

    const int EBINS = 100;
    const double EMIN = 0.0; // [MeV]
    const double EMAX = 10000.0; // [GeV]

    const int QBINS = 1000;
    const double QMIN = 0.0;
    const double QMAX = 150000.0;

    // primary electron
    TH3F* fLongHist3D_electron;
    TH3F* fTranHist3D_electron;
    TH3F* fTranHist3D_electron_1;
    TH3F* fTranHist3D_electron_2;
    TH3F* fTranHist3D_electron_3;
    TH3F* fTranHist3D_electron_4;
    TH3F* fTranHist3D_electron_5;

    TProfile* fLongProfile_electron;
    TProfile* fTranProfile_electron;
    TProfile* fTranProfile_electron_1;
    TProfile* fTranProfile_electron_2;
    TProfile* fTranProfile_electron_3;
    TProfile* fTranProfile_electron_4;
    TProfile* fTranProfile_electron_5;

    // photon
    TH3F* fLongHist3D_photon;
    TH3F* fTranHist3D_photon;
    TH3F* fTranHist3D_photon_1;
    TH3F* fTranHist3D_photon_2;
    TH3F* fTranHist3D_photon_3;
    TH3F* fTranHist3D_photon_4;
    TH3F* fTranHist3D_photon_5;

    TProfile* fLongProfile_photon;
    TProfile* fTranProfile_photon;
    TProfile* fTranProfile_photon_1;
    TProfile* fTranProfile_photon_2;
    TProfile* fTranProfile_photon_3;
    TProfile* fTranProfile_photon_4;
    TProfile* fTranProfile_photon_5;

    // proton
    TH3F* fLongHist3D_proton;
    TH3F* fTranHist3D_proton;
    TH3F* fTranHist3D_proton_1;
    TH3F* fTranHist3D_proton_2;
    TH3F* fTranHist3D_proton_3;
    TH3F* fTranHist3D_proton_4;
    TH3F* fTranHist3D_proton_5;

    TProfile* fLongProfile_proton;
    TProfile* fTranProfile_proton;
    TProfile* fTranProfile_proton_1;
    TProfile* fTranProfile_proton_2;
    TProfile* fTranProfile_proton_3;
    TProfile* fTranProfile_proton_4;
    TProfile* fTranProfile_proton_5;

    // pion
    TH3F* fLongHist3D_pion;
    TH3F* fTranHist3D_pion;
    TH3F* fTranHist3D_pion_1;
    TH3F* fTranHist3D_pion_2;
    TH3F* fTranHist3D_pion_3;
    TH3F* fTranHist3D_pion_4;
    TH3F* fTranHist3D_pion_5;

    TProfile* fLongProfile_pion;
    TProfile* fTranProfile_pion;
    TProfile* fTranProfile_pion_1;
    TProfile* fTranProfile_pion_2;
    TProfile* fTranProfile_pion_3;
    TProfile* fTranProfile_pion_4;
    TProfile* fTranProfile_pion_5;

    // others
    TH3F* fLongHist3D_others;
    TH3F* fTranHist3D_others;
    TH3F* fTranHist3D_others_1;
    TH3F* fTranHist3D_others_2;
    TH3F* fTranHist3D_others_3;
    TH3F* fTranHist3D_others_4;
    TH3F* fTranHist3D_others_5;

    TProfile* fLongProfile_others;
    TProfile* fTranProfile_others;
    TProfile* fTranProfile_others_1;
    TProfile* fTranProfile_others_2;
    TProfile* fTranProfile_others_3;
    TProfile* fTranProfile_others_4;
    TProfile* fTranProfile_others_5;

};


ShowerTemplateMaker::ShowerTemplateMaker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
  // More initializers here.
  fCalorimetryAlg (p.get<fhicl::ParameterSet> ("CaloAlg"))
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
  fFVCutOnRecon       = p.get<bool>("FVCutOnRecon");
  fModBoxA       = p.get<double>("ModBoxA");
  fModBoxB       = p.get<double>("ModBoxB");
  fADCtoE = p.get<std::vector<double>>("ADCtoE");
  fRecombination = p.get<double>("Recombination");
  fSaveOnlyLargestShowerPerSlice  = p.get<bool>("SaveOnlyLargestShowerPerSlice");
  fShowerCutPurity    = p.get<float>("ShowerCutPurity");
  fShowerCutCompleteness = p.get<float>("ShowerCutCompleteness");
  fShowerEnergyCorrection  = p.get<bool>("ShowerEnergyCorrection");
  fShowerEnergyCorrectionSlope = p.get<std::vector<double>>("ShowerEnergyCorrectionSlope");
  fShowerEnergyCorrectionIntercept = p.get<std::vector<double>>("ShowerEnergyCorrectionIntercept");

  this->makeDataProducts();
}

void ShowerTemplateMaker::makeDataProducts() {
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

  fEventTree->Branch("slice_id", &slice_id);

  fEventTree->Branch("shower_id", &shower_id);
  fEventTree->Branch("sh_start_X", &sh_start_X);
  fEventTree->Branch("sh_start_Y", &sh_start_Y);
  fEventTree->Branch("sh_start_Z", &sh_start_Z);
  fEventTree->Branch("sh_direction_X", &sh_direction_X);
  fEventTree->Branch("sh_direction_Y", &sh_direction_Y);
  fEventTree->Branch("sh_direction_Z", &sh_direction_Z);
  fEventTree->Branch("sh_length", &sh_length);
  fEventTree->Branch("showerParticlePDG", &showerParticlePDG);
  fEventTree->Branch("sh_hasPrimary_e", &sh_hasPrimary_e);
  fEventTree->Branch("showerEnergyMC", &showerEnergyMC);
  fEventTree->Branch("sh_purity", &sh_purity);
  fEventTree->Branch("sh_completeness", &sh_completeness);
  fEventTree->Branch("showerEnergyReco", &showerEnergyReco);

  fEventTree->Branch("showerHitLongDist", &showerHitLongDist);
  fEventTree->Branch("showerHitTranDist", &showerHitTranDist);
  fEventTree->Branch("showerHitChargeQ", &showerHitChargeQ);


  // primary electron
  fLongHist3D_electron = tfs->make<TH3F>("fLongHist3D_electron","fLongHist3D_electron; energy (GeV);t (radiation length 14 cm);Q", EBINS, EMIN, EMAX, LBINS, LMIN, LMAX, QBINS, QMIN, QMAX);
  fTranHist3D_electron = tfs->make<TH3F>("fTranHist3D_electron", "fTranHist3D_electron;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_electron_1 = tfs->make<TH3F>("fTranHist3D_electron_1", "fTranHist3D_electron_1;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_electron_2 = tfs->make<TH3F>("fTranHist3D_electron_2", "fTranHist3D_electron_2;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_electron_3 = tfs->make<TH3F>("fTranHist3D_electron_3", "fTranHist3D_electron_3;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_electron_4 = tfs->make<TH3F>("fTranHist3D_electron_4", "fTranHist3D_electron_4;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_electron_5 = tfs->make<TH3F>("fTranHist3D_electron_5", "fTranHist3D_electron_5;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);

  fLongProfile_electron = tfs->make<TProfile>("fLongProfile_electron", "fLongProfile_electron;t;Q", LBINS, LMIN, LMAX);
  fTranProfile_electron = tfs->make<TProfile>("fTranProfile_electron", "fTranProfile_electron;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_electron_1 = tfs->make<TProfile>("fTranProfile_electron_1", "fTranProfile_electron_1;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_electron_2 = tfs->make<TProfile>("fTranProfile_electron_2", "fTranProfile_electron_2;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_electron_3 = tfs->make<TProfile>("fTranProfile_electron_3", "fTranProfile_electron_3;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_electron_4 = tfs->make<TProfile>("fTranProfile_electron_4", "fTranProfile_electron_4;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_electron_5 = tfs->make<TProfile>("fTranProfile_electron_5", "fTranProfile_electron_5;dist;Q", TBINS, TMIN, TMAX);

  // photon
  fLongHist3D_photon = tfs->make<TH3F>("fLongHist3D_photon","fLongHist3D_photon; energy (GeV);t (radiation length 14 cm);Q", EBINS, EMIN, EMAX, LBINS, LMIN, LMAX, QBINS, QMIN, QMAX);
  fTranHist3D_photon = tfs->make<TH3F>("fTranHist3D_photon", "fTranHist3D_photon;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_photon_1 = tfs->make<TH3F>("fTranHist3D_photon_1", "fTranHist3D_photon_1;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_photon_2 = tfs->make<TH3F>("fTranHist3D_photon_2", "fTranHist3D_photon_2;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_photon_3 = tfs->make<TH3F>("fTranHist3D_photon_3", "fTranHist3D_photon_3;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_photon_4 = tfs->make<TH3F>("fTranHist3D_photon_4", "fTranHist3D_photon_4;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_photon_5 = tfs->make<TH3F>("fTranHist3D_photon_5", "fTranHist3D_photon_5;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);

  fLongProfile_photon = tfs->make<TProfile>("fLongProfile_photon", "fLongProfile_photon;t;Q", LBINS, LMIN, LMAX);
  fTranProfile_photon = tfs->make<TProfile>("fTranProfile_photon", "fTranProfile_photon;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_photon_1 = tfs->make<TProfile>("fTranProfile_photon_1", "fTranProfile_photon_1;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_photon_2 = tfs->make<TProfile>("fTranProfile_photon_2", "fTranProfile_photon_2;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_photon_3 = tfs->make<TProfile>("fTranProfile_photon_3", "fTranProfile_photon_3;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_photon_4 = tfs->make<TProfile>("fTranProfile_photon_4", "fTranProfile_photon_4;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_photon_5 = tfs->make<TProfile>("fTranProfile_photon_5", "fTranProfile_photon_5;dist;Q", TBINS, TMIN, TMAX);

  // proton
  fLongHist3D_proton = tfs->make<TH3F>("fLongHist3D_proton","fLongHist3D_proton; energy (GeV);t (radiation length 14 cm);Q", EBINS, EMIN, EMAX, LBINS, LMIN, LMAX, QBINS, QMIN, QMAX);
  fTranHist3D_proton = tfs->make<TH3F>("fTranHist3D_proton", "fTranHist3D_proton;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_proton_1 = tfs->make<TH3F>("fTranHist3D_proton_1", "fTranHist3D_proton_1;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_proton_2 = tfs->make<TH3F>("fTranHist3D_proton_2", "fTranHist3D_proton_2;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_proton_3 = tfs->make<TH3F>("fTranHist3D_proton_3", "fTranHist3D_proton_3;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_proton_4 = tfs->make<TH3F>("fTranHist3D_proton_4", "fTranHist3D_proton_4;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_proton_5 = tfs->make<TH3F>("fTranHist3D_proton_5", "fTranHist3D_proton_5;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);

  fLongProfile_proton = tfs->make<TProfile>("fLongProfile_proton", "fLongProfile_proton;t;Q", LBINS, LMIN, LMAX);
  fTranProfile_proton = tfs->make<TProfile>("fTranProfile_proton", "fTranProfile_proton;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_proton_1 = tfs->make<TProfile>("fTranProfile_proton_1", "fTranProfile_proton_1;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_proton_2 = tfs->make<TProfile>("fTranProfile_proton_2", "fTranProfile_proton_2;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_proton_3 = tfs->make<TProfile>("fTranProfile_proton_3", "fTranProfile_proton_3;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_proton_4 = tfs->make<TProfile>("fTranProfile_proton_4", "fTranProfile_proton_4;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_proton_5 = tfs->make<TProfile>("fTranProfile_proton_5", "fTranProfile_proton_5;dist;Q", TBINS, TMIN, TMAX);

  // pion
  fLongHist3D_pion = tfs->make<TH3F>("fLongHist3D_pion","fLongHist3D_pion; energy (GeV);t (radiation length 14 cm);Q", EBINS, EMIN, EMAX, LBINS, LMIN, LMAX, QBINS, QMIN, QMAX);
  fTranHist3D_pion = tfs->make<TH3F>("fTranHist3D_pion", "fTranHist3D_pion;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_pion_1 = tfs->make<TH3F>("fTranHist3D_pion_1", "fTranHist3D_pion_1;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_pion_2 = tfs->make<TH3F>("fTranHist3D_pion_2", "fTranHist3D_pion_2;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_pion_3 = tfs->make<TH3F>("fTranHist3D_pion_3", "fTranHist3D_pion_3;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_pion_4 = tfs->make<TH3F>("fTranHist3D_pion_4", "fTranHist3D_pion_4;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_pion_5 = tfs->make<TH3F>("fTranHist3D_pion_5", "fTranHist3D_pion_5;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);

  fLongProfile_pion = tfs->make<TProfile>("fLongProfile_pion", "fLongProfile_pion;t;Q", LBINS, LMIN, LMAX);
  fTranProfile_pion = tfs->make<TProfile>("fTranProfile_pion", "fTranProfile_pion;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_pion_1 = tfs->make<TProfile>("fTranProfile_pion_1", "fTranProfile_pion_1;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_pion_2 = tfs->make<TProfile>("fTranProfile_pion_2", "fTranProfile_pion_2;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_pion_3 = tfs->make<TProfile>("fTranProfile_pion_3", "fTranProfile_pion_3;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_pion_4 = tfs->make<TProfile>("fTranProfile_pion_4", "fTranProfile_pion_4;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_pion_5 = tfs->make<TProfile>("fTranProfile_pion_5", "fTranProfile_pion_5;dist;Q", TBINS, TMIN, TMAX);

  // others
  fLongHist3D_others = tfs->make<TH3F>("fLongHist3D_others","fLongHist3D_others; energy (GeV);t (radiation length 14 cm);Q", EBINS, EMIN, EMAX, LBINS, LMIN, LMAX, QBINS, QMIN, QMAX);
  fTranHist3D_others = tfs->make<TH3F>("fTranHist3D_others", "fTranHist3D_others;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_others_1 = tfs->make<TH3F>("fTranHist3D_others_1", "fTranHist3D_others_1;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_others_2 = tfs->make<TH3F>("fTranHist3D_others_2", "fTranHist3D_others_2;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_others_3 = tfs->make<TH3F>("fTranHist3D_others_3", "fTranHist3D_others_3;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_others_4 = tfs->make<TH3F>("fTranHist3D_others_4", "fTranHist3D_others_4;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);
  fTranHist3D_others_5 = tfs->make<TH3F>("fTranHist3D_others_5", "fTranHist3D_others_5;energy (GeV);dist (cm);Q", EBINS, EMIN, EMAX, TBINS, TMIN, TMAX, QBINS, QMIN, QMAX);

  fLongProfile_others = tfs->make<TProfile>("fLongProfile_others", "fLongProfile_others;t;Q", LBINS, LMIN, LMAX);
  fTranProfile_others = tfs->make<TProfile>("fTranProfile_others", "fTranProfile_others;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_others_1 = tfs->make<TProfile>("fTranProfile_others_1", "fTranProfile_others_1;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_others_2 = tfs->make<TProfile>("fTranProfile_others_2", "fTranProfile_others_2;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_others_3 = tfs->make<TProfile>("fTranProfile_others_3", "fTranProfile_others_3;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_others_4 = tfs->make<TProfile>("fTranProfile_others_4", "fTranProfile_others_4;dist;Q", TBINS, TMIN, TMAX);
  fTranProfile_others_5 = tfs->make<TProfile>("fTranProfile_others_5", "fTranProfile_others_5;dist;Q", TBINS, TMIN, TMAX);
}

void ShowerTemplateMaker::beginJob() {
  mf::LogInfo("ShowerTemplateMaker") << ".... begin of job ...." << endl;
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

bool ShowerTemplateMaker::insideFV( double vertex[3] ) {
  double x = vertex[0];
  double y = vertex[1];
  double z = vertex[2];

  if (x >= fFidVolXmin && x <= fFidVolXmax &&
      y >= fFidVolYmin && y <= fFidVolYmax &&
      x >= fFidVolZmin && z <= fFidVolZmax) {
    return true;
  }
  else {
    return false;
  }
}

void ShowerTemplateMaker::saveShowerInformation(art::Ptr<recob::Shower> myShower, art::FindManyP<recob::Hit> myfmshhits, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> myparticles_per_hit, std::vector<std::unordered_map<int, double>> myallhit_trkide, art::Event const& e) {
  // FV cut using recon
  if (fFVCutOnRecon) { // using FV cut on recon shower vertex
    // check if inside Fiducial Volume
    double check_vertex[3];
    check_vertex[0] = myShower->ShowerStart().X();
    check_vertex[1] = myShower->ShowerStart().Y();
    check_vertex[2] = myShower->ShowerStart().Z();

    if (!insideFV(check_vertex)) {
      cout << "\n **** Shower vertex is NOT inside the Fiducial Volume. RETURN ****" << endl;
      return;
    }
   } 

  // may need to consider the SCE correction for position
  shower_id.push_back(myShower->ID());
  sh_start_X.push_back(myShower->ShowerStart().X());
  sh_start_Y.push_back(myShower->ShowerStart().Y());
  sh_start_Z.push_back(myShower->ShowerStart().Z());
  sh_direction_X.push_back(myShower->Direction().X()); // Direction(): direction cosines at the start of the shower
  sh_direction_Y.push_back(myShower->Direction().Y());
  sh_direction_Z.push_back(myShower->Direction().Z());
  sh_length.push_back(myShower->Length()); 

  //cout << "myShower.key(): " << myShower.key() << endl;
  cout << "myShower->ID(): " << myShower->ID() << endl;

  cout << "myShower->ShowerStart().X(): " << myShower->ShowerStart().X() << endl;
  cout << "myShower->ShowerStart().Y(): " << myShower->ShowerStart().Y() << endl;
  cout << "myShower->ShowerStart().Z(): " << myShower->ShowerStart().Z() << endl;
  
  cout << "myShower->Direction().X(): " << myShower->Direction().X() << endl;
  cout << "myShower->Direction().Y(): " << myShower->Direction().Y() << endl;
  cout << "myShower->Direction().Z(): " << myShower->Direction().Z() << endl;

  cout << "myShower->Length(): " << myShower->Length() << endl;



  TVector3 dir;
  dir[0] = myShower->Direction().X();
  dir[1] = myShower->Direction().Y();
  dir[2] = myShower->Direction().Z();
  dir = dir.Unit();

  // shower hits
  std::vector<art::Ptr<recob::Hit>> sh_hits;
  sh_hits = myfmshhits.at(myShower.key());

  // note here we consider hits on all planes to associate a partile to the shower
  std::vector<simb::MCParticle const *> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const *> match_vec;
  std::unordered_map<int, double> shhit_trkide;

  //double sh_total_Q = 0.0;
  //double sh_total_Q_1 = 0.0;
  //double sh_total_Q_0 = 0.0;
  cout << "sh_hits.size(): " << sh_hits.size() << endl;
  //int testCountHitPlane2 = 0;
  //int testCountHitPlane1 = 0;
  //int testCountHitPlane0 = 0;
  std::vector<double> sh_total_Q (3, 0.0);
  std::vector<int> sh_hit_count (3, 0);
  for (size_t ihit = 0; ihit < sh_hits.size(); ++ihit) {
    art::Ptr<recob::Hit> hit = sh_hits[ihit];

    if (hit->WireID().Plane == 2) {
      sh_total_Q[2] += hit->Integral();
      sh_hit_count[2] += 1;
    }
    else if (hit->WireID().Plane == 1) {
      sh_total_Q[1] += hit->Integral();
      sh_hit_count[1] += 1;
    }
    else if (hit->WireID().Plane == 0) {
      sh_total_Q[0] += hit->Integral();
      sh_hit_count[0] += 1;
    }
    particle_vec.clear();
    match_vec.clear();
    myparticles_per_hit.get(hit.key(), particle_vec, match_vec);
    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
      if (!(pi_serv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId()))) {
        cout << "....this hit is from data...." << endl;
        cout << match_vec[i_p]->energy << endl;
      }
      // here, from the trackID of the MCParticle, we can find the EveTrackId (the main mother at very beginning. We should use this to classify a shower)
      shhit_trkide[pi_serv->TrackIdToEveTrackId(particle_vec.at(i_p)->TrackId())] += match_vec[i_p]->energy; // match_vec[i_p]->numElectrons;
    }
  } // end of for sh_hits.size() ihit
  
  /*
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  if (sh_total_Q > 0) {
     double Q_e = fCalorimetryAlg.ElectronsFromADCArea(sh_total_Q, 2);
     cout  << Q_e << endl;
     double rho = detprop->Density();            // LAr density in g/cm^3
     double Wion = 1000./util::kGeVToElectrons;  // 23.6 eV = 1e, Wion in MeV/e
     double E_field = detprop->Efield();        // Electric Field in the drift region in KV/cm
     double Beta = fModBoxB / (rho * E_field);
     double Alpha = fModBoxA;
     double E =  (exp(Beta * Wion * Q_e ) - Alpha) / Beta;
     cout << "E: " << E << endl;
     cout << "E_field: " << E_field << endl;
     cout << "Alpha = " << Alpha << endl;
     cout << "Beta: " << Beta << endl;
     cout << "util::kGeVToElectrons: " << util::kGeVToElectrons << endl; 
     cout << Q_e/0.62*Wion << endl;
     cout <<  sh_total_Q*243.7/0.6* 2.36e-05 << endl;
  }
  */

  int shhit_trackID = 0;
  double shhit_maxenergy = -99999.0;
  // ideally, each plane can measure the shower energy independently. Here, we sum over all planes to calculate the purity and completeness; may use each plane for this purpose.
  double shhit_partial_E = 0.0; // shower energy on all planes from the "shower particle"
  double shhit_total_E = 0.0; // shower energy on all planes
  for (auto const & p_id : shhit_trkide) {
    shhit_total_E += p_id.second;
    if (p_id.second > shhit_maxenergy) {
      shhit_partial_E = p_id.second;
      shhit_trackID = p_id.first;
      shhit_maxenergy = p_id.second;
    }
  }
  
  // primary electron shower
  if (shhit_trackID == MC_lepton_ID) {
    sh_hasPrimary_e.push_back(1);
  }
  else {
    sh_hasPrimary_e.push_back(0);
  }

  // we still keep the shower energy information for each plane
  std::vector<double> sh_energyReco(3);
  sh_energyReco.clear();
  for (int p = 0; p < 3; p++) {
    cout << "sh_total_Q[" << p << "]: " << sh_total_Q[p] << endl;
    cout << "sh_hit_count[" << p << "]: " << sh_hit_count[p] << endl;
    if (fShowerEnergyCorrection) {
      //sh_energyReco.push_back( (myShower->Energy().at(p)-fShowerEnergyCorrectionIntercept[p]) / fShowerEnergyCorrectionSlope[p] );
      sh_energyReco.push_back( (sh_total_Q[p]*fADCtoE[p]/fRecombination*23.6e-6 - fShowerEnergyCorrectionIntercept[p]) / fShowerEnergyCorrectionSlope[p] );
    }
    else {
      //sh_energyReco.push_back(myShower->Energy().at(p));
      sh_energyReco.push_back(sh_total_Q[p]*fADCtoE[p]/fRecombination*23.6e-6);
    }
  }
  showerEnergyReco.push_back(sh_energyReco);

  for (int p = 0; p< 3;p++) {
    cout << "sh_energyReco[" << p << "]: " <<  sh_energyReco[p] << endl;
    cout<<  "showerEnergyReco.back()[" << p << "]: " << showerEnergyReco.back()[p] << endl;
  }


  // to find the shower particle, we usually use hits on all planee; same thing for purity and completeness
  if (pi_serv->TrackIdToMotherParticle_P(shhit_trackID) && sh_total_Q[2] > 0) {
    cout << "shhit_trackID: " << shhit_trackID << endl;
    cout << "pi_serv->TrackIdToMotherParticle_P(shhit_trackID)->PdgCode(): " << pi_serv->TrackIdToMotherParticle_P(shhit_trackID)->PdgCode() << endl;
    const simb::MCParticle* sh_particle;
    sh_particle = pi_serv->TrackIdToParticle_P(shhit_trackID);
    cout << "sh_particle->PdgCode(): " << sh_particle->PdgCode() << endl;
    cout << "sh_particle->Mother(): " << sh_particle->Mother() << endl;

    showerParticlePDG.push_back(sh_particle->PdgCode());
    showerEnergyMC.push_back((sh_particle->E() - sh_particle->Mass())*1000.0); // kinetic energy [MeV]
    
    // here, we consider the energy on three planes, which is about three times of the particle energy deposits.
    double particle_total_E = myallhit_trkide[0][shhit_trackID] + myallhit_trkide[1][shhit_trackID] + myallhit_trkide[2][shhit_trackID];
    if (shhit_total_E > 0) {
      sh_purity.push_back(shhit_partial_E/shhit_total_E);
    }
    else {
      sh_purity.push_back(0.0);
    }
    if (particle_total_E > 0) {
      sh_completeness.push_back(shhit_partial_E/particle_total_E);
    }
    else {
      sh_completeness.push_back(0.0);
    }
    
    // ---------- shower profile ---------
    if (sh_purity.back() > fShowerCutPurity && sh_completeness.back() > fShowerCutCompleteness) { 
      cout << "myShower->Direction().X(): " << myShower->Direction().X() << endl;
      cout << "myShower->Direction().Y(): " << myShower->Direction().Y() << endl;
      cout << "myShower->Direction().Z(): " << myShower->Direction().Z() << endl;
      showerProfile(sh_hits, myShower->ShowerStart(), myShower->Direction(),myShower->Length(), sh_energyReco[2], sh_particle->PdgCode());
    }

  } // end if (pi_serv->TrackIdToMotherParticle_P(shhit_trackID))
  else {
    cout << "....www...." << endl;
    showerParticlePDG.push_back(-99999);
    showerEnergyMC.push_back(-99999.0);
    sh_purity.push_back(-99999.0);
    sh_completeness.push_back(-99999.0);
  }

}

void ShowerTemplateMaker::showerProfile(std::vector< art::Ptr<recob::Hit> >shhits, TVector3 shvertex, TVector3 shdirection, double shlength, double shenergy, int shpdg) {

  TH1F* ltemp = new TH1F("ltemp", "ltemp", LBINS, LMIN, LMAX);
  TH1F* ttemp = new TH1F("ttemp", "ttemp", TBINS, TMIN, TMAX);

  TH1F* ttemp_1 = new TH1F("ttemp_1", "ttemp_1", TBINS, TMIN, TMAX);
  TH1F* ttemp_2 = new TH1F("ttemp_2", "ttemp_2", TBINS, TMIN, TMAX);
  TH1F* ttemp_3 = new TH1F("ttemp_3", "ttemp_3", TBINS, TMIN, TMAX);
  TH1F* ttemp_4 = new TH1F("ttemp_4", "ttemp_4", TBINS, TMIN, TMAX);
  TH1F* ttemp_5 = new TH1F("ttemp_5", "ttemp_5", TBINS, TMIN, TMAX);


  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  double shVertexTime = 0; // in Ticks
  double shVertexWire = 0; // wire index
  cout << "shvertex[0]: " << shvertex[0] << endl;
  cout << "shvertex[1]: " << shvertex[1] << endl;
  cout << "shvertex[2]: " << shvertex[2] << endl;
  convertXYZtoGlobalWireTDC(shvertex, shVertexWire, shVertexTime);
  cout << "shVertexWire: " << shVertexWire << endl;
  cout << "shVertexTime: " << shVertexTime << endl;

  double shTwoTime = 0; // in Ticks
  double shTwoWire = 0; // wire index
  // to avoid a big uncertain in direction, we use direction*factorDir
  // todo: may need remove lenght or rethink how to use it
  // todo: may apply FV cut on the second point
  double factorDir = 5; // ideally, this should be larger, i.e, 20
  convertXYZtoGlobalWireTDC(shvertex+shdirection*factorDir, shTwoWire, shTwoTime);
  cout << "shTwoTime: " << shTwoTime << endl;
  cout << "shTwoWire: " << shTwoWire << endl;

  //auto const* geo = lar::providerFrom<geo::Geometry>();

  double t0 = detprop->TriggerOffset();

  std::vector<double> sh_hitLongDist;
  std::vector<double> sh_hitTranDist;
  std::vector<double> sh_hitChargeQ;

  cout << "shhits.size() : " << shhits.size() << endl;
  // re-define the vertex using very first hit
  // line equation (y2-y1)x - (x2-x1)y + y1x2-y2x1=0
  // let a = y2-y1; b = -(x2-x1); c= y1x2-y2x1
  double a = shTwoTime - shVertexTime;
  double b = -(shTwoWire - shVertexWire);
  double c = shVertexTime*shTwoWire - shTwoTime*shVertexWire;
  for (size_t i = 0; i < shhits.size(); ++i) {
    // collection plane; may need to consider the induction planes, too.
    if (shhits[i]->WireID().Plane != 2) continue;
    unsigned int hitGlobalWire = shhits[i]->WireID().Wire;
    double hitGlobalTDC = shhits[i]->PeakTime();
    
    // for a hit, we find the point on the shower projection line which is closest to the hit point
    // we only consider hits near the default shower vertex, i.e., using the following condition
    // todo: consider the shower direction
    if (shdirection.Z()>=0) {
      if (hitGlobalWire < shTwoWire) {
        double xx = (b*(b*hitGlobalWire - a*hitGlobalTDC) - a*c) / (a*a+b*b);
        double yy = (a*(-b*hitGlobalWire+a*hitGlobalTDC) - b*c) / (a*a+b*b);
        if (xx < shVertexWire) {
          shVertexWire = xx;
          shVertexTime = yy;
        }
      } 
    }
    else {
      if (hitGlobalWire > shTwoWire) {
        double xx = (b*(b*hitGlobalWire - a*hitGlobalTDC) - a*c) / (a*a+b*b);
        double yy = (a*(-b*hitGlobalWire+a*hitGlobalTDC) - b*c) / (a*a+b*b);
        if (xx < shVertexWire) {
          shVertexWire = xx;
          shVertexTime = yy;
        }
      } 
    }
  }


  for (size_t i = 0; i < shhits.size(); ++i) {
    // collection plane; may need to consider the induction planes, too.
    if (shhits[i]->WireID().Plane != 2) continue;
    unsigned int cryostat = shhits[i]->WireID().Cryostat;
    unsigned int tpc = shhits[i]->WireID().TPC;
    unsigned int hitGlobalWire = shhits[i]->WireID().Wire;
    unsigned int hitGlobalPlane = shhits[i]->WireID().Plane;
    double hitGlobalTDC = shhits[i]->PeakTime();
   
    if (tpc != 0 || cryostat !=0) {
      cout << "-----cryostat: " << cryostat << endl;
      cout << "-----tpc: " << tpc << endl;
    }

    // todo: rethink about: sometime the reconstructed vertex is not right, skip upstream hits;
    //if (hitGlobalWire <= shVertexWire) continue;
   
    double wirePitch = geom->WirePitch(shhits[i]->WireID());
    double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()); // [cm/us]
    // detprop->SamplingRate(): the period of the TPC readout electronics clock in [us]
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns

    if (i < 2 ){
      cout << "t0: " << t0 << endl;
      cout << "hitGlobalWire: " << hitGlobalWire << endl;
      cout << "hitGlobalPlane: " << hitGlobalPlane << endl;
      cout << "hitGlobalTDC: " << hitGlobalTDC << endl;
      cout << "wirePitch: " << wirePitch << " [cm]" << endl; 
      cout << "detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()): " << detprop->DriftVelocity(detprop->Efield(),detprop->Temperature()) << " [cm/us]" << endl;
      cout << "detprop->SamplingRate(): " << detprop->SamplingRate() << endl; 
      cout << "tickToDist: " << tickToDist << " [cm]" << endl; 
    }

    // todo: rethink about distance 
    double x_vertex = shVertexTime * tickToDist;
    double y_vertex = shVertexWire * wirePitch;

    double x2_vertex = shTwoTime * tickToDist;
    double y2_vertex = shTwoWire * wirePitch;

    double x_orth = (y2_vertex -y_vertex) + x_vertex;
    double y_orth = -(x2_vertex - x_vertex) + y_vertex;

    double xhit = hitGlobalTDC * tickToDist;
    double yhit = hitGlobalWire * wirePitch;

    double ldist = std::abs((y_orth - y_vertex)*xhit - (x_orth-x_vertex)*yhit + x_orth*y_vertex - y_orth*x_vertex) / std::sqrt( pow((y_orth-y_vertex), 2) + pow((x_orth -x_vertex), 2));
    double tdist = ((y2_vertex - y_vertex)*xhit - (x2_vertex-x_vertex)*yhit + x2_vertex*y_vertex - y2_vertex*x_vertex) / std::sqrt( pow((y2_vertex-y_vertex), 2) + pow((x2_vertex-x_vertex), 2));

    double to3D = 1. / sqrt( pow((x_vertex-x2_vertex)/factorDir,2) + pow((y_vertex-y2_vertex)/factorDir,2) ) ; // distance between two points in 3D space is one; reset by factorDir 
    ldist *= to3D;
    tdist *= to3D;
    double t = ldist/X0;
    // todo: add lifetime correction
    //double Q = shhits[i]->Integral() * fCalorimetryAlg.LifetimeCorrection(shhits[i]->PeakTime()+ t0);
    double Q = shhits[i]->Integral() * fADCtoE[2]; // [e-]
    
    sh_hitLongDist.push_back(t);
    sh_hitTranDist.push_back(tdist);
    sh_hitChargeQ.push_back(Q);

    ltemp->Fill(t, Q);
    ttemp->Fill(tdist, Q);

    if (t < 1) ttemp_1->Fill(tdist, Q);
    else if (t < 2) ttemp_2->Fill(tdist, Q);
    else if (t < 3) ttemp_3->Fill(tdist, Q);
    else if (t < 4) ttemp_4->Fill(tdist, Q);
    else if (t < 5) ttemp_5->Fill(tdist, Q);
  } //  end of loop shhits.size() i
 
  showerHitLongDist.push_back(sh_hitLongDist);
  showerHitTranDist.push_back(sh_hitTranDist);
  showerHitChargeQ.push_back(sh_hitChargeQ);

  cout << "shpdg: " << shpdg << endl;
  
  // primary electron
  if ( abs(shpdg) == 11 && sh_hasPrimary_e.back() == 1 ) {
    cout << "event: " << event << endl;
    cout << "run: " <<run << endl;
    cout << "subrun: " << subrun << endl;
    fillElectronHist(shenergy, ltemp, ttemp, ttemp_1, ttemp_2, ttemp_3, ttemp_4, ttemp_5);
  }
  // photon
  else if (abs(shpdg) == 22) fillPhotonHist(shenergy, ltemp, ttemp, ttemp_1, ttemp_2, ttemp_3, ttemp_4, ttemp_5);
  // proton
  else if (abs(shpdg) == 2212) fillProtonHist(shenergy, ltemp, ttemp, ttemp_1, ttemp_2, ttemp_3, ttemp_4, ttemp_5);
  // pion
  else if (abs(shpdg) == 211) fillPionHist(shenergy, ltemp, ttemp, ttemp_1, ttemp_2, ttemp_3, ttemp_4, ttemp_5);
  // others
  else fillOthersHist(shenergy, ltemp, ttemp, ttemp_1, ttemp_2, ttemp_3, ttemp_4, ttemp_5);
}

void ShowerTemplateMaker::convertXYZtoGlobalWireTDC(TVector3 pointXYZ,  double& pointGlobalWire, double& pointGlobalTDC) {
  // For a coordinate (x,y,z), we need to find the correct cryostat, tpc, wire index. Then we convert that to global wire and tdc.
  // x:0-256 cm, anode at x=0 cm; y = -116 cm - 116 cm; z = 0-1036 cm;
  // todo: rethink if more restriction is needed.
  if (pointXYZ[2] < 0) return;
  //double tpcWidth = 256.0; // [cm]

  int cryostat = 0;
  int tpc = 0;
  int plane = 2; // here, we only consider the collection plane
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;
  auto collectionPlane = geo::PlaneID(cryostat, tpc, plane);

  double pointTDC = detprop->ConvertXToTicks(pointXYZ[0], collectionPlane);
  double pointWire = geom->WireCoordinate(pointXYZ[1],pointXYZ[2], collectionPlane);

  /*
  geo::WireIDIntersection intersection;
  geom->WireIDsIntersect( geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 1)), geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 2)), intersection);
  geom->WireIDsIntersect(geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 0)), geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 2)), intersection);
  cout << intersection.y << endl;
  cout << intersection.z << endl;
  */
  
  cout << "XYZ....www..." << geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 0)) << endl;
  cout << "XYZ....www..." << geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 1)) << endl;
  cout << "XYZ....www..." << geom->WireCoordinate(pointXYZ[1],pointXYZ[2] ,geo::PlaneID(cryostat, tpc, 2)) << endl;
  pointGlobalTDC = pointTDC;
  pointGlobalWire = pointWire;
}

void ShowerTemplateMaker::fillElectronHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5) {
  for (int i = 0; i < LBINS; ++i) {
    //if (ltemp->GetBinContent(i+1) == 0) continue;
    fLongProfile_electron->Fill(ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
    fLongHist3D_electron->Fill(shenergy, ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp->GetBinContent(i+1) == 0) continue;
    fTranProfile_electron->Fill(ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
    fTranHist3D_electron->Fill(shenergy, ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_1->GetBinContent(i+1) == 0) continue;
    fTranProfile_electron_1->Fill(ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
    fTranHist3D_electron_1->Fill(shenergy, ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_2->GetBinContent(i+1) == 0) continue;
    fTranProfile_electron_2->Fill(ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
    fTranHist3D_electron_2->Fill(shenergy, ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_3->GetBinContent(i+1) == 0) continue;
    fTranProfile_electron_3->Fill(ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
    fTranHist3D_electron_3->Fill(shenergy, ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_4->GetBinContent(i+1) == 0) continue;
    fTranProfile_electron_4->Fill(ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
    fTranHist3D_electron_4->Fill(shenergy, ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_5->GetBinContent(i+1) == 0) continue;
    fTranProfile_electron_5->Fill(ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
    fTranHist3D_electron_5->Fill(shenergy, ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
  }

}

void ShowerTemplateMaker::fillPhotonHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5) {
  for (int i = 0; i < LBINS; ++i) {
    //if (ltemp->GetBinContent(i+1) == 0) continue;
    fLongProfile_photon->Fill(ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
    fLongHist3D_photon->Fill(shenergy, ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp->GetBinContent(i+1) == 0) continue;
    fTranProfile_photon->Fill(ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
    fTranHist3D_photon->Fill(shenergy, ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_1->GetBinContent(i+1) == 0) continue;
    fTranProfile_photon_1->Fill(ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
    fTranHist3D_photon_1->Fill(shenergy, ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_2->GetBinContent(i+1) == 0) continue;
    fTranProfile_photon_2->Fill(ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
    fTranHist3D_photon_2->Fill(shenergy, ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_3->GetBinContent(i+1) == 0) continue;
    fTranProfile_photon_3->Fill(ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
    fTranHist3D_photon_3->Fill(shenergy, ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_4->GetBinContent(i+1) == 0) continue;
    fTranProfile_photon_4->Fill(ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
    fTranHist3D_photon_4->Fill(shenergy, ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_5->GetBinContent(i+1) == 0) continue;
    fTranProfile_photon_5->Fill(ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
    fTranHist3D_photon_5->Fill(shenergy, ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
  }
}

void ShowerTemplateMaker::fillProtonHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5) {
  for (int i = 0; i < LBINS; ++i) {
    //if (ltemp->GetBinContent(i+1) == 0) continue;
    fLongProfile_proton->Fill(ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
    fLongHist3D_proton->Fill(shenergy, ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp->GetBinContent(i+1) == 0) continue;
    fTranProfile_proton->Fill(ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
    fTranHist3D_proton->Fill(shenergy, ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_1->GetBinContent(i+1) == 0) continue;
    fTranProfile_proton_1->Fill(ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
    fTranHist3D_proton_1->Fill(shenergy, ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_2->GetBinContent(i+1) == 0) continue;
    fTranProfile_proton_2->Fill(ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
    fTranHist3D_proton_2->Fill(shenergy, ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_3->GetBinContent(i+1) == 0) continue;
    fTranProfile_proton_3->Fill(ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
    fTranHist3D_proton_3->Fill(shenergy, ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_4->GetBinContent(i+1) == 0) continue;
    fTranProfile_proton_4->Fill(ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
    fTranHist3D_proton_4->Fill(shenergy, ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_5->GetBinContent(i+1) == 0) continue;
    fTranProfile_proton_5->Fill(ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
    fTranHist3D_proton_5->Fill(shenergy, ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
  }
}

void ShowerTemplateMaker::fillPionHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5) {
  for (int i = 0; i < LBINS; ++i) {
    //if (ltemp->GetBinContent(i+1) == 0) continue;
    fLongProfile_pion->Fill(ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
    fLongHist3D_pion->Fill(shenergy, ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp->GetBinContent(i+1) == 0) continue;
    fTranProfile_pion->Fill(ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
    fTranHist3D_pion->Fill(shenergy, ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_1->GetBinContent(i+1) == 0) continue;
    fTranProfile_pion_1->Fill(ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
    fTranHist3D_pion_1->Fill(shenergy, ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_2->GetBinContent(i+1) == 0) continue;
    fTranProfile_pion_2->Fill(ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
    fTranHist3D_pion_2->Fill(shenergy, ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_3->GetBinContent(i+1) == 0) continue;
    fTranProfile_pion_3->Fill(ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
    fTranHist3D_pion_3->Fill(shenergy, ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_4->GetBinContent(i+1) == 0) continue;
    fTranProfile_pion_4->Fill(ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
    fTranHist3D_pion_4->Fill(shenergy, ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_5->GetBinContent(i+1) == 0) continue;
    fTranProfile_pion_5->Fill(ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
    fTranHist3D_pion_5->Fill(shenergy, ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
  }
}

void ShowerTemplateMaker::fillOthersHist(double shenergy, TH1F* ltemp, TH1F* ttemp, TH1F* ttemp_1, TH1F* ttemp_2, TH1F* ttemp_3, TH1F* ttemp_4, TH1F* ttemp_5) {
  for (int i = 0; i < LBINS; ++i) {
    //if (ltemp->GetBinContent(i+1) == 0) continue;
    fLongProfile_others->Fill(ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
    fLongHist3D_others->Fill(shenergy, ltemp->GetBinCenter(i+1), ltemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp->GetBinContent(i+1) == 0) continue;
    fTranProfile_others->Fill(ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
    fTranHist3D_others->Fill(shenergy, ttemp->GetBinCenter(i+1), ttemp->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_1->GetBinContent(i+1) == 0) continue;
    fTranProfile_others_1->Fill(ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
    fTranHist3D_others_1->Fill(shenergy, ttemp_1->GetBinCenter(i+1), ttemp_1->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_2->GetBinContent(i+1) == 0) continue;
    fTranProfile_others_2->Fill(ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
    fTranHist3D_others_2->Fill(shenergy, ttemp_2->GetBinCenter(i+1), ttemp_2->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_3->GetBinContent(i+1) == 0) continue;
    fTranProfile_others_3->Fill(ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
    fTranHist3D_others_3->Fill(shenergy, ttemp_3->GetBinCenter(i+1), ttemp_3->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_4->GetBinContent(i+1) == 0) continue;
    fTranProfile_others_4->Fill(ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
    fTranHist3D_others_4->Fill(shenergy, ttemp_4->GetBinCenter(i+1), ttemp_4->GetBinContent(i+1));
  }

  for (int i = 0; i < TBINS; ++i) {
    //if (ttemp_5->GetBinContent(i+1) == 0) continue;
    fTranProfile_others_5->Fill(ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
    fTranHist3D_others_5->Fill(shenergy, ttemp_5->GetBinCenter(i+1), ttemp_5->GetBinContent(i+1));
  }
}

void ShowerTemplateMaker::analyze(art::Event const& e)
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

      if (!fFVCutOnRecon) { // using FV cut on MC interaction vertex
        // check if inside Fiducial Volume
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
  } // end of loop mctruthlist i

  //------------ MCTruth: Geant4 -------------------
  // Note: generator level MCPartilceList is different from Geant4 level MCParticleList.  MCParticleList(Geant4) contains all particles in MCParticleList(generator) but their the indexes (TrackIds) are different.
  simb::MCParticle *MC_lepton = NULL; // Geant4 level
  const sim::ParticleList& plist = pi_serv->ParticleList();
  //ideally one event is corresponding to one mctruth
  cout << "plist.size():" << plist.size() << endl;
  if (plist.size() != 1) {
    cout << "Warning: this event contains more than one mctruth" << endl;
  }
  int testCountPrimarylepton = 0;
  for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar) {
    simb::MCParticle *particle = 0;
    particle = ipar->second; // first is index(TrackId), second is value (point address)
    auto & truth = pi_serv->ParticleToMCTruth_P(particle);
    // beam neutrino only, only look at electron
    //if ( truth->Origin()==simb::kBeamNeutrino && std::abs(particle->PdgCode()) == 11 &&  particle->Mother()==0){ // primary lepton; Mother() = 0 means e^{-} for v+n=e^{-}+p
    if ( truth->Origin()==simb::kBeamNeutrino && particle->Mother()==0 ) {
      
      if (std::abs(particle->PdgCode()) == 11) {// primary lepton; Mother() = 0 means e^{-} for v+n=e^{-}+p
        MC_lepton_ID = particle->TrackId();
        MC_lepton = particle;
        testCountPrimarylepton += 1;
      }
      else if (std::abs(particle->PdgCode()) >= 11 && std::abs(particle->PdgCode()) <= 16) {
        MC_lepton_ID = particle->TrackId();
        MC_lepton = particle;
        testCountPrimarylepton += 1;
      }
    }
  } // end of loop plist ipar

  cout << "testCountPrimarylepton: " << testCountPrimarylepton << endl;
  cout << "MC_lepton_ID: " << MC_lepton_ID << endl;
  cout << "MC_lepton->TrackId(): " << MC_lepton->TrackId() << endl;
  cout << "MC_lepton->PdgCode(): " << MC_lepton->PdgCode() << endl;
  cout << "MC_lepton->E(): " << MC_lepton->E() << endl;
  cout << "leptonEnergy_truth: " << leptonEnergy_truth  << " [GeV]" << endl;

  // -------- hit information --------------
  //find all hits from the primary electron; build a map for total hit charges corresponding to MC_lepton_ID (Geant4)
  std::vector<simb::MCParticle const *> allhit_particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const *> allhit_match_vec;
  std::vector<std::unordered_map<int, double>> allhit_trkide(3);
  
  cout << "hitlist.size(): " << hitlist.size() << endl;
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
    cout << "allhit_trkide[p][MC_lepton_ID]: " << allhit_trkide[p][MC_lepton_ID] << " [MeV]" << endl;
  }

  // --------- shower information -----------
  // ------- shower associated with slice ----------
  cout << "slicelist.size(): " << slicelist.size() << endl;
  int testCountPFPs = 0;
  int testCountShowers = 0;
  for (size_t i = 0; i < slicelist.size(); ++i) {
    cout << "slicelist[i].key(): " << slicelist[i].key() << endl;

    if (fm_pfp_slice_assns.isValid()) {
      std::vector<art::Ptr<recob::PFParticle>> pfs = fm_pfp_slice_assns.at(i);
      cout << "pfs.size():  " << pfs.size() << endl;
      testCountPFPs += pfs.size();
      double maxenergy_sh_energy = -99999.0;
      art::Ptr<recob::Shower> largestShower;
      
      // loop PFParticle for the associated showers
      for (size_t j = 0; j < pfs.size(); ++j) {
        if (fm_pfp_shower_assns.isValid()) {
          std::vector<art::Ptr<recob::Shower>> showers = fm_pfp_shower_assns.at(pfs[j].key()); // here, should use pfs[j].key() to get the right showers.
          if(showers.size() > 0) cout << "showers.size(): " << showers.size() << endl;
          testCountShowers += showers.size();
          for (size_t k = 0; k < showers.size(); ++k) {
            //if (event == 5598 && subrun == 90 && event == 4547) {
              cout << "\n\n" << endl;
              cout << "event: " << event << endl;
              cout << "subrun: " << subrun << endl;
              cout << "run: " << run << endl;
              cout << "showers[k]->ID(): " << showers[k]->ID() << endl;
              cout << "showers[k]->Length(): " << showers[k]->Length() << endl;
              //cout << "showers[k]->Energy().at(0): " << showers[k]->Energy().at(0) << endl;
              //cout << "showers[k]->Energy().at(1): " << showers[k]->Energy().at(1) << endl;
              //cout << "showers[k]->Energy().at(2): " << showers[k]->Energy().at(2) << endl;
            //}
            // use energy at plane 2 to select the largest shower

            // shower hits
            std::vector<art::Ptr<recob::Hit>> sh_hits = fmshhits.at(showers[k].key());
            cout << "sh_hits.size(): " << sh_hits.size() << endl;
            std::vector<double> shHitsVec (3, 0.);
            for (size_t ihit = 0; ihit < sh_hits.size(); ++ihit) {
              art::Ptr<recob::Hit> hit = sh_hits[ihit];
              if (hit->WireID().Plane == 2) shHitsVec[2] += hit->Integral();
              else if (hit->WireID().Plane == 1) shHitsVec[1] += hit->Integral();
              else if (hit->WireID().Plane == 0) shHitsVec[0] += hit->Integral();
            }
            
            // we require hits on collection plane and at least on one of the induction plane
            // todo: rethink about if using "charge == 0" is ok, may switch to number of hits 
            if ( shHitsVec[2] == 0 || (shHitsVec[0]+shHitsVec[1] == 0) ) continue;

            //double k_sh_energy = showers[k]->Energy().at(2);
            double k_sh_energy = shHitsVec[2]*fADCtoE[2]/fRecombination*23.6e-6; //[MeV]
            if (k_sh_energy > maxenergy_sh_energy) {
              maxenergy_sh_energy = k_sh_energy;
              largestShower = showers[k];
            }

            // save all showers in a slice
            if (!fSaveOnlyLargestShowerPerSlice) {
              slice_id.push_back(slicelist[i].key());
              saveShowerInformation(showers[k], fmshhits, particles_per_hit, allhit_trkide, e);
            }
          } // end of for showers.size() k

        } // end of if (fm_pfp_shower_assns.isValid()) 
      } // end of for pfs.size() j

      if (fSaveOnlyLargestShowerPerSlice && largestShower) {
        slice_id.push_back(slicelist[i].key());
        saveShowerInformation(largestShower, fmshhits, particles_per_hit, allhit_trkide, e);
      } 

    } // end of if (fm_pfp_slice_assns.isValid())
  } // end of or slicelist.size() i

  cout << "testCountPFPs: " << testCountPFPs << endl;
  cout << "testCountShowers: " << testCountShowers << endl;

  fEventTree->Fill();
}

void ShowerTemplateMaker::reset() {
  run = -99999;
  subrun = -99999;
  event = -99999;

  // clear the vectors
  slice_id.clear();
  shower_id.clear();
  sh_start_X.clear();
  sh_start_Y.clear();
  sh_start_Z.clear();
  sh_direction_X.clear();
  sh_direction_Y.clear();
  sh_direction_Z.clear();
  sh_length.clear();
  showerParticlePDG.clear();
  sh_hasPrimary_e.clear();
  showerEnergyMC.clear();
  sh_purity.clear();
  sh_completeness.clear();
  showerEnergyReco.clear();
  showerHitLongDist.clear();
  showerHitTranDist.clear();
  showerHitChargeQ.clear();
}

DEFINE_ART_MODULE(ShowerTemplateMaker)

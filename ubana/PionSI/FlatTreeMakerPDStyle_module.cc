////////////////////////////////////////////////////////////////////////
// Class:       FlatTreeMakerPDStyle
// Plugin Type: analyzer (art v3_02_06)
// File:        FlatTreeMakerPDStyle_module.cc
//
// Generated at Wed Jul 17 09:17:30 2019 by Kirsty Duffy using cetskelgen
// from cetlib version v3_07_02.
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
// #include "art_root_io/TFileService.h"

// root includes
#include "TTree.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

class FlatTreeMakerPDStyle;


class FlatTreeMakerPDStyle : public art::EDAnalyzer {
public:
  explicit FlatTreeMakerPDStyle(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlatTreeMakerPDStyle(FlatTreeMakerPDStyle const&) = delete;
  FlatTreeMakerPDStyle(FlatTreeMakerPDStyle&&) = delete;
  FlatTreeMakerPDStyle& operator=(FlatTreeMakerPDStyle const&) = delete;
  FlatTreeMakerPDStyle& operator=(FlatTreeMakerPDStyle&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;
  void beginJob() override;
  void endJob() override;
  void reset();

private:

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  /******************************/
  // Note: this is copied from a ProtoDUNE analyzer and trying to keep the same variable names. In protoDUNE this would be ruth level info of the primary beam particle that generated the event
  // Our equivalent is the neutrino-induced particle that exits the nucleus
  // But for consistency in downstream code, I will continue to use these variable names
  int true_beam_PDG;
  int true_beam_ID;
  std::string true_beam_EndProcess;
  std::vector< double > true_beam_X;
  std::vector< double > true_beam_Y;
  std::vector< double > true_beam_Z;
  std::vector< double > true_beam_PX;
  std::vector< double > true_beam_PY;
  std::vector< double > true_beam_PZ;
  //Truth level info of the daughter MCParticles coming out of the
  //true primary particle
  std::vector< int > true_beam_daughter_PDGs;
  std::vector< int > true_beam_daughter_IDs;
  //FCL pars
  std::string fGeneratorTag;
  bool fVerbose;

};


FlatTreeMakerPDStyle::FlatTreeMakerPDStyle(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}   ,
    fGeneratorTag(p.get<std::string>("GeneratorTag")),
    fVerbose(p.get<bool>("Verbose")){}


void FlatTreeMakerPDStyle::analyze(art::Event const& evt)
{
  //reset containers
  reset();
  run    = evt.run();
  subrun = evt.subRun();
  event  = evt.id().event();

  // // Get various utilities
  // art::ServiceHandle<cheat::BackTrackerService>         bt_serv;
  // art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  // protoana::ProtoDUNETruthUtils                         truthUtil;
  // ////////////////////////////////////////

  // Note: this is copied from a ProtoDUNE analyzer and trying to keep the same variable names. In protoDUNE this function would get the true beam particle that generated the event
  // Our equivalent is the neutrino-induced particle that exits the nucleus
  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  evt.getByLabel(fGeneratorTag, mcp_h);
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  art::fill_ptr_vector(mcp_v, mcp_h);

  // Get true neutrino and find out its MCParticle ID
  auto const mctruth_h = evt.getValidHandle<std::vector<simb::MCTruth>>("generator"); // Get only GENIE MCtruth
  simb::MCTruth top_mctruth = mctruth_h->at(0);
  simb::MCNeutrino top_nu = top_mctruth.GetNeutrino();
  simb::MCParticle top_neutrino = top_nu.Nu();
  int nu_MCPID = top_neutrino.TrackId();

  // Now loop over all MCParticles in the event
  for (auto mcp : mcp_v) {

    // We only care about MCParticles that are true daughters of the neutrino - skip if parent is not neutrino
    if (mcp->Mother() != nu_MCPID) continue;

    // We only want to record weights for true pions and protons, so skip other particles
    true_beam_EndProcess  = mcp->EndProcess();
    true_beam_PDG         = mcp->PdgCode();
    true_beam_ID          = mcp->TrackId();
    if ( !( true_beam_PDG == 211 || true_beam_PDG == 2212 ) ){
      continue;
    }

    if (fVerbose) std::cout << "Looking at neutrino daughter particle with PDG " << true_beam_PDG << std::endl;

    //Here: add in the stepping info from MCTrajectory
    //Get a handle to the Geometry service to look up AuxDetGeos from module numbers
    art::ServiceHandle < geo::Geometry > fGeometryService;
    for( size_t i = 0; i < mcp->NumberTrajectoryPoints(); ++i ){
      double X = mcp->Position(i).X();
      double Y = mcp->Position(i).Y();
      double Z = mcp->Position(i).Z();
      geo::Point_t testpoint1 { X, Y, Z };
      const TGeoMaterial* testmaterial1 = fGeometryService->Material( testpoint1 );
      if( fVerbose ){
        std::cout << "At point (" << X << ", " << Y << ", " << Z << ") " << testmaterial1->GetName()
                  << " " << testmaterial1->GetZ() << " " << testmaterial1->GetA()
                  << " " << testmaterial1->GetDensity() << "\n";
      }
      //For now, just going to reweight the points within the LAr of the TPC
      if ( /*Z > -1. &&*/ !strcmp( testmaterial1->GetName(), "LAr" ) ){
        if(fVerbose) std::cout << "Will save this" << std::endl;
        true_beam_X.push_back( X );
        true_beam_Y.push_back( Y );
        true_beam_Z.push_back( Z );

        true_beam_PX.push_back( mcp->Px(i) );
        true_beam_PY.push_back( mcp->Py(i) );
        true_beam_PZ.push_back( mcp->Pz(i) );
      }

    } // end loop over trajectory points

    // Now find daughters of the MCP
    for( int i = 0; i < mcp->NumberDaughters(); i++ ){
      int daughterID = mcp->Daughter(i);
      if (fVerbose) std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
      for (auto test_mcp : mcp_v){
        if (test_mcp->TrackId() == daughterID){
          int pid = test_mcp->PdgCode();
          true_beam_daughter_PDGs.push_back(pid);
          true_beam_daughter_IDs.push_back( test_mcp->TrackId() );
          break;
        }
      }

    }
    fTree->Fill();
    reset();
  } // loop over MCPs

} // end analyze

void FlatTreeMakerPDStyle::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("FlatTree","");
  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_EndProcess", &true_beam_EndProcess);
  fTree->Branch("true_beam_X", &true_beam_X);
  fTree->Branch("true_beam_Y", &true_beam_Y);
  fTree->Branch("true_beam_Z", &true_beam_Z);
  fTree->Branch("true_beam_PX", &true_beam_PX);
  fTree->Branch("true_beam_PY", &true_beam_PY);
  fTree->Branch("true_beam_PZ", &true_beam_PZ);
  fTree->Branch("true_beam_daughter_PDGs", &true_beam_daughter_PDGs);
  fTree->Branch("true_beam_daughter_IDs", &true_beam_daughter_IDs);

}

void FlatTreeMakerPDStyle::endJob()
{

}

void FlatTreeMakerPDStyle::reset()
{
  true_beam_EndProcess = "";
  true_beam_PDG = -1;
  true_beam_ID = -1;
  true_beam_daughter_PDGs.clear();
  true_beam_daughter_IDs.clear();
  true_beam_X.clear();
  true_beam_Y.clear();
  true_beam_Z.clear();
  true_beam_PX.clear();
  true_beam_PY.clear();
  true_beam_PZ.clear();
}

DEFINE_ART_MODULE(FlatTreeMakerPDStyle)

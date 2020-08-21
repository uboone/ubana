#ifndef MinimalExample_module
#define MinimalExample_module

#include "AnaHelper.h"

// Analyzer class
class MinimalExample : public art::EDAnalyzer
{
public:
  explicit MinimalExample(fhicl::ParameterSet const & pset);
  virtual ~MinimalExample();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Declare fhiclcpp variables
  std::string fMessage;

  // Declare services
  geo::GeometryCore const* fGeometry; // Pointer to the Geometry service

  // Declare trees
  TTree *tDataTree;
  TTree *tMetaTree;

  // Declare analysis variables
  int run, subrun, event;

  // Declare analysis functions
  void ClearData();
}; // End class MinimalExample

MinimalExample::MinimalExample(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMessage(pset.get<std::string>("message"))  
{} // END constructor MinimalExample

MinimalExample::~MinimalExample()
{} // END destructor MinimalExample

void MinimalExample::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tMetaTree = tfs->make<TTree>("Metadata","");
  tMetaTree->Branch("message",&fMessage);
  tMetaTree->Fill();

  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run,"run/I");
  tDataTree->Branch("subrun",&subrun,"subrun/I");
  tDataTree->Branch("event",&event,"event/I");

  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
} // END function beginJob

void MinimalExample::endJob()
{
} // END function endJob

void MinimalExample::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
} // END function ClearData

void MinimalExample::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  printf("\n-------------------------------------------------------\n");
  
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);
  
  // Start performing analysis
  printf("Here's your message: %s\n", fMessage.c_str());

  // Fill tree and finish event loop
  tDataTree->Fill();
  printf("-------------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(MinimalExample)

#endif // END def MinimalExample_module
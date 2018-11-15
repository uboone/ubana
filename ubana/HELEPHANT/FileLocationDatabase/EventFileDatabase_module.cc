#ifndef EventFileDatabase_module
#define EventFileDatabase_module

#include "AnaHelper.h"

// Analyzer class
class EventFileDatabase : public art::EDAnalyzer
{
public:
  explicit EventFileDatabase(fhicl::ParameterSet const & pset);
  virtual ~EventFileDatabase();
  void respondToOpenInputFile(art::FileBlock const& fb);
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Declare trees
  TTree *tDataTree;

  // Declare analysis variables
  int run, subrun, event;
  int f_run, f_subrun, f_event;
  std::string fileName;
  std::string f_fileName;

  // Declare analysis functions
  void ClearData();
}; // End class EventFileDatabase

EventFileDatabase::EventFileDatabase(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor EventFileDatabase

EventFileDatabase::~EventFileDatabase()
{} // END destructor EventFileDatabase

void EventFileDatabase::respondToOpenInputFile(art::FileBlock const& fb)
{
 fileName = fb.fileName();
}

void EventFileDatabase::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tDataTree = tfs->make<TTree>("EventFileDatabase","");
  tDataTree->Branch("run",&f_run);
  tDataTree->Branch("subrun",&f_subrun);
  tDataTree->Branch("event",&f_event);
  tDataTree->Branch("fileName",&f_fileName);
} // END function beginJob

void EventFileDatabase::endJob()
{
} // END function endJob

void EventFileDatabase::ClearData()
{
  f_run = -1;
  f_subrun = -1;
  f_event = -1;
  f_fileName = "";
} // END function ClearData

void EventFileDatabase::analyze(art::Event const & evt)
{
  // Start by clearing all data.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  // Assign values
  f_run = run;
  f_subrun = subrun;
  f_event = event;
  f_fileName = fileName;
  tDataTree->Fill();

} // END function analyze

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(EventFileDatabase)

#endif // END def EventFileDatabase_module
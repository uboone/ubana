#ifndef FindFileWithEvent_module
#define FindFileWithEvent_module

#include "AnaHelper.h"

// Analyzer class
class FindFileWithEvent : public art::EDAnalyzer
{
public:
  explicit FindFileWithEvent(fhicl::ParameterSet const & pset);
  virtual ~FindFileWithEvent();
  void respondToOpenInputFile(art::FileBlock const& fb);
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Declare fhiclcpp variables
  std::vector<int> q_run;
  std::vector<int> q_subrun;
  std::vector<int> q_event;
  bool q_all;

  // Declare trees
  TTree *tDataTree;

  // Declare analysis variables
  int run, subrun, event;
  int f_run, f_subrun, f_event;
  std::string fileName;
  std::string f_fileName;
  int lenQueriedEvents;
  bool isCorrectRun, isCorrectSubrun, isCorrectEvent, isToStore;

  // Declare analysis functions
  void ClearData();
}; // End class FindFileWithEvent

FindFileWithEvent::FindFileWithEvent(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    q_run(pset.get<std::vector<int>>("queriedRun")),
    q_subrun(pset.get<std::vector<int>>("queriedSubrun")),
    q_event(pset.get<std::vector<int>>("queriedEvent")),
    q_all(pset.get<bool>("queryAll"))
{} // END constructor FindFileWithEvent

FindFileWithEvent::~FindFileWithEvent()
{} // END destructor FindFileWithEvent

void FindFileWithEvent::respondToOpenInputFile(art::FileBlock const& fb)
{
 fileName = fb.fileName();
}

void FindFileWithEvent::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&f_run);
  tDataTree->Branch("subrun",&f_subrun);
  tDataTree->Branch("event",&f_event);
  tDataTree->Branch("fileName",&f_fileName);

  lenQueriedEvents = q_event.size();
} // END function beginJob

void FindFileWithEvent::endJob()
{
} // END function endJob

void FindFileWithEvent::ClearData()
{
  f_run = -1;
  f_subrun = -1;
  f_event = -1;
  f_fileName = "";
  isToStore = false;
} // END function ClearData

void FindFileWithEvent::analyze(art::Event const & evt)
{
  // Core analysis. This will be repeated event by event.

  // Start by clearing all data.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  // Check if event is contained in the list of events to be found. If queryAll == true, then information for all events will be stored.
  if (q_all) isToStore = true;
  else
  {
    for (int i=0;i<lenQueriedEvents;i++)
    {
      isCorrectRun = (q_run[i] == run);
      isCorrectSubrun = (q_subrun[i] == subrun);
      isCorrectEvent = (q_event[i] == event);
      if (isCorrectRun && isCorrectSubrun && isCorrectEvent) isToStore = true;
    }
  }

  // If event information is to be stored, fill tree.
  if (isToStore)
  {
    printf("||FOUND EVENT %i [RUN %i, SUBRUN %i] in file %s||\n", event, run, subrun, fileName.c_str());
    f_run = run;
    f_subrun = subrun;
    f_event = event;
    f_fileName = fileName;
    tDataTree->Fill();
    printf("-------------------------------------------------------\n\n");
  }
} // END function analyze

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(FindFileWithEvent)

#endif // END def FindFileWithEvent_module
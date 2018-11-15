#ifndef GetPotCount_module
#define GetPotCount_module

#include "AnaHelper.h"

// Analyzer class
class GetPotCount : public art::EDAnalyzer
{
public:
  explicit GetPotCount(fhicl::ParameterSet const & pset);
  virtual ~GetPotCount();
  void analyze(art::Event const & evt) ;
  void endSubRun(art::SubRun const & sr);
  void respondToOpenInputFile(art::FileBlock const& fb);
  void beginJob();
  void endJob();
private:
  // Declare trees
  TTree *tPotCount;
  //TTree *tEventCount;
  bool fVerbose;
  bool fIsOverlayData;

  // Declare analysis variables
  int run, subrun, nEvents;
  std::vector<int> events;
  bool isRealData;
  double pot, mc_pot, data_pot;
  double mc_potSum = 0;
  double data_potSum = 0;

  // Declare analysis functions
  void ClearData();
}; // End class GetPotCount

GetPotCount::GetPotCount(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fVerbose(pset.get<bool>("verbose")),
    fIsOverlayData(pset.get<bool>("isOverlayData"))
{} // END constructor GetPotCount

GetPotCount::~GetPotCount()
{} // END destructor GetPotCount

void GetPotCount::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;
  tPotCount = tfs->make<TTree>("PotCount","");
  tPotCount->Branch("run",&run,"run/I");
  tPotCount->Branch("subrun",&subrun,"subrun/I");
  tPotCount->Branch("events",&events);
  tPotCount->Branch("nEvents",&nEvents);
  tPotCount->Branch("pot",&pot);
  tPotCount->Branch("mcPot",&mc_pot);
  tPotCount->Branch("dataPot",&data_pot);

} // END function beginJob

void GetPotCount::endJob()
{
} // END function endJob

void GetPotCount::ClearData()
{
  run = -1;
  subrun = -1;
  nEvents = -1;
} // END function ClearData

void GetPotCount::respondToOpenInputFile(art::FileBlock const &fb)
{
    mc_potSum = 0.;
    data_potSum = 0.;
}

void GetPotCount::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  events.push_back(evt.id().event());
  isRealData = evt.isRealData();
} // END function analyze

void GetPotCount::endSubRun(art::SubRun const & sr) 
{ 
  art::Handle<sumdata::POTSummary> potsum_h;

  // For OVERLAY (special case)
  if (fIsOverlayData)
  {
    if (sr.getByLabel("generator", potsum_h)) mc_pot = potsum_h->totpot - mc_potSum;
    else mc_pot = 0.;
    if (sr.getByLabel("beamdata", "bnbETOR860", potsum_h)) data_pot = potsum_h->totpot - data_potSum;
    else data_pot = 0.;

    mc_potSum += mc_pot;
    data_potSum += data_pot;
    pot = mc_potSum;
  }

  // For regular DATA and MONTE CARLO
  else
  {
    if (sr.getByLabel("generator", potsum_h)) mc_pot = potsum_h->totpot;
    else mc_pot = 0.;
    if (sr.getByLabel("beamdata", "bnbETOR860", potsum_h)) data_pot = potsum_h->totpot;
    else data_pot = 0.;
    // For MONTE CARLO
    if (!isRealData)
    {
      pot = mc_pot;
    }
    else
    {
      pot = data_pot;
    }

  } // END if isOverlayData

  std::cout << "----------------------------" << std::endl;
  std::cout << "Total POT / subRun: " << pot << std::endl;
  std::cout << "----------------------------" << std::endl;

  nEvents = events.size();
  tPotCount->Fill();
  events.clear();
} // END function endSubRun

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(GetPotCount)

#endif // END def GetPotCount_module

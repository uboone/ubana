#ifndef GetNormalizationHistograms_module
#define GetNormalizationHistograms_module

#include "AnaHelper.h"

// Analyzer class
class GetNormalizationHistograms : public art::EDAnalyzer
{
public:
  explicit GetNormalizationHistograms(fhicl::ParameterSet const & pset);
  virtual ~GetNormalizationHistograms();
  void analyze(art::Event const & evt) ;
  void beginJob();
  void endJob();
private:
  // Declare trees
  TTree *tNorm;
  std::string fHitLabel, fTrackLabel, fShowerLabel, fBeamFlashLabel, fCosmicFlashLabel;
  bool fExtraInformation, fVerbose;

  // Declare analysis variables
  int run, subrun, event;
  int nHits, nTracks, nShowers, nBeamFlashes, nCosmicFlashes;
  std::vector<double> trackLength, showerStartZ, beamFlashTime, beamFlashTotPE, cosmicFlashTime, cosmicFlashTotPE;

  // Declare analysis functions
  void ClearData();
}; // End class GetNormalizationHistograms

GetNormalizationHistograms::GetNormalizationHistograms(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fTrackLabel(pset.get<std::string>("TrackLabel")),
    fShowerLabel(pset.get<std::string>("ShowerLabel")),
    fBeamFlashLabel(pset.get<std::string>("BeamFlashLabel")),
    fCosmicFlashLabel(pset.get<std::string>("CosmicFlashLabel")),
    fExtraInformation(pset.get<bool>("ExtraInformation")),
    fVerbose(pset.get<bool>("Verbose"))

{} // END constructor GetNormalizationHistograms

GetNormalizationHistograms::~GetNormalizationHistograms()
{} // END destructor GetNormalizationHistograms

void GetNormalizationHistograms::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;
  tNorm = tfs->make<TTree>("Norm","");
  tNorm->Branch("run",&run,"run/I");
  tNorm->Branch("subrun",&subrun,"subrun/I");
  tNorm->Branch("event",&event,"event/I");
  tNorm->Branch("nHits",&nHits);
  tNorm->Branch("nTracks",&nTracks);
  tNorm->Branch("nShowers",&nShowers);
  tNorm->Branch("nBeamFlashes",&nBeamFlashes);
  tNorm->Branch("nCosmicFlashes",&nCosmicFlashes);
  if (fExtraInformation)
  {
    tNorm->Branch("trackLength",&trackLength);
    tNorm->Branch("showerStartZ",&showerStartZ);
    tNorm->Branch("beamFlashTime",&beamFlashTime);
    tNorm->Branch("beamFlashTotPE",&beamFlashTotPE);
    tNorm->Branch("cosmicFlashTime",&cosmicFlashTime);
    tNorm->Branch("cosmicFlashTotPE",&cosmicFlashTotPE);
  }

} // END function beginJob

void GetNormalizationHistograms::endJob()
{
} // END function endJob

void GetNormalizationHistograms::ClearData()
{
  trackLength.clear();
  showerStartZ.clear();
  beamFlashTime.clear();
  beamFlashTotPE.clear();
  cosmicFlashTime.clear();
  cosmicFlashTotPE.clear();
} // END function ClearData

void GetNormalizationHistograms::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.

  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  // Get tags and handles for all data products
  art::InputTag hitTag {fHitLabel};
  art::InputTag trackTag {fTrackLabel};
  art::InputTag showerTag {fShowerLabel};
  art::InputTag beamFlashTag {fBeamFlashLabel};
  art::InputTag cosmicFlashTag {fCosmicFlashLabel};

  const auto& hitHandle = evt.getValidHandle< std::vector<recob::Hit> >(hitTag);
  const auto& trackHandle = evt.getValidHandle< std::vector<recob::Track> >(trackTag);
  const auto& showerHandle = evt.getValidHandle< std::vector<recob::Shower> >(showerTag);
  const auto& beamFlashHandle = evt.getValidHandle< std::vector<recob::OpFlash> >(beamFlashTag);
  const auto& cosmicFlashHandle = evt.getValidHandle< std::vector<recob::OpFlash> >(cosmicFlashTag);

  // Get number of elements in each handle
  nHits = (*hitHandle).size();
  nTracks = (*trackHandle).size();
  nShowers = (*showerHandle).size();
  nBeamFlashes = (*beamFlashHandle).size();
  nCosmicFlashes = (*cosmicFlashHandle).size();

  // Get extra information about individual elements if needed
  if (fExtraInformation)
  {
    for(std::vector<int>::size_type i=0; i!=(*trackHandle).size(); i++)
    {
      art::Ptr<recob::Track> track(trackHandle,i);
      trackLength.push_back(track->Length());
    }
    for(std::vector<int>::size_type i=0; i!=(*showerHandle).size(); i++)
    {
      art::Ptr<recob::Shower> shower(showerHandle,i);
      showerStartZ.push_back(shower->ShowerStart()[2]);
    }
    for(std::vector<int>::size_type i=0; i!=(*beamFlashHandle).size(); i++)
    {
      art::Ptr<recob::OpFlash> flash(beamFlashHandle,i);
      beamFlashTime.push_back(flash->Time());
      beamFlashTotPE.push_back(flash->TotalPE());
    }
    for(std::vector<int>::size_type i=0; i!=(*cosmicFlashHandle).size(); i++)
    {
      art::Ptr<recob::OpFlash> flash(cosmicFlashHandle,i);
      cosmicFlashTime.push_back(flash->Time());
      cosmicFlashTotPE.push_back(flash->TotalPE());
    }
  }

  // Fill the tree
  tNorm->Fill();
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(GetNormalizationHistograms)

#endif // END def GetNormalizationHistograms_module
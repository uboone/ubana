#ifndef QUICKEXPLORER_MODULE
#define QUICKEXPLORER_MODULE

#include "QuickExplorer.h"

QuickExplorer::QuickExplorer(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor QuickExplorer
QuickExplorer::~QuickExplorer() {}
void QuickExplorer::beginJob() {}
void QuickExplorer::endJob() {}

// Example function to get handle and associations
void QuickExplorer::ExampleFunction(art::Event const & evt)
{
  // Get all the handles to objects.
  art::InputTag trackTag {std::string ("pandoraNu")};
  art::ValidHandle< std::vector<recob::Track> > const & trackHandle = evt.getValidHandle< std::vector<recob::Track> >(trackTag);
  art::FindManyP<anab::Calorimetry> caloTrackAss(trackHandle, evt, "pandoraNucalo");
  for(int i=0; i!=int((*trackHandle).size()); i++)
  {
    art::Ptr<recob::Track> track(trackHandle,i);
    std::vector<art::Ptr<anab::Calorimetry>> calos = caloTrackAss.at(i);
    printf("Track %i\n",i);
    printf("|_Length: %.1f\n", track->Length());
    for(int j=0; j!=int(calos.size()); j++)
    {
      printf("  |_CaloObj %i | PlaneID: %i | KE: %.1f | Range: %.1f\n", j, calos[j]->PlaneID().Plane, calos[j]->KineticEnergy(), calos[j]->Range());
    }
  }
}

void QuickExplorer::analyze(art::Event const & evt)
{
  // Determine event ID
  int run = evt.id().run();
  int subrun = evt.id().subRun();
  int event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

  ExampleFunction(evt);
} // END function analyze

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(QuickExplorer)

#endif // END def QuickExplorer_module

#ifndef FAKERAWDIGITS_MODULE
#define FAKERAWDIGITS_MODULE

#include "FakeRawDigits.h"

FakeRawDigits::FakeRawDigits(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor FakeRawDigits

FakeRawDigits::~FakeRawDigits()
{} // END destructor FakeRawDigits

void FakeRawDigits::beginJob()
{
} // END function beginJob

void FakeRawDigits::endJob()
{
} // END function endJob


// MAIN FUNCTION ANALYZE
 void FakeRawDigits::produce(art::Event& evt)
  {

} // END function produce



// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(FakeRawDigits)

#endif // END def FakeRawDigits_module

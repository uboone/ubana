////////////////////////////////////////////////////////////////////////
// Class:       WeightFilter
// Plugin Type: filter (art v3_01_02)
// File:        WeightFilter_module.cc
//
// Generated at Fri Feb 21 05:17:33 2020 by Wesley Ketchum using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <random>
#include <memory>

//include the truth objects
#include "larcoreobj/SummaryData/POTSummary.h"


namespace util {
  class WeightFilter;
}


class util::WeightFilter : public art::EDFilter {
public:
  explicit WeightFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WeightFilter(WeightFilter const&) = delete;
  WeightFilter(WeightFilter&&) = delete;
  WeightFilter& operator=(WeightFilter const&) = delete;
  WeightFilter& operator=(WeightFilter&&) = delete;

  bool filter(art::Event& e) override;
  bool beginSubRun(art::SubRun &) override;
  
private:

  art::InputTag fEventWeightTag;
  art::InputTag fPOTInfoTag;
  
  bool   fApplyPOTWeight;
  double fTargetPOT;
  double fPOTWeight;

  bool   fApplyExternalWeight;
  double fExternalWeight;

  bool   fApplyEventWeight;
  double fEventWeight;

  bool   fVerbose;

  std::default_random_engine fRandomGenerator;
  std::uniform_real_distribution<double> fFlatDistribution;
};


util::WeightFilter::WeightFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}, 
  fApplyPOTWeight(p.get<bool>("ApplyPOTWeight")),
  fTargetPOT(p.get<double>("TargetPOT",10.1E20)),
  fApplyExternalWeight(p.get<bool>("ApplyExternalWeight")),
  fExternalWeight(p.get<double>("ExternalWeight",1.0)),
  fApplyEventWeight(p.get<bool>("ApplyEventWeight")),
  fVerbose(p.get<bool>("Verbose",false)),
  fFlatDistribution(0.0,1.0)
{
    
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  
  if(fApplyEventWeight)
    fEventWeightTag = p.get<art::InputTag>("EventWeightTag");
  
  if(fApplyPOTWeight)
    fPOTInfoTag = p.get<art::InputTag>("POTInfoTag");
}

bool util::WeightFilter::beginSubRun(art::SubRun & s)
{
  fPOTWeight=1.;
  if(fApplyPOTWeight){
    auto const& potsum_handle = s.getValidHandle< sumdata::POTSummary >(fPOTInfoTag);
    fPOTWeight = fTargetPOT / potsum_handle->totpot;
    
    if(fVerbose) std::cout << "\tPOTWeight = " << fTargetPOT << " / " << potsum_handle->totpot 
			   << " = " << fPOTWeight << std::endl;
  }

  return true;
}

bool util::WeightFilter::filter(art::Event& )
{
  //get a random number
  double random = fFlatDistribution(fRandomGenerator);

  double fEventWeight=1.;
  if(fApplyEventWeight){
    //...
  }
  
  //calc total weight
  double weightTot=1.;
  if(fApplyPOTWeight)
    weightTot*=fPOTWeight;
  if(fApplyExternalWeight)
    weightTot*=fExternalWeight;
  if(fApplyEventWeight)
    weightTot*=fEventWeight;


  bool keepEvent = (random < weightTot);

  if(fVerbose) {
    std::cout << "Random number: " << random << std::endl;
    std::cout << "Weights: total = " 
	      << fPOTWeight << " * " << fExternalWeight << " * " << fEventWeight
	      << " = " << weightTot << std::endl;
    std::cout << "Keep event? " << keepEvent << std::endl;
  }

  return keepEvent;
}

DEFINE_ART_MODULE(util::WeightFilter)
    

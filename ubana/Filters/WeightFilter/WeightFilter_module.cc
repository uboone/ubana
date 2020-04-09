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
#include "larsim/EventWeight/Base/MCEventWeight.h"


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
  bool endSubRun(art::SubRun &) override;
  
private:

  std::vector<art::InputTag> fEventWeightTags;
  art::InputTag fPOTInfoTag;
  
  bool   fApplyPOTWeight;
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
  fPOTWeight(p.get<double>("POTWeight",1.0)),
  fApplyExternalWeight(p.get<bool>("ApplyExternalWeight")),
  fExternalWeight(p.get<double>("ExternalWeight",1.0)),
  fApplyEventWeight(p.get<bool>("ApplyEventWeight")),
  fVerbose(p.get<bool>("Verbose",false)),
  fFlatDistribution(0.0,1.0)
{
    
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  this->produces<sumdata::POTSummary,art::InSubRun>();
  
  if(fApplyEventWeight)
    fEventWeightTags = p.get< std::vector<art::InputTag> >("EventWeightTags");
  
  //if(fApplyPOTWeight)
  fPOTInfoTag = p.get<art::InputTag>("POTInfoTag");
}

bool util::WeightFilter::beginSubRun(art::SubRun & s)
{
  /*
  fPOTWeight=1.;
  if(fApplyPOTWeight){
    auto const& potsum_handle = s.getValidHandle< sumdata::POTSummary >(fPOTInfoTag);
    fPOTWeight = fTargetPOT / potsum_handle->totpot;
    
    if(fVerbose) std::cout << "\tPOTWeight = " << fTargetPOT << " / " << potsum_handle->totpot 
			   << " = " << fPOTWeight << std::endl;
  }
  */
  return true;
}

bool util::WeightFilter::endSubRun(art::SubRun & s)
{
  auto const& potsum_handle = s.getValidHandle< sumdata::POTSummary >(fPOTInfoTag);
  
  std::unique_ptr<sumdata::POTSummary> srpot_ptr(new sumdata::POTSummary());
  srpot_ptr->totpot = potsum_handle->totpot*fPOTWeight;
  srpot_ptr->totgoodpot = potsum_handle->totgoodpot*fPOTWeight;
  srpot_ptr->totspills = (int)(potsum_handle->totspills*fPOTWeight);
  srpot_ptr->goodspills = (int)(potsum_handle->goodspills*fPOTWeight);


  if(fVerbose) std::cout << "\tPOTWeight = " << srpot_ptr->totpot << " / " << potsum_handle->totpot 
			 << " = " << fPOTWeight << std::endl;

  s.put(std::move(srpot_ptr));
  
  return true;
}

bool util::WeightFilter::filter(art::Event& e)
{
  //get a random number
  double random = fFlatDistribution(fRandomGenerator);

  double fEventWeight=1.;
  if(fApplyEventWeight){
    //...
    for(auto const& tag : fEventWeightTags){
    
      auto const& evw_handle = e.getValidHandle< std::vector<evwgh::MCEventWeight> >(tag);
      auto const& evw_vec(*evw_handle);
      for(auto const& evw_map : evw_vec){
	for(auto const& evw : evw_map.fWeight){
	  
	  if(evw.second.size()>0) fEventWeight*=evw.second.at(0);	
	  if(fVerbose)
	    std::cout << "\t\tweight " << tag << " ... " << evw.first << " : " << evw.second.at(0) << std::endl;
	  
	}
      }
    }
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
    

////////////////////////////////////////////////////////////////////////
// Class:       OpFlashFilter
// Plugin Type: filter (art v2_05_00)
// File:        OpFlashFilter_module.cc
//
// Generated at Tue Mar 21 19:02:12 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
//
// This is a more generalized "OpFlash" filtering module to look at OpFlash
// collections and filter on them in various ways.
//
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

#include <memory>
#include <limits>
#include <bitset>

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"

namespace op {
  class OpFlashFilter;
}


class op::OpFlashFilter : public art::EDFilter {
public:
  explicit OpFlashFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpFlashFilter(OpFlashFilter const &) = delete;
  OpFlashFilter(OpFlashFilter &&) = delete;
  OpFlashFilter & operator = (OpFlashFilter const &) = delete;
  OpFlashFilter & operator = (OpFlashFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  art::InputTag fOpFlashTag;

  typedef struct {  

    double        MinPE;
    double        MaxPE;
    double        MinTime;
    double        MaxTime;
    
    unsigned int  MinMultiplicity;
    unsigned int  MaxMultiplicity;
    double        MultiplicityMaxPE;
    double        MultiplicityMinPE;

    bool          Invert;
  } FlashCutCollection_t;

  std::vector< std::vector<FlashCutCollection_t> > fCuts;
  bool fVerbose;

  bool PassesCuts(recob::OpFlash const&, FlashCutCollection_t const&, size_t);
  void PrintCuts(FlashCutCollection_t const&);
  
};


op::OpFlashFilter::OpFlashFilter(fhicl::ParameterSet const & p)
{
  this->reconfigure(p);
}

void op::OpFlashFilter::PrintCuts(op::OpFlashFilter::FlashCutCollection_t const& cut)
{
  std::cout << "\tMinPE:\t" << cut.MinPE << std::endl;
  std::cout << "\tMaxPE:\t" << cut.MaxPE << std::endl;
  std::cout << "\tMinTime:\t" << cut.MinTime << std::endl;
  std::cout << "\tMaxTime:\t" << cut.MaxTime << std::endl;
  std::cout << "\tMinMultiplicity:\t" << cut.MinMultiplicity << std::endl;
  std::cout << "\tMaxMultiplicity:\t" << cut.MaxMultiplicity << std::endl;
  std::cout << "\tMultiplicityMinPE:\t" << cut.MultiplicityMinPE << std::endl;
  std::cout << "\tMultiplicityMaxPE:\t" << cut.MultiplicityMaxPE << std::endl;
  std::cout << "\tInvert:\t" << cut.Invert << std::endl;
}

bool op::OpFlashFilter::PassesCuts(recob::OpFlash const& flash,
				   op::OpFlashFilter::FlashCutCollection_t const& cut,
				   size_t N_OPDETS)
{
  bool pass=true;
  
  if(flash.TotalPE() < cut.MinPE) pass=false;
  else if(flash.TotalPE() > cut.MaxPE) pass=false;
  else if(flash.Time() < cut.MinTime) pass=false;
  else if(flash.Time() > cut.MaxTime) pass=false;
  else{
    size_t n_hits=0;
    for(size_t i_opdet=0; i_opdet < N_OPDETS; ++i_opdet)
      if(flash.PE(i_opdet) > cut.MultiplicityMinPE &&
	 flash.PE(i_opdet) < cut.MultiplicityMaxPE) n_hits++;
    if(n_hits < cut.MinMultiplicity) pass=false;
    else if(n_hits > cut.MaxMultiplicity) pass=false;
  }

  if(cut.Invert)
    pass=!pass;
  
  if(pass && fVerbose){
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Passed cut" << std::endl;
    PrintCuts(cut);
  }
  
  return pass;
}

bool op::OpFlashFilter::filter(art::Event & e)
{
  auto const& opflash_handle = e.getValidHandle< std::vector<recob::OpFlash> >(fOpFlashTag);
  auto const& opflashVector(*opflash_handle);

  art::ServiceHandle<geo::Geometry> geo;
  
  
  bool pass_any_collection=false;
  for(auto const& this_cutcollection : fCuts){

    bool pass_all_cuts=true;
    for(auto const& this_cut : this_cutcollection){

      bool any_flash_pass=false;
      for(auto const& flash : opflashVector){
	if(PassesCuts(flash,this_cut,geo->NOpDets())){
	  any_flash_pass=true;
	  break;
	}
      }

      if(any_flash_pass==false){
	pass_all_cuts=false;
	break;
      }

    }

    if(pass_all_cuts==true){
      pass_any_collection=true;
      break;
    }
    
  }

  return pass_any_collection;

}

void op::OpFlashFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fOpFlashTag = p.get<art::InputTag>("OpFlashTag");
  fVerbose = p.get<bool>("Verbose",false);

  auto const& cutcollection_psets = p.get< std::vector< std::vector<fhicl::ParameterSet> > >("FlashCuts");
  fCuts.resize(cutcollection_psets.size());
  for(size_t i_col=0; i_col<fCuts.size(); ++i_col){
    auto const& psets = cutcollection_psets[i_col];
    auto & cutcollection = fCuts[i_col];

    cutcollection.resize(psets.size());
    for(size_t i_c=0; i_c<cutcollection.size(); ++i_c){

      auto & this_cut = cutcollection[i_c];
      auto const& this_pset = psets[i_c];
      
      this_cut.MinPE = this_pset.get<double>("MinPE",std::numeric_limits<double>::lowest());
      this_cut.MaxPE = this_pset.get<double>("MaxPE",std::numeric_limits<double>::max());
      
      this_cut.MinTime = this_pset.get<double>("MinTime",std::numeric_limits<double>::lowest());
      this_cut.MaxTime = this_pset.get<double>("MaxTime",std::numeric_limits<double>::max());
      
      this_cut.MinMultiplicity = this_pset.get<unsigned int>("MinMultiplicity",std::numeric_limits<unsigned int>::lowest());
      this_cut.MaxMultiplicity = this_pset.get<unsigned int>("MaxMultiplicity",std::numeric_limits<unsigned int>::max());
      this_cut.MultiplicityMinPE = this_pset.get<double>("MultiplicityMinPE",std::numeric_limits<double>::lowest());
      this_cut.MultiplicityMaxPE = this_pset.get<double>("MultiplicityMaxPE",std::numeric_limits<double>::max());
      
      this_cut.Invert = this_pset.get<bool>("Invert",false);

      if(fVerbose){
	std::cout << "=========================================" << std::endl;
	PrintCuts(this_cut);
      }
    }
  }

  

}

DEFINE_ART_MODULE(op::OpFlashFilter)

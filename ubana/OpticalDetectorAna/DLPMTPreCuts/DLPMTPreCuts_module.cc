////////////////////////////////////////////////////////////////////////
// Class:       DLPMTPreCuts
// Plugin Type: filter (art v2_06_03)
// File:        DLPMTPreCuts_module.cc
//
// Generated at Thu Apr 13 21:45:28 2017 by Taritree Wongjirad using cetskelgen
// from cetlib version v2_03_00.
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

#include "lardataobj/RecoBase/OpHit.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "larcore/Geometry/Geometry.h"

#include "ubevt/Database/LightYieldService.h"
#include "ubevt/Database/LightYieldProvider.h"
#include "ubevt/Database/UbooneLightYieldProvider.h"

// from dlpmtprecutalgo in UPS
#include "LEEPreCutAlgo.h"

#include <memory>

namespace dl {
  class DLPMTPreCuts;
}


class dl::DLPMTPreCuts : public art::EDFilter {
public:
  explicit DLPMTPreCuts(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DLPMTPreCuts(DLPMTPreCuts const &) = delete;
  DLPMTPreCuts(DLPMTPreCuts &&) = delete;
  DLPMTPreCuts & operator = (DLPMTPreCuts const &) = delete;
  DLPMTPreCuts & operator = (DLPMTPreCuts &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.
  leeprecuts::LEEPreCutAlgo m_algo;
  std::string fOpHitProducer;
  int fBinTickWidth;
  int fWinStartTick;
  int fWinEndTick;
  int fVetoStartTick;
  int fVetoEndTick;
  float fPEThreshold, fPEThresholdScaled;    
  float fPMTMaxFrac;
  bool  fStoreOpticalFilterObj;

  bool fScaleLY;

  std::set<int> fIgnoreChannelList;
};


dl::DLPMTPreCuts::DLPMTPreCuts(fhicl::ParameterSet const & p)
  : EDFilter{p}
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  fStoreOpticalFilterObj = p.get<bool>("StoreOpticalFilterObj",true);
  if(fStoreOpticalFilterObj)
    produces<uboone::UbooneOpticalFilter>();
  
  fOpHitProducer = p.get< std::string >("OpHitProducer");
  fBinTickWidth  = p.get< int >("BinTickWidth",     6);
  fWinStartTick  = p.get< int >("WinStartTick",   190);
  fWinEndTick    = p.get< int >("WinEndTick",     320);
  fPEThreshold   = p.get< float >("PEThreshold", 20.0);
  fVetoStartTick = p.get< int >("VetoStartTick",  60);
  fVetoEndTick   = p.get< int >("VetoEndTick",    190);
  fPMTMaxFrac    = p.get< float > ("PMTMaxFrac",  0.6);

  // Addition to account for LY scaling
  fScaleLY      = p.get< bool > ("ScaleLY",false);
  fPEThresholdScaled = fPEThreshold;


  std::vector<int> tmpvec = p.get< std::vector<int> >("IgnoreChannelList",std::vector<int>());
  fIgnoreChannelList = std::set<int>(tmpvec.begin(),tmpvec.end());

}

bool dl::DLPMTPreCuts::filter(art::Event & e)
{
  
  // if we are to scale by the LY, adjust fPEThreshold appropriately
  if (fScaleLY == true) {
    
    float scaling = 0.;
    int nfactors = 0;

    const ::art::ServiceHandle<geo::Geometry> geo;

    const lariov::LightYieldProvider& ly_provider = art::ServiceHandle<lariov::LightYieldService>()->GetProvider();
    for (unsigned int i=0; i!= geo->NOpDets(); ++i) {
      if (geo->IsValidOpChannel(i) && i<32) {
	scaling += ly_provider.LYScaling(i);
	nfactors += 1;
      }
    }
    
    if (nfactors != 32) std::cout << "ERROR != 32 PMTs " << std::endl;
    
    scaling /= nfactors;

    fPEThresholdScaled = fPEThreshold * scaling;
    
    std::cout << "scaling now is " << fPEThresholdScaled << std::endl;
  }

  // Implementation of required member function here.
  art::Handle< std::vector<recob::OpHit> > ophitHandle;
  e.getByLabel( fOpHitProducer, ophitHandle );
  std::vector<recob::OpHit> ophit_v(*ophitHandle);
  
  std::vector<float> ophit_peaktime_v( ophit_v.size(), 0.0);
  std::vector<float> ophit_pe_v( ophit_v.size(), 0.0);
  std::vector<int>   ophit_femch_v( ophit_v.size(), 0);
  
  for ( size_t i=0; i<ophit_v.size(); i++ ) {
    auto const& ophit = ophit_v.at(i);
    
    if(fIgnoreChannelList.count(ophit.OpChannel())) {
      continue;
    }

    ophit_peaktime_v[i] = ophit.PeakTime();
    ophit_pe_v[i] = ophit.PE();
    ophit_femch_v[i] = ophit.OpChannel();
  }
  
  std::vector<float> flashbins = m_algo.MakeTimeBin( ophit_peaktime_v, ophit_pe_v, fBinTickWidth, fWinStartTick, fWinEndTick );
  std::vector<float> vetobins  = m_algo.MakeTimeBin( ophit_peaktime_v, ophit_pe_v, fBinTickWidth, fVetoStartTick, fVetoEndTick );
  
  std::vector<float> beamPEinfo = m_algo.GetTotalPE( fPEThresholdScaled , flashbins );
  std::vector<float> vetoPEinfo = m_algo.GetTotalPE( fPEThresholdScaled , vetobins );
  
  float maxfrac     = m_algo.PMTMaxFrac( ophit_peaktime_v, ophit_pe_v, ophit_femch_v, beamPEinfo, fBinTickWidth,  fWinStartTick);

  if(fStoreOpticalFilterObj)
    {
      float total_beam_PE=0; for(auto const& pe : flashbins) total_beam_PE+= pe;
      float total_veto_PE=0; for(auto const& pe : vetobins) total_veto_PE+= pe;
      std::unique_ptr<uboone::UbooneOpticalFilter> ubopfilter_obj( new uboone::UbooneOpticalFilter(beamPEinfo[0],vetoPEinfo[0],maxfrac,
												   total_beam_PE,total_veto_PE));
      e.put(std::move(ubopfilter_obj));
    }
  
  if ( beamPEinfo[0]>fPEThresholdScaled && vetoPEinfo[0]<fPEThresholdScaled && maxfrac < fPMTMaxFrac )
    return true;
  
  return false;  
}

DEFINE_ART_MODULE(dl::DLPMTPreCuts)

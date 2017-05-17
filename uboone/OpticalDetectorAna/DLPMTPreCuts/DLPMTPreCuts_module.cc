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
  float fPEThreshold;    
  float fPMTMaxFrac;

};


dl::DLPMTPreCuts::DLPMTPreCuts(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.

  fOpHitProducer = p.get< std::string >("OpHitProducer");
  fBinTickWidth  = p.get< int >("BinTickWidth",     6);
  fWinStartTick  = p.get< int >("WinStartTick",   190);
  fWinEndTick    = p.get< int >("WinEndTick",     320);
  fPEThreshold   = p.get< float >("PEThreshold", 20.0);
  fVetoStartTick = p.get< int >("VetoStartTick",  60);
  fVetoEndTick   = p.get< int >("VetoEndTick",    190);
  fPMTMaxFrac    = p.get< float > ("PMTMaxFrac",  0.6);

}

bool dl::DLPMTPreCuts::filter(art::Event & e)
{
  
  // Implementation of required member function here.
  art::Handle< std::vector<recob::OpHit> > ophitHandle;
  e.getByLabel( fOpHitProducer, ophitHandle );
  std::vector<recob::OpHit> ophit_v(*ophitHandle);
  
  std::vector<float> ophit_peaktime_v( ophit_v.size(), 0.0);
  std::vector<float> ophit_pe_v( ophit_v.size(), 0.0);
  std::vector<int>   ophit_femch_v( ophit_v.size(), 0);
  
  for ( size_t i=0; i<ophit_v.size(); i++ ) {
    auto const& ophit = ophit_v.at(i);
    ophit_peaktime_v[i] = ophit.PeakTime();
    ophit_pe_v[i] = ophit.PE();
    ophit_femch_v[i] = ophit.OpChannel();
  }
  
  std::vector<float> flashbins = m_algo.MakeTimeBin( ophit_peaktime_v, ophit_pe_v, fBinTickWidth, fWinStartTick, fWinEndTick );
  std::vector<float> vetobins  = m_algo.MakeTimeBin( ophit_peaktime_v, ophit_pe_v, fBinTickWidth, fVetoStartTick, fVetoEndTick );
  
  std::vector<float> beamPEinfo = m_algo.GetTotalPE( fPEThreshold , flashbins );
  std::vector<float> vetoPEinfo = m_algo.GetTotalPE( fPEThreshold , vetobins );
  
  float maxfrac     = m_algo.PMTMaxFrac( ophit_peaktime_v, ophit_pe_v, ophit_femch_v, beamPEinfo, fBinTickWidth,  fWinStartTick);
  
  if ( beamPEinfo[0]>fPEThreshold && vetoPEinfo[0]<fPEThreshold && maxfrac < fPMTMaxFrac )
    return true;
  
  return false;  
}

DEFINE_ART_MODULE(dl::DLPMTPreCuts)

#ifndef __SubEventModule__
#define __SubEventModule__

#include <vector>
#include <map>
#include "SubEventModConfig.hh"
#include "ubooneobj/Optical/SubEvent.hh"
#include "ubooneobj/Optical/Flash.hh"
#include "ubooneobj/Optical/FlashList.hh"
#include "ubooneobj/Optical/SubEventList.hh"
#include "WaveformData.hh"

namespace subevent {
  
  // Main Routine ------------------------------------------------------------
  void formSubEvents( WaveformData& wfms, SubEventModConfig& config, std::map< int, double >& pmtspemap, SubEventList& subevents, FlashList& unclaimed_flashes );
  // -------------------------------------------------------------------------

  int findChannelFlash( int ch, std::vector<double>& waveform, SubEventModConfig& config, std::string discrname, Flash& returned_flash );
  int getChannelFlashes( int channel, std::vector< double >& waveform, std::vector< double >& baseline, SubEventModConfig& config, std::string discrname,FlashList& flashes, std::vector<double>& postwfm );
  
  void formFlashes( WaveformData& wfms, SubEventModConfig& config, std::string discrname, FlashList& flashes, WaveformData& postwfms );
  void fillFlashAccumulators( FlashList& flashes, std::map< int, double >& pmtspemap, SubEventModConfig& config, std::vector< double >& peacc, std::vector< double >& hitacc );

  void AnalyzeSubEvents( SubEventList& subevents );
  
}

#endif

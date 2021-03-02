#ifndef OPTICALDETECTORANA_OPDETWAVEFORMANA_CXX
#define OPTICALDETECTORANA_OPDETWAVEFORMANA_CXX

#include "OpDetWaveformAna.h"
#include <limits>
#include <climits>
#include <iostream>

namespace pmtana {

  OpDetWaveformAna::OpDetWaveformAna(const std::string name)
    : _name(name)
    , _hitana_tree   (nullptr)
    , _wfana_tree    (nullptr)
    , _wf_tree       (nullptr)
    , _ev_wf_tree    (nullptr)
    , _ev_hit_tree   (nullptr)
    , _ev_flash_tree (nullptr)
  {
    ClearEvent();
    ClearWaveform();
    _preco_mgr.AddRecoAlgo(&_preco_alg);

    //_preco_mgr.SetPedAlgo(pmtana::kHEAD);
    //_preco_mgr.SePedSampleCosmic (  3 );
    //_preco_mgr.SetPedSampleBeam   ( 10 );
  }

  void OpDetWaveformAna::TickPeriod   ( const double period )
  { _period = period; }
  
  void OpDetWaveformAna::SetEventInfo ( const unsigned int run,
					const unsigned int subrun,
					const unsigned int event  )
  {
    _run    = run;
    _subrun = subrun;
    _event  = event;
  }
  
  void OpDetWaveformAna::AnaWaveform  ( const unsigned int ch,
					const double time_wrt_trigger,
					const std::vector<short>& wf)
  {

    if( !_hitana_tree && !_wfana_tree && !_wf_tree ) return;
    ClearWaveform();
    
    _ch = ch;
    _t_wrt_trigger = time_wrt_trigger;

    if( _wfana_tree || _wf_tree ) {
      for(auto const& adc : wf) {
	if(adc > _max_adc) _max_adc = adc;
	if(adc < _min_adc) _min_adc = adc;
      }
      _wf_size = wf.size();


      if( _wf_tree ) _wf = wf;
    }

    if( _wfana_tree  ) _wfana_tree->Fill();
    if( _wf_tree     ) _wf_tree->Fill();

  }

  void OpDetWaveformAna::AnaEventWaveform  ( const std::vector<raw::OpDetWaveform>& ev_wf_v)
  {

    if( !_ev_wf_tree ) return;
    ClearWaveform();
    
    for (size_t i=0; i < ev_wf_v.size(); i++) {

      auto const wf = ev_wf_v[i];

      auto ch = wf.ChannelNumber();

      if (ch >= 32) continue;

      
      for (size_t n=0; n < wf.size(); n++){
	auto adc = wf[n];
	if(adc > _max_adc_v[ch]) _max_adc_v[ch] = adc;
	if(adc < _min_adc_v[ch]) _min_adc_v[ch] = adc;
	if (n < 1500){
	  _wf_v[ch][n] = ev_wf_v[i][n];
	  _wfsum[n] += wf[n];
	}
      }// for all channels

    }// for all waveforms      
    
    if( _ev_wf_tree ) _ev_wf_tree->Fill();
    
    return;
  }


  void OpDetWaveformAna::AnaEventHit  ( const std::vector<recob::OpHit>& ev_hit_v)
  {

    if( !_ev_hit_tree ) return;
    
    for (size_t i=0; i < ev_hit_v.size(); i++) {

      auto const hit = ev_hit_v[i];

      _hit_ch = hit.OpChannel();
      
      _hit_pe = hit.PE();
      _hit_area = hit.Area();
      _hit_ampl = hit.Amplitude();
      _hit_time = hit.PeakTime();

      if (_ev_hit_tree) _ev_hit_tree->Fill();
      
    }// for all hits
    
    return;
  }
  
  void OpDetWaveformAna::AnaEventFlash  ( const std::vector<recob::OpFlash>& ev_flash_v) {
    
    if( !_ev_flash_tree ) return;
    
    for (size_t i=0; i < ev_flash_v.size(); i++) {
      
      auto const flash = ev_flash_v[i];
      
      _flash_time    = flash.Time();
      _flash_zcenter = flash.ZCenter();
      _flash_zwidth  = flash.ZWidth();
      _flash_ycenter = flash.YCenter();
      _flash_ywidth  = flash.YWidth();
      _flash_petotal = std::accumulate(flash.PEs().begin(),flash.PEs().end(),0);

      if (_ev_flash_tree) _ev_flash_tree->Fill();
      
    }// for all hits
    
    return;
  }
  
  void OpDetWaveformAna::ClearEvent()
  {
    _run = _subrun = _event = std::numeric_limits<unsigned int>::max();
  }

  void OpDetWaveformAna::ClearWaveform()
  {
    _ch = std::numeric_limits<unsigned int>::max();
    _wf_size = 0;

    _ped_mean = _ped_rms = -1;
    _ped_mean_v = std::vector<float>(32,-1);
    _ped_rms_v  = std::vector<float>(32,-1);
    _tstart = _tpeak = _tend = -1;
    _t_wrt_trigger = 0;
    _q = _amp = -1;
    _max_adc = 0;
    _max_adc_v = std::vector<unsigned short>(32,0);
    _min_adc = std::numeric_limits<unsigned short>::max();
    _min_adc_v =  std::vector<unsigned short>(32,std::numeric_limits<unsigned short>::max());
    _wf_size = 0;
    _wf.clear();
    _wfsum = std::vector<int>(1500,0);
    //_wf_v = std::vector<std::vector<short>>(32,std::vector<short>());
  }

  void OpDetWaveformAna::AnaHit       ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }     
    _hitana_tree = ptr;
    _hitana_tree->Branch( "run",      &_run,      "run/i"      );
    _hitana_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _hitana_tree->Branch( "event",    &_event,    "event/i"    );
    _hitana_tree->Branch( "ch",       &_ch,       "ch/i"       );
    _hitana_tree->Branch( "ped_mean", &_ped_mean, "ped_mean/F" );
    _hitana_tree->Branch( "ped_rms",  &_ped_rms,  "ped_rms/F"  );
    _hitana_tree->Branch( "tstart",   &_tstart,   "tstart/D"   );
    _hitana_tree->Branch( "tpeak",    &_tpeak,    "tpeak/D"    );
    _hitana_tree->Branch( "tend",     &_tend,     "tend/D"     );
    _hitana_tree->Branch( "q",        &_q,        "q/D"        );
    _hitana_tree->Branch( "amp",      &_amp,      "amp/D"      );

  }
  void OpDetWaveformAna::AnaWaveform  ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }
    _wfana_tree = ptr;
    _wfana_tree->Branch( "run",       &_run,      "run/i"      );
    _wfana_tree->Branch( "subrun",    &_subrun,   "subrun/i"   );
    _wfana_tree->Branch( "event",     &_event,    "event/i"    );
    _wfana_tree->Branch( "ch",        &_ch,       "ch/i"       );
    _wfana_tree->Branch( "ped_mean",  &_ped_mean, "ped_mean/F" );
    _wfana_tree->Branch( "ped_rms",   &_ped_rms,  "ped_rms/F"  );
    _wfana_tree->Branch( "max_adc",   &_max_adc,  "max_adc/s"  );
    _wfana_tree->Branch( "min_adc",   &_min_adc,  "min_adc/s"  );
    _wfana_tree->Branch( "wf_size",   &_wf_size,  "wf_size/i"  );

  }
  void OpDetWaveformAna::SaveWaveform ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }
    _wf_tree = ptr;
    _wf_tree->Branch( "run",      &_run,      "run/i"      );
    _wf_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _wf_tree->Branch( "event",    &_event,    "event/i"    );
    _wf_tree->Branch( "ch",       &_ch,       "ch/i"       );
    _wf_tree->Branch( "t_wrt_trigger", &_t_wrt_trigger, "t_wrt_trigger/D"   );
    _wf_tree->Branch( "ped_mean", &_ped_mean, "ped_mean/F" );
    _wf_tree->Branch( "ped_rms",  &_ped_rms,  "ped_rms/F"  );
    _wf_tree->Branch( "max_adc",  &_max_adc,  "max_adc/s"  );
    _wf_tree->Branch( "min_adc",  &_min_adc,  "min_adc/s"  );
    _wf_tree->Branch( "wf", "std::vector<short>", &_wf  );
  }

  void OpDetWaveformAna::SaveEvWaveform ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }

    _wf_v = std::vector< std::vector<short> >(32,std::vector<short>(1500,0));

    _ev_wf_tree = ptr;
    _ev_wf_tree->Branch( "run",      &_run,      "run/i"      );
    _ev_wf_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _ev_wf_tree->Branch( "event",    &_event,    "event/i"    );
    _ev_wf_tree->Branch( "ped_mean_v", &_ped_mean_v, "ped_mean_v/F" );
    _ev_wf_tree->Branch( "ped_rms_v",  &_ped_rms_v,  "ped_rms_v/F"  );
    _ev_wf_tree->Branch( "max_adc_v",  &_max_adc_v,  "max_adc_v/s"  );
    _ev_wf_tree->Branch( "min_adc_v",  &_min_adc_v,  "min_adc_v/s"  );
    _ev_wf_tree->Branch( "wfsum", "std::vector<int>", &_wfsum  );
    _ev_wf_tree->Branch( "wf_00", "std::vector<short>", &(_wf_v[0])  );
    _ev_wf_tree->Branch( "wf_01", "std::vector<short>", &(_wf_v[1])  );
    _ev_wf_tree->Branch( "wf_02", "std::vector<short>", &(_wf_v[2])  );
    _ev_wf_tree->Branch( "wf_03", "std::vector<short>", &(_wf_v[3])  );
    _ev_wf_tree->Branch( "wf_04", "std::vector<short>", &(_wf_v[4])  );
    _ev_wf_tree->Branch( "wf_05", "std::vector<short>", &(_wf_v[5])  );
    _ev_wf_tree->Branch( "wf_06", "std::vector<short>", &(_wf_v[6])  );
    _ev_wf_tree->Branch( "wf_07", "std::vector<short>", &(_wf_v[7])  );
    _ev_wf_tree->Branch( "wf_08", "std::vector<short>", &(_wf_v[8])  );
    _ev_wf_tree->Branch( "wf_09", "std::vector<short>", &(_wf_v[9])  );
    _ev_wf_tree->Branch( "wf_10", "std::vector<short>", &(_wf_v[10])  );
    _ev_wf_tree->Branch( "wf_11", "std::vector<short>", &(_wf_v[11])  );
    _ev_wf_tree->Branch( "wf_12", "std::vector<short>", &(_wf_v[12])  );
    _ev_wf_tree->Branch( "wf_13", "std::vector<short>", &(_wf_v[13])  );
    _ev_wf_tree->Branch( "wf_14", "std::vector<short>", &(_wf_v[14])  );
    _ev_wf_tree->Branch( "wf_15", "std::vector<short>", &(_wf_v[15])  );
    _ev_wf_tree->Branch( "wf_16", "std::vector<short>", &(_wf_v[16])  );
    _ev_wf_tree->Branch( "wf_17", "std::vector<short>", &(_wf_v[17])  );
    _ev_wf_tree->Branch( "wf_18", "std::vector<short>", &(_wf_v[18])  );
    _ev_wf_tree->Branch( "wf_19", "std::vector<short>", &(_wf_v[19])  );
    _ev_wf_tree->Branch( "wf_20", "std::vector<short>", &(_wf_v[20])  );
    _ev_wf_tree->Branch( "wf_21", "std::vector<short>", &(_wf_v[21])  );
    _ev_wf_tree->Branch( "wf_22", "std::vector<short>", &(_wf_v[22])  );
    _ev_wf_tree->Branch( "wf_23", "std::vector<short>", &(_wf_v[23])  );
    _ev_wf_tree->Branch( "wf_24", "std::vector<short>", &(_wf_v[24])  );
    _ev_wf_tree->Branch( "wf_25", "std::vector<short>", &(_wf_v[25])  );
    _ev_wf_tree->Branch( "wf_26", "std::vector<short>", &(_wf_v[26])  );
    _ev_wf_tree->Branch( "wf_27", "std::vector<short>", &(_wf_v[27])  );
    _ev_wf_tree->Branch( "wf_28", "std::vector<short>", &(_wf_v[28])  );
    _ev_wf_tree->Branch( "wf_29", "std::vector<short>", &(_wf_v[29])  );
    _ev_wf_tree->Branch( "wf_30", "std::vector<short>", &(_wf_v[30])  );
    _ev_wf_tree->Branch( "wf_31", "std::vector<short>", &(_wf_v[31])  );
  }

  void OpDetWaveformAna::SaveEvHit ( TTree* ptr ) {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }
    
    _ev_hit_tree = ptr;
    _ev_hit_tree->Branch( "run",      &_run,      "run/i"      );
    _ev_hit_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _ev_hit_tree->Branch( "event",    &_event,    "event/i"    );
    _ev_hit_tree->Branch("hit_ch",&_hit_ch, "hit_ch/I");
    _ev_hit_tree->Branch("hit_pe",&_hit_pe,"hit_pe/F");
    _ev_hit_tree->Branch("hit_time",&_hit_time,"hit_time/F");
    _ev_hit_tree->Branch("hit_ampl",&_hit_ampl,"hit_ampl/F");
    _ev_hit_tree->Branch("hit_area",&_hit_area,"hit_area/F");
  }
  
  void OpDetWaveformAna::SaveEvFlash ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }

    _ev_flash_tree = ptr;
    _ev_flash_tree->Branch( "run",      &_run,      "run/i"      );
    _ev_flash_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _ev_flash_tree->Branch( "event",    &_event,    "event/i"    );
    _ev_flash_tree->Branch("flash_time",&_flash_time,"flash_time/F");
    _ev_flash_tree->Branch("flash_petotal",&_flash_petotal,"flash_petotal/F");
    _ev_flash_tree->Branch("flash_zwidth",&_flash_zwidth,"flash_zwidth/F");
    _ev_flash_tree->Branch("flash_zcenter",&_flash_zcenter,"flash_zcenter/F");
    _ev_flash_tree->Branch("flash_ywidth",&_flash_ywidth,"flash_ywidth/F");
    _ev_flash_tree->Branch("flash_ycenter",&_flash_ycenter,"flash_ycenter/F");
  }


}



#endif

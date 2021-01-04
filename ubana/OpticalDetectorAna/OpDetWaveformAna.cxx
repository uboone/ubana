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

  void OpDetWaveformAna::SetEventGain ( const std::vector<double> opch_area_gain_v,
					const std::vector<double> opch_ampl_gain_v)
  {
    _opch_area_gain_v = opch_area_gain_v;
    _opch_ampl_gain_v = opch_ampl_gain_v;
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

      //auto ch = wf.ChannelNumber();

      if (i >= _wf_v.size()) continue;

      
      for (size_t n=0; n < wf.size(); n++){
	auto adc = wf[n];
	if(adc > _max_adc_v[i]) _max_adc_v[i] = adc;
	if(adc < _min_adc_v[i]) _min_adc_v[i] = adc;
	if (n < 1500){
	  _wf_v[i][n] = ev_wf_v[i][n];
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
      
      _hit_area = hit.Area();
      _hit_area_pe = _hit_area / _opch_area_gain_v[_hit_ch%100];
      _hit_ampl = hit.Amplitude();
      _hit_ampl_pe = _hit_ampl / _opch_ampl_gain_v[_hit_ch%100];
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

  void OpDetWaveformAna::SaveEvWaveform ( TTree* ptr, unsigned int nchan, bool savesum, int channelmask )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }

    _wf_v = std::vector< std::vector<short> >(nchan,std::vector<short>(1500,0));

    _ev_wf_tree = ptr;
    _ev_wf_tree->Branch( "run",      &_run,      "run/i"      );
    _ev_wf_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _ev_wf_tree->Branch( "event",    &_event,    "event/i"    );
    _ev_wf_tree->Branch( "ped_mean_v", &_ped_mean_v, "ped_mean_v/F" );
    _ev_wf_tree->Branch( "ped_rms_v",  &_ped_rms_v,  "ped_rms_v/F"  );
    _ev_wf_tree->Branch( "max_adc_v",  &_max_adc_v,  "max_adc_v/s"  );
    _ev_wf_tree->Branch( "min_adc_v",  &_min_adc_v,  "min_adc_v/s"  );

    for (unsigned int n=0; n < nchan; n++){ 
      if ( (channelmask >= 0) && ((unsigned int)channelmask == n) )
	_ev_wf_tree->Branch( TString::Format("wf_%02u",n), "std::vector<short>", &(_wf_v[n])  );
      if (channelmask < 0) // always save
	_ev_wf_tree->Branch( TString::Format("wf_%02u",n), "std::vector<short>", &(_wf_v[n])  );
    }
    if (savesum)
      _ev_wf_tree->Branch( "wfsum", "std::vector<int>", &_wfsum  );
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
    _ev_hit_tree->Branch("hit_area_pe",&_hit_area_pe,"hit_area_pe/F");
    _ev_hit_tree->Branch("hit_ampl_pe",&_hit_ampl_pe,"hit_ampl_pe/F");
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

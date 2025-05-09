////////////////////////////////////////////////////////////////////////
// Class:       OpDigitSaturationCorrection
// Module Type: producer
// File:        OpDigitSaturationCorrection_module.cc
//
// Generated at Sun Oct 25 16:07:45 2015 by David Caratelli using artmod
// from cetpkgsupport v1_08_07.

// DOCUMENTATION FOR THIS MODULE CAN BE FOUND @ DOCDB 4979
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=4979

////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "ubcore/Geometry/UBOpReadoutMap.h"

// data-products
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
//#include "RawData/TriggerData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// C++ includes
#include <memory>
#include <iostream>
#include <cmath>
#include <limits>

// ROOT includes
#include <TTree.h>

// function to get idx of max ADC from a wf
size_t getMaxADC(const std::vector<short>& wf){
  
  size_t max_idx = 0;
  short  max_adc = 0;
  for (size_t i=0; i < wf.size(); i++){
    if (wf.at(i) > max_adc){
      max_adc = wf.at(i);
      max_idx = i;
    }
  }// for all ADCs
  
  return max_idx;
}

class OpDigitSaturationCorrection;

class OpDigitSaturationCorrection : public art::EDProducer {
public:
  explicit OpDigitSaturationCorrection(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpDigitSaturationCorrection(OpDigitSaturationCorrection const &) = delete;
  OpDigitSaturationCorrection(OpDigitSaturationCorrection &&) = delete;
  OpDigitSaturationCorrection & operator = (OpDigitSaturationCorrection const &) = delete;
  OpDigitSaturationCorrection & operator = (OpDigitSaturationCorrection &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Verbosity setter
  void SetVerbose(bool on) { _verbose = on; }

public:
  
  struct LGPulseInfo_t {
    size_t _index;
    double _start;
    double _duration;
    bool IsUseful(double start, double start_diff, double duration=-1, double duration_diff=-1) const {
      if(std::abs(start - _start) > start_diff) return false; // Restrict start time difference
      if(duration>0 && std::abs(duration - _duration) > duration_diff) return false; // restrict duration diff (optional)
      return true;
    }
  };

  typedef std::vector<std::vector<LGPulseInfo_t> > LGChanMap_t;

private:

  // Declare member data here.

  // functions:
  size_t FindMatchingLGPulse(const unsigned int chan, const double time,
			     const LGChanMap_t& LG_ChanMap);

  void MakeTree();
  
  // sets up various mapping to figure out optical channel, detector id, and type mapping
  void SetUpChannelMap();

  // trigger time (usec)
  double _TrigTime;

  // time-tick duration in usec
  double _TDC;

  // time-tick count minimum to define "beamgate"
  double _BEAMGATE_MIN;
  
  // by how much to scale LG so that it matches the HG?
  double _gain_fact;

  // max ADC
  short unsigned int ADC_max;

  // max IDX
  size_t IDX_max;

  // verbosity flag
  bool _verbose;

  // Flag to include beamgate
  bool _include_beamgate;

  // baseline?
  short unsigned int _baseline;

  // vector storing channel number for calibration correction
  std::vector<unsigned int> _gain_chan_v;
  // vector holding the calibration correction factors per PMT
  std::vector<double>       _gain_vals_v;

  // vector storing channel number for baseline values
  std::vector<unsigned int> _baseline_chan_v;
  // vector storing baseline value
  std::vector<double>       _baseline_vals_v;

  // save waveforms in ttree
  bool _save_wfs;

  // flat to allow use LG Beam FEM for HG Cosmic correction 
  // If enabled, _LG_cosmic_* must be empty (else exception)
  // This handle is to coop with UB situation where only 3FEM employeed (09/03/2016)
  bool _use_LG_beam_for_HG_cosmic;

  // Waveform producer name
  std::string _LG_producer;
  std::string _LG_label;
  std::string _HG_producer;
  std::string _HG_label;
  std::string _LG_cosmic_producer;
  std::string _LG_cosmic_label;
  std::string _HG_cosmic_producer;
  std::string _HG_cosmic_label;

  // OpChannel => UBOpChannelType mapping
  // Note these are opdet::Undefined, opdet::HighGain, opdet::LowGain, and opdet::LogicChannel
  std::vector<opdet::UBOpticalChannelType_t> _opch_to_chtype_m; 

  // OpChannel => UBOpChannelCategory mapping
  // Note categories are:
  /*
    opdet::Uncategorized = 0
    opdet::UnspecifiedLogic
    opdet::OpdetCosmicHighGain
    opdet::OpdetCosmicLowGain
    opdet::OpdetBeamHighGain
    opdet::OpdetBeamLowGain
    opdet::BNBLogicPulse
    opdet::NUMILogicPulse
    opdet::FlasherLogicPulse
    opdet::StrobeLogicPulse
  */
  std::vector<opdet::UBOpticalChannelCategory_t> _opch_to_chcategory_m;

  // OpChannel => OpDet mapping
  std::vector<int> _opch_to_opdet_m;

  // OpDet => OpChannel
  std::vector<std::set<int> > _opdet_to_opch_m;

  // OpChannel LG => OpChannel HG on same FEM
  std::map<int,int> _opchLG_to_opchHG_m;
  
  // OpChannel => High/Low gain mapping
  std::vector<bool> _opch_to_highgain_m;

  // Ignore set of channels
  std::set<size_t> _ignore_channels;

  // OpChannel to baseline
  std::vector<double> _opch_to_baseline_m;
  
  // OpChannel to Gain
  std::vector<double> _opch_to_gain_m;

  // Tree variables
  TTree* _tree;
  int    _event;
  int    _subrun;
  int    _run;
  int    _chan_HG;   // OpChannel channel number HG
  int    _chan_LG;   // OpChannel channel number LG
  int    _nticks_HG, _nticks_LG;
  int    _pmt;       // OpDet PMT number
  int    _matched;   // boolean matched HG to LG
  int    _swapped;   // LG corrected swapped in for HG
  double _time_HG, _time_LG;
  int    _max_adc_HG, _max_adc_LG;
  int    _max_tdc_HG, _max_tdc_LG;
  std::vector<short> _wf_HG, _wf_LG, _wf_corr;
  int    _max_adc_corr; // the maximum ADC assigned to the LG corrected WF that replaced the HG saturated WF
};

size_t OpDigitSaturationCorrection::FindMatchingLGPulse(const unsigned int opdet,
							const double start,
							const LGChanMap_t& LG_ChanMap)
{

  // search the _LG_ChanMap for a pulse on this channel at approx. the same time
  if (LG_ChanMap.at(opdet).size() == 0)
    return IDX_max;

  // HG_pulses:
  auto const& LG_pulses = LG_ChanMap.at(opdet);
  
  // go through them and find one with a matching time
  for (auto const& pulse : LG_pulses) {
    if(pulse.IsUseful(start,(2*_TDC))) {
      // we have found a matching pulse!
      if (_verbose)
	std::cout << "found a matching pulse for opdet " << opdet << " @ time " << start << std::endl;
      return pulse._index;
    }
  }// for all HG pulses
  return IDX_max;
}

OpDigitSaturationCorrection::OpDigitSaturationCorrection(fhicl::ParameterSet const & p)
  : EDProducer{p}, _tree(nullptr)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  //produces< std::vector<raw::OpDetWaveform>   >();
  produces< std::vector<raw::OpDetWaveform> >( opdet::UBOpChannelEnumName( opdet::OpdetBeamHighGain ) );
  produces< std::vector<raw::OpDetWaveform> >( opdet::UBOpChannelEnumName( opdet::OpdetCosmicHighGain ) );

  // load per-PMT calibration correction factors
  _gain_chan_v = p.get< std::vector<unsigned int> > ( "GainChannelNumber" );
  _gain_vals_v = p.get< std::vector<double>       > ( "GainValues"        );
  if (_gain_chan_v.size() != _gain_vals_v.size() ){
    std::cerr << "Calibration Correction values and channel-list vectors do not agree in size...";
    throw std::exception();
  }

  // load per-channel baseline values
  _baseline_chan_v = p.get< std::vector<unsigned int> > ("BaselineChannelNumber");
  _baseline_vals_v = p.get< std::vector<double>       > ("BaselineValues"       );
  if (_baseline_chan_v.size() != _baseline_vals_v.size() ){
    std::cerr << "Baseline values and channel-list vectors do not agree in size...";
    throw std::exception();
  }
 
  // Take care of a beamgate if specified
  _include_beamgate = p.get<bool>("IncludeBeamgate");

  // time tick is 15.6 ns (in usec)
  _TDC = 0.0156;

  // beam gate tick count min
  _BEAMGATE_MIN = 800;

  // scale LG to HG
  _gain_fact = 10.;

  // set default baseline
  _baseline = 2048;

  // max ADC
  ADC_max = std::numeric_limits<short unsigned int>::max();

  // max IDX
  IDX_max = std::numeric_limits<size_t>::max();

  // verbosity flag
  _verbose = p.get<bool>("verbose");

  // producer name/labels
  _LG_producer  = p.get<std::string>( "LGProducer"  );
  _LG_label     = p.get<std::string>( "LGLabel","*" );
  _HG_producer = p.get<std::string>( "HGProducer"  );
  _HG_label    = p.get<std::string>( "HGLabel","*" );
  _LG_cosmic_producer = p.get<std::string>( "LGProducerCosmic"  );
  _LG_cosmic_label    = p.get<std::string>( "LGLabelCosmic","*" );
  _HG_cosmic_producer = p.get<std::string>( "HGProducerCosmic"  );
  _HG_cosmic_label    = p.get<std::string>( "HGLabelCosmic","*" );

  _save_wfs = p.get<bool>("SaveWaveformInTTree");
  
  // make tree
  MakeTree();

  _use_LG_beam_for_HG_cosmic = p.get<bool>("UseLGBeamForHGCosmic");
  if(_use_LG_beam_for_HG_cosmic && (!_LG_cosmic_producer.empty() || !_LG_cosmic_label.empty())) {
    std::cerr << "UseLGBeamForHGCosmic enabled but LGProducerCosmic and/or LGLabelCosmic also provided!" << std::endl;
    throw std::exception();
  }

  auto const ignore_channels = p.get<std::vector<size_t> >("IgnoreChannel");
  for(auto const& ch : ignore_channels) _ignore_channels.insert(ch);
}

void OpDigitSaturationCorrection::SetUpChannelMap()
{

  if(!_opch_to_opdet_m.empty()) return;

  // Get generic geometry
  ::art::ServiceHandle<geo::Geometry> geom;
  auto const& channelMap = art::ServiceHandle<geo::WireReadout const>()->Get();

  // Get Optical Readout Channel Map Geometry
  ::art::ServiceHandle<geo::UBOpReadoutMap> ub_geom;

  // Load a full set of Optical Channels
  auto const& readout_channel_set = ub_geom->GetReadoutChannelSet();

  // Loop over channel set, create OpChannel=>OpDet and OpChannel=>ChannelType mapping
  for(auto const& ch : readout_channel_set) {
    
    // Make sure mapping has enough entries
    if(ch >= _opch_to_opdet_m.size()) {
      _opch_to_chtype_m.resize(ch+1,opdet::Undefined);
      _opch_to_chcategory_m.resize(ch+1,opdet::Uncategorized);
      _opch_to_opdet_m.resize(ch+1,-1);
      _opch_to_baseline_m.resize(ch+1,_baseline);
      _opch_to_gain_m.resize(ch+1,1.);
    }

    bool skip=true;
    for(size_t i=0; i<geom->Cryostat().NOpDet(); ++i) {
      try{
        skip = !(channelMap.OpDetFromOpChannel(ch) < geom->Cryostat().NOpDet());
      }catch(...){
	skip = true;
      }
      if(!skip) break;
    }

    _opch_to_chcategory_m[ch] = ub_geom->GetChannelCategory(ch);

    _opch_to_chtype_m[ch] = ub_geom->GetChannelType(ch);

    if(skip) _opch_to_opdet_m[ch] = -1;
    else _opch_to_opdet_m[ch] = channelMap.OpDetFromOpChannel(ch);

    // if this channel is found in the user-defined baseline map, add baseline value
    for (size_t i=0; i < _baseline_chan_v.size(); i++){
      auto const& chan = _baseline_chan_v[i];
      if (chan == ch){
	_opch_to_baseline_m[ ch ] = _baseline_vals_v[ i ];
	break;
      }
    }

    // if this channel is found in the user-defined gain map, add baseline value
    for (size_t i=0; i < _gain_chan_v.size(); i++){
      auto const& chan = _gain_chan_v[i];
      if (chan == ch){
	_opch_to_gain_m[ ch ] = _gain_vals_v[ i ];
	break;
      }
    }
    
  }

  // Reverse-engineer OpDet => set of OpChannel mapping (probably not really useful)
  for(size_t opch=0; opch < _opch_to_opdet_m.size(); ++opch) {
    auto const& opdet = _opch_to_opdet_m[opch];
    if(opdet < 0) continue;
    if(_opdet_to_opch_m.size() <= (size_t)opdet) _opdet_to_opch_m.resize(opdet+1,std::set<int>());
    _opdet_to_opch_m[opdet].insert(opch);
  }

  // loop through _opdet_to_opch_m and for the various entries find which HG and LG channel numbers
  // are associated with the same FEM
  for(auto const& opch_set : _opdet_to_opch_m) {
    
    // map beamgate LG to HG
    for (auto const& opch1 : opch_set) {

      if ( _opch_to_chcategory_m[opch1] == opdet::OpdetBeamHighGain){
	//std::cout<<"mapping high gain channel " << opch1 << " (BeamHighGain) ... " << std::flush;
	for(auto const& opch2 : opch_set) {
	  if ( _opch_to_chcategory_m[opch2] == opdet::OpdetBeamLowGain) {
	    _opchLG_to_opchHG_m[ opch2 ] = opch1;
	    //std::cout<<opch2 << " (BeamLowGain) " << std::endl;
	    break;
	  }
	}
      }// map Beam LG => HG
      if ( _opch_to_chcategory_m[opch1] == opdet::OpdetCosmicHighGain){
	//std::cout<<"mapping high gain channel " << opch1 << " (CosmicHighGain) ... " << std::flush;
	for(auto const& opch2 : opch_set) {
	  if ( !_use_LG_beam_for_HG_cosmic && _opch_to_chcategory_m[opch2] == opdet::OpdetCosmicLowGain) {
	    _opchLG_to_opchHG_m[ opch2 ] = opch1;
	    //std::cout<<opch2 << " (CosmicLowGain) " << std::endl;
	    break;
	  }
	  if ( _use_LG_beam_for_HG_cosmic && _opch_to_chcategory_m[opch2] == opdet::OpdetBeamLowGain) {
	    _opchLG_to_opchHG_m[ opch2 ] = opch1;
	    //std::cout<<opch2 << " (BeamLowGain) " << std::endl;
	    break;
	  } 
	}
      }// map Cosmic LG => HG
    }
    
  }// for all PMTs

  // Report
  if(_verbose) {
    std::cout << "Loaded OpChannel => OpDet ... Channel Type ... Channel Category mapping" << std::endl;
    for(size_t opch=0; opch<_opch_to_opdet_m.size(); ++opch) {
      //if(_opch_to_opdet_m[opch] < 0) continue;
      std::cout << opch << " => " << _opch_to_opdet_m[opch];
      switch(_opch_to_chtype_m[opch]) {
      case opdet::Undefined:    std::cout << " UNDEFINED!"; break;
      case opdet::HighGain:     std::cout << " HG";         break;
      case opdet::LowGain:      std::cout << " LG";         break;
      case opdet::LogicChannel: std::cout << " LOGIC";      break;
      default:
	std::cerr << "Cannot decode the type: " << _opch_to_chtype_m[opch] << std::endl;
	throw std::exception();
      }
      std::cout << " ...";
      switch(_opch_to_chcategory_m[opch]) {
      case opdet::Uncategorized:       std::cout << " UNDEFINED!"        << std::endl; break;
      case opdet::UnspecifiedLogic:    std::cout << " Unspecified Logic" << std::endl; break;
      case opdet::OpdetCosmicHighGain: std::cout << " Cosmic High Gain " << std::endl; break;
      case opdet::OpdetCosmicLowGain:  std::cout << " Cosmic Low Gain " << std::endl; break;
      case opdet::OpdetBeamHighGain:   std::cout << " Beam High Gain " << std::endl; break;
      case opdet::OpdetBeamLowGain:    std::cout << " Beam Low Gain " << std::endl; break;
      case opdet::BNBLogicPulse:       std::cout << " BNB Logic Pulse " << std::endl; break;
      case opdet::NUMILogicPulse:      std::cout << " NuMI Logic Pulse " << std::endl; break;
      case opdet::FlasherLogicPulse:   std::cout << " Flasher Logic Pulse " << std::endl; break;
      case opdet::StrobeLogicPulse:    std::cout << " Strobe Logic Pulse " << std::endl; break;
      default:
	std::cerr << "Cannot decode the category: " << _opch_to_chcategory_m[opch] << std::endl;
	throw std::exception();
      }
    }
    
    std::cout << std::endl;
    std::cout << "Loaded LG OpChannel => HG OpChannel mapping" << std::endl;
    for(auto const& pair : _opchLG_to_opchHG_m)
      std::cout << pair.first << " => " << pair.second << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Loaded OpDet => OpChannel mapping" << std::endl;
    for(size_t opdet=0; opdet<_opdet_to_opch_m.size(); ++opdet) {
      std::cout << opdet << " => ";
      for(auto const& opch : _opdet_to_opch_m[opdet]) { std::cout << opch << " "; }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Loaded OpChan => Baseline mapping" << std::endl;
    for(size_t opch=0; opch < _opch_to_baseline_m.size(); ++opch) 
      std::cout << opch << " => " << _opch_to_baseline_m[ opch ] << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Loaded OpChan => Gain mapping" << std::endl;
    for(size_t opch=0; opch < _opch_to_gain_m.size(); ++opch) 
      std::cout << opch << " => " << _opch_to_gain_m[ opch ] << std::endl;
    std::cout << std::endl;

  }
}

void OpDigitSaturationCorrection::produce(art::Event & e)
{

  if (_verbose)
    std::cout << "Begin new event. Set up channel map." << std::endl;
  
  SetUpChannelMap();

  if (_verbose)
    std::cout << "Start processing data." << std::endl;

  // update tree variables for this event
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // produce OpDetWaveform data-product to be filled within module
  //std::unique_ptr< std::vector<raw::OpDetWaveform> > corrected_wfs(new std::vector<raw::OpDetWaveform>);
  std::unique_ptr< std::vector<raw::OpDetWaveform> > corrected_beam_wfs(new std::vector<raw::OpDetWaveform>);
  std::unique_ptr< std::vector<raw::OpDetWaveform> > corrected_cosmic_wfs(new std::vector<raw::OpDetWaveform>);
  
  // the output waveforms will carry the channel number of the corresponding HG channel they are storing info for
  // if cosmic-discriminated waveforms for OpDet PMT 01 are being read out by OpChan channel numbers 01 for HG and 101 for LG
  // output for this OpDet's cosmic-discriminated data will go to OpChan 01.
  // similary for the beam-gate readout for the same OpDet PMT. If it is read-out by OpChan 201 (HG) and 301 (LG) the data for this
  // stream will be stored in channel number 201.

  // load OpDetWaveform from High Gain raw data: Cosmic Disc. FEM
  art::Handle<std::vector<raw::OpDetWaveform> > opwf_HG_cosmic_v;
  e.getByLabel(_HG_cosmic_producer, _HG_cosmic_label, opwf_HG_cosmic_v);
  // make sure hits look good
  if(!opwf_HG_cosmic_v.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate HG Cosmic Disc. OpDetWf!"<<std::endl;
    throw std::exception();
  }  

  // load OpDetWaveform from High Gain raw data: Beam Window FEM
  art::Handle<std::vector<raw::OpDetWaveform> > opwf_HG_v;
  e.getByLabel(_HG_producer, _HG_label, opwf_HG_v);
  // make sure hits look good
  if(!opwf_HG_v.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate HG OpDetWf!"<<std::endl;
    throw std::exception();
  }

  // load OpDetWaveform from Low Gain raw data: Cosmic Disc. FEM 
  art::Handle<std::vector<raw::OpDetWaveform> > opwf_LG_cosmic_v;
  e.getByLabel( (_use_LG_beam_for_HG_cosmic ? _LG_producer : _LG_cosmic_producer),
		(_use_LG_beam_for_HG_cosmic ? _LG_label    : _LG_cosmic_label   ),
		opwf_LG_cosmic_v);
  // make sure hits look good
  if(!opwf_LG_cosmic_v.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate LG Cosmic Disc. OpDetWf!"<<std::endl;
    throw std::exception();
  }
  
  // load OpDetWaveform from Low Gain raw data: Beam Window FEM
  art::Handle<std::vector<raw::OpDetWaveform> > opwf_LG_v;
  e.getByLabel(_LG_producer, _LG_label, opwf_LG_v);
  // make sure hits look good
  if(!opwf_LG_v.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate LG OpDetWf!"<<std::endl;
    throw std::exception();
  }

  // load trigger data
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  // get the trigger time
  _TrigTime = clockData.TriggerTime();

  if (_verbose)
    std::cout << "Check that channel numbers and mappings are as expected" << std::endl;

  // validate that the channel numbers associated with HG / LG pulses agree with what expected
  // from the channel mapping.
  for (size_t idx = 0; idx < opwf_LG_v->size(); idx++){
    auto const& wf   = opwf_LG_v->at(idx);
    auto const& opch = wf.ChannelNumber();
    auto const& cat  = _opch_to_chcategory_m[opch];
    switch(cat) {
    case opdet::UnspecifiedLogic:    
    case opdet::BNBLogicPulse:       
    case opdet::NUMILogicPulse:      
    case opdet::FlasherLogicPulse:   
    case opdet::StrobeLogicPulse:    
    case opdet::OpdetBeamLowGain: break;
    default:
      std::cerr<<"\033[93m[ERROR]\033[00m ... channel " << opch << " from LG stream associated with producer "
	       << _LG_producer << " is not marked as Beam LG in map produced by SetUpChannelMap function." <<std::endl;
      throw std::exception();
    }
  }
  for (size_t idx = 0; idx < opwf_HG_v->size(); idx++){
    auto const& wf   = opwf_HG_v->at(idx);
    auto const& opch = wf.ChannelNumber();
    auto const& cat  = _opch_to_chcategory_m[opch];
    switch(cat) {
    case opdet::UnspecifiedLogic:    
    case opdet::BNBLogicPulse:       
    case opdet::NUMILogicPulse:      
    case opdet::FlasherLogicPulse:   
    case opdet::StrobeLogicPulse:    
    case opdet::OpdetBeamHighGain: break;
    default:
      std::cerr<<"\033[93m[ERROR]\033[00m ... channel " << opch << " from HG stream associated with producer "
	       << _LG_producer << " is not marked as Beam HG in map produced by SetUpChannelMap function." <<std::endl;
      throw std::exception();
    }
  }
  for (size_t idx = 0; idx < opwf_LG_cosmic_v->size(); idx++){
    auto const& wf   = opwf_LG_cosmic_v->at(idx);
    auto const& opch = wf.ChannelNumber();
    auto const& cat  = _opch_to_chcategory_m[opch];
    bool skip=true;
    switch(cat) {
    case opdet::UnspecifiedLogic:    
    case opdet::BNBLogicPulse:       
    case opdet::NUMILogicPulse:      
    case opdet::FlasherLogicPulse:   
    case opdet::StrobeLogicPulse:    
      skip=true;
    default:
      break;
    }
    if(skip) continue;

    static bool onetime_warning = false;
    if (_use_LG_beam_for_HG_cosmic) {
      if(_opch_to_chcategory_m[opch] != opdet::OpdetBeamLowGain){
	std::cerr<<"\033[93m[ERROR]\033[00m ... channel " << opch << " from LG stream associated with producer "
		 << _LG_producer << " is not marked as Beam LG in map produced by SetUpChannelMap function." <<std::endl;
	throw std::exception();
      }
      if(!onetime_warning) {
	std::cout<<"\033[95m[WARNING]\033[00m ... Using Beam LG to correct Cosmic HG!" << std::endl;
	onetime_warning = true;
      }
      continue;
    }
    if (_opch_to_chcategory_m[opch] != opdet::OpdetCosmicLowGain){
      std::cerr<<"\033[93m[ERROR]\033[00m ... channel " << opch << " from LG stream associated with producer "
	       << _LG_cosmic_producer << " is not marked as Cosmic LG in map produced by SetUpChannelMap function." <<std::endl;
      throw std::exception();
    }
  }
  for (size_t idx = 0; idx < opwf_HG_cosmic_v->size(); idx++){
    auto const& wf   = opwf_HG_cosmic_v->at(idx);
    auto const& opch = wf.ChannelNumber();
    auto const& cat  = _opch_to_chcategory_m[opch];
    switch(cat) {
    case opdet::UnspecifiedLogic:    
    case opdet::BNBLogicPulse:       
    case opdet::NUMILogicPulse:      
    case opdet::FlasherLogicPulse:   
    case opdet::StrobeLogicPulse:    
    case opdet::OpdetCosmicHighGain: break;
    default:
      std::cerr<<"\033[93m[ERROR]\033[00m ... channel " << opch << " from HG stream associated with producer "
	       << _LG_producer << " is not marked as Cosmic HG in map produced by SetUpChannelMap function." <<std::endl;
      throw std::exception();
    }
  }

  if (_verbose)
    std::cout << "Fill LG_BEAM_ChanMap : there are " << opwf_LG_v->size() << " waveforms." << std::endl;

  // fill a vector (for LG BEAM):
  // [ pmt channel ] -> [ idx(wf1), idx(wf2), ...]
  // to quickly find the position of wach LG wf
  //std::vector<std::vector< std::pair<size_t,double> > > LG_BEAM_ChanMap;
  LGChanMap_t LG_BEAM_ChanMap;
  LG_BEAM_ChanMap.clear();
  for (size_t pmt=0; pmt < 32; pmt++)
    LG_BEAM_ChanMap.push_back( std::vector<LGPulseInfo_t>() );
  for (size_t idx = 0; idx < opwf_LG_v->size(); idx++){
    auto const& wf = opwf_LG_v->at(idx);
    // retrieve PMT OpDet number
    auto const& opdet = _opch_to_opdet_m[ wf.ChannelNumber() ];
    // skip if irrelevant channel
    if(_opch_to_chcategory_m[wf.ChannelNumber()] != ::opdet::OpdetBeamLowGain) continue;
    if(_ignore_channels.find(wf.ChannelNumber()) != _ignore_channels.end()) continue;
    // LG indx = ChannelNumber() - 100 (to match HG channel number)
    double wf_time = wf.TimeStamp() - _TrigTime;
    LGPulseInfo_t p;
    p._index = idx;
    p._start = wf_time;
    p._duration = wf.size() * _TDC;
    if(opdet<0) {
      std::cerr << "\033[93m[ERROR]\033[00m ... channel " << wf.ChannelNumber() << " has no corresponding opdet!" << std::endl;
    }
    LG_BEAM_ChanMap.at( opdet ).emplace_back( std::move(p) );
  }

  if (_verbose)
    std::cout << "Fill LG_COSMIC_ChanMap : there are " << opwf_LG_cosmic_v->size() << " waveforms." << std::endl;

  // fill a vector (for LG COSMIC):
  // [ pmt channel ] -> [ idx(wf1), idx(wf2), ...]
  // to quickly find the position of wach LG wf
  LGChanMap_t LG_COSMIC_ChanMap;
  LG_COSMIC_ChanMap.clear();
  for (size_t pmt=0; pmt < 32; pmt++)
    LG_COSMIC_ChanMap.push_back( std::vector<LGPulseInfo_t>() );
  for (size_t idx = 0; idx < opwf_LG_cosmic_v->size(); idx++){
    auto const& wf = opwf_LG_cosmic_v->at(idx);
    // retrieve PMT OpDet number
    auto const& opdet = _opch_to_opdet_m[ wf.ChannelNumber() ];
    // If this "cosmic" is actually "beam", ignore beam window waveforms
    if(_use_LG_beam_for_HG_cosmic && wf.size() > _BEAMGATE_MIN) continue;
    // skip if irrelevant wf channel
    if(_opch_to_chcategory_m[wf.ChannelNumber()] != ::opdet::OpdetCosmicLowGain) continue;
    if(_ignore_channels.find(wf.ChannelNumber()) != _ignore_channels.end()) continue;
    // LG indx = ChannelNumber() - 100 (to match HG channel number)
    double wf_time = wf.TimeStamp() - _TrigTime;
    LGPulseInfo_t p;
    p._index = idx;
    p._start = wf_time;
    p._duration = wf.size() * _TDC;
    if(opdet<0) {
      std::cerr << "\033[93m[ERROR]\033[00m ... channel " << wf.ChannelNumber() << " has no corresponding opdet!" << std::endl;
    }
    LG_COSMIC_ChanMap.at( opdet ).emplace_back( std::move(p) );
  }

  // now loop through HG waveforms
  // if any one saturates
  // find the corresponding LG waveform
  // if available
  // append to the list of new WFs
  // either the un-saturated HG ones
  // or the corrected LG ones

  if (_verbose)
    std::cout << "loop through HG WFs" << std::endl;

  // loop through HG BEAM waveforms
  for (size_t idx = 0; idx < opwf_HG_v->size(); idx++){
    auto const& wf_HG = opwf_HG_v->at(idx);

    if (_verbose)
      std::cout << "Starting a new HG BEAM wf w/ chan num " << wf_HG.ChannelNumber() 
		<< " and " << wf_HG.size() << " ticks" << std::endl;

    // ignore registered channels to be ignored
    if(_ignore_channels.find(wf_HG.ChannelNumber()) != _ignore_channels.end()) continue;
    
    // retrieve OpDet channel number
    auto const& opdet = _opch_to_opdet_m[ wf_HG.ChannelNumber() ];
    if(opdet<0) {
      std::cerr << "\033[93m[ERROR]\033[00m ... channel " << wf_HG.ChannelNumber() << " has no corresponding opdet!" << std::endl;
    }

    if (_verbose)
      std::cout << "corresponding to PMT OpDet number " << opdet << std::endl;
    
    // ignore PMTs with size != cosmic discriminator size
    if (!_include_beamgate && wf_HG.size() > _BEAMGATE_MIN)
      continue;
    
    _time_LG = -1;
    _max_adc_LG = _max_tdc_LG = -1;
    _max_adc_corr = -1;
    
    auto const& max_idx = getMaxADC(wf_HG);
    double wf_time = wf_HG.TimeStamp() - _TrigTime;

    _chan_HG    = wf_HG.ChannelNumber();    
    _time_HG    = wf_time;
    _max_adc_HG = wf_HG.at(max_idx) - _opch_to_baseline_m[ _chan_HG ];
    _max_tdc_HG = max_idx;
    _pmt        = opdet;
    _matched    = 0;
    _swapped    = 0;
    _wf_HG      = wf_HG;
    _nticks_HG  = wf_HG.size();

    // skip if irrelevant channel
    if(_opch_to_chcategory_m[_chan_HG] != ::opdet::OpdetBeamHighGain) continue;

    if (_verbose)
      std::cout << "This HG channel has Ch Num : " << wf_HG.ChannelNumber()
		<< " w/ TimeStamp : " << wf_time
		<< " w/ size : " << wf_HG.size() << std::endl
		<< "\t max ADC @ " << max_idx << " is " << wf_HG.at(max_idx) << std::endl;
    
    // find matching LG channel
    if (_verbose)
      std::cout << "Try and find a matching LG pulse" << std::endl;
    auto const& LG_idx = FindMatchingLGPulse(_pmt, wf_time, LG_BEAM_ChanMap);
    
    // if we did not find a match:
    // add the old HG waveform
    if (LG_idx == IDX_max){
      if (_verbose)
	std::cout << "no match found -> add saturated HG wf" << std::endl;
      corrected_beam_wfs->push_back(wf_HG);
    }
    
    else{
      // if we made it this far -> we found a match in LG information!
      if (_verbose)
	std::cout << "getting corresponding LG pulse w/ idx " << LG_idx << std::endl;
      auto const& wf_LG = opwf_LG_v->at(LG_idx);
      auto const& max_idx_LG = getMaxADC(wf_LG);
      _matched      = 1;
      _chan_LG      = wf_LG.ChannelNumber();
      _time_LG      = wf_LG.TimeStamp() - _TrigTime;
      _max_adc_LG   = wf_LG.at(max_idx_LG) - _opch_to_baseline_m[ _chan_LG ];
      _max_tdc_LG   = max_idx_LG;
      _max_adc_corr = _max_adc_LG * ( _gain_fact * _opch_to_gain_m[ _chan_HG ] );
      _wf_LG        = wf_LG;
      _nticks_LG    = wf_LG.size();
      if (_verbose)
	std::cout << "finished getting corresponding LG pulse w/ idx " << LG_idx << std::endl;
      
      // if the max ADC value for the HG was below saturation -> add the HG wf
      if (_max_adc_HG + 10 < (4095 - _opch_to_baseline_m[ _chan_HG ]) ){
	if (_verbose)
	  std::cout << "HG pulse not saturated -> add HG pulse" << std::endl;
	corrected_beam_wfs->push_back(wf_HG);
      }
      else{
	  if (_verbose)
	    std::cout << "found a match and HG saturates! -> add LG wf w/ correction" << std::endl;
	  // create a new waveform by editing the LG one
	  if (_verbose)
	    std::cout << "corr factor for chan " << wf_LG.ChannelNumber()
		      << " : " << _opch_to_gain_m[ _chan_HG ] << std::endl;
	  std::vector<short unsigned int> adcs;
	  for (size_t n=0; n < wf_LG.size(); n++){
	    int this_ADC_above_baseline = wf_LG.at(n) - _opch_to_baseline_m [ _chan_LG ];
	    // make sure we don't overflow the data-product (not the firmware waveform...)
	    if (this_ADC_above_baseline > (int)( ADC_max / (_gain_fact * _opch_to_gain_m[ _chan_HG ] ) ) ){
	      if (_verbose) std::cout << "new ADC (saturated) : " << ADC_max << std::endl;
	      adcs.push_back(ADC_max);
	    }
	    else{
	      short unsigned int new_ADC = (short unsigned int)( this_ADC_above_baseline * (_gain_fact * _opch_to_gain_m[ _chan_HG ] ) + _opch_to_baseline_m[ _chan_LG ] );
	      if (_verbose) std::cout << "new ADC (not saturated) : " << new_ADC << std::endl;
	      adcs.push_back( new_ADC );
	    }
	  }// for all ADCs
	  raw::OpDetWaveform new_wf(wf_LG.TimeStamp(),
				    wf_HG.ChannelNumber(),
				    //_opchLG_to_opchHG_m [ wf_LG.ChannelNumber() ],
				    adcs);
	  corrected_beam_wfs->push_back(new_wf);
	  _wf_corr = new_wf;
	  _swapped = 1;
	  if (_verbose)
	    std::cout << "saved a corrected waveform w/ " << _wf_corr.size() 
		      << " ticks @ channel number " << new_wf.ChannelNumber() << std::endl;
      }// if this waveform saturates
    }// if LG wf was found
    
    if (_verbose)
      std::cout << "fill tree!" << std::endl;
    if (_tree) _tree->Fill();
    
    if (_verbose)
      std::cout << "move to next wf..." << std::endl;
    
  }// for all WFs

  // loop through HG COSMICS waveforms
  for (size_t idx = 0; idx < opwf_HG_cosmic_v->size(); idx++){
    auto const& wf_HG = opwf_HG_cosmic_v->at(idx);

    if (_verbose)
      std::cout << "Starting a new HG COSMIC wf w/ chan num " << wf_HG.ChannelNumber()
		<< " and " << wf_HG.size() << " ticks" << std::endl;


    // ignore registered channels to be ignored
    if(_ignore_channels.find(wf_HG.ChannelNumber()) != _ignore_channels.end()) continue;
    
    // retrieve OpDet channel number
    auto const& opdet = _opch_to_opdet_m[ wf_HG.ChannelNumber() ];
    if(opdet<0) {
      std::cerr << "\033[93m[ERROR]\033[00m ... channel " << wf_HG.ChannelNumber() << " has no corresponding opdet!" << std::endl;
    }

    if (_verbose)
      std::cout << "corresponding to PMT OpDet number " << opdet << std::endl;
    
    // ignore PMTs with size != cosmic discriminator size
    if (!_include_beamgate && wf_HG.size() > _BEAMGATE_MIN)
      continue;
    
    _time_LG = -1;
    _max_adc_LG = _max_tdc_LG = -1;
    _max_adc_corr = -1;
    
    auto const& max_idx = getMaxADC(wf_HG);
    double wf_time = wf_HG.TimeStamp() - _TrigTime;

    _chan_HG    = wf_HG.ChannelNumber();    
    _time_HG    = wf_time;
    _max_adc_HG = wf_HG.at(max_idx) - _opch_to_baseline_m[ _chan_HG ];
    _max_tdc_HG = max_idx;
    _pmt        = opdet;
    _matched    = 0;
    _swapped    = 0;
    _wf_HG      = wf_HG;
    _nticks_HG  = wf_HG.size();

    if(_opch_to_chcategory_m[_chan_HG] != ::opdet::OpdetCosmicHighGain) continue;

    if (_verbose)
      std::cout << "This HG channel has Ch Num : " << wf_HG.ChannelNumber()
		<< " w/ TimeStamp : " << wf_time
		<< " w/ size : " << wf_HG.size() << std::endl
		<< "\t max ADC @ " << max_idx << " is " << wf_HG.at(max_idx) << std::endl;
    
    // find matching LG channel
    if (_verbose)
      std::cout << "Try and find a matching LG pulse" << std::endl;
    auto const& LG_idx = FindMatchingLGPulse(_pmt, wf_time, LG_BEAM_ChanMap);
    
    // if we did not find a match:
    // add the old HG waveform
    if (LG_idx == IDX_max){
      if (_verbose)
	std::cout << "no match found -> add saturated HG wf" << std::endl;
      corrected_cosmic_wfs->push_back(wf_HG);
    }
    
    else{
      // if we made it this far -> we found a match in LG information!
      if (_verbose)
	std::cout << "getting corresponding LG pulse w/ idx " << LG_idx << std::endl;
      auto const& wf_LG = opwf_LG_cosmic_v->at(LG_idx);
      auto const& max_idx_LG = getMaxADC(wf_LG);
      _matched      = 1;
      _chan_LG      = wf_LG.ChannelNumber();
      _time_LG      = wf_LG.TimeStamp() - _TrigTime;
      _max_adc_LG   = wf_LG.at(max_idx_LG) - _opch_to_baseline_m[ _chan_LG ];
      _max_tdc_LG   = max_idx_LG;
      _max_adc_corr = _max_adc_LG * ( _gain_fact * _opch_to_gain_m[ _chan_HG ] );
      _wf_LG        = wf_LG;
      _nticks_LG    = wf_LG.size();
      if (_verbose)
	std::cout << "finished getting corresponding LG pulse w/ idx " << LG_idx << std::endl;
      
      // if the max ADC value for the HG was below saturation -> add the HG wf
      if (_max_adc_HG + 10 < (4095 - _opch_to_baseline_m[ _chan_HG ]) ){
	if (_verbose)
	  std::cout << "HG pulse not saturated -> add HG pulse" << std::endl;
	corrected_cosmic_wfs->push_back(wf_HG);
      }
      else{
	  if (_verbose)
	    std::cout << "found a match and HG saturates! -> add LG wf w/ correction" << std::endl;
	  // create a new waveform by editing the LG one
	  if (_verbose)
	    std::cout << "corr factor for chan " << wf_LG.ChannelNumber()
		      << " : " << _opch_to_gain_m[ _chan_HG ] << std::endl;
	  std::vector<short unsigned int> adcs;
	  for (size_t n=0; n < wf_LG.size(); n++){
	    int this_ADC_above_baseline = wf_LG.at(n) - _opch_to_baseline_m[ _chan_LG ];
	    // make sure we don't overflow the data-product (not the firmware waveform...)
	    if (this_ADC_above_baseline > (int)( ADC_max / (_gain_fact * _opch_to_gain_m[ _chan_HG ] ) ) ){
	      if (_verbose) std::cout << "new ADC (saturated) : " << ADC_max << std::endl;
	      adcs.push_back(ADC_max);
	    }
	    else{
	      short unsigned int new_ADC = (short unsigned int)( this_ADC_above_baseline * (_gain_fact * _opch_to_gain_m[ _chan_HG ] ) + _opch_to_baseline_m[ _chan_LG ] );
	      if (_verbose) std::cout << "new ADC (not saturated) : " << new_ADC << std::endl;
	      adcs.push_back( new_ADC );
	    }
	  }// for all ADCs
	  raw::OpDetWaveform new_wf(wf_LG.TimeStamp(),
				    wf_HG.ChannelNumber(),
				    //_opchLG_to_opchHG_m [ _chan_LG ],
				    adcs);
	  _wf_corr = new_wf;
	  _swapped = 1;
	  if (_verbose)
	    std::cout << "saved a corrected waveform w/ " << _wf_corr.size() << " ticks" << std::endl;
	  corrected_cosmic_wfs->push_back(new_wf);
      }// if this waveform saturates
    }// if LG wf was found
    
    if (_verbose)
      std::cout << "fill tree!" << std::endl;
    if (_tree) _tree->Fill();
    
    if (_verbose)
      std::cout << "move to next wf..." << std::endl;
    
  }// for all COSMICS WFs

  if (_verbose)
    std::cout << "done merging HG/LG. Save waveforms" << std::endl;
  
  e.put(std::move(corrected_beam_wfs),opdet::UBOpChannelEnumName( opdet::OpdetBeamHighGain ));
  e.put(std::move(corrected_cosmic_wfs),opdet::UBOpChannelEnumName( opdet::OpdetCosmicHighGain ));

  if (_verbose)
    std::cout << "done with this event." << std::endl;
}

void OpDigitSaturationCorrection::MakeTree()
{
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("saturation_tree","Saturation Correction Tree");
  
  // Simple type branch
  _tree->Branch ( "run",        &_run,        "run/i"       );
  _tree->Branch ( "subrun",     &_subrun,     "subrun/i"    );
  _tree->Branch ( "event",      &_event,      "event/i"     );
  _tree->Branch ( "chan_HG",    &_chan_HG,    "chan_HG/i"   );
  _tree->Branch ( "chan_LG",    &_chan_LG,    "chan_LG/i"   );
  _tree->Branch ( "pmt",        &_pmt,        "pmt/i"       );
  _tree->Branch ( "matched",    &_matched,    "matched/I"   ); // both HG and LG wfs present
  _tree->Branch ( "swapped",    &_swapped,    "swapped/I"   ); // LG corrected swapped in for HG 
  _tree->Branch ( "_nticks_LG",   &_nticks_LG,    "nticks_LG/I"   ); // number of ticks in waveform
  _tree->Branch ( "_nticks_HG",   &_nticks_HG,    "nticks_HG/I"   ); // number of ticks in waveform
  _tree->Branch ( "time_HG",    &_time_HG,    "time_HG/D"   );
  _tree->Branch ( "max_adc_HG", &_max_adc_HG, "max_adc_HG/I");
  _tree->Branch ( "max_tdc_HG", &_max_tdc_HG, "max_tdc_HG/I");
  _tree->Branch ( "time_LG",    &_time_LG,    "time_LG/D"   );
  _tree->Branch ( "max_adc_LG", &_max_adc_LG, "max_adc_LG/I");
  _tree->Branch ( "max_tdc_LG", &_max_tdc_LG, "max_tdc_LG/I");
  _tree->Branch ( "max_adc_corr", &_max_adc_corr, "max_adc_corr/I");
  if (_save_wfs == true){
    _tree->Branch ( "_wf_LG","std::vector<short>", &_wf_LG    );
    _tree->Branch ( "_wf_HG","std::vector<short>", &_wf_HG    );
    _tree->Branch ( "_wf_corr","std::vector<short>", &_wf_corr    );
  }
}

DEFINE_ART_MODULE(OpDigitSaturationCorrection)

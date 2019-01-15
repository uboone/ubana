/**
 * \file OpDetWaveformAna
 *
 * \ingroup OpticalDetectorAna
 * 
 * \brief Class def header for a class OpDetWaveformAna
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorAna

    @{*/
#ifndef OPTICALDETECTORANA_OPDETWAVEFORMANA_H
#define OPTICALDETECTORANA_OPDETWAVEFORMANA_H

#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
// #include "larana/OpticalDetector/OpHitFinder/AlgoPedestal.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include <TTree.h>
#include <string>
/**
   \class OpDetWaveformAna
 */
namespace pmtana {
  class OpDetWaveformAna {
    
  public:
    
    /// Default constructor
    OpDetWaveformAna(const std::string name="noname");
    
    /// Default destructor
    virtual ~OpDetWaveformAna() {}

    void ClearEvent    ();
    void ClearWaveform ();
    
    void AnaHit            ( TTree* ptr );
    void AnaWaveform       ( TTree* ptr );
    void SaveWaveform      ( TTree* ptr );
    void SaveEvWaveform    ( TTree* ptr );

    void TickPeriod   ( const double period );    
    void SetEventInfo ( const unsigned int run,
		        const unsigned int subrun,
		        const unsigned int event  );
    void AnaWaveform  ( const unsigned int ch,
			const double time_wrt_trigger,
			const std::vector<short>& wf);
    void AnaEventWaveform  ( const std::vector<raw::OpDetWaveform>& ev_wf_v);

    ::pmtana::PulseRecoManager& GetManager() { return _preco_mgr; };
    ::pmtana::AlgoThreshold&    GetAlgo()    { return _preco_alg; }

  private:

    ::pmtana::PulseRecoManager  _preco_mgr;
    ::pmtana::AlgoThreshold     _preco_alg;

    std::string _name;
    TTree* _hitana_tree;
    TTree* _wfana_tree;
    TTree* _wf_tree;
    TTree* _ev_wf_tree;

    double _period;
    unsigned int _run;
    unsigned int _subrun;
    unsigned int _event;
    unsigned int _ch;
    float  _ped_mean;
    float  _ped_rms;
    std::vector<float>  _ped_mean_v;
    std::vector<float>  _ped_rms_v;
    double _tstart;
    double _tpeak;
    double _tend;
    double _t_wrt_trigger;
    double _q;
    double _amp;
    unsigned short _max_adc;
    unsigned short _min_adc;
    std::vector<unsigned short> _max_adc_v;
    std::vector<unsigned short> _min_adc_v;
    unsigned int   _wf_size;
    std::vector<short> _wf;
    std::vector<int> _wfsum;
    std::vector<std::vector<short>> _wf_v;
  };
}

#endif
/** @} */ // end of doxygen group 

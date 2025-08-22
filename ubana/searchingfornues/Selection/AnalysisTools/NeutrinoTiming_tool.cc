#ifndef ANALYSIS_NEUTRINOTIMING_CXX
#define ANALYSIS_NEUTRINOTIMING_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

// backtracking tools
#include "ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"
#include "ubana/searchingfornues/Selection/CommonDefs/TrackShowerScoreFuncs.h"

//Flux info
#include "dk2nu/tree/dk2nu.h"
#include "ubobj/ReBooNE/BooNEInfo.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       NeutrinoTiming
    // File:        NeutrinoTiming.cc
    //
    //              A basic analysis example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by David Caratelli (dcaratelli@ucsb.edu) on 10/06/2022
    // Code developed by Dante Totani and exported to LArSoft by Erin Yandel
    // Modified by Anyssa Navrer-Agasson March 2025 to include the ns timing 
    // for BNB (MC) and NuMI (data & MC).
    ////////////////////////////////////////////////////////////////////////

  class NeutrinoTiming : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    NeutrinoTiming(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~NeutrinoTiming(){ };
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
    void analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree) override;

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree) override;

    void AddDaughters(const art::Ptr<recob::PFParticle>& pfp,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);
    
  private:

    art::InputTag fPFPproducer;
    art::InputTag fT0producer;
    art::InputTag fFLASHproducer;
    art::InputTag fPMTWFproducer;
    art::InputTag fLogicWFproducer;
    art::InputTag fSpacePointproducer;

    //beam timing plots
    TH1F *H_time;
    TH1F *H_maxH;
    TH1F *H_t0_Beam;
    TH1F *H_SimTime;
    TH1F *H_TruthTime;
    TH1F *H_ns_time;
    TH2F *H_TimeVsPh;

    // beam timing variables
    //int _tickSum;
    //int _tickP[32];
    //double _maxSum; 
    //double _timeSum; 
    double _maxP[32]; 
    double _timeP[32];
    double _RWM_T;
    //double _BeamT0;

    bool fMC;
    bool f_isrun3;
    bool f_isnumi;
    bool fUsePID;
    bool fIsEXT;
    bool fSaveExtraInfo;

    float f_shiftoffset;
  
    float f_ccnd1_a;
    float f_ccnd1_b;
    float f_ccnd2_a;
    float f_ccnd2_b;
    float f_ccnd3_a;
    float f_ccnd3_b;
    float f_ccnd3_c;
    float f_ccnd3_d;
    //float f_ccnd4_a;
    //float f_ccnd4_b;
    int _run;

    //Useful info
    float _pmt_size;
    float _Med_TT3;

    //MC-specific variables
    float _sim_time_offset; // variable stored in TTree
    float _mc_interaction_time;    // variable stored in TTree
    float _interaction_time_modulo; // variable stored in TTree
    float _interaction_time_abs;    // variable stored in TTree
    float _spill_time;
    float _decay_time;
    float _prop_time;
    float _vtxt;
    float _mc_offset;

    //re-calibration variables
    std::vector<float> _PMT_time;
    std::vector<float> _PMT_timeProp; 

    Float_t	f_mcflux_dk2gen; // distance from decay to ray origin
    Float_t	f_mcflux_gen2vtx; // distance from ray origin to event vtx


    //for resimulating beam structure
    double fTimeBetweenBuckets;  //time between buckets
    double fBucketTimeSigma;  //how wide is distribution in bucket
    double fNBucketsPerBatch;   
    double fNFilledBucketsPerBatch;
    double fGlobalOffset;  //always displaced by this (in ns)
    std::vector<double> fbi;
    std::vector<int>    fDisallowedBatchMask;  ///< disallow individual batches
    std::vector<double> fCummulativeBatchPDF;  ///< summed prob for batches
    TRandom* fRndmGen;
    
    double _nuvtx_x, _nuvtx_y, _nuvtx_z;

    Float_t _evtDeltaTimeNS, _evtTimeNS;
    Float_t f_sim_time;
    Float_t f_sim_time_nospill;
    Float_t f_sim_deltatime;
    Float_t	f_truth_nu_pos[4]; // X,Y,Z,T
    double  calib[32];

    //const int samples_numi = 1500;

    std::map<unsigned int, unsigned int> _pfpmap;

    void getBeamWF(std::vector<double> *wf_w_03);
    void getPMTwf(std::vector<double> wfsum,
		  std::vector< std::vector<double> > _wf_v);
    //int& tickSum, int tickP[32],
    //double& maxSum, double& timeSum, double maxP[32], double timeP[32]);
    
    void nsbeamtiming(art::Event const& e, std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v);
    void CalculateCPDF(std::vector<double> bi);
    void get_sim_time(art::Event const& e);
    double TimeOffset();


  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  NeutrinoTiming::NeutrinoTiming(const fhicl::ParameterSet& p)
  {

    fPFPproducer = p.get< art::InputTag >("PFPproducer");
    fT0producer  = p.get< art::InputTag >("T0producer" );
    fPMTWFproducer  = p.get< art::InputTag >("nstimePMTWFproducer", "pmtreadout:OpdetBeamHighGain");
    fLogicWFproducer = p.get< art::InputTag >("LogicWFproducer", "pmtreadout:UnspecifiedLogic");
    fFLASHproducer  = p.get< art::InputTag >("FLASHproducer" );
    fSpacePointproducer = p.get< art::InputTag >("SpacePointproducer");

    f_shiftoffset = p.get<float>("ShiftOffset", 0);
    f_isnumi  = p.get<bool>("isNuMI", false);
    f_isrun3  = p.get<bool>("isRun3", false);
    fMC       = p.get<bool>("isMC", false);
    fUsePID   = p.get<bool>("usePID", false);
    fIsEXT    = p.get<bool>("isEXT", false);
    fSaveExtraInfo   = p.get("SaveExtraInfo", true);

    //Empirical correction
    f_ccnd1_a = p.get<float>("ccnd1_a", 0.529594);
    f_ccnd1_b = p.get<float>("ccnd1_b", 7.13804);
    f_ccnd2_a = p.get<float>("ccnd2_a", 0.068752);
    f_ccnd2_b = p.get<float>("ccnd2_b", 2.32023);
    f_ccnd3_a = p.get<float>("ccnd3_a", 0.4697);
    f_ccnd3_b = p.get<float>("ccnd3_b", 0.004233);
    f_ccnd3_c = p.get<float>("ccnd3_c", 0.000001006);
    f_ccnd3_d = p.get<float>("ccnd3_d", -0.195);
    //f_ccnd4_a = pset.get<float>("ccnd4_a", 0);
    //f_ccnd4_b = pset.get<float>("ccnd4_b", 0);

    //Global MC offset
    _mc_offset = p.get<float>("MCOffset", 0);

    //MC-specific ntuple variables
    f_sim_time           = 0;
    f_sim_time_nospill   = 0;
	  f_sim_deltatime      = 0;
    _mc_interaction_time = 0;
    _sim_time_offset     = 0;




    //Beam parameters
    fTimeBetweenBuckets     = p.get<float>("TimeBetweenBuckets",1e9/53.103e6);
    fBucketTimeSigma        = p.get<float>("BucketTimeSigma",0.750);
    fNBucketsPerBatch       = p.get<int>("NBucketsPerBatch",84);
    fNFilledBucketsPerBatch = p.get<int>("NFilledBucketsPerBatch",81);
    fGlobalOffset           = p.get<double>("GlobalOffset",0);
    

    fbi = p.get<std::vector<double>>("BatchIntensities",{1,1,1,1,1,1}); // all ones would be equal batches, set last to lower to emulate slip stacking
    if(fbi.size()==0){fbi = {1,1,1,1,1,1};}
    fRndmGen = new TRandom3(0);
    CalculateCPDF(fbi);

    //ns beam timing plots for validation
    art::ServiceHandle<art::TFileService> tfs;

    if(f_isnumi){
      H_time= tfs->make<TH1F>("H_time","Time PMT",2000, 0,20000);
      H_maxH= tfs->make<TH1F>("H_maxH","Max amplitude",2400,1800,2100);
      H_t0_Beam= tfs->make<TH1F>("H_t0_Beam","T_0 beam",300,0,150);
      H_TimeVsPh= tfs->make<TH2F>("H_TimeVsPh","H_TimeVsPh",  100, -50,50,  100, 0,1500);
      H_TruthTime = tfs->make<TH1F>("H_Truthtime","H_Truthtime",  100, 2000, 5000);
      H_SimTime   = tfs->make<TH1F>("H_SimTime","H_SimTime",  100, 2000, 7000);
      H_ns_time   = tfs->make<TH1F>("H_ns_time","H_ns_time",  100, 2000, 7000);
    }else{
      H_time      = tfs->make<TH1F>("H_time","Time PMT",500, 0,6000);
      H_maxH      = tfs->make<TH1F>("H_maxH","Max amplitude",800,2000,2100);
      H_t0_Beam   = tfs->make<TH1F>("H_t0_Beam","T_0 beam",800,3000,8450);
      H_TimeVsPh  = tfs->make<TH2F>("H_TimeVsPh","H_TimeVsPh",  100, -50,50,  100, 0,500);
      H_TruthTime = tfs->make<TH1F>("H_Truthtime","H_Truthtime",  100, 2000, 5000);
      H_SimTime   = tfs->make<TH1F>("H_SimTime","H_SimTime",  100, 2000, 7000);
      H_ns_time   = tfs->make<TH1F>("H_ns_time","H_ns_time",  100, 2000, 7000);
    }

  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void NeutrinoTiming::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void NeutrinoTiming::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {

    //std::cout << "[NeutrinoTimingDebug] Start ns timing reconstruction" << std::endl;


    _run = e.run();

    const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
    for (int i=0; i<32; i++){
      calib[i] = gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain");
      //std::cout << "[NeutrinoTimingDebug] calib[i] for i "  << i << " is " << calib[i] << std::endl;
    }

    // load PMT waveforms
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle;
    e.getByLabel( fPMTWFproducer, wf_handle );
    std::cout << "[NeutrinoTimingDebug] wf_handle label: " << fPMTWFproducer << std::endl;

    //art::Handle< std::vector< raw::OpDetWaveform > > wf_handle;
    //if(!fMC) {e.getByLabel( "pmtreadout:OpdetBeamHighGain", wf_handle );}
    //else{e.getByLabel( "mixer:OpdetBeamHighGain", wf_handle );}//MC in this instead

    if(!wf_handle.isValid()) {
      std::cerr<<"BeamTiming: No pmtreadout:OpdetBeamHighGain! Skip neutrino timing."<<std::endl;
      return;
    }

    // load logic waveforms
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle_beam;
    e.getByLabel( fLogicWFproducer, wf_handle_beam );
    
    if(!wf_handle_beam.isValid()) {
      std::cerr<<"\033[93m[WARNING]\033[00m no raw::OpDetWaveform Logic Signal. Skip neutrino timing"<<std::endl;
      return;
    }

    /*
     load waveforms into vectors
    */
    //std::cout << "[NeutrinoTimingDebug] load waveforms into vectors " << std::endl;

    //clear PMT waveform
    auto wfsum = std::vector<double>(1500,0);
    auto _wf_v = std::vector< std::vector<double> >(32,std::vector<double>(1500,0));
    
    for (size_t i=0; i < wf_handle->size(); i++) {
      auto const wf = wf_handle->at(i);
      auto ch = wf.ChannelNumber();
      if (ch >= 32) continue;
      
      for (size_t n=0; n < wf.size(); n++){
        if (n < 1500){
          _wf_v[ch][n] = (wf_handle->at(i))[n];
          wfsum[n] += wf[n];
        }
      }// for all channels
    }// for PMT all waveforms

    //std::cout << "[NeutrinoTimingDebug] calling getPMTwf() " << std::endl;
    getPMTwf(wfsum, _wf_v);
    
    /*
      load Logic waveform
     */
    //std::cout << "[NeutrinoTimingDebug] load logic waveform " << std::endl;

    std::vector<double> *wf_w_03 = new std::vector<double>;

    if(!fMC){ //Load beam waveform if beam on data
      for (size_t i=0; i < wf_handle_beam->size(); i++) {
        auto const wf = wf_handle_beam->at(i);
        auto ch = wf.ChannelNumber();
        //std::cout<<"Channel: "<<ch<<std::endl;
        //if (ch != 39) continue;
        if ( (ch != 39 && !f_isnumi) || (ch != 37 && f_isnumi) ) continue;
        for (size_t n=0; n < wf.size(); n++){
          if (n < 1500){
            wf_w_03->emplace_back((wf_handle_beam->at(i))[n]);
          }
        }// for all channels
      }// for all waveforms
      getBeamWF(wf_w_03);
    }
    
    // load tracks previously created for which T0 reconstruction is requested                                                                                                                                

    //std::cout << "[NeutrinoTimingDebug] load T0 tagged information " << std::endl;

    art::Handle<std::vector<anab::T0> > t0_h;
    e.getByLabel( fT0producer , t0_h );

    // make sure tracks look good                                                                                                                                                                             
    if(!t0_h.isValid()) {
      std::cerr<<"\033[93m[WARNING]\033[00m no anab::T0 for flash-matching. Skip flash-matching"<<std::endl;
      return;
    }
    
    art::ValidHandle<std::vector<recob::PFParticle>> pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    
    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
    auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);  
    art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);
    
    
    // figure out which PFP is the neutrino
    //size_t nupfp = 0;
    for (auto pfp : slice_pfp_v) {
      
      if (pfp->IsPrimary() == false) continue;
      
      if ( (pfp->PdgCode() == 12) || (pfp->PdgCode() == 14) ) {
	
	//nupfp = pfp->Self();
	
	auto vtx = pfp.get<recob::Vertex>();
	if (vtx.size() != 1) {
	  std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	  return;
	}
	std::cout << "Neutrino! \n";
	// save vertex to array
	TVector3 nuvtx;
	Double_t xyz[3] = {};
	vtx.at(0)->XYZ(xyz);
	nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);
	_nuvtx_x = nuvtx.X();
	_nuvtx_y = nuvtx.Y();
	_nuvtx_z = nuvtx.Z();
	
      }// if neutrino
    }// for all PFPs in slice
    
    
    // now build vectors of PFParticles, space-points, and hits for this slice
    // std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit> > > hit_v_v; 
    
    // fill map: pfparticle Self() -> index/key
    _pfpmap.clear();
    for (unsigned int p=0; p < pfp_h->size(); p++) _pfpmap[pfp_h->at(p).Self()] = p;
    
    // loop through all PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {
      
      auto const& pfp = pfp_h->at(p);
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      
      // only primary PFPs have a flash-match score
      if (pfp.IsPrimary() == false) continue;
      
      // here only interested in the neutrino slice
      auto PDG = fabs(pfp.PdgCode());
      if ( (PDG != 12) && (PDG != 14) ) continue;
      
      AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);
      
      // go through these pfparticles and fill info needed for matching
      for (size_t i=0; i < pfp_ptr_v.size(); i++) {    
	      auto key = pfp_ptr_v.at(i).key();
	      // recob::PFParticle ipfp = *pfp_ptr_v.at(i);
	      // pfp_v.push_back(ipfp);  
	      //auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
	      for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	        auto const& spkey = spacepoint_ptr_v.at(sp).key();
	        const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	        for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	          hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
	        }// for all hits associated to this spacepoint
	      }// fpr all spacepoints  
	      spacepoint_v_v.push_back( spacepoint_ptr_v );
	      hit_v_v.push_back( hit_ptr_v );
      }// for all pfp pointers
      //flashmatch::SliceCandidate slice(pfp_ptr_v, spacepoint_v_v, hit_v_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower);
      
    }// for all PFPs

    if (fMC){
      // Generator Info
      auto const& mclist = e.getProduct<std::vector<simb::MCTruth>>("generator");
      if (mclist.size()>0) {

        simb::MCTruth const& mctruth = mclist[0];
        if (mctruth.NeutrinoSet()) {

          simb::MCNeutrino nu = mctruth.GetNeutrino();
          const TLorentzVector& position = nu.Nu().Position(0);
          f_truth_nu_pos[0] = position.X(); // cm
          f_truth_nu_pos[1] = position.Y();
          f_truth_nu_pos[2] = position.Z();
          f_truth_nu_pos[3] = position.T();
        }
      }
   }

    
    std::cout << "[NeutrinoTiming] run ns beam timing " << std::endl;
    nsbeamtiming(e, spacepoint_v_v);

    _interaction_time_modulo = _evtDeltaTimeNS;
    _interaction_time_abs = _evtTimeNS;

    if(fMC){
      _sim_time_offset = f_sim_deltatime;
      _mc_interaction_time = f_sim_time;
    }

    //std::cout << "[NeutrinoTimingDebug] interaction_time_abs: " << _interaction_time_abs << " interaction_time_modulo: "<< _interaction_time_modulo << " mc_interaction_time: " << _mc_interaction_time << " time_offset: "<< _time_offset << std::endl;
  
  return;
  }
  
  void NeutrinoTiming::analyzeEvent(art::Event const &e, bool fData)
  {
 
    return;
  }

  void NeutrinoTiming::setBranches(TTree* _tree)
  {
    _tree->Branch("interaction_time_abs",&_interaction_time_abs,"interaction_time_abs/F");
    _tree->Branch("interaction_time_modulo",&_interaction_time_modulo,"interaction_time_modulo/F");

    if(fSaveExtraInfo){
      _tree->Branch("pmt_size",&_pmt_size,"pmt_size/I");
      _tree->Branch("Med_TT3", &_Med_TT3, "Med_TT3/F");
      _tree->Branch("pmt_timeProp", &_PMT_timeProp, "timeProp/F");
      _tree->Branch("pmt_time", &_PMT_time, "time/F");
    }
    

    if(fMC){
      _tree->Branch("mc_interaction_time",&_mc_interaction_time,"mc_interaction_time/F");
      _tree->Branch("sim_time_offset",&_sim_time_offset,"sim_time_offset/F");

      if(fSaveExtraInfo){
        _tree->Branch("spill_time",&_spill_time,"spill_time/F");
        _tree->Branch("decay_time", &_decay_time, "decay_time/F");
        _tree->Branch("propagation_time", &_prop_time, "propagation_time/F");
      }
    }else{
      if (fSaveExtraInfo) _tree->Branch("beam_t0",&_RWM_T,"beam_t0/D");
    }
    //_tree->Branch("Ph_Tot", &f_Ph_Tot);
  }
  
  void NeutrinoTiming::resetTTree(TTree* _tree)
  {
    _interaction_time_abs    = std::numeric_limits<float>::lowest();
    _interaction_time_modulo = std::numeric_limits<float>::lowest();

    if(fSaveExtraInfo){
      _pmt_size      = std::numeric_limits<int>::lowest();
      _Med_TT3       = std::numeric_limits<float>::lowest();
      _PMT_timeProp.clear();
      _PMT_time.clear();
    }

    if(fMC){
      _mc_interaction_time   = std::numeric_limits<float>::lowest();
      _sim_time_offset       = std::numeric_limits<float>::lowest();
      _spill_time            = std::numeric_limits<float>::lowest();
    }
    else{
      _RWM_T = std::numeric_limits<float>::lowest();
    }
    
  }

  void NeutrinoTiming::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
    auto daughters = pfp_ptr->Daughters();
  
    pfp_v.push_back(pfp_ptr);
  
    std::cout << "\t PFP w/ PdgCode " << pfp_ptr->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  
    for(auto const& daughterid : daughters) {

      if (_pfpmap.find(daughterid) == _pfpmap.end()) {
	        std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
	        continue;
      }
    
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    
      AddDaughters(pfp_ptr, pfp_h, pfp_v);
    
    }// for all daughters


  
    return;
  }


  void NeutrinoTiming::getPMTwf(std::vector<double> wfsum, std::vector< std::vector<double> > _wf_v)
  {

    //analyze waveforms
    int samples_64 = 5;
    int samples = 500;

    if(f_isnumi){
      samples_64 = 23;
      samples = 1500;
    }

    //std::cout << "[NeutrinoTiming:getPMTwf()] samples_64: "<< samples_64 << " samples: "<< samples << std::endl;

 

    
    //=======================================================================================================
    //=======================================================================================================
    //std::vector<std::vector<double>> Help_wf_v(32, std::vector<double> (samples));
    //std::vector<double> x_wf_v(samples), Raw_wf_v(samples), Base_wf_v(samples), Norm_wf_v(samples);
    double Help_wf_v[32][1500];
    double x_wf_v[1500], Raw_wf_v[1500], Base_wf_v[1500], Norm_wf_v[1500];
    double maxZ,max0,base;   int basebinmax,tick;
    double tca,tcb,tcc, TT[32], max[32];
    
    int TF,TB,tickF,tickB,FB, Nss,is;
    double  maxZhelp1,maxZhelp2,maxZhelp3, tickFit1,tickFit2;
    
    //Raw waveform saturation
    int saturation=4094;
    //Parameters for discrete saturated WF reconstruction (1st step)
    double Frac[100]={1.0, 0.951931, 0.93356, 0.838637, 0.719408, 0.701042, 0.565673, 0.44655, 0.39447, 0.352336, 0.28716, 0.245364, 0.216771, 0.194888, 0.178976, 0.16844, 0.161732, 0.157486, 0.154762, 0.153136, 0.152468, 0.152299, 0.152147, 0.151456, 0.150069, 0.148089, 0.145565, 0.142516, 0.139369, 0.136467, 0.133871, 0.13141, 0.128926, 0.126237, 0.12313, 0.119611, 0.115946, 0.112453, 0.111706, 0.109225, 0.106281, 0.103499, 0.100926, 0.0985215, 0.0961512, 0.0938219, 0.0917209, 0.0898817, 0.0882937, 0.0868338, 0.0854753, 0.0841353, 0.0827237, 0.081198, 0.0794796, 0.077565, 0.0756475, 0.0738863, 0.0722821, 0.0708473, 0.0695119, 0.0682504, 0.0672023, 0.0663549, 0.0656337, 0.0649918, 0.0645003, 0.0641535, 0.0638046, 0.0633435, 0.0627506, 0.0621379, 0.0615464, 0.0609178, 0.0601846, 0.0592098, 0.0580465, 0.0568861, 0.0559024, 0.0550731, 0.0541904, 0.0532532, 0.0524181, 0.0517606, 0.0512326, 0.0507392, 0.0502093, 0.0495968, 0.0488915, 0.0480173, 0.0470195, 0.0459744, 0.0448855, 0.0437359, 0.0425199, 0.0412832, 0.0400036, 0.038688, 0.0373173, 0.0358925};
    //double Delay[100]={0.0, 0.0686367, 0.0948574, 0.230854, 0.405159, 0.432604, 0.642806, 0.845613, 0.942555, 1.02616, 1.16775, 1.26923, 1.34523, 1.40802, 1.45675, 1.49068, 1.51306, 1.52757, 1.53702, 1.54272, 1.54507, 1.54567, 1.54621, 1.54865, 1.55359, 1.56069, 1.56984, 1.58104, 1.59279, 1.60379, 1.61377, 1.62337, 1.63319, 1.64398, 1.65665, 1.67129, 1.68688, 1.70208, 1.70537, 1.71644, 1.72981, 1.74271, 1.75486, 1.76644, 1.77806, 1.78969, 1.80036, 1.80987, 1.81819, 1.82594, 1.83324, 1.84054, 1.84831, 1.85684, 1.86658, 1.87763, 1.88891, 1.89947, 1.90926, 1.91815, 1.92656, 1.93462, 1.9414, 1.94695, 1.95171, 1.95598, 1.95928, 1.96162, 1.96398, 1.96711, 1.97117, 1.97539, 1.97951, 1.98391, 1.98909, 1.99605, 2.00448, 2.01302, 2.02038, 2.02665, 2.03342, 2.04069, 2.04726, 2.0525, 2.05674, 2.06073, 2.06505, 2.0701, 2.07597, 2.08334, 2.09188, 2.10099, 2.11065, 2.12106, 2.13232, 2.14403, 2.15646, 2.16957, 2.18362, 2.19868 };
    //Fit function parameters for saturated WF reconstruction (2nd step)
    double pLL[8]={3.93256,4.31002,2.44182,5.12491,0.830928,0.231375,50.9081,-2.69014};
    
    //=====================================================================================================
    //=======================================================================================================
    
    std::vector<float> max_pmt_v(32,0);
    std::vector<int>   max_tick_v(32,0);
    
    for(int q=0; q<32; q++){
      for(int i=0; i<samples; i++){
	      Help_wf_v[q][i]=_wf_v[q][i];
	      //if (i > 100 && i < 110) std::cout << "[NeutrinoTimingDebug::getPMTwf()] PMT " << q << " @ tick " << i << " has ADC " << Help_wf_v[q][i] << std::endl;
	      if (Help_wf_v[q][i] > max_pmt_v[q]) { max_pmt_v[q] = Help_wf_v[q][i]; max_tick_v[q] = i; }
      }
    }

    //for (size_t nn=0; nn < 32; nn++) {
      //std::cout << "[NeutrinoTimingDebug] PMT " << nn << " has maximum tick @ " <<  max_tick_v[nn] << " with value " << max_pmt_v[nn] << "." << std::endl;
    //}

    _wf_v.clear(); _wf_v.shrink_to_fit();
    
    //-----------------------------------------------------------------------------------------------------
    for(int q=0; q<32; q++){
      TT[q]=-9999.;
      max[q]=-9999.;
      maxZ=0.;
      max0=0.;
      base=0.;
      tick=0;
      tickB=0;
      tickF=0;
      TF=0;
      TB=0;
      //Getting raw waveform (Raw_wf_v[i]) only for i<500 since the beam window is between 3 and 5 us -> [i>3*64 && i<5*64]
      for(int i=0; i<samples; i++){
	      x_wf_v[i]=i*1.0;
	      Raw_wf_v[i]=Help_wf_v[q][i];
      }
      //Getting raw wf max amplitude and max amp tick
      for(int i=3*64; i<samples_64*64; i++){if(maxZ<Raw_wf_v[i]){maxZ=Raw_wf_v[i]; tick=i;}}
      //Baseline removal
      TH1F *basehelp= new TH1F("basehelp","basehelp",400, 1900,2200);
      basebinmax=0; for(int i=0; i<3*64; i++){basehelp->Fill(Raw_wf_v[i]);}
      basebinmax=basehelp->GetMaximumBin(); base=basehelp->GetXaxis()->GetBinCenter(basebinmax);
      basehelp->Delete();
      //Getting wf max amp after baseline removal (this is proportional to number of Photons in the rising endge)
      //getting wf baseline subtracted and wf baseline subtracted and normalized for the max amp.
      for(int i=0; i<samples; i++){max0=maxZ-base;
	      Base_wf_v[i]=Raw_wf_v[i]-base; Norm_wf_v[i]=Base_wf_v[i]/max0;}
      //fitting the normalized baseline subtracted wf
      TGraph *gr = new TGraph(samples,x_wf_v,Norm_wf_v);
      TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power(([0]-x)/[1],4))",tick-10, tick);
      fit->SetParameters(tick,2,1);  gr->Fit("fit","Q","",tick-10, tick);
      tca=fit->GetParameter(0);  tcb=fit->GetParameter(1);  tcc=fit->GetParameter(2);
      //timing is the risign edge half height
      TT[q]=(tca-abs(tcb*TMath::Power(-log(0.5/tcc),0.25)))/0.064; max[q]=max0;
      //----------------------------------------------
      //check for saturated wf
      if(maxZ<=saturation){
	      TT[q]=TT[q];
	      max[q]=max[q];
	      //std::cout << "[NeutrinoTimingDebug] not saturated!" << std::endl;
      }
      else if(maxZ>saturation) {
	      //counting the number of ticks above the saturation
	      for(int i=3*64; i<samples_64*64; i++){
	        if(TF==0){if(Raw_wf_v[i+1]>4094 && Raw_wf_v[i]<=4094){tickF=i; TF=1;}}
	        if(TB==0){if(Raw_wf_v[i]>4094 && Raw_wf_v[i+1]<=4094){tickB=i; TB=1;}}
        }
	      FB=tickB-tickF;  if(FB>99){FB=99;}
	      //amplitude discrete correction
	      maxZhelp1=maxZ/Frac[FB]; tick=tickF; Nss=0; is=0;
        for(int i=3*64; i<samples_64*64; i++){if(Raw_wf_v[i]<4095){Nss=Nss+1;}}
      
        double txSS[1500],tySS[1500],txSS2[1500],tySS2[1500];
      
        for(int i=3*64; i<samples_64*64; i++){if(Raw_wf_v[i]<4095){txSS[is]=i*1.0; tySS[is]=Raw_wf_v[i]/maxZhelp1; is=is+1;}}
        TGraph *g1 = new TGraph(Nss,txSS,tySS);
        //for (int uu=0; uu < (7*64); uu++)
        //std::cout << "[NeutrinoTimingDebug] value @ tick " << uu << " is " << tySS[uu] << std::endl;
        TF1 *fitS1 = new TF1("fitS1","[9]*(exp(-TMath::Power(([0]-(x-[8]))/[1],4))*0.5*(TMath::Erf(-(x-[8])-[7])+1.0)+([5]+[4]*exp(-TMath::Power(([2]-(x-[8]))/[3],2)))*exp((-(x-[8]))/[6])*0.5*(TMath::Erf([7]+(x-[8]))+1.0))",tick-30, tick+250);
        fitS1->SetParameters(pLL[0],pLL[1],pLL[2],pLL[3],pLL[4],pLL[5],pLL[6],pLL[7],tick,1.);
        for(int i=0; i<8; i++){fitS1->FixParameter(i,pLL[i]);} g1->Fit("fitS1","Q","",tick-30, tick+250);
        tickFit1=fitS1->GetParameter(8); maxZhelp2=fitS1->GetParameter(9);  maxZhelp3=maxZhelp1/maxZhelp2;
        //amplitude fit correction

        for(int i=0; i<Nss; i++){txSS2[i]=txSS[i]; tySS2[i]=tySS[i]/maxZhelp2;}
        TGraph *g2 = new TGraph(Nss,txSS2,tySS2);

        TF1 *fitS2 = new TF1("fitS2","exp(-TMath::Power(([0]-(x-[8]))/[1],4))*0.5*(TMath::Erf(-(x-[8])-[7])+1.0)+([5]+[4]*exp(-TMath::Power(([2]-(x-[8]))/[3],2)))*exp((-(x-[8]))/[6])*0.5*(TMath::Erf([7]+(x-[8]))+1.0)",tick-30, tick+250);

        fitS2->SetParameters(pLL[0],pLL[1],pLL[2],pLL[3],pLL[4],pLL[5],pLL[6],pLL[7],tickFit1);
        for(int i=0; i<8; i++){fitS2->FixParameter(i,pLL[i]);}

        g2->Fit("fitS2","Q","",tick-30, tick+250);  tickFit2=fitS2->GetParameter(8);
        //timing is the risign edge half height
        TT[q]=tickFit2/0.064; max[q]=maxZhelp3;
      }// if not saturated
      //-------------------------------------------------------------------------------------------------------
      H_time->Fill(TT[q]);
    }
    //-------------------------------------------------------------------------------------------------------

    for(int q=0; q<32; q++){
      _maxP[q]=max[q];
      _timeP[q]=TT[q];
      //std::cout << "[NeutrinoTimingDebug] : PMT "  << q << " has maximum @ tick " << _timeP[q] << " with value " << _maxP[q] << std::endl;
    } //only two variables needed
    
  }
  
  
  void NeutrinoTiming::getBeamWF(std::vector<double> *wf_w_03)
  {
    //=======================================================================================================
    //=======================================================================================================

    double beamBase, wx[500], wy[500], BBmax,pca,pcb,pcc, TT = -99999.0;

    int Btick,tickMax;
    //=======================================================================================================
    //-------baseline calculation-------------------------------------------
    beamBase=0.; BBmax=0.; Btick=0; TT=-9999.;
    for(int i=0; i<4*64; i++){beamBase=beamBase+wf_w_03->at(i);}
    beamBase=beamBase/(4.0*64.0);
    //baseline subtraction
    for(int i=0; i<500; i++){
      wx[i]=i*1.0; 
      wy[i]=wf_w_03->at(i)-beamBase;
      //max amplitude
      if(BBmax<wy[i]){BBmax=wy[i]; Btick=i;}
    }
    H_maxH->Fill(BBmax);
    double BBmax_threshold_l0 = 2000;
    double BBmax_threshold_high = 2100;
    if(f_isnumi) {
      BBmax_threshold_l0   = 1800;
    }
    if(BBmax>BBmax_threshold_l0 && BBmax<BBmax_threshold_high){
    //wf normalization
    for(int i=0; i<500; i++){wy[i]=wy[i]/BBmax;}
    //wf max check
    tickMax=0;
    for(int i=Btick-20; i<Btick+10; i++){if(wy[i-1]<1 && wy[i]==1){tickMax=i;}}
    //wf fit
    TGraph *gr0 = new TGraph(500, &(wx[0]),&(wy[0]));
    TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power((x-[0])/[1],4))",tickMax-6, tickMax);
    fit->SetParameters(tickMax,2,1);     gr0->Fit("fit","Q","",tickMax-6, tickMax);
    pca=fit->GetParameter(0); pcb=fit->GetParameter(1); pcc=fit->GetParameter(2);
    //timing is the rising edge half height
    TT=(pca-abs(pcb*TMath::Power(-log(0.5/pcc),0.25)))/0.064;
    _RWM_T = TT;
    std::cout << "[NeutrinoTimingDebug] RWM time (_RMW_T) : " << _RWM_T << std::endl;
    H_t0_Beam->Fill(TT);
  }
}

  void NeutrinoTiming::nsbeamtiming(art::Event const& e, std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v)
  { 

    Float_t x =	_nuvtx_x;
    Float_t y =	_nuvtx_y;
    Float_t z =	_nuvtx_z;

    //std::cout << "[NeutrinoTimingDebug] vertex @ (" << x << ", " << y << ", " << z << ")" << std::endl;
    
    std::vector<float> *sps_x = new std::vector<float>;
    std::vector<float> *sps_y = new std::vector<float>;
    std::vector<float> *sps_z = new std::vector<float>;

    if(!fUsePID){ //default
      for (size_t s=0; s < spacepoint_v_v.size(); s++) {
        auto sps_v = spacepoint_v_v[s];
        for (size_t ss=0; ss < sps_v.size(); ss++) {
	        auto sps = sps_v[ss];
	        sps_x->push_back( sps->XYZ()[0] );
	        sps_y->push_back( sps->XYZ()[1] );
	        sps_z->push_back( sps->XYZ()[2] );

        }
      }
      //std::cout << "[NeutrinoTimingDebug] there are " << sps_z->size() << " spacepoints" << std::endl;
    }
    
    
    double max[32],time[32];
    double BeamT0 = -99999.;
    for (int ii=0; ii < 32; ii++){
      time[ii] = _timeP[ii];
      max[ii] = _maxP[ii];
    }
    
    //getPMTwf(e,max,time);
    if(!fMC && !fIsEXT) BeamT0 = _RWM_T;
    else BeamT0 = 0;
    
    double PMT0[3]={-11.4545, -28.625, 990.356};  double PMT1[3]={-11.4175, 27.607, 989.712};
    double PMT2[3]={-11.7755, -56.514, 951.865};  double PMT3[3]={-11.6415, 55.313, 951.861};
    double PMT4[3]={-12.0585, -56.309, 911.939};  double PMT5[3]={-11.8345, 55.822, 911.065};
    double PMT6[3]={-12.1765, -0.722, 865.599};   double PMT7[3]={-12.3045, -0.502, 796.208};
    double PMT8[3]={-12.6045, -56.284, 751.905};  double PMT9[3]={-12.5405, 55.625, 751.884};
    double PMT10[3]={-12.6125, -56.408, 711.274}; double PMT11[3]={-12.6615, 55.8, 711.073};
    double PMT12[3]={-12.6245, -0.051, 664.203};  double PMT13[3]={-12.6515, -0.549, 585.284};
    double PMT14[3]={-12.8735, 55.822, 540.929};  double PMT15[3]={-12.6205, -56.205, 540.616};
    double PMT16[3]={-12.5945, -56.323, 500.221}; double PMT17[3]={-12.9835, 55.771, 500.134};
    double PMT18[3]={-12.6185, -0.875, 453.096};  double PMT19[3]={-13.0855, -0.706, 373.839};
    double PMT20[3]={-12.6485, -57.022, 328.341}; double PMT21[3]={-13.1865, 54.693, 328.212};
    double PMT22[3]={-13.4175, 54.646, 287.976};  double PMT23[3]={-13.0075, -56.261, 287.639};
    double PMT24[3]={-13.1505, -0.829, 242.014};  double PMT25[3]={-13.4415, -0.303, 173.743};
    double PMT26[3]={-13.3965, 55.249, 128.354};  double PMT27[3]={-13.2784, -56.203, 128.18};
    double PMT28[3]={-13.2375, -56.615, 87.8695}; double PMT29[3]={-13.5415, 55.249, 87.7605};
    double PMT30[3]={-13.4345, 27.431, 51.1015};  double PMT31[3]={-13.1525, -28.576, 50.4745};
    double PMT[32][3];    for(int j=0; j<3; j++){ PMT[30][j]=PMT30[j]; PMT[31][j]=PMT31[j];
      PMT[0][j]=PMT0[j];   PMT[10][j]=PMT10[j]; PMT[20][j]=PMT20[j]; PMT[1][j]=PMT1[j];   PMT[11][j]=PMT11[j];
      PMT[21][j]=PMT21[j]; PMT[2][j]=PMT2[j];   PMT[12][j]=PMT12[j]; PMT[22][j]=PMT22[j]; PMT[3][j]=PMT3[j];
      PMT[13][j]=PMT13[j]; PMT[23][j]=PMT23[j]; PMT[4][j]=PMT4[j];   PMT[14][j]=PMT14[j]; PMT[24][j]=PMT24[j];
      PMT[5][j]=PMT5[j];   PMT[15][j]=PMT15[j]; PMT[25][j]=PMT25[j]; PMT[6][j]=PMT6[j];   PMT[16][j]=PMT16[j];
      PMT[18][j]=PMT18[j]; PMT[28][j]=PMT28[j]; PMT[9][j]=PMT9[j];   PMT[19][j]=PMT19[j]; PMT[29][j]=PMT29[j];
      PMT[26][j]=PMT26[j]; PMT[7][j]=PMT7[j];   PMT[17][j]=PMT17[j]; PMT[27][j]=PMT27[j]; PMT[8][j]=PMT8[j];}
    double offset[32]={1.03002, -5.18104, -2.11164, -5.99395, -1.25798, 0.633079, 2.87666, 2.21969, 0.885092, 2.35423,
		       -1.63039, -1.83775, -0.859883, 3.4741, 1.84833, 1.58233, -2.71783, 0, 3.18776, 0.982666, 0.728438, 0.280592, -5.27068,
		       -3.27857, -1.41196, 1.59643, 1.41425, -1.62682, -2.55772, 1.49136, -0.522791, 0.974533};

    if(fMC){
      for(int i=0; i<32; i++){offset[i]=0;}//no need to apply the additional pmt calibration to the MC
      for(int j=0; j<3; j++){//do the PMT remapping for the MC
        double temp = PMT[31][j];
        PMT[31][j] = PMT[30][j];
        PMT[30][j] = PMT[29][j];
        PMT[29][j] = PMT[28][j];
        PMT[28][j] = PMT[27][j];
        PMT[27][j] = PMT[26][j];
        PMT[26][j] = temp;
      }
    }
    
    
    //================================================================================================================
    double gap=18.936;
    double MaxLim=2.5;
    std::vector<int> N_pmt;
    double ccnd1, ccnd2,ccnd3;
    double Ph_Tot, RWM_T, nuToF, DPh,DLh, tPhelp,tp;
    //double Med_TT0,Med_TT1,Med_TT2;
    double Med_TT3=-9999.;
    double TT_merged = -9999.;
    //===================================================================================================================
    //===================================================================================================================
    Ph_Tot=0.;
    N_pmt.clear();
    for(int q=0; q<32; q++) {
      max[q]=max[q]/calib[q];
      //std::cout << "[NeutrinoTimingDebug] max[q] for q = " << q << " is " << max[q] << std::endl;
      if((max[q]>MaxLim && q!=17 && q!=28) && ((time[q]>3000.0 && time[q]<5000.0) || f_isnumi ) ){
        N_pmt.push_back(q); 
        Ph_Tot=Ph_Tot+max[q]; 
        //pmtid->push_back(q);
      }
    } 
   //--------------------------------------------------------------------------------------------------------------------
    std::cout << "[NeutrinoTimingDebug] PMT size : " << N_pmt.size() << std::endl;
    
    if(N_pmt.size()>2){
      _pmt_size = N_pmt.size();
      RWM_T=BeamT0;
      //std::cout << "[NeutrinoTimingDebug] RWM_T : " << RWM_T << std::endl;
      double dist = z;
      
      if(f_isnumi) {
        TVector3 target_dir(-0.46, -0.05, -0.885);
        double min_a = -122.86902944472968;  
        double min_b = 80.60659897339974; 
        double min_c = 59.34119182916038;
        dist = ( (min_a-x)*target_dir[0] + (min_b-y)*target_dir[1] + (min_c-z)*target_dir[2] ) / sqrt(target_dir[0]*target_dir[0] + target_dir[1]*target_dir[1] + target_dir[2]*target_dir[2] );
      }
      nuToF=dist*0.033356;
      //std::cout << "[NeutrinoTimingDebug] nuToF : " << nuToF << std::endl;
      std::vector<double> timeProp = std::vector<double>(N_pmt.size(),0);
      for(uint i=0; i<N_pmt.size(); i++){
        tp=5000000000.0;
        for(uint j=0; j<sps_x->size(); j++){
          DPh=abs(sqrt(TMath::Power(x-sps_x->at(j),2)+TMath::Power(y-sps_y->at(j),2)+TMath::Power(z-sps_z->at(j),2)));
          DLh=abs(sqrt(TMath::Power(PMT[N_pmt.at(i)][0]-sps_x->at(j),2)+TMath::Power(PMT[N_pmt.at(i)][1]-sps_y->at(j),2)+TMath::Power(PMT[N_pmt.at(i)][2]-sps_z->at(j),2)));
          tPhelp=(DPh*0.033356)+(DLh*0.0746);
          if(tPhelp<tp){tp=tPhelp;}
        }
        timeProp[i]=tp;
      
        //std::cout << "[NeutrinoTimingDebug] timeProp: " << timeProp[i] << std::endl;
      }
     
      double TT3_array[32];
      if(f_isrun3 && _run>17200 && _run<17400 && !f_isnumi){
        if(RWM_T>5450){ 
          f_shiftoffset=118.3;
        }
      }
      float RWM_offset = 5700.0 - f_shiftoffset;
      if (f_isnumi) RWM_offset = 40.0;
      if(fMC || fIsEXT) RWM_offset = 0;
      if(fMC) { //Force corrections to 0 for MC for now, calibration upcoming
        f_ccnd1_a = 0;
        f_ccnd1_b = 0;
        f_ccnd2_a = 0;
        f_ccnd2_b = 0;
        f_ccnd3_a = 0;
        f_ccnd3_b = 0;
        f_ccnd3_c = 0;
        f_ccnd3_d = 0;
      }

      
      for(uint i=0; i<N_pmt.size(); i++){
        
	      ccnd1= timeProp[i]*(f_ccnd1_a)-(f_ccnd1_b);
        //std::cout << "[NeutrinoTimingDebug] ccnd1 "<< ccnd1 << std::endl;
	      ccnd2= max[N_pmt.at(i)]*(f_ccnd2_a)-(f_ccnd2_b);
        ///std::cout << "[NeutrinoTimingDebug] ccnd2 " << ccnd2 << std::endl;
	      if(Ph_Tot>150){ccnd3=f_ccnd3_a-f_ccnd3_b*Ph_Tot+f_ccnd3_c*Ph_Tot*Ph_Tot;}
	      else{ccnd3=f_ccnd3_d;}
        //std::cout << "[NeutrinoTimingDebug] ccnd3 " << ccnd3 << std::endl;
	
	      //all the corrections
	      TT3_array[i]=(time[N_pmt.at(i)])-RWM_T+RWM_offset-nuToF-timeProp[i]-offset[N_pmt.at(i)]+ccnd1+ccnd2+ccnd3;
        std::cout << "[NeutrinoTimingDebug] timeProp: " << timeProp[i] << std::endl;
        std::cout << "[NeutrinoTimingDebug] PMT time: " << time[N_pmt.at(i)]<< std::endl;
        std::cout << "[NeutrinoTimingDebug] Corrected PMT " << i << " timing " << TT3_array[i] << std::endl;
        if(fSaveExtraInfo){
          _PMT_timeProp.push_back(timeProp[i]);
          _PMT_time.push_back(time[N_pmt.at(i)]);
        }
      }
      Med_TT3=TMath::Median((Long64_t)N_pmt.size(),TT3_array);
      //Fill a 2d histogram with  TT3_array[i] vs max[N_pmt.at(i)] this is usefull to check for any errors
      for(uint i=0; i<N_pmt.size(); i++){
	      H_TimeVsPh->Fill( TT3_array[i]-Med_TT3, max[N_pmt.at(i)]);
      }
    }// if PMT size > 2


    _evtTimeNS = Med_TT3;
    _Med_TT3 = Med_TT3;
    if(fMC) {
      get_sim_time(e);
      _evtTimeNS = Med_TT3 + f_sim_deltatime + _mc_offset;
    }

    std::cout << "[NeutrinoTimingDebug] Med_TT3 "<< _Med_TT3 <<std::endl;
    H_ns_time->Fill(_evtTimeNS);
  
    //Merge Peaks
    double Shift=3166.9;
    double bunches = 81;
    if(f_isnumi){Shift=11567.87; gap=18.83; bunches=503;}

    //_evtTimeNS = TThelp;
    double TThelp= Med_TT3 - Shift + gap * 0.5;
    std::cout << "[NeutrinoTimingDebug] TThelp : "  << TThelp << std::endl;
    //merge peaks
    if(TThelp>=0 && TThelp < gap * bunches){
      TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
    }
    else {TT_merged=-9999;}

    _evtDeltaTimeNS = TT_merged;
    
    std::cout << "[NeutrinoTimingDebug] evtTimeNS: "<< _evtTimeNS <<std::endl;
    std::cout <<  "[NeutrinoTimingDebug] evtDeltaTimeNS: "<< _evtDeltaTimeNS <<std::endl;
    
  }

  void NeutrinoTiming::get_sim_time(art::Event const& e){
    _vtxt = -9999;
    if (auto mcfluxListHandle = e.getHandle<std::vector<simb::MCFlux>>("generator")){
      //std::cout << "[NeutrinoTimingDebug] Found MCFlux object" << std::endl;

      if (mcfluxListHandle->size()>0) {
            simb::MCFlux const& mcflux = mcfluxListHandle->front();
            f_mcflux_dk2gen = mcflux.fdk2gen; // distance from decay to ray origin
            f_mcflux_gen2vtx = mcflux.fgen2vtx; // distance from ray origin to event vtx
      } 

      _spill_time = TimeOffset();
      _prop_time = (f_mcflux_dk2gen+f_mcflux_gen2vtx)*100*0.033356;

      if(!f_isnumi){
        std::cout<<"Get ReBooNE info"<<std::endl;
        art::Handle<std::vector<BooNEInfo>> flux_handle_bnb;
        std::vector<art::Ptr<BooNEInfo> > bnb_flux;
        e.getByLabel("generator", flux_handle_bnb);

        art::fill_ptr_vector(bnb_flux, flux_handle_bnb);
        if(bnb_flux.size()==0) std::cout<<"no ReBooNE found"<<std::endl;
        else {
          //std::cout<<"Get neutrino start time"<<std::endl;
          art::Ptr<BooNEInfo> this_bnb_flux = bnb_flux.at(0);
          _decay_time = this_bnb_flux->nu_startt; //start time of neutrino, convert to ns
          _vtxt = this_bnb_flux->nu_vtx_t;
          //std::cout<<"Neutrino start time: "<< decay_time << std::endl;
        }

      }else{
        std::cout<<"Get ReDk2Nu info"<<std::endl;
        art::Handle<std::vector<bsim::Dk2Nu>> dk2nu_flux_handle;
	      std::vector<art::Ptr<bsim::Dk2Nu> > dk2nu_flux;
        e.getByLabel("generator",dk2nu_flux_handle);
        art::fill_ptr_vector(dk2nu_flux,dk2nu_flux_handle);
	      if(dk2nu_flux.size()==0){
          std::cout<<"no redk2nu found"<<std::endl;
	      }
	      else{
          art::Ptr<bsim::Dk2Nu> this_dk2nu_flux = dk2nu_flux.at(0);
	        std::vector<bsim::Ancestor> ancestors = this_dk2nu_flux->ancestor;
          _decay_time = ancestors.back().startt; //start time of neutrino 
        }
      }
      f_sim_time = _prop_time + _decay_time + _spill_time;
      f_sim_deltatime = f_sim_time - f_truth_nu_pos[3];
      std::cout << "[NeutrinoTimingDebug] Truth_time: "<< f_truth_nu_pos[3] << " f_sim_time: " << f_sim_time << "  f_sim_deltatime: " << f_sim_deltatime << "  spill_time: "<<  _spill_time <<  "  decay_time: "<< _decay_time << "  propagation_time: "<< _prop_time << std::endl; 
      H_SimTime->Fill(f_sim_time);
      H_TruthTime->Fill(f_truth_nu_pos[3]);
    }

  }

  double NeutrinoTiming::TimeOffset(){
    // calculate in small to large

    // pick a time within a bucket
    double offset = fRndmGen->Gaus(0.0,fBucketTimeSigma);

    // pick a bucket within a batch
    // assume all ~ buckets constant in batch until we have another model
    offset +=  fTimeBetweenBuckets * (double)fRndmGen->Integer(fNFilledBucketsPerBatch);

    // pick a bucket
    bool   disallowed = true;
    size_t ibatch = 0;
    size_t nbatch = fCummulativeBatchPDF.size();
    double r = 2;
    while ( disallowed ) {
      r = fRndmGen->Uniform();
      for (ibatch=0; ibatch<nbatch; ++ibatch) {
        if ( r <= fCummulativeBatchPDF[ibatch] ) break;
      }
      disallowed = ( fDisallowedBatchMask[ibatch] != 0 );
    }
    offset += fTimeBetweenBuckets*(double)fNBucketsPerBatch*(double)ibatch;

    // finally the global offset
    return offset + fGlobalOffset;
  }

  void NeutrinoTiming::CalculateCPDF(std::vector<double> bi){
    fCummulativeBatchPDF.clear();
    double sum = 0;
    size_t nbi = bi.size();
    for (size_t i=0; i < nbi; ++i) {
      sum += bi[i];
      fCummulativeBatchPDF.push_back(sum);
    }
    // normalize to unit probability
    for (size_t i=0; i < nbi; ++i) fCummulativeBatchPDF[i] /= sum;
    // make sure the mask vector keeps up (but never make it smaller)
    // allowing all new batches
    if ( nbi > fDisallowedBatchMask.size() )
      fDisallowedBatchMask.resize(nbi,0);
  }
  
  
  DEFINE_ART_CLASS_TOOL(NeutrinoTiming)
} // namespace analysis

#endif

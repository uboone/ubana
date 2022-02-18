////////////////////////////////////////////////////////////////////////
// Class:       GetPOT
// Plugin Type: analyzer (art v3_01_02)
// File:        GetPOT_module.cc
//
// Generated at Fri Dec  6 09:53:59 2019 by Krishan Mistry using cetskelgen
// from cetlib version v3_05_01.
// A module to get the POT information from a art-root file. 
// It will fill a tree with the run, subrun information and POT information
// It will also create a file for you with the run-subrun info this is useful for
// data where you can give this file to zarko's tool
// of course, this module will only do pot counting, so you should merge this code into your own analyser to insert this info in
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Principal/Run.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "canvas/Utilities/InputTag.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "larcoreobj/SummaryData/POTSummary.h" 

#include "TTree.h" 
#include "TBranch.h" 
#include "TFile.h" 

class GetPOT;

class GetPOT : public art::EDAnalyzer {
public:
  explicit GetPOT(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GetPOT(GetPOT const&) = delete;
  GetPOT(GetPOT&&) = delete;
  GetPOT& operator=(GetPOT const&) = delete;
  GetPOT& operator=(GetPOT&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endSubRun(art::SubRun const &sr) override;

private:

  // POT 
  bool run_pot_counting = false;

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot;
  std::ofstream _run_subrun_list_file;

  std::string _potsum_producer_mc;
  std::string _potsum_producer_data;
  std::string _potsum_instance;

  // Other
  int iteration{0}; // index to count number of entries run over
  std::string _mode;
  bool _is_mc{false};
  bool _is_data{false};
  bool _is_overlay{false};
  bool _cosmic_only{false};

};

GetPOT::GetPOT(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
  // Call appropriate consumes&lt;&gt;() for any products to be retrieved by this module.
  _potsum_producer_mc    = p.get&lt;std::string&gt;("POTSummaryProducerMC");
  _potsum_producer_data  = p.get&lt;std::string&gt;("POTSummaryProducerData");
  _potsum_instance       = p.get&lt;std::string&gt;("POTSummaryInstance");

}

void GetPOT::beginJob() {
  // Implementation of optional member function here.

  art::ServiceHandle&lt;art::TFileService&gt; fs;

  _sr_tree = fs-&gt;make&lt;TTree&gt;("pottree","");
  _sr_tree-&gt;Branch("run",       &_sr_run,       "run/I"      );
  _sr_tree-&gt;Branch("subrun",    &_sr_subrun,    "subrun/I"   );
  _sr_tree-&gt;Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree-&gt;Branch("endtime",   &_sr_endtime,   "endtime/D"  );
  _sr_tree-&gt;Branch("pot",       &_sr_pot,       "pot/D"      );
  _run_subrun_list_file.open ("run_subrun_list.txt", std::ofstream::out | std::ofstream::trunc);

}

void GetPOT::analyze(art::Event const& e) {

  if      (_mode == "EXT")     _cosmic_only = true;
  else if (_mode == "Data")    _is_data     = true;
  else if (_mode == "Overlay"){
    _is_mc       = true;
    _is_overlay  = true;
  } 
  else _is_mc = true;

  std::cout &lt;&lt; "[Analyze] ------------------------------------------- [Analyze]" &lt;&lt; std::endl;
  std::cout &lt;&lt; "[Analyze] Running over entry: " &lt;&lt; iteration &lt;&lt; std::endl;
  iteration++;

  if(_cosmic_only == true) std::cout &lt;&lt; "[Analyze] Running in Cosmic Only Configuration! " &lt;&lt; std::endl;

  if(_is_mc == true)       std::cout &lt;&lt; "[Analyze] Running with Monte Carlo " &lt;&lt; std::endl;

  if(_is_data == true)     std::cout &lt;&lt; "[Analyze] Running with Data " &lt;&lt; std::endl;

  if(_is_data == true)     {run_pot_counting = false; std::cout &lt;&lt; "[Analyze] Do Not Count MC POT" &lt;&lt; std::endl;}

  if(_is_mc == true)       {run_pot_counting = true;  std::cout &lt;&lt; "[Analyze] Count MC POT" &lt;&lt; std::endl;}

}

void GetPOT::endSubRun(art::SubRun const & sr) {

  bool _debug = false;

    if (_debug) std::cout &lt;&lt; "[Analysis::endSubRun] Starts" &lt;&lt; std::endl;

    // Saving run and subrun number on file so that we can run Zarko's script easily
    _run_subrun_list_file &lt;&lt; sr.run() &lt;&lt; " " &lt;&lt; sr.subRun() &lt;&lt; std::endl;

    _sr_run       = sr.run();
    _sr_subrun    = sr.subRun();
    _sr_begintime = sr.beginTime().value();
    _sr_endtime   = sr.endTime().value();

    art::Handle&lt;sumdata::POTSummary&gt; potsum_h;

    // MC
    if (_is_mc) {

        if (_debug) std::cout &lt;&lt; "[Analysis::endSubRun] Getting POT for MC" &lt;&lt; std::endl;

        if(sr.getByLabel(_potsum_producer_mc, potsum_h)) {
            if (_debug) std::cout &lt;&lt; "[Analysis::endSubRun] POT are valid" &lt;&lt; std::endl;
            _sr_pot = potsum_h-&gt;totpot;
        }
        else
            _sr_pot = 0.;
    }

    // Data
    if (_is_data) {

        if (_debug) std::cout &lt;&lt; "[Analysis::endSubRun] Getting POT for DATA, producer " &lt;&lt; _potsum_producer_data &lt;&lt; ", instance " &lt;&lt; _potsum_instance &lt;&lt; std::endl;

        if (sr.getByLabel(_potsum_producer_data, _potsum_instance, potsum_h)) {
            if (_debug) std::cout &lt;&lt; "[Analysis::endSubRun] POT are valid" &lt;&lt; std::endl;
            _sr_pot = potsum_h-&gt;totpot;
        }
        else
            _sr_pot = 0;

    }

    _sr_tree-&gt;Fill();

    if (_debug) std::cout &lt;&lt; "[Analysis::endSubRun] Ends" &lt;&lt; std::endl;
}

DEFINE_ART_MODULE(GetPOT)

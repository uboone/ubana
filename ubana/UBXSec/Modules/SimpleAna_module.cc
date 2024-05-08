////////////////////////////////////////////////////////////////////////
// Class:       SimpleAna
// Plugin Type: analyzer (art v2_05_00)
// File:        SimpleAna_module.cc
//
// Generated at Tue Oct 31 16:36:23 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class SimpleAna
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module with simple example to retrieve results
 *
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version analyzer (art v2_05_00)
 *
 * \date 2017/03/10
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Tue Oct 31 16:36:23 2017
 *
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "ubobj/UBXSec/SelectionResult.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoVector.h"
#include "TTree.h"

class SimpleAna;


class SimpleAna : public art::EDAnalyzer {
public:
  explicit SimpleAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleAna(SimpleAna const &) = delete;
  SimpleAna(SimpleAna &&) = delete;
  SimpleAna & operator = (SimpleAna const &) = delete;
  SimpleAna & operator = (SimpleAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.
  TTree* _tree;
  int _run, _subrun, _event;
  bool _status_ccpi0, _status_ccincl;
  double _bnb_correction;
  double _nu_energy;
};


SimpleAna::SimpleAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree","");
  _tree->Branch("run",             &_run,              "run/I");
  _tree->Branch("subrun",          &_subrun,           "subrun/I");
  _tree->Branch("event",           &_event,            "event/I");
  _tree->Branch("bnb_correction",  &_bnb_correction,   "bnb_correction/D");
  _tree->Branch("status_ccincl",   &_status_ccincl,    "status_ccincl/O");
  _tree->Branch("status_ccpi0",    &_status_ccpi0,     "status_ccpi0/O");
  _tree->Branch("nu_energy",       &_nu_energy,        "nu_energy/D");
}

void SimpleAna::analyze(art::Event const & e)
{

  std::cout << "[SimpleAna] Simple UBXSec Validation Module. Starts." << std::endl;

  // Implementation of required member function here.
  auto selection_h = e.getHandle<std::vector<ubana::SelectionResult>>("UBXSec");
  if (!selection_h || selection_h->empty()) {
    std::cout << "[SimpleAna] SelectionResult handle is not valid or empty." << std::endl;
  }
  std::vector<ubana::SelectionResult> const& selection_v = *selection_h;

  if (selection_v[0].GetSelectionStatus()) {
    std::cout << "[SimpleAna] Event is selected" << std::endl;
  } else {
    std::cout << "[SimpleAna] Event is not selected" << std::endl;
    std::cout << "[SimpleAna] Failure reason " << selection_v[0].GetFailureReason()  << std::endl;
  }

  std::map<std::string,bool> failure_map = selection_v[0].GetCutFlowStatus();

  std::cout << "[SimpleAna] Now Printing Cut Flow Status" << std::endl;

  for (auto iter : failure_map) {
    std::cout << "[SimpleAna] Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
  }


  //
  // Trigger study
  //
  art::Handle<art::TriggerResults> trigger_h;
  e.getByLabel("TriggerResults::UBXSec", trigger_h);
  if (!trigger_h.isValid()) {
    std::cout << "[SimpleAna] TriggerResults handle is not valid." << std::endl;
  }
  art::TriggerResults trigger = *trigger_h;
  std::cout << "[SimpleAna] Trigger size " << trigger.size() << std::endl;
  if (trigger.size() != 2) {
    std::cout << "[SimpleAna] Trigger size is not 2. Returning now." << std::endl;
    return;
  }
  std::cout << "[SimpleAna] accept()? " << (trigger.accept() ? "YES" : "NO") << std::endl;
  std::cout << "[SimpleAna] accept(0)? " << (trigger.accept(0) ? "YES" : "NO") << std::endl;
  std::cout << "[SimpleAna] accept(1)? " << (trigger.accept(1) ? "YES" : "NO") << std::endl;

  std::cout << "[SimpleAna] Simple UBXSec Validation Module. Ends." << std::endl;


  //
  // MCTruth
  //

  bool in_tpcactive = false;
  bool energy_range = false;
  bool right_flavour = false;

  auto const& mclist = e.getProduct<std::vector<simb::MCTruth>>("generator");
  int iList = 0;
  ::geoalgo::Vector truth_nu_vtx (mclist[iList].GetNeutrino().Nu().Vx(),mclist[iList].GetNeutrino().Nu().Vy(),mclist[iList].GetNeutrino().Nu().Vz());
  auto const& tpc = art::ServiceHandle<geo::Geometry>{}->TPC();
  ::geoalgo::AABox tpcvol(0, (-1.)*(tpc.HalfHeight()), 0.,
              tpc.HalfWidth()*2, tpc.HalfHeight(), tpc.Length());

  if(tpcvol.Contain(truth_nu_vtx)) in_tpcactive = true;
  else in_tpcactive = false;

  if (mclist[iList].GetNeutrino().Nu().E() > 0.05 && mclist[iList].GetNeutrino().Nu().E() < 1.5)
    energy_range = true;
  else
    energy_range = false;

  if (mclist[iList].GetNeutrino().CCNC() == 0 && mclist[iList].GetNeutrino().Nu().PdgCode() == 12)
    right_flavour = true;
  else
    right_flavour = false;

  _nu_energy = mclist[iList].GetNeutrino().Nu().E();


  //
  // EventWeight
  //

  double bnb_weight = 1.;
  auto eventweight_h = e.getHandle<std::vector<evwgh::MCEventWeight>>("eventweight");
  if(!eventweight_h){
    std::cout << "[SimpleAna] MCEventWeight product not found..." << std::endl;
  }
  std::vector<evwgh::MCEventWeight> const& eventweight_v = *eventweight_h;

  if (eventweight_v.size() > 0) {
    evwgh::MCEventWeight const& evt_wgt = eventweight_v[0];
    for (auto entry : evt_wgt.fWeight) {
      if (entry.first.find("bnbcorrection") != std::string::npos) {
        bnb_weight *= entry.second.at(0);
        std::cout << "[SimpleAna] BNB Correction Weight: " << bnb_weight << std::endl;
      }
    }
  }

  //
  // Save to tree
  //

  if (right_flavour && in_tpcactive && energy_range) {
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
    _bnb_correction = bnb_weight;
    _status_ccincl = trigger.accept(0);
    _status_ccpi0 = trigger.accept(1);
    _tree->Fill();
  }

}

DEFINE_ART_MODULE(SimpleAna)

#ifndef HITCOSMICTAG_COSMICTAGMANAGER_CXX
#define HITCOSMICTAG_COSMICTAGMANAGER_CXX

#include <sstream>
#include <map>
#include <set>
#include "CosmicTagManager.h"

#include "HitCosmicTagException.h"
#include "StartHitFinderFactory.h"
#include "HitOrdererFactory.h"
#include "HitSmootherFactory.h"
#include "DqDsCalculatorFactory.h"
#include "DqDsSmootherFactory.h"
#include "LocalLinearityCalculatorFactory.h"
#include "CustomAlgoFactory.h"


namespace cosmictag {

  CosmicTagManager::CosmicTagManager(const std::string name)
    : _alg_start_hit_finder(nullptr)
    , _alg_hit_orderer(nullptr)
    , _alg_hit_smoother(nullptr)
    , _alg_dqds_calculator(nullptr)
    , _alg_dqds_smoother(nullptr)
    , _alg_linearity_calculator(nullptr)
    , _name(name)
  {
    _configured = false;

    if(_make_csv) _csvfile.open ("stopping_muon_tagger_helper.csv", std::ofstream::out | std::ofstream::trunc);
    if(_make_csv) _csvfile << "n,i,dqdx,dqdx_slider,linearity" << std::endl;
  }

  const std::string& CosmicTagManager::Name() const
  { return _name; }


  void CosmicTagManager::AddCustomAlgo(BaseAlgorithm* alg)
  {
    if(_custom_alg_m.find(alg->AlgorithmName()) != _custom_alg_m.end()) {
      std::stringstream ss;
      ss << "Duplicate name: " << alg->AlgorithmName() << std::endl;
      throw HitCosmicTagException(ss.str());
      }
    _custom_alg_m[alg->AlgorithmName()] = alg;
  }
  
  void CosmicTagManager::Configure(const Config_t& main_cfg) 
  {
   
    auto const& mgr_cfg = main_cfg.get<cosmictag::Config_t>(Name());

    this->set_verbosity((msg::Level_t)(mgr_cfg.get<unsigned int>("Verbosity")));

    auto const start_hit_finder_name      = mgr_cfg.get<std::string>("StartHitFinderAlgo","");
    auto const hit_orderer_name           = mgr_cfg.get<std::string>("HitOrdererAlgo","");
    auto const hit_smoother_name          = mgr_cfg.get<std::string>("HitSmootherAlgo","");
    auto const dqds_calculator_name       = mgr_cfg.get<std::string>("DqDsCalculatorAlgo","");
    auto const dqds_smoother_name         = mgr_cfg.get<std::string>("DqDsSmootherAlgo","");
    auto const linearity_calculator_name  = mgr_cfg.get<std::string>("LocalLinearityCalculatorAlgo","");
    
    std::vector<std::string> custom_algo_v;
    custom_algo_v = mgr_cfg.get<std::vector<std::string> >("CustomAlgo",custom_algo_v);

    if(!start_hit_finder_name.empty()) _alg_start_hit_finder = StartHitFinderFactory::get().create(start_hit_finder_name,start_hit_finder_name);
    if(!hit_orderer_name.empty()) _alg_hit_orderer = HitOrdererFactory::get().create(hit_orderer_name,hit_orderer_name);
    if(!hit_smoother_name.empty()) _alg_hit_smoother = HitSmootherFactory::get().create(hit_smoother_name,hit_smoother_name);
    if(!dqds_calculator_name.empty()) _alg_dqds_calculator = DqDsCalculatorFactory::get().create(dqds_calculator_name,dqds_calculator_name);
    if(!dqds_smoother_name.empty()) _alg_dqds_smoother = DqDsSmootherFactory::get().create(dqds_smoother_name,dqds_smoother_name);
    if(!linearity_calculator_name.empty()) _alg_linearity_calculator = LocalLinearityCalculatorFactory::get().create(linearity_calculator_name,linearity_calculator_name);

    for(auto const& name : custom_algo_v)
      if(!name.empty()) AddCustomAlgo(CustomAlgoFactory::get().create(name,name));


    _make_csv = mgr_cfg.get<bool>("MakeCSV", false);
    

    // Checks
  
    if (_alg_start_hit_finder) {

      _alg_start_hit_finder->Configure(main_cfg.get<cosmictag::Config_t>(_alg_start_hit_finder->AlgorithmName()));
      
    }
    if (_alg_hit_orderer) {

      _alg_hit_orderer->Configure(main_cfg.get<cosmictag::Config_t>(_alg_hit_orderer->AlgorithmName()));
      _alg_hit_orderer->CollectionCoplanar(_collection_coplanar);
    }
    if (_alg_hit_smoother) {

      _alg_hit_smoother->Configure(main_cfg.get<cosmictag::Config_t>(_alg_hit_smoother->AlgorithmName()));
      
    }
    if (_alg_dqds_calculator) {

      _alg_dqds_calculator->Configure(main_cfg.get<cosmictag::Config_t>(_alg_dqds_calculator->AlgorithmName()));
      
    }
    if (_alg_dqds_smoother) {

      _alg_dqds_smoother->Configure(main_cfg.get<cosmictag::Config_t>(_alg_dqds_smoother->AlgorithmName()));
      
    }
    if (_alg_linearity_calculator) {

      _alg_linearity_calculator->Configure(main_cfg.get<cosmictag::Config_t>(_alg_linearity_calculator->AlgorithmName()));
      
    }

     
    // Custom algo config
    for (auto& name_ptr : _custom_alg_m) {

      name_ptr.second->Configure(main_cfg.get<cosmictag::Config_t>(name_ptr.first));

    }

    _configured = true;
  }

  BaseAlgorithm* CosmicTagManager::GetAlgo(cosmictag::AlgoType type)
  {
    if (!_configured)
      CT_WARNING() << "Algorithm may be not configured yet!" << std::endl;

    // Figure out the type of a provided algorithm
    switch (type) {

    // Start Hit Finder
    case kStartHitFinder:
      return _alg_start_hit_finder;

    // Hit Orderer
    case kHitOrderer:
      return _alg_hit_orderer;

    case kHitSmoother:
      return _alg_hit_smoother;

    // dq/ds Calculator
    case kDqDsCalculator:
      return _alg_dqds_calculator;

    // dq/ds Smoother
    case kDqDsSmoother:
      return _alg_dqds_smoother;

    // Local Linearity 
    case kLocalLinearity:
      return _alg_linearity_calculator;

    // Fuck it
    default:
      std::stringstream ss;
      ss << "Unsupported algorithm type: " << type;
      throw HitCosmicTagException(ss.str());
    }
    return nullptr;
  }

  cosmictag::BaseAlgorithm* CosmicTagManager::GetCustomAlgo(std::string name)
  {
    if(_custom_alg_m.find(name) == _custom_alg_m.end()) {
      //CT_ERROR() << Form("Algorithm name %s not found!",name.c_str()) << std::endl;
      throw HitCosmicTagException();
    }
    return _custom_alg_m[name];
  }

  

  void CosmicTagManager::Emplace(cosmictag::SimpleCluster && obj)
  { _cluster = obj; }

  void CosmicTagManager::SetStartHit(cosmictag::SimpleHit && hit)
  { _start_hit = hit; }

  SimpleCluster CosmicTagManager::GetCluster()
  {
    return _cluster;
  }
  
  // CORE FUNCTION
  bool CosmicTagManager::Run()
  {     
    
    if (!_alg_start_hit_finder)
      throw HitCosmicTagException("Start hit finder algorithm is required! (not attached)");
    if (!_alg_hit_orderer)
      throw HitCosmicTagException("Hit orderer algorithm is required! (not attached)");
    if (!_alg_hit_smoother)
      throw HitCosmicTagException("Hit smoother algorithm is required! (not attached)");
    if (!_alg_dqds_calculator)
      throw HitCosmicTagException("dq/ds calculator algorithm is required! (not attached)");
    if (!_alg_dqds_smoother)
      throw HitCosmicTagException("dq/ds smoother algorithm is required! (not attached)");
    if (!_alg_linearity_calculator)
      throw HitCosmicTagException("local linearity algorithm is required! (not attached)");

    if (!_configured)
      throw HitCosmicTagException("Have not configured yet!");

    bool status = true;
   
    CT_DEBUG() << "Running start hit finder now." << std::endl;
    status = _alg_start_hit_finder->FindStartHit(_cluster, _start_hit);
    if (!status) return false;

    CT_DEBUG() << "Running hit orderer algo now." << std::endl;
    _alg_hit_orderer->CollectionCoplanar(_collection_coplanar);
    status = _alg_hit_orderer->OrderHits(_cluster);
    if (!status) return false;

    CT_DEBUG() << "Running hit smoother algo now." << std::endl;
    status = _alg_hit_smoother->Smooth(_cluster);
    if (!status) return false;

    CT_DEBUG() << "Running dq/ds calculator now." << std::endl;
    status = _alg_dqds_calculator->CalculateDqDs(_cluster);
    if (!status) return false;

    CT_DEBUG() << "Running dq/ds smoother now." << std::endl;
    status = _alg_dqds_smoother->SmoothDqDs(_cluster);
    if (!status) return false;

    CT_DEBUG() << "Running linearity calculator now." << std::endl;
    status = _alg_linearity_calculator->CalculateLocalLinearity(_cluster);
    if (!status) return false;

    _ready = true;

    return true;

  }



  bool CosmicTagManager::MakeDecision(std::string name) 
  {

    if (!_configured)
      throw HitCosmicTagException("Have not configured yet!");

    if (!_ready)
      throw HitCosmicTagException("Have not run yet!");


    if(_custom_alg_m.find(name) == _custom_alg_m.end()) {
      //CT_ERROR() << Form("Algorithm name %s not found!",name.c_str()) << std::endl;
      throw HitCosmicTagException();
    }
    //return _custom_alg_m[name];

    return false;
  }



  void CosmicTagManager::PrintConfig() {
    
    std::cout << "---- COSMIC TAG MANAGER PRINTING CONFIG     ----" << std::endl
	      << "_name = " << _name << std::endl 
	      << "_alg_start_hit_finder?" << std::endl;
    if (_alg_start_hit_finder)
      std::cout << "\t" << _alg_start_hit_finder->AlgorithmName() << std::endl;
    std::cout << "_alg_hit_orderer?" << std::endl;
    if (_alg_hit_orderer)
      std::cout << "\t" << _alg_hit_orderer->AlgorithmName() << std::endl;
    std::cout << "_alg_hit_smoother?" << std::endl;
    if (_alg_hit_smoother)
      std::cout << "\t" << _alg_hit_smoother->AlgorithmName() << std::endl;
    std::cout << "_alg_dqds_calculator?" << std::endl;
    if (_alg_dqds_calculator)
      std::cout << "\t" << _alg_dqds_calculator->AlgorithmName() << std::endl;
    std::cout << "_alg_dqds_smoother?" << std::endl;
    if (_alg_dqds_smoother)
      std::cout << "\t" << _alg_dqds_smoother->AlgorithmName() << std::endl;
    std::cout << "_alg_linearity_calculator?" << std::endl;
    if (_alg_linearity_calculator)
      std::cout << "\t" << _alg_linearity_calculator->AlgorithmName() << std::endl;
    std::cout << "_custom_alg_m?" << std::endl;
    for (auto& name_ptr : _custom_alg_m)
      std::cout << "\t" << name_ptr.first << std::endl;
    std::cout << "---- END FLASH MATCH MANAGER PRINTING CONFIG ----" << std::endl;
    
  }

  void CosmicTagManager::PrintClusterStatus() {

    auto & s_hit_v = _cluster._s_hit_v;
    auto & _dqds_slider = _cluster._dqds_slider;
    auto & _linearity_v = _cluster._linearity_v;

    CT_NORMAL() << "Current Cluster Status:" << std::endl;

    int counter = 0;
    for (auto h : s_hit_v) {
      CT_NORMAL() << "index " << counter
                  << ", wire: " << h.wire
                  << ", time: " << h.time*4
                  << ", dqdx_slider: " << _dqds_slider.at(counter)
                  << ", linearity: " << _linearity_v.at(counter) << std::endl;
      counter++;
    }
  }

  void CosmicTagManager::PrintOnFile(int index) {

    if (_cluster._s_hit_v.size() != _cluster._linearity_v.size() 
      || _cluster._s_hit_v.size() != _cluster._dqds_slider.size()) {
      CT_CRITICAL() << "Vectors in cluster have different size!" << std::endl;
      throw HitCosmicTagException();
    }

    if (!_make_csv) {
      CT_CRITICAL() << "_make_csv has to be true to print on file." << std::endl;
      throw HitCosmicTagException();
    }

    for (size_t i = 0; i < _cluster._dqds_slider.size(); i++) {
      _csvfile << index << "," 
               << i << "," 
               << _cluster._dqds_v.at(i) << "," 
               << _cluster._dqds_slider.at(i) << ", "
               << _cluster._linearity_v.at(i)
               << std::endl;
    }
  }


}

#endif

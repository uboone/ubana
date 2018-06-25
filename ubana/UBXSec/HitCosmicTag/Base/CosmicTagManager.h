/**
 * \file CosmicTagManager.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class CosmicTagManager
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Base

    @{*/
#ifndef HITCOSMICTAG_COSMICTAGMANAGER_H
#define HITCOSMICTAG_COSMICTAGMANAGER_H

#include "LoggerFeature.h"
#include "HitCosmicTagFMWKInterface.h"
#include "BaseAlgorithm.h"
#include "BaseStartHitFinderAlgo.h"
#include "BaseHitOrdererAlgo.h"
#include "BaseHitSmootherAlgo.h"
#include "BaseDqDsCalculatorAlgo.h"
#include "BaseDqDsSmootherAlgo.h"
#include "BaseLocalLinearityCalculatorAlgo.h"

#include <fstream>

namespace cosmictag {
  /**
     \class CosmicTagManager
  */
  class CosmicTagManager : public LoggerFeature {

  public:
    
    /// Default constructor
    CosmicTagManager(const std::string name="CosmicTagManager");
    
    /// Default destructor
    ~CosmicTagManager(){}

    /// Name getter
    const std::string& Name() const;

    /// Configuration
    void Configure(const Config_t& cfg);

    /// Algorithm getter
    cosmictag::BaseAlgorithm* GetAlgo(cosmictag::AlgoType type);

    /// Custom algorithm getter
    cosmictag::BaseAlgorithm* GetCustomAlgo(std::string name);
		 
    /// Emplacer of a Cluster
    void Emplace(cosmictag::SimpleCluster && obj);

    void SetStartHit(cosmictag::SimpleHit && hit);

    SimpleCluster GetCluster();

    /**
       CORE FUNCTION: executes algorithms to find a match of TPC object and flash provided by users. \n
       The execution takes following steps:             \n
       0) TPC filter algorithm if provided (optional)   \n
       1) Flash filter algorithm if provided (optional) \n
       3) Flash matching algorithm (required)           \n
       4) Returns match information for created TPC object & flash pair which respects the outcome of 3)
     */
    bool Run();

    bool MakeDecision(std::string);

    /// Clears locally kept TPC object (QClusterArray_t) and flash (FlashArray_t), both provided by a user
    void Reset()
    { /*_tpc_object_v.clear(); _flash_v.clear();*/ _ready = false; _collection_coplanar = false;}

    ///
    void CollectionCoplanar(bool status)
    { _collection_coplanar = status; }

    /// Configuration option: true => allows an assignment of the same flash to multiple TPC objects
    void CanReuseFlash(bool ok=true)
    { /*_allow_reuse_flash = ok;*/ }

    void PrintConfig();

    void PrintClusterStatus();

    void PrintOnFile(int index);


  private:

    void AddCustomAlgo(BaseAlgorithm* alg);

    BaseStartHitFinderAlgo*           _alg_start_hit_finder;          ///< Find start hit in cluster Algorithm
    BaseHitOrdererAlgo*               _alg_hit_orderer;               ///< Order Hits Algorithm
    BaseHitSmootherAlgo*              _alg_hit_smoother;              ///< Smooth Hits Algorithm
    BaseDqDsCalculatorAlgo*           _alg_dqds_calculator;           ///< Calculate dq/ds algorithm
    BaseDqDsSmootherAlgo*             _alg_dqds_smoother;             ///< Smooths dqds algorithm
    BaseLocalLinearityCalculatorAlgo* _alg_linearity_calculator;      ///< Calculate linearity algorithm
   
    /**
       A set of custom algorithms (not to be executed but to be configured)
    */
    std::map<std::string,cosmictag::BaseAlgorithm*> _custom_alg_m;

    /// The cluster to be analysed
    SimpleCluster _cluster;

    /// The start hit of the cluster
    SimpleHit _start_hit;

    /// Collection coplanar flas
    bool _collection_coplanar = false;

    /// Configuration readiness flag
    bool _configured = false;

    /// Readiness flag
    bool _ready = false;

    /// If true makes a CSV file
    bool _make_csv = false;

    /// Configuration file
    std::string _config_file;
    /// Name
    std::string _name;
    /// Request boolean to store full matching result (per Match function call)
    //bool _store_full;   

    std::ofstream _csvfile; 
  };
}

#endif
/** @} */ // end of doxygen group 


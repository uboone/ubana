#ifndef ANALYSIS_BDT_CXX
#define ANALYSIS_BDT_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "ubana/XGBoost/xgboost/c_api.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BDT
    // File:        BDT.cc
    //
    //              A basic analysis example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
    //
    ////////////////////////////////////////////////////////////////////////

  class BDT : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    BDT(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~BDT();
    
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
    
  private:

    art::InputTag fTrigResProducer;
    BoosterHandle booster_nuNCpi0;
    BoosterHandle booster_numuCCpi0;
    BoosterHandle booster_numuCC;
    BoosterHandle booster_ext;
    BoosterHandle booster_cosmic;
    BoosterHandle booster_global;
    float _bdt_nuNCpi0;
    float _bdt_numuCCpi0;
    float _bdt_numuCC;
    float _bdt_ext;
    float _bdt_cosmic;
    float _bdt_global;
    int _pass_antibdt_filter;
    BoosterHandle booster_pi0_np;
    BoosterHandle booster_nonpi0_np;
    BoosterHandle booster_bkg_0p;
    float _bdt_pi0_np;
    float _bdt_nonpi0_np;
    float _bdt_bkg_0p;
    float _anglediff_Y;
    float _anglediff_V;
    float _anglediff_U;
    float _trkpid;
    TTree* _mytree;
    bool fVerbose;
    bool fOnlyLegacy;
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  BDT::BDT(const fhicl::ParameterSet& p)
  {
    fTrigResProducer = p.get<art::InputTag>("TrigResProducer");
    fVerbose = p.get<bool>("Verbose", false);
    fOnlyLegacy = p.get<bool>("OnlyLegacy", false);
    //documentation in xgboost/include/xgboost/c_api.h
    int xgtest = -1;
    xgtest = XGBoosterCreate(NULL, 0, &booster_nuNCpi0);
    xgtest = XGBoosterCreate(NULL, 0, &booster_numuCCpi0);
    xgtest = XGBoosterCreate(NULL, 0, &booster_numuCC);
    xgtest = XGBoosterCreate(NULL, 0, &booster_ext);
    xgtest = XGBoosterCreate(NULL, 0, &booster_cosmic);
    xgtest = XGBoosterCreate(NULL, 0, &booster_global);
    xgtest = XGBoosterCreate(NULL, 0, &booster_pi0_np);
    xgtest = XGBoosterCreate(NULL, 0, &booster_nonpi0_np);
    xgtest = XGBoosterCreate(NULL, 0, &booster_bkg_0p);

    std::string _filename;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file("searchingfornues/booster_nopid_ncpi0.model",_filename);
    xgtest = XGBoosterLoadModel(booster_nuNCpi0, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_ccpi0.model",_filename);
    xgtest = XGBoosterLoadModel(booster_numuCCpi0, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_cc.model",_filename);
    xgtest = XGBoosterLoadModel(booster_numuCC, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_ext.model",_filename);
    xgtest = XGBoosterLoadModel(booster_ext, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_cosmic.model",_filename);
    xgtest = XGBoosterLoadModel(booster_cosmic, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid.model",_filename);
    xgtest = XGBoosterLoadModel(booster_global, _filename.c_str());
    sp.find_file("searchingfornues/booster_pi0_0304_extnumi.model",_filename);
    xgtest = XGBoosterLoadModel(booster_pi0_np, _filename.c_str());
    sp.find_file("searchingfornues/booster_nonpi0_0304_extnumi.model",_filename);
    xgtest = XGBoosterLoadModel(booster_nonpi0_np, _filename.c_str());
    sp.find_file("searchingfornues/booster_bkg_0304_noext.model",_filename);
    xgtest = XGBoosterLoadModel(booster_bkg_0p, _filename.c_str());
    assert(xgtest==0);

    // Need to do something with xgtest in prof build.

    if(xgtest != 0)
      throw cet::exception("BDT_tool") << "xgtest = " << xgtest;
  }

  BDT::~BDT()
  {
    int xgtest = -1;
    xgtest = XGBoosterFree(booster_nuNCpi0);
    xgtest = XGBoosterFree(booster_numuCCpi0);
    xgtest = XGBoosterFree(booster_numuCC);
    xgtest = XGBoosterFree(booster_ext);
    xgtest = XGBoosterFree(booster_cosmic);
    xgtest = XGBoosterFree(booster_global);
    xgtest = XGBoosterFree(booster_pi0_np);
    xgtest = XGBoosterFree(booster_nonpi0_np);
    xgtest = XGBoosterFree(booster_bkg_0p);
    assert(xgtest==0);

    // Need to do something with xgtest in prof build.

    if(xgtest != 0)
      std::cout << "In BDT destructor xgtest = " << xgtest << std::endl;   // Can't throw here.
  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void BDT::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void BDT::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {
    art::Handle<art::TriggerResults> trigRes;
    if (fTrigResProducer!="") {
      e.getByLabel(fTrigResProducer,trigRes);
      fhicl::ParameterSet pset;
      if (!fhicl::ParameterSetRegistry::get(trigRes->parameterSetID(), pset)) { throw cet::exception("PSet Not Found???"); }
      std::vector<std::string> trigger_path_names = pset.get<std::vector<std::string> >("trigger_paths", {});
      if (trigger_path_names.size()!=trigRes->size()) { throw cet::exception("Size mismatch???"); }
      for (size_t itp=0;itp<trigRes->size();itp++) {
	if (fVerbose) std::cout << "filter name=" << trigger_path_names.at(itp) << " decision=" << trigRes->at(itp).accept() << std::endl;
	if (trigger_path_names.at(itp)=="antibdt") _pass_antibdt_filter = trigRes->at(itp).accept();
      }
    }

    std::vector<float> data;
    std::vector<std::string> variables{"shr_dedx_Y", "shr_distance", "trk_distance", "pt", "hits_y",
	"shr_tkfit_dedx_Y", "shr_tkfit_dedx_U", "shr_tkfit_dedx_V", "p",
	"hits_ratio", "shr_dedx_U", "shr_dedx_V", "n_tracks_contained", "n_showers_contained",
	"shr_theta", "trk_len", "trk_score", "shr_score", "shr_energy_tot_cali", "trk_energy_tot",
	"shr_phi", "trk_theta", "trk_phi", "tksh_angle", "tksh_distance", "CosmicIP",
	"shr_pca_2", "shr_pca_1", "shr_pca_0",
	"topological_score", "slpdg"};
    std::vector<std::string> type{"F", "F", "F", "F", "i",
	"F", "F", "F", "F",
	"F", "F", "F", "I", "I",
	"F", "F", "F", "F", "F", "F",
	"F", "F", "F", "F", "F", "F",
	"F", "F", "F",
	"F", "I"};
    //retrieve variables from the tree abd fill data vector accordingly
    for (size_t iv=0;iv<variables.size();++iv) {
      if (type[iv]=="F") {
	float* tmp = (float*) _mytree->GetBranch(variables[iv].c_str())->GetAddress();
	data.push_back(*tmp);
	if (fVerbose) std::cout << type[iv] << " " << variables[iv] << "=" << data.back() << " - " << *tmp << std::endl;
      } else  if (type[iv]=="I") {
	int* tmp = (int*) _mytree->GetBranch(variables[iv].c_str())->GetAddress();
	data.push_back(*tmp);
	if (fVerbose) std::cout << type[iv] << " " << variables[iv] << "=" << data.back() << " - " << *tmp << std::endl;
      } else  if (type[iv]=="i") {
	unsigned int* tmp = (unsigned int*) _mytree->GetBranch(variables[iv].c_str())->GetAddress();
	data.push_back(*tmp);
	if (fVerbose) std::cout << type[iv] << " " << variables[iv] << "=" << data.back() << " - " << *tmp << std::endl;
      } else {exit(1);}
    }
    //get predictions: 5 BDTs against specific backgrounds
    int xgtest = -1;
    DMatrixHandle mat_nuNCpi0;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_nuNCpi0);
    bst_ulong out_len_nuNCpi0 = 0;
    const float* out_result_nuNCpi0   = NULL;
    xgtest = XGBoosterPredict(booster_nuNCpi0  , mat_nuNCpi0, 0, 317, &out_len_nuNCpi0, &out_result_nuNCpi0  );// NB: setting by hand ntree_limit=booster.best_iteration
    _bdt_nuNCpi0   = *out_result_nuNCpi0  ;
    if (fVerbose) std::cout << "bdt_nuNCpi0=" << *out_result_nuNCpi0 << " " << _bdt_nuNCpi0 << std::endl;
    //
    DMatrixHandle mat_numuCCpi0;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_numuCCpi0);
    bst_ulong out_len_numuCCpi0 = 0;
    const float* out_result_numuCCpi0 = NULL;
    xgtest = XGBoosterPredict(booster_numuCCpi0, mat_numuCCpi0, 0,  343, &out_len_numuCCpi0, &out_result_numuCCpi0);
    _bdt_numuCCpi0 = *out_result_numuCCpi0;
    if (fVerbose) std::cout << "bdt_numuCCpi0=" << *out_result_numuCCpi0 << " " << _bdt_numuCCpi0 << std::endl;
    //
    DMatrixHandle mat_numuCC;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_numuCC);
    bst_ulong out_len_numuCC = 0;
    const float* out_result_numuCC    = NULL;
    xgtest = XGBoosterPredict(booster_numuCC   , mat_numuCC, 0, 491, &out_len_numuCC, &out_result_numuCC   );
    _bdt_numuCC    = *out_result_numuCC   ;
    if (fVerbose) std::cout << "bdt_numuCC=" << *out_result_numuCC << " " << _bdt_numuCC << std::endl;
    //
    DMatrixHandle mat_ext;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_ext);
    bst_ulong out_len_ext = 0;
    const float* out_result_ext       = NULL;
    xgtest = XGBoosterPredict(booster_ext      , mat_ext, 0, 265, &out_len_ext, &out_result_ext      );
    _bdt_ext       = *out_result_ext      ;
    if (fVerbose) std::cout << "bdt_ext=" << *out_result_ext << " " << _bdt_ext << std::endl;
    //
    DMatrixHandle mat_cosmic;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_cosmic);
    bst_ulong out_len_cosmic = 0;
    const float* out_result_cosmic    = NULL;
    xgtest = XGBoosterPredict(booster_cosmic   , mat_cosmic, 0, 304, &out_len_cosmic, &out_result_cosmic   );
    _bdt_cosmic    = *out_result_cosmic   ;
    if (fVerbose) std::cout << "bdt_cosmic=" << *out_result_cosmic << " " << _bdt_cosmic << std::endl;
    //
    // const char** dump_result_cosmic       = NULL;
    // xgtest = XGBoosterDumpModel(booster_cosmic,"",0, &out_len_cosmic, &dump_result_cosmic);
    // for (bst_ulong k=0;k<out_len_cosmic;k++){
    //   printf("%s \n",dump_result_cosmic[k]);
    // }
    //

    //get global BDT prediction
    std::vector<float> data_global;
    data_global.push_back(_bdt_ext      );
    data_global.push_back(_bdt_nuNCpi0  );
    data_global.push_back(_bdt_numuCC   );
    data_global.push_back(_bdt_numuCCpi0);
    data_global.push_back(_bdt_cosmic   );
    DMatrixHandle mat_global;
    xgtest = XGDMatrixCreateFromMat(data_global.data(),1,data_global.size(),0,&mat_global);
    bst_ulong out_len_global = 0;
    const float* out_result_global    = NULL;
    xgtest = XGBoosterPredict(booster_global, mat_global, 0, 115, &out_len_global, &out_result_global   );
    _bdt_global    = *out_result_global   ;
    if (fVerbose) std::cout << "bdt_global=" << *out_result_global << " " << _bdt_global << std::endl;

    if (fOnlyLegacy) return;

    //new BDTs
    float* secondshower_Y_dir = (float*) _mytree->GetBranch("secondshower_Y_dir")->GetAddress();
    float* secondshower_V_dir = (float*) _mytree->GetBranch("secondshower_V_dir")->GetAddress();
    float* secondshower_U_dir = (float*) _mytree->GetBranch("secondshower_U_dir")->GetAddress();
    float* shrclusdir2 = (float*) _mytree->GetBranch("shrclusdir2")->GetAddress();
    float* shrclusdir1 = (float*) _mytree->GetBranch("shrclusdir1")->GetAddress();
    float* shrclusdir0 = (float*) _mytree->GetBranch("shrclusdir0")->GetAddress();
    _anglediff_Y = std::abs( (*secondshower_Y_dir)-(*shrclusdir2) );
    _anglediff_V = std::abs( (*secondshower_V_dir)-(*shrclusdir1) );
    _anglediff_U = std::abs( (*secondshower_U_dir)-(*shrclusdir0) );

    std::vector<float>** trk_llr_pid_score_v = (std::vector<float>**) _mytree->GetBranch("trk_llr_pid_score_v")->GetAddress();
    unsigned int* trk_id = (unsigned int*) _mytree->GetBranch("trk_id")->GetAddress();
    unsigned int trk_ilen = (**trk_llr_pid_score_v).size();
    if ( ((*trk_id)-1) < trk_ilen ) _trkpid = (**trk_llr_pid_score_v)[(*trk_id)-1];
    if (fVerbose) std::cout << "trk_id=" << (*trk_id) << " trk_ilen=" << trk_ilen << " trkpid=" <<  _trkpid << std::endl;

    std::vector<float> data_np;
    std::vector<std::string> variables_np{"shr_score","tksh_distance","tksh_angle","shr_tkfit_dedx_max","trkfit","trkpid","subcluster","shrmoliereavg",
	"trkshrhitdist2","hits_ratio","secondshower_Y_nhit","secondshower_Y_vtxdist","secondshower_Y_dot","anglediff_Y","CosmicIPAll3D","CosmicDirAll3D"};
    std::vector<std::string> type_np{"F","F","F","F","F","F","i","f",
	"F","F","I","F","F","F","F","F"};
    //retrieve variables from the tree and fill data vector accordingly
    for (size_t iv=0;iv<variables_np.size();++iv) {
      if (type_np[iv]=="F") {
	float* tmp = (float*) _mytree->GetBranch(variables_np[iv].c_str())->GetAddress();
	data_np.push_back(*tmp);
	if (fVerbose) std::cout << type_np[iv] << " " << variables_np[iv] << "=" << data_np.back() << " - " << *tmp << std::endl;
      } else  if (type_np[iv]=="f") {
	Float16_t* tmp = (Float16_t*) _mytree->GetBranch(variables_np[iv].c_str())->GetAddress();
	data_np.push_back(*tmp);
	if (fVerbose) std::cout << type_np[iv] << " " << variables_np[iv] << "=" << data_np.back() << " - " << *tmp << std::endl;
      } else  if (type_np[iv]=="I") {
	int* tmp = (int*) _mytree->GetBranch(variables_np[iv].c_str())->GetAddress();
	data_np.push_back(*tmp);
	if (fVerbose) std::cout << type_np[iv] << " " << variables_np[iv] << "=" << data_np.back() << " - " << *tmp << std::endl;
      } else  if (type_np[iv]=="i") {
	unsigned int* tmp = (unsigned int*) _mytree->GetBranch(variables_np[iv].c_str())->GetAddress();
	data_np.push_back(*tmp);
	if (fVerbose) std::cout << type_np[iv] << " " << variables_np[iv] << "=" << data_np.back() << " - " << *tmp << std::endl;
      } else {exit(1);}
    }
    //get predictions: pi0 ad nonpi0 BDTs for Np channel
    DMatrixHandle mat_pi0_np;
    xgtest = XGDMatrixCreateFromMat(data_np.data(),1,data_np.size(),0,&mat_pi0_np);
    bst_ulong out_len_pi0_np = 0;
    const float* out_result_pi0_np   = NULL;
    xgtest = XGBoosterPredict(booster_pi0_np, mat_pi0_np, 0, 996, &out_len_pi0_np, &out_result_pi0_np  );// NB: setting by hand ntree_limit=booster.best_iteration
    _bdt_pi0_np = *out_result_pi0_np;
    if (fVerbose) std::cout << "bdt_pi0_np=" << *out_result_pi0_np << " " << _bdt_pi0_np << std::endl;
    //
    DMatrixHandle mat_nonpi0_np;
    xgtest = XGDMatrixCreateFromMat(data_np.data(),1,data_np.size(),0,&mat_nonpi0_np);
    bst_ulong out_len_nonpi0_np = 0;
    const float* out_result_nonpi0_np   = NULL;
    xgtest = XGBoosterPredict(booster_nonpi0_np, mat_nonpi0_np, 0, 999, &out_len_nonpi0_np, &out_result_nonpi0_np  );// NB: setting by hand ntree_limit=booster.best_iteration
    _bdt_nonpi0_np   = *out_result_nonpi0_np  ;
    if (fVerbose) std::cout << "bdt_nonpi0_np=" << *out_result_nonpi0_np << " " << _bdt_nonpi0_np << std::endl;
    //

    if (fVerbose) {
      float* reco_e = (float*) _mytree->GetBranch("reco_e")->GetAddress();
      std::cout << "reco_e=" << (*reco_e) << std::endl;
    }

    std::vector<float> data_0p;
    std::vector<std::string> variables_0p{"shrmoliereavg","shr_score", "trkfit","subcluster","CosmicIPAll3D","CosmicDirAll3D",
            "secondshower_Y_nhit","secondshower_Y_vtxdist","secondshower_Y_dot","anglediff_Y",
            "secondshower_V_nhit","secondshower_V_vtxdist","secondshower_V_dot","anglediff_V",
            "secondshower_U_nhit","secondshower_U_vtxdist","secondshower_U_dot","anglediff_U",
            "shr_tkfit_2cm_dedx_U", "shr_tkfit_2cm_dedx_V", "shr_tkfit_2cm_dedx_Y",
            "shr_tkfit_gap10_dedx_U", "shr_tkfit_gap10_dedx_V", "shr_tkfit_gap10_dedx_Y",
            "shrMCSMom", "DeltaRMS2h", "shrPCA1CMed_5cm", "CylFrac2h_1cm"};
    std::vector<std::string> type_0p{"f","F","F","i","F","F",
	"I","F","F","F",
	"I","F","F","F",
	"I","F","F","F",
	"F","F","F",
	"F","F","F",
	"f","f","f","f"};
    //retrieve variables from the tree and fill data vector accordingly
    for (size_t iv=0;iv<variables_0p.size();++iv) {
      if (type_0p[iv]=="F") {
	float* tmp = (float*) _mytree->GetBranch(variables_0p[iv].c_str())->GetAddress();
	data_0p.push_back(*tmp);
	if (fVerbose) std::cout << type_0p[iv] << " " << variables_0p[iv] << "=" << data_0p.back() << " - " << *tmp << std::endl;
      } else  if (type_0p[iv]=="f") {
	Float16_t* tmp = (Float16_t*) _mytree->GetBranch(variables_0p[iv].c_str())->GetAddress();
	data_0p.push_back(*tmp);
	if (fVerbose) std::cout << type_0p[iv] << " " << variables_0p[iv] << "=" << data_0p.back() << " - " << *tmp << std::endl;
      } else  if (type_0p[iv]=="I") {
	int* tmp = (int*) _mytree->GetBranch(variables_0p[iv].c_str())->GetAddress();
	data_0p.push_back(*tmp);
	if (fVerbose) std::cout << type_0p[iv] << " " << variables_0p[iv] << "=" << data_0p.back() << " - " << *tmp << std::endl;
      } else  if (type_0p[iv]=="i") {
	unsigned int* tmp = (unsigned int*) _mytree->GetBranch(variables_0p[iv].c_str())->GetAddress();
	data_0p.push_back(*tmp);
	if (fVerbose) std::cout << type_0p[iv] << " " << variables_0p[iv] << "=" << data_0p.back() << " - " << *tmp << std::endl;
      } else {exit(1);}
    }
    //get predictions: bkg BDT for 0p channel
    DMatrixHandle mat_bkg_0p;
    xgtest = XGDMatrixCreateFromMat(data_0p.data(),1,data_0p.size(),0,&mat_bkg_0p);
    bst_ulong out_len_bkg_0p = 0;
    const float* out_result_bkg_0p   = NULL;
    xgtest = XGBoosterPredict(booster_bkg_0p, mat_bkg_0p, 0, 741, &out_len_bkg_0p, &out_result_bkg_0p  );// NB: setting by hand ntree_limit=booster.best_iteration
    _bdt_bkg_0p = *out_result_bkg_0p;
    if (fVerbose) std::cout << "bdt_bkg_0p=" << *out_result_bkg_0p << " " << _bdt_bkg_0p << std::endl;
    //

    assert(xgtest==0);

    // Need to do something with xgtest in prof build.

    if(xgtest != 0)
      throw cet::exception("BDT_tool") << "xgtest = " << xgtest;

    return;
  }

  void BDT::analyzeEvent(art::Event const &e, bool fData)
  {
    // std::cout << "analyze event" << std::endl;
  }

  void BDT::setBranches(TTree* _tree)
  {
    _tree->Branch("bdt_nuNCpi0"  ,&_bdt_nuNCpi0  ,"bdt_nuNCpi0/F"  );
    _tree->Branch("bdt_numuCCpi0",&_bdt_numuCCpi0,"bdt_numuCCpi0/F");
    _tree->Branch("bdt_numuCC"   ,&_bdt_numuCC   ,"bdt_numuCC/F"   );
    _tree->Branch("bdt_ext"      ,&_bdt_ext      ,"bdt_ext/F"      );
    _tree->Branch("bdt_cosmic"   ,&_bdt_cosmic   ,"bdt_cosmic/F"   );
    _tree->Branch("bdt_global"   ,&_bdt_global   ,"bdt_global/F"   );
    _tree->Branch("pass_antibdt_filter", &_pass_antibdt_filter,"bdt_global/I");
    _tree->Branch("bdt_pi0_np"   ,&_bdt_pi0_np   ,"bdt_pi0_np/F"   );
    _tree->Branch("bdt_nonpi0_np",&_bdt_nonpi0_np,"bdt_nonpi0_np/F");
    _tree->Branch("bdt_bkg_0p"   ,&_bdt_bkg_0p   ,"bdt_bkg_0p/F"   );
    _tree->Branch("anglediff_Y"   ,&_anglediff_Y   ,"anglediff_Y/F");
    _tree->Branch("anglediff_V"   ,&_anglediff_V   ,"anglediff_V/F");
    _tree->Branch("anglediff_U"   ,&_anglediff_U   ,"anglediff_U/F");
    _tree->Branch("trkpid"        ,&_trkpid        ,"trkpid/F"     );
    _mytree = _tree;//not ideal, be careful...
  }
  
  void BDT::resetTTree(TTree* _tree)
  {
    _bdt_nuNCpi0   = 9999.;
    _bdt_numuCCpi0 = 9999.;
    _bdt_numuCC    = 9999.;
    _bdt_ext       = 9999.;
    _bdt_cosmic    = 9999.;
    _bdt_global    = 9999.;
    _pass_antibdt_filter = -9999;
    _bdt_pi0_np    = -9999.;
    _bdt_nonpi0_np = -9999.;
    _bdt_bkg_0p    = -9999.;
    _anglediff_Y   = -9999.;
    _anglediff_V   = -9999.;
    _anglediff_U   = -9999.;
    _trkpid        = 9999.;
  }

  
  DEFINE_ART_CLASS_TOOL(BDT)
} // namespace analysis

#endif

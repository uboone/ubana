#pragma once
#ifndef ANALYSIS_NEUTRONWEIGHTS_CXX
#define ANALYSIS_NEUTRONWEIGHTS_CXX

#include "AnalysisToolBase.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "canvas/Persistency/Common/TriggerResults.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <map>
#include <random>
#include <numeric>
#include <functional>
#include "TRandom3.h"

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <cassert>
#include <algorithm>

namespace analysis {

class NeutronWeights : public AnalysisToolBase {
public:
  explicit NeutronWeights(fhicl::ParameterSet const& pset);
  ~NeutronWeights() override = default;

  void configure(fhicl::ParameterSet const& pset);
  void setBranches(TTree* t) override;   // sets the tree
  void resetTTree(TTree* t) override;
//  void analyzeEvent(art::Event const& e, bool fData) override;
//  void analyzeSlice(art::Event const&, std::vector<ProxyPfpElem_t>&, bool, bool) override;// override {}

  void analyzeEvent(art::Event const &e, bool fData) override;

  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

private:
  // ---- helpers ----
  static inline bool isInsideTPC(float x, float y, float z) {
    return (x > 0.f && x < 256.35f &&
            y > -116.5f && y < 116.5f &&
            z > 0.f && z < 1036.8f);
  }
  double inDetectorSegmentLength(float sx, float sy, float sz,
                                 float ex, float ey, float ez) const;

template <class T>
T* getBranchAddr(std::string const& name, char const* who) const {
  auto* br = _tree ? _tree->GetBranch(name.c_str()) : nullptr;
  if (!br) {
    throw cet::exception(who) << "Missing branch '" << name << "'";
  }
  void* addr = br->GetAddress();
  if (!addr) {
    throw cet::exception(who) << "Null address for branch '" << name
                              << "'. Was SetBranchAddress called upstream?";
  }
     return *reinterpret_cast<T**>(addr);
     }

/*  template <class T>
  T* getBranchAddr(std::string const& name, char const* who) const {
    auto* br = _tree ? _tree->GetBranch(name.c_str()) : nullptr;
    if (!br) throw cet::exception(who) << "Missing branch '" << name << "'";
    auto* addr = br->GetAddress();
    if (!addr) throw cet::exception(who) << "Null address for '" << name << "'";
    return reinterpret_cast<T*>(addr);
  }
*/
  // ---- config ----
  std::string fXSecFile;
  int         fNUniverses{1000};
  int         fSeedBase{1337};

  // input branch names
  std::string b_all_mc_pdg         {"all_mc_pdg"};
  std::string b_all_mc_vx          {"all_mc_vx"};
  std::string b_all_mc_vy          {"all_mc_vy"};
  std::string b_all_mc_vz          {"all_mc_vz"};
  std::string b_all_mc_endx        {"all_mc_endx"};
  std::string b_all_mc_endy        {"all_mc_endy"};
  std::string b_all_mc_endz        {"all_mc_endz"};
  std::string b_all_mc_E           {"all_mc_E"};
  std::string b_all_mc_end_process {"all_mc_end_process"};

  // ---- tree ----
  TTree* _tree {nullptr};

  // ---- outputs ----
  std::vector<float> _w_neutron_reint;           // event multisim weights (size=NUniverses)
  std::vector<float> _weight_neutron_argon_xsec; // event CV (size=1)
  std::vector<float> _neutron_track_weights;     // per-neutron CV
  std::vector<float> _neutron_track_len;         // per-neutron [cm]
  std::vector<float> _neutron_ke;                // per-neutron [GeV]
  std::vector<float> _neutron_last_w;            // per-neutron, universe 0 (debug)
};

} // namespace analysis

#endif

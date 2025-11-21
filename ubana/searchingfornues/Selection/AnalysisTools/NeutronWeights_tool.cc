#include "NeutronWeights_tool.hh"

// Reweighter
#include "NeutronReweighter.h"
#include "TRandom3.h"
#include "TBranch.h"
#include <numeric>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>

namespace analysis {

NeutronWeights::NeutronWeights(fhicl::ParameterSet const& pset) { configure(pset); }

void NeutronWeights::configure(fhicl::ParameterSet const& p) {

  fXSecFile    = p.get<std::string>("XSecFile",
                  "xsecHistograms/neutron_inelastic_cross_section_mod.root");
  fNUniverses  = p.get<int>("NUniverses", 1000);
  fSeedBase    = p.get<int>("ZSeedBase", 1337);


  // clamp for safety
  if (fNUniverses < 0)    fNUniverses = 0;
  if (fNUniverses > 5000) fNUniverses = 5000;

  //// reading variables for secondary particles /////

  b_all_mc_pdg         = p.get<std::string>("all_mc_pdg",         b_all_mc_pdg);
  b_all_mc_vx          = p.get<std::string>("all_mc_vx",          b_all_mc_vx);
  b_all_mc_vy          = p.get<std::string>("all_mc_vy",          b_all_mc_vy);
  b_all_mc_vz          = p.get<std::string>("Branch_all_mc_vz",          b_all_mc_vz);
  b_all_mc_endx        = p.get<std::string>("Branch_all_mc_endx",        b_all_mc_endx);
  b_all_mc_endy        = p.get<std::string>("Branch_all_mc_endy",        b_all_mc_endy);
  b_all_mc_endz        = p.get<std::string>("Branch_all_mc_endz",        b_all_mc_endz);
  b_all_mc_E           = p.get<std::string>("Branch_all_mc_E",           b_all_mc_E);
  b_all_mc_end_process = p.get<std::string>("Branch_all_mc_end_process", b_all_mc_end_process);

  _w_neutron_reint.assign(static_cast<size_t>(fNUniverses), 1.f);

  std::cout << "[NeutronWeights] XSecFile=" << fXSecFile
            << " NUniverses=" << fNUniverses
            << " Seed=" << fSeedBase << std::endl;
}

void NeutronWeights::setBranches(TTree* t) {

  _tree = t;  // TTree

  // outputs
  _tree->Branch("w_neutron_reint",           &_w_neutron_reint);
  _tree->Branch("weight_neutron_argon_xsec", &_weight_neutron_argon_xsec);
  _tree->Branch("neutron_track_weights",     &_neutron_track_weights);
  _tree->Branch("neutron_track_len",         &_neutron_track_len);
  _tree->Branch("neutron_ke",                &_neutron_ke);
  _tree->Branch("neutron_last_w",            &_neutron_last_w);
}

void NeutronWeights::resetTTree(TTree*) {

  std::fill(_w_neutron_reint.begin(), _w_neutron_reint.end(), 1.f);
  _weight_neutron_argon_xsec.clear();
  _neutron_track_weights.clear();
  _neutron_track_len.clear();
  _neutron_ke.clear();
  _neutron_last_w.clear();

}

double NeutronWeights::inDetectorSegmentLength(float sx, float sy, float sz,
                                               float ex, float ey, float ez) const
{


  const double step = 0.1; // cm
  const double dx = ex - sx, dy = ey - sy, dz = ez - sz;
  const double L = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (!std::isfinite(L) || L <= 0.0) return 0.0;

  const double ux = dx / L, uy = dy / L, uz = dz / L;
  const int nSteps = static_cast<int>(std::ceil(L / step));
  double inside = 0.0;

  float px = sx, py = sy, pz = sz;
  bool prevInside = isInsideTPC(px, py, pz);

  for (int i = 0; i < nSteps; ++i) {
    const float nx = px + ux * step;
    const float ny = py + uy * step;
    const float nz = pz + uz * step;
    const bool nowInside = isInsideTPC(nx, ny, nz);
    if (prevInside && nowInside) inside += step;
    else if (prevInside && !nowInside) inside += step;
    px = nx; py = ny; pz = nz; prevInside = nowInside;
  }
  return inside;
}

void NeutronWeights::analyzeSlice(art::Event const&,std::vector<ProxyPfpElem_t>&,bool, bool){} // changing from slice to event

void NeutronWeights::analyzeEvent(art::Event const &e, bool fData) {

  if (!_tree) {
    throw cet::exception("NeutronWeights")
      << "TTree is null. setBranches(...) must be called before analyzeEvent().";
  }
  resetTTree(_tree);

  // ---- READ inputs ----
  auto* pdg_v     = getBranchAddr<std::vector<int>>(         "all_mc_pdg",         "NeutronWeights");
  auto* vx_v      = getBranchAddr<std::vector<float>>(       "all_mc_vx",          "NeutronWeights");
  auto* vy_v      = getBranchAddr<std::vector<float>>(       "all_mc_vy",          "NeutronWeights");
  auto* vz_v      = getBranchAddr<std::vector<float>>(       "all_mc_vz",          "NeutronWeights");
  auto* endx_v    = getBranchAddr<std::vector<float>>(       "all_mc_endx",        "NeutronWeights");
  auto* endy_v    = getBranchAddr<std::vector<float>>(       "all_mc_endy",        "NeutronWeights");
  auto* endz_v    = getBranchAddr<std::vector<float>>(       "all_mc_endz",        "NeutronWeights");
  auto* E_v       = getBranchAddr<std::vector<float>>(       "all_mc_E",           "NeutronWeights");
  auto* endproc_v = getBranchAddr<std::vector<std::string>>( "all_mc_end_process", "NeutronWeights");

  const size_t N = pdg_v->size();

  if (N == 0 || N > 1000000) return;

  // ---- build per-neutron inputs ----
  std::vector<double> KE; KE.reserve(N);
  std::vector<double> lengths; lengths.reserve(N);
  std::vector<std::string> endproc; endproc.reserve(N);

  _neutron_ke.reserve(N);
  _neutron_track_len.reserve(N);

  constexpr double Mn = 0.939565; // GeV

  for (size_t i = 0; i < N; ++i) {
    if (pdg_v->at(i) != 2112) continue;

    const double KE_GeV = static_cast<double>(E_v->at(i)) - Mn;
    const double Lcm    = inDetectorSegmentLength(vx_v->at(i),  vy_v->at(i),  vz_v->at(i),
                                                  endx_v->at(i), endy_v->at(i), endz_v->at(i));
    _neutron_ke.push_back(static_cast<float>(KE_GeV));
    _neutron_track_len.push_back(static_cast<float>(Lcm));

    if (KE_GeV < 0.1 || Lcm <= 0.0) continue; // ~100 MeV or no in-FV path

    KE.push_back(KE_GeV);
    lengths.push_back(Lcm);
    endproc.push_back(endproc_v->at(i));
  }

  if (KE.empty()) {
	_weight_neutron_argon_xsec.push_back(1.f);
	return;}


// int seed = 1337;

  // ---- reweighter ----
  NeutronReweighter rw;
  rw.SetXsecFileName(fXSecFile);
//  rw.SetRandomSeed(seed);
  rw.Configure();
  rw.ConfigureEvent(KE, lengths, endproc);


  // CV per-neutron & event (universe = -1)
  _neutron_track_weights.reserve(KE.size());
  for (size_t i = 0; i < KE.size(); ++i) {
    const double nom = rw.GetNominalXsec(KE[i]);
    const double cv  = rw.GetUniverseXsec(KE[i], -1);
    const double w   = rw.CalculateSegmentWeight(KE[i], lengths[i], endproc[i], nom, cv);
    _neutron_track_weights.push_back(static_cast<float>(w));
  }
  const float cv_event =
      std::accumulate(_neutron_track_weights.begin(), _neutron_track_weights.end(),
                      1.f, std::multiplies<float>());
  _weight_neutron_argon_xsec.push_back(cv_event);

  // debug per-neutron (universe 0)
  _neutron_last_w.reserve(KE.size());
  for (size_t i = 0; i < KE.size(); ++i) {
    const double nom = rw.GetNominalXsec(KE[i]);
    const double u0  = rw.GetUniverseXsec(KE[i], 0);
    const double w0  = rw.CalculateSegmentWeight(KE[i], lengths[i], endproc[i], nom, u0);
    _neutron_last_w.push_back(static_cast<float>(w0));
  }

  // multisim event weights
  for (int u = 0; u < fNUniverses; ++u) {
    _w_neutron_reint[static_cast<size_t>(u)] = static_cast<float>(rw.GetEventWeight(u));
  }
}
DEFINE_ART_CLASS_TOOL(NeutronWeights)
} // namespace analysis



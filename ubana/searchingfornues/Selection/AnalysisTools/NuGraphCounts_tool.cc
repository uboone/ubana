#ifndef ANALYSIS_NUGRAPHCOUNTS_CXX
#define ANALYSIS_NUGRAPHCOUNTS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardata/Utilities/FindManyInChainP.h"
#include "larcore/Geometry/WireReadout.h"
#include "ubana/searchingfornues/Selection/CommonDefs/ShowerBranchTagger.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       NuGraphCounts
// File:        NuGraphCounts.cc
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

class NuGraphCounts : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  NuGraphCounts(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~NuGraphCounts(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief Save truth info for event associated to neutrino
     */
  void SaveTruth(art::Event const &e);

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

private:
  //
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fNG2producer; // nugraph2 producer
  art::InputTag fNG2FiltProducer; // nugraph2 producer (filter)

  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

  float _wire2cm, _time2cm;

  // TTree variables

  //nu graph slice hit counts for each semantic label (MIP/HIP/shower/Michel/diffuse) and for background
  int slcng2mip;
  int slcng2hip;
  int slcng2shr;
  int slcng2mcl;
  int slcng2dfs;
  int slcng2bkg;

  //nu graph hit count of hits in PFPs (i.e. clustered) for each semantic label (MIP/HIP/shower/Michel/diffuse) and for background
  int clung2mip;
  int clung2hip;
  int clung2shr;
  int clung2mcl;
  int clung2dfs;
  int clung2bkg;

  //nu graph pfp counts
  std::vector<int> pfng2semlabel;  //NG2 semantic label of the PFP based on majority voting of its hits
  std::vector<float> pfng2mipfrac; //Fraction of PFP hits that are labeled by NuGraph as MIP/HIP/shower/Michel/diffuse and background
  std::vector<float> pfng2hipfrac;
  std::vector<float> pfng2shrfrac;
  std::vector<float> pfng2mclfrac;
  std::vector<float> pfng2dfsfrac;
  std::vector<float> pfng2bkgfrac;
  std::vector<float> pfng2mipavrg; //Average MIP/HIP/shower/Michel/diffuse and background score of hits in the PFP
  std::vector<float> pfng2hipavrg;
  std::vector<float> pfng2shravrg;
  std::vector<float> pfng2mclavrg;
  std::vector<float> pfng2dfsavrg;
  std::vector<float> pfng2bkgavrg;

  // HIP activity counts
  int nhits_r1cm;//total number of hits within 1/3/5/10 cm of neutrino vertex
  int nhits_r3cm;
  int nhits_r5cm;
  int nhits_r10cm;
  int ng2hip_r1cm;//number of HIP-labeled hits within 1/3/5/10 cm of neutrino vertex
  int ng2hip_r3cm;
  int ng2hip_r5cm;
  int ng2hip_r10cm;
  int ng2clu_hip_r1cm;//number of HIP-labeled hits within 1/3/5/10 cm of neutrino vertex that are clustered in PFPs
  int ng2clu_hip_r3cm;
  int ng2clu_hip_r5cm;
  int ng2clu_hip_r10cm;
  int ng2clu_hippfp_r1cm;//number of HIP-labeled hits within 1/3/5/10 cm of neutrino vertex that are clustered in a HIP-labeled PFPs
  int ng2clu_hippfp_r3cm;
  int ng2clu_hippfp_r5cm;
  int ng2clu_hippfp_r10cm;
  std::vector<float> pfng2hip_r1cm;//number of HIP-labeled hits within 1/3/5/10 cm of neutrino vertex in each PFP
  std::vector<float> pfng2hip_r3cm;
  std::vector<float> pfng2hip_r5cm;
  std::vector<float> pfng2hip_r10cm;

};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
NuGraphCounts::NuGraphCounts(const fhicl::ParameterSet &p)
{
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fNG2producer = p.get<art::InputTag>("NG2producer");
  fNG2FiltProducer = p.get<art::InputTag>("NG2FiltProducer");

  // get detector specific properties
  auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  _wire2cm = channelMap.Plane(geo::PlaneID{0,0,0}).WirePitch();
  _time2cm = sampling_rate(clockData) / 1000.0 * detProp.DriftVelocity( detProp.Efield(), detProp.Temperature() );
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void NuGraphCounts::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void NuGraphCounts::analyzeEvent(art::Event const &e, bool fData) {}

void NuGraphCounts::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
  // somehow proxies don't work for the slice-hit association, so go back to old assns
  art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

  // auto GNNDescription = e.getHandle<anab::MVADescription<5>>(art::InputTag("NuGraph", "semantic"));

  // grab vertex neutrino vertex
  TVector3 nuvtx;
  for (auto pfp : slice_pfp_v)
  {
    if ( (pfp->PdgCode() == 12) || (pfp->PdgCode() == 14) ) {
      auto vtx = pfp.get<recob::Vertex>();
      if (vtx.size() != 1) {
	std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	return;
      }	
      // save vertex to array
      nuvtx = TVector3(vtx[0]->position().X(),vtx[0]->position().Y(),vtx[0]->position().Z());
    }// if neutrino PFP
  }

  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
    e,
    art::InputTag("gaushit"), //tag of the hit collection we ran the GNN on
    proxy::withParallelData<anab::FeatureVector<1>>(fNG2FiltProducer),
    proxy::withParallelData<anab::FeatureVector<5>>(fNG2producer));

  nhits_r1cm = 0;
  nhits_r3cm = 0;
  nhits_r5cm = 0;
  nhits_r10cm = 0;
  ng2hip_r1cm = 0;
  ng2hip_r3cm = 0;
  ng2hip_r5cm = 0;
  ng2hip_r10cm = 0;
  ng2clu_hip_r1cm = 0;
  ng2clu_hip_r3cm = 0;
  ng2clu_hip_r5cm = 0;
  ng2clu_hip_r10cm = 0;
  ng2clu_hippfp_r1cm = 0;
  ng2clu_hippfp_r3cm = 0;
  ng2clu_hippfp_r5cm = 0;
  ng2clu_hippfp_r10cm = 0;

  std::vector<int> ng2semclucounts(5,0);
  std::vector<int> ng2semslccounts(5,0);
  int ng2bkgclucounts = 0;
  int ng2bkgslccounts = 0;
  for (auto& h : hitsWithScores) {
    auto scores = h.get<anab::FeatureVector<5>>();
    std::vector<float> ng2semscores;
    for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
    unsigned int sem_label = arg_max(ng2semscores);
    ng2semslccounts[sem_label]++;
    //
    auto fscore = h.get<anab::FeatureVector<1>>();
    if (fscore[0]<0.25) ng2bkgslccounts++;
    //
    auto vtxdistance = searchingfornues::HitPtDistance(nuvtx,h,_wire2cm,_time2cm);
    if (vtxdistance < 1.) nhits_r1cm++;
    if (vtxdistance < 3.) nhits_r3cm++;
    if (vtxdistance < 5.) nhits_r5cm++;
    if (vtxdistance < 10.) nhits_r10cm++;
    if (sem_label==1) {
      if (vtxdistance < 1.) ng2hip_r1cm++;
      if (vtxdistance < 3.) ng2hip_r3cm++;
      if (vtxdistance < 5.) ng2hip_r5cm++;
      if (vtxdistance < 10.) ng2hip_r10cm++;
    }
    //
  }

  for (auto pfp : slice_pfp_v)
  {

    int pfp_ng2hip_r1cm = 0;
    int pfp_ng2hip_r3cm = 0;
    int pfp_ng2hip_r5cm = 0;
    int pfp_ng2hip_r10cm = 0;

    //exclude the neutrino pfp
    if (pfp->PdgCode()!=11 && pfp->PdgCode()!=13) continue;

    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit>> hit_v;
    auto clus_pxy_v = pfp.get<recob::Cluster>();
    if (clus_pxy_v.size() != 0)
    {
      for (auto ass_clus : clus_pxy_v)
      {
        // get cluster proxy
        const auto &clus = clus_proxy[ass_clus.key()];
        auto clus_hit_v = clus.get<recob::Hit>();
        for (const auto &hit : clus_hit_v)
          hit_v.push_back(hit);
      } // for all clusters associated to PFP
    }
    //
    if (hit_v.size()>0) {
      std::vector<int> ng2sempfpcounts(5,0);
      std::vector<float> ng2sempfptotscr(5,0);
      int ng2bkgpfpcounts = 0;
      float ng2bkgpfptotscr = 0.;
      for (auto& hit : hit_v) {
	  auto scores = hitsWithScores[hit.key()].get<anab::FeatureVector<5>>();
	  std::vector<float> ng2semscores;
	  for (size_t i=0;i<scores.size();i++) {
	    ng2semscores.push_back(scores[i]);
	    ng2sempfptotscr[i] += scores[i];
	  }
	  unsigned int sem_label = arg_max(ng2semscores);
	  ng2sempfpcounts[sem_label]++;
	  ng2semclucounts[sem_label]++;
	  //
	  if (sem_label==1) {
	    auto vtxdistance = searchingfornues::HitPtDistance(nuvtx,hit,_wire2cm,_time2cm);
	    if (vtxdistance < 1.)  pfp_ng2hip_r1cm++;
	    if (vtxdistance < 3.)  pfp_ng2hip_r3cm++;
	    if (vtxdistance < 5.)  pfp_ng2hip_r5cm++;
	    if (vtxdistance < 10.) pfp_ng2hip_r10cm++;
	  }
	  //
	  auto fscore = hitsWithScores[hit.key()].get<anab::FeatureVector<1>>();
	  if (fscore[0]<0.25) {
	    ng2bkgclucounts++;
	    ng2bkgpfpcounts++;
	  }
	  ng2bkgpfptotscr += fscore[0];
      }
      pfng2semlabel.push_back(arg_max(ng2sempfpcounts));
      pfng2mipfrac.push_back(float(ng2sempfpcounts[0])/hit_v.size());
      pfng2hipfrac.push_back(float(ng2sempfpcounts[1])/hit_v.size());
      pfng2shrfrac.push_back(float(ng2sempfpcounts[2])/hit_v.size());
      pfng2mclfrac.push_back(float(ng2sempfpcounts[3])/hit_v.size());
      pfng2dfsfrac.push_back(float(ng2sempfpcounts[4])/hit_v.size());
      pfng2bkgfrac.push_back(float(ng2bkgpfpcounts)/hit_v.size());
      pfng2mipavrg.push_back(float(ng2sempfptotscr[0])/hit_v.size());
      pfng2hipavrg.push_back(float(ng2sempfptotscr[1])/hit_v.size());
      pfng2shravrg.push_back(float(ng2sempfptotscr[2])/hit_v.size());
      pfng2mclavrg.push_back(float(ng2sempfptotscr[3])/hit_v.size());
      pfng2dfsavrg.push_back(float(ng2sempfptotscr[4])/hit_v.size());
      pfng2bkgavrg.push_back(float(ng2bkgpfptotscr)/hit_v.size());
      pfng2hip_r1cm.push_back(pfp_ng2hip_r1cm);
      pfng2hip_r3cm.push_back(pfp_ng2hip_r3cm);
      pfng2hip_r5cm.push_back(pfp_ng2hip_r5cm);
      pfng2hip_r10cm.push_back(pfp_ng2hip_r10cm);
      ng2clu_hip_r1cm  += pfp_ng2hip_r1cm;
      ng2clu_hip_r3cm  += pfp_ng2hip_r3cm;
      ng2clu_hip_r5cm  += pfp_ng2hip_r5cm;
      ng2clu_hip_r10cm += pfp_ng2hip_r10cm;
      if (pfng2semlabel.back() == 1) {
	ng2clu_hippfp_r1cm  += pfp_ng2hip_r1cm;
	ng2clu_hippfp_r3cm  += pfp_ng2hip_r3cm;
	ng2clu_hippfp_r5cm  += pfp_ng2hip_r5cm;
	ng2clu_hippfp_r10cm += pfp_ng2hip_r10cm;
      }
    } else {
      pfng2semlabel.push_back(-1);
      pfng2mipfrac.push_back(-1);
      pfng2hipfrac.push_back(-1);
      pfng2shrfrac.push_back(-1);
      pfng2mclfrac.push_back(-1);
      pfng2dfsfrac.push_back(-1);
      pfng2bkgfrac.push_back(-1);
      pfng2mipavrg.push_back(-1);
      pfng2hipavrg.push_back(-1);
      pfng2shravrg.push_back(-1);
      pfng2mclavrg.push_back(-1);
      pfng2dfsavrg.push_back(-1);
      pfng2bkgavrg.push_back(-1);
    }
    //
  }
  //
  slcng2mip = ng2semslccounts[0];
  slcng2hip = ng2semslccounts[1];
  slcng2shr = ng2semslccounts[2];
  slcng2mcl = ng2semslccounts[3];
  slcng2dfs = ng2semslccounts[4];
  slcng2bkg = ng2bkgslccounts;
  //
  clung2mip = ng2semclucounts[0];
  clung2hip = ng2semclucounts[1];
  clung2shr = ng2semclucounts[2];
  clung2mcl = ng2semclucounts[3];
  clung2dfs = ng2semclucounts[4];
  clung2bkg = ng2bkgclucounts;

  return;
}

void NuGraphCounts::setBranches(TTree *_tree)
{
  //
  _tree->Branch("slcng2mip", &slcng2mip, "slcng2mip/I");
  _tree->Branch("slcng2hip", &slcng2hip, "slcng2hip/I");
  _tree->Branch("slcng2shr", &slcng2shr, "slcng2shr/I");
  _tree->Branch("slcng2mcl", &slcng2mcl, "slcng2mcl/I");
  _tree->Branch("slcng2dfs", &slcng2dfs, "slcng2dfs/I");
  _tree->Branch("slcng2bkg", &slcng2bkg, "slcng2bkg/I");
  //
  _tree->Branch("clung2mip", &clung2mip, "clung2mip/I");
  _tree->Branch("clung2hip", &clung2hip, "clung2hip/I");
  _tree->Branch("clung2shr", &clung2shr, "clung2shr/I");
  _tree->Branch("clung2mcl", &clung2mcl, "clung2mcl/I");
  _tree->Branch("clung2dfs", &clung2dfs, "clung2dfs/I");
  _tree->Branch("clung2bkg", &clung2bkg, "clung2bkg/I");
  //
  _tree->Branch("pfng2semlabel", &pfng2semlabel);
  _tree->Branch("pfng2mipfrac", &pfng2mipfrac);
  _tree->Branch("pfng2hipfrac", &pfng2hipfrac);
  _tree->Branch("pfng2shrfrac", &pfng2shrfrac);
  _tree->Branch("pfng2mclfrac", &pfng2mclfrac);
  _tree->Branch("pfng2dfsfrac", &pfng2dfsfrac);
  _tree->Branch("pfng2bkgfrac", &pfng2bkgfrac);
  _tree->Branch("pfng2mipavrg", &pfng2mipavrg);
  _tree->Branch("pfng2hipavrg", &pfng2hipavrg);
  _tree->Branch("pfng2shravrg", &pfng2shravrg);
  _tree->Branch("pfng2mclavrg", &pfng2mclavrg);
  _tree->Branch("pfng2dfsavrg", &pfng2dfsavrg);
  _tree->Branch("pfng2bkgavrg", &pfng2bkgavrg);
  //
  _tree->Branch("nhits_r1cm", &nhits_r1cm, "nhits_r1cm/I");
  _tree->Branch("nhits_r3cm", &nhits_r3cm, "nhits_r3cm/I");
  _tree->Branch("nhits_r5cm", &nhits_r5cm, "nhits_r5cm/I");
  _tree->Branch("nhits_r10cm", &nhits_r10cm, "nhits_r10cm/I");
  _tree->Branch("ng2hip_r1cm", &ng2hip_r1cm, "ng2hip_r1cm/I");
  _tree->Branch("ng2hip_r3cm", &ng2hip_r3cm, "ng2hip_r3cm/I");
  _tree->Branch("ng2hip_r5cm", &ng2hip_r5cm, "ng2hip_r5cm/I");
  _tree->Branch("ng2hip_r10cm", &ng2hip_r10cm, "ng2hip_r10cm/I");
  _tree->Branch("ng2clu_hip_r1cm", &ng2clu_hip_r1cm, "ng2clu_hip_r1cm/I");
  _tree->Branch("ng2clu_hip_r3cm", &ng2clu_hip_r3cm, "ng2clu_hip_r3cm/I");
  _tree->Branch("ng2clu_hip_r5cm", &ng2clu_hip_r5cm, "ng2clu_hip_r5cm/I");
  _tree->Branch("ng2clu_hip_r10cm", &ng2clu_hip_r10cm, "ng2clu_hip_r10cm/I");
  _tree->Branch("ng2clu_hippfp_r1cm", &ng2clu_hippfp_r1cm, "ng2clu_hippfp_r1cm/I");
  _tree->Branch("ng2clu_hippfp_r3cm", &ng2clu_hippfp_r3cm, "ng2clu_hippfp_r3cm/I");
  _tree->Branch("ng2clu_hippfp_r5cm", &ng2clu_hippfp_r5cm, "ng2clu_hippfp_r5cm/I");
  _tree->Branch("ng2clu_hippfp_r10cm", &ng2clu_hippfp_r10cm, "ng2clu_hippfp_r10cm/I");
  _tree->Branch("pfng2hip_r1cm", &pfng2hip_r1cm);
  _tree->Branch("pfng2hip_r3cm", &pfng2hip_r3cm);
  _tree->Branch("pfng2hip_r5cm", &pfng2hip_r5cm);
  _tree->Branch("pfng2hip_r10cm", &pfng2hip_r10cm);
}

void NuGraphCounts::resetTTree(TTree *_tree)
{
  //nu graph slice hit counts
  slcng2mip = std::numeric_limits<int>::min();
  slcng2hip = std::numeric_limits<int>::min();
  slcng2shr = std::numeric_limits<int>::min();
  slcng2mcl = std::numeric_limits<int>::min();
  slcng2dfs = std::numeric_limits<int>::min();
  slcng2bkg = std::numeric_limits<int>::min();
  //nu graph clustered hit counts
  clung2mip = std::numeric_limits<int>::min();
  clung2hip = std::numeric_limits<int>::min();
  clung2shr = std::numeric_limits<int>::min();
  clung2mcl = std::numeric_limits<int>::min();
  clung2dfs = std::numeric_limits<int>::min();
  clung2bkg = std::numeric_limits<int>::min();
  //nu graph pfp counts
  pfng2semlabel.clear();
  pfng2mipfrac.clear();
  pfng2hipfrac.clear();
  pfng2shrfrac.clear();
  pfng2mclfrac.clear();
  pfng2dfsfrac.clear();
  pfng2bkgfrac.clear();
  pfng2mipavrg.clear();
  pfng2hipavrg.clear();
  pfng2shravrg.clear();
  pfng2mclavrg.clear();
  pfng2dfsavrg.clear();
  pfng2bkgavrg.clear();
  //
  nhits_r1cm = std::numeric_limits<int>::min();
  nhits_r3cm = std::numeric_limits<int>::min();
  nhits_r5cm = std::numeric_limits<int>::min();
  nhits_r10cm = std::numeric_limits<int>::min();
  ng2hip_r1cm = std::numeric_limits<int>::min();
  ng2hip_r3cm = std::numeric_limits<int>::min();
  ng2hip_r5cm = std::numeric_limits<int>::min();
  ng2hip_r10cm = std::numeric_limits<int>::min();
  ng2clu_hip_r1cm = std::numeric_limits<int>::min();
  ng2clu_hip_r3cm = std::numeric_limits<int>::min();
  ng2clu_hip_r5cm = std::numeric_limits<int>::min();
  ng2clu_hip_r10cm = std::numeric_limits<int>::min();
  ng2clu_hippfp_r1cm = std::numeric_limits<int>::min();
  ng2clu_hippfp_r3cm = std::numeric_limits<int>::min();
  ng2clu_hippfp_r5cm = std::numeric_limits<int>::min();
  ng2clu_hippfp_r10cm = std::numeric_limits<int>::min();
  pfng2hip_r1cm.clear();
  pfng2hip_r3cm.clear();
  pfng2hip_r5cm.clear();
  pfng2hip_r10cm.clear();
}

DEFINE_ART_CLASS_TOOL(NuGraphCounts)
} // namespace analysis

#endif

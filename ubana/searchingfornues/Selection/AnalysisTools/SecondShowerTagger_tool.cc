#ifndef SELECTION_SECONDSHOWERTAGGER_CXX
#define SELECTION_SECONDSHOWERTAGGER_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "larcore/Geometry/WireReadout.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "ubana/searchingfornues/Selection/CommonDefs/TrackShowerScoreFuncs.h"
#include "ubana/searchingfornues/Selection/CommonDefs/TrackFitterFunctions.h"
#include "ubana/searchingfornues/Selection/CommonDefs/ShowerBranchTagger.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       SecondShowerTagger
    // File:        SecondShowerTagger.cc
    //
    //              A basic selection example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by David Caratelli (davidc@fnal.gov) on 01/30/2019
    //
    ////////////////////////////////////////////////////////////////////////
    
  class SecondShowerTagger : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    SecondShowerTagger(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~SecondShowerTagger(){};
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree) override;

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree) override;
    
  private:

    art::InputTag fClusterproducer;
    art::InputTag fPFPproducer;
    art::InputTag fHitproducer;

    float _wire2cm, _time2cm;

    // NTuple variables
    // charge associated to the largest 2D cluster on each plane
    float _secondshower_U_charge, _secondshower_V_charge, _secondshower_Y_charge;
    // nhit associated to the largest 2D cluster on each plane
    int _secondshower_U_nhit, _secondshower_V_nhit, _secondshower_Y_nhit;
    // 2D distance from vertex for the largest 2D cluster on each plane
    float _secondshower_U_vtxdist, _secondshower_V_vtxdist, _secondshower_Y_vtxdist;
    // ratio of eigenvalues of the two sub-clusters
    float _secondshower_U_eigenratio, _secondshower_V_eigenratio, _secondshower_Y_eigenratio;
    // dot product between vtx -> closest hit in cluster and charge-weighted cluster direction w.r.t. cloest hit in cluster
    float _secondshower_U_dot, _secondshower_V_dot, _secondshower_Y_dot;
    // charge-weighted direction of cluster w.r.t. vertex
    float _secondshower_U_dir, _secondshower_V_dir, _secondshower_Y_dir;
    // number of clusters per plane
    int _secondshower_U_03_ncluster, _secondshower_V_03_ncluster, _secondshower_Y_03_ncluster;
    int _secondshower_U_10_ncluster, _secondshower_V_10_ncluster, _secondshower_Y_10_ncluster;
    int _secondshower_U_20_ncluster, _secondshower_V_20_ncluster, _secondshower_Y_20_ncluster;
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  SecondShowerTagger::SecondShowerTagger(const fhicl::ParameterSet& pset)
  {
    configure(pset);

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
  void SecondShowerTagger::configure(fhicl::ParameterSet const & pset)
  {

    fClusterproducer = pset.get< art::InputTag > ("Clusterproducer", "");
    fHitproducer     = pset.get< art::InputTag > ("Hitproducer"    , "");

}

  
void SecondShowerTagger::analyzeEvent(art::Event const &e, bool fData)
{
  
}

  
  //----------------------------------------------------------------------------
  /// selectEvent
  ///
  /// Arguments:
  ///
  /// art::Event
  /// slice track pointer vector
  /// slice shower pointer vector
  ///
  void SecondShowerTagger::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    std::cout << "[SecondShowerTagger] START analyzeSlice" << std::endl;

    // STEP (1) load neutrino vertex
    TVector3 nuvtx;
    Double_t xyz[3] = {};

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (const auto& pfp_pxy : slice_pfp_v) {
      
      auto PDG = fabs(pfp_pxy->PdgCode());

      if ( (PDG == 12) || (PDG == 14) ) {

	// grab vertex
	auto vtx = pfp_pxy.get<recob::Vertex>();
	if (vtx.size() != 1) {
	  std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	  return;
	}
	
	// save vertex to array
	vtx.at(0)->XYZ(xyz);
	nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);

      }// if neutrino PFP


      
    }// for all PFParticles

    // STEP (1.5): find clusters associated with PFParticles in slice
    std::cout << "[SecondShowerTagger] MESSAGE" << std::endl;
    // clusters associated to each PFParticle
    std::vector<std::pair<int, size_t>> used_cluster_idx_v;
    for (const auto& pfp_pxy : slice_pfp_v) {
      std::cout << "[SecondShowerTagger] New PFParticle " << std::endl;
      //
      float trkscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
      if ((trkscore < 0) || (trkscore >= 0.5)) { std::cout << "[SecondShowerTagger] track-shower score is " << trkscore << " -> track-like. Skip..." << std::endl; continue; }
      //
      auto cluster_v = pfp_pxy.get<recob::Cluster>();
      for (size_t c=0; c < cluster_v.size(); c++) {
	std::cout << "[SecondShowerTagger] \t Cluster index is " << (*cluster_v[c]).ID() << std::endl;
	std::pair<int, size_t> used_cluster_idx = {(*cluster_v[c]).ID(), (*cluster_v[c]).Plane().Plane};
        used_cluster_idx_v.push_back( used_cluster_idx );
//        used_cluster_idx_v.push_back( (*cluster_v[c]).ID() );
        std::cout << "cluster_v[c].Plane().Plane" << (*cluster_v[c]).Plane().Plane << std::endl;
        std::cout << "cluster_v[c].NHits()" << (*cluster_v[c]).NHits() << std::endl;
      }// for each associated cluster
    }// for all PFParticles
    
    //const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at( pfp_ptr.key() );

    /*
    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

    // grab clusters associated with PFParticles
    art::FindManyP<recob::Cluster> pfp_cluster_assn_v(pfp_h, e, fPFPproducer);

    // store indices of all clusters associated to pfparticles
    std::vector<size_t> cluster_idx_v;

    // loop over all PFPs 
    for (const auto& pfp_pxy : slice_pfp_v) {

      const std::vector<art::Ptr<recob::Cluster> >& clusters(pfp_cluster_assn_v.at(pfp_pxy.index()));
      // get cluster indices
      for (auto& cl : clusters) shr_clusters.push_back(*cl);
      for (size_t c=0; c < clusters.size(); c++) {
	
	cluster_idx_v.push_back(clusters[c].index());
      }
    }// for all pfparticles
    */

    
    // STEP (2) : find, per plane, largest 2D cluster

    // grab clustnoselectionfilter_run1_overlay_cc0pinp.fclers themselves
    auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fClusterproducer);
    // get hits associated to clusters
    art::FindManyP<recob::Hit> clus_hit_assn_v(cluster_h, e, fClusterproducer);
    // grab hits themselves
    auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitproducer);

    _secondshower_U_03_ncluster = 0;
    _secondshower_V_03_ncluster = 0;
    _secondshower_Y_03_ncluster = 0;

    _secondshower_U_10_ncluster = 0;
    _secondshower_V_10_ncluster = 0;
    _secondshower_Y_10_ncluster = 0;

    _secondshower_U_20_ncluster = 0;
    _secondshower_V_20_ncluster = 0;
    _secondshower_Y_20_ncluster = 0;
    
    for (unsigned int pl=0; pl < 3; pl++) {
    
      // reset plane variables
      float mosthits = 0;
      float charge = 0.;
      float vtxdistancemin = 1e6; // min distance to neutrino vertex

      TVectorD eigenVal(2);
      TMatrixD eigenVec(2,2);

      float dot = 0;
      float angle = 0;
      
      // loop through clusters
      for (size_t c=0; c < cluster_h->size(); c++) {
	//auto clus = cluster_h->at(c);

	// get associated hits
	auto clus_hit_v = clus_hit_assn_v.at( c );
	
	// create vector of gaushits corresponding to new proton hits
	auto gaushit_hit_v = searchingfornues::getGaussHits(clus_hit_v, hit_h);
	
	// no hits associated to the cluster? weird...skip...
	if (gaushit_hit_v.size() == 0) continue;
	
	unsigned int PLANE = (gaushit_hit_v.at(0))->WireID().Plane;
	
	// focus one plane at a time
	if (PLANE != pl) continue;

	// is this cluster already used? -> SKIP
	std::cout << "Plane " << pl  << std::endl;
        // std::vector<size_t> cluster_idx = {static_cast<size_t>(c), pl};
	std::pair<int, size_t> cluster_idx = {cluster_h->at(c).ID(), pl}; // Giuseppe's suggestion
	if (std::find(used_cluster_idx_v.begin(), used_cluster_idx_v.end(), cluster_idx) != used_cluster_idx_v.end()){
	  std::cout << "[SecondShowerTagger] cluster " << c << " is already used in a reco PFParticle" << std::endl;
	  continue;
	}
	std::cout << "[SecondShowerTagger] cluster " << c << " is an un-matched cluster" << std::endl;
	
	int clushits = clus_hit_v.size();

	if (clushits >= 3) {
	  if (pl==0) { _secondshower_U_03_ncluster += 1; }
	  if (pl==1) { _secondshower_V_03_ncluster += 1; }
	  if (pl==2) { _secondshower_Y_03_ncluster += 1; }
	}

	if (clushits >= 10) {
	  if (pl==0) { _secondshower_U_10_ncluster += 1; }
	  if (pl==1) { _secondshower_V_10_ncluster += 1; }
	  if (pl==2) { _secondshower_Y_10_ncluster += 1; }
	}

	if (clushits >= 20) {
	  if (pl==0) { _secondshower_U_20_ncluster += 1; }
	  if (pl==1) { _secondshower_V_20_ncluster += 1; }
	  if (pl==2) { _secondshower_Y_20_ncluster += 1; }
	}

	// we are focusing on the largest cluster per plane
	// continue if not this one
	if (clushits <= mosthits) continue;

	mosthits = clushits;
	
	// reset variables
	charge = 0;
	vtxdistancemin = 1e6;
	
	for (size_t hi=0; hi < gaushit_hit_v.size(); hi++) {
	  auto hit = gaushit_hit_v.at(hi);
	  //gammaWire += hit->WireID().Wire * _wire2cm * hit->Integral();
	  //gammaTime += (hit->PeakTime() - trigger_offset(clockData))  * _time2cm * hit->Integral();
	  charge += hit->Integral();
	  auto vtxdistance = searchingfornues::HitPtDistance(nuvtx,hit,_wire2cm,_time2cm);
	  if (vtxdistance < vtxdistancemin) { vtxdistancemin = vtxdistance; }
	}

	searchingfornues::PCA(clus_hit_v, _wire2cm, _time2cm, eigenVal, eigenVec);
	
	// measure alignment between cluster direction and vtx -> cluster start
	dot = searchingfornues::ClusterVtxAlignment(nuvtx, clus_hit_v,_wire2cm,_time2cm);

	// save charge-weieghted 2D direction of cluster hits to neutrino vertex 
	// this will be used to compare the direction with the candidate shower direction
	angle = searchingfornues::ClusterVtxDirection(nuvtx, clus_hit_v,_wire2cm,_time2cm);
	
      }// for all 2D clusters

      if (pl==0) {
	_secondshower_U_charge  = charge;
	_secondshower_U_nhit    = mosthits;
	_secondshower_U_vtxdist = vtxdistancemin;
	_secondshower_U_eigenratio = eigenVal(1) / eigenVal(0);
	_secondshower_U_dot     = dot;
	_secondshower_U_dir     = angle;
      }
      if (pl==1) {
	_secondshower_V_charge  = charge;
	_secondshower_V_nhit    = mosthits;
	_secondshower_V_vtxdist = vtxdistancemin;
	_secondshower_V_eigenratio = eigenVal(1) / eigenVal(0);
	_secondshower_V_dot     = dot;
	_secondshower_V_dir     = angle;
      }
      if (pl==2) {
	_secondshower_Y_charge  = charge;
	_secondshower_Y_nhit    = mosthits;
	_secondshower_Y_vtxdist = vtxdistancemin;
	_secondshower_Y_eigenratio = eigenVal(1) / eigenVal(0);
	_secondshower_Y_dot     = dot;
	_secondshower_Y_dir     = angle;
      }

    }// for all planes
    
    return;
  }
  
  
  void SecondShowerTagger::setBranches(TTree* _tree) {

    _tree->Branch("secondshower_U_charge"    , &_secondshower_U_charge    , "secondshower_U_charge/F"    );
    _tree->Branch("secondshower_U_nhit"      , &_secondshower_U_nhit      , "secondshower_U_nhit/I"      );
    _tree->Branch("secondshower_U_vtxdist"   , &_secondshower_U_vtxdist   , "secondshower_U_vtxdist/F"   );
    _tree->Branch("secondshower_U_eigenratio", &_secondshower_U_eigenratio, "secondshower_U_eigenratio/F");
    _tree->Branch("secondshower_U_dot"       , &_secondshower_U_dot       , "secondshower_U_dot/F"       );
    _tree->Branch("secondshower_U_dir"       , &_secondshower_U_dir       , "secondshower_U_dir/F"       );

    _tree->Branch("secondshower_V_charge"    , &_secondshower_V_charge    , "secondshower_V_charge/F"    );
    _tree->Branch("secondshower_V_nhit"      , &_secondshower_V_nhit      , "secondshower_V_nhit/I"      );
    _tree->Branch("secondshower_V_vtxdist"   , &_secondshower_V_vtxdist   , "secondshower_V_vtxdist/F"   );
    _tree->Branch("secondshower_V_eigenratio", &_secondshower_V_eigenratio, "secondshower_V_eigenratio/F");
    _tree->Branch("secondshower_V_dot"       , &_secondshower_V_dot       , "secondshower_V_dot/F"       );
    _tree->Branch("secondshower_V_dir"       , &_secondshower_V_dir       , "secondshower_V_dir/F"       );

    _tree->Branch("secondshower_Y_charge"    , &_secondshower_Y_charge    , "secondshower_Y_charge/F"    );
    _tree->Branch("secondshower_Y_nhit"      , &_secondshower_Y_nhit      , "secondshower_Y_nhit/I"      );
    _tree->Branch("secondshower_Y_vtxdist"   , &_secondshower_Y_vtxdist   , "secondshower_Y_vtxdist/F"   );
    _tree->Branch("secondshower_Y_eigenratio", &_secondshower_Y_eigenratio, "secondshower_Y_eigenratio/F");
    _tree->Branch("secondshower_Y_dot"       , &_secondshower_Y_dot       , "secondshower_Y_dot/F"       );
    _tree->Branch("secondshower_Y_dir"       , &_secondshower_Y_dir       , "secondshower_Y_dir/F"       );

    _tree->Branch("secondshower_U_03_ncluster"      , &_secondshower_U_03_ncluster      , "secondshower_U_03_ncluster/I"      );
    _tree->Branch("secondshower_V_03_ncluster"      , &_secondshower_V_03_ncluster      , "secondshower_V_03_ncluster/I"      );
    _tree->Branch("secondshower_Y_03_ncluster"      , &_secondshower_Y_03_ncluster      , "secondshower_Y_03_ncluster/I"      );

    _tree->Branch("secondshower_U_10_ncluster"      , &_secondshower_U_10_ncluster      , "secondshower_U_10_ncluster/I"      );
    _tree->Branch("secondshower_V_10_ncluster"      , &_secondshower_V_10_ncluster      , "secondshower_V_10_ncluster/I"      );
    _tree->Branch("secondshower_Y_10_ncluster"      , &_secondshower_Y_10_ncluster      , "secondshower_Y_10_ncluster/I"      );

    _tree->Branch("secondshower_U_20_ncluster"      , &_secondshower_U_20_ncluster      , "secondshower_U_20_ncluster/I"      );
    _tree->Branch("secondshower_V_20_ncluster"      , &_secondshower_V_20_ncluster      , "secondshower_V_20_ncluster/I"      );
    _tree->Branch("secondshower_Y_20_ncluster"      , &_secondshower_Y_20_ncluster      , "secondshower_Y_20_ncluster/I"      );


    return;
  }

  void SecondShowerTagger::resetTTree(TTree* _tree) {

    _secondshower_U_charge  = std::numeric_limits<float>::lowest();
    _secondshower_U_nhit    = std::numeric_limits<int>::lowest();
    _secondshower_U_vtxdist = std::numeric_limits<float>::max();
    _secondshower_U_eigenratio  = std::numeric_limits<float>::lowest();
    _secondshower_U_dot     = std::numeric_limits<float>::lowest();
    _secondshower_U_dir     = std::numeric_limits<float>::lowest();

    _secondshower_V_charge  = std::numeric_limits<float>::lowest();
    _secondshower_V_nhit    = std::numeric_limits<int>::lowest();
    _secondshower_V_vtxdist = std::numeric_limits<float>::max();
    _secondshower_V_eigenratio  = std::numeric_limits<float>::lowest();
    _secondshower_V_dot     = std::numeric_limits<float>::lowest();
    _secondshower_V_dir     = std::numeric_limits<float>::lowest();

    _secondshower_Y_charge  = std::numeric_limits<float>::lowest();
    _secondshower_Y_nhit    = std::numeric_limits<int>::lowest();
    _secondshower_Y_vtxdist = std::numeric_limits<float>::max();
    _secondshower_Y_eigenratio  = std::numeric_limits<float>::lowest();
    _secondshower_Y_dot     = std::numeric_limits<float>::lowest();
    _secondshower_Y_dir     = std::numeric_limits<float>::lowest();

    _secondshower_U_03_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_V_03_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_Y_03_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_U_10_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_V_10_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_Y_10_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_U_20_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_V_20_ncluster = std::numeric_limits<int>::lowest();
    _secondshower_Y_20_ncluster = std::numeric_limits<int>::lowest();

    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(SecondShowerTagger)
} // namespace selection

#endif

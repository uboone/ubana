#ifndef ANALYSIS_COSMICIP_CXX
#define ANALYSIS_COSMICIP_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

// backtracking tools
#include "ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"
#include "ubana/searchingfornues/Selection/CommonDefs/TrackShowerScoreFuncs.h"

#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       CosmicIP
    // File:        CosmicIP.cc
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

  class CosmicIP : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    CosmicIP(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~CosmicIP(){ };
    
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

    //void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);

    //bool TrackInTime(const art::Ptr<recob::Track>& pfp_track_assn_v);

    art::InputTag fPFPproducer;
    std::vector<art::InputTag> fPFPproducersInit;

    float fTrkShrScore;
    float _CosmicIP;          // 3D distance of shower start from closest spacepoint of primary muon (i.e. cosmic)
    float _CosmicIPAll3D;     // 3D distance of shower start from closest spacepoint of any pfp not in the neutrino slice
    float _CosmicDirAll3D;    // cosine of 3D direction difference between shower and closest pfp not in the neutrino slice
    float _CosmicIPAll2DEnds; // 2D distance of shower cluster endpoints from closest cluster endpoints from any pfp not in the neutrino slice (minimum fron all planes)
    float _CosmicDirAll2DEnds;// cosine of 2D direction difference between shower cluster and closest cluster from pfps not in the neutrino slice
    float _CosmicIPAll2DOvlp; // 2D distance of shower cluster endpoints from closest line connecting cluster endpoints from any pfp not in the neutrino slice (minimum fron all planes)
    float _CosmicDirAll2DOvlp;// cosine of 2D direction difference between shower cluster and cluster whose endpoints for the closest line 
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  CosmicIP::CosmicIP(const fhicl::ParameterSet& p)
  {

    fPFPproducer      = p.get< art::InputTag >("PFPproducer");
    fPFPproducersInit = p.get< std::vector<art::InputTag> >("PFPproducersInit");
    fTrkShrScore      = p.get< float >("TrkShrScore");
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void CosmicIP::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void CosmicIP::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {

    // std::cout << "[NEW EVENT]" << e.event() << std::endl;
    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);
    float wire2cm = channelMap.Plane(geo::PlaneID{0,0,0}).WirePitch();
    float time2cm = sampling_rate(clockData) / 1000.0 * detProp.DriftVelocity( detProp.Efield(), detProp.Temperature() );

    // set defaults
    _CosmicIP = 9999.;
    _CosmicIPAll3D = 9999.;
    _CosmicDirAll3D = 9999.;
    _CosmicIPAll2DEnds = 9999.;
    _CosmicDirAll2DEnds = 9999.;
    _CosmicIPAll2DOvlp = 9999.;
    _CosmicDirAll2DOvlp = 9999.;

    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

    // grab clusters associated with PFParticles
    art::FindManyP<recob::Cluster> pfp_cluster_assn_v(pfp_h, e, fPFPproducer);

    std::vector<size_t> slice_pfp_keys;
    std::vector<recob::Cluster> shr_clusters;

    // first get the candidate shower's start point
    float shr_candidate_energy = 0.;
    TVector3 shrStart;
    TVector3 shrDir;
    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++) {
      
      auto const &pfp_pxy = slice_pfp_v.at(i_pfp);
      slice_pfp_keys.push_back(pfp_pxy.index());
      
      auto PDG = fabs(pfp_pxy->PdgCode());
      
      // skip neutrino PFP
      if ((PDG == 12) || (PDG == 14))
	continue;
      
      // grab shower/track score
      auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
      
      // 1 -> track-like
      if (trkshrscore > fTrkShrScore)
	continue;
      
      if ( pfp_pxy.get<recob::Shower>().size() == 0 )
	continue;
      
      auto const &shr = pfp_pxy.get<recob::Shower>().at(0);
      
      if (shr->Energy()[2] > shr_candidate_energy) {
	shr_candidate_energy = shr->Energy()[2];
	shrStart = shr->ShowerStart();
	shrDir = shr->Direction();
	const std::vector<art::Ptr<recob::Cluster> >& clusters(pfp_cluster_assn_v.at(pfp_pxy.index()));
	for (auto& cl : clusters) shr_clusters.push_back(*cl);
	// for (auto& cl : clusters) std::cout << "shr plane=" << cl->Plane().Plane << " start(w,t,a)=" << cl->StartWire() << ", " << cl->StartTick() << ", " << cl->StartAngle() << " end(w,t,a)=" << cl->EndWire() << ", " << cl->EndTick() << ", " << cl->EndAngle() << std::endl;
      }
    }// for all PFParticles

    if (shr_candidate_energy == 0)
      return;
    

    //loop over input tags (only one should matter), then over slices, and then over pfps
    for (auto tag : fPFPproducersInit) {
      art::Handle<std::vector<recob::PFParticle>> inputPfp;
      if (e.getByLabel(tag,inputPfp)==false) continue;

      //std::cout << "tag=" << tag << std::endl;
      auto const& inputSlice = e.getValidHandle<std::vector<recob::Slice>>(tag);
      art::FindOneP<larpandoraobj::PFParticleMetadata> assocPfpMetadata(inputPfp, e, tag);
      art::FindManyP<recob::PFParticle> assocSlicePfp(inputSlice, e, tag);
      art::FindManyP<recob::Cluster> assocPfpCluster(inputPfp, e, tag);

      std::vector<art::Ptr<recob::Slice> > slice_ptr_v;
      std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
      art::fill_ptr_vector(slice_ptr_v,inputSlice);
      art::fill_ptr_vector(pfp_ptr_v,inputPfp);
      std::map<int, art::Ptr<recob::PFParticle>> pfpmap;
      lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfp_ptr_v,pfpmap);

      // grab spacepoints associated with PFParticles
      art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(inputPfp, e, tag);
      // grab tracks associated with PFParticles
      art::FindManyP<recob::Track> pfp_track_assn_v(inputPfp, e, tag);

      for (auto& slice : slice_ptr_v) {
	const std::vector<art::Ptr<recob::PFParticle> >& nuSlicePFPs(assocSlicePfp.at(slice.key()));
	for (auto& pfp : nuSlicePFPs) {
	  auto pfp_parent = lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpmap,pfp);

	  //skip pfps whose parent is not a clear cosmic
	  auto parentParticleMetadata = assocPfpMetadata.at(pfp_parent.key());
	  auto parentParticlePropertiesMap = parentParticleMetadata->GetPropertiesMap();
	  //std::cout << "ipfp self=" << pfp->Self() << " primary=" << pfp->IsPrimary() << " pdg=" << pfp->PdgCode() << std::endl;
	  bool clearCosmic = false;
	  for (auto it = parentParticlePropertiesMap.begin(); it != parentParticlePropertiesMap.end(); it++) {
	    if ( it->first == "IsClearCosmic" ) {
	      clearCosmic = true;
	      break;
	    }
	  }
	  if (!clearCosmic) continue;

	  //std::cout << "slice key=" << slice.key() << " pfp key=" << pfp.key() << " parent key=" << pfp_parent.key() << std::endl;
	  // get spacepoints associated to PFParticle
	  //const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
	  const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at( pfp.key() );

	  float dmin = 9999.;
      
	  // loop through spacepoints, find closest ditance to shower start point
	  for (size_t s=0; s < spacepoint_ptr_v.size(); s++) {
	    auto sps = spacepoint_ptr_v[s];
	    auto spsx = sps->XYZ()[0];
	    auto spsy = sps->XYZ()[1];
	    auto spsz = sps->XYZ()[2];
	    double d = sqrt(  ((spsx-shrStart.X())*(spsx-shrStart.X())) + 
			      ((spsy-shrStart.Y())*(spsy-shrStart.Y())) +
			      ((spsz-shrStart.Z())*(spsz-shrStart.Z())) );
	    if (d < dmin) {
	      dmin = d;
	    }
	  }// for all spacepoints

	  // start from primary PFParticles
	  if (pfp->IsPrimary() == true) {
	    // cut on PDG code. We do not want the neutrino candidate!
	    if (pfp->PdgCode() == 13) {
	      if (dmin < _CosmicIP) _CosmicIP = dmin;
	    }
	  }

	  // make sure this pfp is NOT in the neutrino slice
	  //if (std::find(slice_pfp_keys.begin(),slice_pfp_keys.end(),p)!=slice_pfp_keys.end()) continue;

	  if (dmin < _CosmicIPAll3D) {
	    _CosmicIPAll3D = dmin;
	    const std::vector<art::Ptr<recob::Track> >& tkv(pfp_track_assn_v.at(pfp.key()));
	    if (tkv.size()>0) {
	      _CosmicDirAll3D = tkv[0]->StartDirection<TVector3>().Dot(shrDir);
	    }
	  }

	  float dmin2dends = 9999.;
	  float dmin2dovlp = 9999.;
	  float dirmin2dends = 9999.;
	  float dirmin2dovlp = 9999.;
	  const std::vector<art::Ptr<recob::Cluster> >& clusters(assocPfpCluster.at(pfp.key()));
	  // loop over pfp clusters
	  for (auto& cl : clusters) {
	    // std::cout << "plane=" << cl->Plane().Plane
	    // 	  << " start(w,t,a)=" << cl->StartWire() << ", " << cl->StartTick() << ", " << cl->StartAngle()
	    // 	  << " end(w,t,a)=" << cl->EndWire() << ", " << cl->EndTick() << ", " << cl->EndAngle() << std::endl;
	    float clsw = cl->StartWire()*wire2cm;
	    float clew = cl->EndWire()*wire2cm;
	    float clst = cl->StartTick()*time2cm;
	    float clet = cl->EndTick()*time2cm;
	    // loop over shower clusters
	    for (auto& scl : shr_clusters) {
	      //make sure we compare clusters on the same plane
	      if (scl.Plane().Plane!=cl->Plane().Plane) continue;
	      float shsw = scl.StartWire()*wire2cm;
	      float shew = scl.EndWire()*wire2cm;
	      float shst = scl.StartTick()*time2cm;
	      float shet = scl.EndTick()*time2cm;
	      //
	      // compute distance between end points, save minimum distance and direction difference
	      float dss = sqrt( (clsw-shsw)*(clsw-shsw) + (clst-shst)*(clst-shst) );
	      float dee = sqrt( (clew-shew)*(clew-shew) + (clet-shet)*(clet-shet) );
	      float dse = sqrt( (clsw-shew)*(clsw-shew) + (clst-shet)*(clst-shet) );
	      float des = sqrt( (clew-shsw)*(clew-shsw) + (clet-shst)*(clet-shst) );
	      // std::cout << "dss=" << dss << " dee=" << dee << " dse=" << dse << " des=" << des << std::endl;
	      if (dss<dmin2dends) {
		dmin2dends = dss;
		dirmin2dends = std::cos(scl.StartAngle())*std::cos(cl->StartAngle()) + std::sin(scl.StartAngle())*std::sin(cl->StartAngle());
	      }
	      if (dee<dmin2dends) {
		dmin2dends = dee;
		dirmin2dends = std::cos(scl.EndAngle())*std::cos(cl->EndAngle()) + std::sin(scl.EndAngle())*std::sin(cl->EndAngle());
	      }
	      if (dse<dmin2dends) {
		dmin2dends = dse;
		dirmin2dends = std::cos(scl.EndAngle())*std::cos(cl->StartAngle()) + std::sin(scl.EndAngle())*std::sin(cl->StartAngle());
	      }
	      if (des<dmin2dends) {
		dmin2dends = des;
		dirmin2dends = std::cos(scl.StartAngle())*std::cos(cl->EndAngle()) + std::sin(scl.StartAngle())*std::sin(cl->EndAngle());
	      }
	      //
	      // consider clusters that are overlapping in time and wires, find the one with minimum distance between shower end points and line connecting endpoints of overlapping cluster, save minimum distance and direction difference
	      bool startIn = (shsw>=std::min(clsw,clew) && shsw<=std::max(clsw,clew) && shst>=std::min(clst,clet) && shst<=std::max(clst,clet));
	      bool endIn = (shew>=std::min(clsw,clew) && shew<=std::max(clsw,clew) && shet>=std::min(clst,clet) && shet<=std::max(clst,clet));
	      // std::cout << "startIn=" << startIn << " endIn=" << endIn << std::endl;
	      if (startIn) {
		float distS = std::abs( (clet-clst)*shsw - (clew-clsw)*shst + clew*clst - clet*clsw )/std::sqrt( (clet-clst)*(clet-clst) + (clew-clsw)*(clew-clsw) );
		float dangleSS = std::cos(scl.StartAngle())*std::cos(cl->StartAngle()) + std::sin(scl.StartAngle())*std::sin(cl->StartAngle());
		float dangleSE = std::cos(scl.StartAngle())*std::cos(cl->EndAngle()) + std::sin(scl.StartAngle())*std::sin(cl->EndAngle());
		// std::cout << "distS=" << distS << " dangleSS=" << dangleSS << " dangleSE=" << dangleSE << std::endl;
		if (distS<dmin2dovlp) {
		  dmin2dovlp = distS;
		  dirmin2dovlp = ((fabs(dangleSS)>fabs(dangleSE)) ? dangleSS : dangleSE);
		}
	      }
	      if (endIn) {
		float distE = std::abs( (clet-clst)*shew - (clew-clsw)*shet + clew*clst - clet*clsw )/std::sqrt( (clet-clst)*(clet-clst) + (clew-clsw)*(clew-clsw) );
		float dangleES = std::cos(scl.EndAngle())*std::cos(cl->StartAngle()) + std::sin(scl.EndAngle())*std::sin(cl->StartAngle());
		float dangleEE = std::cos(scl.EndAngle())*std::cos(cl->EndAngle()) + std::sin(scl.EndAngle())*std::sin(cl->EndAngle());
		// std::cout << "distE=" << distE << " dangleES=" << dangleES << " dangleEE=" << dangleEE << std::endl;
		if (distE<dmin2dovlp) {
		  dmin2dovlp = distE;
		  dirmin2dovlp = ((fabs(dangleES)>fabs(dangleEE)) ? dangleES : dangleEE);
		}
	      }
	      // std::cout << "dmin2dovlp=" << dmin2dovlp << " dirmin2dovlp=" << dirmin2dovlp << std::endl;
	    }
	  }

	  if (dmin2dends < _CosmicIPAll2DEnds) {
	    _CosmicIPAll2DEnds = dmin2dends;
	    _CosmicDirAll2DEnds = dirmin2dends;
	  }
	  if (dmin2dovlp < _CosmicIPAll2DOvlp) {
	    _CosmicIPAll2DOvlp = dmin2dovlp;
	    _CosmicDirAll2DOvlp = dirmin2dovlp;
	  }



	}//pfps
      }//slices
    }//tag
    return;
  }

  void CosmicIP::analyzeEvent(art::Event const &e, bool fData)
  {
    // std::cout << "analyze event" << std::endl;
  }

  void CosmicIP::setBranches(TTree* _tree)
  {
    _tree->Branch("CosmicIP",&_CosmicIP,"CosmicIP/F");
    _tree->Branch("CosmicIPAll3D",&_CosmicIPAll3D,"CosmicIPAll3D/F");
    _tree->Branch("CosmicDirAll3D",&_CosmicDirAll3D,"CosmicDirAll3D/F");
    _tree->Branch("CosmicIPAll2DEnds",&_CosmicIPAll2DEnds,"CosmicIPAll2DEnds/F");
    _tree->Branch("CosmicDirAll2DEnds",&_CosmicDirAll2DEnds,"CosmicDirAll2DEnds/F");
    _tree->Branch("CosmicIPAll2DOvlp",&_CosmicIPAll2DOvlp,"CosmicIPAll2DOvlp/F");
    _tree->Branch("CosmicDirAll2DOvlp",&_CosmicDirAll2DOvlp,"CosmicDirAll2DOvlp/F");
  }
  
  void CosmicIP::resetTTree(TTree* _tree)
  {
    _CosmicIP = 9999.;
    _CosmicIPAll3D = 9999.;
    _CosmicDirAll3D = 9999.;
    _CosmicIPAll2DEnds = 9999.;
    _CosmicDirAll2DEnds = 9999.;
    _CosmicIPAll2DOvlp = 9999.;
    _CosmicDirAll2DOvlp = 9999.;
  }

  
  DEFINE_ART_CLASS_TOOL(CosmicIP)
} // namespace analysis

#endif

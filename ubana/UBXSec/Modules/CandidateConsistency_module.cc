////////////////////////////////////////////////////////////////////////
// Class:       CandidateConsistency
// Plugin Type: producer (art v2_05_00)
// File:        CandidateConsistency_module.cc
//
// Generated at Sun Dec 31 07:59:34 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class CandidateConsistency
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module that tags TPCObjects if thet are not cosntistent with being neutrino candidate
 * 
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version producer (art v2_05_00)
 *
 * \date 2017/03/10
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Sun Dec 31 07:59:34 2017
 *
 */

#include "art/Framework/Core/EDProducer.h"
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <algorithm>
#include <memory>

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "ubobj/UBXSec/FlashMatch.h"
#include "ubobj/UBXSec/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "ubobj/UBXSec/TPCObject.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

//Algorithms include
#include "ubana/UBXSec/Algorithms/UBXSecHelper.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "ubana/UBXSec/Algorithms/StoppingMuonTaggerHelper.h"
#include "ubana/UBXSec/HitCosmicTag/Base/DataTypes.h"
#include "ubana/UBXSec/HitCosmicTag/Base/CosmicTagManager.h"
#include "ubana/UBXSec/HitCosmicTag/Algorithms/StopMuMichel.h"
#include "ubana/UBXSec/HitCosmicTag/Algorithms/StopMuBragg.h"
#include "ubana/UBXSec/HitCosmicTag/Algorithms/CosmicSimpleMIP.h"

#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "ubcore/Geometry/UBOpReadoutMap.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"

const std::vector<float> endPt1 = {-9999., -9999., -9999.};
const std::vector<float> endPt2 = {-9999., -9999., -9999.};

class CandidateConsistency;


class CandidateConsistency : public art::EDProducer {
public:
  explicit CandidateConsistency(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CandidateConsistency(CandidateConsistency const &) = delete;
  CandidateConsistency(CandidateConsistency &&) = delete;
  CandidateConsistency & operator = (CandidateConsistency const &) = delete;
  CandidateConsistency & operator = (CandidateConsistency &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string _tpcobject_producer;
  std::string _shower_producer;
  std::string _track_producer;

  double      _tolerance;
  double      _dqds_threshold;
  double      _distance_cut;
  double      _perc_losen_hit_cut;
  double      _linearity_threshold;
  double      _linearity_threshold_track;
  double      _dqds_average_cut;
  int         _good_ch_status;
  bool        _debug;

  ::flashana::FlashMatchManager       _mgr;

  void ContainPoint(geo::Point_t&);
  std::vector<double> GetFlashHypo(art::Ptr<recob::Track>, double);
  bool IsCathodeCrossing(art::Ptr<recob::Track>, double &);

  ::art::ServiceHandle<geo::Geometry> geo;
  geo::WireReadoutGeom const* _channelMap = &art::ServiceHandle<geo::WireReadout const>()->Get();

  ::cosmictag::CosmicTagManager _ct_manager;

};


CandidateConsistency::CandidateConsistency(fhicl::ParameterSet const & p) : EDProducer{p}
{

  _tpcobject_producer        = p.get<std::string>("TPCObjectProducer",  "TPCObjectMaker::UBXSec");
  _shower_producer           = p.get<std::string>("ShowerProducer",     "pandoraNu::UBXSec");
  _track_producer            = p.get<std::string>("TrackProducer",      "pandoraNu::UBXSec");

  _tolerance                 = p.get<double>("Tolerance", 5.);
  _dqds_threshold            = p.get<double>("DqDsThreshold", 60000);
  _distance_cut              = p.get<double>("DistanceCut", 8);
  _perc_losen_hit_cut        = p.get<double>("PercentageLosenHits", 30);
  _linearity_threshold       = p.get<double>("LinearityThreshold", 0.7);
  _linearity_threshold_track = p.get<double>("LinearityThresholdTrack", 0.9);
  _dqds_average_cut          = p.get<double>("DqDsAverageThreshold", 90000);
  _good_ch_status            = p.get<double>("GoodChannelStatus", 4);

  _debug                     = p.get<bool>("DebugMode", true);

  _ct_manager.Configure(p.get<cosmictag::Config_t>("CosmicTagManager"));

  _mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));

  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<ubana::TPCObject, anab::CosmicTag>>();

}

void CandidateConsistency::produce(art::Event & e)
{

  if(_debug) std::cout << "[CandidateConsistency] Starts." << std::endl;

  // Instantiate the output
  std::unique_ptr<std::vector<anab::CosmicTag>> cosmicTagVector (new std::vector<anab::CosmicTag>);
  std::unique_ptr<art::Assns<ubana::TPCObject, anab::CosmicTag>> assnOutCosmicTagTPCObject(new art::Assns<ubana::TPCObject, anab::CosmicTag>);

  const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  // Construct map wire->channel for collection plane
  std::map<int,int> wire_to_channel;
  for (unsigned int ch = 0; ch < 8256; ch++) {
    for ( auto wire_id : _channelMap->ChannelToWire(ch)) {
      if (wire_id.Plane == 2) {
        wire_to_channel[wire_id.Wire] = ch;
      }
    }
  }

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    throw cet::exception(__PRETTY_FUNCTION__) << "Cannote locate ubana::TPCObject." << std::endl;
  }

  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);

  std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
  art::fill_ptr_vector(tpcobj_v, tpcobj_h);

  // Collect tracks and map tracks->hits
  lar_pandora::TrackVector track_v;
  lar_pandora::TracksToHits tracks_to_hits;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _track_producer, track_v, tracks_to_hits);

  // Collect showers and map showers->hits
  lar_pandora::ShowerVector shower_v;
  lar_pandora::ShowersToHits showers_to_hits;
  lar_pandora::LArPandoraHelper::CollectShowers(e, _shower_producer, shower_v, showers_to_hits);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);

  for (size_t i = 0; i < tpcobj_v.size(); i++) {

    if (_debug) std::cout << "[CandidateConsistency] >>>>> TPCObject " << i << std::endl;

    bool consistency_failed = false;
    double dqds_average = -999;

    art::Ptr<ubana::TPCObject> tpcobj = tpcobj_v.at(i);

    std::vector<art::Ptr<recob::Track>>      tracks  = tpcobjToTrackAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::PFParticle>> pfps    = tpcobjToPFPAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::Shower>>     showers = tpcobjToShowerAssns.at(tpcobj.key());

    recob::Vertex vertex = tpcobj->GetVertex();
    double xyz[3];
    vertex.XYZ(xyz);
    TVector3 vtx_xyz(xyz[0], xyz[1], xyz[2]);


    //
    // Check if 2 PFP and 1 track case (one pfp is always the neutrino pfs)
    //
    bool _1pfp_1track_case = false;
    if (pfps.size() == 2 && tracks.size() == 1) 
      _1pfp_1track_case = true;

    bool _1pfp_1track_case1 = false;
    bool _1pfp_1track_case2 = false;
    bool _1pfp_1track_case3 = false;


    //
    // Collect objects that are close to the vertex
    //
    std::vector<art::Ptr<recob::Track>>  selected_tracks;
    std::vector<geo::Point_t>            selected_tracks_vertex;
    std::vector<art::Ptr<recob::Shower>> selected_showers;
    std::vector<geo::Point_t>            selected_showers_vertex;

    // 1) Tracks
    for (auto t : tracks) {

      if ( (t->Vertex<TVector3>() - vtx_xyz).Mag() < _tolerance ) {
        selected_tracks.push_back(t);
        selected_tracks_vertex.push_back(t->Vertex());
      }

      if ( (t->End<TVector3>() - vtx_xyz).Mag() < _tolerance ) {
        selected_tracks.push_back(t);
        selected_tracks_vertex.push_back(t->End());
      }


    }

    // 2) Showers
    for (auto s : showers) {

      if ( (s->ShowerStart() - vtx_xyz).Mag() < _tolerance ) {
        selected_showers.push_back(s);
        selected_showers_vertex.push_back(geo::vect::toPoint(s->ShowerStart()));
      }

    }



    //
    // Now analyse the selected objects
    //

    // 1) Tracks 

    if(_debug) std::cout << "[CandidateConsistency] Working of Tracks" << std::endl;


    for (size_t i = 0; i < selected_tracks.size(); i++) {

      if(_debug) std::cout << "[CandidateConsistency] Track number " << i << std::endl;

      auto iter = tracks_to_hits.find(selected_tracks.at(i));
      if (iter == tracks_to_hits.end()) {
        std::cout << "[CandidateConsistency] Track in TPCObject not found by Pandora?!" << std::endl;
        continue;
      }
      auto hits = iter->second;

      // Take only collection plane hits
      std::vector<cosmictag::SimpleHit> simple_hit_v;
      for (auto h : hits) {

        if (h->View() != 2) continue;
        
        cosmictag::SimpleHit sh;
        sh.t = detProp.ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,2));
        sh.w = h->WireID().Wire * _channelMap->Plane(geo::PlaneID(0,0,2)).WirePitch();

        sh.plane = h->View();
        sh.integral = h->Integral();

        sh.time = h->PeakTime() / 4;
        sh.wire = h->WireID().Wire;

        simple_hit_v.emplace_back(sh);

      }

      // Emplacing simple hits to the manager
      cosmictag::SimpleCluster sc(simple_hit_v);
      _ct_manager.Emplace(std::move(sc));

      // Creating an approximate start hit given the TPCObject vertex
      auto& vertex = selected_tracks_vertex.at(i);
      this->ContainPoint(vertex);
      geo::PlaneID const plane_2{0, 0, 2};
      double vertex_t = detProp.ConvertXToTicks(vertex.X(), plane_2)/4.;
      int vertex_w    = _channelMap->Plane(plane_2).NearestWireID(vertex).Wire;

      cosmictag::SimpleHit start;
      start.time = vertex_t;
      start.wire = vertex_w;
      start.plane = 2;
      _ct_manager.SetStartHit(std::move(start));

      // Running the cluster analyser
      bool passed = _ct_manager.Run();

      if (passed) {
        _ct_manager.PrintClusterStatus();
      } else {
        std::cout << "not passed" << std::endl;
        continue;
      }

      cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();
      std::vector<double> dqds_v = processed_cluster._dqds_slider;

      // Calulated the trunc median
      double dqds_trunc = UBXSecHelper::GetDqDxTruncatedMean(dqds_v);
      if(_debug) std::cout << "[CandidateConsistency] dqds_trunc is " << dqds_trunc << std::endl;

      if(_1pfp_1track_case) {

        // Flash hypo test
        double x_offset = 0.;
        bool is_cathode_crossing = this->IsCathodeCrossing(selected_tracks.at(i), x_offset);

        if (is_cathode_crossing) {

          if (_debug) std::cout << "Is cathode crossing, x_offset set to: " << x_offset << std::endl;

          std::vector<double> pe_v = this->GetFlashHypo(selected_tracks.at(i), x_offset);

          double _pe_threshold = 8.;
          bool below_threshold = true;
          for (auto pe : pe_v) {
            if (pe > _pe_threshold) {
              below_threshold = false;
              break;
            }
          }
          if (below_threshold) {
            if (_debug) std::cout << "Is cathode crossing and produces flash below threshold." << std::endl;
            _1pfp_1track_case3 = true;
          }
        } else {
          if (_debug) std::cout << "Is not cathode crossing" << std::endl;
          _1pfp_1track_case3 = false;
        }


        std::vector<double> linearity_v = processed_cluster._linearity_v;

        bool is_linear = true;
        for (size_t l = 2; l < linearity_v.size()-2; l++) {
          if (linearity_v.at(l) < _linearity_threshold_track){
            if(_debug) std::cout << "[CandidateConsistency] linearity below threshold at hit n " << l << std::endl;
            is_linear = false;
          }
        }

        // Calculate dqds average of first 5 hits
        double n_hits = 5;
        dqds_average = 0.;
        for (int i = 0; i < n_hits; i++) {
          dqds_average += processed_cluster._dqds_v.at(i);
        }
        dqds_average /= n_hits;
        if(_debug) std::cout << "[CandidateConsistency] dqds average on first " 
                             << n_hits << " is: " << dqds_average << std::endl;


        // Check start hit not close to dead region
        bool close_dead_region = false;
        int start_hit_wire = processed_cluster._s_hit_v.at(0).wire;
        for (int wire = start_hit_wire - 2; wire <= start_hit_wire + 2; wire ++) {
          auto iter = wire_to_channel.find(wire);
          if (iter != wire_to_channel.end()) {
            if (chanFilt.Status(iter->second) != _good_ch_status) {
              if(_debug) std::cout << "[CandidateConsistency] Is close to dead wire" << std::endl;
              close_dead_region = true;
            }
          }
        }



        // Bump Finder
        double ratio_1 = processed_cluster._dqds_v.at(1) / processed_cluster._dqds_v.at(0);
        double ratio_2 = processed_cluster._dqds_v.at(2) / processed_cluster._dqds_v.at(1);

        std::cout << "ratio_1: " << ratio_1 << std::endl;
        std::cout << "ratio_2: " << ratio_2 << std::endl;

        _1pfp_1track_case1 = !(ratio_1 > 1.3 && ratio_2 < 0.95 && is_linear && close_dead_region);
        _1pfp_1track_case2 = !(dqds_average > 100000 && is_linear);

      } // pfp_1track_case ends

    } // selected tracks loop



    // 1) Showers 

    if(_debug) std::cout << "[CandidateConsistency] Working of Showers" << std::endl;
    
    for (size_t i = 0; i < selected_showers.size(); i++) {

      if(_debug) std::cout << "[CandidateConsistency] Shower number " << i << std::endl;

      auto iter = showers_to_hits.find(selected_showers.at(i));
      if (iter == showers_to_hits.end()) {
        std::cout << "[CandidateConsistency] Shower in TPCObject not found by Pandora?!" << std::endl;
        continue;
      }
      auto hits = iter->second;

      // Take only collection plane hits
      std::vector<cosmictag::SimpleHit> simple_hit_v;
      int n_hits_original = 0;
      for (auto h : hits) {

        if (h->View() != 2) continue;
        
        cosmictag::SimpleHit sh;
        sh.t = detProp.ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,2));
        sh.w = h->WireID().Wire * _channelMap->Plane(geo::PlaneID(0,0,2)).WirePitch();

        sh.plane = h->View();
        sh.integral = h->Integral();

        sh.time = h->PeakTime() / 4;
        sh.wire = h->WireID().Wire;

        simple_hit_v.emplace_back(sh);
        n_hits_original ++;

      }

      // Emplacing simple hits to the manager
      cosmictag::SimpleCluster sc(simple_hit_v);
      _ct_manager.Emplace(std::move(sc));

      // Creating an approximate start hit
      auto& vertex = selected_showers_vertex.at(i);
      this->ContainPoint(vertex);
      geo::PlaneID const plane_2{0, 0, 2};
      double vertex_t = detProp.ConvertXToTicks(vertex.X(), plane_2)/4.;
      int vertex_w    = _channelMap->Plane(plane_2).NearestWireID(vertex).Wire;

      cosmictag::SimpleHit start;
      start.time = vertex_t;
      start.wire = vertex_w;
      start.plane = 2;
      _ct_manager.SetStartHit(std::move(start));

      // Running the cluster analyser
      bool passed = _ct_manager.Run();

      if (passed) {
        _ct_manager.PrintClusterStatus();
      } else {
        if(_debug) std::cout << "[CandidateConsistency] Not passed" << std::endl;
        continue;
      }


      cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();
      std::vector<double> dqds_v = processed_cluster._dqds_slider;
      std::vector<cosmictag::SimpleHit> hit_v = processed_cluster._s_hit_v;


      // Check that vertex-start-hit and start hit in cluster are close
      bool are_close = false;
      TVector3 vertex_hit (vertex_w, vertex_t, 0);
      TVector3 start_hit (hit_v.at(0).wire, hit_v.at(0).time, 0);
      if ( (vertex_hit-start_hit).Mag() < _distance_cut )
        are_close = true;
      if(_debug) std::cout << "[CandidateConsistency] Distance is " << (vertex_hit-start_hit).Mag() << std::endl;

      // Check the percentage of non clustered hits
      double n_hit_ratio = (double)hit_v.size() / (double)n_hits_original * 100.;
      bool enough_hits = true;
      if (n_hit_ratio < _perc_losen_hit_cut) {
        if(_debug) std::cout << "[CandidateConsistency] Too many losen hits" << std::endl;
        enough_hits = false;
      }

      // Calulated the trunc median
      double dqds_trunc = UBXSecHelper::GetDqDxTruncatedMean(dqds_v);
      if(_debug) std::cout << "[CandidateConsistency] dqds_trunc is " << dqds_trunc << std::endl;

      // Verify linearity
      std::vector<double> linearity_v = processed_cluster._linearity_v;
      bool linearity_failed = false;
      for (size_t l = 2; l < linearity_v.size()-2; l++) {
        if (linearity_v.at(l) < _linearity_threshold) {
          if(_debug) std::cout << "[CandidateConsistency] Linearity is below threshold: " << linearity_v.at(l) << std::endl;
          linearity_failed = true;
	}
      }

      // If we have a shower, with low dqds, and whose hits
      // are not linear, then this is likely an electron
      if (dqds_trunc < _dqds_threshold && linearity_failed && are_close && enough_hits) {
        consistency_failed = true;
        if(_debug) std::cout << "[CandidateConsistency] Consistency Failed." << std::endl;
      }


    } // selected showers loop

    double cosmicScore = 0.;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;

    if (consistency_failed) {
      cosmicScore = 1.;
      tag_id = anab::CosmicTagID_t::kGeometry_Y; // ... don't have a proper id for these
    }

    if (_1pfp_1track_case1) {
      cosmicScore = 0.5;
      tag_id = anab::CosmicTagID_t::kGeometry_Z; // ... don't have a proper id for these
    }
    if (_1pfp_1track_case2) {
      cosmicScore = 0.6;
      tag_id = anab::CosmicTagID_t::kGeometry_Z; // ... don't have a proper id for these
    }
    if (_1pfp_1track_case3) {
      cosmicScore = 0.7;
      tag_id = anab::CosmicTagID_t::kGeometry_Z; // ... don't have a proper id for these
    }
    //cosmicScore = dqds_average;

    if (_debug) std::cout << "[CandidateConsistency] Cosmic Score is " << cosmicScore << std::endl;
    cosmicTagVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id); 
    util::CreateAssn(*this, e, *cosmicTagVector, tpcobj, *assnOutCosmicTagTPCObject);

  } // tpcobject loop

  e.put(std::move(cosmicTagVector));
  e.put(std::move(assnOutCosmicTagTPCObject));

  if(_debug) std::cout << "[CandidateConsistency] Ends." << std::endl;
}

void CandidateConsistency::ContainPoint(geo::Point_t& point) {

  auto const& tpc = geo->TPC();
  constexpr double e = std::numeric_limits<double>::epsilon();
  point.SetX(std::clamp(point.X(), e, 2.*tpc.HalfWidth() - e));
  point.SetY(std::clamp(point.Y(), -tpc.HalfWidth() + e, tpc.HalfWidth() - e));
  point.SetZ(std::clamp(point.Z(), e, tpc.Length() - e));

}



bool CandidateConsistency::IsCathodeCrossing(art::Ptr<recob::Track> trk_ptr, double &x_offset) {

  x_offset = -1.;

  auto const& tpc = geo->TPC();
  double _BOTTOM = -tpc.HalfHeight() + 25.; //cm
  double _TOP    = tpc.HalfHeight() - 25.; //cm

  std::vector<TVector3> sorted_trk;

  // Sort track points first
  sorted_trk.clear();
                                                                                                            
  auto const&N = trk_ptr->NumberTrajectoryPoints();
  auto const&start = trk_ptr->LocationAtPoint(0);
  auto const&end   = trk_ptr->LocationAtPoint( N - 1 );

  // if points are ordered correctly                                                                                                                                       
  if (start.Y() > end.Y()){
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( trk_ptr->LocationAtPoint<TVector3>(i) );
  }
  
  // otherwise flip order                                                                                                                                                 
  else {
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( trk_ptr->LocationAtPoint<TVector3>( N - i - 1) );
  }

  // Case 1: track exits bottom
  if ( sorted_trk.at( sorted_trk.size() - 1).Y() < _BOTTOM ) {
    auto const& top    = sorted_trk.at(0);
    auto const& bottom = sorted_trk.at( sorted_trk.size() - 1 );

    if (top.X() > bottom.X()){
      // Track crosses 
      x_offset = tpc.HalfWidth()*2 - top.X();
      return true;
    }
  }


  // Case 2: track enters top
   if (sorted_trk.at(0).Y() > _TOP) {
    auto const& top    = sorted_trk.at(0);
    auto const& bottom = sorted_trk.at( sorted_trk.size() - 1 );

    if (top.X() < bottom.X()){
      // Track crosses 
      x_offset = tpc.HalfWidth()*2 - bottom.X();
      return true;
    }
  }

  return false;

}


std::vector<double> CandidateConsistency::GetFlashHypo(art::Ptr<recob::Track> trk_ptr, double x_offset) {

  ::geoalgo::Trajectory track_geotrj;
  track_geotrj.resize(trk_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));

  for (size_t pt_idx=0; pt_idx < trk_ptr->NumberTrajectoryPoints(); ++pt_idx) {
    auto const& pt = trk_ptr->LocationAtPoint(pt_idx);
    track_geotrj[pt_idx][0] = pt.X() + x_offset;
    track_geotrj[pt_idx][1] = pt.Y();
    track_geotrj[pt_idx][2] = pt.Z();
  }        

  flashana::QCluster_t qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(track_geotrj);

  flashana::Flash_t flashHypo;
  flashHypo.pe_v.resize(geo->NOpDets());
  ((flashana::PhotonLibHypothesis*)(_mgr.GetAlgo(flashana::kFlashHypothesis)))->FillEstimate(qcluster,flashHypo);

  if (_debug) {
    for (auto pe : flashHypo.pe_v) {
      std::cout << "pe: " << pe << std::endl;
    }
  }

  return flashHypo.pe_v;
}

DEFINE_ART_MODULE(CandidateConsistency)

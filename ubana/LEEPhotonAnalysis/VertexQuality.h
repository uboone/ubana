

#ifndef VERTEXQUALITY_H
#define VERTEXQUALITY_H

#include "ParticleAssociations.h"
#include "RecoMCMatching.h"

#include "art/Framework/Principal/Event.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TTree.h"


class VertexQuality {

  std::string ftrack_producer;
  std::string fshower_producer;
  
  RecoMCMatching const * frmcm;

  geoalgo::AABox ftpc_volume;

  TTree * fvertex_tree;
  TTree * fvertex_tree_event;

  int ftpc_volume_contained;

  double fdist;
  double fdistx;
  double fdisty;
  double fdistz;

  int ftrue_track_total;
  int ftrue_shower_total;
  int freco_track_total;
  int freco_shower_total;
  int fcorrect_track_total;
  int fcorrect_shower_total;
  
  int freco_vertex_present;

  std::vector<double> ftrack_matching_ratio_v;
  std::vector<int> ftrack_true_pdg_v;
  std::vector<int> ftrack_true_origin_v;

  std::vector<double> fshower_matching_ratio_v;
  std::vector<int> fshower_true_pdg_v;
  std::vector<int> fshower_true_origin_v;

 public:

  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const & e, int const geant_track_id);

  VertexQuality(std::string const & track_producer,
		std::string const & shower_producer,
		RecoMCMatching const & rmcm);
  
  void GetTrueObjects(art::Event const & e,
		      std::vector<size_t> & mctrack_v,
		      std::vector<size_t> & mcshower_v,
		      std::vector<size_t> & mcparticle_v);
  void GetTrueRecoObjects(art::Event const & e,
			  std::vector<size_t> & track_v,
			  std::vector<size_t> & shower_v);
  void GetTrueObjects(art::Event const & e,
		      std::vector<int> & mcparticle_v);
  double GetTrueTotal(std::vector<int> const & mcparticle_v);

  void Reset();
  void FillTree(art::Event const & e,
		TTree * tree, 
		ParticleAssociations const & pas,
		size_t const closest_index,
		geoalgo::Point_t const & true_nu_vtx,
		std::vector<size_t> const & track_v,
		std::vector<size_t> const & shower_v);
  void RunDist(art::Event const & e,
	       ParticleAssociations const & pas);

};


#endif

#ifndef RECOMCMATCHING_H
#define RECOMCMATCHING_H



#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "TTree.h"

#include "SimchInfo.h"



struct RecoMCMatch {

  int mc_type;
  size_t mc_index;
  double max;
  double total;
  double ratio;

  std::unordered_map<int, double> trkide_map;

  RecoMCMatch(){}
  RecoMCMatch(int const imc_type, size_t const imc_index, double const imax, double const itotal, double const iratio);
  
};



class RecoMCMatching {

  std::string fmatch_type;
  std::string fsimch_producer;
  std::string fhit_producer;
  std::string ftrack_producer;
  std::string fshower_producer;
  std::string frmcmassociation_producer;

  bool fconsider_mcparticles;

  SimchInfo simchi;

  std::vector<std::vector<double>> fmcq_vv;
  std::vector<RecoMCMatch> ftrack_matches;
  std::vector<RecoMCMatch> fshower_matches;
  std::vector<double> fmcshower_charge;
  std::vector<double> fmctrack_charge;
  std::vector<double> fmcparticle_charge;

  TTree * fhit_tree;
  int fhit_from_reco_track;
  double fhit_time;
  double fhit_tdc;

  TTree * fsimch_tree;
  int fsimch_tdc;

 public:
  int const fmc_type_shower;
  int const fmc_type_track;
  int const fmc_type_particle;

  RecoMCMatching();

  

  void ConfigureSimch(std::string const & match_type,
		      std::string const & simch_producer,
		      std::string const & hit_producer,
		      std::string const & track_producer,
		      std::string const & shower_producer) {

    fmatch_type = match_type;
    fsimch_producer = simch_producer;
    fhit_producer = hit_producer;
    ftrack_producer = track_producer;
    fshower_producer = shower_producer;

  }

  void Configure(std::string const & hit_producer,
		 std::string const & track_producer,
		 std::string const & shower_producer,
		 std::string const & rmcmassociation_producer) {
    
    fhit_producer = hit_producer;
    ftrack_producer = track_producer;
    fshower_producer = shower_producer;
    frmcmassociation_producer = rmcmassociation_producer;

  }

  void ConsiderMCParticles(bool const consider_mcparticles = true) {fconsider_mcparticles = consider_mcparticles;}

  void FillHitTree(art::Event const & e);
  void FillSimchTree(art::Event const & e);

  //Old simchannel method
  void Match(art::Event const & e,
	     std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  void MatchMCBTAlg(art::Event const & e,
		    std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  void MatchSimchInfo(art::Event const & e,
		      std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  void FillRatio(std::vector<RecoMCMatch> & ratio_v,
		 size_t const ev_mctrack_size,
		 size_t const ev_mcshower_size,
		 size_t const ratio_index,
		 size_t const best_mc_index,
		 double const max,
		 double const total,
		 bool const last);
  void MatchAll(art::Event const & e);
  void CoutMatches(art::Event const & e);

  //New association method
  void MatchWAssociations(art::Event const & e);
  void FillAssociationVector(std::unordered_map<int, size_t> const & tp_map,
			     std::unordered_map<int, size_t> const & sp_map,
			     std::unordered_map<int, size_t> const & mcp_map,
			     art::FindManyP<recob::Hit> const & hits_per_object,
			     art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> const & particles_per_hit,
			     std::vector<RecoMCMatch> & object_matches);
  
  std::vector<RecoMCMatch> const & GetTrackMatches() const {
    return ftrack_matches;
  }
  std::vector<RecoMCMatch> const & GetShowerMatches() const {
    return fshower_matches;
  }      
  std::vector< std::vector<double> > const & GetMatchVV() const {
    return fmcq_vv;
  }
  std::vector<double> const & GetMCShowerCharge() const {
    return fmcshower_charge;
  }
  std::vector<double> const & GetMCTrackCharge() const {
    return fmctrack_charge;
  }
  std::vector<double> const & GetMCParticleCharge() const {
    return fmcparticle_charge;
  }  

};



#endif

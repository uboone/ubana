


#ifndef SIMCHINFO_H
#define SIMCHINFO_H


#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"



class SimchInfo {

  std::vector<size_t> ftrkid_to_index;
  size_t fnum_parts;
  std::map<raw::ChannelID_t, std::map<int, std::vector<double>>> fconsidered_tick_map;
  std::vector<std::vector<double>> fmcq_vv;

  void Reset();
  void Register(std::vector<unsigned int> const & track_id_v);

 public:

  struct WireRange_t;

  void FillSimchInfo(std::vector<sim::SimChannel> const & simch_v,
		     std::vector<std::vector<unsigned int>> const & g4_trackid_v,
		     std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  
  std::vector<double> MCQ(art::Ptr<recob::Hit> const & hit) const;
  std::vector<double> MCQ(std::vector<art::Ptr<recob::Hit>> const & hit_v) const;

  std::vector<std::vector<double>> const & GetMCQVV() const {return fmcq_vv;}

};



struct SimchInfo::WireRange_t {
  unsigned int ch;
  double start, end;
  WireRange_t() {
    ch    = std::numeric_limits<unsigned int>::max();
    start = end = std::numeric_limits<double>::max();
  }
  WireRange_t(unsigned int c,double s, double e)
    { ch = c; start = s; end = e; }		
};



#endif

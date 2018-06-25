


#include "SimchInfo.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"



void SimchInfo::Reset() {

 ftrkid_to_index.clear();
 fnum_parts = 0;
 fconsidered_tick_map.clear();
 fmcq_vv.clear();

}



void SimchInfo::Register(std::vector<unsigned int> const & track_id_v) {

  unsigned int max_id = 0;
  for(auto const & id : track_id_v) if(max_id < id) max_id = id;
  if(ftrkid_to_index.size() <= max_id)
    ftrkid_to_index.resize(max_id+1, SIZE_MAX);
  
  for(auto const & id : track_id_v) {
    if(ftrkid_to_index[id] == SIZE_MAX) { 
      ftrkid_to_index[id] = fnum_parts;
    }
    else {
      std::cout << "ERROR: doubly used TrackID: " << id << "\n";
      exit(1);
    }
  }
  ++fnum_parts;

}



std::vector<double> SimchInfo::MCQ(art::Ptr<recob::Hit> const & hit) const {
  
  std::vector<double> res(fnum_parts, 0);

  auto ctm_it = fconsidered_tick_map.find(hit->Channel());
  if(ctm_it == fconsidered_tick_map.end()) return res;
  std::map<int, std::vector<double>> const & ch_info = ctm_it->second;

  detinfo::DetectorClocks const * ts = lar::providerFrom<detinfo::DetectorClocksService>();
  
  auto itlow = ch_info.lower_bound(int(ts->TPCTick2TDC(hit->PeakTime()-hit->RMS())));
  auto const itup  = ch_info.upper_bound(int(ts->TPCTick2TDC(hit->PeakTime()+hit->RMS())+1));

  while(itlow != ch_info.end() && itlow != itup) {
    std::vector<double> const & edep_info = itlow->second;
    for(size_t part_index = 0; part_index < fnum_parts; ++part_index) {
      res.at(part_index) += edep_info.at(part_index);
    }
    ++itlow;
  }

  return res;

}



std::vector<double> SimchInfo::MCQ(std::vector<art::Ptr<recob::Hit>> const & hit_v) const {

  std::vector<double> res(fnum_parts, 0);

  for(art::Ptr<recob::Hit> const & hit : hit_v) {
    std::vector<double> const tmp_res = MCQ(hit);
    for(size_t i = 0; i < res.size(); ++i) res.at(i) += tmp_res.at(i);
  }

  return res;

}



void SimchInfo::FillSimchInfo(std::vector<sim::SimChannel> const & simch_v,
			      std::vector<std::vector<unsigned int>> const & g4_trackid_v,
			      std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v) {
  
  Reset();
  for(std::vector<unsigned int> const & v : g4_trackid_v) Register(v);
  ++fnum_parts;

  detinfo::DetectorClocks const * ts = lar::providerFrom<detinfo::DetectorClocksService>();
  
  for(std::vector< art::Ptr<recob::Hit> > const & hit_v : reco_to_hit_v) {
    for(art::Ptr<recob::Hit> const & hit : hit_v) {
      auto const ctm_it = fconsidered_tick_map.emplace(hit->Channel(), std::map< int, std::vector<double> >()).first;
      int const start = ts->TPCTick2TDC(hit->PeakTime()-hit->RMS());
      int const end = ts->TPCTick2TDC(hit->PeakTime()+hit->RMS()) + 1;
      for(int i = start; i <= end; ++i) ctm_it->second.emplace(i, std::vector<double>(fnum_parts, 0));
    }
  }

  for(auto const & simch : simch_v) {
      
    auto ctm_it = fconsidered_tick_map.find(simch.Channel());
    if(ctm_it == fconsidered_tick_map.end()) continue;
    std::map<int, std::vector<double>> & ch_info = ctm_it->second;
    
    for(auto const & time_ide : simch.TDCIDEMap()) {
      
      auto chi_it = ch_info.find(time_ide.first);
      if(chi_it == ch_info.end()) continue;
      std::vector<double> & edep_info = chi_it->second;

      for(auto const & ide : time_ide.second) {
	
	size_t index = SIZE_MAX;
	if(ide.trackID >= 0 && ide.trackID < int(ftrkid_to_index.size())){
	  index = ftrkid_to_index.at(ide.trackID);
	}
	if(fnum_parts <= index) {
	  edep_info.back() += ide.numElectrons;
	}
	else {
	  edep_info.at(index) += ide.numElectrons;
	}
      }
      
    }

  }

  fmcq_vv.reserve(reco_to_hit_v.size());
  for(std::vector<art::Ptr<recob::Hit>> const & hit_v : reco_to_hit_v) fmcq_vv.push_back(MCQ(hit_v));

}

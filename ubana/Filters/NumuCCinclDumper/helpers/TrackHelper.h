#ifndef TRACKHELPER_H
#define TRACKHELPER_H

#include "lardataobj/RecoBase/Track.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"


class TrackHelper
{
public:
  explicit TrackHelper();
  ~TrackHelper() = default;

  bool getPID(std::map<std::string, float> &pid_map, 
              const art::Ptr<recob::Track> &this_track,
              const art::FindManyP<anab::ParticleID> &trackPIDAssn);

  void getRangeMomentum(float length, float& mom_proton, float& mom_muon);

private:
  trkf::TrackMomentumCalculator _trkmom;
};

#endif

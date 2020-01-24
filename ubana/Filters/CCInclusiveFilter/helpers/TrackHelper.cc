#include "TrackHelper.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

TrackHelper::TrackHelper()
{
}

bool TrackHelper::getPID(std::map<std::string, float> &pid_map,
                         const art::Ptr<recob::Track> &this_track,
                         const art::FindManyP<anab::ParticleID> &trackPIDAssn)
{
    bool ok = false;
    if (!trackPIDAssn.isValid())
    {
        std::cout << "[ParticleIDValidation] trackPIDAssn.isValid() == false. Skipping track." << std::endl;
        return false;
    }
    else
    {
        std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(this_track.key());
        if (trackPID.size() != 0)
        {
            std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
            for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
            {
                anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
                if (AlgScore.fPlaneMask.none() || AlgScore.fPlaneMask.count() > 1 || (AlgScore.fPlaneMask.count() == 1 && !(AlgScore.fPlaneMask.test(0) || AlgScore.fPlaneMask.test(1) || AlgScore.fPlaneMask.test(2))))
                {
                    std::cout << "[TrackHelper::getPID] Bad AlgScore" << std::endl;
                }
                else
                {
                    int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
                    if (planeid==2 && AlgScore.fAlgName == "Chi2")
                    {
                        if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward)
                        {
                            if (AlgScore.fAssumedPdg == 2212)
                            {
                                pid_map.insert(std::pair<std::string, float>("chi2_proton", AlgScore.fValue));
                                ok = true;
                            }
                            else if (AlgScore.fAssumedPdg == 13)
                            {
                                pid_map.insert(std::pair<std::string, float>("chi2_muon", AlgScore.fValue));
                            }
                        }
                    }
                }
            }
        }
        return ok;
    }
}

void TrackHelper::getRangeMomentum(float length, float &mom_proton, float &mom_muon)
{
    mom_proton = _trkmom.GetTrackMomentum(length, 2212);
    mom_muon = _trkmom.GetTrackMomentum(length, 13);
}
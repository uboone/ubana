////////////////////////////////////////////////////////////////////////
// Class:       PandoraAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        PandoraAnalysis_module.cc
//
// Generated at Thu Nov 16 02:35:54 2023 by Andrew Chappell using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <TTree.h>

namespace pandora {
    class PandoraAnalysis;
}


class pandora::PandoraAnalysis : public art::EDAnalyzer {
    public:
        explicit PandoraAnalysis(fhicl::ParameterSet const& p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        PandoraAnalysis(PandoraAnalysis const&) = delete;
        PandoraAnalysis(PandoraAnalysis&&) = delete;
        PandoraAnalysis& operator=(PandoraAnalysis const&) = delete;
        PandoraAnalysis& operator=(PandoraAnalysis&&) = delete;

        // Required functions.
        void analyze(art::Event const& e) override;

        // Selected optional functions.
        void beginJob() override;
        void endJob() override;
        void respondToOpenInputFile(art::FileBlock const &block) override;

        void reset();
        void fillMCInfo(art::Event const& e);

    private:

        // Declare member data here.
        TTree *fTree;
        art::RunNumber_t fRun;
        art::SubRunNumber_t fSubRun;
        art::EventNumber_t fEvent;
        int fLocalEvent;
        int fIsNC;
        int fNPi0;
        int fNPiC;
        std::vector<double> fTopologicalScore;
        std::vector<double> fNPfps;
        std::vector<int> fPfpPdgVec;
        std::vector<int> fPfpClearCosmicVec;
        std::vector<int> fPfpPrimaryVec;
        std::vector<int> fPfpNumChildrenVec;
        std::vector<int> fPfpNumClustersVec;
        std::vector<int> fPfpNumHitsVec;
        std::vector<int> fPfpNumSpacePointsVec;
        std::vector<double> fPfpTrackLengthVec;
        std::vector<double> fPfpShowerLengthVec;
        std::vector<double> fPfpStartXVec;
        std::vector<double> fPfpStartYVec;
        std::vector<double> fPfpStartZVec;
        std::vector<double> fPfpStartDirectionXVec;
        std::vector<double> fPfpStartDirectionYVec;
        std::vector<double> fPfpStartDirectionZVec;
        std::vector<double> fPfpTrackThetaVec;
        std::vector<double> fPfpTrackPhiVec;
        std::vector<double> fPfpShowerOpeningAngleVec;
        std::vector<double> fMcStartXVec;
        std::vector<double> fMcStartYVec;
        std::vector<double> fMcStartZVec;
        std::vector<double> fMcDirXVec;
        std::vector<double> fMcDirYVec;
        std::vector<double> fMcDirZVec;
        std::map<int, int> fMCToHitMap;
};


pandora::PandoraAnalysis::PandoraAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
{
    // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pandora::PandoraAnalysis::analyze(art::Event const& e)
{
    this->reset();
    this->fillMCInfo(e);
    fRun = e.id().run();
    fSubRun = e.id().subRun();
    fEvent = e.id().event();
    ++fLocalEvent;

    // generator truth
    art::Handle< std::vector<simb::MCTruth>> mcTruthHandle;
    std::vector< art::Ptr<simb::MCTruth>> mcTruthVector;
    e.getByLabel("generator", mcTruthHandle);
    if (!mcTruthHandle.isValid())
        return;
    art::fill_ptr_vector(mcTruthVector, mcTruthHandle);
    if (mcTruthVector.size() != 1)
        return;
    art::Ptr<simb::MCTruth> mcTruth{mcTruthVector.front()};
    const simb::MCNeutrino mcNeutrino{mcTruth->GetNeutrino()};
    const simb::MCParticle mcParticleNu{mcNeutrino.Nu()};
    fIsNC = mcNeutrino.CCNC() == simb::kNC;

    // geant truth
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;
    if (!e.getByLabel("largeant", mcParticleHandle))
        return;
    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);

    // flash match slice selection
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;
    art::InputTag pfpInputTag("pandora", "");
    if (!e.getByLabel(pfpInputTag, pfpHandle))
        return;
    art::fill_ptr_vector(pfpVector, pfpHandle);
    std::cout << "PFPs: " << pfpVector.size() << std::endl;

    fNPfps.emplace_back(pfpVector.size());

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfpVector, neutrinoPFPs);

    art::FindManyP<recob::Cluster> clusterPfpAssn(pfpVector, e, pfpInputTag);
    art::FindManyP<recob::SpacePoint> pfpSpacepointAssn(pfpHandle, e, pfpInputTag);
    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn(pfpHandle, e, pfpInputTag);
    // Currently just looping over all PFPs, should do this by slice and look for nu etc
    for (auto pfp : pfpVector)
    {
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata{metadataAssn.at(pfp.key())};
        if (pfpMetadata.empty())
            continue;
        const auto &properties{pfpMetadata.front()->GetPropertiesMap()};
        if (properties.find("NuScore") != properties.end())
        {
            fTopologicalScore.emplace_back(pfpMetadata.front()->GetPropertiesMap().at("NuScore"));
            fPfpClearCosmicVec.emplace_back(0);
        }
        if (properties.find("IsClearCosmic") != properties.end())
        {
            fPfpClearCosmicVec.emplace_back(pfpMetadata.front()->GetPropertiesMap().at("IsClearCosmic"));
        }
        fPfpPdgVec.emplace_back(pfp->PdgCode());
        fPfpPrimaryVec.emplace_back(pfp->IsPrimary());
        fPfpNumChildrenVec.emplace_back(pfp->NumDaughters());
        std::vector<art::Ptr<recob::Cluster>> pfpClusters{clusterPfpAssn.at(pfp.key())};
        fPfpNumClustersVec.emplace_back(pfpClusters.size());
        std::vector<art::Ptr<recob::SpacePoint>> pfpSpacepoints{pfpSpacepointAssn.at(pfp.key())};
        fPfpNumSpacePointsVec.emplace_back(pfpSpacepoints.size());
    }

    // Tracks
    art::FindManyP<recob::Track> trackPfpAssn(pfpHandle, e, pfpInputTag);
    art::Handle< std::vector<recob::Track>> trackHandle;
    e.getByLabel("pandora", trackHandle);

    art::Handle<std::vector<recob::Hit>> hitHandle;
    e.getByLabel("gaushit", hitHandle);

    if (trackHandle.isValid())
    {
        art::FindManyP<recob::Hit> trackHitAssn(trackHandle, e, pfpInputTag);

        for (auto pfp : pfpVector)
        {
            std::vector<art::Ptr<recob::Track>> pfpTracks{trackPfpAssn.at(pfp.key())};
            if ((pfp->PdgCode() == 13) && !pfpTracks.empty())
            {
                art::Ptr<recob::Track> track{pfpTracks.front()};
                fPfpTrackLengthVec.emplace_back(track->Length());
                fPfpStartXVec.emplace_back(track->Start().X());
                fPfpStartYVec.emplace_back(track->Start().Y());
                fPfpStartZVec.emplace_back(track->Start().Z());
                fPfpStartDirectionXVec.emplace_back(track->StartDirection().X());
                fPfpStartDirectionYVec.emplace_back(track->StartDirection().Y());
                fPfpStartDirectionZVec.emplace_back(track->StartDirection().Z());
                fPfpTrackThetaVec.emplace_back(track->Theta());
                fPfpTrackPhiVec.emplace_back(track->Phi());
                std::vector<art::Ptr<recob::Hit>> pfpHits{trackHitAssn.at(track.key())};
                fPfpNumHitsVec.emplace_back(pfpHits.size());

                art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcHitAssn(hitHandle, e, "gaushitTruthMatch");
                std::map<const art::Ptr<simb::MCParticle>, int> sharedHitsMap;
                for (const art::Ptr<recob::Hit> &hit : pfpHits)
                {
                    const std::vector<art::Ptr<simb::MCParticle>> &matchedMC{mcHitAssn.at(hit.key())};
                    auto data{mcHitAssn.data(hit.key())};
                    int j{0};
                    for (const art::Ptr<simb::MCParticle> &mc : matchedMC)
                    {
                        auto datum{data.at(j)};
                        ++j;
                        if (datum->isMaxIDE != 1)
                            continue;

                        if (sharedHitsMap.find(mc) == sharedHitsMap.end())
                            sharedHitsMap[mc] = 0;
                        sharedHitsMap[mc]++;
                        //sharedHitsMap[mc->TrackId()]++;
                    }
                }
                int mcBestHits{0};
                double mcX{0.}, mcY{0.}, mcZ{0.}, mcDirX{0.}, mcDirY{0.}, mcDirZ{0.};
                for (auto const& [mc, nHits] : sharedHitsMap)
                {
                    if (nHits > mcBestHits)
                    {
                        mcBestHits = nHits; 
                        mcX = mc->Position().X();
                        mcY = mc->Position().Y();
                        mcZ = mc->Position().Z();
                        const TLorentzVector &mom{mc->Momentum()};
                        const TLorentzVector dir{mom * (1. / mom.Mag())};
                        mcDirX = dir.X();
                        mcDirY = dir.Y();
                        mcDirZ = dir.Z();
                    }
                }
                fMcStartXVec.emplace_back(mcX);
                fMcStartYVec.emplace_back(mcY);
                fMcStartZVec.emplace_back(mcZ);
                fMcDirXVec.emplace_back(mcDirX);
                fMcDirYVec.emplace_back(mcDirY);
                fMcDirZVec.emplace_back(mcDirZ);
            }
        }
    }

    // Showers
    art::FindManyP<recob::Shower> showerPfpAssn(pfpHandle, e, pfpInputTag);
    art::Handle< std::vector<recob::Shower>> showerHandle;
    e.getByLabel("pandora", showerHandle);
    if (showerHandle.isValid())
    {
        art::FindManyP<recob::Hit> showerHitAssn(showerHandle, e, pfpInputTag);

        for (auto pfp : pfpVector)
        {
            std::vector<art::Ptr<recob::Shower>> pfpShowers{showerPfpAssn.at(pfp.key())};
            if ((pfp->PdgCode() == 11) && !pfpShowers.empty())
            {
                art::Ptr<recob::Shower> shower{pfpShowers.front()};
                fPfpShowerLengthVec.emplace_back(shower->Length());
                fPfpStartXVec.emplace_back(shower->ShowerStart().X());
                fPfpStartYVec.emplace_back(shower->ShowerStart().Y());
                fPfpStartZVec.emplace_back(shower->ShowerStart().Z());
                fPfpStartDirectionXVec.emplace_back(shower->Direction().X());
                fPfpStartDirectionYVec.emplace_back(shower->Direction().Y());
                fPfpStartDirectionZVec.emplace_back(shower->Direction().Z());
                fPfpShowerOpeningAngleVec.emplace_back(shower->OpenAngle());
                std::vector<art::Ptr<recob::Hit>> pfpHits{showerHitAssn.at(shower.key())};
                fPfpNumHitsVec.emplace_back(pfpHits.size());

                art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcHitAssn(hitHandle, e, "gaushitTruthMatch");
                std::map<const art::Ptr<simb::MCParticle>, int> sharedHitsMap;
                for (const art::Ptr<recob::Hit> &hit : pfpHits)
                {
                    const std::vector<art::Ptr<simb::MCParticle>> &matchedMC{mcHitAssn.at(hit.key())};
                    auto data{mcHitAssn.data(hit.key())};
                    int j{0};
                    for (const art::Ptr<simb::MCParticle> &mc : matchedMC)
                    {
                        auto datum{data.at(j)};
                        ++j;
                        if (datum->isMaxIDE != 1)
                            continue;

                        if (sharedHitsMap.find(mc) == sharedHitsMap.end())
                            sharedHitsMap[mc] = 0;
                        sharedHitsMap[mc]++;
                    }
                }
                int mcBestHits{0};
                double mcX{0.}, mcY{0.}, mcZ{0.}, mcDirX{0.}, mcDirY{0.}, mcDirZ{0.};
                for (auto const& [mc, nHits] : sharedHitsMap)
                {
                    if (nHits > mcBestHits)
                    {
                        mcBestHits = nHits; 
                        mcX = mc->Position().X();
                        mcY = mc->Position().Y();
                        mcZ = mc->Position().Z();
                        const TLorentzVector &mom{mc->Momentum()};
                        const TLorentzVector dir{mom * (1. / mom.Mag())};
                        mcDirX = dir.X();
                        mcDirY = dir.Y();
                        mcDirZ = dir.Z();
                    }
                }
                fMcStartXVec.emplace_back(mcX);
                fMcStartYVec.emplace_back(mcY);
                fMcStartZVec.emplace_back(mcZ);
                fMcDirXVec.emplace_back(mcDirX);
                fMcDirYVec.emplace_back(mcDirY);
                fMcDirZVec.emplace_back(mcDirZ);
            }
        }
    }

    std::cout << "Length: " << fPfpPdgVec.size() << std::endl;

    fTree->Fill();
}

void pandora::PandoraAnalysis::fillMCInfo(art::Event const& e)
{
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    e.getByLabel("gaushit", hitHandle);
    art::fill_ptr_vector(hitVector, hitHandle);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcHitAssn(hitHandle, e, "gaushitTruthMatch");
    
    for (const art::Ptr<recob::Hit> &hit : hitVector)
    {
        const std::vector<art::Ptr<simb::MCParticle>> &matchedMC{mcHitAssn.at(hit.key())};
        auto data{mcHitAssn.data(hit.key())};
        int i{0};
        for (const art::Ptr<simb::MCParticle> &mc : matchedMC)
        {
            auto datum{data.at(i)};
            ++i;
            if (datum->isMaxIDE != 1)
                continue;

            if (fMCToHitMap.find(mc->TrackId()) == fMCToHitMap.end())
                fMCToHitMap[mc->TrackId()] = 0;
            fMCToHitMap[mc->TrackId()]++;
        }
    }
    for (auto const& [mc, nHits] : fMCToHitMap)
        std::cout << mc << ": " << nHits << std::endl;
}

void pandora::PandoraAnalysis::reset()
{
    fPfpPdgVec.clear();
    fNPfps.clear();
    fTopologicalScore.clear();
    fPfpPrimaryVec.clear();
    fPfpNumChildrenVec.clear();
    fPfpNumClustersVec.clear();
    fPfpTrackLengthVec.clear();
    fPfpShowerLengthVec.clear();
    fPfpStartXVec.clear();
    fPfpStartYVec.clear();
    fPfpStartZVec.clear();
    fPfpStartDirectionXVec.clear();
    fPfpStartDirectionYVec.clear();
    fPfpStartDirectionZVec.clear();
    fPfpTrackThetaVec.clear();
    fPfpTrackPhiVec.clear();
    fPfpShowerOpeningAngleVec.clear();
    fPfpNumHitsVec.clear();
    fPfpNumSpacePointsVec.clear();
    fMcStartXVec.clear();
    fMcStartYVec.clear();
    fMcStartZVec.clear();
    fMcDirXVec.clear();
    fMcDirYVec.clear();
    fMcDirZVec.clear();
    fMCToHitMap.clear();
}

void pandora::PandoraAnalysis::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("pandoraOutput", "Pandora Output Tree");
    fTree->Branch("run", &fRun);
    fTree->Branch("subrun", &fSubRun);
    fTree->Branch("event", &fEvent);
    fTree->Branch("localEvent", &fLocalEvent);
    fTree->Branch("isnc", &fIsNC);
    fTree->Branch("fNPfps", &fNPfps);
    fTree->Branch("topScore", &fTopologicalScore);
    fTree->Branch("pfpPdgVec", &fPfpPdgVec);
    fTree->Branch("pfpClearCosmicVec", &fPfpClearCosmicVec);
    fTree->Branch("pfpPrimaryVec", &fPfpPrimaryVec);
    fTree->Branch("pfpNumChildrenVec", &fPfpNumChildrenVec);
    fTree->Branch("fPfpNumClustersVec", &fPfpNumClustersVec);
    fTree->Branch("fPfpTrackLengthVec", &fPfpTrackLengthVec);
    fTree->Branch("fPfpShowerLengthVec", &fPfpShowerLengthVec);
    fTree->Branch("fPfpStartXVec", &fPfpStartXVec);
    fTree->Branch("fPfpStartYVec", &fPfpStartYVec);
    fTree->Branch("fPfpStartZVec", &fPfpStartZVec);
    fTree->Branch("fPfpStartDirectionXVec", &fPfpStartDirectionXVec);
    fTree->Branch("fPfpStartDirectionYVec", &fPfpStartDirectionYVec);
    fTree->Branch("fPfpStartDirectionZVec", &fPfpStartDirectionZVec);
    fTree->Branch("fPfpTrackThetaVec", &fPfpTrackThetaVec);
    fTree->Branch("fPfpTrackPhiVec", &fPfpTrackPhiVec);
    fTree->Branch("fPfpShowerOpeningAngleVec", &fPfpShowerOpeningAngleVec);
    fTree->Branch("fPfpNumHitsVec", &fPfpNumHitsVec);
    fTree->Branch("fPfpNumSpacePointsVec", &fPfpNumSpacePointsVec);
    fTree->Branch("fMcStartXVec", &fMcStartXVec);
    fTree->Branch("fMcStartYVec", &fMcStartYVec);
    fTree->Branch("fMcStartZVec", &fMcStartZVec);
    fTree->Branch("fMcDirXVec", &fMcDirXVec);
    fTree->Branch("fMcDirYVec", &fMcDirYVec);
    fTree->Branch("fMcDirZVec", &fMcDirZVec);
}

void pandora::PandoraAnalysis::endJob()
{
}

void pandora::PandoraAnalysis::respondToOpenInputFile(art::FileBlock const &block)
{
    (void)block;
    fLocalEvent = -1;
}

DEFINE_ART_MODULE(pandora::PandoraAnalysis)

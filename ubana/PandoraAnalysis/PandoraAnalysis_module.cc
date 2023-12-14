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
        double fFlashMatchScore;
        double fTopologicalScore;
        double fNPfps;
        std::vector<int> fPfpPdgVec;
        std::vector<int> fPfpClearCosmicVec;
        std::vector<int> fPfpPrimaryVec;
        std::vector<int> fPfpNumChildrenVec;
        std::vector<int> fPfpNumClustersVec;
        std::vector<double> fPfpTrackLengthVec;
        std::vector<double> fPfpShowerLengthVec;
};


pandora::PandoraAnalysis::PandoraAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
{
    // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pandora::PandoraAnalysis::analyze(art::Event const& e)
{
    this->reset();
    fRun = e.id().run();
    fSubRun = e.id().subRun();
    fEvent = e.id().event();
    ++fLocalEvent;
    fFlashMatchScore = -std::numeric_limits<double>::max();
    fTopologicalScore = -std::numeric_limits<double>::max();

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

    fNPfps = pfpVector.size();

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfpVector, neutrinoPFPs);

    art::FindManyP<recob::Cluster> clusterPfpAssn(pfpVector, e, pfpInputTag);
    art::FindManyP<recob::Track> trackPfpAssn(pfpVector, e, pfpInputTag);
    art::FindManyP<recob::Shower> showerPfpAssn(pfpVector, e, pfpInputTag);
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
            fTopologicalScore = pfpMetadata.front()->GetPropertiesMap().at("NuScore");
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
        std::vector<art::Ptr<recob::Track>> pfpTracks{trackPfpAssn.at(pfp.key())};
        if ((pfp->PdgCode() == 13) && !pfpTracks.empty())
        {
            art::Ptr<recob::Track> track{pfpTracks.front()};
            fPfpTrackLengthVec.emplace_back(track->Length());
        }
        std::vector<art::Ptr<recob::Shower>> pfpShowers{showerPfpAssn.at(pfp.key())};
        if ((pfp->PdgCode() == 11) && !pfpTracks.empty())
        {
            art::Ptr<recob::Shower> shower{pfpShowers.front()};
            fPfpShowerLengthVec.emplace_back(shower->Length());
        }

    }

    // pandora slice selection
/*    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::InputTag sliceInputTag("pandora", "allOutcomes");
    if (!e.getByLabel(sliceInputTag, sliceHandle))
        return;

    art::fill_ptr_vector(sliceVector, sliceHandle);
    art::FindManyP<recob::PFParticle> pfpAssn(sliceHandle, e, sliceInputTag);

    art::InputTag pfpInputTagAll("pandora", "allOutcomes");
    art::Handle<std::vector<recob::PFParticle>> pfpHandleAll;
    if (!e.getByLabel(pfpInputTagAll, pfpHandleAll))
        return;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssnAll(pfpHandleAll, e, pfpInputTagAll);

    for (const art::Ptr<recob::Slice> &slice : sliceVector)
    {
        const std::vector<art::Ptr<recob::PFParticle>> slicePFPs{pfpAssn.at(slice.key())};
        for (const art::Ptr<recob::PFParticle> &pfp : slicePFPs)
        {
            if (!pfp->IsPrimary())
                continue;
            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMeta{metadataAssnAll.at(pfp.key())};
            if (pfpMeta.empty())
                continue;
            const larpandoraobj::PFParticleMetadata::PropertiesMap &pfpProperties{pfpMeta.front()->GetPropertiesMap()};
            if (!pfpProperties.empty() && (pfpProperties.find("NuScore") != pfpProperties.end()))
            {
                const double thisScore{pfpProperties.at("NuScore")};
                if (thisScore > fTopologicalScore)
                {
                    fTopologicalScore = thisScore;
                }
            }
        }
    }*/

    std::cout << "Length: " << fPfpPdgVec.size() << std::endl;

    fTree->Fill();
}

void pandora::PandoraAnalysis::reset()
{
    fPfpPdgVec.clear();
    fPfpPrimaryVec.clear();
    fPfpNumChildrenVec.clear();
    fPfpNumClustersVec.clear();
    fPfpTrackLengthVec.clear();
    fPfpShowerLengthVec.clear();
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
    fTree->Branch("flashScore", &fFlashMatchScore);
    fTree->Branch("topScore", &fTopologicalScore);
    fTree->Branch("nPfps", &fNPfps);
    fTree->Branch("pfpPdgVec", &fPfpPdgVec);
    fTree->Branch("pfpClearCosmicVec", &fPfpClearCosmicVec);
    fTree->Branch("pfpPrimaryVec", &fPfpPrimaryVec);
    fTree->Branch("pfpNumChildrenVec", &fPfpNumChildrenVec);
    fTree->Branch("fPfpNumClustersVec", &fPfpNumClustersVec);
    fTree->Branch("fPfpTrackLengthVec", &fPfpTrackLengthVec);
    fTree->Branch("fPfpShowerLengthVec", &fPfpShowerLengthVec);
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

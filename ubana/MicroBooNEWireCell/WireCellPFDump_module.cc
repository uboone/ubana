////////////////////////////////////////////////////////////////////////
// Class:       WireCellPFDump
// Plugin Type: analyzer (art v3_01_02)
// File:        WireCellPFDump_module.cc
//
// 2021-02-23 modified by Haiwang Yu (hyu@bnl.gov)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "ubobj/WcpPort/NuSelectionContainment.h"
#include "ubobj/WcpPort/NuSelectionMatch.h"
#include "ubobj/WcpPort/NuSelectionTruth.h"
#include "ubobj/WcpPort/NuSelectionCharge.h"
#include "ubobj/WcpPort/NuSelectionSTM.h"
#include "ubobj/WcpPort/NuSelectionBDT.h"
#include "ubobj/WcpPort/NuSelectionKINE.h"

#include "lardataobj/MCBase/MCShower.h"

#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>


class WireCellPFDump;


class WireCellPFDump : public art::EDAnalyzer {
    public:
        explicit WireCellPFDump(fhicl::ParameterSet const& p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        WireCellPFDump(WireCellPFDump const&) = delete;
        WireCellPFDump(WireCellPFDump&&) = delete;
        WireCellPFDump& operator=(WireCellPFDump const&) = delete;
        WireCellPFDump& operator=(WireCellPFDump&&) = delete;

        // Required functions.
        void analyze(art::Event const& e) override;

        // Selected optional functions.
        void endSubRun(art::SubRun const& sr) override;

        // user defined
        void reconfigure(fhicl::ParameterSet const& pset);
        void initOutput();
        void resetOutput();

    private:

        // Declare member data here.

        // fcl config
        bool f_PFDump;
        bool fMC;
        std::string fFileType;
        std::string fPFInputTag; // inputTag -- label:instance:process
        std::string fPFtruthInputTag; // inputTag -- label:instance:process

        // output
        TTree* fTreePFDump;

        enum LIMITS {
            MAX_TRACKS = 30000,
        };

        Int_t           f_run;
        Int_t           f_subRun;
        Int_t           f_event;
        int truth_Ntrack;  // number of tracks in MC
        int truth_id[MAX_TRACKS];  // track id; size == truth_Ntrack
        int truth_pdg[MAX_TRACKS];  // track particle pdg; size == truth_Ntrack
        int truth_process[MAX_TRACKS];  // track generation process code; size == truth_Ntrack
        int truth_mother[MAX_TRACKS];  // mother id of this track; size == truth_Ntrack
        float truth_startXYZT[MAX_TRACKS][4];  // start position of this track; size == truth_Ntrack
        float truth_endXYZT[MAX_TRACKS][4];  // end position of this track; size == truth_Ntrack
        float truth_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == truth_Ntrack
        float truth_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == truth_Ntrack
        std::vector<std::vector<int> > *truth_daughters;  // daughters id of this track; vector

        int reco_Ntrack;  // number of tracks in MC
        int reco_id[MAX_TRACKS];  // track id; size == reco_Ntrack
        int reco_pdg[MAX_TRACKS];  // track particle pdg; size == reco_Ntrack
        int reco_process[MAX_TRACKS];  // track generation process code; size == reco_Ntrack
        int reco_mother[MAX_TRACKS];  // mother id of this track; size == reco_Ntrack
        float reco_startXYZT[MAX_TRACKS][4];  // start position of this track; size == reco_Ntrack
        float reco_endXYZT[MAX_TRACKS][4];  // end position of this track; size == reco_Ntrack
        float reco_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == reco_Ntrack
        float reco_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == reco_Ntrack
        std::vector<std::vector<int> > *reco_daughters;  // daughters id of this track; vector
};


WireCellPFDump::WireCellPFDump(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // ,
    // More initializers here.
{
    // Call appropriate consumes<>() for any products to be retrieved by this module.

    // fcl config
    reconfigure(p);

    initOutput();

}

void WireCellPFDump::reconfigure(fhicl::ParameterSet const& pset)
{
    std::cout<<"------------ WireCellPFDump::reconfigure ----------"<<std::endl;  
    fMC = pset.get<bool>("MC"); // overlay and full mc
    f_PFDump = pset.get<bool>("PFDump", true);

    fFileType = pset.get<std::string>("FileType", "empty");
    fPFInputTag = pset.get<std::string>("PF_inputtag");
    fPFtruthInputTag = pset.get<std::string>("PFtruth_inputtag");
}

void WireCellPFDump::initOutput()
{
    std::cout<<"------------ WireCellPFDump::initOutput ----------"<<std::endl;  

    art::ServiceHandle<art::TFileService> tfs;
    // PFDump
    if(f_PFDump) {
        fTreePFDump = tfs->make<TTree>("T_PFDump", "T_PFDump");
        fTreePFDump->Branch("run", &f_run);
        fTreePFDump->Branch("subrun", &f_subRun);
        fTreePFDump->Branch("event", &f_event);
        fTreePFDump->Branch("file_type", &fFileType);
        fTreePFDump->Branch("truth_Ntrack", &truth_Ntrack);
        fTreePFDump->Branch("truth_id", &truth_id, "truth_id[truth_Ntrack]/I");
        fTreePFDump->Branch("truth_pdg", &truth_pdg, "truth_pdg[truth_Ntrack]/I");
        fTreePFDump->Branch("truth_process", &truth_process, "truth_process[truth_Ntrack]/I");
        fTreePFDump->Branch("truth_mother", &truth_mother, "truth_mother[truth_Ntrack]/I");
        fTreePFDump->Branch("truth_startXYZT", &truth_startXYZT, "truth_startXYZT[truth_Ntrack][4]/F");
        fTreePFDump->Branch("truth_endXYZT", &truth_endXYZT, "truth_endXYZT[truth_Ntrack][4]/F");
        fTreePFDump->Branch("truth_startMomentum", &truth_startMomentum, "truth_startMomentum[truth_Ntrack][4]/F");
        fTreePFDump->Branch("truth_endMomentum", &truth_endMomentum, "truth_endMomentum[truth_Ntrack][4]/F");
        fTreePFDump->Branch("truth_daughters", &truth_daughters);

        fTreePFDump->Branch("reco_Ntrack", &reco_Ntrack);
        fTreePFDump->Branch("reco_id", &reco_id, "reco_id[reco_Ntrack]/I");
        fTreePFDump->Branch("reco_pdg", &reco_pdg, "reco_pdg[reco_Ntrack]/I");
        fTreePFDump->Branch("reco_process", &reco_process, "reco_process[reco_Ntrack]/I");
        fTreePFDump->Branch("reco_mother", &reco_mother, "reco_mother[reco_Ntrack]/I");
        fTreePFDump->Branch("reco_startXYZT", &reco_startXYZT, "reco_startXYZT[reco_Ntrack][4]/F");
        fTreePFDump->Branch("reco_endXYZT", &reco_endXYZT, "reco_endXYZT[reco_Ntrack][4]/F");
        fTreePFDump->Branch("reco_startMomentum", &reco_startMomentum, "reco_startMomentum[reco_Ntrack][4]/F");
        fTreePFDump->Branch("reco_endMomentum", &reco_endMomentum, "reco_endMomentum[reco_Ntrack][4]/F");
        fTreePFDump->Branch("reco_daughters", &reco_daughters);
    }

}

void WireCellPFDump::analyze(art::Event const& e)
{
    // Implementation of required member function here.

    // reset output at the beginning of each event
    resetOutput();

    // read NuMetrics
    std::cout<<" RUN: "<<e.run()<<"\n SUBRUN: "<<e.subRun()<<"\n EVENT: "<<e.id().event()<<std::endl;
    f_run = e.run();
    f_subRun = e.subRun();
    f_event = e.id().event();

    // reco start [nested loop]
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    if (! e.getByLabel(fPFInputTag, particleHandle)) return;
    std::vector< art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);
    std::cout << "particles.size(): " << particles.size() << std::endl;

    if(f_PFDump) {
        reco_daughters->resize(particles.size());
    }
    for (auto const& particle: particles){

        if(f_PFDump) {
            ++reco_Ntrack;
            reco_id[reco_Ntrack-1] = particle->TrackId();
            reco_pdg[reco_Ntrack-1] = particle->PdgCode();
            reco_process[reco_Ntrack-1] = 0;
            reco_mother[reco_Ntrack-1] = particle->Mother();
            auto start_pos = particle->Position();
            reco_startXYZT[reco_Ntrack-1][0] = start_pos.X();
            reco_startXYZT[reco_Ntrack-1][1] = start_pos.Y();
            reco_startXYZT[reco_Ntrack-1][2] = start_pos.Z();
            reco_startXYZT[reco_Ntrack-1][3] = start_pos.T();
            auto end_pos = particle->EndPosition();
            reco_endXYZT[reco_Ntrack-1][0] = end_pos.X();
            reco_endXYZT[reco_Ntrack-1][1] = end_pos.Y();
            reco_endXYZT[reco_Ntrack-1][2] = end_pos.Z();
            reco_endXYZT[reco_Ntrack-1][3] = end_pos.T();
            auto start_mom = particle->Momentum();
            reco_startMomentum[reco_Ntrack-1][0] = start_mom.Px();
            reco_startMomentum[reco_Ntrack-1][1] = start_mom.Py();
            reco_startMomentum[reco_Ntrack-1][2] = start_mom.Pz();
            reco_startMomentum[reco_Ntrack-1][3] = start_mom.E();
            auto end_mom = particle->EndMomentum();
            reco_endMomentum[reco_Ntrack-1][0] = end_mom.Px();
            reco_endMomentum[reco_Ntrack-1][1] = end_mom.Py();
            reco_endMomentum[reco_Ntrack-1][2] = end_mom.Pz();
            reco_endMomentum[reco_Ntrack-1][3] = end_mom.E();
            for (int i=0; i<particle->NumberDaughters(); ++i) {
                reco_daughters->at(reco_Ntrack-1).push_back(particle->Daughter(i));
            }
        }
    }

    if(fMC == true){
        /// truth start
        art::Handle< std::vector<simb::MCParticle> > particleHandle2;
        if (! e.getByLabel(fPFtruthInputTag, particleHandle2)) return;
        std::vector< art::Ptr<simb::MCParticle> > particles2;
        art::fill_ptr_vector(particles2, particleHandle2);
        // Get space charge correction
        auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

        if(f_PFDump) {
            truth_daughters->resize(particles2.size());
        }
        for (auto const& particle: particles2){

            if(f_PFDump) {
                ++truth_Ntrack;
                truth_id[truth_Ntrack-1] = particle->TrackId();
                truth_pdg[truth_Ntrack-1] = particle->PdgCode();
                truth_process[truth_Ntrack-1] = 0;
                truth_mother[truth_Ntrack-1] = particle->Mother();
                auto start_pos = particle->Position();
                auto start_sce_offset = SCE->GetPosOffsets(geo::Point_t(start_pos.X(), start_pos.Y(), start_pos.Z()));
                truth_startXYZT[truth_Ntrack-1][0] = start_pos.X() - start_sce_offset.X();
                truth_startXYZT[truth_Ntrack-1][1] = start_pos.Y() + start_sce_offset.Y();
                truth_startXYZT[truth_Ntrack-1][2] = start_pos.Z() + start_sce_offset.Z();
                truth_startXYZT[truth_Ntrack-1][3] = start_pos.T();
                truth_startXYZT[truth_Ntrack-1][0] = (truth_startXYZT[truth_Ntrack-1][0] + 0.6)*1.101/1.098 + start_pos.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
                auto end_pos = particle->EndPosition();
                auto end_sce_offset = SCE->GetPosOffsets(geo::Point_t(end_pos.X(), end_pos.Y(), end_pos.Z()));
                truth_endXYZT[truth_Ntrack-1][0] = end_pos.X() - end_sce_offset.X();
                truth_endXYZT[truth_Ntrack-1][1] = end_pos.Y() + end_sce_offset.Y();
                truth_endXYZT[truth_Ntrack-1][2] = end_pos.Z() + end_sce_offset.Z();
                truth_endXYZT[truth_Ntrack-1][3] = end_pos.T();
                truth_endXYZT[truth_Ntrack-1][0] = (truth_endXYZT[truth_Ntrack-1][0] + 0.6)*1.101/1.098 + end_pos.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
                auto start_mom = particle->Momentum();
                truth_startMomentum[truth_Ntrack-1][0] = start_mom.Px();
                truth_startMomentum[truth_Ntrack-1][1] = start_mom.Py();
                truth_startMomentum[truth_Ntrack-1][2] = start_mom.Pz();
                truth_startMomentum[truth_Ntrack-1][3] = start_mom.E();
                auto end_mom = particle->EndMomentum();
                truth_endMomentum[truth_Ntrack-1][0] = end_mom.Px();
                truth_endMomentum[truth_Ntrack-1][1] = end_mom.Py();
                truth_endMomentum[truth_Ntrack-1][2] = end_mom.Pz();
                truth_endMomentum[truth_Ntrack-1][3] = end_mom.E();
                for (int i=0; i<particle->NumberDaughters(); ++i) {
                    truth_daughters->at(truth_Ntrack-1).push_back(particle->Daughter(i));
                }
            }

        }
    }

    if(f_PFDump) fTreePFDump->Fill();
}

void WireCellPFDump::endSubRun(art::SubRun const& sr)
{
    // Implementation of optional member function here.
}

void WireCellPFDump::resetOutput()
{
    // live period within each event
    // maybe redundant here
    if(f_PFDump){
        truth_Ntrack = 0;
        truth_daughters->clear();
        reco_Ntrack = 0;
        reco_daughters->clear();
    }
}

DEFINE_ART_MODULE(WireCellPFDump)

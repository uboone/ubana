////////////////////////////////////////////////////////////////////////
// Class:       SingleShowerFilter
// Plugin Type: filter (art v2_05_01)
// File:        SingleShowerFilter_module.cc
//
////////////////////////////////////////////////////////////////////////

// Base includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "Pandora/PdgTable.h"
#include <string>
#include <map>


#include "TTree.h"

#include <memory>

class SingleShowerFilter;


class SingleShowerFilter : public art::EDFilter {
    public:
        explicit SingleShowerFilter(fhicl::ParameterSet const & p);

        SingleShowerFilter(SingleShowerFilter const &) = delete;
        SingleShowerFilter(SingleShowerFilter &&) = delete;
        SingleShowerFilter & operator = (SingleShowerFilter const &) = delete;
        SingleShowerFilter & operator = (SingleShowerFilter &&) = delete;

        bool filter(art::Event & e) override;
        void endJob() override;

        void quickFinalStateParticles(const std::map< size_t, art::Ptr<recob::PFParticle>> &pfParticleMap, std::vector< art::Ptr<recob::PFParticle> > &crParticles, std::vector< art::Ptr<recob::PFParticle> > &nuParticles );


    private:
        std::string m_pfp_producer;
        double m_ShowerCut;
        int m_NumShowerCut;
};


SingleShowerFilter::SingleShowerFilter(fhicl::ParameterSet const & pset) :
    EDFilter(pset)
{
    //producer datalabels
    m_pfp_producer = pset.get<std::string>("PFP_producer", "pandora");
    //cut variables
    m_ShowerCut = pset.get<double>("ShowerScoreCut", 0.5);
    m_NumShowerCut = pset.get<double>("NumShowerCut", 1);
}

bool SingleShowerFilter::filter(art::Event & evt)
{

    //Collect the PFParticles from the event. This is the core!
    art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
    std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
    art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

    //This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;
    // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
    std::map< size_t, art::Ptr<recob::PFParticle>> pfParticleMap;
    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
        {
            const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
            if (!pfParticleMap.insert(std::map< size_t, art::Ptr<recob::PFParticle>>::value_type(pParticle->Self(), pParticle)).second)
            {
                throw cet::exception("SingleShowerFilter") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
            }
        }


    //--------------------------------
    // Produce two PFParticle vectors containing final-state particles:
    // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
    // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
    std::vector< art::Ptr<recob::PFParticle> > crParticles;
    std::vector< art::Ptr<recob::PFParticle> > nuParticles;
    quickFinalStateParticles(pfParticleMap, crParticles, nuParticles);
    
    //add the associaton between PFP and metadata, this is important to look at the slices and scores
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  m_pfp_producer);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
        pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
    }

    int num_showers = 0;
    for(size_t p=0; p< nuParticles.size();p++){
        auto pfp = nuParticles[p];
        float score = 0;

        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatalist = pfParticleToMetadataMap[pfp];

        //will be empty in circumstances outlined above, not every PFP has stored metadata
        if (!metadatalist.empty()){

            //std::cout<<"metadatalist not empty for pfp with index and pdg code: "<<pfp->Self()<<"/"<<pfp->PdgCode()<<", primary = "<<pfp->IsPrimary()<<std::endl;

            for(art::Ptr<larpandoraobj::PFParticleMetadata> data:metadatalist){

                //get the metadata properties
                std::map<std::string, float> propertiesmap  = data->GetPropertiesMap();
                //std::cout<<"the number of items in the metadata properties map is "<<propertiesmap.size()<<std::endl;
                //
                //
                if( propertiesmap.count("TrackScore")>0){
                    score  = propertiesmap["TrackScore"];
                }else{ 
                    score = 1;

                }
            }
        }
        
        if(score <= m_ShowerCut) num_showers++;
        if(num_showers>m_NumShowerCut) return false;    
    }

    if(num_showers==m_NumShowerCut){
        return true;
    }else{
        return false;
    }

}



    void SingleShowerFilter::quickFinalStateParticles(const std::map< size_t, art::Ptr<recob::PFParticle>> &pfParticleMap, std::vector< art::Ptr<recob::PFParticle> > &crParticles, std::vector< art::Ptr<recob::PFParticle> > &nuParticles )
    {

        int found = 0;
        int primaries = 0;
        int full = 0;
        for (std::map< size_t, art::Ptr<recob::PFParticle>>::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
        {
            const art::Ptr<recob::PFParticle> pParticle(it->second);

            full++;
            // Only look for primary particles
            if (!pParticle->IsPrimary()) continue;

            // Check if this particle is identified as the neutrino
            const int pdg(pParticle->PdgCode());
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);


            primaries++;
            if(isNeutrino){
                found++;
            }

            // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
            if (!isNeutrino)
            {
                crParticles.push_back(pParticle);
                continue;
            }

            // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
            //       If this is not the case please handle accordingly
            if (!nuParticles.empty())
            {
                throw cet::exception("SinglePhoton") << "  This event contains multiple reconstructed neutrinos!";
            }

            // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
            for (const size_t daughterId : pParticle->Daughters())
            {
                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";

                nuParticles.push_back(pfParticleMap.at(daughterId));
            }
        }


    }




void SingleShowerFilter::endJob(){
    // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SingleShowerFilter)

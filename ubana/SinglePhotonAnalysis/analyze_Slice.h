#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearSlices(){
        m_reco_slice_num = 0;
        m_reco_slice_nuscore.clear();
    }

    void SinglePhoton::ResizeSlices(size_t size){
        m_reco_slice_nuscore.resize(size);
    }


    void SinglePhoton::CreateSliceBranches(){
        vertex_tree->Branch("reco_slice_nuscore",&m_reco_slice_nuscore);
        vertex_tree->Branch("reco_slice_num",&m_reco_slice_num);
    }

    void SinglePhoton::AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap, std::vector<art::Ptr<recob::PFParticle>>  &pfParticleVector ){
        std::vector<double> nuscore_slices;



        /*
         * Grabbed this info from Giuseppe:
         * Here's the structure of these Metadata
         * Primary PfParticles are either
         * 1) IsClearCosmic = 1 (for unambiguous cosmics)
         * 2) NuScore = 0.108586, SliceIndex = 1 (for cosmic slices)
         * 3) IsNeutrino = 1, NuScore = 0.170914, SliceIndex = 2 (for the nu slice)
         * Then, for PfParticles that are daughter of the nu, the track score is saved, e.g.:
         * 4) TrackScore = 0.671488
         * PfParticles that are not primary and that are not daugthers of the neutrino have empty Metadata
         */


        //print out some stuff about the metadata
        for (auto pair: pfParticleToMetadataMap){
            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatalist= pair.second;
            art::Ptr<recob::PFParticle> pfp = pair.first;

            //will be empty in circumstances outlined above
            if (!metadatalist.empty()){

                std::cout<<"metadatalist not empty for pfp with index and pdg code: "<<pfp->Self()<<"/"<<pfp->PdgCode()<<", primary = "<<pfp->IsPrimary()<<std::endl;

                for(art::Ptr<larpandoraobj::PFParticleMetadata> data:metadatalist){

                    //const pandora::PropertiesMap &pfParticlePropertiesMap(metadata->GetPropertiesMap()); 
                    std::map<std::string, float> propertiesmap  = data->GetPropertiesMap();
                    std::cout<<"the number of items in the metadata properties map is "<<propertiesmap.size()<<std::endl;
                    for (auto it:propertiesmap ){
                        std::cout << "  - " << it.first << " = " << it.second << std::endl;
                        if (it.first == "NuScore"){
                            nuscore_slices.push_back(it.second);
                        }
                    }
                }

            }

        }

        m_reco_slice_num = nuscore_slices.size();
        this->ResizeSlices(m_reco_slice_num); 
        m_reco_slice_nuscore = nuscore_slices;

    }
}

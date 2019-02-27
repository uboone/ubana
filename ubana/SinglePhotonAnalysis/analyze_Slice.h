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

    //called once per event to get all the slice info
    //fills a map between the neutrino score for the slice and the primary reco PFP
    void SinglePhoton::AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap, std::vector<art::Ptr<recob::PFParticle>>  &pfParticleVector){
        std::vector<double> nuscore_slices; //this is a temporary vector to store neutrino score per slice for this event
        std::vector<art::Ptr<recob::PFParticle>> primary_pfps; //store the primary PFP for each slice

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


        //for all PFP metadata in the event
        for (auto pair: pfParticleToMetadataMap){
            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatalist= pair.second; //get the metadata
            art::Ptr<recob::PFParticle> pfp = pair.first; //get the corresponding PFP

            //will be empty in circumstances outlined above, not every PFP has stored metadata
            if (!metadatalist.empty()){

                //std::cout<<"metadatalist not empty for pfp with index and pdg code: "<<pfp->Self()<<"/"<<pfp->PdgCode()<<", primary = "<<pfp->IsPrimary()<<std::endl;

                //for each PFP per event
                for(art::Ptr<larpandoraobj::PFParticleMetadata> data:metadatalist){

                    //const pandora::PropertiesMap &pfParticlePropertiesMap(metadata->GetPropertiesMap()); 

                    //get the metadata properties
                    std::map<std::string, float> propertiesmap  = data->GetPropertiesMap();
                    //std::cout<<"the number of items in the metadata properties map is "<<propertiesmap.size()<<std::endl;

                    //for each of the things in the list
                    for (auto it:propertiesmap ){
                        //std::cout << "  - " << it.first << " = " << it.second << std::endl;

                        //store the neutrino score for each slice
                        if (it.first == "NuScore"){
                            nuscore_slices.push_back(it.second);

                            //if it's the neutrino score it also means it's the primary particle so save the pfp
                            primary_pfps.push_back(pfp);
                            pfParticleToNuScoreMap[pfp] = it.second;
                        }
                    }//for each item in properties map
                }//for each PFP/metadata

            }//if the list isn't empty

        }//for all PFP/metadata in event

       //store stuff in the output tree
        m_reco_slice_num = nuscore_slices.size();//the number of slices also corresponds to the number of neutrino scores

        std::cout<<"saving the info for "<<m_reco_slice_num<<" slices"<<std::endl;
        this->ResizeSlices(m_reco_slice_num); 
        m_reco_slice_nuscore = nuscore_slices;

    }

    /*
      //here we put in some reco truth matching thing where given an true interaction, find the corresponding reco objects
        //then for the given recob::shower/track(s), match it back to the primary pfp for slice
        //can do slice score, completeness, etc. as function of shower energy, conversion length, etc. 
    void SinglePhoton::GetSliceMatchInteraction(std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap, std::vector<art::Ptr<recob::Shower>>& this_shower, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>>& showerToPFParticleMap){
        //for the moment just doing a test with a single shower, will expand to more particles
        //for each shower, get the corresponding PFP
            //look up the PFP chain to get the primary
            //store all the primaries
              

    }
    */
       
}

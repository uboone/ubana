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
    //loops over all PFP's to find the primary and then associate to a slice
    void SinglePhoton::AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap,  PFParticleIdMap &pfParticleMap){
        std::vector<std::pair<art::Ptr<recob::PFParticle>, int>> primarySliceIdVec; //maps a primary PFP to a slice index
        std::map<int, double> sliceIdToNuScoreMap; //maps a slice index to the associated neutrino score
        std::vector<art::Ptr<recob::PFParticle>> clearCosmicPFP;

        std::vector<double> nuscore_slices; //this is a temporary vector to store neutrino score per slice for this event
        //std::vector<art::Ptr<recob::PFParticle>> primary_pfps; //store the primary PFP for each slice

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
                    int temp_ind = -1;
                    double temp_score = -1.0;
                    int clear_cosmic = -1;
                    //for each of the things in the list
                    for (auto it:propertiesmap ){
                        // std::cout << "  - " << it.first << " = " << it.second << std::endl;
                        if (it.first == "SliceIndex"){
                            temp_ind = it.second;
                            // std::cout << "  - " << it.first << " = " << it.second << std::endl;
                        }
                        //store the neutrino score for each slice
                        if (it.first == "NuScore"){
                            nuscore_slices.push_back(it.second);
                            temp_score = it.second;
                            //std::cout << "  - " << it.first << " = " << it.second << std::endl;
                            //if it's the neutrino score it also means it's the primary particle so save the pfp
                            //primary_pfps.push_back(pfp);
                            //pfParticleToNuScoreMap[pfp] = it.second;

                        
                     }
                         if (it.first == "IsClearCosmic"){
                             clear_cosmic = 1;
                         }
                    }//for each item in properties map

                    //if there is a neutrino score it's the primary PFP, so save the score+slice info
                    if(temp_score != -1.0){
                        primarySliceIdVec.push_back(std::pair(pfp, temp_ind));
                        //primaryToSliceIdMap[pfp] = temp_ind;
                        sliceIdToNuScoreMap[temp_ind] = temp_score;
                        std::cout<<"this primary PFP at index "<<pfp->Self()<<std::endl;

                    }
                    if( clear_cosmic> 0){
                        //std::cout<<"clear cosmic with PFP id "<<pfp->Self()<<std::endl;
                        clearCosmicPFP.push_back(pfp);
                    }

                }//for each PFP/metadata

            }//if the list isn't empty

        }//for all PFP/metadata in event

        //now we have all the primary pfp's and the corresponding slices+scores
        //the next step is to look at all the pfp's in the event, find the primary, and then store the slice ind
        std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind

        //for all pfp's in the event
        //std::cout<<"looking at all PFP's"<<std::endl;
        //for (unsigned int i = 0; i< pfParticleVector.size(); i++){
        for(auto item: pfParticleMap){
            art::Ptr<recob::PFParticle> start_pfp = item.second;
            //skipping pfp's that are clear cosmics
            auto iter = find(clearCosmicPFP.begin(), clearCosmicPFP.end(), start_pfp);
            if(iter != clearCosmicPFP.end()) continue;

            //std::cout<<"START: looking for match for pfp - track id/pdg code "<<start_pfp->Self()<<"/"<<start_pfp->PdgCode()<<std::endl; 

            art::Ptr<recob::PFParticle> this_pfp = start_pfp;
            art::Ptr<recob::PFParticle> parent_pfp ;

            //if this is a primary particle, skip next part
            if(!this_pfp->IsPrimary()){
                parent_pfp = pfParticleMap[this_pfp->Parent()];
                //std::cout<<"parent of start particle is "<<parent_pfp->Self()<<"/"<<parent_pfp->PdgCode()<<std::endl;   
                //if not primary, iterate up parent chain
                while(!parent_pfp->IsPrimary()){
                  //  std::cout<<"not primary - track id/pdg code "<<parent_pfp->Self()<<"/"<<parent_pfp->PdgCode()<<std::endl; 
                    //std::cout<<"iterating, current particle has index "<<parent_pfp->Self()<<std::endl;
                    parent_pfp = pfParticleMap[this_pfp->Parent()];
                    this_pfp = parent_pfp;
                }//while not primary, iterate up chain
            } else{
                parent_pfp = start_pfp;
            }

            //std::cout<<"this particle was primary at track id/pdg code "<<parent_pfp->Self()<<"/"<<parent_pfp->PdgCode()<<std::endl; 

            //get slice id for this primary
            int slice_id = -1;//initialize to invalid

            //for each primary pfp
            for(unsigned int j = 0; j <primarySliceIdVec.size(); j++){
                //if the parent primary matches to a primary in the list, save the slice id
                if (parent_pfp == primarySliceIdVec[j].first){
                    slice_id = primarySliceIdVec[j].second;
                    break;
                }
            }//for all primary PFP's in event

            //store original pfp and it's slice id
            if(slice_id < 0 ){
                std::cout<<"no matching slice found for this PFP with primary id "<<parent_pfp->Self()<<std::endl;
            } else {
                // allPFPSliceIdVec[i] = std::pair(start_pfp,slice_id);
                allPFPSliceIdVec.push_back(std::pair(start_pfp,slice_id));
            }
        }//for all pfp's in the event 


        /*
         * store stuff in the output tree
         */
        m_reco_slice_num = nuscore_slices.size();//the number of slices also corresponds to the number of neutrino scores

        //currently this is junk, just a placeholder
        std::cout<<"saving the info for "<<m_reco_slice_num<<" slices"<<std::endl;
        this->ResizeSlices(m_reco_slice_num); 
        m_reco_slice_nuscore = nuscore_slices;

        std::cout<<"the number of clear cosmic PFP's in the event is "<<clearCosmicPFP.size()<<std::endl;
        std::cout<<"the number of items in the primary to slice vector is "<<primarySliceIdVec.size()<<" and in the slice to score map is "<<sliceIdToNuScoreMap.size()<<std::endl;
        std::cout<<"the number of PFP's matched to a slice is "<< allPFPSliceIdVec.size()<<"/"<< pfParticleMap.size()<<std::endl;
        if ((clearCosmicPFP.size() +  allPFPSliceIdVec.size())!= pfParticleMap.size()){
            std::cout<<"BIG ERROR, UNACCOUNTED FOR PFP's, (clearCosmicPFP.size() +  allPFPSliceIdVec.size())!= pfParticleMap.size())"<<std::endl;
        }
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

#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearSlices(){
        m_reco_slice_num = 0;
        m_reco_slice_nuscore.clear();
        m_reco_slice_shower_num_matched_signal = -999;
        m_reco_slice_track_num_matched_signal = -999;
        m_reco_slice_shower_matched_sliceId.clear();
        m_reco_slice_track_matched_sliceId.clear();
        m_reco_slice_shower_matched_energy.clear();
        m_reco_slice_track_matched_energy.clear();
        m_reco_slice_shower_matched_conversion.clear();
        m_reco_slice_shower_matched_overlay_frac.clear();
    }

    //resizes the branches that are filled for every slice int the event
    void SinglePhoton::ResizeSlices(size_t size){
        m_reco_slice_nuscore.resize(size);
    }

    //resize the branches that are filled for matched track and shower objects
    void SinglePhoton::ResizeMatchedSlices(size_t size_shower ,size_t size_track){
        m_reco_slice_shower_matched_sliceId.resize(size_shower);
        m_reco_slice_track_matched_sliceId.resize( size_track);
        m_reco_slice_shower_matched_energy.resize(size_shower);
        m_reco_slice_track_matched_energy.resize( size_track);
        m_reco_slice_shower_matched_conversion.resize(size_shower);
        m_reco_slice_shower_matched_overlay_frac.resize(size_shower);

    }



    void SinglePhoton::CreateSliceBranches(){
        vertex_tree->Branch("reco_slice_nuscore",&m_reco_slice_nuscore);
        vertex_tree->Branch("reco_slice_num",&m_reco_slice_num);
        vertex_tree->Branch("reco_slice_shower_num_matched_signal",& m_reco_slice_shower_num_matched_signal);
        vertex_tree->Branch("reco_slice_track_num_matched_signal",& m_reco_slice_track_num_matched_signal);
    }

    void SinglePhoton::CreateMatchedSliceBranches(){
        vertex_tree->Branch("reco_slice_shower_matched_sliceId", & m_reco_slice_shower_matched_sliceId);
        vertex_tree->Branch("reco_slice_track_matched_sliceId", & m_reco_slice_track_matched_sliceId);
        vertex_tree->Branch("reco_slice_shower_matched_energy",& m_reco_slice_shower_matched_energy);
        vertex_tree->Branch("reco_slice_track_matched_energy",& m_reco_slice_track_matched_energy);
        vertex_tree->Branch("reco_slice_shower_matched_conversion",&  m_reco_slice_shower_matched_conversion); 
        vertex_tree->Branch("reco_slice_shower_matched_overlay_frac",&  m_reco_slice_shower_matched_overlay_frac); 
    }

    //called once per event to get all the slice info
    //fills a map between the neutrino score for the slice and the primary reco PFP
    //loops over all PFP's to find the primary and then associate to a slice
    void SinglePhoton::AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap,
            PFParticleIdMap &pfParticleMap,
            std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> &primaryPFPSliceIdVec,
            std::map<int, double> &sliceIdToNuScoreMap,
            std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
            std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap){


        //std::vector<std::pair<art::Ptr<recob::PFParticle>, int>> primaryPFPSliceIdVec; //maps a primary PFP to a slice index
        // std::map<int, double> sliceIdToNuScoreMap; //maps a slice index to the associated neutrino score
        std::vector<art::Ptr<recob::PFParticle>> clearCosmicPFP;
        std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind


        std::vector<double> nuscore_slices; //this is a temporary vector to store neutrino score per slice for this event
        //std::vector<art::Ptr<recob::PFParticle>> primary_pfps; //store the primary PFP for each slice
        // sliceIdToPFPMap.clear(); //clear between events

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
                        primaryPFPSliceIdVec.push_back(std::pair(pfp, temp_ind));
                        //primaryToSliceIdMap[pfp] = temp_ind;
                        sliceIdToNuScoreMap[temp_ind] = temp_score;
                        if(m_is_verbose)std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t found primary PFP at index "<<pfp->Self()<<std::endl;

                    }
                    if( clear_cosmic> 0){
                        //std::cout<<"clear cosmic with PFP id "<<pfp->Self()<<std::endl;
                        clearCosmicPFP.push_back(pfp);
                        PFPToClearCosmicMap[pfp] = true;
                    }else{
                        PFPToClearCosmicMap[pfp] = false;

                    }

                }//for each PFP/metadata

            }//if the list isn't empty

        }//for all PFP/metadata in event

        //now we have all the primary pfp's and the corresponding slices+scores
        //the next step is to look at all the pfp's in the event, find the primary, and then store the slice ind
        // std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind

        //for all pfp's in the event
        //std::cout<<"looking at all PFP's"<<std::endl;
        //std::map<int, std::vector<art::Ptr<recob::PFParticle>>> sliceIdToPFPMap;

        //for (unsigned int i = 0; i< pfParticleVector.size(); i++){
        for(auto item: pfParticleMap){
            art::Ptr<recob::PFParticle> start_pfp = item.second;
            //no longer skipping pfp's that are clear cosmics
            //auto iter = find(clearCosmicPFP.begin(), clearCosmicPFP.end(), start_pfp);
            //if(iter != clearCosmicPFP.end()) continue;

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
            for(unsigned int j = 0; j <primaryPFPSliceIdVec.size(); j++){
                //if the parent primary matches to a primary in the list, save the slice id
                if (parent_pfp == primaryPFPSliceIdVec[j].first){
                    slice_id = primaryPFPSliceIdVec[j].second;
                    break;
                }
            }//for all primary PFP's in event

            //store original pfp and it's slice id
            if(slice_id < 0 ){
                if(m_is_verbose)std::cout<<"no matching slice found for this PFP with primary id "<<parent_pfp->Self()<<std::endl;
            } else {
                // allPFPSliceIdVec[i] = std::pair(start_pfp,slice_id);
                allPFPSliceIdVec.push_back(std::pair(start_pfp,slice_id));
                // PFPToSliceIdMap[start_pfp] = slice_id; 
            }
            PFPToSliceIdMap[start_pfp] = slice_id; 

            // sliceIdToPFPMap[slice_id].push_back(start_pfp);
        }//for all pfp's in the event

        //for (auto pair: sliceIdToPFPMap){
        //     std::cout<<"in slice ID "<<pair.first<<" there are "<<pair.second.size()<<" PFP's"<<std::endl;
        //} 


        /*
         * store stuff in the output tree
         */
        m_reco_slice_num = nuscore_slices.size();//the number of slices also corresponds to the number of neutrino scores

        //currently this is junk, just a placeholder
        //std::cout<<"saving the info for "<<m_reco_slice_num<<" slices"<<std::endl;
        this->ResizeSlices(m_reco_slice_num); 
        m_reco_slice_nuscore = nuscore_slices;

        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of clear cosmic PFP's in the event is "<<clearCosmicPFP.size()<<std::endl;
        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of items in the primary to slice vector is "<<primaryPFPSliceIdVec.size()<<" and in the slice to score map is "<<sliceIdToNuScoreMap.size()<<std::endl;
        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of PFP's matched to a slice is "<< allPFPSliceIdVec.size()<<"/"<< pfParticleMap.size()<<std::endl;
        if ((clearCosmicPFP.size() +  allPFPSliceIdVec.size())!= pfParticleMap.size()){
            std::cout<<"BIG ERROR, UNACCOUNTED FOR PFP's, (clearCosmicPFP.size() +  allPFPSliceIdVec.size())!= pfParticleMap.size())"<<std::endl;
        }

        if (PFPToSliceIdMap.size() != pfParticleMap.size()){
            std::cout<<"BIG ERROR, PFPToSliceIdMap.size() != pfParticleMap.size(): "<<std::endl;
            std::cout<<PFPToSliceIdMap.size()<<"/"<<pfParticleMap.size()<<std::endl;
        }

        if (PFPToClearCosmicMap.size() != pfParticleToMetadataMap.size()){
            std::cout<<"BIG ERROR, PFPToClearCosmicMap.size() != pfParticleToMetadataMap.size(): "<<std::endl;
            std::cout<<PFPToClearCosmicMap.size()<<"/"<<pfParticleToMetadataMap.size()<<std::endl;

        } 

        //std::vector<std::vector<int>> pfp_pdg_slice; 
        //for(auto item: allPFPSliceIdVec){
        //std::cout<<"the pfp with id "<<item.first->Self()<<" is associated to slice "<<item.second<<std::endl;
        //      pfp_pdg_slice[item.second].push_back(item.first->PdgCode());
        // }
    }


    //here we put in some reco truth matching thing where given an true interaction, find the corresponding reco objects
    //then for the given recob::shower/track(s), match it back to the primary pfp for slice
    //can do slice score, completeness, etc. as function of shower energy, conversion lengt
    int SinglePhoton::GetShowerSlice(art::Ptr<recob::Shower>& this_shower, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>>& showerToPFParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec){
        //for the shower, get the associated PFP 
        art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[this_shower];

        int slice = -1;
        //for the pfp, get the slice
        for(auto pair: allPFPSliceIdVec){
            if(pair.first == pfp){
                slice =  pair.second;
            }
        }

        //return the slice or -1 if there isn't an associated slice - this means clear cosmic
        return slice;
    }

    int SinglePhoton::GetTrackSlice(art::Ptr<recob::Track>& this_track, std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>& trackToPFParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec){
        //for the track, get the associated PFP 
        art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[this_track];

        int slice = -1;
        //for the pfp, get the slice
        for(auto pair: allPFPSliceIdVec){
            if(pair.first == pfp){
                slice =  pair.second;
            }
        }

        //return the slice or -1 if there isn't an associated slice - this means clear cosmic
        return slice;
    }

    //for a given signal def, finds the MCParticles in event
    //loops over association between reco tracks/showers to get associated slice(s)
    //can also look at things like shower energy, conversion length, etc.
    void SinglePhoton::FindSignalSlice(std::string signal_def, std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,
            std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle> > & showerToPFParticleMap, 
            std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec, 
            std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
            std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle> > & trackToNuPFParticleMap,
            std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap){
        // std::vector<recob::Shower> shower_from_truth; //stores the recob::Showers which are matched to the signal MCP's
        // std::vector<recob::Track> track_from_truth; //stores the recob::Tracks matched to the signal MCP's

        //std::vector<simb::MCParticle>;

        std::vector<double> matched_reco_slice_shower_overlay_fraction;
        //std::vector<int> matched_reco_slice_shower_matched;
        std::vector<art::Ptr<simb::MCParticle>> matched_reco_slice_shower_MCP;
        //std::vector<double> matched_reco_slice_track_overlay_fraction;
        //std::vector<int> matched_reco_slice_track_matched;
        std::vector<art::Ptr<simb::MCParticle>> matched_reco_slice_track_MCP;


        //first check if in the event there's a match to a given signal
        if(signal_def == "ncdelta"){
            //std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t looking for signal def "<<signal_def<<", m_mctruth_is_delta_radiative = "<<m_mctruth_is_delta_radiative<<std::endl; 
            if(m_mctruth_is_delta_radiative== true){
                std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t looking for signal def "<<signal_def<<", m_mctruth_is_delta_radiative = "<<m_mctruth_is_delta_radiative<<std::endl; 
                //collect primary particles
                //first look for sim showers
                for (unsigned int j = 0; j< m_sim_shower_parent_pdg.size(); j++){
                    int parent= m_sim_shower_parent_pdg[j];
                    int pdg =  m_sim_shower_pdg[j];

                    //std::cout<<"found sim photon shower with pdg "<<pdg<<" and parent pdg "<<parent<<std::endl;

                    //if this sim shower is a photon and it's primary (parent pdg is -1)
                    if(parent == -1 && pdg ==22){
                        //first check that this particle isn't alread saved
                        //use map from track ID to get MCP
                        int id = m_sim_shower_trackID[j];
                        art::Ptr<simb::MCParticle> mcp = MCParticleToTrackIDMap[id];
                        // std::cout<<"found sim photon shower with track ID "<<id<<std::endl;
                        if (std::find(matched_reco_slice_shower_MCP.begin(),matched_reco_slice_shower_MCP.end(),mcp) ==  matched_reco_slice_shower_MCP.end()){
                            //if this shower is matched to a recob:shower
                            if (m_sim_shower_matched[j] > 0){
                                //matched_reco_slice_shower_matched.push_back(m_sim_shower_matched[j]);
                                matched_reco_slice_shower_MCP.push_back(mcp);

                                //save the overlay fraction and whether it's matched
                                matched_reco_slice_shower_overlay_fraction.push_back(m_sim_shower_overlay_fraction[j]);
                                //std::cout<<"saving sim photon shower with track ID "<<id<<std::endl;
                            }
                        } //if the particle isn't already stored

                    }//if it's a photon from the neutrino interaction
                }//for all sim showers

                if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of sim showers-MCP matches associated to the true ncdelta is "<<matched_reco_slice_shower_MCP.size()<<std::endl;

                m_reco_slice_shower_num_matched_signal = matched_reco_slice_shower_MCP.size();
                m_reco_slice_shower_matched_overlay_frac = matched_reco_slice_shower_overlay_fraction;

                // std::cout<<"checking sim tracks"<<std::endl;
                //then repeat for sim tracks
                for (unsigned int k = 0; k< m_sim_track_parent_pdg.size(); k++){
                    int parent= m_sim_track_parent_pdg[k];
                    int pdg =  m_sim_track_pdg[k];

                    //std::cout<<"for this track at trackID "<<m_sim_track_trackID[k]<< " and pdg "<<pdg <<" the parent is "<<parent<<std::endl;
                    //if this sim track is a photon and it's primary (parent pdg is -1)
                    if((parent == -1 ||parent == 12 || parent ==14 ) && pdg == 2212){
                        //use map from track ID to get MCP
                        int id = m_sim_track_trackID[k];
                        art::Ptr<simb::MCParticle> mcp = MCParticleToTrackIDMap[id];

                        if (m_sim_track_matched[k] > 0){
                            if(std::find(matched_reco_slice_track_MCP.begin(),matched_reco_slice_track_MCP.end(),mcp) ==  matched_reco_slice_track_MCP.end()){

                                matched_reco_slice_track_MCP.push_back(mcp);

                                //  std::cout<<"found a candiate proton track"<<std::endl;
                                //save the overlay fraction and whether it's matched
                                //matched_reco_slice_track_overlay_fraction.push_back(m_reco_slice_track_overlay_fraction[j]);
                                //std::cout<<"found sim proton track with track ID "<<id<<std::endl;
                            }//if not already stored
                            //else{
                            //    std::cout<<"not matched to recob track"<<std::endl;
                            //}

                        }//if matched
                    }//if proton from neutrino interaction
                }//for all sim tracks
                m_reco_slice_track_num_matched_signal = matched_reco_slice_track_MCP.size();


                this->ResizeMatchedSlices(m_reco_slice_shower_num_matched_signal ,m_reco_slice_track_num_matched_signal);  


                std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of sim tracks-MCP matches associated to the true ncdelta is "<<matched_reco_slice_track_MCP.size()<<std::endl;
                //check if either 1g1p or 1g0p topology
                std::cout<<"Status - ";
                //if it's a true 1g1p event
                if (m_mctruth_delta_radiative_1g1p_or_1g1n == 0){
                    //std::cout<<"Status - ";
                    if ( m_reco_slice_track_num_matched_signal == 1 && m_reco_slice_shower_num_matched_signal ==1){
                        std::cout<<"there is one shower one track for 1g1p"<<std::endl;
                    }else if (m_reco_slice_track_num_matched_signal == 1 && m_reco_slice_shower_num_matched_signal == 0){
                        std::cout<<"error, true 1g1p but no shower, 1 track"<<std::endl;
                    } else if (m_reco_slice_track_num_matched_signal == 0 && m_reco_slice_shower_num_matched_signal == 1){
                        std::cout<<"true 1g1p with 1 shower, no track (looks like 1g0p)"<<std::endl;
                    } else {
                        std::cout<<"true 1g1p but neither track nor shower were reconstructed"<<std::endl;
                    }            
                } else{
                    if ( m_reco_slice_track_num_matched_signal == 1 && m_reco_slice_shower_num_matched_signal ==1){
                        std::cout<<"there is one shower one track for 1g0p"<<std::endl;
                    }else if (m_reco_slice_track_num_matched_signal == 1 && m_reco_slice_shower_num_matched_signal == 0){
                        std::cout<<"error, true 1g0p but no shower, 1 track"<<std::endl;
                    } else if (m_reco_slice_track_num_matched_signal == 0 && m_reco_slice_shower_num_matched_signal == 1){
                        std::cout<<"correct, true 1g0p with 1 shower, no track"<<std::endl;
                    } else {
                        std::cout<<"true 1g0p but neither track nor shower were reconstructed"<<std::endl;
                    } 
                }

            }//for nc delta signal
        }//for signal def

        //std::cout<< showerToMCParticleMap.size()<<std::endl;


        //for each MCP, associated to a reco track/shower, get the corresponding PFP and slice
        int i_shr =0;
        for(art::Ptr<simb::MCParticle> mcp: matched_reco_slice_shower_MCP){
            //get the recob shower
            art::Ptr<recob::Shower> this_shr;
            for(auto pair: showerToMCParticleMap){
                if (pair.second == mcp){
                    this_shr = pair.first;    
                }
            }
            art::Ptr<recob::PFParticle> this_pfp;
            if(!this_shr.isNull()){
                this_pfp = showerToPFParticleMap[this_shr];
            }

            //get the slice
            if(!this_pfp.isNull()){
                for(auto pair :allPFPSliceIdVec){
                    art::Ptr<recob::PFParticle> pfp = pair.first;
                    if (this_pfp == pfp){
                if(m_is_verbose)        std::cout<<"found recob shower - MCP at track id "<<mcp->TrackId()<<" in slice "<<pair.second <<std::endl;
                        m_reco_slice_shower_matched_sliceId[i_shr] = pair.second; 
                        m_reco_slice_shower_matched_energy[i_shr] = mcp->E();
                    }
                }
            } else{
              if(m_is_verbose)  std::cout<<"no corresponding slice found for recob shower - MCP at track id "<<mcp->TrackId()<<std::endl;
            }
            i_shr++;
        } 

        int i_trk= 0;
        for(art::Ptr<simb::MCParticle> mcp: matched_reco_slice_track_MCP){
            //get the recob track
            art::Ptr<recob::Track> this_trk;
            for(auto pair: trackToMCParticleMap){
                if (pair.second == mcp){
                    this_trk = pair.first;    
                }
            }
            art::Ptr<recob::PFParticle> this_pfp;
            if(!this_trk.isNull()){
                this_pfp = trackToNuPFParticleMap[this_trk];
            }

            //get the slice
            if(!this_pfp.isNull()){
                for(auto pair :allPFPSliceIdVec){
                    art::Ptr<recob::PFParticle> pfp = pair.first;
                    if (this_pfp == pfp){
                if(m_is_verbose)        std::cout<<"found recob track - MCP at track id "<<mcp->TrackId()<<" in slice "<<pair.second <<std::endl;
                        m_reco_slice_track_matched_sliceId[i_trk] = pair.second;
                        m_reco_slice_track_matched_energy[i_trk]= mcp->E();
                    }
                }
            } else{
                if(m_is_verbose)std::cout<<"no corresponding slice found for recob track - MCP at track id "<<mcp->TrackId()<<std::endl;
            }
            i_trk++;
        }


        //I need to fill these
        /*
           m_reco_slice_shower_matched_conversion 
           */ 

        //if there is a match, store vectors of the showers/tracks and get info like shower energy, conversion distance
        //if there's a partial match then some of the simb things were not reconstructed given the purity/completeness requirements?
        //if there's no match then ??


        //for the given showers/tracks, search for the corresponding slices
        //do the pandora calc for slice correctness
        //good purity:
        //correct = all simb particles in one slice + no additional
        //incorrect = all simb particles in one slice but some additional
        //incorrect = not all simb particles in one slice but good
        //bad purity or not matched at all:
        //incorrect = remaining good simb particles in one slice
        //incorrect = remaining good simb particles in different slices


    }//findslice




    }

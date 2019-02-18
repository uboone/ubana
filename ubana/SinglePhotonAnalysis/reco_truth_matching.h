#include "SinglePhoton_module.h"
#include <climits>
//#include <typeinfo>

namespace single_photon
{

    //same as recoMCmatching but specifically for recob::showers
    std::vector<double> showerRecoMCmatching(std::vector<art::Ptr<recob::Shower>>& objectVector,
            std::map<art::Ptr<recob::Shower>,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
            std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
            std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
            std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap ){


        std::vector<double> vec_fraction_matched;

        //for each recob::track/shower in the event
        for(size_t i=0; i<objectVector.size();++i){
            auto object = objectVector[i];

            //get the associated reco PFP
            const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];

            //check to see where the PFP is in the chain
            // If the particle is primary, it doesn't have a parent
            if (pfp->IsPrimary()){std::cout<<"this is the primary particle"<<std::endl;}
            //else get the parent
            else{
                const auto parentIterator = pfParticleIdMap.find(pfp->Parent());
                if (parentIterator == pfParticleIdMap.end()){
                    std::cout<<"error: no parent but not primary"<<std::endl;
                }

                const int parentPDG = parentIterator->second->PdgCode();
                std::cout<<"the parent pdg code is "<<parentPDG<<std::endl;

                auto daughers_vec = pfp->Daughters();
                std::cout<<"the number of daugter particles is "<<daughers_vec.size() <<std::endl;

                for (auto const daughterpfp: daughers_vec){
                    const auto daughterIterator = pfParticleIdMap.find(daughterpfp);
                    if (daughterIterator == pfParticleIdMap.end()){
                        std::cout<<"error: didn't find that daughter"<<std::endl;
                    } else {
                        art::Ptr<recob::PFParticle> daughters = daughterIterator->second;
                        std::cout<<"the daughter pdg code is "<<daughters->PdgCode()<<std::endl;
                    }
                }               

                //const int parentPDG = parentIterator->second->PdgCode();
                //std::cout<<"the parent pdg code is "<<parentPDG<<std::endl;
            }

            //putting in the PFP pdg code as a check
            int pdg = pfp->PdgCode();

            //and get the hits associated to the reco PFP
            std::vector< art::Ptr<recob::Hit> > obj_hits_ptrs = pfParticleToHitsMap[pfp];

            std::unordered_map<int,double> objide; //map between the MCParticle track ID and the backtracker energy

            //energy for an MCParticle that comprises the most energy when sum over associated hits in PFP
            //total energy of the reco PFP taken from the sum of the hits associated to an MCParticle
            double maxe=-1, tote=0;                

            //simb::MCParticle const * best_matched_mcparticle = NULL; //pointer for the particle match we will calculate
            art::Ptr<simb::MCParticle> best_matched_mcparticle; //pointer for the MCParticle match we will calculate

            //    std::vector<simb::MCParticle const *> particle_vec;
            //    std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

            std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
            std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

            bool found_a_match = false;
            int n_associated_mcparticle_hits = 0;
            std::vector<art::Ptr<simb::MCParticle>> asso_mcparticles_vec;


            std::cout<<"REC: This object has "<<obj_hits_ptrs.size()<<" hits associated with it"<<std::endl;

            //loop only over hits associated to this reco PFP
            for(size_t i_h=0; i_h < obj_hits_ptrs.size(); ++i_h){

                particle_vec.clear(); match_vec.clear(); //only store per hit

                //for the hit, fill the backtracker info 
                mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(), particle_vec, match_vec);

                //mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //the .key() gives us the index in the original collection
                //std::cout<<"REC: hit "<<i_h<<" has "<<particle_vec.size()<<" MCparticles assocaied: "<<std::endl;

                //if there is an MCParticle associated to this hit
                if(particle_vec.size()>0) n_associated_mcparticle_hits++;

                //loop over MCparticles finding which is the MCparticle with most "energy" matched correctly
                //for each MCParticle associated with this hit
                for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                    //add the energy of the back tracked hit for this MCParticle to the track id for the MCParticle in the map
                    objide[ particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy; //store energy per track id

                    //if the id isn't already in the map, store it in the vector of all associated MCParticles
                    if(std::find(asso_mcparticles_vec.begin(), asso_mcparticles_vec.end(),  particle_vec[i_p]) == asso_mcparticles_vec.end()){
                        asso_mcparticles_vec.push_back(particle_vec[i_p]);
                    }

                    //add the energy of the back tracked hit to the total energy for the PFP
                    tote += match_vec[i_p]->energy; //calculate total energy deposited

                    //want the MCParticle with the max total energy summed from the back tracker hit energy from hits in PFP
                    if( objide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                        maxe = objide[ particle_vec[i_p]->TrackId() ];
                        best_matched_mcparticle = particle_vec[i_p];
                        found_a_match = true;
                    }
                }//end loop over particles per hit
            }
            if(n_associated_mcparticle_hits == 0){
                //This will only occur if the whole recob::PFParticle is associated with an overlay object

            }//for each recob::track/shower in the event

            std::cout << "SinglePhoton::recoMC()\t||\t the number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       
            std::cout<<"SinglePhoton::recoMC()\t||\t the stored number of assocaited MCParticles is "<<asso_mcparticles_vec.size()<<std::endl;

            std::map<int, art::Ptr<simb::MCParticle>> mother_MCP_map;
            //art::Ptr<simb::MCParticle> this_mcp;
            int this_mcp_id = -1;
            //bool isValid = true;
            int last_mcp_id = -1;
            for(auto mcp:asso_mcparticles_vec){
                //isValid = true;
                //if (asso_mcparticles_vec.size() < 20){
                std::cout<<"looking at an MCP with pdg code "<<mcp->PdgCode()<<" and status code "<<mcp->StatusCode()<<std::endl;
                //std::cout<<"the mother of this MCP is track id "<<mcp->Mother()<<" and there are "<<mcp->NumberDaughters()<<" daughters"<<std::endl;
                //art::Ptr<simb::MCParticle> this_mother;
                //if this particle has a mother (track id != -1), keep going up
                //if no mother, store the mother
                this_mcp_id = mcp->TrackId();
                last_mcp_id = this_mcp_id;
                //this_mcp_mother_id = mcp->Mother();
                while(this_mcp_id >= 0 ){
                    //while ((MCParticleToTrackIdMap[this_mcp_id]->Mother()) > 0){
                    //move up the chain to the mother
                    art::Ptr<simb::MCParticle> this_mcp = MCParticleToTrackIdMap[this_mcp_id];
                    std::cout<<"going up the tree got mother particle"<<std::endl;

                    if (this_mcp.isNull()){
                        std::cout<<"null pointer at id "<<this_mcp_id<<std::endl;
                        //isValid = false;
                        this_mcp_id = last_mcp_id;
                        break;
                    }

                    std::cout<<"going up the tree at an MCP with track id  "<<this_mcp_id<<", pdg code "<<this_mcp->PdgCode()<<", and status code "<<this_mcp->StatusCode()<<std::endl;

                    last_mcp_id = this_mcp_id;
                    this_mcp_id =  this_mcp->Mother();


                    std::cout<<"-------------------"<<std::endl;

                    //this_mcp_id =  MCParticleToTrackIdMap[this_mcp_id]->Mother();
                    //std::cout<<"going up the tree at an MCP with track id  "<<this_mcp_id<<", pdg code "<<this_mcp->PdgCode()<<", and status code "<<this_mcp->StatusCode()<<std::endl;

                    //if( this_mcp_id == 88){
                    //    std::cout<<"the mother id for particle with id 88 = "<< MCParticleToTrackIdMap[this_mcp_id]->Mother()<<std::endl;
                    //} 
                }

                std::cout<<"reached the end of the loop with track ID "<<this_mcp_id<<std::endl;

                if (this_mcp_id >= 0){
                    std::cout<<"storing the mother mother particle with track id "<<this_mcp_id<<" and pdg code "<<MCParticleToTrackIdMap[this_mcp_id]->PdgCode()<<" and status code "<<MCParticleToTrackIdMap[this_mcp_id]->StatusCode()<<std::endl;
                    mother_MCP_map[this_mcp_id] = MCParticleToTrackIdMap[this_mcp_id];
                } else{
                    std::cout<<"error, the mother mother id was "<<this_mcp_id <<std::endl;
                }
                }//for each MCParticle that's associated to a the recob::Shower

                //}

                std::cout<<"SinglePhoton::recoMC()\t||\t the number of source mother particles is "<<mother_MCP_map.size()<<std::endl;

                
                //for each source mother particle, if it's a photon, follow the chain and sum the hits
                for(auto pair: mother_MCP_map){
                    art::Ptr<simb::MCParticle> particle = pair.second;
                    int numDaughters = -1;
                    if (particle->PdgCode()==22){
                        std::cout<<"starting from photon mother particle at head of chain"<<std::endl;
                        numDaughters = particle->NumberDaughters();
                        std::cout<<"this particle has "<<numDaughters<<" daughters"<<std::endl;
                        //loop over the daughters
                        for (int i = 0; i < numDaughters; i++){
                            int daughterId = particle->Daughter(i);
                            std::cout<<"the daughter track Id from mother mother particle is"<<daughterId<<std::endl;
                        } 

                    }
                }

                //for (auto )   


                if(found_a_match){
                    mcParticleVector.push_back(best_matched_mcparticle);
                    objectToMCParticleMap[object] = mcParticleVector.back();
                }else{
                    // mcParticleVector.push_back(0);
                }
                vec_fraction_matched.push_back(maxe/tote);
                // if(m_is_verbose){
                //     std::cout << "SinglePhoton::recoMC()\t||\t the fracrion matched is "<<maxe/tote<<std::endl;
                // }


                if(!found_a_match){
                    std::cout << "SinglePhoton::recoMC()\t||\t NO MATCH NO MATCH (from my loop) for PFP with pdg  "<<pdg<<std::endl;
                    std::cout<<" count "<<objectToMCParticleMap.count(object)<<std::endl;
                }else{
                    std::cout << "SinglePhoton::recoMC()\t||\t Final Match (from my loop) for PFP with pdg "<<pdg<<" is " << best_matched_mcparticle->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
                        << " pdg=" << best_matched_mcparticle->PdgCode()
                        << " trkid=" << best_matched_mcparticle->TrackId()
                        << " ke=" << best_matched_mcparticle->E()-best_matched_mcparticle->Mass()<< "\n";
                }

            }//end vector loop.
            return vec_fraction_matched;
        }



        //Typenamed for recob::Track and recob::Shower
        template<typename T>
            std::vector<double> recoMCmatching(std::vector<T>& objectVector,
                    std::map<T,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
                    std::map<T,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
                    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
                    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector){

                std::vector<double> vec_fraction_matched;

                //for each recob::track/shower in the event
                for(size_t i=0; i<objectVector.size();++i){
                    auto object = objectVector[i];

                    //get the associated reco PFP
                    const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];

                    //putting in the PFP pdg code as a check
                    int pdg = pfp->PdgCode();

                    //and get the hits associated to the reco PFP
                    std::vector< art::Ptr<recob::Hit> > obj_hits_ptrs = pfParticleToHitsMap[pfp];

                    std::unordered_map<int,double> objide; //map between the MCParticle track ID and the backtracker energy

                    //energy for an MCParticle that comprises the most energy when sum over associated hits in PFP
                    //total energy of the reco PFP taken from the sum of the hits associated to an MCParticle
                    double maxe=-1, tote=0;                

                    //simb::MCParticle const * best_matched_mcparticle = NULL; //pointer for the particle match we will calculate
                    art::Ptr<simb::MCParticle> best_matched_mcparticle; //pointer for the MCParticle match we will calculate

                    //    std::vector<simb::MCParticle const *> particle_vec;
                    //    std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

                    std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
                    std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

                    bool found_a_match = false;
                    int n_associated_mcparticle_hits = 0;

                    std::cout<<"REC: This object has "<<obj_hits_ptrs.size()<<" hits associated with it"<<std::endl;

                    //loop only over hits associated to this reco PFP
                    for(size_t i_h=0; i_h < obj_hits_ptrs.size(); ++i_h){

                        particle_vec.clear(); match_vec.clear(); //only store per hit

                        //for the hit, fill the backtracker info 
                        mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(), particle_vec, match_vec);

                        //mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
                        //the .key() gives us the index in the original collection
                        //std::cout<<"REC: hit "<<i_h<<" has "<<particle_vec.size()<<" MCparticles assocaied: "<<std::endl;

                        //if there is an MCParticle associated to this hit
                        if(particle_vec.size()>0) n_associated_mcparticle_hits++;

                        //loop over MCparticles finding which is the MCparticle with most "energy" matched correctly
                        //for each MCParticle associated with this hit
                        for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                            //add the energy of the back tracked hit for this MCParticle to the track id for the MCParticle in the map
                            objide[ particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy; //store energy per track id

                            //add the energy of the back tracked hit to the total energy for the PFP
                            tote += match_vec[i_p]->energy; //calculate total energy deposited

                            //want the MCParticle with the max total energy summed from the back tracker hit energy from hits in PFP
                            if( objide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                                maxe = objide[ particle_vec[i_p]->TrackId() ];
                                best_matched_mcparticle = particle_vec[i_p];
                                found_a_match = true;
                            }
                        }//end loop over particles per hit
                    }
                    if(n_associated_mcparticle_hits == 0){
                        //This will only occur if the whole recob::PFParticle is associated with an overlay object

                    }//for each recob::track/shower in the event

                    std::cout << "SinglePhoton::recoMC()\t||\t the number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       

                    if(found_a_match){
                        mcParticleVector.push_back(best_matched_mcparticle);
                        objectToMCParticleMap[object] = mcParticleVector.back();
                    }else{
                        // mcParticleVector.push_back(0);
                    }
                    vec_fraction_matched.push_back(maxe/tote);
                    // if(m_is_verbose){
                    //     std::cout << "SinglePhoton::recoMC()\t||\t the fracrion matched is "<<maxe/tote<<std::endl;
                    // }


                    if(!found_a_match){
                        std::cout << "SinglePhoton::recoMC()\t||\t NO MATCH NO MATCH (from my loop) for PFP with pdg  "<<pdg<<std::endl;
                        std::cout<<" count "<<objectToMCParticleMap.count(object)<<std::endl;
                    }else{
                        std::cout << "SinglePhoton::recoMC()\t||\t Final Match (from my loop) for PFP with pdg "<<pdg<<" is " << best_matched_mcparticle->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
                            << " pdg=" << best_matched_mcparticle->PdgCode()
                            << " trkid=" << best_matched_mcparticle->TrackId()
                            << " ke=" << best_matched_mcparticle->E()-best_matched_mcparticle->Mass()<< "\n";
                    }

                }//end vector loop.
                return vec_fraction_matched;
            }


        //Typenamed for simb::MCTrack and sim::MCShower
        template<typename T>
            void perfectRecoMatching(
                    std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
                    std::vector<T>& mcObjectVector,
                    std::map<art::Ptr<simb::MCParticle>,T>& mcParticleToMCObjectMap
                    ){

                for(size_t ip=0; ip<mcParticleVector.size(); ++ip){
                    const art::Ptr<simb::MCParticle> particle = mcParticleVector[ip];
                    int particle_trackID = particle->TrackId();

                    std::vector<int> id_matches;
                    std::vector<int> mother_id_matches;
                    std::vector<int> ancestor_id_matches;

                    for(size_t io=0; io<mcObjectVector.size(); ++io){

                        const T object = mcObjectVector[io];
                        int object_trackID =object->TrackID(); 
                        int object_mother_trackID =object->MotherTrackID(); 
                        int object_ancestor_trackID =object->AncestorTrackID();

                        if(object_trackID == particle_trackID ) id_matches.push_back(io);
                        if(object_mother_trackID == particle_trackID ) mother_id_matches.push_back(io);
                        if(object_ancestor_trackID == particle_trackID ) ancestor_id_matches.push_back(io);
                    }

                    int num_id_matches=id_matches.size();
                    int num_mother_id_matches=mother_id_matches.size();
                    int num_ancestor_id_matches=ancestor_id_matches.size();

                    //So im not sure how this works but something like this
                    if(num_id_matches > 1){
                        std::cout<<"Well hot saussage.. more than 1 id match "<<num_id_matches<<std::endl;
                    }else if(num_id_matches == 1){
                        //We have a direct match?
                        mcParticleToMCObjectMap[particle] = mcObjectVector[id_matches.front()];
                    }else if(num_mother_id_matches == 1){
                        //We have a mother match? I guess this is like "Photon?->e+e-"
                        //mcParticleToMCObjectMap[particle] = mcObjectVector[mother_id_matches.front()];
                    }else if(num_ancestor_id_matches == 1){
                        //We have a mother match? I guess this is like Neutron->photon->e+e-"
                        //mcParticleToMCObjectMap[particle] = mcObjectVector[ancestor_id_matches.front()];
                    }else if(num_ancestor_id_matches >1 || num_mother_id_matches >1){
                        std::cout<<"Well hot saussage.. more than 1 mother or ancestor. Actually thats very reasonable hehe."<<num_mother_id_matches<<" "<<num_ancestor_id_matches<<std::endl;
                    }else if(num_id_matches == 0 && num_ancestor_id_matches == 0 && num_mother_id_matches ==0){
                        std::cout<<"NO matches for trackid, mother trackid or ancestor trackid. Hmm"<<num_mother_id_matches<<" "<<num_ancestor_id_matches<<std::endl;
                    }

                    //What if multiple mothers matches?! no idea.


                }//MCParticleLoop

                return;
            }


        void testbed( std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector, const art::Event &evt){

            std::map<int,art::Ptr<simb::MCParticle> > crap_map;
            for(size_t j=0;j< mcParticleVector.size();j++){
                const art::Ptr<simb::MCParticle> mcp = mcParticleVector[j];
                //std::cout<<"PARG: "<<j<<" trackid: "<<mcp->TrackId()<<" key: "<<mcp.key()<<std::endl;
                crap_map[mcp->TrackId()] = mcParticleVector[mcp.key()];
            }
            art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= evt.getValidHandle<std::vector<simb::MCParticle>>("largeant");
            art::FindManyP<simb::MCTruth,sim::GeneratedParticleInfo> genieMCTruthHandle(mcParticleHandle, evt, "largeant");

            std::vector<art::Ptr<simb::MCTruth>> GenieMCTruth;
            std::vector<sim::GeneratedParticleInfo const *> geninfo;



            for(size_t i=0; i< mcParticleVector.size();i++){


                art::Ptr<simb::MCParticle> nth_mother = mcParticleVector[i];

                //if(nth_mother->PdgCode() != 22 && nth_mother->PdgCode() != 11) continue;

                std::cout<<"----------------------------------------------------------"<<std::endl;
                std::cout<<"SinglePhoton::testbed()\t||\t On Particle (trackid: "<<nth_mother->TrackId()<<") pdg: "<<nth_mother->PdgCode()<<", status_code "<<nth_mother->StatusCode()<<"  MotherTrackID: "<<nth_mother->Mother()<<std::endl;

                int n_generation = 1;

                while(nth_mother->Mother() != 0){

                    nth_mother = crap_map[nth_mother->Mother()]; 
                    std::cout<<"SinglePhoton::testbed()\t||\t -- and "<<n_generation<<"-mother trackid "<<nth_mother->TrackId()<<" is a pdg: "<<nth_mother->PdgCode()<<" and status_code "<<nth_mother->StatusCode()<<std::endl;
                    n_generation++;
                }



                GenieMCTruth.clear(); geninfo.clear(); //tidy up this loop
                genieMCTruthHandle.get(mcParticleVector[i].key(),GenieMCTruth,geninfo);

                std::cout<<"SinglePhoton::testbed()\t||\t "<<" GenieMCTruth.size() "<<GenieMCTruth.size()<<" geninfo.size() "<<geninfo.size()<<std::endl;
                for(size_t k=0; k< GenieMCTruth.size(); k++){
                    std::cout<<"SinglePhoton::testbed()\t||\t -- "<<k<<": has "<<GenieMCTruth[k]->NParticles()<<" particlesand geninfo_index: "<<geninfo[k]->generatedParticleIndex()<<std::endl;
                    if((int)geninfo[k]->generatedParticleIndex() > GenieMCTruth[k]->NParticles() || geninfo[k]->generatedParticleIndex()==ULONG_MAX){
                        std::cout<<"SinglePhoton::testbed()\t||\t -- Thats way wrong.."<<std::endl;
                    }else{
                        const simb::MCParticle mp = GenieMCTruth[k]->GetParticle(geninfo[k]->generatedParticleIndex());
                        std::cout<<"SinglePhoton::testbed()\t||\t -- is a pdg: "<<mp.PdgCode()<<"  with statuscode:"<<mp.StatusCode()<<std::endl;
                    }
                    //std::cout<<"SinglePhoton::testbed()\t||\t "<<" "<<GenieMCTruth[0]->NParticles()<<" "<<geninfo[0]->generatedParticleIndex()<<std::endl;
                }


            }//particleloop
        }
    }//namespace end

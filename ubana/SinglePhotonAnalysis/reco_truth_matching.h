#include "SinglePhoton_module.h"


namespace single_photon
{

    //Typenamed for recob::Track and recob::Shower
    template<typename T>
        std::vector<double> recoMCmatching(std::vector<T>& objectVector,
                std::map<T,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
                std::map<T,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
                std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
                art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
                std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector){

            std::vector<double> vec_fraction_matched;

            for(size_t i=0; i<objectVector.size();++i){
                auto object = objectVector[i];
                const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];

                std::vector< art::Ptr<recob::Hit> > obj_hits_ptrs = pfParticleToHitsMap[pfp];

                std::unordered_map<int,double> objide;
                double maxe=-1, tote=0;
                //simb::MCParticle const * best_matched_mcparticle = NULL; //pointer for the particle match we will calculate
                art::Ptr<simb::MCParticle> best_matched_mcparticle; //pointer for the particle match we will calculate

                //    std::vector<simb::MCParticle const *> particle_vec;
                //    std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

                std::vector<art::Ptr<simb::MCParticle>> particle_vec;
                std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

                bool found_a_match = false;
                int n_associated_mcparticle_hits = 0;

                std::cout<<"REC: This object has "<<obj_hits_ptrs.size()<<" hits associated with it"<<std::endl;

                //loop only over our hits
                for(size_t i_h=0; i_h < obj_hits_ptrs.size(); ++i_h){

                    particle_vec.clear(); match_vec.clear(); //tidy up this loop

                    mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(), particle_vec, match_vec);
                    //mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
                    //the .key() gives us the index in the original collection
                    //std::cout<<"REC: hit "<<i_h<<" has "<<particle_vec.size()<<" MCparticles assocaied: "<<std::endl;
                    if(particle_vec.size()>0) n_associated_mcparticle_hits++;


                    //loop over MCparticles finding which is the MCparticle with most "energy" matched correctly
                    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                        objide[ particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy; //store energy per track id
                        tote += match_vec[i_p]->energy; //calculate total energy deposited
                        if( objide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                            maxe = objide[ particle_vec[i_p]->TrackId() ];
                            best_matched_mcparticle = particle_vec[i_p];
                            found_a_match = true;
                        }
                    }//end loop over particles per hit
                }
                if(n_associated_mcparticle_hits == 0){
                    //This will only occur if the whole recob::PFParticle is associated with an overlay object

                }


                if(found_a_match){
                    mcParticleVector.push_back(best_matched_mcparticle);
                    objectToMCParticleMap[object] = mcParticleVector.back();
                }else{
                   // mcParticleVector.push_back(0);
                }
                vec_fraction_matched.push_back(maxe/tote);



                if(!found_a_match){
                    std::cout << "SinglePhoton::recoMC()\t||\t NO MATCH NO MATCH (from my loop)  "<<std::endl;
                    std::cout<<" count "<<objectToMCParticleMap.count(object)<<std::endl;
                }else{
                    std::cout << "SinglePhoton::recoMC()\t||\t Final Match (from my loop) is " << best_matched_mcparticle->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
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
}//namespace end

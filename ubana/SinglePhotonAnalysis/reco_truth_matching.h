#include "SinglePhoton_module.h"


namespace single_photon
{


    template<typename T>
    void recoMCmatching(std::vector<T>& objectVector,
                        std::map<T,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
                        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
                        art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit ){

        for(size_t i=0; i<objectVector.size();++i){
            const T object = objectVector[i];
            const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];

            std::vector< art::Ptr<recob::Hit> > obj_hits_ptrs = pfParticleToHitsMap[pfp];

            std::unordered_map<int,double> objide;
            double maxe=-1, tote=0;
            simb::MCParticle const* maxp_me = NULL; //pointer for the particle match we will calculate

            //std::vector<art::Ptr<simb::MCParticle>> particle_vec; //will hold all MCparticles that "Share" a recob::Hit
            //std::vector<art::Ptr<anab::BackTrackerHitMatchingData>> match_vec;//hold matching data, that I dont think I use at the moment
            std::vector<simb::MCParticle const *> particle_vec;
            std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

            //loop only over our hits
            for(size_t i_h=0; i_h<obj_hits_ptrs.size(); ++i_h){

                particle_vec.clear(); match_vec.clear(); //tidy up this loop
                mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //the .key() gives us the index in the original collection

     //           std::cout << "PLARGGGG \t\tThere are " << particle_vec.size() << " particles matched to hit " << i_h << "\n";

                //loop over MCparticles finding which is the MCparticle with most "energy" matched correctly
                for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                    objide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
                    tote += match_vec[i_p]->energy; //calculate total energy deposited
                    if( objide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                        maxe = objide[ particle_vec[i_p]->TrackId() ];
                        maxp_me = particle_vec[i_p];
                    }
                }//end loop over particles per hit

            }


            std::cout << "SinglePhoton::recoMC()\t||\t Final Match (from my loop) is " << maxp_me->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
                << " pdg=" << maxp_me->PdgCode()
                << " trkid=" << maxp_me->TrackId()
                << " ke=" << maxp_me->E()-maxp_me->Mass()
                << "\nSinglePhoton::recoMC()\t||\t start (x,y,z)=(" << maxp_me->Vx()
                << "," << maxp_me->Vy()
                << "," << maxp_me->Vz()
                << ")\tend (x,y,z)=(" << maxp_me->EndX()
                << "," << maxp_me->EndY()
                << "," << maxp_me->EndZ() << ")" << "\n";



        }
    }

}//namespace end

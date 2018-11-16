#include "SinglePhoton_module.h"


namespace single_photon
{

    //Typenamed for recob::Track and recob::Shower
    template<typename T>
        int badChannelMatching(std::vector<int>& badchannels, 
                std::vector<T>& objects, 
                std::map< T, art::Ptr<recob::PFParticle> > & objectToPFParticleMap,
                std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,
                const geo::GeometryCore * geom,
                std::vector<std::pair<int,int>> & bad_channel_list_fixed_mcc9){




            for(size_t i_trk=0; i_trk<objects.size(); i_trk++){


                const T object = objects[i_trk];
                const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];
                const std::vector<art::Ptr<recob::Hit>> hits = pfParticleToHitsMap[pfp];

                int min_dist_from_bad_channel = 99999;

                for(size_t h=0; h<hits.size(); h++){
                const art::Ptr<recob::Hit> hit = hits[h];
                    
            
                    //int nch = (int)badchannels.size()/3;
                    //for(int i=0; i<nch; i++){
                    for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                      //  const int offset = 3*i;
                      //  int bc = badchannels[offset+0];                                       
                        int bc = bad_channel_list_fixed_mcc9[i].first;                                       
                        int ok = bad_channel_list_fixed_mcc9[i].second;       
                        if(ok>1)continue;
                        int dist =hit->Channel()-bc;
                        auto hs = geom->ChannelToWire(bc);
                        //std::cout<<"AG: "<<hs.size()<<"  BC("<<bc<<"): "<<hs[0]<<" ours: ("<<hit->Channel()<<"): "<<hit->WireID()<<std::endl;
                        //this is the right format for my plotting routine
                        //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<< badchannels[offset+1]<<" "<<badchannels[offset+2]<<std::endl;
                        std::vector<double> start(3);
                        std::vector<double> end(3);
                        auto result = geom->WireEndPoints(hs[0]);
                        
                //        std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 

                        

                        if(fabs(dist) < min_dist_from_bad_channel) min_dist_from_bad_channel = fabs(dist);
                   
            //            std::cout<<"BD: "<<badchannels[offset+0]<<" "<<badchannels[offset+1]<<" "<<badchannels[offset+2]<<" "<< fabs(badchannels[offset+1]-badchannels[offset+2])<<" "<<dist<<std::endl;
                    }
                }//hitloop
            }

            return 0;
        }

}//namespace end

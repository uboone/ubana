#include "SinglePhoton_module.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLine.h"

namespace single_photon
{
    void SinglePhoton::ClearSecondShowers(){
        m_sss_num_unassociated_hits=0;
        m_sss_num_associated_hits=0;

    }

    void SinglePhoton::ResizeSecondShowers(size_t size){

    }


    void SinglePhoton::CreateSecondShowerBranches(){
        vertex_tree->Branch("sss_num_unassociated_hits",&m_sss_num_unassociated_hits,"sss_num_unassociated_hits/I");
        vertex_tree->Branch("sss_num_associated_hits",&m_sss_num_associated_hits,"sss_num_associated_hits/I");

    }

    void SinglePhoton::SecondShowerSearch(
            const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
            const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
            const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
            const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap){


        int total_track_hits =0;
        int total_shower_hits =0;
        int nu_slice_id = -999;

        std::vector<art::Ptr<recob::Hit>> associated_hits;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits_plane0;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits_plane1;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits_plane2;


        for(size_t t =0; t< tracks.size(); t++){
            art::Ptr<recob::Track> track = tracks[t];
            art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];
            int sliceid = pfParticleToSliceIDMap.at(pfp);
            auto slicehits = sliceIDToHitsMap.at(sliceid);
            auto trackhits = pfParticleToHitsMap.at(pfp);

            std::cout<<"SSS: track "<<t<<" is in slice "<<sliceid<<" which has "<<slicehits.size()<<" hits. This track has  "<<trackhits.size()<<" of them. "<<std::endl;
            total_track_hits+=trackhits.size();
            if(nu_slice_id !=  sliceid && nu_slice_id != -999){
                std::cout<<"ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
                exit(EXIT_FAILURE);
            }   
            nu_slice_id = sliceid;


            for(auto &h: trackhits){
                associated_hits.push_back(h);
            }

        }

        for(size_t s =0; s< showers.size(); s++){
            art::Ptr<recob::Shower> shower = showers[s];
            art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(shower);
            int sliceid = pfParticleToSliceIDMap.at(pfp);
            auto slicehits = sliceIDToHitsMap.at(sliceid);
            auto showerhits = pfParticleToHitsMap.at(pfp);

            std::cout<<"SSS: shower "<<s<<" is in slice "<<sliceid<<" which has "<<slicehits.size()<<" hits. This shower has  "<<showerhits.size()<<" of them. "<<std::endl;
            total_shower_hits+=showerhits.size();
            if(nu_slice_id !=  sliceid && nu_slice_id!=-999){
                std::cout<<"ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
                exit(EXIT_FAILURE);
            }
            nu_slice_id = sliceid;

            for(auto &h: showerhits){
                associated_hits.push_back(h);
            }


        }

        std::cout<<"SSS: So in total we have "<<total_shower_hits<<" shower hits and "<<total_track_hits<<" track hits"<<" assocaedHits total: "<<associated_hits.size()<<std::endl;
        m_sss_num_associated_hits = total_shower_hits+total_track_hits;

        if(nu_slice_id>=0){
            std::cout<<"SSS: So that leaves "<<sliceIDToHitsMap.at(nu_slice_id).size()-total_shower_hits-total_track_hits<<" hits not included in tracks and showers"<<std::endl;
            m_sss_num_unassociated_hits = sliceIDToHitsMap.at(nu_slice_id).size()-total_shower_hits-total_track_hits;

            if(m_sss_num_unassociated_hits<0){
                std::cout<<"ERROR!! Number of unassociated hits is negative, i.e: num_associated: "<<m_sss_num_associated_hits<<" and total slice hits: "<<sliceIDToHitsMap.at(nu_slice_id).size()<<std::endl;
                exit(EXIT_FAILURE);
            }




            std::vector<art::Ptr<recob::Hit>> slicehits = sliceIDToHitsMap.at(nu_slice_id);
            for(auto &h: slicehits){

                bool is_associated = false;
                for(auto &a: associated_hits){
                    if(h==a){
                        is_associated = true;
                        break;
                    }
                }

                if(!is_associated){
                    unassociated_hits.push_back(h);
                    auto plane_view = h->View();
                    switch((int)plane_view){
                        case (0) :
                            unassociated_hits_plane0.push_back(h);
                            break;
                        case (1) :
                            unassociated_hits_plane1.push_back(h);
                            break;
                        case (2) :
                            unassociated_hits_plane2.push_back(h);
                            break;
                    }

                }

            }

            std::vector<std::vector<art::Ptr<recob::Hit>>> unassociated_hits_all = {unassociated_hits_plane0,unassociated_hits_plane1,unassociated_hits_plane2};

            std::cout<<" associated_hits.size() "<<associated_hits.size()<<" unassociated_hits.size() "<<unassociated_hits.size()<<" p0: "<<unassociated_hits_plane0.size()<<" p1:  "<<unassociated_hits_plane1.size()<<" p2: "<<unassociated_hits_plane2.size()<<std::endl;


            if(bool_make_sss_plots && showers.size()>0){

                //TFile *f = new TFile("t.root","recreate");
                //f->cd();

                std::string print_name = "sss_"+std::to_string(m_run_number)+"_"+std::to_string(m_subrun_number)+"_"+std::to_string(m_event_number);
                TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),2400,1600);
                can->Divide(3,2,0,0.1);

                double tick_max = 0;
                double tick_min = 1e10;
                std::vector<double> chan_max(3,0);
                std::vector<double> chan_min(3,1e10);


                //First grab all shower clusters
                std::vector<std::vector<TGraph *>> pts_shr( showers.size(), std::vector<TGraph *>(3)  );

                for(size_t s =0; s< showers.size(); s++){
                    art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(showers[s]);
                    auto showerhits = pfParticleToHitsMap.at(pfp);

                    std::vector<TGraph*> t_pts(3);
                    std::vector<std::vector<double>>  vec_t(3);
                    std::vector<std::vector<double>>  vec_c(3);

                    for(auto &h: showerhits){
                        double wire = (double)h->WireID().Wire;
                        vec_c[(int)h->View()].push_back(wire);
                        //vec_c[(int)h->View()].push_back((double)h->Channel());
                        vec_t[(int)h->View()].push_back((double)h->PeakTime());
                        tick_max = std::max(tick_max, (double)h->PeakTime());
                        tick_min = std::min(tick_min, (double)h->PeakTime());
                        chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                        chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);
                    }
                    t_pts[0] = new TGraph(vec_c[0].size(), &(vec_c[0])[0], &(vec_t[0])[0]);
                    t_pts[1] = new TGraph(vec_c[1].size(), &(vec_c[1])[0], &(vec_t[1])[0]);
                    t_pts[2] = new TGraph(vec_c[2].size(), &(vec_c[2])[0], &(vec_t[2])[0]);

                    pts_shr[s] = t_pts;
                }


                std::vector<std::vector<TGraph *>> pts_trk( tracks.size(), std::vector<TGraph *>(3)  );
                for(size_t t =0; t< tracks.size(); t++){
                    art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap.at(tracks[t]);
                    auto trackhits = pfParticleToHitsMap.at(pfp);

                    std::vector<TGraph*> t_pts(3);
                    std::vector<std::vector<double>>  vec_t(3);
                    std::vector<std::vector<double>>  vec_c(3);

                    for(auto &h: trackhits){
                        double wire = (double)h->WireID().Wire;
                        vec_c[(int)h->View()].push_back(wire);
                        vec_t[(int)h->View()].push_back((double)h->PeakTime());
                        tick_max = std::max(tick_max, (double)h->PeakTime());
                        tick_min = std::min(tick_min, (double)h->PeakTime());
                        chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                        chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

                    }
                    t_pts[0] = new TGraph(vec_c[0].size(), &(vec_c[0])[0], &(vec_t[0])[0]);
                    t_pts[1] = new TGraph(vec_c[1].size(), &(vec_c[1])[0], &(vec_t[1])[0]);
                    t_pts[2] = new TGraph(vec_c[2].size(), &(vec_c[2])[0], &(vec_t[2])[0]);

                    pts_trk[t] = t_pts;
                }
                //Now the "Unassociated Hits"

                std::vector<TGraph *> g_unass(3);
                std::vector<std::vector<std::vector<double>>> pts_to_recluster(3);
                for(int i=0; i<3; i++){

                    std::vector<double> vec_t;
                    std::vector<double> vec_c;

                    for(auto &h: unassociated_hits_all[i]){
                        double wire = (double)h->WireID().Wire;
                        vec_c.push_back(wire);
                        vec_t.push_back((double)h->PeakTime());
                        tick_max = std::max(tick_max, (double)h->PeakTime());
                        tick_min = std::min(tick_min, (double)h->PeakTime());
                        chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                        chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

                        //for reclustering
                        std::vector<double> pt = {wire,vec_t.back()};
                        pts_to_recluster[(int)h->View()].push_back(pt);

                    }

                    g_unass[i] = new TGraph(vec_c.size(), &vec_c[0], &vec_t[0]);


                }
                //Plotting now


                std::vector<int> tcols = {kRed-7, kBlue-7, kGreen-3, kOrange-3, kCyan-3, kMagenta-3, kGreen+1 , kRed+1};
                int used_col=0;
                if(showers.size()+tracks.size() > tcols.size()){
                    for(int i =0; i< (int)(showers.size()+tracks.size() - tcols.size() +2); i++){
                        tcols.push_back((int)rangen->Uniform(400,900));
                    }
                }

                std::cout<<"Tick Min: "<<tick_min<<" Max: "<<tick_max<<std::endl;
                auto const TPC = (*geom).begin_TPC();
                auto ID = TPC.ID();
                int fCryostat = ID.Cryostat;
                int fTPC = ID.TPC;
                std::cout<<TPC.ID()<<"= the beginning TPC ID" <<std::endl;
                std::cout<<"the cryostat id = "<<fCryostat<<std::endl;  
                std::cout<<"the tpc id = "<<fTPC<<std::endl;  


                //Plotting the vertex position on the plot.
                std::vector<TGraph*> g_vertex(3);
                for(int i=0; i<3; i++){
                    can->cd(i+1);

                    std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
                    std::vector<double>  time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, *theDetector)};

                    chan_max[i] = std::max( chan_max[i],wire[0]);
                    chan_min[i] = std::min( chan_min[i],wire[0]);

                    g_vertex[i] = new TGraph(1,&wire[0],&time[0]);
                    g_vertex[i]->SetMarkerStyle(29);
                    g_vertex[i]->SetMarkerSize(4);
                    g_vertex[i]->SetMarkerColor(kMagenta-3);
                    g_vertex[i]->GetYaxis()->SetRangeUser(tick_min*0.9,tick_max*1.1);
                    g_vertex[i]->GetXaxis()->SetLimits(chan_min[i]*0.9,chan_max[i]*1.1);
                    g_vertex[i]->SetTitle(("Plane " +std::to_string(i)).c_str());
                    g_vertex[i]->Draw("ap");

                    if(i>0){
                        g_vertex[i]->GetYaxis()->SetLabelOffset(999);
                        g_vertex[i]->GetYaxis()->SetLabelSize(0);
                    }


                    can->cd(i+4);
                    g_vertex[i]->Draw("ap");

                }


                for(size_t t=0; t< pts_trk.size(); t++){
                    int tcol = tcols[used_col];
                    used_col++;

                    for(int i=0; i<3; i++){
                        can->cd(i+1);
                        if(pts_trk[t][i]->GetN()>0){//need a check in case this track has no hits on this plane.
                            pts_trk[t][i]->Draw("p same"); 
                            pts_trk[t][i]->SetMarkerColor(tcol);
                            pts_trk[t][i]->SetMarkerStyle(20);
                            pts_trk[t][i]->SetMarkerSize(0.75);
                        }
                    }
                }

                for(size_t t=0; t< pts_shr.size(); t++){
                    int tcol = tcols[used_col];
                    used_col++;

                    for(int i=0; i<3; i++){
                        can->cd(i+1);
                        if(pts_shr[t][i]->GetN()>0){
                            pts_shr[t][i]->Draw("p same"); //used in the vertex
                            pts_shr[t][i]->SetMarkerColor(tcol);
                            pts_shr[t][i]->SetMarkerStyle(20);
                            pts_shr[t][i]->SetMarkerSize(0.75);
                        }
                    }
                }


                for(int i=0; i<3; i++){
                    can->cd(i+1);

                    if(g_unass[i]->GetN()>0){
                        g_unass[i]->Draw("p same");
                        g_unass[i]->SetMarkerColor(kBlack);
                        g_unass[i]->SetMarkerStyle(20);
                        g_unass[i]->SetMarkerSize(0.75);
                    }
                    g_vertex[i]->Draw("p same");
                }



                //******************************** DeadWireRegions********************************************
                for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                    int badchan = bad_channel_list_fixed_mcc9[i].first;                                       
                    int ok = bad_channel_list_fixed_mcc9[i].second;       

                    if(ok>1){
                        auto hs = geom->ChannelToWire(badchan);

                        //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 
                        //                    std::cout<<wireids[0].Plane<<" "<<result.start().X()<<std::endl;
                        int thisp = (int)hs[0].Plane;
                        double bc = hs[0].Wire;


                        if(chan_min[thisp]*0.9> bc && bc < chan_max[thisp]*1.1 ){
                            can->cd(thisp+1);
                            //TLine *l = new TLine(bc,tick_min*0.9,bc,tick_max*1.1);
                            //l->Draw("same");
                        }
                    }
                }



                //*****************************DBSCAN***********************************
                int min_pts = 5;
                double eps = 100.0;
                std::vector<int> num_clusters(3,0);

                std::vector<std::vector<TGraph*>> g_clusters(3);
                std::vector<std::vector<int>> cluster_labels(3);
                for(int i=0; i<3; i++){

                    std::cout<<"Starting to run DBSCAN for plane: "<<i<<" has "<<pts_to_recluster[i].size()<<" pts to do using eps: "<<eps<<" and min_pts: "<<min_pts<<std::endl; 
                    DBSCAN ReCluster(eps,min_pts);
                    cluster_labels[i] =  ReCluster.Scan2D(pts_to_recluster[i]);

                    for(auto &c: cluster_labels[i]){
                        num_clusters[i] = std::max(c,num_clusters[i]);
                    }
                    std::cout<<"On this plane, DBSCAN found: "<<num_clusters[i]<<" clusters"<<std::endl;
                }

                //Step next, loop over and make plots again
                for(int i=0; i<3; i++){
                    std::vector<std::vector<double>> vec_time(num_clusters[i]+1);
                    std::vector<std::vector<double>> vec_wire(num_clusters[i]+1);
                    std::vector<TGraph*> tmp_g_clusters(num_clusters[i]+1);

                    if(cluster_labels[i].size() != pts_to_recluster[i].size()){
                        std::cout<<"ERROR!! someting amiss cluster labels of size "<<cluster_labels[i].size()<<" and pts  in this plane "<<pts_to_recluster[i].size()<<std::endl;
                    }  

                    for(size_t k=0; k< pts_to_recluster[i].size(); k++){
                        //std::cout<<vec_wire.size()<<" "<<cluster_labels[i][k]<<std::endl;
                        vec_wire[cluster_labels[i][k]].push_back(pts_to_recluster[i][k][0]); 
                        vec_time[cluster_labels[i][k]].push_back(pts_to_recluster[i][k][1]); 

                    }

                    for(int c=0; c< num_clusters[i]+1; c++){
                        int tcol = kBlack;
                        if(c>0) tcol = rangen->Uniform(400,900);//tcols[tcols.size()-c];
                        tmp_g_clusters[c] = new TGraph(vec_wire[c].size(),&(vec_wire[c])[0],&(vec_time[c])[0] );
                        can->cd(i+4);
                        if(tmp_g_clusters[c]->GetN()>0){
                            tmp_g_clusters[c]->Draw("p same");
                            tmp_g_clusters[c]->SetMarkerColor(tcol);
                            tmp_g_clusters[c]->SetMarkerStyle(20);
                            tmp_g_clusters[c]->SetMarkerSize(0.75);
                        }
                    }
                    g_clusters[i] = tmp_g_clusters;
                }

                //********** Some Error Checking ********************//

                /*for(int i=0; i<3; i++){

                    std::cout<<"Plane "<<i<<" Vertex pts "<<g_vertex[i]->GetN()<<std::endl;
                    for(size_t s=0; s< pts_shr.size(); s++){
                        std::cout<<"Plane "<<i<<" Shower "<<s<<" pts "<<pts_shr[s][i]->GetN()<<std::endl;
                    }
                    for(size_t t=0;t < pts_trk.size(); t++){
                        std::cout<<"Plane "<<i<<" Track "<<t<<" pts "<<pts_trk[t][i]->GetN()<<std::endl;
                    }
                    for(int c=0; c< num_clusters[i]+1; c++){
                        std::cout<<"Plane "<<i<<" Cluster "<<c<<" pts "<<g_clusters[i][c]->GetN()<<std::endl;
                    }
                }*/



                std::cout<<"Done Plotting clusters"<<std::endl;
                can->Update();
                //can->Write();
                can->SaveAs((print_name+".pdf").c_str(),"pdf");
                //f->Close();
                std::cout<<"PRINTING"<<std::endl;
                //bool_make_sss_plots=false;
                delete can;





            }

        }else{

            std::cout<<" No Neutrino Slice Found/No Showers & tracks"<<std::endl;

        }


        return;
    }
}

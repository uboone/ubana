#include "SinglePhoton_module.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
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

        std::vector< std::map<size_t, std::vector<art::Ptr<recob::Hit>>>> v_newClusterToHitsMap(3);//one for each plane

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


            if(bool_make_sss_plots && showers.size()==1 && tracks.size()==1){

                //TFile *f = new TFile("t.root","recreate");
                //f->cd();

                std::string print_name = "sss_"+std::to_string(m_run_number)+"_"+std::to_string(m_subrun_number)+"_"+std::to_string(m_event_number);
                TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,1600);
                can->Divide(4,2,0,0.1);

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
                std::vector<std::vector<std::vector<double>>> pts_to_recluster(3); //plane, point, {wire,tic}
            //make a map tp actual hits here I guess.
            std::vector<std::map<int, art::Ptr<recob::Hit>>> mapPointIndexToHit(3);

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
                    mapPointIndexToHit[(int)h->View()][pts_to_recluster[(int)h->View()].size()-1] = h;
                }

                g_unass[i] = new TGraph(vec_c.size(), &vec_c[0], &vec_t[0]);


            }
            //Plotting now
            double plot_point_size=0.6;        

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
            
            std::vector<double> vertex_time(3); 
            std::vector<double> vertex_wire(3); 
            
            
            std::vector<TGraph*> g_vertex(3);
            for(int i=0; i<3; i++){
                TPad * pader = (TPad*)can->cd(i+1);

                if(i==0 || i ==4) pader->SetLeftMargin(0.1);

                std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
                std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, *theDetector)};

                vertex_time[i] = time[0];
                vertex_wire[i] = wire[0];

                chan_max[i] = std::max( chan_max[i],wire[0]);
                chan_min[i] = std::min( chan_min[i],wire[0]);

                g_vertex[i] = new TGraph(1,&wire[0],&time[0]);
                g_vertex[i]->SetMarkerStyle(29);
                g_vertex[i]->SetMarkerSize(4);
                g_vertex[i]->SetMarkerColor(kMagenta-3);
                g_vertex[i]->GetYaxis()->SetRangeUser(tick_min*0.98,tick_max*1.02);
                g_vertex[i]->GetXaxis()->SetLimits(chan_min[i]*0.98,chan_max[i]*1.02);
                g_vertex[i]->SetTitle(("Plane " +std::to_string(i)).c_str());
                g_vertex[i]->GetYaxis()->SetTitle("Peak Hit Time Tick");
                g_vertex[i]->GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
                g_vertex[i]->Draw("ap");

                if(i>0){
                    g_vertex[i]->GetYaxis()->SetLabelOffset(999);
                    g_vertex[i]->GetYaxis()->SetLabelSize(0);
                }


                can->cd(i+5);
                g_vertex[i]->Draw("ap");

            }

            //******************************** DeadWireRegions********************************************
            for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                int badchan = bad_channel_list_fixed_mcc9[i].first;                                       
                int ok = bad_channel_list_fixed_mcc9[i].second;       

                if(ok>1)continue;
                auto hs = geom->ChannelToWire(badchan);

                //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 
                int thisp = (int)hs[0].Plane;
                double bc = hs[0].Wire;

                //                        std::cout<<"WIRE "<<thisp<<" "<<bc<<" "<<hs.size()<<std::endl;

                if(chan_min[thisp]*0.98 < bc && bc < chan_max[thisp]*1.02 ){
                    can->cd(thisp+1);
                    TLine *l = new TLine(bc,tick_min*0.98,bc,tick_max*1.02);
                    l->SetLineColor(kGray+1);
                    l->Draw("same");
                    can->cd(thisp+5);
                    l->Draw("same");
                }
            }

            for(size_t t=0; t< pts_trk.size(); t++){
                int tcol = tcols[used_col];
                used_col++;

                for(int i=0; i<3; i++){
                    can->cd(i+1);
                    if(pts_trk[t][i]->GetN()>0){//need a check in case this track has no hits on this plane.
                        pts_trk[t][i]->Draw("p same"); 
                        pts_trk[t][i]->SetMarkerColor(tcol);
                        pts_trk[t][i]->SetFillColor(tcol);
                        pts_trk[t][i]->SetMarkerStyle(20);
                        pts_trk[t][i]->SetMarkerSize(plot_point_size);
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
                        pts_shr[t][i]->SetFillColor(tcol);
                        pts_shr[t][i]->SetMarkerStyle(20);
                        pts_shr[t][i]->SetMarkerSize(plot_point_size);
                    }
                }
            }


            for(int i=0; i<3; i++){
                can->cd(i+1);

                if(g_unass[i]->GetN()>0){
                    g_unass[i]->Draw("p same");
                    g_unass[i]->SetMarkerColor(kBlack);
                    g_unass[i]->SetMarkerStyle(20);
                    g_unass[i]->SetMarkerSize(plot_point_size);
                }
                g_vertex[i]->Draw("p same");
            }






            //*****************************DBSCAN***********************************
            int min_pts = 10;
            double eps = 5.0;
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


            //Create final maps;
            for(size_t i=0; i<3; i++){
                for(size_t p=0; p<pts_to_recluster[i].size(); p++){

                    art::Ptr<recob::Hit> h = mapPointIndexToHit[i].at(p);// Get the hit
                    size_t cluster_label = cluster_labels[i][p];//get the cluster index, 0 = noise                                                

                    //std::cout<<i<<" "<<p<<" "<<cluster_label<<" "<<v_newClusterToHitsMap[i].count(cluster_label)<<std::endl;

                    if(v_newClusterToHitsMap[i].count(cluster_label)<1){
                        std::vector<art::Ptr<recob::Hit>> t_hs= {h};
                        v_newClusterToHitsMap[i][cluster_label] = t_hs;
                    }else{
                        v_newClusterToHitsMap[i].at(cluster_label).push_back(h);//add it to list
                    }
                }
            }
        
            
            //for(size_t i=0; i<3; i++){
              //  for(int c=0; c<num_clusters[i]+1; c++){
                    //auto v_newClusterToHitsMap[i][c];
               // }
           // }


            int max_n_clusters = std::max(num_clusters[0]+1, std::max(num_clusters[1]+1,num_clusters[2]+1))+2;
            std::vector<int> cluster_colors(max_n_clusters,0);
            std::vector<int> base_col = {632,416, 600, 400, 616,  432,  800, 820,  840, 860,  880, 900};

            for(int j=0; j< max_n_clusters; j++){
                int b = (int)rangen->Uniform(0,11);
                int mod = (int)rangen->Uniform(-10,+3);

                cluster_colors[j] = base_col[b]+mod;
            }
            int c_offset = 0;

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
                    if(c>0) tcol = cluster_colors[c+c_offset];
                    tmp_g_clusters[c] = new TGraph(vec_wire[c].size(),&(vec_wire[c])[0],&(vec_time[c])[0] );
                    can->cd(i+5);
                    if(
                            tmp_g_clusters[c]->GetN()>0){
                        tmp_g_clusters[c]->Draw("p same");
                        tmp_g_clusters[c]->SetMarkerColor(tcol);
                        tmp_g_clusters[c]->SetFillColor(tcol);
                        tmp_g_clusters[c]->SetMarkerStyle(20);
                        tmp_g_clusters[c]->SetMarkerSize(plot_point_size);
                    }
                }
                g_clusters[i] = tmp_g_clusters;
                c_offset += num_clusters[i];
            }

            //******************* INFO Plotting *******************************

            TPad *p_top_info = (TPad*)can->cd(4);
            p_top_info->cd();

            TLatex pottex;
            pottex.SetTextSize(0.045);
            pottex.SetTextAlign(13);  //align at top
            pottex.SetNDC();
            std::string pot_draw = "Run: "+std::to_string(m_run_number)+" SubRun: "+std::to_string(m_subrun_number)+" Event: "+std::to_string(m_event_number);
            pottex.DrawLatex(.1,.94, pot_draw.c_str());

            TLegend * l_top = new TLegend(0.5,0.5,0.85,0.85);

            for(size_t t=0; t< pts_shr.size(); t++){
                std::string sname = "Shower "+std::to_string(t);
                if(pts_shr[t][0]->GetN()>0){
                    l_top->AddEntry(pts_shr[t][0],sname.c_str(),"f");
                }else if(pts_shr[t][1]->GetN()>0){
                    l_top->AddEntry(pts_shr[t][1],sname.c_str(),"f");
                }else if(pts_shr[t][2]->GetN()>0){
                    l_top->AddEntry(pts_shr[t][2],sname.c_str(),"f");
                }
            }

            for(size_t t=0; t< pts_trk.size(); t++){
                std::string sname = "Track "+std::to_string(t);
                if(pts_trk[t][0]->GetN()>0){
                    l_top->AddEntry(pts_trk[t][0],sname.c_str(),"f");
                }else if(pts_trk[t][1]->GetN()>0){
                    l_top->AddEntry(pts_trk[t][1],sname.c_str(),"f");
                }else if(pts_trk[t][2]->GetN()>0){
                    l_top->AddEntry(pts_trk[t][2],sname.c_str(),"f");
                }
            }
            l_top->SetLineWidth(0);
            l_top->SetLineColor(kWhite);
            l_top->Draw("same");

            //Clusters

            can->cd(8);
            for(int i=0; i<3; i++){
                TLegend * l_bot = new TLegend(0.1+i*0.25,0.1,0.1+i*0.25+0.25,0.89);

                for(int c=0; c< num_clusters[i]+1; c++){
                    if(c==0)continue;

                    int num_hits_in_cluster = v_newClusterToHitsMap[i][c].size();
                    auto hitz = v_newClusterToHitsMap[i][c];
                    auto ssscorz = this->ScoreCluster(i,c, hitz ,vertex_wire[i], vertex_time[i], showers[0]);


                    std::string sname = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Hits: "+std::to_string(num_hits_in_cluster)+"}{#splitline{PCA "+std::to_string(ssscorz.pca_0)+"}{Theta:" +std::to_string(ssscorz.pca_theta)+"}}}";
                    l_bot->AddEntry(g_clusters[i][c],sname.c_str(),"f");

                }
                l_bot->SetLineWidth(0);
                l_bot->SetLineColor(kWhite);
                l_bot->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
                l_bot->Draw("same");

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




 sss_score SinglePhoton::ScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower){
            sss_score score(p,cl);
            score.n_hits = hits.size();
            std::cout<<"SSS || Starting Score Calcultaion "<<std::endl;

            std::vector<double> t_wires;
            std::vector<double> t_ticks;

            // ************* Some simple metrics relative to study point (usually vertex) ***************
            score.max_dist_tick = 0;
            score.min_dist_tick = 1e10;
            score.mean_dist_tick = 0;

            score.max_dist_wire = 0;
            score.min_dist_wire = 1e10;
            score.mean_dist_wire = 0;
           


            for(auto &h: hits){
                double h_tick = (double)h->PeakTime();
                double h_wire = (double)h->WireID().Wire;

                score.max_dist_tick = std::max(score.max_dist_tick, fabs(h_tick-vertex_tick));
                score.min_dist_tick = std::min(score.min_dist_tick, fabs(h_tick-vertex_tick));

                score.max_dist_wire = std::max(score.max_dist_wire, fabs(h_wire-vertex_wire));
                score.min_dist_wire = std::min(score.min_dist_wire, fabs(h_wire-vertex_wire));

                score.mean_dist_tick += fabs(h_tick-vertex_tick);
                score.mean_dist_wire += fabs(h_wire-vertex_wire);

                t_wires.push_back(h_wire);
                t_ticks.push_back(h_tick);
            }

//            TGraph * g_pts = new TGraph(t_wires.size(),&t_ticks[0],&t_wires[0]);

            score.mean_dist_tick = score.mean_dist_tick/(double)score.n_hits;
            score.mean_dist_wire = score.mean_dist_wire/(double)score.n_hits;

            // **************** Metrics of Pointing: Does this cluster "point" back to the vertex? *************************
            // **************** First off, PCA
        
            TPrincipal* principal = new TPrincipal(2,"D");

            for(int i = 0; i < score.n_hits; i++){
                std::vector<double> tmp_pts = {t_wires[i], t_ticks[i]};
                principal->AddRow(&tmp_pts[0]);
            }
            principal->MakePrincipals();
            principal->Print();

            TVectorD * eigenval = (TVectorD*) principal->GetEigenValues();
            TMatrixD * eigenvec = (TMatrixD*) principal->GetEigenVectors();
            TMatrixD * covar = (TMatrixD*) principal->GetCovarianceMatrix();

            score.pca_0 = (*eigenval)(0);
            score.pca_1 = (*eigenval)(1);

            (*eigenvec).Print();
            (*covar).Print();
            std::cout<<" SSS | Eigen: "<<score.pca_0<<" "<<score.pca_1<<std::endl;

            score.pca_theta = atan((*covar)[0][0]/(*covar)[0][1])*180.0/3.14159;

            delete principal;

    

            return score;
 }




}

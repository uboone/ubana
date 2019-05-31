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
#include "TEllipse.h"
namespace single_photon
{



    void SinglePhoton::ClearSecondShowers(){
        m_sss_num_unassociated_hits=0;
        m_sss_num_associated_hits=0;

        m_sss_num_candidates = 0;

        m_sss_candidate_num_hits.clear();
        m_sss_candidate_num_wires.clear();
        m_sss_candidate_num_ticks.clear();
        m_sss_candidate_plane.clear();
        m_sss_candidate_PCA.clear();
        m_sss_candidate_impact_parameter.clear();
        m_sss_candidate_fit_slope.clear();
        m_sss_candidate_fit_constant.clear();
        m_sss_candidate_mean_tick.clear();
        m_sss_candidate_max_tick.clear();
        m_sss_candidate_min_tick.clear();
        m_sss_candidate_min_wire.clear();
        m_sss_candidate_max_wire.clear();
        m_sss_candidate_mean_wire.clear();
        m_sss_candidate_min_dist.clear();
    }

    void SinglePhoton::ResizeSecondShowers(size_t size){

    }


    void SinglePhoton::CreateSecondShowerBranches(){
        vertex_tree->Branch("sss_num_unassociated_hits",&m_sss_num_unassociated_hits,"sss_num_unassociated_hits/I");
        vertex_tree->Branch("sss_num_associated_hits",&m_sss_num_associated_hits,"sss_num_associated_hits/I");

        vertex_tree->Branch("sss_num_candidates",&m_sss_num_candidates,"sss_num_candidates/I");
        vertex_tree->Branch("sss_candidate_num_hits",&m_sss_candidate_num_hits);
        vertex_tree->Branch("sss_candidate_num_wires",&m_sss_candidate_num_wires);
        vertex_tree->Branch("sss_candidate_num_ticks",&m_sss_candidate_num_ticks);
        vertex_tree->Branch("sss_candidate_plane",&m_sss_candidate_plane);
        vertex_tree->Branch("sss_candidate_PCA",&m_sss_candidate_PCA);
        vertex_tree->Branch("sss_candidate_impact_parameter",&m_sss_candidate_impact_parameter); 
        vertex_tree->Branch("sss_candidate_fit_slope",&m_sss_candidate_fit_slope);
        vertex_tree->Branch("sss_candidate_fit_constant",&m_sss_candidate_fit_constant);
        vertex_tree->Branch("sss_candidate_mean_tick",&m_sss_candidate_mean_tick);
        vertex_tree->Branch("sss_candidate_max_tick",&m_sss_candidate_max_tick);
        vertex_tree->Branch("sss_candidate_min_tick",&m_sss_candidate_min_tick);
        vertex_tree->Branch("sss_candidate_mean_wire",&m_sss_candidate_mean_wire);
        vertex_tree->Branch("sss_candidate_max_wire",&m_sss_candidate_max_wire);
        vertex_tree->Branch("sss_candidate_min_wire",&m_sss_candidate_min_wire);
        vertex_tree->Branch("sss_candidate_min_dist",&m_sss_candidate_min_dist);



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
                TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,2400);
                can->Divide(4,3,0,0.1);

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
                    tcols.push_back(tcols[(int)rangen->Uniform(0,7)]+(int)rangen->Uniform(-5,5));
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

                if(i==0 || i ==4 || i == 8) pader->SetLeftMargin(0.1);

                std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
                std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, *theDetector)};

                vertex_time[i] = time[0];
                vertex_wire[i] = wire[0];

                if(i==0) m_vertex_pos_wire_p0 = wire[0];
                if(i==1) m_vertex_pos_wire_p1 = wire[0];
                if(i==2) m_vertex_pos_wire_p2 = wire[0];
                m_vertex_pos_tick = time[0];

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

                can->cd(i+9);
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
                    can->cd(thisp+9);
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
                can->cd(i+5);
                g_vertex[i]->Draw("p same");
                can->cd(i+9);
                g_vertex[i]->Draw("p same");

                double rad_cm = 12.0;
                TEllipse * ell_p = new TEllipse(vertex_wire[i],vertex_time[i],rad_cm/0.3,rad_cm*25);
                ell_p->SetLineColor(kRed);
                ell_p->SetFillStyle(0);
                ell_p->Draw("same");

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
                
                TLegend * l_bot2 = new TLegend(0.1+i*0.25,0.1,0.1+i*0.25+0.25,0.89);

                for(int c=0; c< num_clusters[i]+1; c++){
                    if(c==0)continue;

                    int num_hits_in_cluster = v_newClusterToHitsMap[i][c].size();
                    auto hitz = v_newClusterToHitsMap[i][c];
                    auto ssscorz = this->ScoreCluster(i,c, hitz ,vertex_wire[i], vertex_time[i], showers[0]);
                    int is_in_shower = this->CompareToShowers(i,c, hitz ,vertex_wire[i], vertex_time[i], showers, showerToPFParticleMap, pfParticleToHitsMap,eps);

                    // std::string sname = makeSplitlineString({"Cluster: ","Hits: ","PCA: ","Theta: "},{c,num_hits_in_cluster});

                    std::string sname = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Hits: "+std::to_string(num_hits_in_cluster)+"}{#splitline{PCA "+std::to_string(ssscorz.pca_0)+"}{#splitline{Theta:" +std::to_string(ssscorz.pca_theta)+"}{#splitline{Wires: "+std::to_string(ssscorz.n_wires)+ "}{#splitline{Ticks: "+std::to_string(ssscorz.n_ticks)+"}{#splitline{ReMerged: "+std::to_string(is_in_shower)+"}{}}}}}}}";
                    l_bot->AddEntry(g_clusters[i][c],sname.c_str(),"f");


                    //Here we will only plot those that pass in bottom:
                    if(ssscorz.pass && is_in_shower ==-1){
                        can->cd(i+9);

                        if(g_clusters[i][c]->GetN()>0){
                            TGraph * tmp = (TGraph*)g_clusters[i][c]->Clone(("tmp_"+std::to_string(i)+std::to_string(c)).c_str());


                            int Npts = 20;
                            TGraph * core  = (TGraph*)this->GetNearestNpts(i,c,hitz,vertex_wire[i],vertex_time[i],Npts);

                            core->Draw("p same");
                            tmp->Draw("p same");
                            core->Fit("pol1","Q","same",chan_min[i],chan_max[i]);
                            core->GetFunction("pol1")->SetLineWidth(1); 
                            core->GetFunction("pol1")->SetLineStyle(3); 
                            core->GetFunction("pol1")->SetLineColor(g_clusters[i][c]->GetMarkerColor()); 
                            double con = core->GetFunction("pol1")->GetParameter(0);
                            double slope = core->GetFunction("pol1")->GetParameter(1);

                            //lets map (wire,tick) to a rudamentary (cm,cm);
                            //double slope2 = slope*25*0.3;
                            //double con2 = con*25;
                                
                            double impact_parameter = 1e10;// fabs(slope*vertex_wire[i] +vertex_time[i]+con)/sqrt(slope*slope+1.0*1.0);

                            //rudimentary!
                            for(double k=chan_min[i]; k< chan_max[i];k++){
                                double y = slope*k+con;
                                double dist = sqrt(pow(k*0.3-vertex_wire[i]*0.3,2)+pow(y/25.0-vertex_time[i]/25.0,2));
                                impact_parameter = std::min(impact_parameter,dist);
                            }


                            m_sss_num_candidates++;

                            m_sss_candidate_num_hits.push_back(num_hits_in_cluster);
                            m_sss_candidate_num_wires.push_back((int)ssscorz.n_wires);
                            m_sss_candidate_num_ticks.push_back((int)ssscorz.n_ticks);
                            m_sss_candidate_plane.push_back((int)i);
                            m_sss_candidate_PCA.push_back(ssscorz.pca_0);
                            m_sss_candidate_impact_parameter.push_back(impact_parameter);
                            m_sss_candidate_fit_slope.push_back(slope);
                            m_sss_candidate_fit_constant.push_back(con);
                            m_sss_candidate_mean_tick.push_back(ssscorz.mean_tick);
                            m_sss_candidate_max_tick.push_back(ssscorz.max_tick);
                            m_sss_candidate_min_tick.push_back(ssscorz.min_tick);
                            m_sss_candidate_min_wire.push_back(ssscorz.min_wire);
                            m_sss_candidate_max_wire.push_back(ssscorz.max_wire);
                            m_sss_candidate_mean_wire.push_back(ssscorz.mean_wire);
                            m_sss_candidate_min_dist.push_back(ssscorz.min_dist);
                            
                            std::string sname2 = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Impact: "+std::to_string(impact_parameter)+"}{MinDist: "+std::to_string(ssscorz.min_dist)+"}}";
                            l_bot2->AddEntry(tmp,sname2.c_str(),"f");
                        }
                    }          

                }
                can->cd(8);
                l_bot->SetLineWidth(0);
                l_bot->SetLineColor(kWhite);
                l_bot->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
                l_bot->Draw("same");

                can->cd(12);
                l_bot2->SetLineWidth(0);
                l_bot2->SetLineColor(kWhite);
                l_bot2->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
                l_bot2->Draw("same");


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



    TGraph* SinglePhoton::GetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts){

         std::vector<double>t_wire;
         std::vector<double>t_tick;
        // std::vector<double>t_dist;

         std::vector<double>all_wire;
         std::vector<double>all_tick;
         std::vector<double>all_dist;


        for(size_t h = 0; h< hitz.size(); h++){
                auto hit = hitz[h];
                double h_wire = (double)hit->WireID().Wire;
                double h_tick = (double)hit->PeakTime();
 
                double dd =sqrt(pow(h_wire*0.3-vertex_wire*0.3,2)+pow(h_tick/25.0- vertex_tick/25.0,2));
                all_wire.push_back(h_wire);   
                all_tick.push_back(h_tick);   
                all_dist.push_back(dd);
        }

        std::vector<size_t> sorted_in = sort_indexes(all_dist);
        size_t max_e = std::min((size_t)Npts,hitz.size());

        for(size_t i =0; i<max_e; i++){
            t_wire.push_back(all_wire[sorted_in[hitz.size()-1-i]]);
            t_tick.push_back(all_tick[sorted_in[hitz.size()-1-i]]);
        }

        return new TGraph(t_wire.size(),&t_wire[0],&t_tick[0]);
    }

    sss_score SinglePhoton::ScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower){
        sss_score score(p,cl);
        score.n_hits = hits.size();
        std::cout<<"SSS || Starting Score Calcultaion "<<std::endl;

        std::vector<double> t_wires;
        std::vector<double> t_ticks;

        // 
        int n_min_ticks = 4;
        int n_min_wires = 3;
        double n_max_pca = 0.9999;

        score.pass = true;

        // ************* Some simple metrics relative to study point (usually vertex) ***************
        score.max_dist_tick = 0;
        score.min_dist_tick = 1e10;
        score.mean_dist_tick = 0;

        score.max_dist_wire = 0;
        score.min_dist_wire = 1e10;
        score.mean_dist_wire = 0;

        score.max_dist = 0;
        score.min_dist = 1e10;
        score.mean_dist = 0;

        score.mean_tick =0;
        score.max_tick =0;
        score.min_tick =1e10;
        
        score.mean_wire =0;
        score.max_wire =0;
        score.min_wire =1e10;

        score.n_wires = 0;
        score.n_ticks = 0;

        score.impact_parameter = -99;

        std::map<int,bool> wire_count;
        std::map<int,bool> tick_count;

        for(auto &h: hits){
            double h_tick = (double)h->PeakTime();
            double h_wire = (double)h->WireID().Wire;

            score.mean_wire += h_wire;
            score.mean_tick += h_tick;

            score.max_wire = std::max(score.max_wire, h_wire);
            score.min_wire = std::min(score.min_wire, h_wire);

            score.max_tick = std::max(score.max_tick, h_tick);
            score.min_tick = std::min(score.min_tick, h_tick);

            score.max_dist_tick = std::max(score.max_dist_tick, fabs(h_tick-vertex_tick));
            score.min_dist_tick = std::min(score.min_dist_tick, fabs(h_tick-vertex_tick));

            score.max_dist_wire = std::max(score.max_dist_wire, fabs(h_wire-vertex_wire));
            score.min_dist_wire = std::min(score.min_dist_wire, fabs(h_wire-vertex_wire));

            score.mean_dist_tick += fabs(h_tick-vertex_tick);
            score.mean_dist_wire += fabs(h_wire-vertex_wire);
            
            //wierd dits
            double dd =sqrt(pow(h_wire*0.3-vertex_wire*0.3,2)+pow(h_tick/25.0- vertex_tick/25.0,2));
            score.mean_dist += dd;
            score.max_dist = std::max(dd,score.max_dist);
            score.min_dist = std::min(dd,score.min_dist);



            t_wires.push_back(h_wire);
            t_ticks.push_back(h_tick);

            if(wire_count.count((int)h_wire)<1){
                wire_count[((int)h_wire)] = true;
                score.n_wires++;
            }
            if(tick_count.count((int)h_tick)<1){
                tick_count[((int)h_tick)] = true;
                score.n_ticks++;
            }

        }

        //            TGraph * g_pts = new TGraph(t_wires.size(),&t_ticks[0],&t_wires[0]);

        score.mean_tick = score.mean_tick/(double)score.n_hits;
        score.mean_wire = score.mean_wire/(double)score.n_hits;
    
        score.mean_dist = score.mean_dist/(double)score.n_hits;

        score.mean_dist_tick = score.mean_dist_tick/(double)score.n_hits;
        score.mean_dist_wire = score.mean_dist_wire/(double)score.n_hits;

        // **************** Metrics of Pointing: Does this cluster "point" back to the vertex? *************************
        // **************** First off, PCA

        TPrincipal* principal = new TPrincipal(2,"D");
        double mod_wire = 1.0;
        double mod_tick = 1.0;

        for(int i = 0; i < score.n_hits; i++){
            std::vector<double> tmp_pts = {t_wires[i]*mod_wire, t_ticks[i]/mod_tick};
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


        double slope = ((*covar)[0][0]/(*covar)[0][1]);
        double c = score.mean_tick*mod_wire - slope*score.mean_wire/mod_tick;
        score.impact_parameter = fabs(slope*vertex_wire*mod_wire +vertex_tick/mod_tick+c)/sqrt(slope*slope+1.0*1.0);


        if(score.n_wires < n_min_wires || score.n_ticks < n_min_ticks || score.pca_0 >= n_max_pca){

            score.pass = false;
        }



        delete principal;

        return score;
    }

    int SinglePhoton::CompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,
            const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      const   std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,                    double eps){


        for(size_t s =0; s< showers.size(); s++){
            art::Ptr<recob::Shower> shower = showers[s];
            art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(shower);
            std::vector<art::Ptr<recob::Hit>> showerhits = pfParticleToHitsMap.at(pfp);

            bool in_primary_shower = false;
            for(size_t h = 0; h< hitz.size(); h++){
                auto hit = hitz[h];
                double h_wire = (double)hit->WireID().Wire;
                double h_tick = (double)hit->PeakTime();


                for(auto &sh: showerhits){

                    if(sh->View() != hit->View()) continue;

                    double sh_wire = (double)sh->WireID().Wire;
                    double sh_tick = (double)sh->PeakTime();


                    double dist = sqrt(pow(sh_wire*0.3-h_wire*0.3,2)+pow(sh_tick/25.0-h_tick/25.0,2));

                    if(dist<=eps){
                        in_primary_shower = true;
                        return (int)s;
                    }

                }

            }

            if(in_primary_shower){
                return (int)s;
            }
        }


        return -1;
    }

}

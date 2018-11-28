#include "SinglePhoton_module.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TruncMean.h"


namespace single_photon
{


    void SinglePhoton::ClearTracks(){
        m_reco_asso_tracks=0;
        m_reco_track_length.clear();
        m_reco_track_dirx.clear();
        m_reco_track_diry.clear();
        m_reco_track_dirz.clear();
        m_reco_track_startx.clear();
        m_reco_track_starty.clear();
        m_reco_track_startz.clear();
        m_reco_track_endx.clear();
        m_reco_track_endy.clear();
        m_reco_track_endz.clear();

        m_reco_track_theta_yz.clear();
        m_reco_track_phi_yx.clear();

        m_reco_track_num_trajpoints.clear();
        m_reco_track_num_spacepoints.clear();
        m_reco_track_proton_kinetic_energy.clear();

        m_reco_track_spacepoint_principal0.clear();
        m_reco_track_spacepoint_principal1.clear();
        m_reco_track_spacepoint_principal2.clear();

        m_reco_track_spacepoint_chi.clear();
        m_reco_track_spacepoint_max_dist.clear();

        m_reco_track_mean_dEdx.clear();
        m_reco_track_mean_dEdx_start_half.clear();
        m_reco_track_mean_dEdx_end_half.clear();
        m_reco_track_good_calo.clear();
        m_reco_track_mean_trunc_dEdx.clear();
        m_reco_track_mean_trunc_dEdx_start_half.clear();
        m_reco_track_mean_trunc_dEdx_end_half.clear();
        m_reco_track_trunc_PIDA.clear();

        m_sim_track_matched.clear();
        m_sim_track_energy.clear();
        m_sim_track_kinetic_energy.clear();
        m_sim_track_mass.clear();
        m_sim_track_pdg.clear();
        m_sim_track_origin.clear();
        m_sim_track_process.clear();
        m_sim_track_startx.clear();
        m_sim_track_starty.clear();
        m_sim_track_startz.clear();

        m_reco_track_pid_bragg_likelihood_plane2.clear();
        m_reco_track_pid_pida_plane2.clear();
        m_reco_track_pid_chi_plane2.clear();
        m_reco_track_end_to_nearest_dead_wire_plane0.clear();
        m_reco_track_end_to_nearest_dead_wire_plane1.clear();
        m_reco_track_end_to_nearest_dead_wire_plane2.clear();

    }

    void SinglePhoton::ResizeTracks(size_t size){
        m_reco_track_length.resize(size);
        m_reco_track_dirx.resize(size);
        m_reco_track_diry.resize(size);
        m_reco_track_dirz.resize(size);
        m_reco_track_endx.resize(size);
        m_reco_track_endy.resize(size);
        m_reco_track_endz.resize(size);

        m_reco_track_startx.resize(size);
        m_reco_track_starty.resize(size);
        m_reco_track_startz.resize(size);
        m_reco_track_num_trajpoints.resize(size);
        m_reco_track_num_spacepoints.resize(size);
        m_reco_track_proton_kinetic_energy.resize(size);


        m_reco_track_spacepoint_principal0.resize(size);
        m_reco_track_spacepoint_principal1.resize(size);
        m_reco_track_spacepoint_principal2.resize(size);

        m_reco_track_spacepoint_chi.resize(size);
        m_reco_track_spacepoint_max_dist.resize(size);

        m_reco_track_theta_yz.resize(size);
        m_reco_track_phi_yx.resize(size);

        m_reco_track_mean_dEdx.resize(size);
        m_reco_track_mean_dEdx_start_half.resize(size);
        m_reco_track_mean_dEdx_end_half.resize(size);
        m_reco_track_good_calo.resize(size);
        m_reco_track_mean_trunc_dEdx.resize(size);
        m_reco_track_mean_trunc_dEdx_start_half.resize(size);
        m_reco_track_mean_trunc_dEdx_end_half.resize(size);
        m_reco_track_trunc_PIDA.resize(size);

        m_sim_track_matched.resize(size);
        m_sim_track_energy.resize(size);
        m_sim_track_mass.resize(size);
        m_sim_track_kinetic_energy.resize(size);
        m_sim_track_pdg.resize(size);
        m_sim_track_origin.resize(size);
        m_sim_track_process.resize(size);
        m_sim_track_startx.resize(size);
        m_sim_track_starty.resize(size);
        m_sim_track_startz.resize(size);

        m_reco_track_pid_bragg_likelihood_plane2.resize(size);
        m_reco_track_pid_pida_plane2.resize(size);
        m_reco_track_pid_chi_plane2.resize(size);


        m_reco_track_end_to_nearest_dead_wire_plane0.resize(size);
        m_reco_track_end_to_nearest_dead_wire_plane1.resize(size);
        m_reco_track_end_to_nearest_dead_wire_plane2.resize(size);
    }

    void SinglePhoton::CreateTrackBranches(){
        vertex_tree->Branch("reco_asso_tracks",&m_reco_asso_tracks,"reco_asso_tracks/I");
        vertex_tree->Branch("reco_track_displacement", &m_reco_track_length);
        vertex_tree->Branch("reco_track_dirx", &m_reco_track_dirx);
        vertex_tree->Branch("reco_track_diry", &m_reco_track_diry);
        vertex_tree->Branch("reco_track_dirz", &m_reco_track_dirz);
        vertex_tree->Branch("reco_track_startx", &m_reco_track_startx);
        vertex_tree->Branch("reco_track_starty", &m_reco_track_starty);
        vertex_tree->Branch("reco_track_startz", &m_reco_track_startz);
        vertex_tree->Branch("reco_track_endx", &m_reco_track_endx);
        vertex_tree->Branch("reco_track_endy", &m_reco_track_endy);
        vertex_tree->Branch("reco_track_endz", &m_reco_track_endz);

        vertex_tree->Branch("reco_track_theta_yz", &m_reco_track_theta_yz);
        vertex_tree->Branch("reco_track_phi_yx", &m_reco_track_phi_yx);

        vertex_tree->Branch("reco_track_num_trajpoints", &m_reco_track_num_trajpoints);
        vertex_tree->Branch("reco_track_num_spacepoints", &m_reco_track_num_spacepoints);
        vertex_tree->Branch("reco_track_proton_kinetic_energy", &m_reco_track_proton_kinetic_energy);

        vertex_tree->Branch("reco_track_spacepoint_principal0",&m_reco_track_spacepoint_principal0);
        vertex_tree->Branch("reco_track_spacepoint_principal1",&m_reco_track_spacepoint_principal1);
        vertex_tree->Branch("reco_track_spacepoint_principal2",&m_reco_track_spacepoint_principal2);



        vertex_tree->Branch("reco_track_spacepoint_chi",&m_reco_track_spacepoint_chi);
        vertex_tree->Branch("reco_track_spacepoint_max_dist",&m_reco_track_spacepoint_max_dist);

        vertex_tree->Branch("reco_track_mean_dEdx",&m_reco_track_mean_dEdx);
        vertex_tree->Branch("reco_track_mean_dEdx_start_half",&m_reco_track_mean_dEdx_end_half);
        vertex_tree->Branch("reco_track_mean_dEdx_end_half",&m_reco_track_mean_dEdx_start_half);
        vertex_tree->Branch("reco_track_good_calo",&m_reco_track_good_calo);
        vertex_tree->Branch("reco_track_mean_trunc_dEdx",&m_reco_track_mean_trunc_dEdx);
        vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half",&m_reco_track_mean_trunc_dEdx_end_half);
        vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half",&m_reco_track_mean_trunc_dEdx_start_half);
        vertex_tree->Branch("reco_track_trunc_PIDA",&m_reco_track_trunc_PIDA);

        vertex_tree->Branch("reco_track_pid_bragg_likelihood_plane2",&m_reco_track_pid_bragg_likelihood_plane2);
        vertex_tree->Branch("reco_track_pid_pida_plane2",&m_reco_track_pid_pida_plane2);
        vertex_tree->Branch("reco_track_pid_chi_plane2",&m_reco_track_pid_chi_plane2);

        vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane0",&m_reco_track_end_to_nearest_dead_wire_plane0);
        vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane1",&m_reco_track_end_to_nearest_dead_wire_plane1);
        vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane2",&m_reco_track_end_to_nearest_dead_wire_plane2);

        vertex_tree->Branch("sim_track_matched",&m_sim_track_matched);
        vertex_tree->Branch("sim_track_energy",&m_sim_track_energy);
        vertex_tree->Branch("sim_track_kinetic_energy",&m_sim_track_kinetic_energy);
        vertex_tree->Branch("sim_track_mass",&m_sim_track_mass);
        vertex_tree->Branch("sim_track_pdg",&m_sim_track_pdg);
        vertex_tree->Branch("sim_track_origin",&m_sim_track_origin);
        vertex_tree->Branch("sim_track_process",&m_sim_track_process);
        vertex_tree->Branch("sim_track_startx",&m_sim_track_startx);
        vertex_tree->Branch("sim_track_starty",&m_sim_track_starty);
        vertex_tree->Branch("sim_track_startz",&m_sim_track_startz);
    }






    void SinglePhoton::AnalyzeTracks(const std::vector<art::Ptr<recob::Track>>& tracks,
            std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToNuPFParticleMap,
            std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfParticleToSpacePointsMap){


        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Starting recob::Track analysis"<<std::endl;;

        m_reco_asso_tracks = tracks.size();
        int i_trk=0;

        this->ResizeTracks(m_reco_asso_tracks);

        //const double adc2eU(5.1e-3);
        //const double adc2eV(5.2e-3);
        // const double adc2eW(5.4e-3);


        //    const double tau(theDetector->ElectronLifetime());


        for (TrackVector::const_iterator iter = tracks.begin(), iterEnd = tracks.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Track> track = *iter;
            const art::Ptr<recob::PFParticle> pfp = trackToNuPFParticleMap[track];
            const std::vector< art::Ptr<recob::SpacePoint> > trk_spacepoints = pfParticleToSpacePointsMap[pfp];


            int m_trkid = track->ID();
            double m_length = track->Length();
            auto m_trk_dir = track->Direction();

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t On Track: "<<i_trk<<" with TrackID: "<<m_trkid<<" and length: "<<m_length<<""<<std::endl;;

            m_reco_track_num_spacepoints[i_trk] = (int)trk_spacepoints.size();


            m_reco_track_startx[i_trk]= track->Start().X();   
            m_reco_track_starty[i_trk]= track->Start().Y();   
            m_reco_track_startz[i_trk]= track->Start().Z();   

            m_reco_track_length[i_trk] =m_length;
            m_reco_track_dirx[i_trk] = m_trk_dir.first.X();   
            m_reco_track_diry[i_trk] = m_trk_dir.first.Y();   
            m_reco_track_dirz[i_trk] = m_trk_dir.first.Z();   

            m_reco_track_endx[i_trk] = track->End().X();   
            m_reco_track_endy[i_trk]= track->End().Y();   
            m_reco_track_endz[i_trk]= track->End().Z();   

            m_reco_track_theta_yz[i_trk] = atan2(m_reco_track_diry[i_trk],m_reco_track_dirz[i_trk]);
            m_reco_track_phi_yx[i_trk] = atan2(m_reco_track_diry[i_trk],m_reco_track_dirx[i_trk]);

            std::vector<double> tmp_trk_start = {m_reco_track_startx[i_trk],m_reco_track_starty[i_trk],m_reco_track_startz[i_trk]};
            std::vector<double> tmp_trk_end = {m_reco_track_endx[i_trk],m_reco_track_endy[i_trk],m_reco_track_endz[i_trk]};
            double max_dist_from_line = -9999999;

            m_reco_track_spacepoint_chi[i_trk] = 0.0;
            //Principal componant analysis of SPACEPOINTS and not trajectory points. Will add a few things here in future.

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Beginining PCA analysis of track"<<std::endl;;

            TPrincipal* principal = new TPrincipal(3,"ND");
            for(int x = 0; x < m_reco_track_num_spacepoints[i_trk]; x++){
                std::vector<double> tmp_spacepoints = {trk_spacepoints[x]->XYZ()[0],trk_spacepoints[x]->XYZ()[1] , trk_spacepoints[x]->XYZ()[2]};
                principal->AddRow(&tmp_spacepoints[0]);

                double dist = dist_line_point(tmp_trk_start,tmp_trk_end,tmp_spacepoints);
                if(dist> max_dist_from_line) max_dist_from_line = dist;
                m_reco_track_spacepoint_chi[i_trk] += dist*dist;
            }
            m_reco_track_spacepoint_max_dist[i_trk]= max_dist_from_line;

            principal->MakePrincipals();
            TVectorD * eigen = (TVectorD*) principal->GetEigenValues();

            m_reco_track_spacepoint_principal0[i_trk]=(*eigen)(0);
            m_reco_track_spacepoint_principal1[i_trk]=(*eigen)(1);
            m_reco_track_spacepoint_principal2[i_trk]=(*eigen)(2);

            delete principal;

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Completed PCA analysis of track. Primary componant: "<<m_reco_track_spacepoint_principal0.back()<<""<<std::endl;;

            //range based energy calculation assuming
            m_reco_track_proton_kinetic_energy[i_trk] = proton_length2energy_tgraph.Eval(m_length)/1000.0; 
            if(m_length == 0.0) m_reco_track_proton_kinetic_energy[i_trk]=0.0;

            // Dead Wire Approximity
            m_reco_track_end_to_nearest_dead_wire_plane0[i_trk] = distanceToNearestDeadWire(0, m_reco_track_endy[i_trk], m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);
            m_reco_track_end_to_nearest_dead_wire_plane1[i_trk] = distanceToNearestDeadWire(1, m_reco_track_endy[i_trk], m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);
            m_reco_track_end_to_nearest_dead_wire_plane2[i_trk] = distanceToNearestDeadWire(2, m_reco_track_endy[i_trk], m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);


            //A loop over the trajectory points
            size_t const traj_size = track->CountValidPoints();
            m_reco_track_num_trajpoints[i_trk] = (int)traj_size;

            for(unsigned int  p = 0; p < traj_size; ++p) {
                //recob::Track::TrajectoryPoint_t const & trajp = track->TrajectoryPoint(j);
                //recob::Track::Point_t const & pos = trajp.position;
                //recob::Track::Vector_t const & mom = trajp.momentum;

            }


            i_trk++;

        }


        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Finished."<<std::endl;;
    }

    void SinglePhoton::RecoMCTracks(const std::vector<art::Ptr<recob::Track>>& tracks,  
            std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap, 
            std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,
            std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap){


        if(m_is_verbose) std::cout<<"SinglePhoton::RecoMCTracks()\t||\t Begininning recob::Track Reco-MC suite"<<std::endl;;

        int i_trk = 0;
        //for (TrackVector::const_iterator iter = tracks.begin(), iterEnd = tracks.end(); iter != iterEnd; ++iter)
        for(size_t k =0; k< tracks.size();++k)    

        {

         //   const art::Ptr<recob::Track> track = *iter;
            const art::Ptr<recob::Track> track = tracks[k];
            m_sim_track_matched[i_trk] = 0;
            std::cout<<"INSIDE : "<<trackToMCParticleMap.count(track)<<std::endl;
                
            if(trackToMCParticleMap.count(track)>0){
                
                const art::Ptr<simb::MCParticle> mcparticle = trackToMCParticleMap[track];
                std::cout<<"count2: "<<MCParticleToMCTruthMap.count(mcparticle)<<std::endl;
                const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruthMap[mcparticle];

                std::vector<double> corrected(3);
                this->spacecharge_correction(mcparticle, corrected);

                m_sim_track_matched[i_trk] = 1;
                m_sim_track_energy[i_trk] = mcparticle->E();
                m_sim_track_mass[i_trk] = mcparticle->Mass();
                m_sim_track_kinetic_energy[i_trk] = m_sim_track_energy[i_trk]-m_sim_track_mass[i_trk];
                m_sim_track_pdg[i_trk] = mcparticle->PdgCode();
                m_sim_track_process[i_trk] = mcparticle->Process();
                m_sim_track_startx[i_trk] = corrected[0];
                m_sim_track_starty[i_trk] = corrected[1];
                m_sim_track_startz[i_trk] = corrected[2];
                m_sim_track_origin[i_trk] = mctruth->Origin();
            }
            i_trk++;
        }

    }





    void SinglePhoton::AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::map<art::Ptr<recob::Track>,art::Ptr<anab::Calorimetry>> &trackToCaloMap)
    {

        if(m_is_verbose) std::cout<<"SinglePhoton::CollectCalo(recob::Track)\t||\t Starting calo module for recob::Track"<<std::endl;;

        for(size_t i_trk = 0; i_trk<tracks.size(); ++i_trk){
            const art::Ptr<recob::Track>      track = tracks[i_trk];
            const art::Ptr<anab::Calorimetry> calo  = trackToCaloMap[track];

            size_t calo_length = calo->dEdx().size();

            TruncMean tm;
            std::vector<double> trunc_dEdx;
            std::vector<double> res_range_good;
            std::vector<double> dEdx_good;

            m_reco_track_good_calo[i_trk] =  0;
            m_reco_track_mean_dEdx[i_trk] =  0.0;
            m_reco_track_mean_dEdx_start_half[i_trk] =  0.0;
            m_reco_track_mean_dEdx_end_half[i_trk] =  0.0;
            m_reco_track_mean_trunc_dEdx[i_trk] =  0.0;
            m_reco_track_mean_trunc_dEdx_start_half[i_trk] =  0.0;
            m_reco_track_mean_trunc_dEdx_end_half[i_trk] =  0.0;
            m_reco_track_trunc_PIDA[i_trk] =  0.0;


            //First off look over ALL points
            for (size_t k = 0; k < calo_length; ++k) {
                double res_range =    calo->ResidualRange()[k];
                double dEdx =         calo->dEdx()[k];

                m_reco_track_mean_dEdx[i_trk] += dEdx;
                if(k <= calo_length/2){ 
                    m_reco_track_mean_dEdx_start_half[i_trk]+=dEdx;
                }else{
                    m_reco_track_mean_dEdx_end_half[i_trk]+=dEdx;
                }

                bool is_sensible = dEdx < m_track_calo_max_dEdx; 
                bool is_nan =dEdx != dEdx; 
                bool is_inf = std::isinf(dEdx);
                bool is_nonzero = dEdx> m_track_calo_min_dEdx;

                if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length-1){
                    res_range_good.push_back(res_range);
                    dEdx_good.push_back(dEdx);
                }

                //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
            }// End of first loop.

            m_reco_track_good_calo[i_trk] = 0;
            if(res_range_good.size() >= m_track_calo_min_dEdx_hits){
                m_reco_track_good_calo[i_trk] = res_range_good.size();

                //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimum_dx in the res_range
                double tenth_track = std::max(res_range_good.front(), res_range_good.back())/m_track_calo_trunc_fraction;
                double min_dx = 999;
                for(int j = res_range_good.size()-1; j>1; j--){
                    double dx = fabs(res_range_good[j]-res_range_good[j-1]);
                    if(dx < min_dx) min_dx = dx;
                }
                double rad = std::max( min_dx*2, tenth_track); 

                //Calculate the residual range
                tm.setRadius(rad);
                tm.CalcTruncMeanProfile(res_range_good,dEdx_good, trunc_dEdx);			

                double pida_sum_trunc=0.0;
                //Calculate the mean truncated mean dEdx
                for(size_t k=0; k< trunc_dEdx.size(); k++){
                    double dEdx = trunc_dEdx[k];
                    m_reco_track_mean_trunc_dEdx[i_trk] += dEdx;
                    if(k <= trunc_dEdx.size()/2){ 
                        m_reco_track_mean_trunc_dEdx_start_half[i_trk]+=dEdx;
                    }else{
                        m_reco_track_mean_trunc_dEdx_end_half[i_trk]+=dEdx;
                    }


                    if(trunc_dEdx[k] != trunc_dEdx[k] || std::isinf(trunc_dEdx[k]) || trunc_dEdx[k]<0){
                        std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx[k]<<std::endl;
                        std::cout<<"Vector Length : "<<trunc_dEdx.size()<<std::endl;
                        std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
                        //for(int m=0; m<trunc_dEdx.size(); m++){
                        //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
                        //}
                        std::cout<<"Using Radius: "<<rad<<std::endl;
                        //exit(EXIT_FAILURE);
                        m_reco_track_good_calo[i_trk] = 0; 
                    }

                    pida_sum_trunc += trunc_dEdx[k]/(pow(res_range_good[k],-0.42));
                }
                m_reco_track_trunc_PIDA[i_trk] = pida_sum_trunc;           
            }

            m_reco_track_mean_dEdx[i_trk]            *=1.0/((double)calo_length);
            m_reco_track_mean_dEdx_start_half[i_trk] *=2.0/((double)calo_length);
            m_reco_track_mean_dEdx_end_half[i_trk]   *=2.0/((double)calo_length);
            m_reco_track_mean_trunc_dEdx[i_trk]            *=1.0/((double)trunc_dEdx.size());
            m_reco_track_mean_trunc_dEdx_start_half[i_trk] *=2.0/((double)trunc_dEdx.size());
            m_reco_track_mean_trunc_dEdx_end_half[i_trk]   *=2.0/((double)trunc_dEdx.size());
            m_reco_track_trunc_PIDA[i_trk]  *=1.0/((double)trunc_dEdx.size());
        }
    }



    void SinglePhoton::CollectPID( std::vector<art::Ptr<recob::Track>> & tracks,
            std::map< art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> & trackToPIDMap){

        for(size_t i_trk=0; i_trk<tracks.size(); ++i_trk){
            art::Ptr<recob::Track> track = tracks[i_trk];
            art::Ptr<anab::ParticleID> pid = trackToPIDMap[track];
            if (!pid) {
                std::cout << "[analyze_Tracks] bad PID object" << std::endl;
                continue;
            }

            // For each PID object, create vector of PID scores for each algorithm
            // Loop over this and get scores for algorithm of choice
            // But first, prepare garbage values, just in case
            std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pid->ParticleIDAlgScores();
            double pidScore_BL_plane2 = -999;
            double pidScore_PIDA_plane2 = -999;
            double pidScore_Chi_plane2 = -999;

            int planeid = 2;
            for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++) {
                anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
                //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
                ///////////******// HARD-CODING ALERT PLZ FIX LATER K THX BYE //////////////////*******
                if (planeid != 2){
                    std::cout << "[ParticleIDValidation] Not using information for plane " 
                        << planeid << " (using plane 2 calorimetry only)" << std::endl;
                    continue;
                }
                if (AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 13){
                    pidScore_BL_plane2 = AlgScore.fValue;
                }
                if (AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && TMath::Abs(AlgScore.fAssumedPdg) == 13){
                    pidScore_Chi_plane2 = AlgScore.fValue;
                }
                if (AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA && TMath::Abs(AlgScore.fAssumedPdg) == 13){
                    pidScore_PIDA_plane2 = AlgScore.fValue;
                }



            }
            std::cout << "Setting pid score " << pidScore_BL_plane2 << std::endl;
            m_reco_track_pid_bragg_likelihood_plane2[i_trk] = pidScore_BL_plane2;
            m_reco_track_pid_pida_plane2[i_trk] = pidScore_PIDA_plane2;
            m_reco_track_pid_chi_plane2[i_trk] = pidScore_Chi_plane2;

        }
        return;
    }



}

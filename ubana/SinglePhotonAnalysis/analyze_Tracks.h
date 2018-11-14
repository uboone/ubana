#include "SinglePhoton_module.h"
#include "TPrincipal.h"
#include "TVectorD.h"

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
    
        m_sim_track_energy.clear();
        m_sim_track_pdg.clear();
        m_sim_track_origin.clear();
        m_sim_track_process.clear();
        m_sim_track_startx.clear();
        m_sim_track_starty.clear();
        m_sim_track_startz.clear();


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

        m_sim_track_energy.resize(size);
        m_sim_track_pdg.resize(size);
        m_sim_track_origin.resize(size);
        m_sim_track_process.resize(size);
        m_sim_track_startx.resize(size);
        m_sim_track_starty.resize(size);
        m_sim_track_startz.resize(size);

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


        vertex_tree->Branch("sim_track_energy",&m_sim_track_energy);
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



        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Starting recob::Track analysis\n";

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

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t On Track: "<<i_trk<<" with TrackID: "<<m_trkid<<" and length: "<<m_length<<"\n";

            m_reco_track_num_spacepoints[i_trk] = (int)trk_spacepoints.size();


            m_reco_track_startx[i_trk] = track->Start().X();   
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
            
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Beginining PCA analysis of track\n";
           
            TPrincipal* principal = new TPrincipal(3,"ND");
            for(int x = 0; x < m_reco_track_num_spacepoints[i_trk]; x++){
                std::vector<double> tmp_spacepoints = {trk_spacepoints[x]->XYZ()[0],trk_spacepoints[x]->XYZ()[1] , trk_spacepoints[x]->XYZ()[2]};
                principal->AddRow(&tmp_spacepoints[0]);

                double dist = this->dist_line_point(tmp_trk_start,tmp_trk_end,tmp_spacepoints);
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

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Completed PCA analysis of track. Primary componant: "<<m_reco_track_spacepoint_principal0.back()<<"\n";

            //range based energy calculation assuming
            m_reco_track_proton_kinetic_energy[i_trk] = proton_length2energy_tgraph.Eval(m_length)/1000.0; 
            if(m_length == 0.0) m_reco_track_proton_kinetic_energy[i_trk]=0.0;


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


        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Finished.\n";
    }

    void SinglePhoton::RecoMCTracks(const std::vector<art::Ptr<recob::Track>>& tracks,  
            std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap, 
            std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,
            std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap){


        auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

        if(m_is_verbose) std::cout<<"SinglePhoton::RecoMCTracks()\t||\t Begininning recob::Track Reco-MC suite\n";

        int i_trk = 0;
        for (TrackVector::const_iterator iter = tracks.begin(), iterEnd = tracks.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Track> track = *iter;
            const art::Ptr<simb::MCParticle> mcparticle = trackToMCParticleMap[track];
            const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruthMap[mcparticle];

            double kx = mcparticle->Position().X();
            double ky = mcparticle->Position().Y();
            double kz = mcparticle->Position().Z();


            m_sim_track_energy[i_trk] = mcparticle->E();
            m_sim_track_pdg[i_trk] = mcparticle->PdgCode();
            m_sim_track_process[i_trk] = mcparticle->Process();
            m_sim_track_startx[i_trk] = kx+SCE->GetPosOffsets(geo::Point_t(kx,ky,kz)).X();
            m_sim_track_startx[i_trk] = ky+SCE->GetPosOffsets(geo::Point_t(kx,ky,kz)).Y();
            m_sim_track_startx[i_trk] = kz+SCE->GetPosOffsets(geo::Point_t(kx,ky,kz)).Z();
            m_sim_track_origin[i_trk] = mctruth->Origin();

            i_trk++;
        }

    }





    void SinglePhoton::CollectCalo(const art::Event &evt, const art::Ptr<recob::Track> &track)
    {

        if(m_is_verbose) std::cout<<"SinglePhoton::CollectCalo(recob::Track)\t||\t Starting calo module for recob::Track\n";

        art::ValidHandle< std::vector<recob::Track>> theTracks = evt.getValidHandle<std::vector<recob::Track>>(m_trackLabel);
        //art::Handle< std::vector<anab::Calorimetry> > theCalo;
        //evt.getByLabel(label, theCalo);

        if (!theTracks.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find Track Information... " << "\n";
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theTracks->size() << " CaloInfo " << "\n";
        }

        art::FindManyP<anab::Calorimetry> theCaloAssns(theTracks, evt, m_caloLabel);

        int i = track->ID();
        const std::vector< art::Ptr<anab::Calorimetry> > calo = theCaloAssns.at(i);

        for (unsigned int j=0; j<calo.size(); ++j)
        {
            //std::cout<<"Plane: "<<j<<"\n";
            for (size_t iTrkHit = 0; iTrkHit < calo[j]->dEdx().size(); ++iTrkHit) {
 //               std::cout<<"\t"<<iTrkHit<<" "<<calo[j]->dEdx()[iTrkHit]<<" "<<calo[j]->ResidualRange()[iTrkHit]<<"\n";
            }
        }
    }

    double SinglePhoton::dist_line_point( std::vector<double>&X1, std::vector<double>& X2, std::vector<double>& X0){
        double x1 =X1.at(0);
        double y1 =X1.at(1);
        double z1 =X1.at(2);

        double x2 =X2.at(0);
        double y2 =X2.at(1);
        double z2 =X2.at(2);

        double x0 =X0.at(0);
        double y0 =X0.at(1);
        double z0 =X0.at(2);

        double x10 = x1-x0;
        double y10 = y1-y0;
        double z10 = z1-z0;

        double x21 = x2-x1;
        double y21 = y2-y1;
        double z21 = z2-z1;

        double t = -(x10*x21+y10*y21+z10*z21)/fabs(x21*x21+y21*y21+z21*z21 );

        double d2 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0)+(z2-z1)*(z1-z0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));


        return sqrt(d2);

    }


}

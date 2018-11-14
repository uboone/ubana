#include "SinglePhoton_module.h"


namespace single_photon
{
    void SinglePhoton::ClearShowers(){
        m_reco_asso_showers=0;
        m_reco_shower_startx.clear();
        m_reco_shower_starty.clear();
        m_reco_shower_startz.clear();
        m_reco_shower_dirx.clear();
        m_reco_shower_diry.clear();
        m_reco_shower_dirz.clear();
        m_reco_shower_theta_yz.clear();
        m_reco_shower_phi_yx.clear();
        m_reco_shower_conversion_distance.clear();
        m_reco_shower_openingangle.clear();
        m_reco_shower_length.clear();
        m_reco_shower_delaunay_num_triangles_plane0.clear();
        m_reco_shower_delaunay_num_triangles_plane1.clear();
        m_reco_shower_delaunay_num_triangles_plane2.clear();
        m_reco_shower_num_hits_plane0.clear();
        m_reco_shower_num_hits_plane1.clear();
        m_reco_shower_num_hits_plane2.clear();
        m_reco_shower_delaunay_area_plane0.clear();
        m_reco_shower_delaunay_area_plane1.clear();
        m_reco_shower_delaunay_area_plane2.clear();

        m_sim_shower_energy.clear();
        m_sim_shower_pdg.clear();
        m_sim_shower_origin.clear();
        m_sim_shower_process.clear();
        m_sim_shower_startx.clear();
        m_sim_shower_starty.clear();
        m_sim_shower_startz.clear();

    }

    void SinglePhoton::ResizeShowers(size_t size){
        m_reco_shower_startx.resize(size);
        m_reco_shower_starty.resize(size);
        m_reco_shower_startz.resize(size);
        m_reco_shower_dirx.resize(size);
        m_reco_shower_diry.resize(size);
        m_reco_shower_dirz.resize(size);
        m_reco_shower_theta_yz.resize(size);
        m_reco_shower_phi_yx.resize(size);
        m_reco_shower_conversion_distance.resize(size);
        m_reco_shower_openingangle.resize(size);
        m_reco_shower_length.resize(size);
        m_reco_shower_delaunay_num_triangles_plane0.resize(size);
        m_reco_shower_delaunay_num_triangles_plane1.resize(size);
        m_reco_shower_delaunay_num_triangles_plane2.resize(size);
        m_reco_shower_num_hits_plane0.resize(size);
        m_reco_shower_num_hits_plane1.resize(size);
        m_reco_shower_num_hits_plane2.resize(size);
        m_reco_shower_delaunay_area_plane0.resize(size);
        m_reco_shower_delaunay_area_plane1.resize(size);
        m_reco_shower_delaunay_area_plane2.resize(size);


        m_sim_shower_energy.resize(size);
        m_sim_shower_pdg.resize(size);
        m_sim_shower_origin.resize(size);
        m_sim_shower_process.resize(size);
        m_sim_shower_startx.resize(size);
        m_sim_shower_starty.resize(size);
        m_sim_shower_startz.resize(size);
    }

    void SinglePhoton::CreateShowerBranches(){
        vertex_tree->Branch("reco_asso_showers",&m_reco_asso_showers,"reco_asso_showers/I");
        vertex_tree->Branch("reco_shower_length", &m_reco_shower_length);
        vertex_tree->Branch("reco_shower_opening_angle", &m_reco_shower_openingangle);
        vertex_tree->Branch("reco_shower_dirx", &m_reco_shower_dirx);
        vertex_tree->Branch("reco_shower_diry", &m_reco_shower_diry);
        vertex_tree->Branch("reco_shower_dirz", &m_reco_shower_dirz);
        vertex_tree->Branch("reco_shower_startx", &m_reco_shower_startx);
        vertex_tree->Branch("reco_shower_starty", &m_reco_shower_starty);
        vertex_tree->Branch("reco_shower_startz", &m_reco_shower_startz);
        vertex_tree->Branch("reco_shower_theta_yz",&m_reco_shower_theta_yz);
        vertex_tree->Branch("reco_shower_phi_yx",&m_reco_shower_phi_yx);
        vertex_tree->Branch("reco_shower_conversion_distance",& m_reco_shower_conversion_distance);
        vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane0",&m_reco_shower_delaunay_num_triangles_plane0);
        vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane1",&m_reco_shower_delaunay_num_triangles_plane1);
        vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane2",&m_reco_shower_delaunay_num_triangles_plane2);
        vertex_tree->Branch("reco_shower_num_hits_plane0",&m_reco_shower_num_hits_plane0);
        vertex_tree->Branch("reco_shower_num_hits_plane1",&m_reco_shower_num_hits_plane1);
        vertex_tree->Branch("reco_shower_num_hits_plane2",&m_reco_shower_num_hits_plane2);
        vertex_tree->Branch("reco_shower_delaunay_area_plane0",&m_reco_shower_delaunay_area_plane0);
        vertex_tree->Branch("reco_shower_delaunay_area_plane1",&m_reco_shower_delaunay_area_plane1);
        vertex_tree->Branch("reco_shower_delaunay_area_plane2",&m_reco_shower_delaunay_area_plane2);


        vertex_tree->Branch("sim_shower_energy",&m_sim_shower_energy);
        vertex_tree->Branch("sim_shower_pdg",&m_sim_shower_pdg);
        vertex_tree->Branch("sim_shower_origin",&m_sim_shower_origin);
        vertex_tree->Branch("sim_shower_process",&m_sim_shower_process);
        vertex_tree->Branch("sim_shower_startx",&m_sim_shower_startx);
        vertex_tree->Branch("sim_shower_starty",&m_sim_shower_starty);
        vertex_tree->Branch("sim_shower_startz",&m_sim_shower_startz);
    }

    void SinglePhoton::AnalyzeShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap){

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Begininning recob::Shower analysis suite\n";

        m_reco_asso_showers=showers.size();
        int i_shr = 0;
        this->ResizeShowers(m_reco_asso_showers);

        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Shower> shower = *iter;
            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];
            const std::vector<art::Ptr<recob::Hit>> hits =  pfParticleToHitMap[pfp];

            //int m_shrid = shower->ID(); This is an used variable, always -999
            double m_length = shower->Length();
            double m_open_angle = shower->OpenAngle();

            TVector3 shr_start = shower->ShowerStart();
            TVector3 shr_dir = shower->Direction();

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t On Shower: "<<i_shr<<" which has length: "<<m_length<<"\n";

            m_reco_shower_startx[i_shr] = shr_start.X();
            m_reco_shower_starty[i_shr] = shr_start.Y();
            m_reco_shower_startz[i_shr] = shr_start.Z();

            m_reco_shower_dirx[i_shr] = shr_dir.X();
            m_reco_shower_diry[i_shr] = shr_dir.Y();
            m_reco_shower_dirz[i_shr] = shr_dir.Z();

            m_reco_shower_length[i_shr] = m_length;
            m_reco_shower_openingangle[i_shr] = m_open_angle;

            m_reco_shower_conversion_distance[i_shr] = sqrt( pow(shr_start.X()-m_vertex_pos_x,2)+pow(shr_start.Y()-m_vertex_pos_y,2)+ pow(shr_start.Z()-m_vertex_pos_z,2)  );

            m_reco_shower_theta_yz[i_shr] = atan2(m_reco_shower_diry[i_shr],m_reco_shower_dirz[i_shr]);
            m_reco_shower_phi_yx[i_shr] = atan2(m_reco_shower_diry[i_shr],m_reco_shower_dirx[i_shr]);


            std::vector<int> t_num(3,0);
            std::vector<int> t_numhits(3,0);
            std::vector<double> t_area(3,0.0);
            //Right, this basically loops over all hits in all planes and for each plane forms the Delaunay triangilization of it and calculates the 2D area inscribed by the convex hull

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Starting Delaunay Triangleization\n";

            auto start = std::chrono::high_resolution_clock::now();
            this->delaunay_hit_wrapper(hits, t_numhits, t_num, t_area);

            auto finish = std::chrono::high_resolution_clock::now();
            auto microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Finished Delaunay Triangleization. It took "<< microseconds.count() << "ms and found "<<t_num[0]+t_num[1]+t_num[2]<<" triangles\n";

            m_reco_shower_delaunay_num_triangles_plane0[i_shr] = t_num[0];
            m_reco_shower_delaunay_num_triangles_plane1[i_shr] = t_num[1];
            m_reco_shower_delaunay_num_triangles_plane2[i_shr] = t_num[2];

            m_reco_shower_delaunay_area_plane0[i_shr] = t_area[0];
            m_reco_shower_delaunay_area_plane1[i_shr] = t_area[1];
            m_reco_shower_delaunay_area_plane2[i_shr] = t_area[2];

            m_reco_shower_num_hits_plane0[i_shr] = t_numhits[0];
            m_reco_shower_num_hits_plane1[i_shr] = t_numhits[1];
            m_reco_shower_num_hits_plane2[i_shr] = t_numhits[2];


            //-------------- Flashes : Was there a flash in the beam_time and if so was it near in Z? --------------------
            double zmin = m_reco_shower_startz[i_shr];
            double zmax = zmin + m_reco_shower_dirz[i_shr]*m_reco_shower_length[i_shr];
            if(zmin > zmax) std::swap(zmin, zmax);

            double ymin = m_reco_shower_starty[i_shr];
            double ymax = zmin + m_reco_shower_diry[i_shr]*m_reco_shower_length[i_shr];
            if(ymin > ymax) std::swap(ymin, ymax);

            //Now loop over all flashes (only in beamtime) and find SOMETHING?

            //---------------- Reco-Truth Matching to MCParticles -------------------//



            i_shr++;
        }


        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Finished.\n";
    }


    //-----------------------------------------------------------------------------------------------------------------------------------------
    void SinglePhoton::RecoMCShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  
            std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap, 
            std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
            std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap){

        if(m_is_verbose) std::cout<<"SinglePhoton::RecoMCShowers()\t||\t Begininning recob::Shower Reco-MC suite\n";

        int i_shr = 0;
        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Shower> shower = *iter;
            const art::Ptr<simb::MCParticle> mcparticle = showerToMCParticleMap[shower];
            const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruthMap[mcparticle];

            std::vector<double> corrected(3);
            this->spacecharge_correction(mcparticle, corrected);

            m_sim_shower_energy[i_shr] = mcparticle->E();
            m_sim_shower_pdg[i_shr] = mcparticle->PdgCode();
            m_sim_shower_process[i_shr] = mcparticle->Process();
            m_sim_shower_startx[i_shr] = mcparticle->Position().X()+corrected[0];
            m_sim_shower_starty[i_shr] = mcparticle->Position().Y()+corrected[1];
            m_sim_shower_startz[i_shr] =mcparticle->Position().Z()+corrected[2];
            m_sim_shower_origin[i_shr] = mctruth->Origin();

            i_shr++;
        }

    }


    void SinglePhoton::CollectCalo(const art::Event &evt, const art::Ptr<recob::Shower> &shower)
    {

    }


}

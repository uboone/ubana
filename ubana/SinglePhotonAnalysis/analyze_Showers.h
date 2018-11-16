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
	//track shower matching
	m_reco_shower_flash_shortest_distz.clear();
	m_reco_shower_flash_shortest_disty.clear();
	m_reco_shower_flash_shortest_distyz.clear();
 
	m_reco_shower_flash_shortest_index_z.clear();
	m_reco_shower_flash_shortest_index_y.clear();
	m_reco_shower_flash_shortest_index_yz.clear();

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
	m_reco_shower_flash_shortest_distz.resize(size);
	m_reco_shower_flash_shortest_disty.resize(size);
	m_reco_shower_flash_shortest_distyz.resize(size);
 
	m_reco_shower_flash_shortest_index_z.resize(size);
	m_reco_shower_flash_shortest_index_y.resize(size);
	m_reco_shower_flash_shortest_index_yz.resize(size);

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
	vertex_tree->Branch("reco_shower_flash_shortest_distz",&m_reco_shower_flash_shortest_distz);
	vertex_tree->Branch("reco_shower_flash_shortest_disty",&m_reco_shower_flash_shortest_disty);
	vertex_tree->Branch("reco_shower_flash_shortest_distyz",&m_reco_shower_flash_shortest_disty);
	vertex_tree->Branch("reco_shower_flash_shortest_index_z", &m_reco_shower_flash_shortest_index_z);
	vertex_tree->Branch("reco_shower_flash_shortest_index_y", &m_reco_shower_flash_shortest_index_y);
	vertex_tree->Branch("reco_shower_flash_shortest_index_yz", &m_reco_shower_flash_shortest_index_yz);
	
	

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
	    if(m_is_verbose) std::cout << "shower start z_plane: " <<zmin << "\n";
            double zmax = zmin + m_reco_shower_dirz[i_shr]*m_reco_shower_length[i_shr];
	    if(m_is_verbose) std::cout << "shower end z_plane: " <<zmax << "\n";
            if(zmin > zmax) std::swap(zmin, zmax);

            double ymin = m_reco_shower_starty[i_shr];
	    if(m_is_verbose) std::cout << "shower start y_plane: " <<ymin << "\n";
            double ymax = zmin + m_reco_shower_diry[i_shr]*m_reco_shower_length[i_shr];
	    if(m_is_verbose) std::cout << "shower end y_plane: " <<ymax << "\n";
            if(ymin > ymax) std::swap(ymin, ymax);

            //Now loop over all flashes (only in beamtime) and find SOMETHING?
	    //Code property of Gray Yarbrough (all rights reserved)
	    //int optical_flash_in_beamgate_counter=0;
	    double shortest_dist_to_flash_z=DBL_MAX;
	    double shortest_dist_to_flash_y=DBL_MAX;
	    double shortest_dist_to_flash_yz=DBL_MAX;
	    //-999 my nonsenese int can change
	    int shortest_dist_to_flash_index_z=-999;
	    int shortest_dist_to_flash_index_y=-999;
	    int shortest_dist_to_flash_index_yz=-999;
	    if(m_is_verbose) std::cout<<" number of flashes: "<< m_reco_num_flashes<< "\n";
	    for(int i_flash = 0; i_flash < m_reco_num_flashes; ++i_flash) {
	      //double const flash_time=m_reco_flash_time[i_flash];
	      double const zcenter=m_reco_flash_zcenter[i_flash];
	      if(m_is_verbose) std::cout<< "flash z center:" <<m_reco_flash_zcenter[i_flash]<< "\n";
	      double const ycenter=m_reco_flash_ycenter[i_flash];
	      if(m_is_verbose) std::cout<< "flash y center:" <<m_reco_flash_ycenter[i_flash]<< "\n";
	      //if((m_is_verbose) std::cout<<"optical flash time"<< m_reco_flash_time[i_flash] << " (opf.Time()>3.2 && opf.Time<4.8)\n

	      //z plane
	      double dist_z=DBL_MAX;
	       if(zcenter < zmin) {
		 dist_z = zmin - zcenter;
		 if(m_is_verbose) std::cout << "z plane flash center - zmin dist: " << dist_z << "\n";
	       }
	       else if(zcenter > zmax) {
		 dist_z = zcenter - zmax;
		 if(m_is_verbose) std::cout << " z plane flash center - zmax dist: " << dist_z << "\n";
	       }
	       else {
		 dist_z = 0;
		 if(m_is_verbose) std::cout << "z plane flash center inside shower\n";
	       }	    
	       if(dist_z < shortest_dist_to_flash_z){
		 shortest_dist_to_flash_z = dist_z;
		 shortest_dist_to_flash_index_z=i_flash;
	       }
	      
	   
	    //y plane

	       double dist_y=DBL_MAX;
	       if(ycenter < ymin) {
		 dist_y = ymin - ycenter;
		 if(m_is_verbose) std::cout << "y plane flash center - zmin dist: " << dist_y << "\n";
	       }
	       else if(ycenter > ymax) {
		 dist_y = ycenter - ymax;
		 if(m_is_verbose) std::cout << " y plane flash center - zmax dist: " << dist_y << "\n";
	       }
	       else {
		 dist_y= 0;
		 if(m_is_verbose) std::cout << "y plane flash center inside shower\n";
	       }	    
	       if(dist_y < shortest_dist_to_flash_y){
		 shortest_dist_to_flash_y = dist_y;
		 shortest_dist_to_flash_index_y=i_flash;
	       }
	       double dist_yz=DBL_MAX;
	       dist_yz=std::sqrt(dist_y*dist_y+dist_z*dist_z);
	       if(dist_yz<shortest_dist_to_flash_yz){
		 shortest_dist_to_flash_yz = dist_yz;
		 shortest_dist_to_flash_index_yz=i_flash;
	       }
	       
	    }

	    
	    //assume setting to nonsense value
	    if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_z = -2;
	    m_reco_shower_flash_shortest_distz[i_shr]=shortest_dist_to_flash_z;
	    if(m_is_verbose) std::cout << "the shortest distance in z plane between a flash and the shower: " << shortest_dist_to_flash_z << "\n";
	    m_reco_shower_flash_shortest_index_z[i_shr]=shortest_dist_to_flash_index_z;
	    if(m_is_verbose) std::cout << "the index closest flash to shower z_plane: " << shortest_dist_to_flash_index_z << "\n";
	     
	    if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_y = -2;
	    m_reco_shower_flash_shortest_disty[i_shr]=shortest_dist_to_flash_y;
	    if(m_is_verbose) std::cout <<"the shortest distance in y plane between a flash and the shower: " << shortest_dist_to_flash_y << "\n";
	    m_reco_shower_flash_shortest_index_y[i_shr]=shortest_dist_to_flash_index_y;
	    if(m_is_verbose) std::cout << "the index closest flash to shower y_plane: " << shortest_dist_to_flash_index_y << "\n";

	    if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_yz = -2;
	    m_reco_shower_flash_shortest_distyz[i_shr]=shortest_dist_to_flash_yz;
	    if(m_is_verbose) std::cout <<"the shortest distance in yz between a flash and the shower: " << shortest_dist_to_flash_yz << "\n";
	    m_reco_shower_flash_shortest_index_yz[i_shr]=shortest_dist_to_flash_index_yz;
	    if(m_is_verbose) std::cout << "the index closest flash to shower yz: " << shortest_dist_to_flash_index_yz << "\n";

	
	    //end optical flash code
            i_shr++;
	}


           if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Finished.\n";
}

    //-----------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::CollectCalo(const art::Event &evt, const art::Ptr<recob::Shower> &shower)
    {

        art::ValidHandle< std::vector<recob::Shower>> theShowers = evt.getValidHandle<std::vector<recob::Shower>>(m_showerLabel);
        if (!theShowers.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find Shower Information... " << "\n";
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theShowers->size() << " Showers " << "\n";
        }
        art::FindManyP<anab::Calorimetry> theCaloAssns(theShowers, evt, m_caloLabel);



        for (unsigned int i = 0; i < theShowers->size(); ++i)
        {
            const art::Ptr<recob::Shower> tmp_shower(theShowers, i);
            const std::vector< art::Ptr<anab::Calorimetry> > calo = theCaloAssns.at(i);

            if(tmp_shower == shower){

                break;
            }
        }

    }


}

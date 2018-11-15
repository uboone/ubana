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

	m_reco_shower_energy.clear();
        m_reco_shower_dQdx_plane0.clear();
        m_reco_shower_dQdx_plane2.clear();
        m_reco_shower_dQdx_plane2.clear();
        m_reco_shower_dEdx_plane0.clear();
        m_reco_shower_dEdx_plane2.clear();
        m_reco_shower_dEdx_plane2.clear();

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

	m_reco_shower_energy.resize(size);
	m_reco_shower_dQdx_plane0.resize(size);
	m_reco_shower_dQdx_plane1.resize(size);
	m_reco_shower_dQdx_plane2.resize(size);
	m_reco_shower_dEdx_plane0.resize(size);
	m_reco_shower_dEdx_plane1.resize(size);
	m_reco_shower_dEdx_plane2.resize(size);


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

	//the calorimetry info
	vertex_tree->Branch("reco_shower_energy",&m_reco_shower_energy);
	vertex_tree->Branch("reco_shower_dQdx_plane0",&m_reco_shower_dQdx_plane0);
	vertex_tree->Branch("reco_shower_dQdx_plane1",&m_reco_shower_dQdx_plane1);
	vertex_tree->Branch("reco_shower_dQdx_plane2",&m_reco_shower_dQdx_plane2);
	vertex_tree->Branch("reco_shower_dEdx_plane0",&m_reco_shower_dEdx_plane0);
	vertex_tree->Branch("reco_shower_dEdx_plane1",&m_reco_shower_dEdx_plane1);
	vertex_tree->Branch("reco_shower_dEdx_plane2",&m_reco_shower_dEdx_plane2);
	


    }

    void SinglePhoton::AnalyzeShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap, std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > & pfParticleToClusterMap,std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> >  & clusterToHitMap ){

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Begininning recob::Shower analysis suite\n";
            
        m_reco_asso_showers=showers.size();
        int i_shr = 0;
        this->ResizeShowers(m_reco_asso_showers);

        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Shower> shower = *iter;
            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];
            const std::vector<art::Ptr<recob::Hit>> hits =  pfParticleToHitMap[pfp];
	    const std::vector<art::Ptr<recob::Cluster>> clusters = pfParticleToClusterMap[pfp];
            
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

	  //-------------- Calorimetry --------------------

    
	    m_reco_shower_energy[i_shr] = CalcEShower(hits);
	    if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 0"<<std::endl;
	    m_reco_shower_dQdx_plane0[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 0 ); 
	    if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 1"<<std::endl;
	    m_reco_shower_dQdx_plane1[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 1 ); 
	    if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 2"<<std::endl; 
	    m_reco_shower_dQdx_plane2[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 2 ); 


            //-------------- Flashes : Was there a flash in the beam_time and if so was it near in Z? --------------------
            double zmin = m_reco_shower_startz[i_shr];
            double zmax = zmin + m_reco_shower_dirz[i_shr]*m_reco_shower_length[i_shr];
            if(zmin > zmax) std::swap(zmin, zmax);

            double ymin = m_reco_shower_starty[i_shr];
            double ymax = zmin + m_reco_shower_diry[i_shr]*m_reco_shower_length[i_shr];
            if(ymin > ymax) std::swap(ymin, ymax);

            //Now loop over all flashes (only in beamtime) and find SOMETHING?
            

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
       // art::FindManyP<anab::Calorimetry> theCaloAssns(theShowers, evt, m_caloLabel);



        for (unsigned int i = 0; i < theShowers->size(); ++i)
        {
         //   const art::Ptr<recob::Shower> tmp_shower(theShowers, i);
           // const std::vector< art::Ptr<anab::Calorimetry> > calo = theCaloAssns.at(i);

            //if(tmp_shower == shower){

              //  break;
           // }
        }

    }


   double SinglePhoton::CalcEShower(std::vector<art::Ptr<recob::Hit>> hits){    
     double energy[3] = {0., 0., 0.};
   
     //for each hit in the shower
     for (art::Ptr<recob::Hit> thishitptr : hits){
	//check the plane
	int plane= thishitptr->View();

	//skip invalid planes     	
	if (plane > 2 || plane < 0)	continue;

	//calc the energy of the hit
	double E = SinglePhoton::QtoEConversion(thishitptr);	

	//add the energy to the plane
	energy[plane] += E;
    }//for each hit
    
    //find the max energy on a single plane
    double max = energy[0];
    for (double en: energy){
	if( en > max){
		max = en;
		}
	}
    if (m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The energy on each plane for this shower is "<<energy[0]<<", "<<energy[1]<<", "<<energy[2]<<std::endl;

    //return the highest energy on any of the planes
    return max;

   }

   double SinglePhoton::QtoEConversion(art::Ptr<recob::Hit> thishitptr){
	double Q = thishitptr->Integral();
	//return the energy value converted to MeV (the factor of 1e-6)
	return Q* m_gain* m_work_function *1e-6 /m_recombination_factor;

   }

  std::vector<double> SinglePhoton::CalcdQdxShower(const art::Ptr<recob::Shower> shower, const std::vector<art::Ptr<recob::Cluster>> clusters, std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane){
	//if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The number of clusters in this shower is "<<clusters.size()<<std::endl;
	std::vector<double> dqdx;

	//get the 3D shower direction
	//note: in previous versions of the code there was a check to fix cases where the shower direction was inverted - this hasn't been implemented
	TVector3 shower_dir(shower->Direction().X(), shower->Direction().Y(),shower->Direction().Z());
	 
	//calculate the pitch for this plane
	double pitch = getPitch(shower_dir, plane);	
	//if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The pitch between the shower and plane "<<plane<<" is "<<pitch<<std::endl;

	//for all the clusters in the shower
	for (art::Ptr<recob::Cluster> thiscluster: clusters){
		//keep only clusters on the plane
		if(thiscluster->View() != plane) continue;

		//calculate the cluster direction
		std::vector<double> cluster_axis = {cos(thiscluster->StartAngle()), sin(thiscluster->StartAngle())};		

		//get the cluster start and and in CM
		std::vector<double> cluster_start = {thiscluster->StartWire() * m_wire_spacing,theDetector->ConvertTicksToX( thiscluster->StartTick(), plane,m_TPC ,m_Cryostat)};
      		//if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t the start position of the cluster is "<<cluster_start[0]<<", "<<cluster_start[1]<<std::endl;
		std::vector<double> cluster_end = {thiscluster->EndWire() * m_wire_spacing,theDetector->ConvertTicksToX( thiscluster->EndTick(), plane,m_TPC ,m_Cryostat) };

		//check that the cluster has non-zero length
		double length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) + pow(cluster_end[1] - cluster_start[1], 2));
		//if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The cluster length is "<<length<<std::endl;
		if (length <= 0) continue;

		//draw a rectangle around the cluster axis 
	  	std::vector<std::vector<double>> rectangle = buildRectangle(cluster_start, cluster_axis, m_width_dqdx_box, m_length_dqdx_box);	

		//get all the hits for this cluster
		std::vector<art::Ptr<recob::Hit>> hits =  clusterToHitMap[thiscluster];
	
		//for each hit in the cluster
		std::cout<<"the number of hits in this cluster is "<<hits.size()<<std::endl;
		int n = 0;
		for (art::Ptr<recob::Hit> thishit: hits){	
			std::cout<<"starting hit "<<++n<<std::endl;
			//get the hit position in cm
			std::vector<double> thishit_pos = {thishit->WireID().Wire * m_wire_spacing, theDetector->ConvertTicksToX( thishit->PeakTime(), plane,m_TPC ,m_Cryostat)};
			std::cout<<"the position of this hit in the cluster is "<<thishit_pos[0]<<", "<<thishit_pos[1]<<std::endl;

			//check if inside the box
			if (insideBox(thishit_pos, rectangle)){
				std::cout<<"the position of this hit inside of the rectangle is "<<thishit_pos[0]<<", "<<thishit_pos[1]<<std::endl;
				double q = thishit->Integral() * m_gain;
				std::cout<<"the q for this hit is "<<q<<std::endl;
				double this_dqdx = q/pitch; 
				dqdx.push_back(this_dqdx);
				std::cout<<"the calculated dq/dx for this hit is "<<this_dqdx<<std::endl;
			}//if hit falls inside the box
		//	else{
		//		std::cout<<"hit not inside box"<<std::endl;
		//	}
		}//for each hit inthe cluster
	}//for each cluster
	//std::cout<<"the number of dq/dx points on this plane is "<<dqdx.size()<<std::endl;
	//for (double mydqdx: dqdx){
	//	std::cout<<"final dqdx for plane "<<plane<<" = "<<mydqdx<<std::endl;
	//}
	return dqdx;
  }

  double SinglePhoton::getPitch(TVector3 shower_dir, int plane){
	//geo::TPCGeo const& tpc = geom->TPC(0);
	
	//get the wire direction for this plane
	TVector3 wire_dir =  geom->TPC(m_TPC).Plane(plane).FirstWire().Direction();
	//if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The wire_pitch for plane "<<plane<<" is "<<wire_direction[0]<<", "<<wire_direction[1]<<", "<<wire_direction[2]<<std::endl;	
	
	//take the dot product between the wire direction and the shower direction
	double cos = wire_dir.Dot(shower_dir);

	//want only positive values so take abs, normalize by the lengths of the shower and wire
	cos = abs(cos)/(wire_dir.Mag() * shower_dir.Mag());	

	//If the cos is 0 shower is perpendicular and therefore get infinite distance 
	if (cos == 0){ return std::numeric_limits<double>::max(); }

	//output is always >= the wire spacing
	return m_wire_spacing/cos;
  }


  std::vector<std::vector<double>> SinglePhoton::buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length){
	std::vector<std::vector<double>> corners;

	//get the axis perpedicular to the cluster axis
	double perp_axis[2] = {-cluster_axis[1], cluster_axis[0]};
	
	//create a vector for each corner of the rectangle on the plane
	//c1 = bottom left corner
	std::vector<double> c1 = {cluster_start[0] + perp_axis[0] * width / 2,  cluster_start[1] + perp_axis[1] * width / 2};
	//c2 = top left corner
	std::vector<double> c2 = {c1[0] + cluster_axis[0] * length, c1[1] + cluster_axis[1] * length};
	//c3 = bottom right corner
	std::vector<double> c3 = {cluster_start[0] - perp_axis[0] * width / 2, cluster_start[1] - perp_axis[1] * width / 2};
	//c4 = top right corner
	std::vector<double> c4 ={c3[0] + cluster_axis[0] * length, c3[1] + cluster_axis[1] * length}; 
	
	//save each of the vectors
	corners.push_back(c1);
	corners.push_back(c2);
	corners.push_back(c3);
	corners.push_back(c4);

	return corners;
  }

  bool SinglePhoton::insideBox(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle){
	std::cout<<"checking if inside"<<std::endl;
	//for a rectangle this is a known value but this is the most general
	int n_vertices = (int)rectangle.size();

	bool inside = false;
	int i, j = 0;
	//for each pair of vertices
 	for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
		//std::cout<<"checking if inside at i, j = "<<i<<", "<<j<<std::endl;
		//if the hit y coordinate is between the y and x coordinates of two vertices
   		if ( ((rectangle[i][1]> thishit_pos[1]) != (rectangle[j][1]>thishit_pos[1])) 
				&&(thishit_pos[0] < (rectangle[j][0]-rectangle[i][0]) * (thishit_pos[1]-rectangle[i][1]) / (rectangle[j][1]-rectangle[i][1]) + rectangle[i][0]) ){     			inside = true;
		}
  	}
	//std::cout<<"done checking, inside = "<<inside<<std::endl;
 	return inside;
  }

}

#ifndef VERTEXBUILDER_H
#define VERTEXBUILDER_H

#include "ParticleAssociations.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TTree.h"
#include <fstream>//read and write txt file.

/**********************
 *stuct VertexBuilderTree: it store values of variables extracted 
	from the artObjects? (values from the raw root file)
 *
 *class VertexBuilder: it contains functions for determining the vertex of objects
	assigned by the class ParticleAssociations_all.
 *
 *********************/

using namespace std;

namespace single_photon{
struct VertexBuilderTree {

	TTree * ftree;
	int frun_number;
	int fsubrun_number;
	int fevent_number;
	int ftrack_number;
	int fshower_number;
//	int fassociation_track_number; //CHECK,number of track associations;
//	int fassociation_shower_number; //CHECK,number of shower associations;
	int fassociation_final_number; //number of total associations;

	VertexBuilderTree() :
		ftree(nullptr){}
/*********************************************************************
 *	Setup() - Obtain the addresses of event variables 
 *		(indices for particular events, showers, and tracks).
 *********************************************************************/
	void Setup() {
		art::ServiceHandle< art::TFileService > tfs;
		ftree = tfs->make<TTree>("VertexBuilder", "");
		ftree->Branch("run_number", &frun_number, "run_number/I");
		ftree->Branch("subrun_number", &fsubrun_number, "subrun_number/I");
		ftree->Branch("event_number", &fevent_number, "event_number/I");  
		ftree->Branch("track_number", &ftrack_number, "track_number/I");
		ftree->Branch("shower_number", &fshower_number, "shower_number/I");
//		ftree->Branch("association_track_number", &fassociation_track_number, "association_track_number/I");
//		ftree->Branch("association_shower_number", &fassociation_shower_number, "association_shower_number/I");
		ftree->Branch("association_final_number", &fassociation_final_number, "association_final_number/I");
	}

};


class VertexBuilder {

	geoalgo::GeoAlgo const falgo;

	size_t fobject_id;

	double fstart_prox;//the maximum track proximity threshold (in cm)
	double fshower_prox;//the maximum shower proximity threshold (in cm)
	double fcpoa_vert_prox;//the maximum distance btw shower start to a track end point(in cm); 
	double fcpoa_trackend_prox;//the maximum distance btw mid point of impact parameter of shower to the vertex (in cm);

	DetectorObjects_all const * fdetos;

	VertexBuilderTree * fvbt;

	bool fverbose;

	void CheckSetVariables();// exit if the criteria are not set.

	void Erase(std::multimap<size_t, geoalgo::Point_t const *> & pn,
			std::multimap<size_t, geoalgo::Point_t const *>::iterator const best_it,
			geoalgo::Point_t const & sv);

	/*****************************
	 *
	 * AssociateTracks() - this associates tracks to either end points of other tracks
	 *	whose are nearby.
	 *
	 * ***************************/
	void AssociateTracks(ParticleAssociations_all & pas);

	/*
	 * FindClosestApproach() - find the center points of two points along 
	 *		the backward projection that give minimum impact parameter; this
	 *		also returns the distance of two of such points.
	 *shr1, shr2 - two showers with start points and direction
	 *vtx - a vertex to be updated.
	 */
	double FindClosestApproach(const geoalgo::HalfLine_t & shr1,
			const geoalgo::HalfLine_t & shr2,
			geoalgo::Point_t & vtx) const;

	/*****************************
	 *
	 * AssociateShowers() - Run after AssociateTracks(); 
	 *	this associates showers to either end points of a track
	 *	that is considered to be from the same vertex.
	 *
	 * ***************************/
	void AssociateShowers(ParticleAssociations_all & pas);


	void AddLoneTracks(ParticleAssociations_all & pas);
	void AddLoneShowers(ParticleAssociations_all & pas);


	/*****************************
	 *
	 * AssociateShowers() - Run after AssociateTracks(); 
	 *	this associates showers to either end points of a track
	 *	that is considered to be from the same vertex.
	 *
	 * ***************************/
	void FillVBT(ParticleAssociations_all & pas);

	public:
	vector <double> f_dist_tt;//track&track, start_prox
	vector <double> f_dist_sx;//shower&anything, shower_prox
	vector <double> f_dist_st;//shower&track, cpoa_vert_prox
	vector <double> f_dist_sst;//shower&shower&track, cpoa_trackend_prox


	VertexBuilder();

	void SetVerbose(bool const verbose = true) {//allow output info. to the terminal
		fverbose = verbose;
	}

	void SetParameters (std::vector< double > p){

		fstart_prox			= p[0];
		fshower_prox		= p[1];
		fcpoa_vert_prox		= p[2];
		fcpoa_trackend_prox = p[3];

	 f_dist_tt.push_back(999);//track&track, start_prox
	 f_dist_sx.push_back(999);//shower&anything, shower_prox
	 f_dist_st.push_back(999);//shower&track, cpoa_vert_prox
	 f_dist_sst.push_back(999);//shower&shower&track, cpoa_trackend_prox


		if(fverbose){
			cout<<right<<setw(82)<<" Max. track proximity threshold (t_max)= "<<fstart_prox<<endl;
			cout<<right<<setw(82)<<" Max. shower proximity threshold (s_max)= "<<fshower_prox<<endl;
			cout<<right<<setw(82)<<" Max. distance btw shower start & cloest approach (dp_max)= "<<fcpoa_vert_prox<<endl;
			cout<<right<<setw(82)<<" Max. distance btw midway point of impact parameter to a precandidate vertex (a_max)= "<<fcpoa_trackend_prox<<endl;
		}
	}



//	void SetVBT(VertexBuilderTree * vbt) {fvbt = vbt;}

	//  void AddTracks(art::ValidHandle<std::vector<recob::Track>> const & ev_t);
	//  void AddShowers(art::ValidHandle<std::vector<recob::Shower>> const & ev_s);
//	void AddTracks(std::vector<art::Ptr<recob::Track>> const & ev_t);
//	void AddShowers(std::vector<art::Ptr<recob::Shower>> const & ev_t);
	
	// The main code is here; it calls up all the necessary process.
	void Run(ParticleAssociations_all & pas);

};















//HEADER FILE ARE ABOVE














VertexBuilder::VertexBuilder() :
  fobject_id(0),
  fstart_prox(-1),
  fshower_prox(-1),
  fcpoa_vert_prox(-1),
  fcpoa_trackend_prox(-1),
  fdetos(nullptr),
  fvbt(nullptr),
  fverbose(true) {}


void VertexBuilder::CheckSetVariables() { //this was initialized 

  if(fstart_prox == -1) {//It start_prox was not previously defined, exit.
    std::cout << "fstart_prox, the maximum track proximity threshold, is not set\n";
    exit(1);
  }

  if(fshower_prox == -1) {
    std::cout << "fshower_prox, the maximum shower proximity threshold, is  not set\n";
    exit(1);
  }

  if(fcpoa_vert_prox == -1) {
    std::cout << "fcpoa_vert_prox, the maximum distance btw shower start to the cloest approach, is not set\n";
    exit(1);
  }

  if(fcpoa_trackend_prox == -1) {
    std::cout << "fcpoa_trackend_prox, the maximum distance btw mid point of impact parameter of shower to the vertex, is not set\n";
    exit(1);
  }

}


void VertexBuilder::Erase(std::multimap<size_t, geoalgo::Point_t const *> & pn,
			  std::multimap<size_t, geoalgo::Point_t const *>::iterator const best_it,
			  geoalgo::Point_t const & sv) {//Remove candidates that have been analyzed.

  size_t const index = best_it->first;
  pn.erase(best_it);

  if(fdetos->GetRecoType(index) == fdetos->ftrack_reco_type) {
    auto const pn_it = pn.find(index);
    if(fdetos->GetTrack(index).ftrajectory.Length() < fstart_prox ||
       (pn_it != pn.end() && pn_it->second->Dist(sv) < fstart_prox)) {
      pn.erase(pn_it);
    }
  }

}


void VertexBuilder::AssociateTracks(ParticleAssociations_all & pas) {//Group tracks

	//Ingredients: Key for tracks, Point coordinates along the track (recob:: Trajectory)
	int screen_width = 86;
	bool first_time = true;

	std::multimap<size_t, geoalgo::Point_t const *> pn;//This maps an # to a pointer to the track start/end points
	/* Structure of pn:
	 *	Iterator Key	Value	
	 *	0		Key[0]	point1
	 *	1		Key[0]	point2
	 *	2		Key[1]	pointa
	 *	3		Key[1]	pointb
	 *	.
	 *	example of calling first row:
	 *	pn.begin(), pn[0].first, pn[0].second
	 */

	//#1 load up tracks to pn: 1 id for 2 end points;
	for(size_t const i : fdetos->GetTrackIndices()) {

		Track const & t = fdetos->GetTrack(i);
		geoalgo::Vector ref_point(1,1,1);//use this to pick valid end points; i.e. kill (-999,-999,-999)

		if(t.ftrajectory.front().Dot(ref_point) > -2996){
			pn.emplace(t.fid, &t.ftrajectory.front());//track ID and point;
		}
		if(t.ftrajectory.back().Dot(ref_point) > -2996){
			pn.emplace(t.fid, &t.ftrajectory.back());
		}
		/* appended map of one of this loop:
		 * < track_id, {track beginning point, track ending point}>
		 */
	}

	if(fverbose) {
		std::cout << "STEP I: Load " << pn.size() << " track end points for evaluation.\n\n";

		if(pn.size()) {//Print the title
			cout<<setw(9)<<right<<"Track ID ";
			cout<<setw(15)<<left<<"| Track Length ";
			cout<<setw(37)<<left<<"| Track End Point Coordinates (x,y,z) [cm] from Pandora"<<endl;
			for(int i = 0; i<screen_width; i++) 
				cout<<"-";
			cout<<endl;
			for(auto p : pn){
				cout<<setw(9)<<left<< p.first<<"| ";//Index
				cout<<setw(13)<<left<<fdetos->GetTrack(p.first).ftrajectory.Length();//Track Length
				//ftrajectory is the obejct in art library
				cout<<"| "<<*p.second<<endl;//Coordinates
			}
			for(int i = 0; i<screen_width; i++) 
				cout<<"-";
			cout<<"\n\nSTEP II: Look for vertex candidate among tracks:"<<endl;
		}
	}

	//#2 evaluate all loaded end points of tracks; 
	while(pn.size() > 1) {//while there are tracks to be evaluated, loop them and try to make candidate vertex.

		auto best_match = pn.end(); //best match (iterator) is assumed to be the iterator of the last track;
		auto best_candidate = best_match; //best candidate;

		double best_dist = fstart_prox;//Max. track proximity threshold (t_max)
		if(fverbose)
			std::cout <<  "\tRemain "<< pn.size() <<" track end points for evaluation," << std::endl;

		for(auto match_this = pn.begin(); match_this != pn.end(); ++match_this) {//iterator to indicate the point we are evaluating

			if(std::next(match_this) == pn.end()) break;//stop when there is only one track end point left;
			if(fverbose)
				std::cout << "\t\tEnd point from the track (ID: " << match_this->first <<"):"<< std::endl;


			for(auto compare_this = std::next(match_this); compare_this != pn.end(); ++compare_this) {
				if(fverbose)
					std::cout << "\t\t\tCompare to a end point from a track (ID: " << compare_this->first<<"); ";
				if(match_this->first == compare_this->first){//if two end points from the same track, skip this round.
					if(fverbose) std::cout<<"\t\t\t\tSame track, skip!"<<endl;
					continue;
				}

				double const dist = compare_this->second->Dist(*match_this->second);//distance between match_this point to compare_this point;
				if(first_time){
					f_dist_tt.erase(f_dist_tt.begin());
					first_time = false;
				}
				f_dist_tt.push_back(dist);

				if(fverbose)
					std::cout << "distance btw them: " << dist << ", the shortest: "
						<< best_dist << std::endl;
				//				if(dist==0) std::cout<<" Reject this point."<<endl;
				if(dist < best_dist) {
					if(fverbose)
						std::cout << "\t\t\t\t>>Take this track end point as a better candidate vertex!\n\n";	
					best_match = match_this;
					best_candidate = compare_this;
					best_dist = dist;
				}
			}
		}

		if(best_match == pn.end() || pn.size()<2) {
			if(fverbose)
				std::cout << "\tNah, tracks are all single. No vertex candidate for tracks.\n";
			return;
		}
		//at this moment, we have:
		//	best_match - the iterator of an end-point of one track
		//	best_candidate - the iterator of an end-point of another track
		//	best_dist - the cloest distance between these two points

		std::vector<size_t> track_IDs;
		track_IDs.push_back(best_match->first);
		track_IDs.push_back(best_candidate->first);

		std::vector<geoalgo::Point_t> points;
		points.push_back(*best_match->second);
		points.push_back(*best_candidate->second);

		geoalgo::Sphere sphere(falgo.boundingSphere(points));//falgo is an GeoAlgo object, this sets the mid point of two track endpoints as new candidate verex.

		Erase(pn, best_match, sphere.Center());//remove the points that has been added as candidate vertex.
		Erase(pn, best_candidate, sphere.Center());

		//#3 add more tracks around the vertex candidate
		if(fverbose){
			std::cout << "\n\nSTEP III: Add tracks around the best track (ID: "<< best_match->first<<") end point." <<endl;
		}

		do {//have to break to leave.
			auto best_o = pn.end();//best object to add;
			double sbest_dist = fstart_prox;//Max. track proximity threshold (t_max)

			geoalgo::Point_t const & sphere_center = sphere.Center();//this sets the candidate vertex;

			if(fverbose)
				std::cout << "\tFind the distance between vertices and sphere centre; number of track end points to be considered: " << pn.size() << std::endl;

			for(auto this_iterator = pn.begin(); this_iterator != pn.end(); ++this_iterator) {

				double const dist = this_iterator->second->Dist(sphere_center);
				if(fverbose){
					std::cout << "\t\tTrack (ID: " << this_iterator->first<< ") end point ";
					std::cout << "is at distant: " << dist << " sbest_dist: " << sbest_dist << std::endl;
				}

				if(std::find(track_IDs.begin(), track_IDs.end(), this_iterator->first) != track_IDs.end()) {
					if(fverbose)
						std::cout << " This is not the last track?(CHECK) Skip\n";
					continue;
				}

				if(dist < sbest_dist) {//when find a close track end point, attch it to the candidate vertex.
					if(fverbose)
						std::cout << "\t\t\tTake this track. Update the best distance!\n";
					best_o = this_iterator;
					sbest_dist = dist;
				}
			}

			if(best_o == pn.end()) {//exit the loop only when 
				if(fverbose)
					std::cout << "\tNo more tracks can be added.\n";
				break;
			}

			track_IDs.push_back(best_o->first);
			points.push_back(*best_o->second);

			Erase(pn, best_o, *best_o->second);

			//s = algo.boundingSphere(points);

		} while(true);//if false, then it only loops once;
		//points - point objects that contains positions;
		//boundingSphere(points) is a minimum sphere that contains points.
		pas.AddAssociation(track_IDs, points, sphere.Center(), falgo.boundingSphere(points).Radius());

		if(fverbose){//no worry, the code would not run this part is there were no candidate vertex.
			cout<<"\nSummary: add "<<track_IDs.size();
			cout<<" tracks (IDs: ";
			for(size_t i = 0; i < track_IDs.size(); i++){
				cout<<track_IDs[i]<<" ";
			}
			cout<<") to the candidate vertex: (";
			cout<<falgo.boundingSphere(points).Center().at(0)<<",";
			cout<<falgo.boundingSphere(points).Center().at(1)<<",";
			cout<<falgo.boundingSphere(points).Center().at(2)<<")"<<endl;
		}
		//move to next track end point, but actually, only one vertex candidate is considered here, so no more end point here;
		//Keng, current version takes only the best candidate vertex, maybe take all under consideration is better?
	}
	//  return;
}

double VertexBuilder::FindClosestApproach(
		const geoalgo::HalfLine_t & shr1,
		const geoalgo::HalfLine_t & shr2,
		geoalgo::Point_t & vtx) const {
  // Find mininum impact parameter between a two showers
  // flip their directions and look backwards

  // Create a half-line pointing backwards from the shower
  geoalgo::HalfLine_t shr1Bkwd(shr1.Start(), shr1.Dir() * (-1));
  geoalgo::HalfLine_t shr2Bkwd(shr2.Start(), shr2.Dir() * (-1));

  // coordinates for closest points on the two objects
  geoalgo::Point_t PtShr1(3);
  geoalgo::Point_t PtShr2(3);
  double IP = falgo.SqDist(shr1Bkwd, shr2Bkwd, PtShr1, PtShr2);//give the square of the cloesest distance in 3D
  //see https://nusoft.fnal.gov/larsoft/doxsvn/html/classgeoalgo_1_1GeoAlgo.html#a605d7e7a2736237727cad00c80c27eb8 for more;

  // build candidate vertex
  vtx = (PtShr1 + PtShr2) / 2.;

  return sqrt(IP);
}


void VertexBuilder::AssociateShowers(ParticleAssociations_all & pas) {
		
	//for creterias, cr1=10, cr2, cr3..
	//#1 load showers
	//#2 find precandidate
	//- distance shower_i to shower_j < cr1
	//- distance shower_i to track_k
	// ---> next shower_i
	//
	int screen_width = 86;
	bool first_time = true;
	bool first_time_2 = true;
	bool first_time_3 = true;

	std::map<size_t, Shower const *> shower_map;

	//#1 load showers
	for(size_t const i : fdetos->GetShowerIndices()) { 
		shower_map.emplace(i, &pas.GetDetectorObjects().GetShower(i));
	}

	std::vector<ParticleAssociation> const & associations = pas.GetAssociations();//Load vertex candidates, that contains vertices, radius of the vertex (bounding phere)

	while(shower_map.size()) {
		if(fverbose) {
			std::cout << "STEP I: Load " << shower_map.size() << " showers for evaluation.\n\n";

			if(shower_map.size()) {
				//Print the title
				cout<<setw(10)<<right<<"Shower ID ";
				cout<<"| Shower Start Point Coordinates (x,y,z) [cm] from Pandora ";
				cout<<"| Shower Direction from Pandora"<<endl;
				for(int i = 0; i<screen_width; i++) 
					cout<<"-";
				cout<<endl;

				for(auto const &c : shower_map){
					cout<<setw(10)<<left<< c.first<<"| ";//Index
					cout<<left<<c.second->fcone.Start();//Shower Coordinates
					cout<<" | "<<c.second->fcone.Dir()<<endl;//Direction
				}
				for(int i = 0; i<screen_width; i++) 
					cout<<"-";
				cout<<"\n\n STEP II:Look for vertex candidates among showers."<<endl;
			}
		}

  //#2 find precandidate among shower;
		size_t best_shower_id = SIZE_MAX;
		size_t best_other_id = SIZE_MAX;
		size_t index = SIZE_MAX;//updated when the vertex candidate is associated as the last action.
		bool addedToAnPas = false;//true means a shower is attached to a vertex candidate.

		Double_t best_dist = fshower_prox;//save the fshower_prox
		geoalgo::Point_t best_vert(2000, 2000, 2000);
		
		//#1 find precandidate vertex for any showers (shower or track)
//		for(auto const & c : shower_map) {//look at each shower candidate, 
		for(auto this_iterator = shower_map.begin(); this_iterator != shower_map.end(); ++ this_iterator) {//look at each shower candidate, 

			if(fverbose)
				std::cout << "\t\tLook at the shower (ID: " << this_iterator->first <<")"<< std::endl;

			geoalgo::Point_t const & c_start = this_iterator->second->fcone.Start();
			geoalgo::Vector_t const & c_dir = this_iterator->second->fcone.Dir();

		for(auto that_iterator = std::next(this_iterator); that_iterator != shower_map.end(); ++ that_iterator) {//look at each shower candidate, 
//			for(auto const & c2 : shower_map) {//compare showers
				if(fverbose) std::cout << "\t\t\tCompare to a shower (ID: " << that_iterator->first <<")"<< std::endl;

				if( that_iterator->first == this_iterator->first) {
					if(fverbose) std::cout << "\t\t\t\tSame shower, skip\n";
					continue;
				}
				geoalgo::Point_t temp_vert;

				double dist = FindClosestApproach( that_iterator->second->fcone, this_iterator->second->fcone, temp_vert);//minimum impact parameter
				if(first_time){
					f_dist_sx.erase(f_dist_sx.begin());
					first_time = false;
				}
				f_dist_sx.push_back(dist);
				//this is the shortest distance between two backward projection in 3D.
				//update temp_vert as the mid point of these two points that gives the shortest distance.

				if(fverbose)
					std::cout << "\t\t\tCompare dist of shower bkw-projections, " << dist << ", and the shortest-dist so far is "
						<< best_dist << ".\n";

				double temp_dist = best_dist;//fshower_prox now is temp_dist

				if(dist < temp_dist) {//min. impact parameter < max shower proximity threshold

					if(fverbose) std::cout << "\t\t\t\tUpdate the shortest-dist and the precandidate vertex!\n";

					best_shower_id = this_iterator->first;
					best_other_id = that_iterator->first;
					best_vert = temp_vert;
					best_dist = dist;
					index = SIZE_MAX;
					addedToAnPas = false;
				}
			}

			if(fverbose) std::cout << "\t\t >>Finish comparing showers to the shower with ID: "<<this_iterator->first<<endl;

			for(size_t const i : fdetos->GetTrackIndices()) {//compare tracks

				if(fverbose)
					std::cout << "\t\t\tCompare to a track with ID: " << i << std::endl;

				Track const & t = fdetos->GetTrack(i);

				geoalgo::Point_t dont_care;
				geoalgo::Point_t temp_vert;

				Double_t dist =
					sqrt(falgo.SqDist(t.ftrajectory, 
								geoalgo::HalfLine(c_start,
									c_dir*-1),
								temp_vert,//this point will be changed along the track; id does not matter what value it is.
								dont_care));//this is along the backward projection
				if(first_time){
					f_dist_sx.erase(f_dist_sx.begin());
					first_time = false;
				}
				f_dist_sx.push_back(dist);

				if(fverbose)
					std::cout << "\t\t\tCompare the impact parameter to the track, " << dist << ", and the shortest_dist so far is "
						<< best_dist << ".\n";

				if(dist < best_dist) {

					if(fverbose) std::cout << "\t\t\t\tUpdate the shortest-dist and the precandidate vertex!\n";

					best_shower_id = this_iterator->first;
					best_other_id = i;
					best_vert = temp_vert;
					best_dist = dist;
					index  = SIZE_MAX;
					addedToAnPas = false;
				}
			}

			if(fverbose) std::cout << "\t\t >>Finish comparing tracks to the shower with ID: "<<this_iterator->first<<endl;

			for(size_t i = 0; i < associations.size(); ++i) {//compare vertex candidates

				if(fverbose)
					std::cout << "\t\t\tCompare the impact parameter to the vertex candidate with ID: "
						<< i << std::endl;

				ParticleAssociation const & pa = associations.at(i);

				double dist = sqrt(falgo.SqDist(
								pa.GetRecoVertex(),
								geoalgo::HalfLine(c_start, c_dir*-1)));

				if(first_time){
					f_dist_sx.erase(f_dist_sx.begin());
					first_time = false;
				}
				f_dist_sx.push_back(dist);

				if(fverbose)
					std::cout << "\t\t\tCompare the impact parameter to the vertex candidate, " << dist << ", and the shortesest-dist so far, "
						<< best_dist << ".\n";

				if(dist < best_dist) {

					if(fverbose) std::cout << "\t\t\t\tUpdate the shortest-dist and the precandidate vertex!\n";

					best_shower_id = this_iterator->first;
					best_other_id = SIZE_MAX;//shower matches no shower nor track.
					best_vert = falgo.ClosestPt(pa.GetRecoVertex(), this_iterator->second->fcone);//c.second->fcone is a Trajectory_t
					best_dist = dist;
					index = i;//the ith association
					addedToAnPas = true;
				}
			}

			if(fverbose) std::cout << "\t\t >>Finish comparing vertex candidates to the shower with ID: "<<this_iterator->first<<endl;

		}// >>Finish looking at all showers with other showers, tracks, and vertex candidates.
		//it normally ends up with a shower with a small distance to another shower/track/vertex candidate


		if(fverbose) std::cout << "\n\tA brief report after evaluating all showers:\n"
			<< "\tThe shortest_dist to other objects: " << best_dist << ", which is required to be not further than "
				<< fshower_prox << "\n";

		if(best_dist >= fshower_prox) {//ds>=s_{max}
			if(fverbose) std::cout << "\t\tFail to satisfy the above condition, NO shower is going to be considered.\n";
			return;
		} else{
			if(fverbose) std::cout << "\t\tNext is to see if the found precandidate vertex qualified to be a vertex candidate.\n";
		}

//		if(fverbose) std::cout << "\tbest_shower_id: " << best_shower_id
//			<< " best_dist: " << best_dist
//				<< "\n\tindex: " << index << " == -1 ?\n";
		
		//#2 see if the precandidate be qualified as a vertex candidate
		if(addedToAnPas) { //the cloest object to the shower is a vertex candidate, so includes it into that vertex candidate (pas)

			pas.AddShower(index, best_shower_id, best_vert);

		} else{//the cloest object to the shower is a shower or a track, rather than a vertex candidate; 
		//the following determine whether to add this shower to a current vertex candidate with or without other objects,
		//	or just create a new vertex candidate for it.

			size_t association_index = SIZE_MAX;//a large number;
			double best_association_dist = fcpoa_vert_prox;

			if(fverbose) {
				std::cout << "\t\tCheck if we can connect the to an existent vertex candidate. ";
				std::cout << "Number of vertex candidates to be considered: " << associations.size() << std::endl;
			}

			for(size_t i = 0; i < associations.size(); ++i) {//find the cloest candidate vertex to the precandidate vertex


				double const dist = best_vert.Dist(associations.at(i).GetRecoVertex());
				if(first_time_2){
					f_dist_st.erase(f_dist_st.begin());
					first_time_2 = false;
				}
				f_dist_st.push_back(dist);

				if(fverbose){
					std::cout << "\t\t\tVertex candidate with ID: " << i << std::endl;
					std::cout << ", whose distance to the new precandidate vertex : " << dist
						<< "; compare this distance to the shortest_association_dist: "
						<< best_association_dist << std::endl;
				}
				if(dist < best_association_dist) {//<fcpoa_vert_prox, Max. distance btw shower start & cloest approach (dp_max);

					if(fverbose)
						std::cout << "\t\t\t\tThe shortest_association_dist is updated.\n";

					association_index = i;//get the index of association that gives the best association distance.
					best_association_dist = dist;

				}
			}
			
			//get the type (shower/track) of object that is cloest to the shower that we are currently looking at
			int const reco_type = fdetos->GetRecoType(best_other_id); //2 for track or 1 for shower;

			if(fverbose) std::cout << "\t\tThe precandidate is obtained from the shower with:\n";

			if(reco_type == fdetos->fshower_reco_type) {//the case of the shower associated to a shower

				if(fverbose)
					std::cout << "\t\t\ta shower\n"
						<< "\t\t\tCompare the shortest_association_dist, "<< best_association_dist
						<< ", which is required to be less than "
						<< fcpoa_vert_prox << ".\n";//max. shower projection distance;

				if(best_association_dist < fcpoa_vert_prox) {//this is a yes when the precandidate vertex is close to a candidate vertex
					if(fverbose) std::cout << "\t\t\t\tPass; promote the precandidate connecting two showers as a vertex candidate: "
						<< association_index << std::endl;
					pas.AddShower(association_index, best_shower_id, best_vert);
					pas.AddShower(association_index, best_other_id, best_vert);

				} else{//no candidate vertex around the precandidate vertex, then look for track endpoints instead.

					size_t best_track = SIZE_MAX;
					geoalgo::Point_t const * best_tp = nullptr;
					geoalgo::Point_t const * best_other_tp = nullptr;
					double best_showerp_dist = fcpoa_vert_prox;

					for(size_t const i : fdetos->GetTrackIndices()) {//look at all tracks;

						geoalgo::Trajectory const & t = fdetos->GetTrack(i).ftrajectory;

						//look at two end points of a track, and see if one of them is close to the precandidate vertex
						double const trackend_dist = t.back().Dist(best_vert);
						if(trackend_dist < best_showerp_dist) {
							best_track = i;
							best_tp = &t.back();
							best_other_tp = &t.front();
							best_showerp_dist = trackend_dist;
						}
						double const trackstart_dist = t.front().Dist(best_vert);
						if(trackstart_dist < best_showerp_dist) {
							best_track = i;
							best_tp = &t.front();
							best_other_tp = &t.back();
							best_showerp_dist = trackstart_dist;
						}

					}

					if(best_showerp_dist < fcpoa_vert_prox) {//it means at least one track end points has been used to connect the precandidate vertex

						std::vector<size_t> const & index_positions =
							pas.GetAssociationIndicesFromObject(best_track);//iterator that gives the best_track inside pas;
						switch(index_positions.size()){
							case 0://not an end point is associated to a candidate vertex;
								{//this is possible. OMG CHECK, this might mean the track is not an associated one.
									if(fverbose) {
										std::cout << "\t\t\t\tFail, create a new candiate vertex\n";
										if(best_tp == nullptr) std::cout << "The track-end point does not exist!\n";//tp = track point
									}

									std::vector<size_t> objects;
									objects.push_back(best_shower_id);
									objects.push_back(best_other_id);
									objects.push_back(best_track);
									std::vector<geoalgo::Point_t> verts(2, best_vert);//mid-point of two showers;

									verts.push_back(*best_tp);
									pas.AddAssociation(objects,
											verts,
											geoalgo::Sphere(*best_tp, best_dist).Center(),//bounding sphere
											best_dist);//"goodness"
								}
								break;
							case 1://1 end point is associated to a candidate vertex;
								{
									size_t const index =
										pas.GetAssociationIndices().at(index_positions.front());

									geoalgo::Point_t const & added_point =
										associations.at(index).GetVertexFromNode(best_track);//add the vertex candidate (obtained from two track end points) from pas;

									double const point_dist = added_point.Dist(*best_tp);//bounding sphere radius, distance btw track end point and vertex candidate;
									double const otherpoint_dist = added_point.Dist(*best_other_tp);//same as above, but with the other track end point;

									if(otherpoint_dist < point_dist) {//another end point is cloeser than the chosen one;

										if(associations.at(index).GetRecoVertex().
												Dist(*best_tp) < fstart_prox) {//mid point of both showers, track end point, and candidate are nearby.
												//if the vertex candidate is close enough to the track end point,
												//		then we add both showers (the track was already associated) to this candidate vertex;
											pas.AddShower(index, best_shower_id, best_vert);
											pas.AddShower(index, best_other_id, best_vert);

										} else {//the candidate vertex is far away from the track end point, 
												//		then make a new candidate vertex for these 2 showers and the track;
											if(fverbose) {
												std::cout << "\t\t\t\tFail, create a new candiate vertex\n";
												if(best_tp == nullptr) std::cout << "The track-end point does not exist!\n";//tp = track point
											}

											std::vector<size_t> showers;
											showers.push_back(best_shower_id);
											showers.push_back(best_other_id);
											showers.push_back(best_track);
											std::vector<geoalgo::Point_t> verts(2, best_vert);
											verts.push_back(*best_tp);
											pas.AddAssociation(showers,
													verts,
													geoalgo::Sphere(*best_tp, best_dist).Center(),
													best_dist);
										}
									} else {//this is cloesd enought; associate both showers to the candidate vertex
										pas.AddShower(index, best_shower_id, best_vert);
										pas.AddShower(index, best_other_id, best_vert);
									}
								}
								break;
							case 2://2 end points are associated to a candidate vertex;
								{//for the track end point that is closer to the candidate vertex, 
								//		add associate both showers respecting to that end point to the candidate vertex.
									size_t const indexa =
										pas.GetAssociationIndices().at(index_positions.front());
									double dista = 
										associations.at(indexa).GetRecoVertex().Dist(best_vert);

									size_t const indexb =
										pas.GetAssociationIndices().at(index_positions.back());   
									double distb = 
										associations.at(indexb).GetRecoVertex().Dist(best_vert);

									if(dista < distb) {
										pas.AddShower(indexa, best_shower_id, best_vert);
										pas.AddShower(indexa, best_other_id, best_vert);
									} else {
										pas.AddShower(indexb, best_shower_id, best_vert);      
										pas.AddShower(indexb, best_other_id, best_vert);
									}
								}
								break;
							default: {}
						}
					} else {//no track end points are near the precandidate vertex
						//Only add the (smallest impact parameter) point along shower bkw projection to the vertex candidate.
						std::vector<size_t> showers;
						showers.push_back(best_shower_id);
						showers.push_back(best_other_id);
						std::vector<geoalgo::Point_t> verts(2, best_vert);
						pas.AddAssociation(showers,
								verts,
								geoalgo::Sphere(best_vert, best_dist).Center(),
								best_dist);
					}
				}

				if(fverbose)
					std::cout << "\t\t\t >>Finish comparing to the shower with ID: " << best_other_id << std::endl;

				shower_map.erase(best_other_id);//When a shower is analyzed, eliminate it;

			} else if(reco_type == fdetos->ftrack_reco_type) {//the case that the precandidate vertex associates the shower to a track
			//if(fverbose) std::cout << "\t\tThe precandidate is obtained from the shower with:\n";

				if(fverbose) std::cout << "\t\t\ta track\n";

				geoalgo::Trajectory const & t = fdetos->GetTrack(best_other_id).ftrajectory;

				double best_trackend_dist = t.front().Dist(best_vert);
				geoalgo::Point_t const * point = &t.front();
				geoalgo::Point_t const * otherpoint = &t.back();

				double const trackend_dist = t.back().Dist(best_vert);

				if(first_time_3){
					f_dist_sst.erase(f_dist_sst.begin());
					first_time_3 = false;
				}
				f_dist_sst.push_back(trackend_dist);
				if(trackend_dist < best_trackend_dist) {//compare edges of q track to best_vert point distance
				//update the best if t.back() is the cloest to the best_vert
					best_trackend_dist = trackend_dist;
					point = &t.back();
					otherpoint = &t.front();
				}

				if(fverbose)
					std::cout << "\t\t\tbest_trackend_dist: "
						<< best_trackend_dist
						<< " < fcpoa_vert_prox: "
						<< fcpoa_vert_prox << " ?\n";

				if(best_trackend_dist < fcpoa_trackend_prox) {

					std::vector<size_t> const index_positions =
						pas.GetAssociationIndicesFromObject(best_other_id);
					//this helps to find the object with best_other_id as value
					// in the pas.

					if(fverbose)
						std::cout << "\t\t\t\tyes\n"
							<< "\t\t\t\tindex_positions.size(): "
							<< index_positions.size() << std::endl;

					std::vector<size_t> objects;
					std::vector<geoalgo::Point_t> verts;

					switch(index_positions.size()){
						case 0://the track has not been added to the pas;
							{
								if(fverbose)
									std::cout << "\t\t\t\t\tsize 0\n";

								objects.push_back(best_shower_id);
								objects.push_back(best_other_id);
								verts.push_back(best_vert);
								verts.push_back(*point);
								pas.AddAssociation(objects,
										verts,
										geoalgo::Sphere(*point, best_dist).Center(),
										best_dist);
							}
							break;
						case 1://the track has 1 end point added to the pas;
							{
								if(fverbose)
									std::cout << "\t\t\t\t\tsize 1\n";

								size_t const index =
									pas.GetAssociationIndices().at(index_positions.front());

								geoalgo::Point_t const & added_point =
									associations.at(index).GetVertexFromNode(best_other_id);

								double const point_dist = added_point.Dist(*point);
								double const otherpoint_dist = added_point.Dist(*otherpoint);

								if(fverbose)
									std::cout << "\t\t\t\t\totherpoint_dist: "
										<< otherpoint_dist
										<< " < point_dist: "
										<< point_dist << " ?\n";

								if(otherpoint_dist < point_dist) {

									if(fverbose)
										std::cout << "\t\t\t\t\t\tyes\n"
											<< "\t\t\t\t\t\tcenter_point_dist: "
											<< associations.at(index).GetRecoVertex().Dist(*point)
											<< " < fstart_prox: " << fstart_prox << " ?\n";

									if(associations.at(index).GetRecoVertex().
											Dist(*point) < fstart_prox) {
										if(fverbose) std::cout << "\t\t\t\t\t\t\tyes\n";
										pas.AddShower(index, best_shower_id, best_vert);
									} else {

										objects.push_back(best_shower_id);
										objects.push_back(best_other_id);
										//		std::vector<geoalgo::Point_t> verts;
										verts.push_back(best_vert);
										verts.push_back(*point);
										pas.AddAssociation(objects,
												verts,
												geoalgo::Sphere(*point, best_dist).Center(),
												best_dist);      
									}
								} else{

									if(fverbose)
										std::cout << "\t\t\t\t\t\tno\n"
											<< "\t\t\t\t\t\tadd id: " << best_shower_id
											<< " to association: " << index << std::endl;

									pas.AddShower(index, best_shower_id, best_vert); 
								}
							}
							break;
						case 2://the track is inside the pas twice
							{
								if(fverbose)
									std::cout << "\t\t\t\t\tsize 2\n";

								size_t const indexa =
									pas.GetAssociationIndices().at(index_positions.front());
								double dista = 
									associations.at(indexa).GetRecoVertex().Dist(best_vert);

								size_t const indexb = pas.GetAssociationIndices().at(index_positions.back());   
								double distb = 
									associations.at(indexb).GetRecoVertex().Dist(best_vert);

								if(dista < distb) {
									pas.AddShower(indexa,
											best_shower_id,
											best_vert);
								} else {
									pas.AddShower(indexb,
											best_shower_id,
											best_vert);
								}
							}
							break;
						default:
							std::cout << "Warning: more than two indices found, node: "
								<< best_other_id << std::endl;
					}
				} else {

					if(fverbose) std::cout << "\t\t\t\tno\n";

					pas.GetDetectorObjects().SetAssociated(best_shower_id);
				}
			}
		}

		shower_map.erase(best_shower_id);

	}

}


void VertexBuilder::AddLoneTracks(ParticleAssociations_all & pas) {
	//fdetos, detectorobjects_all
  for(size_t const gn : fdetos->GetTrackIndices()) {
  //goes over recob::track index, gn
	
	//Track is a struct defined in DetectorObjects.h
    Track const & t = fdetos->GetTrack(gn);//t is track

    if(t.fis_associated) continue;//this will jump to the end of the loop and goes with next gn; triger this when the track is associated.

    geoalgo::Point_t const * track_end = nullptr;
    double zmin = 2000;//the unit .. 2000 cm??? it doesnt matter though; it changes after first if-else

    geoalgo::Point_t const & front = t.ftrajectory.front();//first point along the trajectory.
    if(front.at(2) < zmin) {//first point's z coordinate compared to zmin
      track_end = &front;
      zmin = front.at(2);
    }

    geoalgo::Point_t const & back = t.ftrajectory.back();
    if(back.at(2) < zmin) {
      track_end = &back;
      zmin = back.at(2);
    }

    if(track_end) {

      pas.AddAssociation(std::vector<size_t>(1, gn),
			 std::vector<geoalgo::Point_t>(1, *track_end),
			 geoalgo::Sphere(*track_end, 0).Center(),
			 0);

	} else{
		std::cout << "Warning: No track end pointer\n";
	}
  }

}


void VertexBuilder::AddLoneShowers(ParticleAssociations_all & pas) {

  for(size_t const gn : fdetos->GetShowerIndices()) {

    Shower const & s = fdetos->GetShower(gn);

    if(s.fis_associated) continue;

    geoalgo::Point_t const p = s.fcone.Start();

    pas.AddAssociation(std::vector<size_t>(1, gn),
		       std::vector<geoalgo::Point_t>(1, p),
		       geoalgo::Sphere(p, 0).Center(),
		       0);

  }

}

// Fill the TTree with vertex info.
void VertexBuilder::FillVBT(ParticleAssociations_all & pas) {

  fvbt->ftrack_number = fdetos->GetTrackIndices().size();
  fvbt->fshower_number = fdetos->GetShowerIndices().size();
  fvbt->ftree->Fill();

}


void VertexBuilder::Run(ParticleAssociations_all & pas) {//Analysis the tracks & showers objects.

  CheckSetVariables();//Variables are set in SinglePhoton_module.cc.

  fdetos = &pas.GetDetectorObjects();//initialize an empty object;

//Make two associations, for tracks and showers.
  if(fverbose) std::cout << ">>>>> Associate tracks\n";

  AssociateTracks(pas);//Gives candidate vertex from tracks.

  //fvbt is the object, it means.. if the object is not empty, do something.
//  if(fvbt) fvbt->fassociation_track_number = pas.GetAssociations().size();

  if(fverbose) std::cout << "\n>>>> Associate showers\n";

  AssociateShowers(pas);//select associated candidates for showers

//if(fvbt) fvbt->fassociation_shower_number = pas.GetAssociations().size();

//CHECK, repeat for long objects? Add them into consideration even if they are not associated to a vertex?
  if(fverbose) std::cout << ">>>> Add lone tracks\n";
  AddLoneTracks(pas);
  if(fverbose) std::cout << ">>>> Add lone showers\n";
  AddLoneShowers(pas);

  if(fverbose) std::cout << ">>>> Get shower associations\n";
  pas.GetShowerAssociations();//CHECK

  if(fvbt) fvbt->fassociation_final_number = pas.GetSelectedAssociations().size();
  cout<<"\n\n # of Vertex found: "<<pas.GetSelectedAssociations().size()<<endl;

  if(fvbt) {//after the association is finished (found the vertex), fill in in the tree.
	  if(fverbose) std::cout << "Fill VBT\n";
	  FillVBT(pas);
  }

  pas.NodeCheck();

}
}


#endif

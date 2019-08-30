

#ifndef PARTICLEASSOCIATIONS_H
#define PARTICLEASSOCIATIONS_H

#include "../SinglePhoton_module.h"
#include "DetectorObjects.h"
#include <fstream>//read and write txt file.

using namespace std;

//what does f in fxxxx  mean ????

namespace single_photon{
class ParticleAssociations;


class ParticleAssociation{//Yes, it is with and without s!
	//private   
  friend ParticleAssociations;

  // this is a collection of indices and vertices,
  // appended in the AddShower1() function in this class.
  // Their elements come from AddShower() function in the class ParticleAssociations
  std::vector<size_t> findices;
  std::vector<geoalgo::Point_t> fvertices;


  geoalgo::Point_t fvertex;
  double fgoodness;
  
  std::multimap<size_t, size_t> fconnected_associations;
  
public:
//Constructor, for internal used only, I guess.
// indices (input) - to label 
  ParticleAssociation(std::vector<size_t> const & indices,
		      std::vector<geoalgo::Point_t> const & vertices,
		      geoalgo::Point_t const & vertex,
		      double const goodness);

  void AddConnection(size_t const i, size_t const n) {
    fconnected_associations.emplace(i, n);
  }
	//Note: this is not for adding pandora_object.. 
	// 2-argument function;
  void AddShower1(size_t const n, geoalgo::Point_t const & vert) {
    findices.push_back(n);
    fvertices.push_back(vert);
  }
  
  std::vector<size_t> const & GetObjectIndices() const {
    return findices;
  }
  
  std::vector<geoalgo::Point_t> const & GetVertices() const {
    return fvertices;
  }
  
  geoalgo::Point_t const & GetVertexFromNode(size_t const n) const {
    
    auto const nv_itb = findices.begin();
    auto const nv_ite = findices.end();
    
    return fvertices.at(std::find(nv_itb, nv_ite, n) - nv_itb);
    
  }
  
  geoalgo::Point_t const & GetRecoVertex() const {
    return fvertex;
  }
  
  double GetGoodness() const {
    return fgoodness;
  }
  
  std::multimap<size_t, size_t> const & GetConnections() const {
    return fconnected_associations;
  }
  
  void PrintAssociation() const;
  
};


class ParticleAssociations{

  friend class SinglePhoton;//OK~ now SinglePhoton can use things in this class!
  DetectorObjects fdetos;//here introduces the DetectorObjects.h

  std::vector<ParticleAssociation> fassociations;
  std::vector<size_t> fobject_index_v;
  std::vector<size_t> fassociation_index_vec;
  std::vector<size_t> fignore_association_vec;
  std::vector<size_t> fselected_associations;//Indices of selected associations (candidate vertice). Done under GetShowerAssociations()

  bool fverbose;

public:

  ParticleAssociations();

  DetectorObjects & GetDetectorObjects() {return fdetos;}//gives out the DetectorObjects.
  DetectorObjects const & GetDetectorObjects() const {return fdetos;}

  void Reset();

  void AddAssociation(std::vector<size_t> const & nodes,
		      std::vector<geoalgo::Point_t> const & vertices,
		      geoalgo::Point_t const & vertex,
		      double const goodness = 0);
    

  /*************************
   *
   * AddShower() - 
   * index - this is only for the fassociation_index_vec; 
   *	I guess, this is to append the nth association.
   * n - this coud be the index (ID) for a pandora_object.
   * vert - a pandora vertex.
   *
   * ***********************/
  void AddShower(size_t const index,
		 size_t const n,
		 geoalgo::Point_t const & vert);//used 

  std::vector<ParticleAssociation> const & GetAssociations() const {
    return fassociations;
  }

  std::vector<size_t> const & GetObjectIndices() const {
    return fobject_index_v;
  }
  
  std::vector<size_t> const & GetAssociationIndices() const {
    return fassociation_index_vec;
  }

  std::vector<size_t> GetAssociationIndicesFromObject(size_t const n);

  void IgnoreAssociation(size_t const n) {//CHECK, when do we use IgnoreAssociation??
	cout<<"\n\n\n\n CHECK Here is where I use the IgnoreAssociation"<<endl;
    fignore_association_vec.push_back(n);
  }
//  bool Ignore(size_t const i) const;
  //bool HandleLoop(std::vector<size_t> const & previously_considered);
  void IgnoreThis(size_t const to_ignore, size_t const connected_index, std::vector<size_t> & previously_considered);


  /*********************
   *
   * ******************************/
  void IgnoreAssociationsConnectedTo(size_t const i);

  void ClearIgnored() {
    fignore_association_vec.clear();
  }

  void PrintAssociations() const;
  void PrintAssociations(std::vector<size_t> const & associations_to_print) const;

  void PrintNodes() const;

  void Test() const;

  void NodeCheck();

  void GetShowerAssociations();

  std::vector<size_t> & GetSelectedAssociations() {return fselected_associations;}
  std::vector<size_t> const & GetSelectedAssociations() const {return fselected_associations;}

  //void DeleteAssociation(size_t const s);

  void SetVerbose(bool const verbose = true) {fverbose = verbose;}

};






//THE HEADER FILE IS THE ABOVE










ParticleAssociation::ParticleAssociation(std::vector<size_t> const & indices,
					 std::vector<geoalgo::Point_t> const & vertices,
					 geoalgo::Point_t const & vertex,
					 Double_t const goodness) :
  findices(indices),
  fvertices(vertices),
  fvertex(vertex),
  fgoodness(goodness){}


void ParticleAssociation::PrintAssociation() const {
  
  for(size_t i = 0; i < findices.size(); ++i)
    std::cout << "Object index: " << findices.at(i) << std::endl;
  
  std::cout << "\nVertex: " << fvertex << ", goodness: "
	    << fgoodness << std::endl;
  
  std::cout << "\nConnections:\n";
  for(std::pair<size_t, size_t> const & pair : fconnected_associations) 
    std::cout << "Association index: " << pair.first << " object index: " << pair.second << std::endl;
  
}




ParticleAssociations::ParticleAssociations() :
  fverbose(false){}//DAMN, : fverbose(false) was misunderstood as the constructor in SinglePhoton..


void ParticleAssociations::Reset() {

  fassociations.clear();
  fobject_index_v.clear();
  fassociation_index_vec.clear();

}


void ParticleAssociations::AddAssociation(std::vector<size_t> const & nodes,
					  std::vector<geoalgo::Point_t> const & vertices,
					  geoalgo::Point_t const & vertex,
					  Double_t const goodness) {
  
  fassociations.push_back(ParticleAssociation(nodes, vertices, vertex, goodness));
  
  size_t const pos = fassociations.size() - 1;
  
  for(size_t const n : nodes) {
    
    auto const nv_itb = fobject_index_v.begin();
    auto const nv_ite = fobject_index_v.end();
    
    for(auto nv_it = nv_itb; nv_it != nv_ite; ++nv_it) {
      
      nv_it = std::find(nv_it, nv_ite, n);
      
      if(nv_it != nv_ite) {
	
	size_t const index = fassociation_index_vec.at(nv_it - nv_itb);
	
	fassociations.at(index).AddConnection(pos, n);
	fassociations.at(pos).AddConnection(index, n);
	
      }
      
      else break;
      
    }
    
    fobject_index_v.push_back(n);
    fassociation_index_vec.push_back(pos);
    
  }      

  for(size_t const i : nodes) {  
    fdetos.SetAssociated(i);
  }

  
}


void ParticleAssociations::AddShower(size_t const index,
				     size_t const n,
				     geoalgo::Point_t const & vert) {
  
  if(index > fassociations.size() - 1) {
    std::cout << "No association with index: " << index << std::endl;
    return;
  }
  
  //Dont know why consruct a ParticleAssociation here..
  ParticleAssociation & association = fassociations.at(index);
  std::vector<size_t> const & group = association.GetObjectIndices();
  
  if(std::find(group.begin(), group.end(), n) != group.end()) {
    std::cout << "Object ID: " << n << " already added\n";
    return;
  }
  
  fassociation_index_vec.push_back(index);

  fobject_index_v.push_back(n);
  
  association.AddShower1(n, vert);//this is for single association; Only Here, the 2nd and 3rd argument from the AddShower();

  fdetos.SetAssociated(n);
  
}


std::vector<size_t> ParticleAssociations::GetAssociationIndicesFromObject(size_t const n) {

  auto const nv_itb = fobject_index_v.begin();
  auto const nv_ite = fobject_index_v.end();
  auto nv_it = nv_itb;
  
  std::vector<size_t> indices;
  
  while(nv_it != nv_ite) {
    
    nv_it = std::find(nv_it, nv_ite, n);
    
    if(nv_it != nv_ite) {
      indices.push_back(nv_it - nv_itb);
      ++nv_it;
    }
    
  }
  
  return indices;
  
}

/*CHECK
bool ParticleAssociations::Ignore(size_t const i) const {
  auto const sk_it =
    std::find(fignore_association_vec.begin(), fignore_association_vec.end(), i);//return the iterator who points to the number i.
  if(sk_it == fignore_association_vec.end()) {
		cout<<" Ok.. Ignore() is going to return false CHECK)"<<endl;
		//Null==Null??
	  return false; //it seems like this will never be truth.
  }
  else{
	  return true;
  }
}
*/

/*
bool ParticleAssociations::HandleLoop(std::vector<size_t> const & previously_considered) {

  if(fverbose) std::cout << "Handle loop\n";

  size_t delete_association = SIZE_MAX;
  double largest_sphere = -1;
  for(size_t const s : previously_considered) {
    ParticleAssociation const & pa = fassociations.at(s);
    if(pa.GetObjectIndices().size() == 2 && pa.GetConnections().size() == 2 && pa.GetGoodness() > largest_sphere) {
      delete_association = s;
      largest_sphere = pa.GetGoodness();
    }
  }

  if(delete_association == SIZE_MAX) return false;
  DeleteAssociation(delete_association);
  return true;

}
*/


void ParticleAssociations::IgnoreThis(size_t const to_ignore, size_t const connected_index, std::vector<size_t> & previously_considered) {

  if(std::find(previously_considered.begin(), previously_considered.end(), to_ignore) != previously_considered.end()) return;
  
  previously_considered.push_back(to_ignore);
  ParticleAssociation const & pa = fassociations.at(to_ignore);
  fignore_association_vec.push_back(to_ignore);
  
  for(std::pair<size_t, size_t> const & p : pa.GetConnections()) {
    if(p.first == connected_index) continue;
    IgnoreThis(p.first, to_ignore, previously_considered);
  }

}



void ParticleAssociations::IgnoreAssociationsConnectedTo(size_t const i) {
  
  ParticleAssociation const & pa = fassociations.at(i);
  
  for(std::pair<size_t, size_t> const & p : pa.GetConnections()) {
    std::vector<size_t> previously_considered;  
    IgnoreThis(p.first, i, previously_considered);
  }
  
}


void ParticleAssociations::PrintAssociations() const {

  std::cout << "----------------------------------------\n\n";
  
  for(size_t i = 0; i < fassociations.size(); ++i) {
 
    std::cout << "Association: " << i << std::endl;
    
    fassociations.at(i).PrintAssociation();
    
    std::cout << std::endl;
    
  }
  
  std::cout << std::endl;
  
}


void ParticleAssociations::PrintAssociations(std::vector<size_t> const & associations_to_print) const {

  std::cout << "----------------------------------------\n\n";
  
  for(size_t i = 0; i < fassociations.size(); ++i) {
    
    if(std::find(associations_to_print.begin(), associations_to_print.end(), i) == associations_to_print.end()) continue;

    std::cout << "Association: " << i << std::endl;
    
    fassociations.at(i).PrintAssociation();
    
    std::cout << std::endl;
    
  }
  
  std::cout << std::endl;
  
}


void ParticleAssociations::PrintNodes() const  {
  
  std::cout << "Nodes:\n";
  
  for(size_t i = 0; i < fobject_index_v.size(); ++i)
    std::cout << "Nodes: " << fobject_index_v.at(i) << " Index: "
	      << fassociation_index_vec.at(i) << std::endl;
  
}


void ParticleAssociations::Test() const {

  for(ParticleAssociation const & pae : fassociations) {
    
    for(std::pair<size_t, size_t> const & pair : pae.GetConnections()) {
      
      size_t const i = pair.first;
      
      std::vector<size_t> const & c = fassociations.at(i).GetObjectIndices();
      auto c_itb = c.begin();
      auto c_ite = c.end();
      
      auto c_it = c_ite;
      
      for(size_t const n : pae.GetObjectIndices()) {
	
	c_it = std::find(c_itb, c_ite, n);
	
	if(c_it != c_ite) break;
	
      }
      
      if(c_it == c_ite) std::cout << "Complain\n";
      
    }
    
  }
  
}


void ParticleAssociations::NodeCheck() {

  std::vector<size_t> nodes;

  for(size_t const n : fobject_index_v) {
    if(nodes.size() < n+1) nodes.resize(n+1, 0);    
    ++nodes.at(n);
  }

  for(size_t i = 0; i < nodes.size(); ++i) {

    size_t const s = nodes.at(i);
    if(s == 0) continue;

    DetectorObject const & deto = fdetos.GetDetectorObject(i);

    if(deto.freco_type == fdetos.ftrack_reco_type && s > 2)
      std::cout << "track > 2 entries: " << s << std::endl;
    if(deto.freco_type == fdetos.fshower_reco_type && s > 1)
      std::cout << "shower > 1 entries: " << s << std::endl;

  }
  
}


void ParticleAssociations::GetShowerAssociations() {
	if(fverbose) std::cout << "GetShowerAssociations\n";

	std::multimap<double, size_t> pa_map;//Partical Association Map; it maps vertex index to the vertex's z-coordinate respectiely.

//	int temp_num_showers = 0;
//	int temp_num_tracks = 0;

	if(fverbose) std::cout << "Number of particle associations: " << fassociations.size() << "\n";

	for(size_t i = 0; i < fassociations.size(); ++i) {
		if(fverbose) std::cout << "\tAssociation: " << i << "\n";
		ParticleAssociation const & pa = fassociations.at(i); //the ith associated particle (vertex)
		bool asso_shower = false;
		for(size_t const s : pa.GetObjectIndices()) {
			if(fverbose) std::cout << "\t\tObject index: " << s << "\n";
			if(fdetos.GetRecoType(s) == fdetos.fshower_reco_type) {//fdetos is detector objects
//				temp_num_showers++;//Keng
				if(fverbose) std::cout << "\t\tis shower\n";


				asso_shower = true;
				break;
			}
		}
		if(asso_shower) {
			pa_map.emplace(pa.GetRecoVertex().at(2), i);//
		}
	}

//	std::ofstream save_shower("temp_num_showertrack.txt");//recored limits of boundary;
//
////	m_bobbyshowers = temp_num_showers;
//
//	save_shower<< temp_num_showers <<" "<< temp_num_tracks<<std::endl;
//	save_shower.close();

/* the old if-else condition;
bool ParticleAssociations::Ignore(size_t const i) const {
  auto const sk_it =
    std::find(fignore_association_vec.begin(), fignore_association_vec.end(), i);//return the iterator who points to the number i.
  if(sk_it == fignore_association_vec.end()) {
		cout<<" Ok.. Ignore() is going to return false CHECK)"<<endl;
		//Null==Null??
	  return false; //it seems like this will never be truth.
  }
  else{
	  return true;
  }
}
 */
	if(fverbose) std::cout << "Loop over filled map\n";

	for(std::pair<double, size_t> const & p : pa_map) {
		if(fverbose) std::cout << "\tAssociation: " << p.second << " " << "z-position: " << p.first << "\n"; 

		auto const sk_it = std::find(fignore_association_vec.begin(), fignore_association_vec.end(), p.second);//return the iterator who points to the number i.
		if( sk_it != fignore_association_vec.end() ){ //Ignore(the association index)
			cout<<"Ignore(p.second) is true"<<endl;
			continue;//proceed when the Ignore() gives false; it means sk_it == fignore_association_vec.end();
		}
		//if sk_it is not the last fignore_associatino_vec, then do the follwoing
		cout<<"CHECK I do Have a list to ignore "<<fignore_association_vec.size()<<endl;
		IgnoreAssociationsConnectedTo(p.second);
		fselected_associations.push_back(p.second);//CHECK, here is where we select the vertex.
	}

	if(fverbose) std::cout << "ClearIgnored\n";

	ClearIgnored();

}


/*
void ParticleAssociations::DeleteAssociation(size_t const s) {
  
  if(fverbose) {
    std::cout << "Delete association " << s << "\nDelete other association connections to this association\n";
  }

  std::vector<size_t> considered_connected_associations;
  for(std::pair<size_t, size_t> const & p : fassociations.at(s).GetConnections()) {  
    if(std::find(considered_connected_associations.begin(), considered_connected_associations.end(), p.first) != considered_connected_associations.end()) continue;
    std::multimap<size_t, size_t> & connections = fassociations.at(p.first).fconnected_associations;
    auto c_it = connections.find(s);
    if(c_it == connections.end()) {
      std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nNo connections found\n";
      exit(1);
    }
    while(c_it != connections.end()) {
      connections.erase(c_it);
      c_it = connections.find(s);
    }
    considered_connected_associations.push_back(p.first);
  }

  if(fverbose) std::cout << "Delete association in ParticleAssociations vector\n";

  fassociations.erase(fassociations.begin()+s);

  if(fverbose) std::cout << "Delete ParticleAssociations bookeeping of association and object indices\n";

  auto av_it = std::find(fassociation_index_vec.begin(), fassociation_index_vec.end(), s);
  while(av_it != fassociation_index_vec.end()) {
    fobject_index_v.erase(fobject_index_v.begin()+size_t(av_it-fassociation_index_vec.begin()));
    fassociation_index_vec.erase(av_it);
  }

  if(fverbose) {
    std::cout << "End delete association\n";
  }
    
}
*/
}


#endif

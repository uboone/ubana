////////////////////////////////////////////////////////////////////////
// Class:       PhotoNuclear
// Plugin Type: filter (art v3_01_02)
// File:        PhotoNuclear_module.cc
//
// Generated at Mon Feb 10 12:04:26 2020 by Yeon-jae Jwa using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

// LArSoft includes                                                                                                                                                      

#include "canvas/Persistency/Common/FindManyP.h"
 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"


class PhotoNuclear;


class PhotoNuclear : public art::EDFilter {
public:
  explicit PhotoNuclear(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotoNuclear(PhotoNuclear const&) = delete;
  PhotoNuclear(PhotoNuclear&&) = delete;
  PhotoNuclear& operator=(PhotoNuclear const&) = delete;
  PhotoNuclear& operator=(PhotoNuclear&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  //void beginJob() override;
  //void endJob() override;

private:

  // Declare member data here.

};


PhotoNuclear::PhotoNuclear(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool PhotoNuclear::filter(art::Event& e)
{
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel("generator", truthHandle);

  const art::FindManyP<simb::MCParticle> truthParticles(truthHandle, e, "largeant");
  assert(truthParticles.isValid());

  bool photonu_absorbed = false;

  for (size_t itruth=0; itruth<truthParticles.size(); itruth++) {

    auto const& mcparticles = truthParticles.at(itruth);
    int photon_count = 0;
    int primary_photon_count = 0;
    std::vector<int> pi0_trackIDs;

    // collecting pi0s in mcparticles
    for (size_t i=0; i<mcparticles.size(); i++) {

      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();
      if (pdg ==111){
	//std::cout <<"pi0 :: EndProcess : " << p.EndProcess() << " , Statuscode : " << p.StatusCode() << " , TrackID : " << p.TrackId() << " , P : " << p.P() << std::endl;
        pi0_trackIDs.push_back(p.TrackId());
      }
    }
    //std::cout << "size of pi0 set : " << pi0_trackIDs.size() << std::endl;

    for (size_t i=0; i<mcparticles.size(); i++) {
      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();
      bool is_from_pi0=false;
      if (pdg == 22){
	//    std::cout <<itruth<<"-th truth, "<< photon_count<<"-th photon :: EndProcess : " << p.EndProcess() << " , Mother: " << p.Mother() << " , Statuscode : " << p.StatusCode() << " , TrackID : " << p.TrackId() << " , P : " << p.P() << std::endl;
        photon_count++;
        if (std::find(pi0_trackIDs.begin(), pi0_trackIDs.end(),p.Mother())!=pi0_trackIDs.end()){
          is_from_pi0 = true;
          primary_photon_count++;
	  //	  std::cout <<itruth<<"-th truth, "<< primary_photon_count<<"-th primary photon :: EndProcess : " << p.EndProcess() << " , Mother: " << p.Mother() << " , Statuscode : " << p.StatusCode() << " , TrackID : " << p.TrackId() << " , P : " << p.P() << " , primary? :" <<is_from_pi0<<std::endl;
        }
      }
      if (!is_from_pi0) continue;

      std::string endProc = p.EndProcess();
      photonu_absorbed = (endProc.find("photonNuclear") != std::string::npos);
      if (photonu_absorbed) {
	std::cout << "gonna return true" <<std::endl;
	return true;
      }
    }
  }
  std::cout << "gonna return false" << std::endl;
  return photonu_absorbed;
  
}


DEFINE_ART_MODULE(PhotoNuclear)

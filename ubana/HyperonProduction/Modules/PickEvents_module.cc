////////////////////////////////////////////////////////////////////////
// Class:       PickEvents
// Plugin Type: filter (art v3_01_02)
// File:        PickEvents_module.cc
//
// Generated at Mon Nov 23 05:39:28 2020 by Christopher Thorpe using cetskelgen
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"

//local includes

//#include "../Alg/FV.h" //fiducial volume checker


namespace hyperon {
	class PickEvents;
}


class hyperon::PickEvents : public art::EDFilter {
	public:
		explicit PickEvents(fhicl::ParameterSet const& p);
		// The compiler-generated destructor is fine for non-base
		// classes without bare pointers or other resource use.

		// Plugins should not be copied or assigned.
		PickEvents(PickEvents const&) = delete;
		PickEvents(PickEvents&&) = delete;
		PickEvents& operator=(PickEvents const&) = delete;
		PickEvents& operator=(PickEvents&&) = delete;

		// Required functions.
		bool filter(art::Event& e) override;

		// Selected optional functions.
		void beginJob() override;
		void endJob() override;

	private:

		std::string eventlistfilename;
		std::ifstream* fInputEventList;

		int run,subrun,event;

		int nfound;
		int nprocessed;

		std::vector<int> runToUse;
		std::vector<int> subrunToUse;
		std::vector<int> eventToUse;

};


hyperon::PickEvents::PickEvents(fhicl::ParameterSet const& p)
	: EDFilter{p}
	, eventlistfilename{p.get<std::string>("EventList")}  // ,
	// More initializers here.
{
	
}

bool hyperon::PickEvents::filter(art::Event& e)
{

	//get the run/subrun/event numbers
	run = e.run();
	subrun = e.subRun();
	event = e.event();

	//true if event was in the list 
	bool found=false;

	//seach lists of corresponding run/sub/event
	for(size_t i=0;i<runToUse.size();i++){

		if( run == runToUse.at(i) && subrun == subrunToUse.at(i) && event == eventToUse.at(i) ) found = true;

	}

	if(found){ nfound++; std::cout << "Got " << nfound << " events of " << runToUse.size() << std::endl; }

	nprocessed++;
	if( nprocessed % 1000 == 0 ) std::cout << "Processed " << nprocessed << " events" << std::endl;

	return found;
}

void hyperon::PickEvents::beginJob()
{
	
	nfound=0;
	nprocessed=0;
	
	// Implementation of optional member function here.
	fInputEventList = new std::ifstream(eventlistfilename,std::fstream::in);	

	//read all the data into vectors
	std::cout << "Reading runs,subruns and events from text file" << std::endl;	
	int thisRun,thisSubrun,thisEvent;
	
	std::string line;
	std::istringstream inputLine;


	while( std::getline(*fInputEventList,line) ){

		inputLine.clear();
		inputLine.str(line);


	inputLine >> thisRun >> thisSubrun >> thisEvent;
	std::cout << thisRun << " " << thisSubrun << " " << thisEvent << std::endl;

	runToUse.push_back(thisRun);
	subrunToUse.push_back(thisSubrun);
	eventToUse.push_back(thisEvent);

	if( runToUse.size() >= 500 ){
	std::cout << "Limit of 500 events stored, not adding any more" << std::endl;
	break;
	}

	}

		
	

}

void hyperon::PickEvents::endJob()
{
	// Implementation of optional member function here.
}

DEFINE_ART_MODULE(hyperon::PickEvents)

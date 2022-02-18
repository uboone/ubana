////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        MyAnalysis_module.cc
//
// Generated at Wed Jul 14 15:15:04 2021 by Ohana Rodrigues using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <map>
#include <string>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialCutsCouple.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"
#include "larcore/Geometry/Geometry.h"



class MyAnalysis;

class MyAnalysis : public art::EDAnalyzer {
public:
  explicit MyAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyAnalysis(MyAnalysis const&) = delete;
  MyAnalysis(MyAnalysis&&) = delete;
  MyAnalysis& operator=(MyAnalysis const&) = delete;
  MyAnalysis& operator=(MyAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  int p_PDG;
  std::string fMCParticleProducer;
  art::ServiceHandle<geo::Geometry> fGeometryService;
  TH1F* fHist;  //!< Output histogram

};


MyAnalysis::MyAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  
  //fMCParticleProducer = p.get<std::string>("MCParticleProducer", "largeant");
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void MyAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
   // Get MCParticles for each MCTruth in this event
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel("generator", truthHandle);
  const art::FindManyP<simb::MCParticle> truthParticles(truthHandle, e, fMCParticleProducer);
  assert(truthParticles.isValid());

  // Loop over sets of MCTruth-associated particles
  for (size_t itruth=0; itruth<truthParticles.size(); itruth++) {
    std::cout << "New truthParticle---" << std::endl;
    // Loop over MCParticles in the event
    auto const& mcparticles = truthParticles.at(itruth);
    for (size_t i=0; i<mcparticles.size(); i++) {
      // Get the MC Particle and get some information about the true particle
      const simb::MCParticle& p = *mcparticles.at(i);
      p_PDG = p.PdgCode();
      // int mcpID = p.TrackId(); unused variable
      std::string EndProcess  = p.EndProcess();
      // double mass = 0.; 
      //      if( TMath::Abs(p_PDG) == 211 ) mass = 139.57;
      //          else if( p_PDG == 2212 ) mass = 938.28;
      //Get the list of processes from the true trajectory
      const std::vector< std::pair< size_t, unsigned char > > & processes = p.Trajectory().TrajectoryProcesses();
      std::map< size_t, std::string > process_map;
      for( auto it = processes.begin(); it != processes.end(); ++it ){
        process_map[ it->first ] = p.Trajectory().KeyToProcess( it->second ); 
      }
      for( size_t j = 0; j < p.NumberTrajectoryPoints(); ++j ){
        double X = p.Position(j).X();
        double Y = p.Position(j).Y();
        double Z = p.Position(j).Z();
        geo::Point_t testpoint1 { X, Y, Z };
        const TGeoMaterial* testmaterial1 = fGeometryService->Material( testpoint1 );
        //For now, just going to reweight the points within the LAr of the TPC
        // TODO check if this is right
        if ( !strcmp( testmaterial1->GetName(), "LAr" ) ){
          auto itProc = process_map.find(i);
          std::cout << "Process " << itProc->second << std::endl;
          if( itProc != process_map.end() && itProc->second == "hadElastic" ){
            //Push back the index relative to the start of the reweightable steps
            std::cout<<"Process"<<itProc->second<<"is hadElastic";
            // if (fDebug) std::cout << "Elastic index: " << trajpoint_X.size() - 1 << std::endl;
          }
        }
      } // end loop over trajectory points
    }// end loop over mcparticles
  }// end loop over truthParticles
}// close analyze function

void MyAnalysis::beginJob()
{
  // Implementation of optional member function here.
}

void MyAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MyAnalysis)

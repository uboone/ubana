////////////////////////////////////////////////////////////////////////
// Class:       MCRecoCheck
// Plugin Type: analyzer (art v3_00_00)
// File:        MCRecoCheck_module.cc
//
// Generated at Wed Jan  2 09:53:19 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TTree.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/MCBase/MCShower.h" 
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"

class MCRecoCheck;


class MCRecoCheck : public art::EDAnalyzer {
public:
  explicit MCRecoCheck(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCRecoCheck(MCRecoCheck const&) = delete;
  MCRecoCheck(MCRecoCheck&&) = delete;
  MCRecoCheck& operator=(MCRecoCheck const&) = delete;
  MCRecoCheck& operator=(MCRecoCheck&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _tree;
  int    _nmcshr; // number of mcshowers in the event
  double _xtimeoffset;
  double _xsceoffset_start, _ysceoffset_start, _zsceoffset_start;
  double _xsceoffset_edep, _ysceoffset_edep, _zsceoffset_edep;
  double _mcpart_x_start, _mcpart_y_start, _mcpart_z_start, _mcpart_e_start;
  double _mcpart_x_edep, _mcpart_y_edep, _mcpart_z_edep, _mcpart_e_edep, _mcpart_q_edep;
  double _mcshr_x_start, _mcshr_y_start, _mcshr_z_start, _mcshr_e_start;
  double _mcshr_x_edep, _mcshr_y_edep, _mcshr_z_edep, _mcshr_e_edep;
  double _mcshr_q0_edep, _mcshr_q1_edep, _mcshr_q2_edep;

};


MCRecoCheck::MCRecoCheck(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","tree");
  _tree->Branch("_nmcshr",&_nmcshr,"nmcshr/I");
  _tree->Branch("_xtimeoffset",&_xtimeoffset,"xtimeoffset/D");
  _tree->Branch("_xsceoffset_start",&_xsceoffset_start,"xsceoffset_start/D");
  _tree->Branch("_ysceoffset_start",&_ysceoffset_start,"ysceoffset_start/D");
  _tree->Branch("_zsceoffset_start",&_zsceoffset_start,"zsceoffset_start/D");
  _tree->Branch("_xsceoffset_edep",&_xsceoffset_edep,"xsceoffset_edep/D");
  _tree->Branch("_ysceoffset_edep",&_ysceoffset_edep,"ysceoffset_edep/D");
  _tree->Branch("_zsceoffset_edep",&_zsceoffset_edep,"zsceoffset_edep/D");
  _tree->Branch("_mcpart_x_start",&_mcpart_x_start,"mcpart_x_start/D");
  _tree->Branch("_mcpart_y_start",&_mcpart_y_start,"mcpart_y_start/D");
  _tree->Branch("_mcpart_z_start",&_mcpart_z_start,"mcpart_z_start/D");
  _tree->Branch("_mcpart_e_start",&_mcpart_e_start,"mcpart_e_start/D");
  _tree->Branch("_mcpart_x_edep",&_mcpart_x_edep,"mcpart_x_edep/D");
  _tree->Branch("_mcpart_y_edep",&_mcpart_y_edep,"mcpart_y_edep/D");
  _tree->Branch("_mcpart_z_edep",&_mcpart_z_edep,"mcpart_z_edep/D");
  _tree->Branch("_mcpart_e_edep",&_mcpart_e_edep,"mcpart_e_edep/D");
  _tree->Branch("_mcpart_q_edep",&_mcpart_q_edep,"mcpart_q_edep/D");
  _tree->Branch("_mcshr_x_start",&_mcshr_x_start,"mcshr_x_start/D");
  _tree->Branch("_mcshr_y_start",&_mcshr_y_start,"mcshr_y_start/D");
  _tree->Branch("_mcshr_z_start",&_mcshr_z_start,"mcshr_z_start/D");
  _tree->Branch("_mcshr_e_start",&_mcshr_e_start,"mcshr_e_start/D");
  _tree->Branch("_mcshr_x_edep",&_mcshr_x_edep,"mcshr_x_edep/D");
  _tree->Branch("_mcshr_y_edep",&_mcshr_y_edep,"mcshr_y_edep/D");
  _tree->Branch("_mcshr_z_edep",&_mcshr_z_edep,"mcshr_z_edep/D");
  _tree->Branch("_mcshr_e_edep",&_mcshr_e_edep,"mcshr_e_edep/D");
  _tree->Branch("_mcshr_q0_edep",&_mcshr_q0_edep,"mcshr_q0_edep/D");
  _tree->Branch("_mcshr_q1_edep",&_mcshr_q1_edep,"mcshr_q1_edep/D");
  _tree->Branch("_mcshr_q2_edep",&_mcshr_q2_edep,"mcshr_q2_edep/D");

}

void MCRecoCheck::analyze(art::Event const& e)
{

  // load MCShowers
  auto const& mcshr_h = e.getValidHandle<std::vector<sim::MCShower> >("mcreco");
  _nmcshr = mcshr_h->size();
  if (mcshr_h->size() != 1) {
    std::cout << "ERROR : number of MCShowers is " << _nmcshr << std::endl;
    _tree->Fill();
    return;
  }

  auto mcshower = mcshr_h->at(0);
  
  _mcshr_x_start = mcshower.Start().X();
  _mcshr_y_start = mcshower.Start().Y();
  _mcshr_z_start = mcshower.Start().Z();
  _mcshr_e_start = mcshower.Start().E();
  _mcshr_x_edep  = mcshower.DetProfile().X();
  _mcshr_y_edep  = mcshower.DetProfile().Y();
  _mcshr_z_edep  = mcshower.DetProfile().Z();
  _mcshr_e_edep  = mcshower.DetProfile().E();
  _mcshr_q0_edep  = mcshower.Charge(0);
  _mcshr_q1_edep  = mcshower.Charge(1);
  _mcshr_q2_edep  = mcshower.Charge(2);

  // load mcparticles
  auto const& mcpart_h = e.getValidHandle<std::vector<simb::MCParticle> > ("largeant");
  // find primary mcparticle [ for single particle files ]
  size_t showerpartidx = 0;
  size_t showertrackid = 0;
  size_t nprimary = 0;
  for (size_t p=0; p < mcpart_h->size(); p++) {
    auto mcpart = mcpart_h->at(p);
    //showertrackid = mcpart.TrackId();
    //auto mother  = mcpart.Mother();
    auto process = mcpart.Process();
    //std::cout << "MCParticle " << p << " has trackID : " << trackid << " and mother : " << mother << " and status code : " << process << std::endl;
    if (process == "primary") { 
      showerpartidx = p; 
      showertrackid = mcpart.TrackId();
      nprimary += 1; 
    }
  }// for all mcparticles
  if (nprimary != 1) {
    std::cout << "ERROR : there are " << nprimary << " mcparticles. Quit " << std::endl;
    return; 
  }

  auto showerpart = mcpart_h->at(showerpartidx);

  // simulation/drift time correction
  auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();

  double g4Ticks = detClocks->TPCG4Time2Tick(showerpart.Position(0).T()) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  _mcpart_x_start = showerpart.Position(0).X();
  _mcpart_y_start = showerpart.Position(0).Y();
  _mcpart_z_start = showerpart.Position(0).Z();
  _mcpart_e_start = showerpart.Momentum(0).E() * 1000.;
  std::cout << "vtx @ [" << _mcpart_x_start << ", " << _mcpart_y_start << ", " << _mcpart_z_start << " ]" << std::endl;
  _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
  
  
  // SCE corrections
  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset_start = SCE->GetPosOffsets(geo::Point_t(_mcpart_x_start,_mcpart_y_start,_mcpart_z_start));
  std::cout << "offset_start : " << offset_start.X() << ", " << offset_start.Y() << ", " << offset_start.Z() << std::endl;
  _xsceoffset_start = offset_start.X();
  _ysceoffset_start = offset_start.Y();
  _zsceoffset_start = offset_start.Z();


  _mcpart_x_edep = -999;
  _mcpart_y_edep = -999;
  _mcpart_z_edep = -999;

  // get conversion point for photons
  for(unsigned int j = 0; j < mcpart_h->size(); j++){
    auto mcpart = mcpart_h->at(j);
    //auto mcp_pdg =  mcpart.PdgCode();
    //only want photons
    //if(mcp_pdg != 22) continue;
    
    size_t mctrkid = mcpart.TrackId();
    // std::cout<<"trackid = "<<trackid<<" and mctrkid = "<<(unsigned int)mctrkid<<std::endl;
    // if the track id for the MCParticle matches the MCShower
    if ( showertrackid == mctrkid ){
      //std::cout<<"trackid = "<<trackid<<" and mctrkid = "<<mctrkid<<", "<<(unsigned int)mctrkid<<std::endl;
      
      //get the trajectory for the MCP
      simb::MCTrajectory trj = mcpart.Trajectory();
      auto npoints = trj.size();
      
      //only considering photons with at least two points in the trajectory
      if(npoints <2) continue;
      
      //std::cout<<"the traj point at "<<0<<" is "<<trj.Position(0).X()<<", "<<trj.Position(0).Y()<<", "<<trj.Position(0).Z()<<std::endl;  
      
      //save the conversion point of the photon (point at position 1 in traj of the MCP) ->note: need to put purity cut >0 to get good match MCP/MCS
      _mcpart_x_edep = trj.Position(1).X();
      _mcpart_y_edep = trj.Position(1).Y();
      _mcpart_z_edep = trj.Position(1).Z();
      
    }//for MCP that matches MCS
  }//for each MCP in event

  // load SimChannel
  auto const& simch_h = e.getValidHandle<std::vector<sim::SimChannel> >("largeant"); //("driftWC:simpleSC");
  
  // integrate simchannel energy deposits to be able to validate mcshower.DetProfile().E()
  _mcpart_e_edep = 0;
  _mcpart_q_edep = 0;
  for (unsigned int s = 0; s < simch_h->size(); s++) {
    auto simch = simch_h->at(s);
    const auto &sch_map(simch.TDCIDEMap());
      // in MicroBooNE the collection plane starts at channel 4800
      if (simch.Channel() < 4800) continue;
      // Loop over ticks
      for(auto tdc_iter = sch_map.begin(); tdc_iter!=sch_map.end(); ++tdc_iter) {
        // for c2: hit_time is unused
	//unsigned short hit_time = (*tdc_iter).first;
	// Loop over IDEs
	for(auto const &ide : (*tdc_iter).second) {
	  _mcpart_e_edep += ide.energy;
	  _mcpart_q_edep += ide.numElectrons;
	}// for all IDEs in SimChannel
      }// for all TDCs
  }// for all simchannels  

  auto offset_edep = SCE->GetPosOffsets(geo::Point_t(_mcpart_x_edep,_mcpart_y_edep,_mcpart_z_edep));
  std::cout << "offset_edep : " << offset_edep.X() << ", " << offset_edep.Y() << ", " << offset_edep.Z() << std::endl;
  _xsceoffset_edep = offset_edep.X();
  _ysceoffset_edep = offset_edep.Y();
  _zsceoffset_edep = offset_edep.Z();
  
  _tree->Fill();

  return;
}

void MCRecoCheck::beginJob()
{
  // Implementation of optional member function here.
}

void MCRecoCheck::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MCRecoCheck)

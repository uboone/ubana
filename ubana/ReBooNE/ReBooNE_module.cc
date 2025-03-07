////////////////////////////////////////////////////////////////////////
// Class:       ReBooNE2
// Plugin Type: producer (art v3_01_02)
// File:        ReBooNE2_module.cc
//
// Generated at Fri Dec 20 09:24:25 2024 by Anyssa Navrer-Agasson using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Additional framework includes
//#include "art_root_io/TFileService.h"
//#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "nusimdata/SimulationBase/MCFlux.h"
//#include "GENIE/Framework/ParticleData/PDGUtils.h"

#include <memory>

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TObject.h"

//Additional includes
#include "ubana/ReBooNE/Utils/BooNENtuple.h"
#include "ubana/ReBooNE/Utils/BeamNtuple.h"
#include "ubana/ReBooNE/Utils/BooNEInfo.h"

namespace reboone {
  class ReBooNE;
}


class reboone::ReBooNE : public art::EDProducer {
public:
  explicit ReBooNE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ReBooNE(ReBooNE const&) = delete;
  ReBooNE(ReBooNE&&) = delete;
  ReBooNE& operator=(ReBooNE const&) = delete;
  ReBooNE& operator=(ReBooNE&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::string fGenLabel;
  std::string fLocTemplate;
  std::string fBooNETree;

  // run, event -> list of entries in file
  std::map<int, std::map<int, std::vector<long>>> fMapOfFluxTreeEntries;
  std::map<int, TTree*> fTrees;

  BooNENtuple fBooneNtp;
  BeamNtuple fBeamNtp;
  BooNEInfo fBooNEInfo;

  TFile *fCurrFile;
  
  std::string fTmploc;

  TFile *fCachedEntriesFile;

};


reboone::ReBooNE::ReBooNE(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fGenLabel(p.get<std::string>("truth_label")),
  fLocTemplate(p.get<std::string>("flux_location_template")), 
  fBooNETree(p.get<std::string>("BooNETree"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<BooNEInfo> >();
  produces< art::Assns<simb::MCTruth, BooNEInfo> >();

  // Call appropriate consumes<>() functions here.
  //Consumes the MCTruth and MCFlux data products as well as their association
  consumes< std::vector<simb::MCTruth> >(fGenLabel);
  consumes< std::vector<simb::MCFlux>  >(fGenLabel);
  consumes< art::Assns<simb::MCTruth, simb::MCFlux> >(fGenLabel);
}

void reboone::ReBooNE::produce(art::Event& e)
{
  // Implementation of required member function here.
  art::Handle<std::vector<simb::MCFlux>> fluxHandle;
  art::Handle<std::vector<simb::MCTruth>> truthHandle;

  
  std::unique_ptr<std::vector<BooNEInfo>> boonecol(new std::vector<BooNEInfo>);
  std::unique_ptr< art::Assns<simb::MCTruth, BooNEInfo> > tdkassn(new art::Assns<simb::MCTruth, BooNEInfo>);
  art::PtrMaker<BooNEInfo> makeBooNEPtr(e);

  if(e.getByLabel(fGenLabel, fluxHandle) && e.getByLabel(fGenLabel, truthHandle)) {
    art::FindManyP<simb::MCTruth> flux2truth(fluxHandle, e, fGenLabel);
    for(size_t iflux = 0; iflux < fluxHandle->size(); ++iflux) {
      auto const& f = fluxHandle->at(iflux); //f is the current MCFlux object
      const int run = f.frun;
      const int nu_pdg = f.fntype;
      const double nu_ene = f.fnenergyn;
      //const double decay_vertex_x = f.fvx;
      //const double decay_vertex_y = f.fvy;
      //const double decay_vertex_z = f.fvz;

      float windowbase_corr[3] = {620.963 , 500.563 , -1001.35};

      std::cout << "Flux job: " << run << std::endl;
      
      auto load_trees = [&](const std::string& fn) -> std::vector<TTree*> {
        
        std::vector<TTree*> boone_trees;
        fCurrFile = new TFile(fn.c_str());
        std::cout << "Opening file: " << fn << std::endl;
        TTree *tBoone = dynamic_cast<TTree*>(fCurrFile->Get(fBooNETree.c_str()));
        if(!tBoone) {
          throw cet::exception("LogicError") << "Cannot find BooNE ntuple tree " << fBooNETree << " in " << fn <<  std::endl;
        }
        tBoone->SetBranchAddress("beamwgt",&fBooneNtp.beamwgt);
        tBoone->SetBranchAddress("ntp",&fBooneNtp.ntp);
        tBoone->SetBranchAddress("npart",&fBooneNtp.npart);
        tBoone->SetBranchAddress("id",fBooneNtp.id);
        tBoone->SetBranchAddress("ini_pos",&fBooneNtp.ini_pos[0][0]);
        tBoone->SetBranchAddress("ini_mom",&fBooneNtp.ini_mom[0][0]);
        tBoone->SetBranchAddress("ini_eng",fBooneNtp.ini_eng);
        tBoone->SetBranchAddress("ini_t",fBooneNtp.ini_t);
        tBoone->SetBranchAddress("fin_mom",&fBooneNtp.fin_mom[0][0]);
        tBoone->SetBranchAddress("fin_pol",&fBooneNtp.fin_pol[0][0]);
        boone_trees.push_back(tBoone);

        TTree *tWindow = dynamic_cast<TTree*>(fCurrFile->Get("h220"));
        if(!tWindow) {
          throw cet::exception("LogicError") << "Cannot find Beam ntuple tree in "<< fn <<  std::endl;
        }
        tWindow->SetBranchAddress("pot", &fBeamNtp.pot);
        tWindow->SetBranchAddress("tank_pos_beam",&fBeamNtp.tank_pos_beam[0]);
        tWindow->SetBranchAddress("targ_pos_beam",&fBeamNtp.targ_pos_beam[0]);
        tWindow->SetBranchAddress("windowbase",&fBeamNtp.windowbase[0]);
        tWindow->SetBranchAddress("windowdir1",&fBeamNtp.windowdir1[0]);
        tWindow->SetBranchAddress("windowdir2",&fBeamNtp.windowdir2[0]);
        boone_trees.push_back(tWindow);


        return boone_trees;
      };
      
      if(fTrees.find(run) == fTrees.end()) { //if tree number not in tree map
        std::string fname = Form(fLocTemplate.c_str(), run);
        std::vector<TTree*> boone_trees = load_trees(fname);
        std::cout << "Trees loaded!" << std::endl;
        std::cout << "Entries in tBoone tree: " << boone_trees[0]->GetEntries() << std::endl;

        boone_trees[1]->GetEntry(0);
        std::cout << "POT count: " << fBeamNtp.pot << std::endl;
        for(long i = 0; i < boone_trees[0]->GetEntries(); ++i) {
            boone_trees[0]->GetEntry(i);
            fMapOfFluxTreeEntries[run][fBeamNtp.pot].push_back(i);
          }
        fTrees[run] = boone_trees[0];
      }

      std::vector<BooNEInfo> found_boones;
        found_boones.clear();
        found_boones.reserve(1);

      double targ_x = fBeamNtp.targ_pos_beam[0]/100;
      double targ_y = fBeamNtp.targ_pos_beam[1]/100;
      double targ_z = fBeamNtp.targ_pos_beam[2]/100;
      double tank_z = (fBeamNtp.tank_pos_beam[2]+windowbase_corr[2])/100;
      //const int event = f.fevtno;

      //std::cout << "Looking for event: " << event << std::endl;
      //std::cout << "At vertex: " << decay_vertex_x << " " << decay_vertex_y << " " << decay_vertex_z << std::endl;
      //std::cout << "With energy: " << nu_ene << std::endl;
      //std::cout << "Beam info: " << fBeamNtp.targ_pos_beam[2] << " " << fBeamNtp.windowbase[2]  << " " << windowbase_corr[2] << std::endl;
      //std::cout << "Target pos: " << targ_x << " " << targ_y << " " << targ_z << std::endl;

      

      for(int entry = 0; entry < fTrees[run]->GetEntries(); ++entry ){

        
        fTrees[run]->GetEntry(entry);
            
        if ( fBooneNtp.ntp == 1 ){
	        fBooNEInfo.pdg = 12; //nue
        }
        else if ( fBooneNtp.ntp == 2 ){
	        fBooNEInfo.pdg  = -12; //nuebar
        }
        else if ( fBooneNtp.ntp == 3 ){
	        fBooNEInfo.pdg = 14; //numu
        }
        else if ( fBooneNtp.ntp == 4 ){
	        fBooNEInfo.pdg = -14; //numubar
        }
        else{
	        std::cout << "Neutrino type not recognized!!! ntp = " << fBooneNtp.ntp << std::endl;
        }


        double nu_x = fBooneNtp.ini_pos[0][0]/100;
        double nu_y = fBooneNtp.ini_pos[0][1]/100;
        double nu_z = fBooneNtp.ini_pos[0][2]/100;
      
        double nu_momx = fBooneNtp.ini_mom[0][0];
        double nu_momy = fBooneNtp.ini_mom[0][1];
        double nu_momz = fBooneNtp.ini_mom[0][2];
            
        fBooNEInfo.nu_vtx_x    = nu_x + (nu_momx/nu_momz)*(tank_z-nu_z) + targ_x;
        fBooNEInfo.nu_vtx_y    = nu_y + (nu_momy/nu_momz)*(tank_z-nu_z) + targ_y;
        fBooNEInfo.nu_vtx_z    = tank_z+targ_z; 

        fBooNEInfo.nu_startt   = fBooneNtp.ini_t[0]; 
        fBooNEInfo.nu_startx   = fBooneNtp.ini_pos[0][0]/100;
        fBooNEInfo.nu_starty   = fBooneNtp.ini_pos[0][1]/100;
        fBooNEInfo.nu_startz   = fBooneNtp.ini_pos[0][2]/100;




        if(nu_pdg == fBooNEInfo.pdg && nu_ene == fBooneNtp.ini_eng[0]) 
        {
          
          double dist = sqrt( pow(nu_x-fBooNEInfo.nu_vtx_x,2) +
			                          pow(nu_y-fBooNEInfo.nu_vtx_y,2) +
			                          pow(nu_z-tank_z,2) );

          fBooNEInfo.nu_vtx_t    = fBooneNtp.ini_t[0]/1e9 + dist/TMath::C(); //s
          //std::cout << "Found a corresponding entry!" << std::endl;
          //std::cout << "i : " << entry << std::endl;
          //std::cout << "BooNE energy: " <<  fBooneNtp.ini_eng[0] << std::endl;
          //std::cout << "BooNE vertex: " << fBooNEInfo.nu_vtx_x << " " << fBooNEInfo.nu_vtx_y << " " << fBooNEInfo.nu_vtx_z << std::endl;
          //std::cout << "Nu momentum: " << nu_momx << " " << nu_momy << " " << nu_momz << std::endl;
          //std::cout << "Nu position: " << nu_x << " " << nu_y << " " << nu_z << std::endl;
          //std::cout << "Dist: " << dist << std::endl;
          //std::cout << "Time: " << fBooneNtp.ini_t[0] << std::endl;
          //std::cout << "BooNE time: " << fBooNEInfo.nu_vtx_t << std::endl;
          found_boones.push_back(fBooNEInfo);  
        }
      }
      
      
      //std::cout << "Found " << found_boones.size() << " corresponding entries!" << std::endl; 
      
      boonecol->push_back(found_boones.front());
      auto boonePtr = makeBooNEPtr(boonecol->size()-1);
      for(auto const& mctruthPtr : flux2truth.at(iflux)) {
        tdkassn->addSingle(mctruthPtr, boonePtr);
      }
    }


  }
  else {
    throw cet::exception("Configuration") << "there should be a truth generator product in the event" <<  std::endl;
  }
  e.put(std::move(boonecol));
  e.put(std::move(tdkassn));
}

DEFINE_ART_MODULE(reboone::ReBooNE)

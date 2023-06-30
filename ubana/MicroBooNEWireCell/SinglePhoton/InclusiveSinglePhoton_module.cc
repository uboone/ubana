////////////////////////////////////////////////////////////////////////
// Class:       InclusiveSinglePhoton
// Plugin Type: filter (art v2_05_00)
// File:        InclusiveSinglePhoton_module.cc
//
// Generated at Fri Jun 23 10:33:44 2017 by Robert Murrells using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"


#include "TTree.h"

#include <memory>

class InclusiveSinglePhoton : public art::EDFilter {

  TTree * ftree;

  int frun;
  int fsubrun;
  int fevent;


  Int_t           f_truth_nuPdg;
  Bool_t          f_truth_isNC;

  Float_t	f_truth_corr_showervtxX; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
  Float_t	f_truth_corr_showervtxY;
  Float_t	f_truth_corr_showervtxZ;
  Float_t f_truth_showerKE;
  Float_t	f_truth_showerMomentum[4];
  Int_t   f_truth_showerPdg;
  Int_t   f_truth_showerMother;
  std::string   f_truth_showerProcess;
  Float_t	f_truth_corr_muonvtxX; // primary muon, vtx = nu vtx
  Float_t	f_truth_corr_muonvtxY;
  Float_t	f_truth_corr_muonvtxZ;
  Float_t	f_truth_muonvtxX;
  Float_t	f_truth_muonvtxY;
  Float_t	f_truth_muonvtxZ;
  Float_t	f_truth_muonendX; // may not in TPC active
  Float_t	f_truth_muonendY;
  Float_t	f_truth_muonendZ;
  Float_t	f_truth_muonMomentum[4];
  Int_t   f_truth_Npi0;
  Int_t 	f_truth_NCDelta; // Radiative Delta label (both CC and NC)
  Int_t   f_truth_single_photon;
  Float_t f_truth_photon_angle;
  Float_t f_truth_photon_dis;

public:
  explicit InclusiveSinglePhoton(fhicl::ParameterSet const & p);

  InclusiveSinglePhoton(InclusiveSinglePhoton const &) = delete;
  InclusiveSinglePhoton(InclusiveSinglePhoton &&) = delete;
  InclusiveSinglePhoton & operator = (InclusiveSinglePhoton const &) = delete;
  InclusiveSinglePhoton & operator = (InclusiveSinglePhoton &&) = delete;

  void cout_stuff(art::Event & e, bool passed);
  void FillTree(art::Event & e);
  void Reset();
  bool filter(art::Event & e) override;

};


InclusiveSinglePhoton::InclusiveSinglePhoton(fhicl::ParameterSet const & p) :
  EDFilter(p),
  ftree(nullptr) {

  if(true) {

    art::ServiceHandle<art::TFileService> tfs;
    ftree = tfs->make<TTree>("InclusiveSinglePhotonFilterTree", "");

    ftree->Branch("run", &frun, "run/I");
    ftree->Branch("subrun", &fsubrun, "subrun/I");
    ftree->Branch("event", &fevent, "event/I");
    ftree->Branch("truth_nuPdg", 		&f_truth_nuPdg);
    ftree->Branch("truth_isNC", 		&f_truth_isNC);
    ftree->Branch("truth_corr_showervtxX", 	&f_truth_corr_showervtxX);
    ftree->Branch("truth_corr_showervtxY", 	&f_truth_corr_showervtxY);
    ftree->Branch("truth_corr_showervtxZ", 	&f_truth_corr_showervtxZ);
    ftree->Branch("truth_showerKE", 		&f_truth_showerKE);
    ftree->Branch("truth_showerMomentum", 	&f_truth_showerMomentum, "truth_showerMomentum[4]/F");
    ftree->Branch("truth_showerPdg", 		&f_truth_showerPdg);
    ftree->Branch("truth_showerMother", 		&f_truth_showerMother);
    ftree->Branch("truth_showerProcess", 		&f_truth_showerProcess);
    ftree->Branch("truth_corr_muonvtxX",	&f_truth_corr_muonvtxX);
    ftree->Branch("truth_corr_muonvtxY", 	&f_truth_corr_muonvtxY);
    ftree->Branch("truth_corr_muonvtxZ", 	&f_truth_corr_muonvtxZ);
    ftree->Branch("truth_muonvtxX",		&f_truth_muonvtxX);
    ftree->Branch("truth_muonvtxY", 		&f_truth_muonvtxY);
    ftree->Branch("truth_muonvtxZ", 		&f_truth_muonvtxZ);
    ftree->Branch("truth_muonendX",		&f_truth_muonendX);
    ftree->Branch("truth_muonendY", 		&f_truth_muonendY);
    ftree->Branch("truth_muonendZ", 		&f_truth_muonendZ);
    ftree->Branch("truth_muonMomentum", 	&f_truth_muonMomentum, "truth_muonMomentum[4]/F");
    ftree->Branch("truth_Npi0",		&f_truth_Npi0);
    ftree->Branch("truth_NCDelta",            	&f_truth_NCDelta);
    ftree->Branch("truth_single_photon",        &f_truth_single_photon);
    ftree->Branch("truth_photon_angle", &f_truth_photon_angle);
    ftree->Branch("truth_photon_dis", &f_truth_photon_dis);

  }

}


void InclusiveSinglePhoton::cout_stuff(art::Event & e, bool passed = false) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  std::cout << passed << "\n"
	    << "==========================\n";
  for(simb::MCTruth const & mct : *ev_mct) {
    std::cout << "----------------------------\n";
    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      std::cout <<"FULL: "<< mcp.TrackId() << " " << mcp.PdgCode() << " " << mcp.Mother() << " " << mcp.StatusCode() << "\n";
    }
  }

}


void InclusiveSinglePhoton::Reset() {

  frun = -1;
  fsubrun = -1;
  fevent = -1;

  f_truth_nuPdg = -1;
  f_truth_isNC = -1;

	f_truth_corr_showervtxX = -1; // truth -(SCE)-> SED -(nu time offset)-> reco [trigger offset O(10) ns ignored]
	f_truth_corr_showervtxY = -1;
	f_truth_corr_showervtxZ = -1;
	f_truth_showerKE = -1;
	f_truth_showerMomentum[0] = -1;
	f_truth_showerMomentum[1] = -1;
	f_truth_showerMomentum[2] = -1;
	f_truth_showerMomentum[3] = -1;
  f_truth_showerPdg = -1;
  f_truth_showerMother = -1.;
  f_truth_showerProcess = "";
	f_truth_corr_muonvtxX = -1; //
	f_truth_corr_muonvtxY = -1;
	f_truth_corr_muonvtxZ = -1;
	f_truth_muonvtxX = -1; //
	f_truth_muonvtxY = -1;
	f_truth_muonvtxZ = -1;
	f_truth_muonendX = -1; //
	f_truth_muonendY = -1;
	f_truth_muonendZ = -1;
	f_truth_muonMomentum[0] = -1;
	f_truth_muonMomentum[1] = -1;
	f_truth_muonMomentum[2] = -1;
	f_truth_muonMomentum[3] = -1;
  f_truth_Npi0 = 0;
  f_truth_NCDelta = -1;
  f_truth_single_photon = -1;
  f_truth_photon_angle = -1;
  f_truth_photon_dis = -1;

}


void InclusiveSinglePhoton::FillTree(art::Event & e) {

  frun = e.id().run();
  fsubrun = e.id().subRun();
  fevent = e.id().event();

  /*art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  simb::MCNeutrino const & mcn = ev_mct->at(mct_index).GetNeutrino();

  fnu_pdg = mcn.Nu().PdgCode();
  fccnc = mcn.CCNC();
  fmode = mcn.Mode();
  finteraction_type = mcn.InteractionType();

  fis_nc_delta_radiative = is_nc_delta_radiative;
  if(parent_index != SIZE_MAX) {
    fparent_status_code = ev_mct->at(mct_index).GetParticle(parent_index).StatusCode();
    fparent_pdg = ev_mct->at(mct_index).GetParticle(parent_index).PdgCode();
  }*/

  ftree->Fill();

}


bool InclusiveSinglePhoton::filter(art::Event & e) {

  Reset();

  //art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  cout_stuff(e,true);

  // Generator Info
   art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
   e.getByLabel("generator",mctruthListHandle);
   std::vector<art::Ptr<simb::MCTruth> > mclist;
   art::fill_ptr_vector(mclist, mctruthListHandle);
   art::Ptr<simb::MCTruth> mctruth;
   if (mclist.size()>0) {
       mctruth = mclist.at(0);
       if (mctruth->NeutrinoSet()) {
           simb::MCNeutrino nu = mctruth->GetNeutrino();
           f_truth_nuPdg = nu.Nu().PdgCode();
           f_truth_isNC = nu.CCNC();
        }
    }

    bool isnue = false;
    if (f_truth_isNC!=1 && abs(f_truth_nuPdg)==12 ){
      isnue = true;
      std::cout<<"InclusiveSinglePhotonFilter: Failed, NueCC"<<std::endl;
      return false;
    }


	art::Handle< std::vector<simb::MCParticle> > particleHandle2;
	if (! e.getByLabel("largeant", particleHandle2)) {
    std::cout<<"InclusiveSinglePhotonFilter: Failed, no largeant"<<std::endl;
    return false;
  }
    	std::vector< art::Ptr<simb::MCParticle> > particles2;
    	art::fill_ptr_vector(particles2, particleHandle2);
	// Get space charge correction
	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

        /**
         *  Construct track-id to MCParticle index mapping
         */

         // load MCParticles
         auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");

         // map mcp track-id to index
         std::map<int,int> MCPmap;
         // map mcp track-id to daughter track-ids
         for (size_t p=0; p < mcp_h->size(); p++) {
           auto mcp = mcp_h->at(p);
           MCPmap[mcp.TrackId()] = p;
         }


         int true_photons = 0;
         bool truth_dalitz = false;

         f_truth_photon_angle = -1;
         f_truth_photon_dis = -1;
         TVector3 dir_1;
         TVector3 start_1;
         int mother_2 = -1;

	for (auto const& particle: particles2){
//for(size_t i = 0; i < ev_mct->size(); ++i) {

  //simb::MCTruth const & mct = ev_mct->at(0);

//  for(int i = 0; i < mct.NParticles(); ++i) {
//    simb::MCParticle const & particle = mct.GetParticle(i);

		if( particle->Mother()==0 && (particle->PdgCode() == 13 || particle->PdgCode() == -13) ){
			const TLorentzVector& position = particle->Position(0);
			f_truth_muonvtxX = position.X();
			f_truth_muonvtxY = position.Y();
			f_truth_muonvtxZ = position.Z();
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_muonvtxX = position.X() - sce_offset.X();
			f_truth_corr_muonvtxY = position.Y() + sce_offset.Y();
			f_truth_corr_muonvtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_muonvtxX = (f_truth_corr_muonvtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us

			const TLorentzVector& endposition = particle->EndPosition();
			f_truth_muonendX = endposition.X();
			f_truth_muonendY = endposition.Y();
			f_truth_muonendZ = endposition.Z();

			const TLorentzVector& momentum = particle->Momentum(0);
			f_truth_muonMomentum[0] = momentum.Px();
			f_truth_muonMomentum[1] = momentum.Py();
			f_truth_muonMomentum[2] = momentum.Pz();
			f_truth_muonMomentum[3] = momentum.E();
		}

    //if (particle->PdgCode() == 111){f_truth_Npi0++;}

    if( particle->PdgCode() == 22 && particle->Process()!="eBrem" && particle->Process()!="annihil" && particle->StatusCode() == 1 &&
        !(f_truth_isNC!=1 && f_truth_nuPdg==12)){
      f_truth_showerPdg = particle->PdgCode();
			const TLorentzVector& position = particle->Position(0);
			//f_truth_corr_showervtxX = position.X(); // units: cm inherit from larsoft
			//f_truth_corr_showervtxY = position.Y(); // units: cm inherit from larsoft
			//f_truth_corr_showervtxZ = position.Z(); // units: cm inherit from larsoft
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_showervtxX = position.X() - sce_offset.X();
			f_truth_corr_showervtxY = position.Y() + sce_offset.Y();
			f_truth_corr_showervtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_showervtxX = (f_truth_corr_showervtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us
			const TLorentzVector& showerMom = particle->Momentum(0);
			float f_truth_showerKE_temp = showerMom.E() - showerMom.M();


      float ex = particle->EndX();
      float ey = particle->EndY();
      float ez = particle->EndZ();
      auto end_sce_offset = SCE->GetPosOffsets(geo::Point_t(ex, ey, ez));
      ex -= end_sce_offset.X();
      ey += end_sce_offset.Y();
      ez += end_sce_offset.Z();
      if (f_truth_showerKE_temp  > 0.02 &&
          ex > 3.0 && ex < 253.0 && ey > -113.0 && ey < 114.0 && ez > 3.0 && ez < 1034.0){
            true_photons+=1;
            f_truth_showerKE = f_truth_showerKE_temp;
            f_truth_showerMomentum[0] = showerMom.Px();
      			f_truth_showerMomentum[1] = showerMom.Py();
      			f_truth_showerMomentum[2] = showerMom.Pz();
      			f_truth_showerMomentum[3] = showerMom.E();
            f_truth_showerProcess = particle->Process();
            if (true_photons==1){
              TVector3 dir(showerMom.Px(),showerMom.Py(),showerMom.Pz());
              TVector3 dis(ex,ey,ez);
              dir_1 = dir;
              start_1 = dis;
              f_truth_photon_dis = start_1.Mag();
              mother_2 = particle->Mother();
            }

            if (true_photons==2 && particle->Mother()==mother_2){
              TVector3 dir_2(showerMom.Px(),showerMom.Py(),showerMom.Pz());
              TVector3 start_2(ex,ey,ez);
              f_truth_photon_angle = (dir_1.Angle(dir_2))*(180./TMath::Pi());
              f_truth_photon_dis = (start_1-start_2).Mag();
            }

            if (particle->Mother()!=-1){
              simb::MCParticle mother = mcp_h->at( MCPmap[particle->Mother()] );//mct.GetParticle(particle->Mother());
              if (mother.PdgCode()==22){
                if (mother.Mother()!=-1.){
                  simb::MCParticle mothermother = mcp_h->at( MCPmap[mother.Mother()] );//mct.GetParticle(mother.Mother() );
                  if (mother.PdgCode()==22){
                    if (mother.Mother()!=-1.){
                      simb::MCParticle mothermothermother = mcp_h->at( MCPmap[mothermother.Mother()] );//mct.GetParticle(mother.Mother() );
                      f_truth_showerMother = mothermothermother.PdgCode();
                    }
                  }else{
                  f_truth_showerMother = mothermother.PdgCode();
                }
              }
              }else{
                f_truth_showerMother = mother.PdgCode();
              }
            }
      }

		}

    int ndaughters = particle->NumberDaughters();
      if ( particle->PdgCode() == 111 ) {
        f_truth_Npi0++;
      for (int i_d=0; i_d<ndaughters; i_d++){
        simb::MCParticle daughter = mcp_h->at( MCPmap[particle->Daughter(i_d)] );
        if (ndaughters>2 && abs(daughter.PdgCode())==11){ truth_dalitz = true;}
      }
  	}

    /*if( abs(particle->PdgCode()) == 11 && particle->Process()!="eBrem" && particle->Process()!="annihil"
      && particle->StatusCode() == 1 && particle->Momentum(0).E() > 0.02){
      truth_dalitz = true;
      return false;
    }*/

	}
//}

	if (mclist.size()>0) {
		mctruth = mclist.at(0);
		if (mctruth->NeutrinoSet()) {
			simb::MCNeutrino nu = mctruth->GetNeutrino();
			// one can access more neutrino GENIE info

			// NC Delta radiate events
			// Delta (top 2 BRs: Delta+ and Delta0)
			// Type 1: Delta --> Gamma
			// Type 2: Delta --> Gamma (internal) --> Gamma (final state)
			f_truth_NCDelta = 0;
			for(int k=0; k<mctruth->NParticles(); k++){
				simb::MCParticle const & mcp = mctruth->GetParticle(k);
				if(mcp.PdgCode()==22 && (mcp.StatusCode()==1 || mcp.StatusCode()==14)){
					const simb::MCParticle mother = mctruth->GetParticle(mcp.Mother());
          f_truth_showerMother = mother.PdgCode();
          if (mother.PdgCode()==22){
            const simb::MCParticle mothermother = mctruth->GetParticle(mother.Mother());
            f_truth_showerMother = mothermother.PdgCode();
          }
					//if(abs(mother.PdgCode()) == 2114 || abs(mother.PdgCode()) == 2214 || abs(mother.PdgCode()) == 2224 || abs(mother.PdgCode()) == 1114) {
          if(abs(f_truth_showerMother) == 2114 || abs(f_truth_showerMother) == 2214 || abs(f_truth_showerMother) == 2224 || abs(f_truth_showerMother) == 1114) {
          f_truth_NCDelta = 1;
					break;
					}
				}
			}

		}

	}




  bool isnucc = false;
  if (f_truth_isNC!=1 && abs(f_truth_nuPdg)==14 && f_truth_muonMomentum[3]- 0.105658 > 0.1){
    isnucc = true;
    std::cout<<"InclusiveSinglePhotonFilter: Failed, NumuCC"<<std::endl;
    return false;
  }

  //if(f_truth_Npi0>0){FillTree(e);}

  if (true_photons==1 && !truth_dalitz && !isnucc && !isnue){
    f_truth_single_photon = 1;
    if(ftree) FillTree(e);
    std::cout<<"InclusiveSinglePhotonFilter: Passed, one photon"<<std::endl;
    return true;
  }else if (true_photons==2 && !truth_dalitz && !isnucc && f_truth_photon_angle>-1 && f_truth_photon_angle<20 && !isnue){
    f_truth_single_photon = 1;
    if(ftree) FillTree(e);
    std::cout<<"InclusiveSinglePhotonFilter: Passed, overlapping photons"<<std::endl;
    return true;
  }

	//

	/// truth end


  std::cout<<"InclusiveSinglePhotonFilter: Failed, other"<<std::endl;
  return false;

}


DEFINE_ART_MODULE(InclusiveSinglePhoton)

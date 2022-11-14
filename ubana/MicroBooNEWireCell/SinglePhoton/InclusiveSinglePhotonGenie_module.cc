////////////////////////////////////////////////////////////////////////
// Class:       InclusiveSinglePhotonGenie
// Plugin Type: filter (art v2_05_00)
// File:        InclusiveSinglePhotonGenie_module.cc
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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"


#include "TTree.h"

#include <memory>

class InclusiveSinglePhotonGenie : public art::EDFilter {

  TTree * ftree;

  int frun;
  int fsubrun;
  int fevent;


  Int_t           f_truth_nuPdg;
  Bool_t          f_truth_isNC;

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

public:
  explicit InclusiveSinglePhotonGenie(fhicl::ParameterSet const & p);

  InclusiveSinglePhotonGenie(InclusiveSinglePhotonGenie const &) = delete;
  InclusiveSinglePhotonGenie(InclusiveSinglePhotonGenie &&) = delete;
  InclusiveSinglePhotonGenie & operator = (InclusiveSinglePhotonGenie const &) = delete;
  InclusiveSinglePhotonGenie & operator = (InclusiveSinglePhotonGenie &&) = delete;

  void cout_stuff(art::Event & e, bool passed);
  void FillTree(art::Event & e);
  void Reset();
  bool filter(art::Event & e) override;

};


InclusiveSinglePhotonGenie::InclusiveSinglePhotonGenie(fhicl::ParameterSet const & p) :
  ftree(nullptr) {

  if(true) {

    art::ServiceHandle<art::TFileService> tfs;
    ftree = tfs->make<TTree>("InclusiveSinglePhotonGenieFilterTree", "");

    ftree->Branch("run", &frun, "run/I");
    ftree->Branch("subrun", &fsubrun, "subrun/I");
    ftree->Branch("event", &fevent, "event/I");
    ftree->Branch("truth_nuPdg", 		&f_truth_nuPdg);
    ftree->Branch("truth_isNC", 		&f_truth_isNC);
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

  }

}


void InclusiveSinglePhotonGenie::cout_stuff(art::Event & e, bool passed = false) {

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


void InclusiveSinglePhotonGenie::Reset() {

  frun = -1;
  fsubrun = -1;
  fevent = -1;

  f_truth_nuPdg = -1;
  f_truth_isNC = -1;

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

}


void InclusiveSinglePhotonGenie::FillTree(art::Event & e) {

  frun = e.id().run();
  fsubrun = e.id().subRun();
  fevent = e.id().event();

  ftree->Fill();

}


bool InclusiveSinglePhotonGenie::filter(art::Event & e) {

  Reset();

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  //cout_stuff(e,true);

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

    if (f_truth_isNC!=1 && f_truth_nuPdg==12 ){
      std::cout<<"InclusiveSinglePhotonGenieFilter: Failed, NueCC"<<std::endl;
      return false;
    }


	// Get space charge correction
	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();




for(size_t i = 0; i < ev_mct->size(); ++i) {

  simb::MCTruth const & mct = ev_mct->at(0);

  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & particle = mct.GetParticle(i);

		if( particle.Mother()==0 && (particle.PdgCode() == 13 || particle.PdgCode() == -13) ){
			const TLorentzVector& position = particle.Position(0);
			f_truth_muonvtxX = position.X();
			f_truth_muonvtxY = position.Y();
			f_truth_muonvtxZ = position.Z();
			auto sce_offset = SCE->GetPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
			f_truth_corr_muonvtxX = position.X() - sce_offset.X();
			f_truth_corr_muonvtxY = position.Y() + sce_offset.Y();
			f_truth_corr_muonvtxZ = position.Z() + sce_offset.Z();
			f_truth_corr_muonvtxX = (f_truth_corr_muonvtxX + 0.6)*1.101/1.098 + position.T()*1e-3*1.101*0.1; //T: ns; 1.101 mm/us

			const TLorentzVector& endposition = particle.EndPosition();
			f_truth_muonendX = endposition.X();
			f_truth_muonendY = endposition.Y();
			f_truth_muonendZ = endposition.Z();

			const TLorentzVector& momentum = particle.Momentum(0);
			f_truth_muonMomentum[0] = momentum.Px();
			f_truth_muonMomentum[1] = momentum.Py();
			f_truth_muonMomentum[2] = momentum.Pz();
			f_truth_muonMomentum[3] = momentum.E();
		}

	}
}


  if (f_truth_isNC!=1 && abs(f_truth_nuPdg)==14 && f_truth_muonMomentum[3]- 0.105658 > 0.1){
    std::cout<<"InclusiveSinglePhotonGenieGenieFilter: Failed, NumuCC"<<std::endl;
    return false;
  }



	/// truth end


  std::cout<<"InclusiveSinglePhotonGenieGenieFilter: Passed"<<std::endl;
  if(ftree) FillTree(e);
  return true;

}


DEFINE_ART_MODULE(InclusiveSinglePhotonGenie)

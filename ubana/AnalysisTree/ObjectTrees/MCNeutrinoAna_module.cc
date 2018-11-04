////////////////////////////////////////////////////////////////////////
// Class:       MCNeutrinoAna
// Plugin Type: analyzer (art v2_11_03)
// File:        MCNeutrinoAna_module.cc
//
// Generated at Mon Oct 22 18:29:07 2018 by Wesley Ketchum using cetskelgen
// from cetlib version v3_03_01.
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

//include for the TFileService/ROOT
#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

//include the truth objects
#include "nusimdata/SimulationBase/MCTruth.h"

namespace ana {
  class MCNeutrinoAna;
}


class ana::MCNeutrinoAna : public art::EDAnalyzer {
public:
  explicit MCNeutrinoAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCNeutrinoAna(MCNeutrinoAna const &) = delete;
  MCNeutrinoAna(MCNeutrinoAna &&) = delete;
  MCNeutrinoAna & operator = (MCNeutrinoAna const &) = delete;
  MCNeutrinoAna & operator = (MCNeutrinoAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fInputTag;

  TTree* fMCTruthTree;
  TTree* fParticleTree;

  unsigned int fRun;
  unsigned int fSubrun;
  unsigned int fEvent;
  
  //truth info
  unsigned int fTruthIndex;
  int fNParticles;
  int fOrigin;

  //neutrino info
  int fMode;
  int fInteractionType;
  int fCCNC;
  int fTarget;
  int fHitNucl;
  int fHitQuark;
  float fHadronicMass;
  float fInteractionX;
  float fInteractionY;
  float fInteractionQ2;
  float fInteractionPt;
  float fInteractionTheta;

  //particle info
  unsigned int fParticleIndex;
  int fStatus;
  int fTrackId;
  int fPdgCode;
  int fMother;
  std::string fProcess;
  std::string fEndProcess;
  float fWeight;
  int fRescatter;

  float fStart_x;
  float fStart_y;
  float fStart_z;
  float fStart_t;
  float fEnd_x;
  float fEnd_y;
  float fEnd_z;
  float fEnd_t;
  
  float fPx;
  float fPy;
  float fPz;
  float fE;
  float fP;
  float fPt;
  float fMass;
  float fEnd_Px;
  float fEnd_Py;
  float fEnd_Pz;
  float fEnd_E;

  void SetupTruthTree();
  void SetupParticleTree();

  void FillEventInfo(art::Event const&);
  void FillMCTruthInfo(simb::MCTruth const&);
  void FillParticleInfo(simb::MCParticle const&);

};

void ana::MCNeutrinoAna::SetupTruthTree()
{
  fMCTruthTree->Branch("run",&fRun,"run/i");
  fMCTruthTree->Branch("subrun",&fSubrun,"subrun/i");
  fMCTruthTree->Branch("event",&fEvent,"event/i");
  fMCTruthTree->Branch("truth_index",&fTruthIndex,"truth_index/i");
  fMCTruthTree->Branch("n_particles",&fNParticles,"n_particles/I");
  fMCTruthTree->Branch("origin",&fOrigin,"origin/I");
  fMCTruthTree->Branch("mode",&fMode,"mode/I");
  fMCTruthTree->Branch("interaction_type",&fInteractionType,"interaction_type/I");
  fMCTruthTree->Branch("ccnc",&fCCNC,"ccnc/I");
  fMCTruthTree->Branch("target",&fTarget,"target/I");
  fMCTruthTree->Branch("hit_nucl",&fHitNucl,"hit_nucl/I");
  fMCTruthTree->Branch("hit_quark",&fHitQuark,"hit_quark/I");
  fMCTruthTree->Branch("hadronic_mass",&fHadronicMass,"hadronic_mass/F");
  fMCTruthTree->Branch("interaction_x",&fInteractionX,"interaction_x/F");
  fMCTruthTree->Branch("interaction_y",&fInteractionY,"interaction_y/F");
  fMCTruthTree->Branch("interaction_q2",&fInteractionQ2,"interaction_q2/F");
  fMCTruthTree->Branch("interaction_pt",&fInteractionPt,"interaction_pt/F");
  fMCTruthTree->Branch("interaction_theta",&fInteractionTheta,"interaction_theta/F");
}

void ana::MCNeutrinoAna::SetupParticleTree()
{
  fParticleTree->Branch("run",&fRun,"run/i");
  fParticleTree->Branch("subrun",&fSubrun,"subrun/i");
  fParticleTree->Branch("event",&fEvent,"event/i");
  fParticleTree->Branch("truth_index",&fTruthIndex,"truth_index/i");
  fParticleTree->Branch("p_index",&fParticleIndex,"p_index/i");
  fParticleTree->Branch("status",&fStatus,"status/I");
  fParticleTree->Branch("trackid",&fTrackId,"trackid/I");
  fParticleTree->Branch("pdgcode",&fPdgCode,"pdgcode/I");
  fParticleTree->Branch("mother",&fMother,"mother/I");
  fParticleTree->Branch("process",&fProcess);
  fParticleTree->Branch("endprocess",&fEndProcess);
  fParticleTree->Branch("weight",&fWeight,"weight/F");
  fParticleTree->Branch("rescatter",&fRescatter,"rescatter/I");
  
  fParticleTree->Branch("start_x",&fStart_x,"start_x/F");
  fParticleTree->Branch("start_y",&fStart_y,"start_y/F");
  fParticleTree->Branch("start_z",&fStart_z,"start_z/F");
  fParticleTree->Branch("start_t",&fStart_t,"start_t/F");
  
  fParticleTree->Branch("end_x",&fEnd_x,"end_x/F");
  fParticleTree->Branch("end_y",&fEnd_y,"end_y/F");
  fParticleTree->Branch("end_z",&fEnd_z,"end_z/F");
  fParticleTree->Branch("end_t",&fEnd_t,"end_t/F");
  
  fParticleTree->Branch("px",&fPx,"px/F");
  fParticleTree->Branch("py",&fPy,"py/F");
  fParticleTree->Branch("pz",&fPz,"pz/F");
  fParticleTree->Branch("p",&fP,"p/F");
  fParticleTree->Branch("pt",&fPt,"pt/F");
  fParticleTree->Branch("e",&fE,"e/F");
  fParticleTree->Branch("mass",&fMass,"mass/F");
  
  fParticleTree->Branch("end_px",&fEnd_Px,"end_px/F");
  fParticleTree->Branch("end_py",&fEnd_Py,"end_py/F");
  fParticleTree->Branch("end_pz",&fEnd_Pz,"end_pz/F");
  fParticleTree->Branch("end_e",&fEnd_E,"end_e/F");
}

void ana::MCNeutrinoAna::FillEventInfo(art::Event const& e)
{
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.event();
}

void ana::MCNeutrinoAna::FillMCTruthInfo(simb::MCTruth const& truth)
{
  fNParticles = truth.NParticles();
  fOrigin = truth.Origin();

  auto const& nu = truth.GetNeutrino();
  fMode = nu.Mode();
  fInteractionType = nu.InteractionType();
  fCCNC = nu.CCNC();
  fTarget = nu.Target();
  fHitNucl = nu.HitNuc();
  fHitQuark = nu.HitQuark();
  fHadronicMass = nu.W();
  fInteractionX = nu.X();
  fInteractionY = nu.Y();
  fInteractionQ2 = nu.QSqr();
  fInteractionPt = nu.Pt();
  fInteractionTheta = nu.Theta();
}

void ana::MCNeutrinoAna::FillParticleInfo(simb::MCParticle const& part)
{
  fStatus = part.StatusCode();
  fTrackId = part.TrackId();
  fPdgCode = part.PdgCode();
  fMother = part.Mother();
    
  fProcess = part.Process();
  fEndProcess = part.EndProcess();
    
  fWeight = part.Weight();
  fRescatter = part.Rescatter();
    
  fStart_x = part.Vx();
  fStart_y = part.Vy();
  fStart_z = part.Vz();
  fStart_t = part.T();
    
  fEnd_x = part.EndX();
  fEnd_y = part.EndY();
  fEnd_z = part.EndZ();
  fEnd_t = part.EndT();
  
  fPx = part.Px();
  fPy = part.Py();
  fPz = part.Pz();
  fP = part.P();
  fPt = part.Pt();
  fE = part.E();
  fMass = part.Mass();
    
  fEnd_Px = part.EndPx();
  fEnd_Py = part.EndPy();
  fEnd_Pz = part.EndPz();
  fEnd_E = part.EndE();
}

ana::MCNeutrinoAna::MCNeutrinoAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
{
  fInputTag = p.get<art::InputTag>("InputTag");
}

void ana::MCNeutrinoAna::analyze(art::Event const & e)
{
  //first fill the event-level information
  FillEventInfo(e);

  //grab the MCTruth object
  auto const& mctruth_handle = e.getValidHandle< std::vector<simb::MCTruth> >(fInputTag);
  auto const& mctruth_vec(*mctruth_handle);

  //loop over mctruths
  for(size_t i_t=0; i_t<mctruth_vec.size(); ++i_t){
    fTruthIndex = i_t;
    FillMCTruthInfo(mctruth_vec[i_t]);
    fMCTruthTree->Fill();

    //loop over particles
    for(size_t i_p=0; (int)i_p<fNParticles; ++i_p){
      fParticleIndex = i_p;
      FillParticleInfo(mctruth_vec[i_t].GetParticle(i_p));
      fParticleTree->Fill();
    }//end loop over particles

  }//end loop over mctruth

}//end analyze event

void ana::MCNeutrinoAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fMCTruthTree = tfs->make<TTree>("mctruth_tree","MCTruth Tree");
  SetupTruthTree();

  fParticleTree = tfs->make<TTree>("particle_tree","MCParticle Tree");
  SetupParticleTree();
}

DEFINE_ART_MODULE(ana::MCNeutrinoAna)

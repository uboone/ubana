////////////////////////////////////////////////////////////////////////
// Class:       SingleProton
// Plugin Type: analyzer (art v3_01_02)
// File:        SingleProton_module.cc
//
// Generated at Tue Jun  4 14:26:57 2019 by Lu Ren using cetskelgen
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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "ubana/NCElastic/Algorithms/BackTrackerTruthMatch.h"
#include "ubana/NCElastic/Algorithms/TrackFeatures.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "TTree.h"

class SingleProton;


class SingleProton : public art::EDAnalyzer {
public:
  explicit SingleProton(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SingleProton(SingleProton const&) = delete;
  SingleProton(SingleProton&&) = delete;
  SingleProton& operator=(SingleProton const&) = delete;
  SingleProton& operator=(SingleProton&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endSubRun(art::SubRun const& sr) override;
  void ClearLocalData();
private:

  // Declare member data here.
  std::string m_pfp_producer;
  std::string m_hit_producer;
  std::string m_geant_producer;
  std::string m_hit_mcp_producer;
  std::string m_genie_producer;
  std::string _hitassoclabel;
  std::string _tpcobject_producer;
  std::string m_potsum_producer;
  std::string m_potsum_instance;
  std::string fMCParticleModuleLabel;
  std::string fHitModuleLabel;
  std::string fMCParticleToMCTruthAssModuleLabel;
  bool _debug=true;
  bool m_isOverlay;
  bool m_isData;
  
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::ShowerVector pfshowers;
  lar_pandora::TrackVector pftracks,trackVector2;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::TracksToHits tracksToHits;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::PFParticlesToClusters particlesToClusters;
  lar_pandora::PFParticlesToShowers particlesToShowers;
  lar_pandora::PFParticlesToTracks particlesToTracks;
  lar_pandora::PFParticlesToHits pfParticlesToHits;
  lar_pandora::HitsToPFParticles hitsToPfParticles;
  lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints;
  lar_pandora::ClustersToHits clustersToHits;
  lar_pandora::HitsToSpacePoints hitsToSpacePoints;
  lar_pandora::SpacePointsToHits spacePointsToHits;
  
  TTree* _tree1;
  TTree* _sr_tree;
  int _run, _subrun, _event;
  float _sr_begintime, _sr_endtime;
  float _sr_pot;
  
  int _mc_ccnc=-9999, _mc_mode=-9999, _mc_interactiontype=-9999, _mc_pdg=-9999, _mc_primary=-9999, _mc_origin=-9999, _mc_hitnuc=-9999, _mc_nupdg=-9999;
  float _mc_q2=-9999., _mc_nu_vtxx, _mc_nu_vtxy, _mc_nu_vtxz, _mc_nu_vtxx_sce, _mc_nu_vtxy_sce, _mc_nu_vtxz_sce, _mc_enu;
  float _mc_length, _mc_start_x,_mc_start_y, _mc_start_z, _mc_end_x, _mc_end_y, _mc_end_z, _mc_theta, _mc_phi, _mc_ke, _mc_mom;

  int _n_pfp=-9999;
  float _track_score=-1;
  float _dislen_ratio=-1;
 float _KE_len=-1;
  float _reco_q2=-9999., _reco_nu_vtxx, _reco_nu_vtxy, _reco_nu_vtxz, _reco_length=-9999., _reco_start_x=-9999., _reco_start_y=-9999., _reco_start_z=-9999., _reco_end_x=-9999., _reco_end_y=-9999., _reco_end_z=-9999., _reco_theta=-9999., _reco_phi=-9999., _reco_ke=-9999., _reco_mom=-9999.;
  int _nhits_0=-9999, _nhits_1=-9999, _nhits_2=-9999;
  float _chi2_p_0=-9999., _chi2_p_1=-9999., _chi2_p_2=-9999.;
 float _chi2_mu_0=-9999., _chi2_mu_1=-9999., _chi2_mu_2=-9999.;
 float _LL_0=-9999., _LL_1=-9999., _LL_2=-9999.;
 float _TM_dedx_0=-9999., _TM_dedx_1=-9999., _TM_dedx_2=-9999.;
 float _PIDA_0=-9999., _PIDA_1=-9999., _PIDA_2=-9999.;
 float _start_dedx_0=-9999., _start_dedx_1=-9999., _start_dedx_2=-9999.;
 float _ratio_dedx_0=-9999., _ratio_dedx_1=-9999., _ratio_dedx_2=-9999.;
 float _avg_dedx_0=-9999., _avg_dedx_1=-9999., _avg_dedx_2=-9999.;
 float _total_dedx_0=-9999., _total_dedx_1=-9999., _total_dedx_2=-9999.;
  float _KE_calo_0=-9999., _KE_calo_1=-9999., _KE_calo_2=-9999.;

};


SingleProton::SingleProton(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraFlashEventBuilding");
  m_hit_producer = p.get<std::string>("hit_producer", "gaushit");
  m_geant_producer = p.get<std::string>("geant_producer", "largeant");
  m_hit_mcp_producer = p.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");
  m_isOverlay    = p.get<bool>("Overlay", false);
  m_genie_producer = p.get<std::string>("genie_producer","generator");
  m_potsum_producer                = p.get<std::string>("POTSummaryProducer","generator");
  //  m_potsum_producer= p.get<std::string>("POTSummaryProducer");
  //  m_potsum_instance = p.get<std::string>("POTSummaryInstance");
  art::ServiceHandle<art::TFileService> tfs;
  _tree1 = tfs->make<TTree>("tree","");
  _tree1->Branch("run",    &_run,    "run/I");
  _tree1->Branch("subrun", &_subrun, "subrun/I");
  _tree1->Branch("event",  &_event,  "event/I");
  _tree1->Branch("mc_ccnc",&_mc_ccnc,"mc_ccnc/I");
  _tree1->Branch("mc_mode",&_mc_mode,"mc_mode/I");
  _tree1->Branch("mc_interactiontype",&_mc_interactiontype,"mc_interactiontype/I");
  _tree1->Branch("mc_hitnuc",&_mc_hitnuc,"mc_hitnuc/I");
  _tree1->Branch("mc_q2",&_mc_q2,"mc_q2/F");
  _tree1->Branch("mc_nu_vtxx",&_mc_nu_vtxx,"mc_nu_vtxx/F");
  _tree1->Branch("mc_nu_vtxy",&_mc_nu_vtxy,"mc_nu_vtxy/F");
  _tree1->Branch("mc_nu_vtxz",&_mc_nu_vtxz,"mc_nu_vtxz/F");
  _tree1->Branch("mc_nu_vtxx_sce",&_mc_nu_vtxx_sce,"mc_nu_vtxx_sce/F");
  _tree1->Branch("mc_nu_vtxy_sce",&_mc_nu_vtxy_sce,"mc_nu_vtxy_sce/F");
  _tree1->Branch("mc_nu_vtxz_sce",&_mc_nu_vtxz_sce,"mc_nu_vtxz_sce/F");
  _tree1->Branch("mc_enu",&_mc_enu,"mc_enu/F");
  _tree1->Branch("mc_pdg",&_mc_pdg,"mc_pdg/I");
  _tree1->Branch("mc_nupdg",&_mc_pdg,"mc_nupdg/I");
  _tree1->Branch("mc_primary",&_mc_primary,"mc_primary/I");
  _tree1->Branch("mc_origin",&_mc_origin,"mc_origin/I");
  _tree1->Branch("mc_length",&_mc_length,"mc_length/F");
  _tree1->Branch("mc_start_x",&_mc_start_x,"mc_start_x/F");
  _tree1->Branch("mc_start_y",&_mc_start_y,"mc_start_y/F");
  _tree1->Branch("mc_start_z",&_mc_start_z,"mc_start_z/F");
  _tree1->Branch("mc_end_x",&_mc_end_x,"mc_end_x/F");
  _tree1->Branch("mc_end_y",&_mc_end_y,"mc_end_y/F");
  _tree1->Branch("mc_end_z",&_mc_end_z,"mc_end_z/F");
  _tree1->Branch("mc_theta",&_mc_theta,"mc_theta/F");
  _tree1->Branch("mc_phi",&_mc_phi,"mc_phi/F");
  _tree1->Branch("mc_ke",&_mc_ke,"mc_ke/F");
  _tree1->Branch("mc_mom",&_mc_mom,"mc_mom/F");
  _tree1->Branch("n_pfp",&_n_pfp,"n_pfp/I");
  _tree1->Branch("trk_score",&_track_score,"trk_score/F");
  _tree1->Branch("KE_len",&_KE_len,"KE_len/F");
  _tree1->Branch("dislen_ratio",&_dislen_ratio,"dislen_ratio/F");

  _tree1->Branch("reco_q2",&_reco_q2,"reco_q2/F");
  _tree1->Branch("reco_nu_vtxx",&_reco_nu_vtxx,"reco_nu_vtxx/F");
  _tree1->Branch("reco_nu_vtxy",&_reco_nu_vtxy,"reco_nu_vtxy/F");
  _tree1->Branch("reco_nu_vtxz",&_reco_nu_vtxz,"reco_nu_vtxz/F");
  _tree1->Branch("reco_length",&_reco_length,"reco_length/F");
  _tree1->Branch("reco_start_x",&_reco_start_x,"reco_start_x/F");
  _tree1->Branch("reco_start_y",&_reco_start_y,"reco_start_y/F");
  _tree1->Branch("reco_start_z",&_reco_start_z,"reco_start_z/F");
  _tree1->Branch("reco_end_x",&_reco_end_x,"reco_end_x/F");
  _tree1->Branch("reco_end_y",&_reco_end_y,"reco_end_y/F");
  _tree1->Branch("reco_end_z",&_reco_end_z,"reco_end_z/F");
  _tree1->Branch("reco_theta",&_reco_theta,"reco_theta/F");
  _tree1->Branch("reco_phi",&_reco_phi,"reco_phi/F");
  _tree1->Branch("reco_ke",&_reco_ke,"reco_ke/F");
  _tree1->Branch("reco_mom",&_reco_mom,"reco_mom/F");
  _tree1->Branch("nhits_0",&_nhits_0,"nhits_0/I");
  _tree1->Branch("nhits_1",&_nhits_1,"nhits_1/I");
  _tree1->Branch("nhits_2",&_nhits_2,"nhits_2/I");
  _tree1->Branch("chi2_p_0",&_chi2_p_0,"chi2_p_0/F");
  _tree1->Branch("chi2_p_1",&_chi2_p_1,"chi2_p_1/F");
  _tree1->Branch("chi2_p_2",&_chi2_p_2,"chi2_p_2/F");
  _tree1->Branch("chi2_mu_0",&_chi2_mu_0,"chi2_mu_0/F");
  _tree1->Branch("chi2_mu_1",&_chi2_mu_1,"chi2_mu_1/F");
  _tree1->Branch("chi2_mu_2",&_chi2_mu_2,"chi2_mu_2/F");

  _tree1->Branch("LL_0",&_LL_0,"LL_0/F");
  _tree1->Branch("LL_1",&_LL_1,"LL_1/F");
  _tree1->Branch("LL_2",&_LL_2,"LL_2/F");
 _tree1->Branch("TM_dedx_0",&_TM_dedx_0,"TM_dedx_0/F");
  _tree1->Branch("TM_dedx_1",&_TM_dedx_1,"TM_dedx_1/F");
  _tree1->Branch("TM_dedx_2",&_TM_dedx_2,"TM_dedx_2/F");
 _tree1->Branch("PIDA_0",&_PIDA_0,"PIDA_0/F");
  _tree1->Branch("PIDA_1",&_PIDA_1,"PIDA_1/F");
  _tree1->Branch("PIDA_2",&_PIDA_2,"PIDA_2/F");
 _tree1->Branch("start_dedx_0",&_start_dedx_0,"start_dedx_0/F");
  _tree1->Branch("start_dedx_1",&_start_dedx_1,"start_dedx_1/F");
  _tree1->Branch("start_dedx_2",&_start_dedx_2,"start_dedx_2/F");
 _tree1->Branch("ratio_dedx_0",&_ratio_dedx_0,"ratio_dedx_0/F");
  _tree1->Branch("ratio_dedx_1",&_ratio_dedx_1,"ratio_dedx_1/F");
  _tree1->Branch("ratio_dedx_2",&_ratio_dedx_2,"ratio_dedx_2/F");
 _tree1->Branch("avg_dedx_0",&_avg_dedx_0,"avg_dedx_0/F");
  _tree1->Branch("avg_dedx_1",&_avg_dedx_1,"avg_dedx_1/F");
  _tree1->Branch("avg_dedx_2",&_avg_dedx_2,"avg_dedx_2/F");
 _tree1->Branch("total_dedx_0",&_total_dedx_0,"total_dedx_0/F");
  _tree1->Branch("total_dedx_1",&_total_dedx_1,"total_dedx_1/F");
  _tree1->Branch("total_dedx_2",&_total_dedx_2,"total_dedx_2/F");
 _tree1->Branch("KE_calo_0",&_KE_calo_0,"KE_calo_0/F");
  _tree1->Branch("KE_calo_1",&_KE_calo_1,"KE_calo_1/F");
  _tree1->Branch("KE_calo_2",&_KE_calo_2,"KE_calo_2/F");


  _sr_tree = tfs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                &_run,                "run/I");
  _sr_tree->Branch("subrun",             &_subrun,             "subrun/I");
  _sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/F");
  _sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/F");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/F");
}

void SingleProton::analyze(art::Event const& e)
{
  // Implementation of required member function here.
 
  

    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
    m_isData = e.isRealData();
    if(_debug) std::cout << "=======Start of Single Proton Module for  Run/Subrun/Event: "<<_run<<" / "<<_subrun<<" / "<<_event<<"======="<<std::endl;
    if(m_isOverlay)
      {
	art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
	std::vector<art::Ptr<simb::MCTruth> > mclist;
	if (e.getByLabel(m_genie_producer,mctruthListHandle))
	  {
	    art::fill_ptr_vector(mclist, mctruthListHandle);
	  }
	if(!mclist.empty())
	  {
	    if(mclist[0]->NeutrinoSet())
	      {
		art::Ptr<simb::MCTruth> mctruth = mclist[0];
		_mc_ccnc   = mctruth->GetNeutrino().CCNC();
		_mc_mode   = mctruth->GetNeutrino().Mode();
		_mc_hitnuc = mctruth->GetNeutrino().HitNuc();
		_mc_q2     = mctruth->GetNeutrino().QSqr();
		_mc_nu_vtxx    = mctruth->GetNeutrino().Nu().Vx();
		_mc_nu_vtxy    = mctruth->GetNeutrino().Nu().Vy();
		_mc_nu_vtxz    = mctruth->GetNeutrino().Nu().Vz();
		_mc_enu    = mctruth->GetNeutrino().Nu().E();
		_mc_interactiontype= mctruth->GetNeutrino().InteractionType();
		if(_debug)
		  {
		    std::cout<<"[Single Proton] CCNC: "<<_mc_ccnc<<std::endl;
		    std::cout<<"[Single Proton] Mode: "<<_mc_mode<<std::endl;
		    std::cout<<"[Single Proton] Hitnuc: "<<_mc_hitnuc<<std::endl;
		    std::cout<<"[Single Proton] Q2: "<<_mc_q2<<std::endl;
		    std::cout<<"[Single Proton] nu (x,y,z):  ("<<_mc_nu_vtxx<<", "<<_mc_nu_vtxy<<", "<<_mc_nu_vtxz<<")"<<std::endl;
		    std::cout<<"[Single Proton] Ev:  "<<_mc_enu<<std::endl;
		  }
	      }
	  }
	
	
      }
    larpandora.CollectPFParticleMetadata(e, m_pfp_producer, pfparticles, particlesToMetadata);
    larpandora.BuildPFParticleMap(pfparticles, particleMap);
    larpandora.CollectTracks(e,m_pfp_producer, trackVector2, tracksToHits);

    art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
    // const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle = evt.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraMCSMu");
    
     art::FindManyP<anab::Calorimetry> trackCaloAssn(trackHandle, e, "pandoraFMLucaliSCE");
	
   

    const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, "pandoraFMLucalipidSCE");
    if (!trackPIDAssn.isValid()){
      std::cout << "[NuCC::FillReconstructed] trackPIDAssn.isValid() == false" << std::endl;
    }
    if (pfparticles.size() == 0)
      {
	if(_debug) std::cout << "[Single Proton] No reconstructed PFParticles in event." << std::endl;
      }
    else
      { 
	larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
	if (pfneutrinos.size() != 1)
	  {
	    if(_debug) std::cout << "[Single Proton] Number of Neutrino Slices found: " << pfneutrinos.size() << std::endl;
	  }
	else // main code for everything starts here
	  { 
	    if(_debug) std::cout << "[Single Proton] There is a neutrino slice! " << pfneutrinos.size() << std::endl;
	    std::vector<art::Ptr<recob::PFParticle> > finalStatePFParticles;
	    //  std::unordered_map<size_t, art::Ptr<recob::PFParticle> > pfParticleIdMap;
	     for (unsigned int n = 0; n < pfparticles.size(); ++n)
	       {
		 const art::Ptr<recob::PFParticle> particle = pfparticles.at(n);
		 if(particle->IsPrimary()) continue;
		 const auto parentIterator = particleMap.find(particle->Parent());
		 if (parentIterator == particleMap.end())
		   throw cet::exception("WorkshopAnalyzer") << "PFParticle parent not found" << std::endl;
		 if (!parentIterator->second->IsPrimary())
		   continue;
		 const int parentPDG = std::abs(parentIterator->second->PdgCode());
		 if (abs(parentPDG) != 14)
		   continue;
		 finalStatePFParticles.push_back(particle);
	       }
	     _n_pfp=finalStatePFParticles.size();
	     if(_n_pfp!=1) 
	       {
		 std::cout<<"[Single Proton] not exact one particle in final state??? "<<finalStatePFParticles.size()<<std::endl;	 
	       }
	     else{
	       art::FindManyP<recob::Track> pfPartToTrackAssoc(finalStatePFParticles, e, m_pfp_producer);
	       const auto track = pfPartToTrackAssoc.at(0);
	       const auto &trackVtxPosition = track.front()->Vertex();
	       const auto &trackEndPosition = track.front()->End();
	       const auto &track_direction  = track.front()->StartDirection();
	       _reco_start_x=trackVtxPosition.x();
	       _reco_start_y=trackVtxPosition.y();
	       _reco_start_z=trackVtxPosition.z();
	       _reco_end_x=trackEndPosition.x();
	       _reco_end_y=trackEndPosition.y();
	       _reco_end_z=trackEndPosition.z();
	       _reco_length= track.front()->Length();
	       _reco_theta=track_direction.Theta();
	       _reco_phi=track_direction.Phi();
	       _dislen_ratio=TMath::Sqrt( TMath::Power(_reco_end_x-_reco_start_x,2) + TMath::Power(_reco_end_y-_reco_start_y,2)  + TMath::Power(_reco_end_z-_reco_start_z,2))/ _reco_length;
	       double A=17.;
	       double bp1=0.58;
	       _KE_len=A/bp1* TMath::Power(_reco_length,bp1);
	       //  _reco_ke, _reco_mom;
	       if(_debug)
		 {
		   	
		   std::cout<<"[Single Proton] reco vtx x, y, z: ("<<trackVtxPosition.x()<<", "<<trackVtxPosition.y()<<", "<<trackVtxPosition.z()<<" )"<<std::endl;
		   std::cout<<"[Single Proton] reco end x, y, z: ("<<trackEndPosition.x()<<", "<<trackEndPosition.y()<<", "<<trackEndPosition.z()<<" )"<<std::endl;
		   std::cout<<"[Single Proton] reco length: "<< _reco_length<<std::endl;
		   std::cout<<"[Single Proton] reco dis/len: "<<_dislen_ratio<<std::endl;
		   
                   std::cout<<"[Single Proton] reco K.E.: "<< _KE_len<<std::endl;
		   std::cout<<"[Single Proton] reco theta: "<< _reco_theta<<std::endl;
		   std::cout<<"[Single Proton] reco phi: "<< _reco_phi<<std::endl;
		 }
	       // pid stuffs for the track
	       std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track.front().key());
	       if (trackPID.size() != 0)
		 {
		   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
		   //  
		   //   if (planenum<0 || planenum>2) continue;
		   for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
		     {
		       anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
		       int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
			  
		       if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward)
			 {
			   if (AlgScore.fAlgName == "Chi2"&&TMath::Abs(AlgScore.fAssumedPdg) == 2212)
			     {
			       if(planenum==0)    _chi2_p_0=AlgScore.fValue;
			       if(planenum==1)    _chi2_p_1=AlgScore.fValue;
			       if(planenum==2)    _chi2_p_2=AlgScore.fValue;
			     }
			    if (AlgScore.fAlgName == "Chi2"&&TMath::Abs(AlgScore.fAssumedPdg) == 13)
			     {
			       if(planenum==0)    _chi2_mu_0=AlgScore.fValue;
			       if(planenum==1)    _chi2_mu_1=AlgScore.fValue;
			       if(planenum==2)    _chi2_mu_2=AlgScore.fValue;
			     }
			    
			 }
		       if (AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA)
			 {
			   if(planenum==0)    _PIDA_0= AlgScore.fValue;
			   if(planenum==1)    _PIDA_1= AlgScore.fValue;
			   if(planenum==2)    _PIDA_2= AlgScore.fValue;
			   
			 }
       
		     }
		 }//end of pid info
	       //Calorimetry, dedx, etc
               if(trackCaloAssn.isValid())
		{
		std::vector<art::Ptr<anab::Calorimetry> > calos = trackCaloAssn.at(track.front().key());
		_nhits_0=calos.at(0)->dEdx().size();
		_nhits_1=calos.at(1)->dEdx().size();
		_nhits_2=calos.at(2)->dEdx().size();
		_KE_calo_0=calos.at(0)->KineticEnergy();
		_KE_calo_1=calos.at(1)->KineticEnergy();
		_KE_calo_2=calos.at(2)->KineticEnergy();
		_total_dedx_0=0;
		_total_dedx_1=0;
		_total_dedx_2=0;
		_start_dedx_0=0;
		_start_dedx_1=0;
		_start_dedx_2=0;
		float _end_dedx_0=0;
		float _end_dedx_1=0;
		float _end_dedx_2=0;
		int ii0=-1;
		int jj0=_nhits_0;
		int ii1=-1;
		int jj1=_nhits_1;
		int ii2=-1;
		int jj2=_nhits_2;
		
		if(_nhits_0>0)
		  {
		   for( auto& de : calos.at(0)->dEdx())
		     {
		       _total_dedx_0 += de;
		       ii0++;
		       jj0--;
		       if(ii0<6){_start_dedx_0+=de;}
		       if(_nhits_0-jj0<6){_end_dedx_0+=de;}
		     }
		   _avg_dedx_0 = (float)_total_dedx_0 / (float)_nhits_0;
		   if(_start_dedx_0>0)
		     {
		       _ratio_dedx_0=_end_dedx_0/_start_dedx_0;
		     }
		  }//plane 0
		if(_nhits_1>0)
		  {
		    for( auto& de : calos.at(1)->dEdx())
		      {
			_total_dedx_1 += de;
			ii1++;
			jj1--;
			if(ii1<6){_start_dedx_1+=de;}
			if(_nhits_1-jj1<6){_end_dedx_1+=de;}
		      }
		    _avg_dedx_1 = (float)_total_dedx_1 / (float)_nhits_1;
		    if(_start_dedx_1>0)
		      {
			_ratio_dedx_1=_end_dedx_1/_start_dedx_1;
		      }
		  }//plane 1
		if(_nhits_2>0)
		  {
		    for( auto& de : calos.at(2)->dEdx())
		      {
			_total_dedx_2 += de;
			ii2++;
			jj2--;
			if(ii2<6){_start_dedx_2+=de;}
			if(_nhits_2-jj2<6){_end_dedx_2+=de;}
		      }
		    _avg_dedx_2 = (float)_total_dedx_2 / (float)_nhits_2;
		    if(_start_dedx_2>0)
		      {
			_ratio_dedx_2=_end_dedx_2/_start_dedx_2;
		      }
		  }//plane 1
	       if (_debug) std::cout<<"[Single Proton] calos.size() "<<calos.size()<<std::endl;
	       if (_debug) std::cout<<"[Single Proton] KE 0 "<<_KE_calo_0<<std::endl;
	       if (_debug) std::cout<<"[Single Proton] KE 1 "<<_KE_calo_1<<std::endl;
	       if (_debug) std::cout<<"[Single Proton] KE 2 "<<_KE_calo_2<<std::endl;
	    	if (_debug) std::cout<<"[Single Proton] NHits 0 "<<_nhits_0<<std::endl;
		if (_debug) std::cout<<"[Single Proton] NHits 1 "<<_nhits_1<<std::endl;
		if (_debug) std::cout<<"[Single Proton] NHits 2 "<<_nhits_2<<std::endl;



		}//end of calo-related stuffs
	       if (_debug) std::cout<<"[Single Proton] chi proton 0, 1, 2: "<<_chi2_p_0<<", "<<_chi2_p_1<<", "<<_chi2_p_2<<std::endl;
	       if (_debug) std::cout<<"[Single Proton] chi muon 0, 1, 2: "<<_chi2_mu_0<<", "<<_chi2_mu_1<<", "<<_chi2_mu_2<<std::endl;
	       //truth information for the track
	       if(m_isOverlay)
		 {
		   // larpandora.BuildPFParticleHitMaps(e, m_pfp_producer, pfParticlesToHits, hitsToPfParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
		   // lar_pandora::TracksToHits tracksToHits;
		   lar_pandora::TracksToHits::const_iterator trkIter2 = tracksToHits.find(track.front());
		   //if (tracksToHits.end() != trkIter2)
		   //  m_trackhits = trkIter2->second.size();
		   std::vector<art::Ptr<recob::Hit>> hitVec;
		   hitVec.clear();
		   hitVec= trkIter2->second;
		   if (_debug)   std::cout << "[Single Proton] DEBUGGING 0 hitvec size "<<hitVec.size()<< std::endl;
		   art::Handle<std::vector<simb::MCParticle>> mcparticle_handle;
		   std::vector<art::Ptr<simb::MCParticle> > largeant_vec;
		   if (e.getByLabel(m_geant_producer,mcparticle_handle)) { art::fill_ptr_vector(largeant_vec,mcparticle_handle); }
		   art::Handle<std::vector<recob::Hit > > hit_handle;
		   std::vector<art::Ptr<recob::Hit> > hit_vec;
		   hit_vec.clear();
		   if(e.getByLabel(m_hit_producer,hit_handle)) { art::fill_ptr_vector(hit_vec,hit_handle); }
		   	   if (_debug)   std::cout << "[Single Proton] DEBUGGING 0 hitVec size "<<hit_vec.size()<< std::endl;
		   	   std::cout<<"[Single Proton] DEBUGGING 0 "<<std::endl;
		   art::FindOneP<simb::MCTruth> MCParticleToMCTruth(mcparticle_handle,e,m_geant_producer);
		   std::cout<<"[Single Proton] DEBUGGING 1.1 "<<std::endl;
		   BackTrackerTruthMatch backtrackertruthmatch;
		   backtrackertruthmatch.MatchToMCParticle(hit_handle,e,hitVec);
		   std::cout<<"[Single Proton] DEBUGGING 1.2 "<<std::endl;
		   art::Ptr< simb::MCParticle > maxp_me = backtrackertruthmatch.ReturnMCParticle();
		   std::cout<<"[Single Proton] DEBUGGING 2 "<<std::endl;
		   if(!(maxp_me.isNull()))
		     {
		       if (_debug) std::cout<<"[Single Proton] DEBUGGING 3 "<<std::endl;
		       if (_debug) std::cout << "Track matched by the BackTracker" << std::endl;
		       const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruth.at(maxp_me.key());
		       if (_debug)   std::cout <<mctruth->Origin()<<" "<<maxp_me.key()<<std::endl;
		       //       std::cout<<mcparticle->PdgCode()<<" "<<maxp_me->PdgCode()<<std::endl;
		       _mc_pdg     = maxp_me->PdgCode();
		       _mc_primary = float(maxp_me->Process() == "primary");
		       _mc_origin  = mctruth->Origin();
		       _mc_length  = maxp_me->Trajectory().TotalLength();
		       _mc_start_x  = maxp_me->Vx();
		       _mc_start_y  = maxp_me->Vy();
		       _mc_start_z  = maxp_me->Vz();
		       _mc_end_x    = maxp_me->EndX();
		       _mc_end_y    = maxp_me->EndY();
		       _mc_end_z    = maxp_me->EndZ();
		       _mc_theta   = maxp_me->Momentum().Theta();
		       _mc_phi     = maxp_me->Momentum().Phi();
		       _mc_ke = maxp_me->E() - maxp_me->Mass();
		       
		       if(_debug)
			 {
			   std::cout<<"[Single Proton] true PDG: "<< _mc_pdg<<std::endl;
			   std::cout<<"[Single Proton] true K.E.: "<< _mc_ke<<std::endl;
			   std::cout<<"[Single Proton] true start x, y, z: ("<<_mc_start_x<<", "<<_mc_start_y<<", "<<_mc_start_z<<" )"<<std::endl;
			   std::cout<<"[Single Proton] true end x, y, z: ("<< _mc_end_x<<", "<< _mc_end_y<<", "<< _mc_end_z<<" )"<<std::endl;
			   std::cout<<"[Single Proton] true length: "<< _mc_length<<std::endl;
			   std::cout<<"[Single Proton] true theta: "<< _mc_theta<<std::endl;
			   std::cout<<"[Single Proton] true phi: "<< _mc_phi<<std::endl;
			 }
		     }
		   else
		     {
		       if (_debug) std::cout<<"[Single Proton] DEBUGGING 4 "<<std::endl;
		       if (_debug) std::cout << "Track not matched by the BackTracker- Cosmic!" << std::endl;
		       _mc_pdg     = -9999;
		       _mc_primary = -9999;
		       _mc_origin  = 2;
		       _mc_length  = -9999;
		       _mc_start_x  = -9999;
		       _mc_start_y  = -9999;
		       _mc_start_z  =-9999;
		       _mc_end_x    = -9999;
		       _mc_end_y    =-9999;
		       _mc_end_z    = -9999;
		       _mc_theta   =-9999;
		       _mc_phi     = -9999;
			 _mc_ke =-9999;
		     }
		   hitVec.clear();
		   largeant_vec.clear();
		   //  maxp_me.clear();
		 }
	     }	     
	  }//end of main code for everything
      }
    _tree1->Fill();
    ClearLocalData();
}

void SingleProton::beginJob()
{
  // Implementation of optional member function here.
}

void SingleProton::endSubRun(art::SubRun const& sr)
{
  // Implementation of optional member function here.
  if (_debug) std::cout << "[UBXSec::endSubRun] Starts" << std::endl;

  _run       = sr.run();
  _subrun    = sr.subRun();
 
  art::Handle<sumdata::POTSummary> potsum_h;

  // MC
  if (m_isOverlay) {
    if (_debug) std::cout << "[NCE::endSubRun] Getting POT for MC" << std::endl;
    if(sr.getByLabel(m_potsum_producer, potsum_h)) {
      if (_debug) std::cout << "[NCE::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot/1.0e16;
//std::cout <<"MC pot:  "<<_sr_pot << std::endl;
    }
    else
      _sr_pot = 0.;
 std::cout <<"MC pot:  "<<_sr_pot << std::endl;
 }
_sr_tree->Fill();
}
void SingleProton::ClearLocalData()
{
  pfparticles.clear();
  pfneutrinos.clear();
  pfdaughters.clear();
  pfshowers.clear();
  pftracks.clear();
  particleMap.clear();
  particlesToMetadata.clear();
  particlesToVertices.clear();
  particlesToClusters.clear();
  particlesToSpacePoints.clear();
  particlesToShowers.clear();
  particlesToTracks.clear();
  clustersToHits.clear();
  hitsToSpacePoints.clear();
  spacePointsToHits.clear();
  tracksToHits.clear();
 _mc_ccnc=-9999; 
 _mc_mode=-9999; 
 _mc_interactiontype=-9999; 
 _mc_pdg=-9999; _mc_primary=-9999; _mc_origin=-9999; _mc_hitnuc=-9999; _mc_nupdg=-9999;
 _mc_q2=-9999.; 
 /*_mc_nu_vtxx; _mc_nu_vtxy; _mc_nu_vtxz; _mc_nu_vtxx_sce; _mc_nu_vtxy_sce; _mc_nu_vtxz_sce; _mc_enu;  _track_score; _reco_nu_vtxx; _reco_nu_vtxy; _reco_nu_vtxz;  _nhits_0; _nhits_1; _nhits_2;
 _mc_length; _mc_start_x;_mc_start_y; _mc_start_z; _mc_end_x; _mc_end_y; _mc_end_z; _mc_theta; _mc_phi; _mc_ke; _mc_mom;
 */
 _n_pfp=-9999;
 _reco_q2=-9999.; _reco_length=-9999.; _reco_start_x=-9999.; _reco_start_y=-9999.; _reco_start_z=-9999.; _reco_end_x=-9999.; _reco_end_y=-9999.; _reco_end_z=-9999.; _reco_theta=-9999.; _reco_phi=-9999.; _reco_ke=-9999.; _reco_mom=-9999.;
 _chi2_p_0=-9999.; _chi2_p_1=-9999.; _chi2_p_2=-9999.;
 _chi2_mu_0=-9999.; _chi2_mu_1=-9999.; _chi2_mu_2=-9999.;
}

DEFINE_ART_MODULE(SingleProton)

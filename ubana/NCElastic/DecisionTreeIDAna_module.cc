////////////////////////////////////////////////////////////////////////
// Class:       DecisionTreeIDAna
// Plugin Type: analyzer (art v2_05_00)
// File:        DecisionTreeIDAna_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "TString.h"
#include "TTree.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

//#include "ubana/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
//#include "ubana/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "ubana/NCElastic/Algorithms/BackTrackerTruthMatch.h"
#include "ubana/NCElastic/Algorithms/TrackFeatures.h"

#include <memory>
#include <iostream>

class DecisionTreeIDAna;

class DecisionTreeIDAna : public art::EDAnalyzer {
public:
  explicit DecisionTreeIDAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DecisionTreeIDAna(DecisionTreeIDAna const &) = delete;
  DecisionTreeIDAna(DecisionTreeIDAna &&) = delete;
  DecisionTreeIDAna & operator = (DecisionTreeIDAna const &) = delete;
  DecisionTreeIDAna & operator = (DecisionTreeIDAna &&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const & e) override;
  void endSubRun(const art::SubRun &sr);
private:

  // Pointer to the feature algorithm
  ::dtfeatures::TrackFeatures fTrackFeatures;

  std::string _trackmodulelabel;
  //std::string _showermodulelabel;

  bool        _getprediction;
  std::string _dtassoclabel;
  bool        _getmctruth;
  bool        _getmcparticle;
  bool        _overlay;
  std::string _geniegenmodulelabel;
  std::string _hitassoclabel;
  std::string _tpcobject_producer;
  std::string _potsum_producer;
  std::string _potsum_instance;
  std::string _fcsvname;
  bool        _writecsv;
  std::string fMCParticleModuleLabel;
  std::string fHitModuleLabel;                    
  std::string fMCParticleToMCTruthAssModuleLabel; 
  fhicl::ParameterSet _truthparams;
 bool _debug=true;
  
  TTree* _tree1;
  int _run, _subrun, _event, _trackid;
  float _nhits;
  float _length,_startx,_starty,_startz,_endx,_endy,_endz;
  float _theta,_phi,_distlenratio;
  float _startdedx,_dedxratio;
  float _trtotaldedx,_traveragededx;
  float _cosmicscore,_coscontscore;
  float _predict_p,_predict_mu,_predict_pi,_predict_em,_predict_cos;
  float   _mc_ccnc=-9999,_mc_mode=-9999,_mc_hitnuc=-9999;
  float _mc_q2=-9999,_mc_nux=-9999,_mc_nuy=-9999,_mc_nuz=-9999;
  float _mc_enu=-9999;
  int   _mcpdg=-9999,_mcprimary=-9999,_mcorigin=-9999;
  float _mclength=-9999,_mcstartx=-9999,_mcstarty=-9999,_mcstartz=-9999;
  float _mcendx=-9999,_mcendy=-9999,_mcendz=-9999,_mctheta=-9999,_mcphi=-9999,_mckinetic=-9999;

  std::ofstream fout;

 // bool _debug = false;   
  bool _is_data, _is_mc;
  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot;

  const anab::CosmicTagID_t TAGID_P  = anab::CosmicTagID_t::kGeometry_YY;
  const anab::CosmicTagID_t TAGID_MU = anab::CosmicTagID_t::kGeometry_YZ;
  const anab::CosmicTagID_t TAGID_PI = anab::CosmicTagID_t::kGeometry_ZZ;
  const anab::CosmicTagID_t TAGID_EM = anab::CosmicTagID_t::kGeometry_XX;
  const anab::CosmicTagID_t TAGID_CS = anab::CosmicTagID_t::kGeometry_XY;
    
};


DecisionTreeIDAna::DecisionTreeIDAna(fhicl::ParameterSet const & p)
  : EDAnalyzer(p) 
{
  _trackmodulelabel = p.get<std::string>("TrackModuleLabel");
  // _showermodulelabel = p.get<std::string>("ShowerModuleLabel","pandora");
  _getprediction    = p.get<bool>("GetDTPrediction", false);
  _dtassoclabel     = p.get<std::string>("DTAssocLabel", "decisiontreeid");
  _getmctruth       = p.get<bool>("GetMCTruth", false);
  _getmcparticle    = p.get<bool>("GetMCParticle", false);
  _overlay    = p.get<bool>("Overlay", false);
  _geniegenmodulelabel = p.get<std::string>("GenieGenModuleLabel","generator");
  _hitassoclabel    = p.get<std::string>("HitAssocLabel","pmtrack");
  _fcsvname         = p.get<std::string>("CSVFileOut","test.csv");
  _writecsv         = p.get<bool>("WriteCSV", false);
  _truthparams      = p.get<fhicl::ParameterSet>("MCTruthMatching");
  //_tpcobject_producer             = p.get<std::string>("TPCObjectProducer");   
  _potsum_producer                = p.get<std::string>("POTSummaryProducer");
  _potsum_instance                = p.get<std::string>("POTSummaryInstance");
  fTrackFeatures.Configure(p.get<fhicl::ParameterSet>("FeaturesConfig"));
  fMCParticleModuleLabel= p.get<std::string>("MCParticleModuleLabel","largeant");
fMCParticleToMCTruthAssModuleLabel= p.get<std::string>("MCParticleToMCTruthAssModuleLabel","largeant");
 fHitModuleLabel=p.get<std::string>("HitModuleLabel","gaushit");

   art::ServiceHandle<art::TFileService> tfs;
  _tree1 = tfs->make<TTree>("tree","");
  _tree1->Branch("run",    &_run,    "run/I");
  _tree1->Branch("subrun", &_subrun, "subrun/I");
  _tree1->Branch("event",  &_event,  "event/I");
  _tree1->Branch("trackid",&_trackid,"trackid/I");
  _tree1->Branch("nhits",&_nhits,"nhits/F");
  _tree1->Branch("length",&_length,"length/F");
  _tree1->Branch("startx",&_startx,"startx/F");
  _tree1->Branch("starty",&_starty,"starty/F");
  _tree1->Branch("startz",&_startz,"startz/F");
  _tree1->Branch("endx",&_endx,"endx/F");
  _tree1->Branch("endy",&_endy,"endy/F");
  _tree1->Branch("endz",&_endz,"endz/F");
  _tree1->Branch("theta",&_theta,"theta/F");
  _tree1->Branch("phi",&_phi,"phi/F");
  _tree1->Branch("distlenratio",&_distlenratio,"distlenratio/F");
  _tree1->Branch("startdedx",&_startdedx,"startdedx/F");
  _tree1->Branch("dedxratio",&_dedxratio,"dedxratio/F");
  _tree1->Branch("trtotaldedx",&_trtotaldedx,"trtotaldedx/F");
  _tree1->Branch("traveragededx",&_traveragededx,"traveragededx/F");
  _tree1->Branch("cosmicscore",&_cosmicscore,"cosmicscore/F");
  _tree1->Branch("coscontscore",&_coscontscore,"coscontscore/F");
  _tree1->Branch("predict_p",&_predict_p,"predict_p/F");
  _tree1->Branch("predict_mu",&_predict_mu,"predict_mu/F");
  _tree1->Branch("predict_pi",&_predict_pi,"predict_pi/F");
  _tree1->Branch("predict_em",&_predict_em,"predict_em/F");
  _tree1->Branch("predict_cos",&_predict_cos,"predict_cos/F");
  _tree1->Branch("mc_ccnc",&_mc_ccnc,"mc_ccnc/F");
  _tree1->Branch("mc_mode",&_mc_mode,"mc_mode/F");
  _tree1->Branch("mc_hitnuc",&_mc_hitnuc,"mc_hitnuc/F");
  _tree1->Branch("mc_q2",&_mc_q2,"mc_q2/F");
  _tree1->Branch("mc_nux",&_mc_nux,"mc_nux/F");
  _tree1->Branch("mc_nuy",&_mc_nuy,"mc_nuy/F");
  _tree1->Branch("mc_nuz",&_mc_nuz,"mc_nuz/F");
  _tree1->Branch("mc_enu",&_mc_enu,"mc_enu/F");
  _tree1->Branch("mcpdg",&_mcpdg,"mcpdg/I");
  _tree1->Branch("mcprimary",&_mcprimary,"mcprimary/I");
  _tree1->Branch("mcorigin",&_mcorigin,"mcorigin/I");
  _tree1->Branch("mclength",&_mclength,"mclength/F");
  _tree1->Branch("mcstartx",&_mcstartx,"mcstartx/F");
  _tree1->Branch("mcstarty",&_mcstarty,"mcstarty/F");
  _tree1->Branch("mcstartz",&_mcstartz,"mcstartz/F");
  _tree1->Branch("mcendx",&_mcendx,"mcendx/F");
  _tree1->Branch("mcendy",&_mcendy,"mcendy/F");
  _tree1->Branch("mcendz",&_mcendz,"mcendz/F");
  _tree1->Branch("mctheta",&_mctheta,"mctheta/F");
  _tree1->Branch("mcphi",&_mcphi,"mcphi/F");
  _tree1->Branch("mckinetic",&_mckinetic,"mckinetic/F");

  _sr_tree = tfs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                &_sr_run,                "run/I");
  _sr_tree->Branch("subrun",             &_sr_subrun,             "subrun/I");
  _sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/D");
  _sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/D");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/D");


}

void DecisionTreeIDAna::analyze(art::Event const & e)
{

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  _is_data = e.isRealData();
  _is_mc   = !_is_data;
  if(_overlay) _is_mc=true;
  // open csv for writing 
  if(_writecsv) fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);
  
  // recover handle for tracks that we want to analyze
  art::Handle< std::vector<recob::Track> > trackVecHandle;
  e.getByLabel(_trackmodulelabel, trackVecHandle);
  //art::Handle< std::vector<recob::Shower> > showerVecHandle;
  //e.getByLabel(_showermodulelabel, showerVecHandle);

  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  std::vector< std::vector<float> > mcparticlevec;

  art::Handle< std::vector<simb::GTruth> > gTruthHandle;
  std::vector<art::Ptr<simb::GTruth > > glist;
  if(_getmctruth&&_is_mc)
  {
    if (e.getByLabel(_geniegenmodulelabel,mctruthListHandle))
    {
      art::fill_ptr_vector(mclist, mctruthListHandle);
    }
    if (e.getByLabel(_geniegenmodulelabel,gTruthHandle))
    {
      art::fill_ptr_vector(glist, gTruthHandle);
    }
  }
  
  art::Handle<std::vector<simb::MCParticle>> mcparticle_handle;
  std::vector<art::Ptr<simb::MCParticle> > largeant_vec;
  if (e.getByLabel(fMCParticleModuleLabel,mcparticle_handle)) { art::fill_ptr_vector(largeant_vec,mcparticle_handle); } 
 art::Handle<std::vector<recob::Hit > > hit_handle;
    std::vector<art::Ptr<recob::Hit> > hit_vec;
     if (e.getByLabel(fHitModuleLabel,hit_handle)) { art::fill_ptr_vector(hit_vec,hit_handle); }

     if (_debug)  std::cout<<"Track Module: "<<_trackmodulelabel<<std::endl;
   if (_debug)  std::cout<<"Hit Module: "<<_hitassoclabel<<std::endl;
    
   //  std::cout<<"Shower Module pandora, number of showers: "<<showerVecHandle->size()<<std::endl;
  if(trackVecHandle.isValid())
  {
    art::FindManyP<anab::CosmicTag> dtAssns(trackVecHandle, e, _dtassoclabel); 
    art::FindManyP<recob::Hit> hitAssns(trackVecHandle, e, _hitassoclabel);
    std::vector< std::vector<float> > evtdata = fTrackFeatures.CreateFeatures(e, trackVecHandle);
    
    // Loop over input tracks
    for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
      {
	art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);
	if (_debug) std::cout<<"number of tracks: "<<trackVecHandle->size()<<std::endl;
	std::vector<float> trkdata = evtdata.at(trkIdx);
	
	_predict_p   = -9999.;
	_predict_mu  = -9999.;
	_predict_pi  = -9999.;
	_predict_em  = -9999.;
	_predict_cos = -9999.;
	
	if(_getprediction && dtAssns.isValid())
	  {
	    
	    std::vector< art::Ptr<anab::CosmicTag> > dtVec = dtAssns.at(trk.key());
	    for(auto const& dttag : dtVec)
	      {
		if(dttag->CosmicType() == TAGID_P)  _predict_p   = dttag->CosmicScore();
		else if(dttag->CosmicType() == TAGID_MU) _predict_mu  = dttag->CosmicScore();
		else if(dttag->CosmicType() == TAGID_PI) _predict_pi  = dttag->CosmicScore();
		else if(dttag->CosmicType() == TAGID_EM) _predict_em  = dttag->CosmicScore();
		else if(dttag->CosmicType() == TAGID_CS) _predict_cos = dttag->CosmicScore();
	      }
	  }
	
	if(TMath::Sign(1,trkdata[3] - trkdata[5]) 
	   == TMath::Sign(1,trk->Vertex().Z() - trk->End().Z()))
	  {
	    _startx = trk->Vertex().X();
	    _endx = trk->End().X();
	  }
	else
	  {
	    _startx = trk->End().X();
	    _endx = trk->Vertex().X();
	  }
	if(_getmctruth&&_is_mc)
	  {
	    if(!mclist.empty())
	      {
		if(mclist[0]->NeutrinoSet())
		  {
		    art::Ptr<simb::MCTruth> mctruth = mclist[0];
		    _mc_ccnc   = mctruth->GetNeutrino().CCNC();
		    _mc_mode   = mctruth->GetNeutrino().Mode();
		    _mc_hitnuc = mctruth->GetNeutrino().HitNuc();
		    _mc_q2     = mctruth->GetNeutrino().QSqr();
		    _mc_nux    = mctruth->GetNeutrino().Nu().Vx();
		    _mc_nuy    = mctruth->GetNeutrino().Nu().Vy();
		    _mc_nuz    = mctruth->GetNeutrino().Nu().Vz();
		    _mc_enu    = mctruth->GetNeutrino().Nu().E();
		  }
	      }
	  }
	
	if(_getmcparticle&& _is_mc && hitAssns.isValid())
	  {
	    
	    
	    std::vector<art::Ptr<recob::Hit>> hitVec = hitAssns.at(trk.key());
	    ////
	    // BackTrackerTruthMatch backtrackertruthmatch;
	    // backtrackertruthmatch.MatchToMCParticle(hitAssns,e,hitVec);
	    // art::Ptr< simb::MCParticle > maxp_me = backtrackertruthmatch.ReturnMCParticle();
	    // std::cout<<"maxpme: "<<maxp_me.isNull()<<std::endl;
	    // ////
	    // get set of track ids associated with hitvec
	    // const std::set<int> g4ids = bt->GetSetOfTrackIDs(hitVec);
	    // find most g4track with highest hit purity
	    // float hitpurity = 0.;
	    // int puridx = -1;
	    // float tmppur;
	    // std::set<int> tmpg4set;
	    //for(auto const& g4id : g4ids)
	    //{
	    // tmpg4set.clear();
	    //tmpg4set.emplace(g4id);
	    // tmppur = bt->HitCollectionPurity(tmpg4set,hitVec);
	    //if(tmppur >= hitpurity)
	    // {
	    //  hitpurity = tmppur;
	    // puridx = g4id;
	    //  / }
	    // }
	    //	std::cout<<"puridx: "<<puridx<<std::endl;
	    if(1)
	      {//std::cout<<"puridx: "<<puridx<<std::endl;
		////
		art::FindOneP<simb::MCTruth> MCParticleToMCTruth(mcparticle_handle,e,fMCParticleToMCTruthAssModuleLabel);
		BackTrackerTruthMatch backtrackertruthmatch;
		backtrackertruthmatch.MatchToMCParticle(hit_handle,e,hitVec);
		art::Ptr< simb::MCParticle > maxp_me = backtrackertruthmatch.ReturnMCParticle();
		///
		
		// const simb::MCParticle* mcparticle = bt->TrackIDToParticle( puridx );
		//art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth( puridx );
		//          if(mcparticle && !(mctruth.isNull()))
		if(!(maxp_me.isNull()))
		  
		  
		  {
		    //double purity = backtrackertruthmatch.ReturnPurity();
		    // double completeness = backtrackertruthmatch.ReturnCompleteness();
		    //int mcparticleid = backtrackertruthmatch.ReturnMCParticleID();
		    
		    if (_debug)   std::cout << "Track matched by the BackTracker" << std::endl;
		    const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruth.at(maxp_me.key());
		    if (_debug)   std::cout <<mctruth->Origin()<<" "<<maxp_me.key()<<std::endl;
		    //	     std::cout<<mcparticle->PdgCode()<<" "<<maxp_me->PdgCode()<<std::endl;
		    _mcpdg     = maxp_me->PdgCode();
		    _mcprimary = float(maxp_me->Process() == "primary");
		    _mcorigin  = mctruth->Origin();
		    _mclength  = maxp_me->Trajectory().TotalLength();
		    _mcstartx  = maxp_me->Vx();
		    _mcstarty  = maxp_me->Vy();
		    _mcstartz  = maxp_me->Vz();
		    _mcendx    = maxp_me->EndX();
		    _mcendy    = maxp_me->EndY();
		    _mcendz    = maxp_me->EndZ();
		    _mctheta   = maxp_me->Momentum().Theta();
		    _mcphi     = maxp_me->Momentum().Phi();
		    _mckinetic = maxp_me->E() - maxp_me->Mass();
		  }
		else
		  {
		    if (_debug) std::cout << "Track not matched by the BackTracker- Cosmic!" << std::endl;
		    _mcpdg     = -9999;
		    _mcprimary = -9999;
		    _mcorigin  = 2;
		    _mclength  = -9999;
		    _mcstartx  = -9999;
		    _mcstarty  = -9999;
		    _mcstartz  =-9999;
		    _mcendx    = -9999;
		    _mcendy    =-9999;
		    _mcendz    = -9999;
		    _mctheta   =-9999;
		    _mcphi     = -9999;
		    _mckinetic =-9999;
		  }
	      }
      }
	   
	

      _trackid      = trk->ID();
      _nhits        = trkdata[0];
      _length       = trkdata[1];
      _starty       = trkdata[2];
      _startz       = trkdata[3];
      _endy         = trkdata[4];
      _endz         = trkdata[5];
      _theta        = trkdata[6];
      _phi          = trkdata[7];
      _distlenratio = trkdata[8];
      _startdedx    = trkdata[9];
      _dedxratio    = trkdata[10];
      _trtotaldedx  = trkdata[11];
      _traveragededx= trkdata[12];
      _cosmicscore  = trkdata[13];
      _coscontscore = trkdata[14];

      _tree1->Fill();

    }
}
  if(_writecsv) fout.close();
}

void DecisionTreeIDAna::beginJob()
{
  if(_writecsv)
  {
    // open csv and add header
    fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);

    fout << "nhits,length,starty,startz,endy,endz,theta,phi,"
            "distlenratio,startdedx,dedxratio,"
            "trtotaldedx,traveragededx,"
            "cosmicscore,coscontscore,startx,endx,"
            "run,subrun,event";

    if(_getprediction)
    {
      fout << ",predict_p,predict_mu,predict_pi,predict_em,predict_cos";
    }
    if(_getmctruth&&_is_mc)
    {
      fout << ",mc_ccnc,mc_mode,mc_hitnuc,mc_q2,mc_nux,mc_nuy,mc_nuz,mc_enu";
    }
    if(_getmcparticle&&_is_mc)
    {
      fout << ",mcpdg,mcprimary,mcorigin,"
              "mclength,mcstartx,mcstarty,mcstartz,"
              "mcendx,mcendy,mcendz,mctheta,mcphi,mckinetic";
    }
    fout << std::endl;
    fout.close();
  }
}

void DecisionTreeIDAna::endSubRun(const art::SubRun& sr) {

  if (_debug) std::cout << "[UBXSec::endSubRun] Starts" << std::endl;

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;

  // MC
  if (_is_mc) {
    if (_debug) std::cout << "[NCE::endSubRun] Getting POT for MC" << std::endl;
    if(sr.getByLabel(_potsum_producer, potsum_h)) {
      if (_debug) std::cout << "[NCE::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot/1.0e16;
//std::cout <<"MC pot:  "<<_sr_pot << std::endl;
    }
    else
      _sr_pot = 0.;
 std::cout <<"MC pot:  "<<_sr_pot << std::endl; 
 }

  // Data - Use Zarko's script instead
  if (_is_data) {
    if (_debug) std::cout << "[NCE::endSubRun] Getting POT for DATA, producer " << _potsum_producer << ", instance " << _potsum_instance << std::endl;
    if (sr.getByLabel(_potsum_producer, _potsum_instance, potsum_h)){
      if (_debug) std::cout << "[NCE::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot/1.0e16;
    }
    else
      _sr_pot = 0;
 std::cout <<"data pot:  "<<_sr_pot << std::endl; 
 }
  if (_debug) std::cout <<"pot:  "<<_sr_pot << std::endl;
  _sr_tree->Fill();

  if (_debug) std::cout << "[NCE::endSubRun] Ends" << std::endl;
}



DEFINE_ART_MODULE(DecisionTreeIDAna)

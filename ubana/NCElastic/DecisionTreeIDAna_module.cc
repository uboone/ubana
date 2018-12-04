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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "ubobj/UBXSec/TPCObject.h"
#include "TString.h"
#include "TTree.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

//#include "ubana/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
//#include "ubana/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"

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

private:

  // Pointer to the feature algorithm
  ::dtfeatures::TrackFeatures fTrackFeatures;

  std::string _trackmodulelabel;

  bool        _getprediction;
  std::string _dtassoclabel;
  bool        _getmctruth;
  bool        _getmcparticle;
  std::string _geniegenmodulelabel;
  std::string _hitassoclabel;
  std::string _tpcobject_producer;

  std::string _fcsvname;
  bool        _writecsv;

  fhicl::ParameterSet _truthparams;

  TTree* _tree1;
  int _run, _subrun, _event, _trackid;
  float _nhits,_length,_startx,_starty,_startz,_endx,_endy,_endz;
  float _theta,_phi,_distlenratio;
  float _startdedx,_dedxratio;
  float _trtotaldedx,_traveragededx;
  float _cosmicscore,_coscontscore;
  float _predict_p,_predict_mu,_predict_pi,_predict_em,_predict_cos;
  float _mc_ccnc,_mc_mode,_mc_hitnuc,_mc_q2,_mc_nux,_mc_nuy,_mc_nuz;
  float _mc_enu;
  float _mcpdg,_mcprimary,_mcorigin;
  float _mclength,_mcstartx,_mcstarty,_mcstartz;
  float _mcendx,_mcendy,_mcendz,_mctheta,_mcphi,_mckinetic;

  std::ofstream fout;

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
  _getprediction    = p.get<bool>("GetDTPrediction", false);
  _dtassoclabel     = p.get<std::string>("DTAssocLabel", "decisiontreeid");
  _getmctruth       = p.get<bool>("GetMCTruth", false);
  _getmcparticle    = p.get<bool>("GetMCParticle", false);
  _geniegenmodulelabel = p.get<std::string>("GenieGenModuleLabel","generator");
  _hitassoclabel    = p.get<std::string>("HitAssocLabel","pandora");
  _fcsvname         = p.get<std::string>("CSVFileOut","test.csv");
  _writecsv         = p.get<bool>("WriteCSV", false);
  _truthparams      = p.get<fhicl::ParameterSet>("MCTruthMatching");
  _tpcobject_producer             = p.get<std::string>("TPCObjectProducer");   
  fTrackFeatures.Configure(p.get<fhicl::ParameterSet>("FeaturesConfig"));

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
  _tree1->Branch("mcpdg",&_mcpdg,"mcpdg/F");
  _tree1->Branch("mcprimary",&_mcprimary,"mcprimary/F");
  _tree1->Branch("mcorigin",&_mcorigin,"mcorigin/F");
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

}

void DecisionTreeIDAna::analyze(art::Event const & e)
{

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // open csv for writing 
  if(_writecsv) fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);

  // recover handle for tracks that we want to analyze
  art::Handle< std::vector<recob::Track> > trackVecHandle;
  e.getByLabel(_trackmodulelabel, trackVecHandle);

  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  std::vector< std::vector<float> > mcparticlevec;

  art::Handle< std::vector<simb::GTruth> > gTruthHandle;
  std::vector<art::Ptr<simb::GTruth > > glist;
  if(_getmctruth)
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
  

  if(trackVecHandle.isValid())
  {
    art::FindManyP<anab::CosmicTag> dtAssns(trackVecHandle, e, _dtassoclabel); 
    art::FindManyP<recob::Hit> hitAssns(trackVecHandle, e, _hitassoclabel);

    std::vector< std::vector<float> > evtdata = fTrackFeatures.CreateFeatures(e, trackVecHandle);

    std::unique_ptr<truth::IMCTruthMatching> bt;

bt= art::make_tool<truth::IMCTruthMatching>( _truthparams);
    if(_getmcparticle)
    {
  //    bt = std::unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(_truthparams));
      bt->Rebuild(e);
    }

    art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
    e.getByLabel(_tpcobject_producer, tpcobj_h);
    if (!tpcobj_h.isValid()) {
      std::cout << "[NCE TPCObject] Cannote locate ubana::TPCObject." << std::endl;
    }
    else 
      std::cout << "[NCE TPCObject] Found ubana::TPCObject." << std::endl;
    
    // Loop over input tracks
    for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
    {
      art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);

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
      if(_getmctruth)
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
      if(_getmcparticle & hitAssns.isValid())
      {
        art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);
        std::vector<art::Ptr<recob::Hit>> hitVec = hitAssns.at(trk.key());

        // get set of track ids associated with hitvec
        const std::set<int> g4ids = bt->GetSetOfTrackIDs(hitVec);
        // find most g4track with highest hit purity
        float hitpurity = 0.;
        int puridx = -1;
        float tmppur;
        std::set<int> tmpg4set;
        for(auto const& g4id : g4ids)
        {
          tmpg4set.clear();
          tmpg4set.emplace(g4id);
          tmppur = bt->HitCollectionPurity(tmpg4set,hitVec);
          if(tmppur >= hitpurity)
          {
            hitpurity = tmppur;
            puridx = g4id;
          }
        }
        if(puridx != -1)
        {
          const simb::MCParticle* mcparticle = bt->TrackIDToParticle( puridx );
          art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth( puridx );
          if(mcparticle && !(mctruth.isNull()))
          {
            _mcpdg     = mcparticle->PdgCode();
            _mcprimary = float(mcparticle->Process() == "primary");
            _mcorigin  = mctruth->Origin();
            _mclength  = mcparticle->Trajectory().TotalLength();
            _mcstartx  = mcparticle->Vx();
            _mcstarty  = mcparticle->Vy();
            _mcstartz  = mcparticle->Vz();
            _mcendx    = mcparticle->EndX();
            _mcendy    = mcparticle->EndY();
            _mcendz    = mcparticle->EndZ();
            _mctheta   = mcparticle->Momentum().Theta();
            _mcphi     = mcparticle->Momentum().Phi();
            _mckinetic = mcparticle->E() - mcparticle->Mass();
          }
        }
      }

      if(_writecsv)
      {
        for(auto it = trkdata.begin(); it != trkdata.end(); ++it)
        {
          if(std::next(it) != trkdata.end()) fout << *it << ",";
          else fout << *it;
        }
        fout << "," << _startx << "," << _endx;
        fout << "," << _run << "," << _subrun << "," << _event;
        if(_getprediction) fout << "," << _predict_p  << ","
                                       << _predict_mu << ","
                                       << _predict_pi << ","
                                       << _predict_em << ","
                                       << _predict_cos;
        if(_getmctruth)
        {
          fout << "," << _mc_ccnc << "," 
                      << _mc_mode << "," 
                      << _mc_hitnuc << "," 
                      << _mc_q2 << "," 
                      << _mc_nux << "," 
                      << _mc_nuy << "," 
                      << _mc_nuz << "," 
                      << _mc_enu;
        }
        if(_getmcparticle)
        {
          fout << "," << _mcpdg << ","
                      << _mcprimary << ","
                      << _mcorigin << ","
                      << _mclength << ","
                      << _mcstartx << ","
                      << _mcstarty << ","
                      << _mcstartz << ","
                      << _mcendx << ","
                      << _mcendy << ","
                      << _mcendz << ","
                      << _mctheta << ","
                      << _mcphi << ","
                      << _mckinetic;
        }

        fout << std::endl;
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
    if(_getmctruth)
    {
      fout << ",mc_ccnc,mc_mode,mc_hitnuc,mc_q2,mc_nux,mc_nuy,mc_nuz,mc_enu";
    }
    if(_getmcparticle)
    {
      fout << ",mcpdg,mcprimary,mcorigin,"
              "mclength,mcstartx,mcstarty,mcstartz,"
              "mcendx,mcendy,mcendz,mctheta,mcphi,mckinetic";
    }
    fout << std::endl;
    fout.close();
  }
}
DEFINE_ART_MODULE(DecisionTreeIDAna)

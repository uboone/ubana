////////////////////////////////////////////////////////////////////////
// Class:       QEEventSelection
// Plugin Type: analyze (art v2_05_00)
// File:        QEEventSelection_module.cc
//
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
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "TString.h"
#include "TTree.h"

#include "ubana/NCElastic/Algorithms/LREvtReconstruction.h"

#include <memory>
#include <cmath>
#include <stdexcept>

class QEEventSelection;

class QEEventSelection : public art::EDAnalyzer {
public:
  explicit QEEventSelection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  QEEventSelection(QEEventSelection const &) = delete;
  QEEventSelection(QEEventSelection &&) = delete;
  QEEventSelection & operator = (QEEventSelection const &) = delete;
  QEEventSelection & operator = (QEEventSelection &&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const & e) override;

private:

  // Pointer to the feature algorithm
  ::qeselection::LREvtReconstruction fEvtReconstruction;

  std::string _GenieModuleLabel;
  bool        _getmctruth;

  bool        _saveall;
  float       _minLR;

  std::string _fcsvname;
  bool        _writecsv;

  TTree* _tree1;
  int _run, _subrun, _event;
  float _reco_q2,_reco_Tp,_reco_Rp,_reco_LR;
  int _reco_TID;
  float _reco_TH,_reco_PH;
  float _reco_nux,_reco_nuy,_reco_nuz;
  float _mc_ccnc,_mc_mode,_mc_hitnuc,_mc_q2,_mc_enu,_mc_enunucrf;
  float _mc_nux,_mc_nuy,_mc_nuz;
  int _mc_int,_mc_nuPDG;
  float _mc_dxsec;
  double _bnb_correction;

  std::ofstream fout;

};


QEEventSelection::QEEventSelection(fhicl::ParameterSet const & p)
  : EDAnalyzer(p) 
{
  _GenieModuleLabel = p.get< std::string > ("GenieModuleLabel");
  _getmctruth       = p.get<bool>("GetMCTruth", false);
  _saveall          = p.get<bool>("SaveAllEvents", false);
  _minLR            = p.get<float>("MinimumLRScore",0.9);
  _fcsvname         = p.get<std::string>("CSVFileOut");
  _writecsv         = p.get<bool>("WriteCSV", false);

  fEvtReconstruction.Configure(p.get<fhicl::ParameterSet>("SelectionConfig"));
   
  art::ServiceHandle<art::TFileService> tfs;
  _tree1 = tfs->make<TTree>("tree","");
  _tree1->Branch("run",    &_run,    "run/I");
  _tree1->Branch("subrun", &_subrun, "subrun/I");
  _tree1->Branch("event",  &_event,  "event/I");
  _tree1->Branch("reco_q2",&_reco_q2,"reco_q2/F");
  _tree1->Branch("reco_Tp",&_reco_Tp,"reco_Tp/F");
  _tree1->Branch("reco_Rp",&_reco_Rp,"reco_Rp/F");
  _tree1->Branch("reco_LR",&_reco_LR,"reco_LR/F");
  _tree1->Branch("reco_TID",&_reco_TID,"reco_TID/I");
  _tree1->Branch("reco_TH",&_reco_TH,"reco_TH/F");
  _tree1->Branch("reco_PH",&_reco_PH,"reco_PH/F");
  _tree1->Branch("reco_nux",&_reco_nux,"reco_nux/F");
  _tree1->Branch("reco_nuy",&_reco_nuy,"reco_nuy/F");
  _tree1->Branch("reco_nuz",&_reco_nuz,"reco_nuz/F");
  _tree1->Branch("mc_ccnc",&_mc_ccnc,"mc_ccnc/F");
  _tree1->Branch("mc_mode",&_mc_mode,"mc_mode/F");
  _tree1->Branch("mc_hitnuc",&_mc_hitnuc,"mc_hitnuc/F");
  _tree1->Branch("mc_nuPDG",&_mc_nuPDG,"mc_nuPDG/I");
  _tree1->Branch("mc_q2",&_mc_q2,"mc_q2/F");
  _tree1->Branch("mc_enu",&_mc_enu,"mc_enu/F");
  _tree1->Branch("mc_enunucrf",&_mc_enunucrf,"mc_enunucrf/F");
  _tree1->Branch("mc_int",&_mc_int,"mc_int/I");
  _tree1->Branch("mc_dxsec",&_mc_dxsec,"mc_dxsec/F");
  _tree1->Branch("bnbcorrection",&_bnb_correction,"bnbcorrection/D");


}

void QEEventSelection::analyze(art::Event const & e)
{

  //get the MC generator information out of the event       
  //these are all handles to mc information.
  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;  
  art::Handle< std::vector<simb::GTruth> > gTruthHandle;
  art::Handle< std::vector<evwgh::MCEventWeight> > evtWeights;

  std::vector<art::Ptr<simb::MCTruth> > mclist;
  std::vector<art::Ptr<simb::GTruth > > glist;

  if(_getmctruth)
  {
    e.getByLabel(_GenieModuleLabel,mcTruthHandle);
    e.getByLabel(_GenieModuleLabel,gTruthHandle);

    art::fill_ptr_vector(mclist, mcTruthHandle);
    art::fill_ptr_vector(glist, gTruthHandle);
  }


  std::vector<float> recovec = fEvtReconstruction.ReconstructEvent(e);

  _reco_q2   = recovec.at(0);
  _reco_Tp   = recovec.at(1);
  _reco_Rp   = recovec.at(2);
  _reco_LR   = recovec.at(3);
  _reco_TID  = recovec.at(4);
  _reco_TH   = recovec.at(5);
  _reco_PH   = recovec.at(6);
  _reco_nux  = recovec.at(7);
  _reco_nuy  = recovec.at(8);
  _reco_nuz  = recovec.at(9);


  if(_saveall || _reco_LR >= _minLR)
  {
    
    // open csv for writing 
    if(_writecsv) fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);
      
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
    
    if(_getmctruth)
    {
      _mc_ccnc     = mclist[0]->GetNeutrino().CCNC();
      _mc_mode     = mclist[0]->GetNeutrino().Mode();
      if(_mc_mode == 10) _mc_hitnuc = mclist[0]->GetNeutrino().HitNuc() - 2000000000;
      else _mc_hitnuc = mclist[0]->GetNeutrino().HitNuc();
      _mc_nuPDG    = mclist[0]->GetNeutrino().Nu().PdgCode();
      _mc_q2       = mclist[0]->GetNeutrino().QSqr();
      _mc_enu      = mclist[0]->GetNeutrino().Nu().E();
      _mc_int      = mclist[0]->GetNeutrino().InteractionType();
      
      // get GTruth info
      TLorentzVector probep4 = glist[0]->fProbeP4;
      _mc_enunucrf = probep4.Energy();
      
      // get default cross sections
      _mc_dxsec = glist[0]->fDiffXsec;
      
      _mc_nux = mclist[0]->GetNeutrino().Nu().EndX();
      _mc_nuy = mclist[0]->GetNeutrino().Nu().EndY();
      _mc_nuz = mclist[0]->GetNeutrino().Nu().EndZ();
      
      // get BNB correction
      _bnb_correction = -1.;
      if(e.getByLabel("eventweight",evtWeights))
      {
        const std::vector< evwgh::MCEventWeight > * evtwgt_vec = evtWeights.product();
        evwgh::MCEventWeight evtwgt = evtwgt_vec->at(0); // just for the first neutrino interaction
        std::map<std::string, std::vector<double>> evtwgt_map = evtwgt.fWeight;
        
        for(std::map<std::string, std::vector<double>>::iterator it = evtwgt_map.begin(); it != evtwgt_map.end(); ++it) {
          if(it->first == "bnbcorrection_FluxHist") _bnb_correction = (it->second).at(0);
        }
      }
    }
    
    if(_writecsv)
    {
      fout << _run << "," << _subrun << "," << _event << ","
           << _reco_q2 << "," << _reco_Tp << "," << _reco_Rp << "," << _reco_LR << ","
           << _reco_TID << "," << _reco_TH << "," << _reco_PH << ","
           << _reco_nux << "," << _reco_nuy << "," << _reco_nuz;
      if(_getmctruth)
      {
        fout << "," << _mc_ccnc << "," << _mc_mode << "," << _mc_hitnuc << "," 
             << _mc_nuPDG << "," << _mc_q2 << "," << _mc_enu << "," << _mc_enunucrf << "," 
             << _mc_int << "," << _mc_dxsec << "," << _mc_nux << "," 
             << _mc_nuy << "," << _mc_nuz << "," << _bnb_correction;
      }
      fout << std::endl;
    }
    _tree1->Fill();
  }
  if(_writecsv) fout.close();
}


void QEEventSelection::beginJob()
{
  if(_writecsv)
  {
    // open csv and add header
    fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);

    fout << "run,subrun,event";
    fout << ",reco_Q2,reco_Tp,reco_Rp,reco_LR";
    fout << ",reco_TID,reco_theta,reco_phi";
    fout << ",reco_nux,reco_nuy,reco_nuz";
    if(_getmctruth)
    {
      fout << ",mc_ccnc,mc_mode,mc_hitnuc,mc_nuPDG";
      fout << ",mc_Q2,mc_Enu,mc_Enu_nucRF,mc_int";
      fout << ",mc_dxsec,mc_inTPC,mc_nux,mc_nuy,mc_nuz";
      fout << ",bnbcorrection";
    }
    fout << std::endl;

    fout.close();
  }
}

DEFINE_ART_MODULE(QEEventSelection)

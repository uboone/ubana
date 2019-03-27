#ifndef CELLTREEUB_MODULE
#define CELLTREEUB_MODULE

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcore/Geometry/Geometry.h"
#include "ubcore/Geometry/UBOpReadoutMap.h"
#include "ubcore/Geometry/UBOpChannelTypes.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "ubobj/MuCS/MuCSData.h"
#include "ubobj/MuCS/MuCSRecoData.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TTimeStamp.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH1S.h"
#include "TLorentzVector.h"
#include "TObjArray.h"

#include <map>
#include <vector>
#include <fstream>
#include <iostream>

#define MAX_TRACKS 30000

using namespace std;

namespace wc{
  class CellTreeUB : public art::EDAnalyzer{
  public:
    explicit CellTreeUB(fhicl::ParameterSet const& pset);
    virtual ~CellTreeUB();
    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void beginSubRun(const art::SubRun& subrun);
    void endSubRun(const art::SubRun& subrun);
    void analyze(const art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    void initOutput();
    void reset();
    void InitProcessMap();
    void processTPC_raw(const art::Event& evt);
    void processTPC_noiseFiltered(const art::Event& evt);
    void processTPC_deconWiener(const art::Event& evt);
    void processTPC_deconGaussian(const art::Event& evt);
    void processTPC_badChannelList(const art::Event& evt);
    void processTPC_channelThreshold(const art::Event& evt);
    void processPMT_raw(const art::Event& evt);
    void processPMT_overlay(const art::Event& evt);
    void processPMT_wfmSaturation(const art::Event& evt);
    void processPMT_reco(const art::Event& evt);
    void processPMT_gainData(const art::Event& evt);
    void processPMT_gainMC(const art::Event& evt);
    void processHWTrigger(const art::Event& evt);
    void processSWTrigger(const art::Event& evt);
    void processBeam(const art::Event& evt);
    void processPOT(const art::SubRun& subrun);
    void processMuCS_data(const art::Event& evt);
    void processMuCS_recoData(const art::Event& evt);
    void processMCParticle(const art::Event& evt);
    void processSimEnergyDepositInside(const art::Event& evt);
    void processSimEnergyDepositOutside(const art::Event& evt);
    void processSimChannel(const art::Event& evt);
  private:
    // --- fhicl parameters ---
    std::string fTPC_rawLabel;
    std::string fTPC_rawProducer;
    std::string fTPC_noiseFilteredLabel;
    std::string fTPC_noiseFilteredProducer;
    std::string fTPC_deconWienerLabel;
    std::string fTPC_deconGaussianLabel;
    std::string fTPC_deconWienerProducer;
    std::string fTPC_deconGaussianProducer;
    std::string fTPC_badChannelListLabel;
    std::string fTPC_badChannelListProducer;
    std::string fTPC_channelThresholdLabel;
    std::string fTPC_channelThresholdProducer;
    std::string fPMT_HG_cosmicLabel;
    std::string fPMT_LG_cosmicLabel;
    std::string fPMT_HG_beamLabel;
    std::string fPMT_LG_beamLabel;
    std::string fPMT_HG_cosmicProducer;
    std::string fPMT_LG_cosmicProducer;
    std::string fPMT_HG_beamProducer;
    std::string fPMT_LG_beamProducer;
    std::string fPMT_wfmSaturationLabel;
    std::string fPMT_overlay_detsim_producer;
    std::string fPMT_overlay_mixer_producer;
    std::string fPMT_overlay_saturation_producer;
    std::string fPMT_recoLabel;
    std::string fHWTrigger_label;
    std::string fSWTrigger_label;
    std::string fBeam_label;
    std::string fPOT_bnb1label;
    std::string fPOT_bnb2label;
    std::string fPOT_numilabel;
    std::string fPOT_producer;
    std::string fMuCS_dataLabel;
    std::string fMuCS_recoDataLabel;
    std::string fSimEnergyDepositInside_label;
    std::string fSimEnergyDepositInside_producer;
    std::string fSimEnergyDepositOutside_label;
    std::string fSimEnergyDepositOutside_producer;
    std::string fSimChannel_label;
    std::string fSimChannel_producer;
    bool _use_LG_beam_for_HG_cosmic;
    bool fSaveTPC_raw;
    bool fSaveTPC_noiseFiltered;
    bool fSaveTPC_deconWiener;
    bool fSaveTPC_deconGaussian;
    bool fSaveTPC_badChannelList;
    bool fSaveTPC_channelThreshold;
    bool fSavePMT_raw;
    bool fSavePMT_overlay;
    bool fSavePMT_wfmSaturation;
    bool fSavePMT_reco;
    bool fSavePMT_gainData;
    bool fSavePMT_gainMC;
    bool fSaveHWTrigger;
    bool fSaveSWTrigger;
    bool fSaveBeam;
    bool fSavePOT;
    bool fSaveMuCSData;
    bool fSaveMuCSRecoData;
    bool fSaveMCParticle;
    bool fSaveSimEnergyDepositInside;
    bool fSaveSimEnergyDepositOutside;
    bool fSaveSimChannel;
    int deconRebin;
    float flashMultPEThreshold;
    bool saveVYshorted;
    short scaleRawdigit;
    short shiftRawdigit;
    std::string mcOption;
    bool fSaveMCTrackPoints;

    TTree *fEventTree;
    int fEvent;
    int fRun;
    int fSubRun;
    double fEventTime;
    // --- TPC raw ---
    int fTPCraw_nChannel;
    vector<int> fTPCraw_channelId;
    TClonesArray *fTPCraw_wf;
    // --- TPC noise filtered ---
    int fTPCnoiseFiltered_nChannel;
    vector<int> fTPCnoiseFiltered_channelId;
    TClonesArray *fTPCnoiseFiltered_wf;
    short fScaleRawdigit;
    short fShiftRawdigit;
    // --- TPC deconvolution Wiener filter ---
    int fTPCdeconWiener_nChannel;
    vector<int> fTPCdeconWiener_channelId;
    TClonesArray *fTPCdeconWiener_wf;
    // --- TPC deconvolution Gaussian filter ---
    int fTPCdeconGaussian_nChannel;
    vector<int> fTPCdeconGaussian_channelId;
    TClonesArray *fTPCdeconGaussian_wf;
    // --- TPC bad channel list ---
    vector<int> fBadChannel;
    vector<int> fBadBegin;
    vector<int> fBadEnd;
    // --- TPC channel thresholds ---
    vector<double> fChannelThreshold;
    // --- PMT raw ---
    vector<short> fOp_cosmic_hg_opch;
    vector<short> fOp_cosmic_lg_opch;
    vector<short> fOp_beam_hg_opch;
    vector<short> fOp_beam_lg_opch;
    vector<double> fOp_cosmic_hg_timestamp;
    vector<double> fOp_cosmic_lg_timestamp;
    vector<double> fOp_beam_hg_timestamp;
    vector<double> fOp_beam_lg_timestamp;
    TClonesArray *fOp_cosmic_hg_wf;
    TClonesArray *fOp_cosmic_lg_wf;
    TClonesArray *fOp_beam_hg_wf;
    TClonesArray *fOp_beam_lg_wf;
    // --- PMT overlay --- ### detsim
    vector<short> fOp_detsim_cosmic_hg_opch;
    vector<short> fOp_detsim_cosmic_lg_opch;
    vector<short> fOp_detsim_beam_hg_opch;
    vector<short> fOp_detsim_beam_lg_opch;
    vector<double> fOp_detsim_cosmic_hg_timestamp;
    vector<double> fOp_detsim_cosmic_lg_timestamp;
    vector<double> fOp_detsim_beam_hg_timestamp;
    vector<double> fOp_detsim_beam_lg_timestamp;
    TClonesArray *fOp_detsim_cosmic_hg_wf;
    TClonesArray *fOp_detsim_cosmic_lg_wf;
    TClonesArray *fOp_detsim_beam_hg_wf;
    TClonesArray *fOp_detsim_beam_lg_wf;
    // --- PMT overlay --- ### mixer
    vector<short> fOp_mixer_beam_hg_opch;
    vector<short> fOp_mixer_beam_lg_opch;
    vector<double> fOp_mixer_beam_hg_timestamp;
    vector<double> fOp_mixer_beam_lg_timestamp;
    TClonesArray *fOp_mixer_beam_hg_wf;
    TClonesArray *fOp_mixer_beam_lg_wf;
    // --- PMT overlay --- ### saturation
    vector<short> fOp_saturation_beam_hg_opch;
    vector<short> fOp_saturation_beam_lg_opch;
    vector<double> fOp_saturation_beam_hg_timestamp;
    vector<double> fOp_saturation_beam_lg_timestamp;
    TClonesArray *fOp_saturation_beam_hg_wf;
    TClonesArray *fOp_saturation_beam_lg_wf;
    // --- PMT waveform saturation ---
    TClonesArray *fOp_wf;
    vector<short> fOp_femch;
    vector<double> fOp_timestamp;
    // --- PMT reconstruction ---
    int oh_nHits;
    vector<short> oh_channel;
    vector<double> oh_bgtime;
    vector<double> oh_trigtime;
    vector<double> oh_pe;
    int of_nFlash;
    vector<float> of_t;
    vector<float> of_peTotal;
    vector<short> of_multiplicity;
    TClonesArray *fPEperOpDet;
    // --- PMT gain ---
    vector<float> fOp_gain;
    vector<float> fOp_gainerror;
    // --- HW Trigger ---
    unsigned int fTriggernumber;
    double fTriggertime;
    double fBeamgatetime;
    unsigned int fTriggerbits;
    // --- SW Trigger ---
    vector<pair<string,bool> > fpassAlgo;
    vector<bool> fpassPrescale;
    vector<uint32_t> fPHMAX;
    vector<uint32_t> fmultiplicity;
    vector<uint32_t> ftriggerTick;
    vector<double> ftriggerTime;
    vector<float> fprescale_weight;
    // --- Beam info ---
    double ftor101;
    double ftortgt;
    double ftrtgtd;
    long long int ft_ms;
    uint8_t fRecordType;
    uint32_t fSeconds;
    uint16_t fMilliSeconds;
    uint16_t fNumberOfDevices;
    map<string,vector<double> > fDataMap;
    // --- POT info ---
    double ftotpot_bnbETOR875;
    double ftotpot_numiETORTGT;
    double ftotpot_bnbETOR860;
    double ftotgoodpot_bnbETOR875;
    double ftotgoodpot_numiETORTGT;
    double ftotgoodpot_bnbETOR860;
    int ftotspills_bnbETOR875;
    int ftotspills_numiETORTGT;
    int ftotspills_bnbETOR860;
    int fgoodspills_bnbETOR875;
    int fgoodspills_numiETORTGT;
    int fgoodspills_bnbETOR860;
    // --- MuCS data ---
    float ft0;
    // --- MuCS reco data ---
    float ftheta_xy;
    float ftheta_xy_rms;
    float fx;
    float fx_rms;
    float ftheta_yz;
    float ftheta_yz_rms;
    float fz;
    float fz_rms;
    float fy;
    float ftheta;
    float fphi;
    int fxmatches;
    int fzmatches;
    // --- MC ---
    std::map<std::string, int> processMap;
    std::map<int, int> savedMCTrackIdMap; // key: id, value: pdg
    int mc_Ntrack;
    int mc_id[MAX_TRACKS];
    int mc_pdg[MAX_TRACKS];
    int mc_process[MAX_TRACKS];
    int mc_mother[MAX_TRACKS];
    float mc_startXYZT[MAX_TRACKS][4];
    float mc_endXYZT[MAX_TRACKS][4];
    float mc_startMomentum[MAX_TRACKS][4];
    float mc_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > mc_daughters;
    TObjArray *fMC_trackPosition;
    int mc_isnu;
    int mc_nGeniePrimaries;
    int mc_nu_pdg;
    int mc_nu_ccnc;
    int mc_nu_mode;
    int mc_nu_intType;
    int mc_nu_target;
    int mc_hitnuc;
    int mc_hitquark;
    double mc_nu_Q2;
    double mc_nu_W;
    double mc_nu_X;
    double mc_nu_Y;
    double mc_nu_Pt;
    double mc_nu_Theta;
    float mc_nu_pos[4];
    float mc_nu_mom[4];
    // --- SimEnergyDeposit (inside) ---
    std::vector<double> sedi_nelectrons;
    std::vector<double> sedi_nphotons;
    std::vector<double> sedi_time_start;
    std::vector<double> sedi_time_end;
    std::vector<double> sedi_x_start;
    std::vector<double> sedi_x_end;
    std::vector<double> sedi_y_start;
    std::vector<double> sedi_y_end;
    std::vector<double> sedi_z_start;
    std::vector<double> sedi_z_end;
    std::vector<int> sedi_pdg;
    std::vector<int> sedi_trackId;
    std::vector<double> sedi_energy;
    // --- SimEnergyDeposit (outside) ---
    std::vector<double> sedo_nelectrons;
    std::vector<double> sedo_nphotons;
    std::vector<double> sedo_time_start;
    std::vector<double> sedo_time_end;
    std::vector<double> sedo_x_start;
    std::vector<double> sedo_x_end;
    std::vector<double> sedo_y_start;
    std::vector<double> sedo_y_end;
    std::vector<double> sedo_z_start;
    std::vector<double> sedo_z_end;
    std::vector<int> sedo_pdg;
    std::vector<int> sedo_trackId;
    std::vector<double> sedo_energy;
    // --- SimChannel ---
    std::vector<float> sc_x;
    std::vector<float> sc_y;
    std::vector<float> sc_z;
    std::vector<float> sc_q;
    std::vector<int> sc_chan;
    std::vector<unsigned int> sc_tick;
    std::vector<unsigned int> sc_tdc;
    float sc_drift_speed=1.098; // mm/us

  }; // class CellTreeUB

  CellTreeUB::CellTreeUB(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset),
      fEventTree(nullptr){
    reconfigure(pset);
    InitProcessMap();
    initOutput();
  }

  CellTreeUB::~CellTreeUB(){
  }

  void CellTreeUB::reconfigure(fhicl::ParameterSet const& pset){
    fTPC_rawLabel = pset.get<std::string>("TPC_rawLabel");
    fTPC_rawProducer = pset.get<std::string>("TPC_rawProducer");
    fTPC_noiseFilteredLabel = pset.get<std::string>("TPC_noiseFilteredLabel");
    fTPC_noiseFilteredProducer = pset.get<std::string>("TPC_noiseFilteredProducer");
    fTPC_deconWienerLabel = pset.get<std::string>("TPC_deconWienerLabel");
    fTPC_deconGaussianLabel = pset.get<std::string>("TPC_deconGaussianLabel");
    fTPC_deconWienerProducer = pset.get<std::string>("TPC_deconWienerProducer");
    fTPC_deconGaussianProducer = pset.get<std::string>("TPC_deconGaussianProducer");
    fTPC_badChannelListLabel = pset.get<std::string>("TPC_badChannelListLabel");
    fTPC_badChannelListProducer = pset.get<std::string>("TPC_badChannelListProducer");
    fTPC_channelThresholdLabel = pset.get<std::string>("TPC_channelThresholdLabel");
    fTPC_channelThresholdProducer = pset.get<std::string>("TPC_channelThresholdProducer");
    fPMT_HG_cosmicLabel = pset.get<std::string>("PMT_HG_cosmicLabel");
    fPMT_LG_cosmicLabel = pset.get<std::string>("PMT_LG_cosmicLabel");
    fPMT_HG_beamLabel = pset.get<std::string>("PMT_HG_beamLabel");
    fPMT_LG_beamLabel = pset.get<std::string>("PMT_LG_beamLabel");
    fPMT_HG_cosmicProducer = pset.get<std::string>("PMT_HG_cosmicProducer");
    fPMT_LG_cosmicProducer = pset.get<std::string>("PMT_LG_cosmicProducer");
    fPMT_HG_beamProducer = pset.get<std::string>("PMT_HG_beamProducer");
    fPMT_LG_beamProducer = pset.get<std::string>("PMT_LG_beamProducer");
    fPMT_overlay_detsim_producer = pset.get<std::string>("PMT_overlay_detsim_producer");
    fPMT_overlay_mixer_producer = pset.get<std::string>("PMT_overlay_mixer_producer");
    fPMT_overlay_saturation_producer = pset.get<std::string>("PMT_overlay_saturation_producer");
    fPMT_wfmSaturationLabel = pset.get<std::string>("PMT_wfmSaturationLabel");
    fPMT_recoLabel = pset.get<std::string>("PMT_recoLabel");
    fHWTrigger_label = pset.get<std::string>("HWTrigger_label");
    fSWTrigger_label = pset.get<std::string>("SWTrigger_label");
    fBeam_label = pset.get<std::string>("Beam_label");
    fPOT_bnb1label = pset.get<std::string>("POT_bnb1label");
    fPOT_bnb2label = pset.get<std::string>("POT_bnb2label");
    fPOT_numilabel = pset.get<std::string>("POT_numilabel");
    fPOT_producer = pset.get<std::string>("POT_producer");
    fMuCS_dataLabel = pset.get<std::string>("MuCS_dataLabel");
    fMuCS_recoDataLabel = pset.get<std::string>("MuCS_recoDataLabel");
    fSimEnergyDepositInside_label = pset.get<std::string>("SimEnergyDepositInside_label");
    fSimEnergyDepositInside_producer = pset.get<std::string>("SimEnergyDepositInside_producer");
    fSimEnergyDepositOutside_label = pset.get<std::string>("SimEnergyDepositOutside_label");
    fSimEnergyDepositOutside_producer = pset.get<std::string>("SimEnergyDepositOutside_producer");
    fSimChannel_label = pset.get<std::string>("SimChannel_label");
    fSimChannel_producer = pset.get<std::string>("SimChannel_producer");
    fSaveTPC_raw = pset.get<bool>("SaveTPC_raw");
    fSaveTPC_noiseFiltered = pset.get<bool>("SaveTPC_noiseFiltered");
    fSaveTPC_deconWiener = pset.get<bool>("SaveTPC_deconWiener");
    fSaveTPC_deconGaussian = pset.get<bool>("SaveTPC_deconGaussian");
    fSaveTPC_badChannelList = pset.get<bool>("SaveTPC_badChannelList");
    fSaveTPC_channelThreshold = pset.get<bool>("SaveTPC_channelThreshold");
    fSavePMT_raw = pset.get<bool>("SavePMT_raw");
    fSavePMT_overlay = pset.get<bool>("SavePMT_overlay");
    fSavePMT_wfmSaturation = pset.get<bool>("SavePMT_wfmSaturation");
    fSavePMT_reco = pset.get<bool>("SavePMT_reco");
    fSavePMT_gainData = pset.get<bool>("SavePMT_gainData");
    fSavePMT_gainMC = pset.get<bool>("SavePMT_gainMC");
    fSaveHWTrigger = pset.get<bool>("SaveHWTrigger");
    fSaveSWTrigger = pset.get<bool>("SaveSWTrigger");
    fSaveBeam = pset.get<bool>("SaveBeam");
    fSavePOT = pset.get<bool>("SavePOT");
    fSaveMuCSData = pset.get<bool>("SaveMuCSData");
    fSaveMuCSRecoData = pset.get<bool>("SaveMuCSRecoData");
    fSaveMCParticle = pset.get<bool>("SaveMCParticle");
    fSaveSimEnergyDepositInside = pset.get<bool>("SaveSimEnergyDepositInside");
    fSaveSimEnergyDepositOutside = pset.get<bool>("SaveSimEnergyDepositOutside");
    fSaveSimChannel = pset.get<bool>("SaveSimChannel");
    _use_LG_beam_for_HG_cosmic = pset.get<bool>("UseLGBeamForHGCosmic");
    flashMultPEThreshold = pset.get<float>("FlashMultPEThreshold");
    deconRebin = pset.get<int>("DeconRebin");
    saveVYshorted = pset.get<bool>("SaveVYshorted");
    scaleRawdigit = pset.get<short>("ScaleRawdigit");
    shiftRawdigit = pset.get<short>("ShiftRawdigit");
    mcOption = pset.get<std::string>("mcOption");
    fSaveMCTrackPoints = pset.get<bool>("SaveMCTrackPoints");
  }

  void CellTreeUB::initOutput(){
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Sim","");
    fEventTree->Branch("eventNo", &fEvent);
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);
    fEventTree->Branch("eventTime", &fEventTime);
    if(fSaveTPC_raw){
      fEventTree->Branch("raw_nChannel", &fTPCraw_nChannel);
      fEventTree->Branch("raw_channelId",&fTPCraw_channelId);
      fTPCraw_wf = new TClonesArray("TH1S");
      fEventTree->Branch("raw_wf", &fTPCraw_wf, 256000, 0);
    }
    if(fSaveTPC_noiseFiltered){
      fEventTree->Branch("nf_nChannel", &fTPCnoiseFiltered_nChannel);
      fEventTree->Branch("nf_channelId",&fTPCnoiseFiltered_channelId);
      fTPCnoiseFiltered_wf = new TClonesArray("TH1S");
      fEventTree->Branch("nf_wf", &fTPCnoiseFiltered_wf, 256000, 0);
      fEventTree->Branch("nf_scale", &fScaleRawdigit);
      fEventTree->Branch("nf_shift",&fShiftRawdigit);
    }
    if(fSaveTPC_deconWiener){
    fEventTree->Branch("calibWiener_nChannel", &fTPCdeconWiener_nChannel);
    fEventTree->Branch("calibWiener_channelId", &fTPCdeconWiener_channelId);
    fTPCdeconWiener_wf = new TClonesArray("TH1F");
    fEventTree->Branch("calibWiener_wf", &fTPCdeconWiener_wf, 256000, 0);
    }
    if(fSaveTPC_deconGaussian){
      fEventTree->Branch("calibGaussian_nChannel", &fTPCdeconGaussian_nChannel);
      fEventTree->Branch("calibGaussian_channelId", &fTPCdeconGaussian_channelId);
      fTPCdeconGaussian_wf = new TClonesArray("TH1F");
      fEventTree->Branch("calibGaussian_wf", &fTPCdeconGaussian_wf, 256000, 0);
    }
    if(fSaveTPC_badChannelList){
      fEventTree->Branch("badChannel", &fBadChannel);
      fEventTree->Branch("badBegin", &fBadBegin);
      fEventTree->Branch("badEnd", &fBadEnd);
    }
    if(fSaveTPC_channelThreshold){
      fEventTree->Branch("channelThreshold", &fChannelThreshold);
    }
    if(fSavePMT_raw){
      fEventTree->Branch("cosmic_hg_opch", &fOp_cosmic_hg_opch);
      fEventTree->Branch("cosmic_lg_opch", &fOp_cosmic_lg_opch);
      fEventTree->Branch("beam_hg_opch", &fOp_beam_hg_opch);
      fEventTree->Branch("beam_lg_opch", &fOp_beam_lg_opch);
      fEventTree->Branch("cosmic_hg_timestamp", &fOp_cosmic_hg_timestamp);
      fEventTree->Branch("cosmic_lg_timestamp", &fOp_cosmic_lg_timestamp);
      fEventTree->Branch("beam_hg_timestamp", &fOp_beam_hg_timestamp);
      fEventTree->Branch("beam_lg_timestamp", &fOp_beam_lg_timestamp);
      fOp_cosmic_hg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("cosmic_hg_wf", &fOp_cosmic_hg_wf, 256000, 0);
      fOp_cosmic_lg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("cosmic_lg_wf", &fOp_cosmic_lg_wf, 256000, 0);
      fOp_beam_hg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("beam_hg_wf", &fOp_beam_hg_wf, 256000, 0);
      fOp_beam_lg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("beam_lg_wf", &fOp_beam_lg_wf, 256000, 0);
    }
    if(fSavePMT_overlay){
      fEventTree->Branch("detsim_cosmic_hg_opch", &fOp_detsim_cosmic_hg_opch);
      fEventTree->Branch("detsim_cosmic_lg_opch", &fOp_detsim_cosmic_lg_opch);
      fEventTree->Branch("detsim_beam_hg_opch", &fOp_detsim_beam_hg_opch);
      fEventTree->Branch("detsim_beam_lg_opch", &fOp_detsim_beam_lg_opch);
      fEventTree->Branch("detsim_cosmic_hg_timestamp", &fOp_detsim_cosmic_hg_timestamp);
      fEventTree->Branch("detsim_cosmic_lg_timestamp", &fOp_detsim_cosmic_lg_timestamp);
      fEventTree->Branch("detsim_beam_hg_timestamp", &fOp_detsim_beam_hg_timestamp);
      fEventTree->Branch("detsim_beam_lg_timestamp", &fOp_detsim_beam_lg_timestamp);
      fOp_detsim_cosmic_hg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("detsim_cosmic_hg_wf", &fOp_detsim_cosmic_hg_wf, 256000, 0);
      fOp_detsim_cosmic_lg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("detsim_cosmic_lg_wf", &fOp_detsim_cosmic_lg_wf, 256000, 0);
      fOp_detsim_beam_hg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("detsim_beam_hg_wf", &fOp_detsim_beam_hg_wf, 256000, 0);
      fOp_detsim_beam_lg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("detsim_beam_lg_wf", &fOp_detsim_beam_lg_wf, 256000, 0);
      fEventTree->Branch("mixer_beam_hg_opch", &fOp_mixer_beam_hg_opch);
      fEventTree->Branch("mixer_beam_lg_opch", &fOp_mixer_beam_lg_opch);
      fEventTree->Branch("mixer_beam_hg_timestamp", &fOp_mixer_beam_hg_timestamp);
      fEventTree->Branch("mixer_beam_lg_timestamp", &fOp_mixer_beam_lg_timestamp);
      fOp_mixer_beam_hg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("mixer_beam_hg_wf", &fOp_mixer_beam_hg_wf, 256000, 0);
      fOp_mixer_beam_lg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("mixer_beam_lg_wf", &fOp_mixer_beam_lg_wf, 256000, 0);
      fEventTree->Branch("saturation_beam_hg_opch", &fOp_saturation_beam_hg_opch);
      fEventTree->Branch("saturation_beam_lg_opch", &fOp_saturation_beam_lg_opch);
      fEventTree->Branch("saturation_beam_hg_timestamp", &fOp_saturation_beam_hg_timestamp);
      fEventTree->Branch("saturation_beam_lg_timestamp", &fOp_saturation_beam_lg_timestamp);
      fOp_saturation_beam_hg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("saturation_beam_hg_wf", &fOp_saturation_beam_hg_wf, 256000, 0);
      fOp_saturation_beam_lg_wf = new TClonesArray("TH1S");
      fEventTree->Branch("saturation_beam_lg_wf", &fOp_saturation_beam_lg_wf, 256000, 0);
    }
    if(fSavePMT_wfmSaturation){
      fOp_wf = new TClonesArray("TH1S");
      fEventTree->Branch("op_wf", &fOp_wf, 256000, 0);
      fEventTree->Branch("op_femch", &fOp_femch);
      fEventTree->Branch("op_timestamp", &fOp_timestamp);
    }
    if(fSavePMT_reco){
      fEventTree->Branch("oh_nHits", &oh_nHits);
      fEventTree->Branch("oh_channel", &oh_channel);
      fEventTree->Branch("oh_bgtime", &oh_bgtime);
      fEventTree->Branch("oh_trigtime", &oh_trigtime); 
      fEventTree->Branch("oh_pe", &oh_pe); 
      fEventTree->Branch("of_nFlash", &of_nFlash);
      fEventTree->Branch("of_t", &of_t);
      fEventTree->Branch("of_peTotal", &of_peTotal); 
      fEventTree->Branch("of_multiplicity", &of_multiplicity);
      fPEperOpDet = new TClonesArray("TH1F");
      fEventTree->Branch("pe_opdet", &fPEperOpDet, 256000, 0);
    }
    if(fSavePMT_gainData){
      fEventTree->Branch("op_gain", &fOp_gain);
      fEventTree->Branch("op_gainerror", &fOp_gainerror);
    }
    if(fSavePMT_gainMC){
      fEventTree->Branch("op_gain", &fOp_gain);
      fEventTree->Branch("op_gainerror", &fOp_gainerror);
    }
    if(fSaveHWTrigger){    
      fEventTree->Branch("triggerNo", &fTriggernumber);
      fEventTree->Branch("triggerTime", &fTriggertime);
      fEventTree->Branch("beamgateTime", &fBeamgatetime);
      fEventTree->Branch("triggerBits", &fTriggerbits);
    }
    if(fSaveSWTrigger){
      fEventTree->Branch("passAlgo", &fpassAlgo);
      fEventTree->Branch("passPrescale", &fpassPrescale);
      fEventTree->Branch("PHMAX", &fPHMAX);
      fEventTree->Branch("multiplicity", &fmultiplicity);
      fEventTree->Branch("triggerTick", &ftriggerTick);
      fEventTree->Branch("swtriggerTime", &ftriggerTime);
      fEventTree->Branch("prescale_weight", &fprescale_weight);
    }
    if(fSaveBeam){
      fEventTree->Branch("tor101", &ftor101);
      fEventTree->Branch("tortgt", &ftortgt);
      fEventTree->Branch("trtgtd", &ftrtgtd);
      fEventTree->Branch("t_ms", &ft_ms);
      fEventTree->Branch("recordType", &fRecordType);
      fEventTree->Branch("seconds", &fSeconds);
      fEventTree->Branch("milliseconds", &fMilliSeconds);
      fEventTree->Branch("numDevices", &fNumberOfDevices);
      fEventTree->Branch("dataMap", &fDataMap);
    }
    if(fSavePOT){
      fEventTree->Branch("totpot_bnbETOR875", &ftotpot_bnbETOR875);
      fEventTree->Branch("totpot_numiETORTGT", &ftotpot_numiETORTGT);
      fEventTree->Branch("totpot_bnbETOR860", &ftotpot_bnbETOR860);
      fEventTree->Branch("totgoodpot_bnbETOR875", &ftotgoodpot_bnbETOR875);
      fEventTree->Branch("totgoodpot_numiETORTGT", &ftotgoodpot_numiETORTGT);
      fEventTree->Branch("totgoodpot_bnbETOR860", &ftotgoodpot_bnbETOR860);
      fEventTree->Branch("totspills_bnbETOR875", &ftotspills_bnbETOR875);
      fEventTree->Branch("totspills_numiETORTGT", &ftotspills_numiETORTGT);
      fEventTree->Branch("totspills_bnbETOR860", &ftotspills_bnbETOR860);
      fEventTree->Branch("goodspills_bnbETOR875", &fgoodspills_bnbETOR875);
      fEventTree->Branch("goodspills_numiETORTGT", &fgoodspills_numiETORTGT);
      fEventTree->Branch("goodspills_bnbETOR860", &fgoodspills_bnbETOR860);
    }
    if(fSaveMuCSData){
      fEventTree->Branch("t0",&ft0);
    }
    if(fSaveMuCSRecoData){
      fEventTree->Branch("theta_xy",&ftheta_xy);
      fEventTree->Branch("theta_xy_rms",&ftheta_xy_rms);
      fEventTree->Branch("x",&fx);
      fEventTree->Branch("x_rms",&fx_rms);
      fEventTree->Branch("theta_yz",&ftheta_yz);
      fEventTree->Branch("theta_yz_rms",&ftheta_yz_rms);
      fEventTree->Branch("z",&fz);
      fEventTree->Branch("z_rms",&fz_rms);
      fEventTree->Branch("y",&fy);
      fEventTree->Branch("theta",&ftheta);
      fEventTree->Branch("phi",&fphi);
      fEventTree->Branch("xmatches",&fxmatches);
      fEventTree->Branch("zmatches",&fzmatches);
    }
    if(fSaveMCParticle){
      fEventTree->Branch("mc_Ntrack", &mc_Ntrack);
      fEventTree->Branch("mc_id", &mc_id, "mc_id[mc_Ntrack]/I");      
      fEventTree->Branch("mc_pdg", &mc_pdg, "mc_pdg[mc_Ntrack]/I");      
      fEventTree->Branch("mc_process", &mc_process, "mc_process[mc_Ntrack]/I");
      fEventTree->Branch("mc_mother", &mc_mother, "mc_mother[mc_Ntrack]/I");  
      fEventTree->Branch("mc_daughters", &mc_daughters);   
      fEventTree->Branch("mc_startXYZT", &mc_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");      
      fEventTree->Branch("mc_endXYZT", &mc_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");      
      fEventTree->Branch("mc_startMomentum", &mc_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");      
      fEventTree->Branch("mc_endMomentum", &mc_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");      
      fMC_trackPosition = new TObjArray();
      fMC_trackPosition->SetOwner(kTRUE);
      fEventTree->Branch("mc_trackPosition", &fMC_trackPosition);      
      fEventTree->Branch("mc_isnu", &mc_isnu);      
      fEventTree->Branch("mc_nGeniePrimaries", &mc_nGeniePrimaries);      
      fEventTree->Branch("mc_nu_pdg", &mc_nu_pdg);      
      fEventTree->Branch("mc_nu_ccnc", &mc_nu_ccnc);      
      fEventTree->Branch("mc_nu_mode", &mc_nu_mode);      
      fEventTree->Branch("mc_nu_intType", &mc_nu_intType);      
      fEventTree->Branch("mc_nu_target", &mc_nu_target);      
      fEventTree->Branch("mc_hitnuc", &mc_hitnuc);      
      fEventTree->Branch("mc_hitquark", &mc_hitquark);      
      fEventTree->Branch("mc_nu_Q2", &mc_nu_Q2);      
      fEventTree->Branch("mc_nu_W", &mc_nu_W);      
      fEventTree->Branch("mc_nu_X", &mc_nu_X);      
      fEventTree->Branch("mc_nu_Y", &mc_nu_Y);      
      fEventTree->Branch("mc_nu_Pt", &mc_nu_Pt);      
      fEventTree->Branch("mc_nu_Theta", &mc_nu_Theta);      
      fEventTree->Branch("mc_nu_pos", &mc_nu_pos, "mc_nu_pos[4]/F");      
      fEventTree->Branch("mc_nu_mom", &mc_nu_mom, "mc_nu_mom[4]/F");      
    }
    if(fSaveSimEnergyDepositInside){
      fEventTree->Branch("sedi_nelectrons",&sedi_nelectrons);
      fEventTree->Branch("sedi_nphotons",&sedi_nphotons);
      fEventTree->Branch("sedi_time_start",&sedi_time_start);
      fEventTree->Branch("sedi_time_end",&sedi_time_end);
      fEventTree->Branch("sedi_x_start",&sedi_x_start);
      fEventTree->Branch("sedi_x_end",&sedi_x_end);
      fEventTree->Branch("sedi_y_start",&sedi_y_start);
      fEventTree->Branch("sedi_y_end",&sedi_y_end);
      fEventTree->Branch("sedi_z_start",&sedi_z_start);
      fEventTree->Branch("sedi_z_end",&sedi_z_end);
      fEventTree->Branch("sedi_pdg",&sedi_pdg);
      fEventTree->Branch("sedi_trackId",&sedi_trackId);
      fEventTree->Branch("sedi_energy",&sedi_energy);
    }
    if(fSaveSimEnergyDepositOutside){
      fEventTree->Branch("sedo_nelectrons",&sedo_nelectrons);
      fEventTree->Branch("sedo_nphotons",&sedo_nphotons);
      fEventTree->Branch("sedo_time_start",&sedo_time_start);
      fEventTree->Branch("sedo_time_end",&sedo_time_end);
      fEventTree->Branch("sedo_x_start",&sedo_x_start);
      fEventTree->Branch("sedo_x_end",&sedo_x_end);
      fEventTree->Branch("sedo_y_start",&sedo_y_start);
      fEventTree->Branch("sedo_y_end",&sedo_y_end);
      fEventTree->Branch("sedo_z_start",&sedo_z_start);
      fEventTree->Branch("sedo_z_end",&sedo_z_end);
      fEventTree->Branch("sedo_pdg",&sedo_pdg);
      fEventTree->Branch("sedo_trackId",&sedo_trackId);
      fEventTree->Branch("sedo_energy",&sedo_energy);
    }
    if(fSaveSimChannel){
      fEventTree->Branch("sc_x",&sc_x);
      fEventTree->Branch("sc_y",&sc_y);
      fEventTree->Branch("sc_z",&sc_z);
      fEventTree->Branch("sc_q",&sc_q);
      fEventTree->Branch("sc_chan",&sc_chan);
      fEventTree->Branch("sc_tdc",&sc_tdc);
      fEventTree->Branch("sc_tick",&sc_tick);
      fEventTree->Branch("sc_drift_speed",&sc_drift_speed);
    }
  }

  void CellTreeUB::beginJob(){
  }

  void CellTreeUB::endJob(){
  }

  void CellTreeUB::beginRun(const art::Run& ){
    mf::LogInfo("CellTreeUB") << "begin run";
  }

  void CellTreeUB::beginSubRun(const art::SubRun& subrun){
    mf::LogInfo("CellTreeUB") << "begin sub run";
  }

  void CellTreeUB::endSubRun(const art::SubRun& subrun){
    if(fSavePOT) processPOT(subrun); 
  }

  void CellTreeUB::analyze(const art::Event& evt){
    reset();
    fEvent = evt.id().event();
    fRun = evt.run();
    fSubRun = evt.subRun();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    fEventTime = tts.AsDouble();
    if(fSaveTPC_raw) processTPC_raw(evt);
    if(fSaveTPC_noiseFiltered) processTPC_noiseFiltered(evt);
    if(fSaveTPC_deconWiener) processTPC_deconWiener(evt);
    if(fSaveTPC_deconGaussian) processTPC_deconGaussian(evt);
    if(fSaveTPC_badChannelList) processTPC_badChannelList(evt);
    if(fSaveTPC_channelThreshold) processTPC_channelThreshold(evt);
    if(fSavePMT_raw) processPMT_raw(evt);
    if(fSavePMT_overlay) processPMT_overlay(evt);
    if(fSavePMT_wfmSaturation) processPMT_wfmSaturation(evt);
    if(fSavePMT_reco) processPMT_reco(evt);
    if(fSavePMT_gainData) processPMT_gainData(evt);
    if(fSavePMT_gainMC) processPMT_gainMC(evt);
    if(fSaveHWTrigger) processHWTrigger(evt);
    if(fSaveSWTrigger) processSWTrigger(evt);
    if(fSaveBeam) processBeam(evt);
    if(fSaveMuCSData) processMuCS_data(evt);
    if(fSaveMuCSRecoData) processMuCS_recoData(evt);
    if(fSaveMCParticle) processMCParticle(evt);
    if(fSaveSimEnergyDepositInside) processSimEnergyDepositInside(evt);
    if(fSaveSimEnergyDepositOutside) processSimEnergyDepositOutside(evt);
    if(fSaveSimChannel) processSimChannel(evt);
    fEventTree->Fill();
  }

  void CellTreeUB::reset(){
    if(fSaveTPC_raw){
      fTPCraw_channelId.clear();
      fTPCraw_wf->Delete();
    }
    if(fSaveTPC_noiseFiltered){
      fTPCnoiseFiltered_channelId.clear();
      fTPCnoiseFiltered_wf->Delete();
    }
    if(fSaveTPC_deconWiener){
      fTPCdeconWiener_channelId.clear();
      fTPCdeconWiener_wf->Delete();
    }
    if(fSaveTPC_deconGaussian){
      fTPCdeconGaussian_channelId.clear();
      fTPCdeconGaussian_wf->Delete();
    }
    if(fSaveTPC_badChannelList){
      fBadChannel.clear();
      fBadBegin.clear();
      fBadEnd.clear();
    }
    if(fSaveTPC_channelThreshold){
      fChannelThreshold.clear();
    }
    if(fSavePMT_raw){
      fOp_cosmic_hg_opch.clear();
      fOp_cosmic_lg_opch.clear();
      fOp_beam_hg_opch.clear();
      fOp_beam_lg_opch.clear();
      fOp_cosmic_hg_timestamp.clear();
      fOp_cosmic_lg_timestamp.clear();
      fOp_beam_hg_timestamp.clear();
      fOp_beam_lg_timestamp.clear();
      fOp_cosmic_hg_wf->Delete();
      fOp_cosmic_lg_wf->Delete();
      fOp_beam_hg_wf->Delete();
      fOp_beam_lg_wf->Delete();
    }
    if(fSavePMT_overlay){
      fOp_detsim_cosmic_hg_opch.clear();
      fOp_detsim_cosmic_lg_opch.clear();
      fOp_detsim_beam_hg_opch.clear();
      fOp_detsim_beam_lg_opch.clear();
      fOp_detsim_cosmic_hg_timestamp.clear();
      fOp_detsim_cosmic_lg_timestamp.clear();
      fOp_detsim_beam_hg_timestamp.clear();
      fOp_detsim_beam_lg_timestamp.clear();
      fOp_detsim_cosmic_hg_wf->Delete();
      fOp_detsim_cosmic_lg_wf->Delete();
      fOp_detsim_beam_hg_wf->Delete();
      fOp_detsim_beam_lg_wf->Delete();
      fOp_mixer_beam_hg_opch.clear();
      fOp_mixer_beam_lg_opch.clear();
      fOp_mixer_beam_hg_timestamp.clear();
      fOp_mixer_beam_lg_timestamp.clear();
      fOp_mixer_beam_hg_wf->Delete();
      fOp_mixer_beam_lg_wf->Delete();
      fOp_saturation_beam_hg_opch.clear();
      fOp_saturation_beam_lg_opch.clear();
      fOp_saturation_beam_hg_timestamp.clear();
      fOp_saturation_beam_lg_timestamp.clear();
      fOp_saturation_beam_hg_wf->Delete();
      fOp_saturation_beam_lg_wf->Delete();
    }
    if(fSavePMT_wfmSaturation){
      fOp_wf->Delete();
      fOp_femch.clear();
      fOp_timestamp.clear();
    }
    if(fSavePMT_reco){
      oh_channel.clear();
      oh_bgtime.clear();
      oh_trigtime.clear();
      oh_pe.clear();
      of_t.clear();
      of_peTotal.clear();
      of_multiplicity.clear();
      fPEperOpDet->Delete();
    }
    if(fSavePMT_gainData){
      fOp_gain.clear();
      fOp_gainerror.clear();
    }
    if(fSavePMT_gainMC){
      fOp_gain.clear();
      fOp_gainerror.clear();
    }
    if(fSaveSWTrigger){
      fpassAlgo.clear();
      fpassPrescale.clear();
      fPHMAX.clear();
      fmultiplicity.clear();
      ftriggerTick.clear();
      ftriggerTime.clear();
      fprescale_weight.clear();
    }
    if(fSaveMCParticle){
      mc_Ntrack = 0;
      for(int i=0; i<MAX_TRACKS; i++){
	mc_id[i] = 0;
	mc_pdg[i] = 0;
	mc_mother[i] = 0;
	for(int j=0; j<4; j++){
	  mc_startXYZT[i][j] = 0;
	  mc_endXYZT[i][j] = 0;
	  mc_startMomentum[i][j] = 0;
	  mc_endMomentum[i][j] = 0;
	}
      }
      mc_daughters.clear();
      savedMCTrackIdMap.clear();
      fMC_trackPosition->Delete();
      mc_isnu = 0;
      mc_nGeniePrimaries = -1;
      mc_nu_pdg = -1;
      mc_nu_ccnc = -1;
      mc_nu_mode = -1;
      mc_nu_intType = -1;
      mc_nu_target = -1;
      mc_hitnuc = -1;
      mc_hitquark = -1;
      mc_nu_Q2 = -1;
      mc_nu_W = -1;
      mc_nu_X = -1;
      mc_nu_Y = -1;
      mc_nu_Pt = -1;
      mc_nu_Theta = -1;
      for(int i=0; i<4; i++){
	mc_nu_pos[i] = 0;
	mc_nu_mom[i] = 0;
      }
    }
    if(fSaveSimEnergyDepositInside){
      sedi_nelectrons.clear();
      sedi_nphotons.clear();
      sedi_time_start.clear();
      sedi_time_end.clear();
      sedi_x_start.clear();
      sedi_x_end.clear();
      sedi_y_start.clear();
      sedi_y_end.clear();
      sedi_z_start.clear();
      sedi_z_end.clear();
      sedi_pdg.clear();
      sedi_trackId.clear();
      sedi_energy.clear();
    }
    if(fSaveSimEnergyDepositOutside){
      sedo_nelectrons.clear();
      sedo_nphotons.clear();
      sedo_time_start.clear();
      sedo_time_end.clear();
      sedo_x_start.clear();
      sedo_x_end.clear();
      sedo_y_start.clear();
      sedo_y_end.clear();
      sedo_z_start.clear();
      sedo_z_end.clear();
      sedo_pdg.clear();
      sedo_trackId.clear();
      sedo_energy.clear();
    }
    if(fSaveSimChannel){
      sc_x.clear();
      sc_y.clear();
      sc_z.clear();
      sc_tdc.clear();
    }
  }

  void CellTreeUB::InitProcessMap(){
    processMap["unknown"]                   = 0;
    processMap["primary"]                   = 1;
    processMap["compt"]                     = 2;
    processMap["phot"]                      = 3;
    processMap["annihil"]                   = 4;
    processMap["eIoni"]                     = 5;
    processMap["eBrem"]                     = 6;
    processMap["conv"]                      = 7;
    processMap["muIoni"]                    = 8;
    processMap["muMinusCaptureAtRest"]      = 9;
    processMap["NeutronInelastic"]          = 10;
    processMap["nCapture"]                  = 11;
    processMap["hadElastic"]                = 12;
    processMap["Decay"]                     = 13;
    processMap["CoulombScat"]               = 14;
    processMap["muPairProd"]                = 15;
    processMap["muBrems"]                   = 16;
    processMap["muPairProd"]                = 17;
    processMap["PhotonInelastic"]           = 18;
    processMap["hIoni"]                     = 19;
    processMap["ProtonInelastic"]           = 20;
    processMap["PionPlusInelastic"]         = 21;
    processMap["CHIPSNuclearCaptureAtRest"] = 22;
    processMap["PionMinusInelastic"]        = 23;
  }

  void CellTreeUB::processTPC_raw(const art::Event& evt){
    art::Handle<std::vector<raw::RawDigit> > rawdigit;
    if(! evt.getByLabel(fTPC_rawProducer, fTPC_rawLabel,  rawdigit)){
      cout << "WARNING: no raw::RawDigit producer" << fTPC_rawProducer 
	   << " or label " << fTPC_rawLabel << endl;
      return;
    }
    std::vector<art::Ptr<raw::RawDigit> > rd_v;
    art::fill_ptr_vector(rd_v, rawdigit);
    fTPCraw_nChannel = rd_v.size();
    int i=0;
    for(auto const& rd : rd_v){
      int refChan = rd->Channel();
      if(saveVYshorted == true && (refChan<3566 || refChan>4305)) continue;
      fTPCraw_channelId.push_back(rd->Channel());
      int nSamples = rd->Samples();
      std::vector<short> uncompressed(nSamples);
      raw::Uncompress(rd->ADCs(), uncompressed, rd->Compression());
      TH1S *h = new((*fTPCraw_wf)[i]) TH1S("","",9600,0,9600);
      for(int j=1; j<=nSamples; j++){ h->SetBinContent(j,uncompressed[j-1]); }
      i++;
    } 
  }

  void CellTreeUB::processTPC_noiseFiltered(const art::Event& evt){ // for old pre-MCC9 L1SP within imaging
    art::Handle<std::vector<recob::Wire> > wire;
    if(! evt.getByLabel(fTPC_noiseFilteredProducer, fTPC_noiseFilteredLabel, wire)){
      cout << "WARNING: no recob::Wire producer " << fTPC_noiseFilteredProducer 
	   << " or label " << fTPC_noiseFilteredLabel << endl;
      return;
    }
    fScaleRawdigit = scaleRawdigit;
    fShiftRawdigit = shiftRawdigit;

    std::vector<art::Ptr<recob::Wire> > w_v;
    art::fill_ptr_vector(w_v, wire);
    fTPCnoiseFiltered_nChannel = w_v.size();
    int i=0;
    for(auto const& w : w_v){
      int refChan = w->Channel();
      if(saveVYshorted == true && (refChan<3566 || refChan>4305)) continue;
      fTPCnoiseFiltered_channelId.push_back(w->Channel());
      std::vector<float> wf = w->Signal();
      int nbin = (int)wf.size();
      TH1S *h = new((*fTPCnoiseFiltered_wf)[i]) TH1S("","",nbin,0,nbin);
      for(int j=1; j<=nbin; j++){ 
	short temp = short(wf[j]*fScaleRawdigit+0.5-fShiftRawdigit);
	h->SetBinContent(j, temp); 
      }
      i++;
    }
  }
  

  void CellTreeUB::processTPC_deconWiener(const art::Event& evt){
    art::Handle<std::vector<recob::Wire> > wire;
    if(! evt.getByLabel(fTPC_deconWienerProducer, fTPC_deconWienerLabel, wire)){
      cout << "WARNING: no recob::Wire producer " << fTPC_deconWienerProducer 
	   << " or label " << fTPC_deconWienerLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::Wire> > w_v;
    art::fill_ptr_vector(w_v, wire);
    fTPCdeconWiener_nChannel = w_v.size();
    int i=0;
    for(auto const& w : w_v){
      fTPCdeconWiener_channelId.push_back(w->Channel());
      std::vector<float> wf = w->Signal();
      int nbin = (int)wf.size()/deconRebin;
      TH1F *h = new((*fTPCdeconWiener_wf)[i]) TH1F("","",nbin,0,nbin);
      int ctr = 0, binCtr = 1;
      float content = 0;
      int nSamples = (int)wf.size();
      for(int j=0; j<nSamples; j++){ 
	ctr++;
	content += wf[j];
	if(ctr %  deconRebin == 0){
	  h->SetBinContent(binCtr, content); 
	  content = 0;
	  ctr = 0;
	  binCtr++;
	}
      }
      i++;
    }
  }

  void CellTreeUB::processTPC_deconGaussian(const art::Event& evt){
    art::Handle<std::vector<recob::Wire> > wire;
    if(! evt.getByLabel(fTPC_deconGaussianProducer, fTPC_deconGaussianLabel, wire)){
      cout << "WARNING: no recob::Wire producer " << fTPC_deconGaussianProducer 
	   << " or no label " << fTPC_deconGaussianLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::Wire> > w_v;
    art::fill_ptr_vector(w_v, wire);
    fTPCdeconGaussian_nChannel = w_v.size();
    int i=0;
    for(auto const& w : w_v){
      fTPCdeconGaussian_channelId.push_back(w->Channel());
      std::vector<float> wf = w->Signal();
      int nbin = (int)wf.size()/deconRebin;
      TH1F *h = new((*fTPCdeconGaussian_wf)[i]) TH1F("","",nbin,0,nbin);
      int ctr = 0, binCtr = 1;
      float content = 0;
      int nSamples = (int)wf.size();
      for(int j=0; j<nSamples; j++){ 
	ctr++;
	content += wf[j];
	if(ctr % deconRebin == 0){
	  h->SetBinContent(binCtr, content);
	  content = 0;
	  ctr = 0;
	  binCtr++;
	} 
      }
      i++;
    }
  }

  void CellTreeUB::processTPC_badChannelList(const art::Event& evt){
    art::Handle<std::vector<int> > bad;
    if(! evt.getByLabel(fTPC_badChannelListProducer, fTPC_badChannelListLabel, bad)){
      cout << "WARNING: no bad channel list with producer " << fTPC_badChannelListProducer 
	   << " or label " << fTPC_badChannelListLabel << endl;
      return;
    }
    vector<int> const& bad_v(*bad);
    int nch = (int)bad_v.size()/3;
    for(int i=0; i<nch; i++){
      const int offset = 3*i;
      fBadChannel.push_back(bad_v[offset+0]);
      fBadBegin.push_back(bad_v[offset+1]);
      fBadEnd.push_back(bad_v[offset+2]);
    } 
  }

  void CellTreeUB::processTPC_channelThreshold(const art::Event& evt){
    art::Handle<std::vector<double> > threshold;
    if(! evt.getByLabel(fTPC_channelThresholdProducer, fTPC_channelThresholdLabel, threshold)){
      cout << "WARNING: no channel thresholds with producer " << fTPC_channelThresholdProducer 
	   << " or label " << fTPC_channelThresholdLabel << endl;
      return;
    }
    vector<double> const& threshold_v(*threshold);
    fChannelThreshold = threshold_v;
  }

  void CellTreeUB::processPMT_raw(const art::Event& evt){
    art::Handle<std::vector<raw::OpDetWaveform> > opwf_cosmic_hg;
    evt.getByLabel(fPMT_HG_cosmicProducer, fPMT_HG_cosmicLabel, opwf_cosmic_hg);
    if(opwf_cosmic_hg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*opwf_cosmic_hg);
      int i=0;
      for(auto const& op : opwf_v){
	fOp_cosmic_hg_opch.push_back(op.ChannelNumber());
	fOp_cosmic_hg_timestamp.push_back(op.TimeStamp());
	int nbins = (int)op.size();
	TH1S *h = new ((*fOp_cosmic_hg_wf)[i]) TH1S("","",nbins,0,nbins);
	for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
	i++;
      }
    }
    if(!opwf_cosmic_hg.isValid()){
      cout << "WARNINg: no raw::OpDetWaveform producer " << fPMT_HG_cosmicProducer
	   << " or no label " << fPMT_HG_cosmicLabel << endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > opwf_beam_hg;
    evt.getByLabel(fPMT_HG_beamProducer, fPMT_HG_beamLabel, opwf_beam_hg);
    if(opwf_beam_hg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*opwf_beam_hg);
      int i=0;
      for(auto const& op : opwf_v){
	fOp_beam_hg_opch.push_back(op.ChannelNumber());
	fOp_beam_hg_timestamp.push_back(op.TimeStamp());
	int nbins = (int)op.size();
	TH1S *h = new ((*fOp_beam_hg_wf)[i]) TH1S("","",nbins,0,nbins);
	for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
	i++;
      }
    }
    if(!opwf_beam_hg.isValid()){
      cout << "WARNINg: no raw::OpDetWaveform producer " << fPMT_HG_beamProducer
	   << " or no label " << fPMT_HG_beamLabel << endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > opwf_cosmic_lg;
    evt.getByLabel( (_use_LG_beam_for_HG_cosmic ? fPMT_LG_beamProducer : fPMT_LG_cosmicProducer),
		    (_use_LG_beam_for_HG_cosmic ? fPMT_LG_beamLabel    : fPMT_LG_cosmicLabel   ),
		    opwf_cosmic_lg);
    if(opwf_cosmic_lg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*opwf_cosmic_lg);
      int i=0;
      for(auto const& op : opwf_v){
	fOp_cosmic_lg_opch.push_back(op.ChannelNumber());
	fOp_cosmic_lg_timestamp.push_back(op.TimeStamp());
	int nbins = (int)op.size();
	TH1S *h = new ((*fOp_cosmic_lg_wf)[i]) TH1S("","",nbins,0,nbins);
	for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
	i++;
      }
    }
    if(!opwf_cosmic_lg.isValid()){
      cout << "WARNINg: no raw::OpDetWaveform producer " << fPMT_LG_cosmicProducer
	   << " or no label " << fPMT_LG_cosmicLabel << endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > opwf_beam_lg;
    evt.getByLabel(fPMT_LG_beamProducer, fPMT_LG_beamLabel, opwf_beam_lg);
    if(opwf_beam_lg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*opwf_beam_lg);
      int i=0;
      for(auto const& op : opwf_v){
	fOp_beam_lg_opch.push_back(op.ChannelNumber());
	fOp_beam_lg_timestamp.push_back(op.TimeStamp());
	int nbins = (int)op.size();
	TH1S *h = new ((*fOp_beam_lg_wf)[i]) TH1S("","",nbins,0,nbins);
	for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
	i++;
      }
    }
    if(!opwf_beam_lg.isValid()){
      cout << "WARNINg: no raw::OpDetWaveform producer " << fPMT_LG_beamProducer
	   << " or no label " << fPMT_LG_beamLabel << endl;
    }
  }

  void CellTreeUB::processPMT_overlay(const art::Event& evt){
    art::Handle<std::vector<raw::OpDetWaveform> > detsim_opwf_cosmic_hg;
    evt.getByLabel(fPMT_overlay_detsim_producer, "OpdetCosmicHighGain", detsim_opwf_cosmic_hg);
    if(detsim_opwf_cosmic_hg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*detsim_opwf_cosmic_hg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_detsim_cosmic_hg_opch.push_back(op.ChannelNumber());
        fOp_detsim_cosmic_hg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_detsim_cosmic_hg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!detsim_opwf_cosmic_hg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_detsim_producer
           << " or no label OpdetCosmicHighGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > detsim_opwf_cosmic_lg;
    evt.getByLabel(fPMT_overlay_detsim_producer, "OpdetCosmicLowGain", detsim_opwf_cosmic_lg);
    if(detsim_opwf_cosmic_lg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*detsim_opwf_cosmic_lg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_detsim_cosmic_lg_opch.push_back(op.ChannelNumber());
        fOp_detsim_cosmic_lg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_detsim_cosmic_lg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!detsim_opwf_cosmic_lg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_detsim_producer
           << " or no label OpdetCosmicLowGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > detsim_opwf_beam_hg;
    evt.getByLabel(fPMT_overlay_detsim_producer, "OpdetBeamHighGain", detsim_opwf_beam_hg);
    if(detsim_opwf_beam_hg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*detsim_opwf_beam_hg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_detsim_beam_hg_opch.push_back(op.ChannelNumber());
        fOp_detsim_beam_hg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_detsim_beam_hg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!detsim_opwf_beam_hg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_detsim_producer
           << " or no label OpdetBeamHighGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > detsim_opwf_beam_lg;
    evt.getByLabel(fPMT_overlay_detsim_producer, "OpdetBeamLowGain", detsim_opwf_beam_lg);
    if(detsim_opwf_beam_lg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*detsim_opwf_beam_lg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_detsim_beam_lg_opch.push_back(op.ChannelNumber());
        fOp_detsim_beam_lg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_detsim_beam_lg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!detsim_opwf_beam_lg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_detsim_producer
           << " or no label OpdetBeamLowGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > mixer_opwf_beam_hg;
    evt.getByLabel(fPMT_overlay_mixer_producer, "OpdetBeamHighGain", mixer_opwf_beam_hg);
    if(mixer_opwf_beam_hg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*mixer_opwf_beam_hg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_mixer_beam_hg_opch.push_back(op.ChannelNumber());
        fOp_mixer_beam_hg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_mixer_beam_hg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!mixer_opwf_beam_hg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_mixer_producer
           << " or no label OpdetBeamHighGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > mixer_opwf_beam_lg;
    evt.getByLabel(fPMT_overlay_mixer_producer, "OpdetBeamLowGain", mixer_opwf_beam_lg);
    if(mixer_opwf_beam_lg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*mixer_opwf_beam_lg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_mixer_beam_lg_opch.push_back(op.ChannelNumber());
        fOp_mixer_beam_lg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_mixer_beam_lg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!mixer_opwf_beam_lg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_mixer_producer
           << " or no label OpdetBeamLowGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > saturation_opwf_beam_hg;
    evt.getByLabel(fPMT_overlay_saturation_producer, "OpdetBeamHighGain", saturation_opwf_beam_hg);
    if(saturation_opwf_beam_hg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*saturation_opwf_beam_hg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_saturation_beam_hg_opch.push_back(op.ChannelNumber());
        fOp_saturation_beam_hg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_saturation_beam_hg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!saturation_opwf_beam_hg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_saturation_producer
           << " or no label OpdetBeamHighGain"<< endl;
    }
    art::Handle<std::vector<raw::OpDetWaveform> > saturation_opwf_beam_lg;
    evt.getByLabel(fPMT_overlay_saturation_producer, "OpdetBeamLowGain", saturation_opwf_beam_lg);
    if(saturation_opwf_beam_lg.isValid()){
      std::vector<raw::OpDetWaveform> const& opwf_v(*saturation_opwf_beam_lg);
      int i=0;
      for(auto const& op : opwf_v){
        fOp_saturation_beam_lg_opch.push_back(op.ChannelNumber());
        fOp_saturation_beam_lg_timestamp.push_back(op.TimeStamp());
        int nbins = (int)op.size();
        TH1S *h = new ((*fOp_saturation_beam_lg_wf)[i]) TH1S("","",nbins,0,nbins);
        for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j-1]);}
        i++;
      }
    }
    if(!saturation_opwf_beam_lg.isValid()){
      cout << "WARNING: no raw::OpDetWaveform producer " << fPMT_overlay_saturation_producer
           << " or no label OpdetBeamLowGain"<< endl;
    }
  }

  void CellTreeUB::processPMT_wfmSaturation(const art::Event& evt){
    art::Handle<std::vector<raw::OpDetWaveform> > opwf;
    if(! evt.getByLabel(fPMT_wfmSaturationLabel, opwf)){
      cout << "WARNING no raw::OpDetWaveform label " << fPMT_wfmSaturationLabel << endl;
      return;
    }
    art::ServiceHandle<geo::UBOpReadoutMap> ub_pmt_map;
    std::vector<raw::OpDetWaveform> const& opwf_v(*opwf);
    int i=0;
    for(auto const& op : opwf_v){
      unsigned int chan = (unsigned int)op.ChannelNumber();
      unsigned int crate, slot, femch;
      ub_pmt_map->GetCrateSlotFEMChFromReadoutChannel(chan, crate, slot, femch);
      if(i == 32) break;
      if(i != (int)femch && op.size() != 1501){ continue; }
      fOp_femch.push_back((short)femch);
      fOp_timestamp.push_back(op.TimeStamp());
      int nbins = (int)op.size();
      TH1S *h = new((*fOp_wf)[i]) TH1S("","",nbins,0,nbins);
      for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j]);}
      i++;
    }
    for(auto const& op : opwf_v){
      if(op.size() == 1501) continue;
      unsigned int chan = (unsigned int)op.ChannelNumber();
      unsigned int crate, slot, femch;
      ub_pmt_map->GetCrateSlotFEMChFromReadoutChannel(chan, crate, slot, femch);
      fOp_femch.push_back((short)femch);
      fOp_timestamp.push_back(op.TimeStamp());
      int nbins = (int)op.size();
      TH1S *h = new((*fOp_wf)[i]) TH1S("","",nbins,0,nbins);
      for(int j=1; j<=nbins; j++){h->SetBinContent(j,op[j]);}
      i++;
    }
  }

  void CellTreeUB::processPMT_reco(const art::Event& evt){
    art::Handle<std::vector<recob::OpHit> > ophit;
    if(! evt.getByLabel(fPMT_recoLabel, ophit)){
      cout << "WARNING: no reco::OpHit label " << fPMT_recoLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::OpHit> > oh_v;
    art::fill_ptr_vector(oh_v, ophit);
    oh_nHits = (int)oh_v.size();
    for(auto const& oh : oh_v){
      oh_channel.push_back(oh->OpChannel());
      oh_bgtime.push_back(oh->PeakTime());
      oh_trigtime.push_back(oh->PeakTimeAbs());
      oh_pe.push_back(oh->PE());
    }
    art::Handle<std::vector<recob::OpFlash> > opflash;
    if(! evt.getByLabel(fPMT_recoLabel, opflash)){
      cout << "WARNING: no reco::OpFlash label " << fPMT_recoLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::OpFlash> > of_v;
    art::fill_ptr_vector(of_v, opflash);
    of_nFlash = (int)of_v.size();
    int i=0;
    art::ServiceHandle<geo::Geometry> fGeometry;
    int nOpDet = fGeometry->NOpDets();
    for(auto const& of : of_v){
      of_t.push_back(of->Time());
      of_peTotal.push_back(of->TotalPE());
      TH1F *h = new ((*fPEperOpDet)[i]) TH1F("","",nOpDet,0,nOpDet);
      int mult = 0;
      for(int j=0; j<nOpDet; j++){
	if(of->PE(j) >= flashMultPEThreshold){
	  mult++;
	}
	h->SetBinContent(j,of->PE(j));
      }
      of_multiplicity.push_back(mult);
      i++;
    }
  }

  void CellTreeUB::processPMT_gainData(const art::Event& evt){
    art::ServiceHandle<geo::Geometry> geo;
    const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
    for(unsigned int a = 0; a != geo->NOpDets(); ++a){
      if(geo->IsValidOpChannel(a)){
	fOp_gain.push_back(gain_provider.Gain(a));
	fOp_gainerror.push_back(gain_provider.GainErr(a));
      }
    }
  }

  void CellTreeUB::processPMT_gainMC(const art::Event& evt){
    art::ServiceHandle<geo::Geometry> geo;
    for(unsigned int a = 0; a != geo->NOpDets(); ++a){
      if(geo->IsValidOpChannel(a)){
	fOp_gain.push_back(120.0);
	fOp_gainerror.push_back(0.05);
      }
    }
  }

  void CellTreeUB::processHWTrigger(const art::Event& evt){
    art::Handle<std::vector<raw::Trigger> > trigger;
    if(! evt.getByLabel(fHWTrigger_label, trigger)){
      cout << "WARNING: no raw::Trigger label " << fHWTrigger_label << endl;
      return;
    }
    std::vector<art::Ptr<raw::Trigger> > t_v;
    art::fill_ptr_vector(t_v, trigger);
    if(t_v.size()){
      fTriggernumber = t_v[0]->TriggerNumber();
      fTriggertime = t_v[0]->TriggerTime();
      fBeamgatetime = t_v[0]->BeamGateTime();
      fTriggerbits = t_v[0]->TriggerBits();
    }
    else{
      fTriggernumber = 0;
      fTriggertime = 0;
      fBeamgatetime = 0;
      fTriggerbits = 0;
    }
  }

  void CellTreeUB::processSWTrigger(const art::Event& evt){
    art::Handle<raw::ubdaqSoftwareTriggerData> swtrigger;
    if(! evt.getByLabel(fSWTrigger_label, swtrigger)){
      cout <<"WARNING: no raw::ubdaqSoftwareTriggerData label" << fSWTrigger_label << endl;
      return;
    }
    std::vector<std::string> algos = swtrigger->getListOfAlgorithms();
    for(auto const& algo : algos) {
      std::pair<std::string,bool> tempPair;
      tempPair.first = algo;
      tempPair.second = swtrigger->passedAlgo(algo);
      fpassAlgo.push_back(tempPair);
      fpassPrescale.push_back(swtrigger->passedPrescaleAlgo(algo));
      fPHMAX.push_back(swtrigger->getPhmax(algo));
      fmultiplicity.push_back(swtrigger->getMultiplicity(algo));
      ftriggerTick.push_back(swtrigger->getTriggerTick(algo));
      ftriggerTime.push_back(swtrigger->getTimeSinceTrigger(algo));
      fprescale_weight.push_back(swtrigger->getPrescale(algo));
    }
  }

  void CellTreeUB::processBeam(const art::Event& evt){
    art::Handle<std::vector<raw::BeamInfo> > beaminfo;
    if(! evt.getByLabel(fBeam_label, beaminfo)){
      cout << "WARNING: no raw::BeamInfo label " << fBeam_label << endl;
      return;
    }
    std::vector<art::Ptr<raw::BeamInfo> > bi_v;
    art::fill_ptr_vector(bi_v, beaminfo);
    if(bi_v.size()){
      ftor101 = bi_v[0]->get_tor101();
      ftortgt = bi_v[0]->get_tortgt();
      ftrtgtd = bi_v[0]->get_trtgtd();
      ft_ms = bi_v[0]->get_t_ms();
      fRecordType = bi_v[0]->GetRecordType();
      fSeconds = bi_v[0]->GetSeconds();
      fMilliSeconds = bi_v[0]->GetMilliSeconds();
      fNumberOfDevices = bi_v[0]->GetNumberOfDevices();
      fDataMap = bi_v[0]->GetDataMap();
    }
  }

  void CellTreeUB::processPOT(const art::SubRun& subrun){
    art::Handle<sumdata::POTSummary> pots1;
    if(! subrun.getByLabel(fPOT_producer, fPOT_bnb1label, pots1)){
	cout << "WARNING:  no sumdata::POTSummary producer " << fPOT_producer
	     << " or label " << fPOT_bnb1label << endl;
	return;
      }
    sumdata::POTSummary const& p1(*pots1);
    ftotpot_bnbETOR875 = p1.totpot;
    ftotgoodpot_bnbETOR875 = p1.totgoodpot;
    ftotspills_bnbETOR875 = p1.totspills;
    fgoodspills_bnbETOR875 = p1.goodspills;

    art::Handle<sumdata::POTSummary> pots2;
    if(! subrun.getByLabel(fPOT_producer, fPOT_bnb2label, pots2)){
	cout << "WARNING:  no sumdata::POTSummary producer " << fPOT_producer
	     << " or label " << fPOT_bnb2label << endl;
	return;
      }
    sumdata::POTSummary const& p2(*pots2);
    ftotpot_bnbETOR860 = p2.totpot;
    ftotgoodpot_bnbETOR860 = p2.totgoodpot;
    ftotspills_bnbETOR860 = p2.totspills;
    fgoodspills_bnbETOR860 = p2.goodspills;

    art::Handle<sumdata::POTSummary> pots3;
    if(! subrun.getByLabel(fPOT_producer, fPOT_numilabel, pots3)){
	cout << "WARNING:  no sumdata::POTSummary producer " << fPOT_producer
	     << " or label " << fPOT_numilabel << endl;
	return;
      }
    sumdata::POTSummary const& p3(*pots3);
    ftotpot_numiETORTGT = p3.totpot;
    ftotgoodpot_numiETORTGT = p3.totgoodpot;
    ftotspills_numiETORTGT = p3.totspills;
    fgoodspills_numiETORTGT = p3.goodspills;
  }

  void CellTreeUB::processMuCS_data(const art::Event& evt){
    art::Handle<std::vector<MuCS::MuCSData> > mdata;
    if(!evt.getByLabel(fMuCS_dataLabel, mdata)){
      cout << "WARNING: no MuCS::MuCSData label "<<fMuCS_dataLabel<<endl;
      return;
    }
    std::vector<art::Ptr<MuCS::MuCSData> > data_v;
    art::fill_ptr_vector(data_v, mdata);
    for(auto const& d : data_v){ ft0 = d->T0(); }
  }

  void CellTreeUB::processMuCS_recoData(const art::Event& evt){
    art::Handle<std::vector<MuCS::MuCSRecoData> > mreco;
    if(!evt.getByLabel(fMuCS_recoDataLabel, mreco)){
      cout<<"WARNING: no MuCS::MuCSRecoData label "<<fMuCS_recoDataLabel<<endl;
      return;
    }
    std::vector<art::Ptr<MuCS::MuCSRecoData> > reco_v;
    art::fill_ptr_vector(reco_v, mreco);
    for(auto const& r : reco_v){
      ftheta_xy = r->theta_xy();
      ftheta_xy_rms = r->theta_xy_rms();
      fx = r->x();
      fx_rms = r->x_rms();
      ftheta_yz = r->theta_yz();
      ftheta_yz_rms = r->theta_yz_rms();
      fz = r->z();
      fz_rms = r->z_rms();
      fy = r->y();
      ftheta = r->theta();
      fphi = r->phi();
      fxmatches = r->xmatches();
      fzmatches = r->zmatches();
    }
  }

  void CellTreeUB::processMCParticle(const art::Event& evt){
    art::Handle<std::vector<simb::MCParticle> > particleHandle;
    if(!evt.getByLabel("largeant", particleHandle)) return;
    std::vector<art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);
    art::FindOneP<simb::MCTruth> fo(particleHandle, evt, "largeant");

    int i=0; // track index in saved MCParticles
    int i_all=0; //track index in all MCParticles
    for(auto const& particle : particles){
      art::Ptr<simb::MCTruth> mctruth = fo.at(i_all);
      i_all++;
      if(mcOption == "nuOnly"){
	if( !(mctruth->Origin() == 1 && particle->Mother()==0) ){ continue; }
      }

      mc_process[i] = processMap[particle->Process()];
      mc_id[i] = particle->TrackId();
      mc_pdg[i] = particle->PdgCode();
      mc_mother[i] = particle->Mother();
      savedMCTrackIdMap[mc_id[i]] = mc_pdg[i];

      int Ndaughters = particle->NumberDaughters();
      vector<int> daughters;
      for(int i=0; i<Ndaughters; i++){ daughters.push_back(particle->Daughter(i)); }
      mc_daughters.push_back(daughters);

      size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
      int last = numberTrajectoryPoints - 1;
      const TLorentzVector& positionStart = particle->Position(0);
      const TLorentzVector& positionEnd = particle->Position(last);
      const TLorentzVector& momentumStart = particle->Momentum(0);
      const TLorentzVector& momentumEnd = particle->Momentum(last);
      positionStart.GetXYZT(mc_startXYZT[i]);
      positionEnd.GetXYZT(mc_endXYZT[i]);
      momentumStart.GetXYZT(mc_startMomentum[i]);
      momentumEnd.GetXYZT(mc_endMomentum[i]);

      if(fSaveMCTrackPoints){
	TClonesArray *Lposition = new TClonesArray("TLorentzVector",numberTrajectoryPoints);
	for(unsigned int j=0; j<numberTrajectoryPoints; j++){
	  new ((*Lposition)[j]) TLorentzVector(particle->Position(j));
	}
	fMC_trackPosition->Add(Lposition);
      }
      i++;
      if(i==MAX_TRACKS){ break; }
    }
    mc_Ntrack = i;

    art::Handle<std::vector<simb::MCTruth> > mctruthListHandle;
    evt.getByLabel("generator",mctruthListHandle);
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mctruthListHandle);
    art::Ptr<simb::MCTruth> mctruth;

    if(mclist.size()>0){
      mctruth = mclist.at(0);
      if(mctruth->NeutrinoSet()){
	simb::MCNeutrino nu = mctruth->GetNeutrino();
	mc_isnu = 1;
	mc_nGeniePrimaries = mctruth->NParticles();
	mc_nu_pdg = nu.Nu().PdgCode();
	mc_nu_ccnc = nu.CCNC();
	mc_nu_mode = nu.Mode();
	mc_nu_intType = nu.InteractionType();
	mc_nu_target = nu.Target();
	mc_hitnuc = nu.HitNuc();
	mc_hitquark = nu.HitQuark();
	mc_nu_Q2 = nu.QSqr();
	mc_nu_W = nu.W();
	mc_nu_X = nu.X();
	mc_nu_Y = nu.Y();
	mc_nu_Pt = nu.Pt();
	mc_nu_Theta = nu.Theta();

	const TLorentzVector& position = nu.Nu().Position(0);
	const TLorentzVector& momentum = nu.Nu().Momentum(0);
	position.GetXYZT(mc_nu_pos);
	momentum.GetXYZT(mc_nu_mom);
      }
    }
  }

  void CellTreeUB::processSimEnergyDepositInside(const art::Event& evt){
    art::Handle<std::vector<sim::SimEnergyDeposit> > sed;
    if(! evt.getByLabel(fSimEnergyDepositInside_producer, fSimEnergyDepositInside_label, sed)){
      cout << "WARNING: no sim::SimEnergyDeposit (inside) with producer " << fSimEnergyDepositInside_producer
           << " or label " << fSimEnergyDepositInside_label << endl;
      return;
    }
    vector<sim::SimEnergyDeposit> const& sed_v(*sed);
    for(size_t i=0; i<sed_v.size(); i++){
      sim::SimEnergyDeposit s = sed_v[i];
      sedi_nelectrons.push_back(s.NumElectrons());
      sedi_nphotons.push_back(s.NumPhotons());
      sedi_time_start.push_back(s.T0());
      sedi_time_end.push_back(s.T1());
      sedi_x_start.push_back(s.Start().X());
      sedi_x_end.push_back(s.End().X());
      sedi_y_start.push_back(s.Start().Y());
      sedi_y_end.push_back(s.End().Y());
      sedi_z_start.push_back(s.Start().Z());
      sedi_z_end.push_back(s.End().Z());
      sedi_pdg.push_back(s.PdgCode());
      sedi_trackId.push_back(s.TrackID());
      sedi_energy.push_back(s.Energy());
    }
  }

  void CellTreeUB::processSimEnergyDepositOutside(const art::Event& evt){
    art::InputTag sedTag(fSimEnergyDepositOutside_producer, fSimEnergyDepositOutside_label);
    art::Handle<std::vector<sim::SimEnergyDeposit> > sed;
    if(! evt.getByLabel(fSimEnergyDepositOutside_producer, fSimEnergyDepositOutside_label, sed)){
      cout << "WARNING: no sim::SimEnergyDeposit (outside) with producer " << fSimEnergyDepositOutside_producer
           << " or label " << fSimEnergyDepositOutside_label << endl;
      return;
    }
    vector<sim::SimEnergyDeposit> const& sed_v(*sed);
    for(size_t i=0; i<sed_v.size(); i++){
      sim::SimEnergyDeposit s = sed_v[i];
      sedo_nelectrons.push_back(s.NumElectrons());
      sedo_nphotons.push_back(s.NumPhotons());
      sedo_time_start.push_back(s.T0());
      sedo_time_end.push_back(s.T1());
      sedo_x_start.push_back(s.Start().X());
      sedo_x_end.push_back(s.End().X());
      sedo_y_start.push_back(s.Start().Y());
      sedo_y_end.push_back(s.End().Y());
      sedo_z_start.push_back(s.Start().Z());
      sedo_z_end.push_back(s.End().Z());
      sedo_pdg.push_back(s.PdgCode());
      sedo_trackId.push_back(s.TrackID());
      sedo_energy.push_back(s.Energy());
    }
  }

  void CellTreeUB::processSimChannel(const art::Event& evt){
    art::InputTag scTag(fSimChannel_producer,fSimChannel_label);
    auto const& scHandle = evt.getValidHandle<std::vector<sim::SimChannel>>(scTag);
    auto const& sc_vec(*scHandle);
    for(size_t i=0; i<sc_vec.size(); i++){
      sim::SimChannel sc = sc_vec[i];
      auto const timeSlices = sc.TDCIDEMap();
      for(auto timeSlice : timeSlices){
	auto const& energyDeposits = timeSlice.second;
	for(auto energyDeposit : energyDeposits){
	  sc_x.push_back(energyDeposit.x);
	  sc_y.push_back(energyDeposit.y);
	  sc_z.push_back(energyDeposit.z);
	  sc_q.push_back(energyDeposit.numElectrons);
	  sc_chan.push_back(sc.Channel());
	  sc_tdc.push_back(timeSlice.first);
	  unsigned int tick = (unsigned int) timeSlice.first - ((4050.0 + 1600.0)/0.5);
	  sc_tick.push_back(tick);
	}
      }
    }
  }

  DEFINE_ART_MODULE(CellTreeUB)
} // namespace wc
#endif

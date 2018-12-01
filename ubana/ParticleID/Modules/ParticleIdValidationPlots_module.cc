////////////////////////////////////////////////////////////////////////
// Class:       ParticleIdValidationPlots
// Plugin Type: analyzer (art v2_05_01)
// File:        ParticleIdValidationPlots_module.cc
//
// Generated at Thu Mar 15 14:40:38 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 *
 * \class ParticleIdValidationPlots
 *
 * \brief ParticleId analyzer module
 *
 * \author Kirsty Duffy (kduffy@fnal.gov), Adam Lister (a.lister1@lancaster.ac.uk)
 *
 * \date 2018/04/18
 *
 * \notes
 * - There are some checks in here to make sure that the plane id is between
 *   0 and 2. This is because the calorimetry module will return a value of
 *   -1 when there are no hits in a specific plane.
 *
 */

// Art
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// UBXSec
#include "ubobj/UBXSec/SelectionResult.h"
#include "ubobj/UBXSec/TPCObject.h"

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

// Algorithms
#include "ubana/ubana/ParticleID/Algorithms/FiducialVolume.h"
#include "ubana/ubana/ParticleID/Algorithms/dQdxSeparatorMarco.h"
#include "ubana/ubana/ParticleID/Algorithms/Bragg_Likelihood_Estimator.h"
#include "ubana/ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

// cpp
#include <vector>

class ParticleIdValidationPlots;


class ParticleIdValidationPlots : public art::EDAnalyzer {
  public:
    explicit ParticleIdValidationPlots(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ParticleIdValidationPlots(ParticleIdValidationPlots const &) = delete;
    ParticleIdValidationPlots(ParticleIdValidationPlots &&) = delete;
    ParticleIdValidationPlots & operator = (ParticleIdValidationPlots const &) = delete;
    ParticleIdValidationPlots & operator = (ParticleIdValidationPlots &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endSubRun(art::SubRun const &sr) override;

  private:

    art::ServiceHandle<art::TFileService> tfs;
    fidvol::fiducialVolume fid;

    std::vector<double> fv;

    bool fIsDataPlots;
    bool fIsUBXSecSelected;
    std::string fTrackLabel;
    std::string fHitLabel;
    std::string fHitTrackAssns;
    std::string fCaloTrackAssns;
    std::string fHitTruthAssns;
    std::string fPIDLabel;
    //std::string fPIDLabelChi2;
    int fNHitsForTrackDirection;

    bool isData;

    /** Setup root trees  */
    TTree *potTree;
    double sr_pot = 0;
    int sr_run = 0;
    int sr_sub_run = 0;

    TTree *pidTree;

    int run = -999;
    int sub_run = -999;
    int event = -999;
    int true_PDG = -999;
    double true_purity = -999;
    double true_start_momentum = -999;
    double true_start_x = -999;
    double true_start_y = -999;
    double true_start_z = -999;
    double true_end_momentum = -999;
    double true_end_x = -999;
    double true_end_y = -999;
    double true_end_z = -999;
    int track_id = -999;
    double track_start_x;
    double track_start_y;
    double track_start_z;
    double track_end_x;
    double track_end_y;
    double track_end_z;
    std::vector<double> track_likelihood_fwd_mu;
    std::vector<double> track_likelihood_fwd_p;
    std::vector<double> track_likelihood_fwd_pi;
    std::vector<double> track_likelihood_fwd_k;
    std::vector<double> track_likelihood_fwd_mip;
    std::vector<double> track_likelihood_bwd_mu;
    std::vector<double> track_likelihood_bwd_p;
    std::vector<double> track_likelihood_bwd_pi;
    std::vector<double> track_likelihood_bwd_k;
    std::vector<double> track_likelihood_shift_fwd_mu;
    std::vector<double> track_likelihood_shift_fwd_p;
    std::vector<double> track_likelihood_shift_fwd_pi;
    std::vector<double> track_likelihood_shift_fwd_k;
    std::vector<double> track_likelihood_shift_bwd_mu;
    std::vector<double> track_likelihood_shift_bwd_p;
    std::vector<double> track_likelihood_shift_bwd_pi;
    std::vector<double> track_likelihood_shift_bwd_k;
    std::vector<double> track_PIDA_mean;
    std::vector<double> track_PIDA_median;
    std::vector<double> track_PIDA_kde;
    std::vector<double> track_dEdx;
    std::vector<double> track_dQdx;
    std::vector<double> track_depE;
    std::vector<double> track_nhits;
    std::vector<double> track_Chi2Proton;
    std::vector<double> track_Chi2Kaon;
    std::vector<double> track_Chi2Pion;
    std::vector<double> track_Chi2Muon;
    double track_length;
    double track_theta;
    double track_phi;
    double track_rangeE_mu;
    double track_rangeE_p;
    bool track_dQdxtruncmeanvslength_isMuon;
    std::vector<float> track_dEdx_perhit_u;
    std::vector<float> track_dEdx_perhit_v;
    std::vector<float> track_dEdx_perhit_y;
    std::vector<float> track_resrange_perhit_u;
    std::vector<float> track_resrange_perhit_v;
    std::vector<float> track_resrange_perhit_y;
    std::vector<std::vector<float>> track_Lmip_perhit;
    std::vector<std::vector<float>> dEdx;
    std::vector<std::vector<float>> resRange;


    /** Histograms for all tracks, i.e. can be used by data */
    TH2F *All_chargeEndOverStart_sm0_5_dEdxrr;
    TH2F *All_chargeEndOverStart_gr2_dEdxrr;
    TH2F *All_chargeEndOverStart_0_5to2_dEdxrr;

    /** Histograms for tracks expected to have Bragg peaks */
    TH1F *TrueBragg_chargeEndOverStart_directionCorrect;
    TH1F *TrueBragg_chargeEndOverStart_directionIncorrect;
    TH2F *TrueBragg_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect;

    TH2F *TrueBragg_chargeEndOverStart_sm0_5_dEdxrr;
    TH2F *TrueBragg_chargeEndOverStart_gr2_dEdxrr;
    TH2F *TrueBragg_chargeEndOverStart_0_5to2_dEdxrr;

    TH2F *TrueBragg_correctdirection;
    TH2F *TrueBragg_incorrectdirection;
    TH1F *TrueBragg_PIDdir;

    TH2F *TrueBragg_correctdirection_PIDdir;
    TH2F *TrueBragg_incorrectdirection_PIDdir;

    /** All tracks which are matched to an MCParticle */
    TH1F *All_chargeEndOverStart_directionCorrect;
    TH1F *All_chargeEndOverStart_directionIncorrect;
    TH2F *All_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *All_chargeEndOverStartVersusNHits_directionIncorrect;

    TH1F *Contained_chargeEndOverStart_directionCorrect;
    TH1F *Contained_chargeEndOverStart_directionIncorrect;
    TH2F *Contained_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *Contained_chargeEndOverStartVersusNHits_directionIncorrect;

    TH2F *All_correctdirection;
    TH2F *All_incorrectdirection;

    // for likelihood-based PID
    particleid::Bragg_Likelihood_Estimator braggcalc;
};

ParticleIdValidationPlots::ParticleIdValidationPlots(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolume");
  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet const p_bragg  = p.get<fhicl::ParameterSet>("BraggAlgo");

  fIsDataPlots = p.get<bool>("IsDataPlotsOnly", "false");
  fIsUBXSecSelected = p.get<bool>("IsUBXSecSelected", "false");
  fTrackLabel = p_labels.get<std::string>("TrackLabel","myTracks::PandoraWorkshopTrackShower");
  fHitLabel = p_labels.get<std::string>("HitLabel","gaushit::McRecoStage1");
  //This should just be the neutrino slice........
  fHitTrackAssns = p_labels.get<std::string>("HitTrackAssn","myTracks::PandoraWorkshopTrackShower");
  fCaloTrackAssns = p_labels.get<std::string>("CaloTrackAssn", "pandoracalo::McRecoStage2");
  fHitTruthAssns = p_labels.get<std::string>("HitTruthAssn","gaushitTruthMatch::McRecoStage1");
  fPIDLabel = p_labels.get<std::string>("ParticleIdLabel");
  //fPIDLabelChi2 = p_labels.get<std::string>("ParticleIdChi2Label");
  fNHitsForTrackDirection = p.get<int>("NHitsForTrackDirection");

  fv = fid.setFiducialVolume(fv, p_fv);
  fid.printFiducialVolume(fv);

  braggcalc.configure(p_bragg);
  braggcalc.checkRange=false;
  braggcalc.nHitsToDrop = 0;
  braggcalc.printConfiguration();

  std::cout << "[ParticleIdValidation] >> Use only UBXSec TPCObj-tracks? " << fIsUBXSecSelected << std::endl;
  std::cout << "[ParticleIDValidation] >> Track label: " << fTrackLabel << std::endl;
  std::cout << "[ParticleIDValidation] >> Hit label: " << fHitLabel << std::endl;
  std::cout << "[ParticleIDValidation] >> Hit-track assns: " << fHitTrackAssns << std::endl;
  std::cout << "[ParticleIDValidation] >> Hit-truth assns: " << fHitTruthAssns << std::endl;
  std::cout << "[ParticleIDValidation] >> Calo-track assns: " << fCaloTrackAssns << std::endl;
  std::cout << "[ParticleIDValidation] >> ParticleID label: " << fPIDLabel << std::endl;
  //std::cout << "[ParticleIDValidation] >> ParticleID Chi2 label: " << fPIDLabelChi2 << std::endl;
  std::cout << "[ParticleIDValidation] >> NHits for track direction: " << fNHitsForTrackDirection << std::endl;

}

void ParticleIdValidationPlots::analyze(art::Event const & e)
{
  run = e.run();
  sub_run = e.subRun();
  event = e.event();

  isData = e.isRealData();

  if (!isData) std::cout << "[ParticleIDValidation] Running simulated data." << std::endl;
  else std::cout << "[ParticleIDValidation] Running physics data." << std::endl;

  /**
   * Get handles to needed information
   */
  //std::cout << "hello0!" << std::endl;
  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector<art::Ptr<recob::Track>> trackCollection;

  //std::cout << "hello1!" << std::endl;
  /**
   * Two options for trackCollection, if isUBXSecSelected then
   * only use tracks actually in the selected TPC object,
   * if not then use all tracks which pass the cosmic filter
   */

  // first fill ptr vector
  art::fill_ptr_vector(trackCollection, trackHandle);
  std::vector<art::Ptr<recob::Track>> trackPtrVector;

  if (fIsUBXSecSelected == true){

    art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
    e.getByLabel("UBXSec", selectionHandle);

    art::FindManyP<ubana::TPCObject> selectedTPCObjects(selectionHandle, e, "UBXSec");
    art::Ptr<ubana::TPCObject> selectedTPCObject;

    if (selectedTPCObjects.at(0).size() == 1){

      selectedTPCObject = selectedTPCObjects.at(0).at(0);

      const std::vector<recob::Track>& selectedTracks = selectedTPCObject->GetTracks();

      for (size_t i = 0; i < selectedTracks.size(); i++){

        recob::Track selectedTrack = selectedTracks.at(i);
        int selectedTrackId = selectedTrack.ID();

        // loop all tracks in the collection and find this track, then add the
        // pointer to that track to the trackPtrVector

        for (auto& track : trackCollection){

          int testTrackId = track->ID();

          std::cout << "[ParticleIDValidation] Checking track " << testTrackId << std::endl;

          if (testTrackId == selectedTrackId){

            std::cout << "[ParticleIDValidation] Found track ID match with ID " << testTrackId << std::endl;
            trackPtrVector.push_back(track);

          }

        }

      }

    }

  }
  else {
    // if fIsUBXSecSelected == false, then we just want all tracks
    // and so the trackPtrVector is just the trackCollection from
    // earlier
    trackPtrVector = trackCollection;
  }

  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  art::FindManyP<recob::Hit> hits_from_tracks(trackHandle, e, fHitTrackAssns);
  art::FindManyP<anab::Calorimetry> calo_from_tracks(trackHandle, e, fCaloTrackAssns);

  /**
   * Variables which need to have scope throughout the code
   */
  TVector3 true_start;
  TVector3 true_end;

  for (auto& track : trackPtrVector){
    //std::cout << "found track" << std::endl;

    /** reset default values */
    dEdx.resize(3);
    resRange.resize(3);
    track_Lmip_perhit.resize(3);
    track_likelihood_fwd_mu.resize(3,-999.);
    track_likelihood_fwd_p.resize(3,-999.);
    track_likelihood_fwd_pi.resize(3,-999.);
    track_likelihood_fwd_k.resize(3,-999.);
    track_likelihood_fwd_mip.resize(3,-999.);
    track_likelihood_bwd_mu.resize(3,-999.);
    track_likelihood_bwd_p.resize(3,-999.);
    track_likelihood_bwd_pi.resize(3,-999.);
    track_likelihood_bwd_k.resize(3,-999.);
    track_likelihood_shift_fwd_mu.resize(3,-999.);
    track_likelihood_shift_fwd_p.resize(3,-999.);
    track_likelihood_shift_fwd_pi.resize(3,-999.);
    track_likelihood_shift_fwd_k.resize(3,-999.);
    track_likelihood_shift_bwd_mu.resize(3,-999.);
    track_likelihood_shift_bwd_p.resize(3,-999.);
    track_likelihood_shift_bwd_pi.resize(3,-999.);
    track_likelihood_shift_bwd_k.resize(3,-999.);
    track_Chi2Muon.resize(3,-999.);
    track_Chi2Proton.resize(3, -999.);
    track_Chi2Kaon.resize(3, -999.);
    track_Chi2Pion.resize(3,-999.);
    track_PIDA_mean.resize(3,-999.);
    track_PIDA_median.resize(3,-999.);
    track_PIDA_kde.resize(3,-999.);
    track_dEdx.resize(3,-999.);
    track_dQdx.resize(3,-999.);
    track_depE.resize(3,-999.);
    track_nhits.resize(3,-999);

    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = calo_from_tracks.at(track->ID());

    track_id = track->ID();
    track_length = track->Length();
    track_theta = track->Theta();
    track_phi = track->Phi();
    track_start_x = track->Start().X();
    track_start_y = track->Start().Y();
    track_start_z = track->Start().Z();
    track_end_x = track->End().X();
    track_end_y = track->End().Y();
    track_end_z = track->End().Z();

    bool TrueBragg = false;
    true_PDG = 0;
    simb::MCParticle const* matched_mcparticle = NULL;
  
    if (!fIsDataPlots){

      /**
       * Get true PDG from associations.
       * We do this by using hit <-> MCParticle associations, looping over the
       * hits and finding the MCParticle which contributed the most charge
       * to each hit.
       *
       * .key() is used to get the index in the original collection
       */

      std::unordered_map<int,double> trkide;
      double maxe=-1, tote=0;
      double purity = -999;
  
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
   
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track->ID());
      
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitHandle,e,fHitTruthAssns);
      //std::cout << "LISTEN HERE!" << particles_per_hit.size()<< std::endl;
      
      for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
        //std::cout << "hello*!, 431" << std::endl;
	particle_vec.clear(); match_vec.clear();
	//std::cout << "hello*!, 433" << std::endl;
        particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
	//std::cout << "hello*!, 435" << std::endl;
      
        for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
          trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
          tote += match_vec[i_p]->energy;
          if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
            maxe = trkide[ particle_vec[i_p]->TrackId() ];
            matched_mcparticle = particle_vec[i_p];
          }
        }//end loop over particles per hit

        purity = maxe/tote;

      }

      true_purity = purity;
      true_PDG = matched_mcparticle->PdgCode();
      true_start_momentum = matched_mcparticle->P();
      true_start_x = matched_mcparticle->Vx();
      true_start_y = matched_mcparticle->Vy();
      true_start_z = matched_mcparticle->Vz();
      true_end_momentum =
        std::sqrt( std::pow(matched_mcparticle->EndMomentum().X(),2)
            + std::pow(matched_mcparticle->EndMomentum().Y(),2)
            + std::pow(matched_mcparticle->EndMomentum().Z(),2));
      true_end_x = matched_mcparticle->EndX();
      true_end_y = matched_mcparticle->EndY();
      true_end_z = matched_mcparticle->EndZ();

      true_start = {true_start_x, true_start_y, true_start_z};
      true_end = {true_end_x, true_end_y, true_end_z};

      // Check if true particle should have a Bragg peak
      TrueBragg = false;
      if (matched_mcparticle->EndPx() == 0
          && matched_mcparticle->EndPy() == 0
          && matched_mcparticle->EndPz() == 0
          && fid.isInFiducialVolume(true_end, fv)){
        TrueBragg = true;
      }

    } // end if(!fIsDataPlots)
   
    std::vector<std::vector<float>> Lmip_perhit(3);

    art::Ptr< anab:: Calorimetry > calo;
    int planenum = -1;
    int hitsToUse = 0;
    for (auto c : caloFromTrack){
      planenum = c->PlaneID().Plane;
      calo = c;
      if (planenum < 0 || planenum > 2){
        std::cout << "[ParticleIDValidation] No calorimetry information for plane "
          << planenum << std::endl;
        continue;
      }
      dEdx.at(planenum).clear();
      resRange.at(planenum).clear();
      dEdx.at(planenum) = calo->dEdx();
      resRange.at(planenum) = calo->ResidualRange();

      // Get MIP likelihood per hit (and store in a vector)
      // Loop through hits (entries in dEdx and resrange vector)
      // For each hit, make a new dEdx vector for that hit only and use it to get the likelihood for that hit. That way we can average the likelihoods over the number of hits we care about later.
      std::vector<float> dEdx_dummy = {0.};
      std::vector<float> rr_dummy = {0.};
      for (size_t i_hit=0; i_hit < dEdx.at(planenum).size(); i_hit++){
        double Lmip = -9999.;
        dEdx_dummy.at(0) = dEdx.at(planenum).at(i_hit);
        rr_dummy.at(0) = resRange.at(planenum).at(i_hit);

        Lmip = braggcalc.getLikelihood(dEdx_dummy,rr_dummy,0,1,planenum);
        Lmip_perhit.at(planenum).push_back(Lmip);
      }
      
      /**
       * Get hit charge of first and final 5 hits of track to try and find the
       * direction of the track. If there are fewer than 10 hits then take half
       * of the total hits.
       */

      double nhits = resRange.at(planenum).size();
      // find how many hits to use
      if (nhits >= 2*fNHitsForTrackDirection)
        hitsToUse = fNHitsForTrackDirection;
      else{
        hitsToUse = std::floor((double)nhits/2.0);
      }

      track_nhits.at(planenum) = nhits;

    }
    double averagedEdxTrackStart=0;
    double averagedEdxTrackEnd=0;

    // loop dEdx and take average of first n hits and last n hits

    TVector3 RecoTrackDir, TrueTrackDir;

    if (dEdx.at(2).size() > 0) {

      for (int i = 0; i < (int)dEdx.at(2).size(); i++){

        if (i < hitsToUse) averagedEdxTrackStart+=dEdx.at(2).at(i);
        else if (i > track_nhits.at(2) - hitsToUse -1) averagedEdxTrackEnd+=dEdx.at(2).at(i);

      }

      averagedEdxTrackStart = averagedEdxTrackStart/hitsToUse;
      averagedEdxTrackEnd   = averagedEdxTrackEnd/hitsToUse;

      double dEdxStartEndRatio = averagedEdxTrackEnd/averagedEdxTrackStart;

      /**
       * Now check the reconstructed direction we found and compare it against
       * the true direction.
       */

      if (!fIsDataPlots){

        RecoTrackDir = {track->End().X()-track->Start().X(),
          track->End().Y()-track->Start().Y(),
          track->End().Z()-track->Start().Z()};
        RecoTrackDir = RecoTrackDir.Unit();

        TrueTrackDir = true_end - true_start;
        TrueTrackDir = TrueTrackDir.Unit();

        // If dot product < 0, track direction is wrong
        if (RecoTrackDir.Dot(TrueTrackDir) < 0){

          All_chargeEndOverStart_directionIncorrect->Fill(dEdxStartEndRatio);
          All_chargeEndOverStartVersusNHits_directionIncorrect->Fill(dEdxStartEndRatio, hitsToUse);

          if (TrueBragg){
            TrueBragg_chargeEndOverStart_directionIncorrect->Fill(dEdxStartEndRatio);
            TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect->Fill(dEdxStartEndRatio, hitsToUse);
          }

          if (fid.isInFiducialVolume(true_start, fv) && fid.isInFiducialVolume(true_end, fv)){
            Contained_chargeEndOverStart_directionIncorrect->Fill(dEdxStartEndRatio);
            Contained_chargeEndOverStartVersusNHits_directionIncorrect->Fill(dEdxStartEndRatio, hitsToUse);
          }

        }
        // else, track direction is correct
        else{

          All_chargeEndOverStart_directionCorrect->Fill(dEdxStartEndRatio);
          All_chargeEndOverStartVersusNHits_directionCorrect->Fill(dEdxStartEndRatio, hitsToUse);

          if (TrueBragg){
            TrueBragg_chargeEndOverStart_directionCorrect->Fill(dEdxStartEndRatio);
            TrueBragg_chargeEndOverStartVersusNHits_directionCorrect->Fill(dEdxStartEndRatio, hitsToUse);
          }

          if (fid.isInFiducialVolume(true_start, fv) && fid.isInFiducialVolume(true_end, fv)){
            Contained_chargeEndOverStart_directionCorrect->Fill(dEdxStartEndRatio);
            Contained_chargeEndOverStartVersusNHits_directionCorrect->Fill(dEdxStartEndRatio, hitsToUse);
          }

        }

      }

      /**
       * Make plots of dEdx versus residual range for different track start/end
       * ratios.
       *
       * Note that here we don't need to check if fIsDataPlots because TrueBragg is
       * false by construction if fIsDataPlots = true
       */

      if (dEdxStartEndRatio < 0.5){
        for (int i=0; i < (int)resRange.at(2).size(); i++){
          All_chargeEndOverStart_sm0_5_dEdxrr->Fill(resRange.at(2).at(i),dEdx.at(2).at(i));
          if (TrueBragg){
            TrueBragg_chargeEndOverStart_sm0_5_dEdxrr->Fill(resRange.at(2).at(i),dEdx.at(2).at(i));
          }
        }
      }
      else if (dEdxStartEndRatio > 2.0){
        for (int i=0; i < (int)resRange.at(2).size(); i++){
          All_chargeEndOverStart_gr2_dEdxrr->Fill(resRange.at(2).at(i),dEdx.at(2).at(i));
          if (TrueBragg){
            TrueBragg_chargeEndOverStart_gr2_dEdxrr->Fill(resRange.at(2).at(i),dEdx.at(2).at(i));
          }
        }
      }
      else {
        for (int i=0; i < (int)resRange.at(2).size(); i++){
          All_chargeEndOverStart_0_5to2_dEdxrr->Fill(resRange.at(2).at(i),dEdx.at(2).at(i));
          if (TrueBragg){
            TrueBragg_chargeEndOverStart_0_5to2_dEdxrr->Fill(resRange.at(2).at(i),dEdx.at(2).at(i));
          }
        }
      }
    }

    /**
     * Calculate PID variables
     */

    art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, fPIDLabel);
    if (!trackPIDAssn.isValid()){
      std::cout << "[ParticleIDValidation] trackPIDAssn.isValid() == false. Skipping track." << std::endl;
      continue;
    }

    double trklen = -999;
    double rangeE_mu = -999;
    double rangeE_p = -999;

    std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track->ID());
    if (trackPID.size() == 0){
      std::cout << "[ParticleIDValidation] No track-PID association found for trackID " << track->ID() << ". Skipping track." << std::endl;
      continue;
    }
    
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

    // Loop through AlgScoresVec and find the variables we want
    for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){

      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);

      if (planeid < 0 || planeid > 2){
        std::cout << "[ParticleIDValidation] No information for planeid " << planeid << std::endl;
        continue;
      }


      /**
       * Algorithm 1: BraggLikelihood
       */
      if (AlgScore.fAlgName == "BraggPeakLLH"){

        if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
          if (AlgScore.fAssumedPdg == 13)   track_likelihood_fwd_mu.at(planeid) = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 2212) track_likelihood_fwd_p.at(planeid) =  AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 211)  track_likelihood_fwd_pi.at(planeid) = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 321)  track_likelihood_fwd_k.at(planeid)  = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 0)    track_likelihood_fwd_mip.at(planeid) = AlgScore.fValue;
        }
        else if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward){
          if (AlgScore.fAssumedPdg == 13)   track_likelihood_bwd_mu.at(planeid) = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 2212) track_likelihood_bwd_p.at(planeid) =  AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 211)  track_likelihood_bwd_pi.at(planeid) = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 321)  track_likelihood_bwd_k.at(planeid) =  AlgScore.fValue;
        }

      } // if fAlgName = BraggPeakLLH

     if (AlgScore.fAlgName == "BraggPeakLLH_shift"){
       if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
         if (AlgScore.fAssumedPdg == 13)   track_likelihood_shift_fwd_mu.at(planeid) = AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 2212) track_likelihood_shift_fwd_p.at(planeid) =  AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 211)  track_likelihood_shift_fwd_pi.at(planeid) = AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 321)  track_likelihood_shift_fwd_k.at(planeid)  = AlgScore.fValue;
       }
       else if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward){
         if (AlgScore.fAssumedPdg == 13)   track_likelihood_shift_bwd_mu.at(planeid) = AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 2212) track_likelihood_shift_bwd_p.at(planeid) =  AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 211)  track_likelihood_shift_bwd_pi.at(planeid) = AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 321)  track_likelihood_shift_bwd_k.at(planeid) =  AlgScore.fValue;
       }

     } // if fAlgName = BraggPeakLLH

      /**
       * Algorithm 2: Chi2
       */

      if(AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
          if (AlgScore.fAssumedPdg == 13)   track_Chi2Muon.at(planeid) = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 2212) track_Chi2Proton.at(planeid) =  AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 211)  track_Chi2Pion.at(planeid) = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 321)  track_Chi2Kaon.at(planeid) =  AlgScore.fValue;
      }

      /**
       * Algorithm 3: PIDA
       */

      if (AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA)
        track_PIDA_mean.at(planeid) = AlgScore.fValue;

      if (AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA)
        track_PIDA_median.at(planeid) = AlgScore.fValue;

      if (AlgScore.fAlgName == "PIDA_kde" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA)
        track_PIDA_kde.at(planeid) = AlgScore.fValue;

      /**
       * Algorithm 4: truncated dE/dx versus residual range
       */
       // dQdx needs to be multiplied by a constant factor to convert from ADC to e/cm
       // Multiply MC by 198 and data by 243
      //FOR MCC9 TAG1 DATA IS NOW 242 on plane 2 and Tag2 MC is 178 on plane 2
       double dQdxcalibval = 178.;
       if (isData){
         dQdxcalibval = 243.;
       }

      if (AlgScore.fAlgName == "TruncatedMean"){
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) track_dEdx.at(planeid) = AlgScore.fValue;
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kdQdxtruncmean) track_dQdx.at(planeid) = AlgScore.fValue * dQdxcalibval;
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kTrackLength) trklen = AlgScore.fValue;
      }
      //std::cout << "planeid = " << planeid << ", trklen = " << trklen << ", dQdx = " << track_dQdx.at(planeid) << std::endl;

      /**
       * Algorithm 5: Deposited energy versus residual range
       */

      if (AlgScore.fAlgName == "DepEvsRangeE"){
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kEdeposited) track_depE.at(planeid) = AlgScore.fValue;
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kEbyRange){
          if (AlgScore.fAssumedPdg == 13) rangeE_mu = AlgScore.fValue;
          if (AlgScore.fAssumedPdg == 2212) rangeE_p = AlgScore.fValue;
        }
      }
    } // Loop over AlgScoresVec


    // Evaluate cut: is the track classed as a muon?
    // Only evaluate for plane 2
    track_dQdxtruncmeanvslength_isMuon=false;
    std::cout << "plane 2 trklen = " << trklen << ", dQdx = " << track_dQdx.at(2) << std::endl;
    // First: all long tracks are MIPs
    int l = std::round(trklen);
    if (l<0){
      track_dQdxtruncmeanvslength_isMuon=false;
    }
    else{
      if (l >= 1000){
        track_dQdxtruncmeanvslength_isMuon=true;
      }
      else{
        // Now evaluate Marco's cut
        double dqdx_cut = _dqdx_cutvals.at(l);

        std::cout << "dqdx_cut = " << dqdx_cut << std::endl;

        if (track_dQdx.at(2) <= dqdx_cut && track_dQdx.at(2)!=-999){
          track_dQdxtruncmeanvslength_isMuon = true;
        }

        std::cout << "[ParticleIDValidation >> MIPConsistencyCheck_Marco] Track length is " << trklen << ", dqdx_cut is " << dqdx_cut << std::endl;
      }

      std::cout << "[ParticleIDValidation >> MIPConsistencyCheck_Marco] \t Truncated mean dQ/dx for candidate is: " << track_dQdx.at(2) << std::endl;
      std::cout << "[ParticleIDValidation >> MIPConsistencyCheck_Marco] \t MIP consistent ? : " << (track_dQdxtruncmeanvslength_isMuon ? "YES" : "NO") << std::endl;
    }

    // Now time to set some variables!
    track_length = trklen;
    track_rangeE_mu = rangeE_mu;
    track_rangeE_p = rangeE_p;
    track_dEdx_perhit_u = dEdx.at(0);
    track_dEdx_perhit_v = dEdx.at(1);
    track_dEdx_perhit_y = dEdx.at(2);
    track_resrange_perhit_u = resRange.at(0);
    track_resrange_perhit_v = resRange.at(1);
    track_resrange_perhit_y = resRange.at(2);
    track_Lmip_perhit = Lmip_perhit;
    //std::cout << "done" << std::endl;
    /**
     * Testing to see whether we can predict the direction of the track from
     * the fwd/bwd likelihood variables
     */
    //std::cout << "hello2!" << std::endl;
    bool PID_fwd = false;
    double Bragg_smallest = std::min({track_likelihood_fwd_mu.at(2), track_likelihood_fwd_p.at(2), track_likelihood_fwd_pi.at(2), track_likelihood_fwd_k.at(2), track_likelihood_fwd_mip.at(2), track_likelihood_bwd_mu.at(2), track_likelihood_bwd_p.at(2), track_likelihood_bwd_pi.at(2), track_likelihood_bwd_k.at(2)});
    if (Bragg_smallest == track_likelihood_fwd_mu.at(2) || Bragg_smallest == track_likelihood_fwd_p.at(2) || Bragg_smallest == track_likelihood_fwd_pi.at(2) || Bragg_smallest == track_likelihood_fwd_k.at(2) || Bragg_smallest == track_likelihood_fwd_mip.at(2)) PID_fwd = true;

    //std::cout << "hello3!" << std::endl;
    // Histogram time
    if (!fIsDataPlots){
      if (TrueBragg){
        // Well-reconstructed truth matching
        if (TMath::Abs(true_PDG) == 13){ // True muons
          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(0.5,1.5);
            else TrueBragg_correctdirection->Fill(0.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(0.5,1.5);
            else TrueBragg_incorrectdirection->Fill(0.5,0.5);
          }
        }

        else if (TMath::Abs(true_PDG) == 2212){ // True protons
          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(1.5,1.5);
            else TrueBragg_correctdirection->Fill(1.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(1.5,1.5);
            else TrueBragg_incorrectdirection->Fill(1.5,0.5);
          }
        }

        else if (TMath::Abs(true_PDG) == 211){ // True pions
          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(2.5,1.5);
            else TrueBragg_correctdirection->Fill(2.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(2.5,1.5);
            else TrueBragg_incorrectdirection->Fill(2.5,0.5);
          }
        }

        else if (TMath::Abs(true_PDG) == 321){ // True kaons
          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(3.5,1.5);
            else TrueBragg_correctdirection->Fill(3.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(3.5,1.5);
            else TrueBragg_incorrectdirection->Fill(3.5,0.5);
          }
        }
      } // end if(TrueBragg)

      // All particles
      if (TMath::Abs(true_PDG) == 13){ // True muons
        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(0.5,1.5);
          else All_correctdirection->Fill(0.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(0.5,1.5);
          else All_incorrectdirection->Fill(0.5,0.5);
        }
      }

      else if (TMath::Abs(true_PDG) == 2212){ // True protons
        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(1.5,1.5);
          else All_correctdirection->Fill(1.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(1.5,1.5);
          else All_incorrectdirection->Fill(1.5,0.5);
        }
      }

      else if (TMath::Abs(true_PDG) == 211){ // True pions
        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(20.5,1.5);
          else All_correctdirection->Fill(2.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(2.5,1.5);
          else All_incorrectdirection->Fill(2.5,0.5);
        }
      }

      else if (TMath::Abs(true_PDG) == 321){ // True kaons
        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(3.5,1.5);
          else All_correctdirection->Fill(3.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(3.5,1.5);
          else All_incorrectdirection->Fill(3.5,0.5);
        }
      }
    }// end !fIsDataPlots

    // Finally, fill the tree
    std::cout << "[ParticleIDValidation] Filling tree. " << std::endl;
    pidTree->Fill();

  } // Loop over tracks


}

void ParticleIdValidationPlots::beginJob(){

  potTree = tfs->make<TTree>("potTree","potTree");
  potTree->Branch("sr_pot", &sr_pot, "sr_pot/D");
  potTree->Branch("sr_run", &sr_run, "sr_run/I");
  potTree->Branch("sr_sub_run", &sr_sub_run, "sr_sub_run/I");

  pidTree = tfs->make<TTree>("pidTree" , "pidTree");

  pidTree->Branch( "run"                      , &run                 ) ;
  pidTree->Branch( "sub_run"                  , &sub_run             ) ;
  pidTree->Branch( "event"                    , &event               ) ;
  pidTree->Branch( "true_PDG"                 , &true_PDG            ) ;
  pidTree->Branch( "true_purity"              , &true_purity         ) ;
  pidTree->Branch( "true_start_momentum"      , &true_start_momentum ) ;
  pidTree->Branch( "true_start_x"             , &true_start_x        ) ;
  pidTree->Branch( "true_start_y"             , &true_start_y        ) ;
  pidTree->Branch( "true_start_z"             , &true_start_z        ) ;
  pidTree->Branch( "true_end_momentum"        , &true_end_momentum   ) ;
  pidTree->Branch( "true_end_x"               , &true_end_x          ) ;
  pidTree->Branch( "true_end_y"               , &true_end_y          ) ;
  pidTree->Branch( "true_end_z"               , &true_end_z          ) ;
  pidTree->Branch( "track_id"                 , &track_id            ) ;
  pidTree->Branch( "track_start_x"            , &track_start_x       ) ;
  pidTree->Branch( "track_start_y"            , &track_start_y       ) ;
  pidTree->Branch( "track_start_z"            , &track_start_z       ) ;
  pidTree->Branch( "track_end_x"              , &track_end_x         ) ;
  pidTree->Branch( "track_end_y"              , &track_end_y         ) ;
  pidTree->Branch( "track_end_z"              , &track_end_z         ) ;
  pidTree->Branch( "track_likelihood_fwd_mu"  , &track_likelihood_fwd_mu    ) ;
  pidTree->Branch( "track_likelihood_fwd_p"   , &track_likelihood_fwd_p     ) ;
  pidTree->Branch( "track_likelihood_fwd_pi"  , &track_likelihood_fwd_pi    ) ;
  pidTree->Branch( "track_likelihood_fwd_k"   , &track_likelihood_fwd_k     ) ;
  pidTree->Branch( "track_likelihood_fwd_mip" , &track_likelihood_fwd_mip   ) ;
  pidTree->Branch( "track_likelihood_bwd_mu"  , &track_likelihood_bwd_mu    ) ;
  pidTree->Branch( "track_likelihood_bwd_p"   , &track_likelihood_bwd_p     ) ;
  pidTree->Branch( "track_likelihood_bwd_pi"  , &track_likelihood_bwd_pi    ) ;
  pidTree->Branch( "track_likelihood_bwd_k"   , &track_likelihood_bwd_k     ) ;
  pidTree->Branch( "track_likelihood_shift_fwd_mu"  , &track_likelihood_shift_fwd_mu    ) ;
  pidTree->Branch( "track_likelihood_shift_fwd_p"   , &track_likelihood_shift_fwd_p     ) ;
  pidTree->Branch( "track_likelihood_shift_fwd_pi"  , &track_likelihood_shift_fwd_pi    ) ;
  pidTree->Branch( "track_likelihood_shift_fwd_k"   , &track_likelihood_shift_fwd_k     ) ;
  pidTree->Branch( "track_likelihood_shift_bwd_mu"  , &track_likelihood_shift_bwd_mu    ) ;
  pidTree->Branch( "track_likelihood_shift_bwd_p"   , &track_likelihood_shift_bwd_p     ) ;
  pidTree->Branch( "track_likelihood_shift_bwd_pi"  , &track_likelihood_shift_bwd_pi    ) ;
  pidTree->Branch( "track_likelihood_shift_bwd_k"   , &track_likelihood_shift_bwd_k     ) ;
  pidTree->Branch( "track_PIDA_mean"          , &track_PIDA_mean          ) ;
  pidTree->Branch( "track_PIDA_median"        , &track_PIDA_median          ) ;
  pidTree->Branch( "track_PIDA_kde"           , &track_PIDA_kde          ) ;
  pidTree->Branch( "track_Chi2Proton"         , &track_Chi2Proton    ) ;
  pidTree->Branch( "track_Chi2Pion"           , &track_Chi2Pion      ) ;
  pidTree->Branch( "track_Chi2Kaon"           , &track_Chi2Kaon      ) ;
  pidTree->Branch( "track_Chi2Muon"           , &track_Chi2Muon      ) ;
  pidTree->Branch( "track_length"             , &track_length        ) ;
  pidTree->Branch( "track_dEdx"               , &track_dEdx          ) ;
  pidTree->Branch( "track_dQdx"               , &track_dQdx          ) ;
  pidTree->Branch( "track_dQdxtruncmeanvslength_isMuon"  , &track_dQdxtruncmeanvslength_isMuon ) ;
  pidTree->Branch( "track_theta"              , &track_theta         ) ;
  pidTree->Branch( "track_phi"                , &track_phi           ) ;
  pidTree->Branch( "track_nhits"              , &track_nhits         ) ;

  pidTree->Branch( "track_depE"              , &track_depE          ) ;
  pidTree->Branch( "track_rangeE_mu"         , &track_rangeE_mu     ) ;
  pidTree->Branch( "track_rangeE_p"          , &track_rangeE_p      ) ;
  pidTree->Branch( "track_dEdx_perhit_u"       , &track_dEdx_perhit_u   ) ;
  pidTree->Branch( "track_dEdx_perhit_v"       , &track_dEdx_perhit_v   ) ;
  pidTree->Branch( "track_dEdx_perhit_y"       , &track_dEdx_perhit_y   ) ;
  pidTree->Branch( "track_resrange_perhit_u"   , &track_resrange_perhit_u ) ;
  pidTree->Branch( "track_resrange_perhit_v"   , &track_resrange_perhit_v ) ;
  pidTree->Branch( "track_resrange_perhit_y"   , &track_resrange_perhit_y ) ;
  pidTree->Branch( "track_Lmip_perhit"         , &track_Lmip_perhit ) ;




  /**
   * Define array of labels for different particle species. We're going to use this
   * to produce a plot showing how often we identify each particles specied
   */
  const char* particles[5] = {"#mu", "p", "#pi", "K", "MIP"};
  All_chargeEndOverStart_sm0_5_dEdxrr  = tfs->make<TH2F>("All_chargeEndOverStart_sm0_5_dEdxrr"  , "All tracks (end/start average charge < 0.5);Residual range (cm); dEdx"                                                          , 150 , 0 , 30 , 400 , 0 , 50);
  All_chargeEndOverStart_gr2_dEdxrr    = tfs->make<TH2F>("All_chargeEndOverStart_gr2_dEdxrr"    , "All tracks (end/start average charge > 2);Residual range (cm); dEdx"                                                            , 150 , 0 , 30 , 400 , 0 , 50);
  All_chargeEndOverStart_0_5to2_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_0_5to2_dEdxrr" , "All tracks (end/start average charge                                                      = 0.5 - 2);Residual range (cm); dEdx" , 150 , 0 , 30 , 400 , 0 , 50);

  if (!fIsDataPlots){
    // ---- True Bragg peak
    TrueBragg_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("TrueBragg_chargeEndOverStart_directionCorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    TrueBragg_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("TrueBragg_chargeEndOverStart_directionIncorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    ;
    TrueBragg_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("TrueBragg_chargeEndOverStartVersusNHits_directionCorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 400, 0, 10, 6, 0, 6);
    TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 400, 0, 10, 6, 0, 6);
    ;

    TrueBragg_chargeEndOverStart_sm0_5_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_sm0_5_dEdxrr","Tracks with true p=0 at end (end/start average charge < 0.5);Residual range (cm); dEdx",150,0,30,400,0,50);
    TrueBragg_chargeEndOverStart_gr2_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_gr2_dEdxrr","Tracks with true p=0 at end (end/start average charge > 2);Residual range (cm); dEdx",150,0,30,400,0,50);
    TrueBragg_chargeEndOverStart_0_5to2_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_0_5to2_dEdxrr","Tracks with true p=0 at end (end/start average charge = 0.5 - 2);Residual range (cm); dEdx",150,0,30,400,0,50);

    TrueBragg_correctdirection = tfs->make<TH2F>("TrueBragg_correctdirection","Tracks with true p=0 at end, reconstructed correct direction;true particle;PID direction correct?",4,0,4,2,0,2);
    TrueBragg_incorrectdirection = tfs->make<TH2F>("TrueBragg_incorrectdirection","Tracks with true p=0 at end, reconstructed incorrect direction;true particle;PID direction correct?",4,0,4,2,0,2);
    for (size_t i=1; i<=4; i++){
      TrueBragg_correctdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
      TrueBragg_incorrectdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
    }
    TrueBragg_correctdirection->GetYaxis()->SetBinLabel(1,"False");
    TrueBragg_correctdirection->GetYaxis()->SetBinLabel(2,"True");
    TrueBragg_incorrectdirection->GetYaxis()->SetBinLabel(1,"False");
    TrueBragg_incorrectdirection->GetYaxis()->SetBinLabel(2,"True");

    All_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("All_chargeEndOverStart_directionCorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    All_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("All_chargeEndOverStart_directionIncorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    ;
    All_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("All_chargeEndOverStartVersusNHits_directionCorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 200, 0, 10, 6, 0, 6);
    All_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("All_chargeEndOverStartVersusNHits_directionIncorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};", 200, 0, 10, 6, 0, 6);
    ;

    Contained_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("Contained_chargeEndOverStart_directionCorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    Contained_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("Contained_chargeEndOverStart_directionIncorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    ;
    Contained_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("Contained_chargeEndOverStartVersusNHits_directionCorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 200, 0, 10, 6, 0, 6);
    Contained_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("Contained_chargeEndOverStartVersusNHits_directionIncorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};", 200, 0, 10, 6, 0, 6);
    ;

    All_correctdirection = tfs->make<TH2F>("All_correctdirection","All tracks, reconstructed correct direction;true particle;PID direction correct?",4,0,4,2,0,2);
    All_incorrectdirection = tfs->make<TH2F>("All_incorrectdirection","All tracks, reconstructed incorrect direction;true particle;PID direction correct?",4,0,4,2,0,2);
    for (size_t i=1; i<=4; i++){
      All_correctdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
      All_incorrectdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
    }
    All_correctdirection->GetYaxis()->SetBinLabel(1,"False");
    All_correctdirection->GetYaxis()->SetBinLabel(2,"True");
    All_incorrectdirection->GetYaxis()->SetBinLabel(1,"False");
    All_incorrectdirection->GetYaxis()->SetBinLabel(2,"True");
  }

}


// endSubRun function for MC POT counting
void ParticleIdValidationPlots::endSubRun(art::SubRun const &sr) {
  // Note: the entire subrun's POT is recorded in the tree for every event.
  // You must only add it once per subrun to get the correct number.

  art::Handle<sumdata::POTSummary> potsum_h;

  if (!isData) { // MC only (data is dealt with using Zarko's script)
    if(sr.getByLabel("generator", potsum_h)) {
      sr_pot = potsum_h->totpot;
    }
  }

  sr_run = sr.run();
  sr_sub_run = sr.subRun();

  potTree->Fill();

}


DEFINE_ART_MODULE(ParticleIdValidationPlots)

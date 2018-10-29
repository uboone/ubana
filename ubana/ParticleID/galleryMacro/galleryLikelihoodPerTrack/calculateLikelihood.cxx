#include "calculateLikelihood.h"

/**
 * calculate likelihood for each track compared to predicted likelihood
 */

const simb::MCParticle* getMatchedParticle(art::FindManyP<recob::Hit> trackHitAssn, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit, art::InputTag hitMatcherTag, size_t it){

  /* get MCParticle from a given track... */

  std::vector< art::Ptr<recob::Hit> > trkHitPtrs = trackHitAssn.at(it);
  std::unordered_map<int, double> trkide;
  double maxe=-1, tote=0;
  const simb::MCParticle* mcp;

  std::vector<simb::MCParticle const*> particleVec;
  std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

  // loop hits
  for (size_t ih = 0; ih < trkHitPtrs.size(); ih++){

    particleVec.clear();
    matchVec.clear();
    particlesPerHit.get(trkHitPtrs[ih].key(), particleVec, matchVec);

    // loop over particles
    for (size_t ip = 0; ip < particleVec.size(); ip++){

      trkide[particleVec[ip]->TrackId()]+=matchVec[ip]->energy;
      tote += matchVec[ip]->energy;
      if (trkide[particleVec[ip]->TrackId()] > maxe){
        maxe = trkide[ particleVec[ip]->TrackId()];
        mcp = particleVec[ip];

      }
    }
  }

  if (maxe < 0.0*tote)
    mcp = 0;

  return mcp;
}


bool isInTPC(double x, double y, double z){

  double DET_X_HIGH = 256.0;
  double DET_X_LOW = 0.0;
  double DET_Y_HIGH = 116.5;
  double DET_Y_LOW = -116.5;
  double DET_Z_HIGH = 1040;
  double DET_Z_LOW = 0.0;

  if (x > DET_X_LOW && x < DET_X_HIGH &&
      y > DET_Y_LOW && y < DET_Y_HIGH &&
      z > DET_Z_LOW && z < DET_Z_HIGH)
    return true;
  else return false;

}

bool isWellReconstructed(recob::Track track, simb::MCParticle mcp){

  double tsy = track.Start().Y();
  double tey = track.End().Y();
  double tsz = track.Start().Z();
  double tez = track.End().Z();

  double msy = mcp.Vy();
  double mey = mcp.EndY();
  double msz = mcp.Vz();
  double mez = mcp.EndZ();

  double twoDStartRes = std::sqrt(std::pow(msy-tsy,2)+std::pow(msz-tsz,2));
  double twoDStartResFlip = std::sqrt(std::pow(msy-tey,2)+std::pow(msz-tez,2));
  double twoDEndRes = std::sqrt(std::pow(mey-tey,2)+std::pow(mez-tez,2));
  double twoDEndResFlip = std::sqrt(std::pow(mey-tsy,2)+std::pow(mez-tsz,2));

  if ((twoDStartRes < 2.0 && twoDEndRes < 2.0) ||
      (twoDStartResFlip < 2.0 && twoDEndResFlip < 2.0)){
    std::cout << "[CALCLIKELIHOOD] Found a match!" << std::endl;
    std::cout << "[CALCLIKELIHOOD] Start Res Fwd : " << twoDStartRes << std::endl;
    std::cout << "[CALCLIKELIHOOD] End Res Fwd   : " << twoDEndRes << std::endl;
    std::cout << "[CALCLIKELIHOOD] Start Res Bwd : " << twoDStartResFlip << std::endl;
    std::cout << "[CALCLIKELIHOOD] End Res Bwd   : " << twoDEndResFlip << std::endl;
    return true;
  }
  else {
    /*    std::cout << "[CALCLIKELIHOOD] Start Res Fwd : " << twoDStartRes << std::endl;
          std::cout << "[CALCLIKELIHOOD] End Res Fwd   : " << twoDEndRes << std::endl;
          std::cout << "[CALCLIKELIHOOD] Start Res Bwd : " << twoDStartResFlip << std::endl;
          std::cout << "[CALCLIKELIHOOD] End Res Bwd   : " << twoDEndResFlip << std::endl;
          */    return false;
  }
}

int main(int argv, char** argc){

  std::vector<std::string> filenames;

  // extract input filename from input textfile
  std::string file_name;
  std::ifstream input_file(argc[1]);
  while (getline(input_file,file_name))
    filenames.push_back(file_name);

  double nevents;
  if (argc[2])
    nevents = std::atoi(argc[2]);
  else nevents = std::numeric_limits<int>::max();


  /**
   * This is just a placeholder parameter set to pass to the Bragg algorithm.
   * It doesn't need to be filled with anything, the algorithm configures
   * variables to be some default value.
   */
  fhicl::ParameterSet const p;
  particleid::Bragg_negLogL_Estimator algorithm;
  algorithm.configure(p);

  std::vector<double> fv;
  fidvol::fiducialVolume fid;
  fv = fid.setFiducialVolume(fv, p);
  fid.printFiducialVolume(fv);

  art::InputTag trackTag { "pandoraNu::McRecoStage2" };
  art::InputTag hitTag { "pandoraCosmicHitRemoval::McRecoStage2" };
  art::InputTag hitMatcherTag { "crHitRemovalTruthMatch::McRecoStage2" };
  art::InputTag caloTag { "pandoraNucali" }; // calibrated

  TFile *fBadLikelihood = new TFile("badLikelihood.root", "RECREATE");
  TFile *fGoodLikelihood = new TFile("goodLikelihood.root", "RECREATE");

  particleid::Theory_dEdx_resrange blah;
  TGraph* protonTheory = (TGraph*)blah.g_ThdEdxRR_Proton;
  TGraph* muonTheory   = (TGraph*)blah.g_ThdEdxRR_Muon;
  fBadLikelihood->cd();
  protonTheory->Write();
  muonTheory->Write();

  fGoodLikelihood->cd();
  protonTheory->Write();
  muonTheory->Write();

  int n_correct = 0;
  int n_total = 0;
  int eventCounter = 0;

  // begin event loop
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

    if (eventCounter > nevents) throw; 

    std::cout << "------ Processing "
      << "Run " << ev.eventAuxiliary().run() << ", "
      << "Event " << ev.eventAuxiliary().event() << ", "
      << "Time " << ev.eventAuxiliary().time().timeHigh() << "------" << std::endl;

    const auto& trackHandle = ev.getValidHandle< std::vector<recob::Track> >(trackTag);
    const auto& trackVec(*trackHandle);
    //const auto& mcpHandle = ev.getValidHandle< std::vector< simb::MCParticle > >("largeant");
    const auto& hitHandle = ev.getValidHandle< std::vector<recob::Hit> >(hitTag);

    art::FindManyP<recob::Hit> trackHitAssn(trackHandle, ev, trackTag);
    art::FindMany<anab::Calorimetry> trackCaloAssn(trackHandle, ev, caloTag);

    for (size_t it = 0; it < trackVec.size(); it++){

      // get tracks and mcparticles
      art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit(hitHandle, ev, hitMatcherTag);
      recob::Track track = trackVec.at(it);
      const simb::MCParticle* mcp = getMatchedParticle(trackHitAssn, particlesPerHit, hitMatcherTag, it);

      if (!mcp) continue;

      //if (!isInTPC(mcp.EndX(), mcp.EndY(), mcp.EndZ())) continue;

      std::vector< const anab::Calorimetry* > calos;
      trackCaloAssn.get(it, calos);

      int pdgCode = std::abs(mcp->PdgCode());

      std::cout << "[CALCLIKELIHOOD] Number of Calorimetry objects: " << calos.size() << std::endl;
      for (unsigned int i = 0; i < calos.size(); i++){

        const anab::Calorimetry* calo = calos.at(i);

        if (calo->PlaneID().Plane != 2) continue;

        const std::vector<double>& dedx = calo->dEdx();
        const std::vector<double>& resrg = calo->ResidualRange();

        double totalProtonLogLikelihoodBwd = algorithm.getNegLogL(dedx, resrg, 2212, false); 
        double totalMuonLogLikelihoodBwd   = algorithm.getNegLogL(dedx, resrg, 13, false);
        double totalPionLogLikelihoodBwd   = algorithm.getNegLogL(dedx, resrg, 211, false);
        double totalProtonLogLikelihoodFwd = algorithm.getNegLogL(dedx, resrg, 2212, true);
        double totalMuonLogLikelihoodFwd   = algorithm.getNegLogL(dedx, resrg, 13, true);
        double totalPionLogLikelihoodFwd   = algorithm.getNegLogL(dedx, resrg, 211, true);
        double totalMuonLogLikelihoodNoBragg = algorithm.getNegLogL(dedx, resrg, 0, true);

        double Bragg_mu = std::min(totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd);
        double Bragg_p  = std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd);
        double Bragg_pi = std::min(totalPionLogLikelihoodFwd, totalPionLogLikelihoodBwd);
        double Bragg_MipNoBragg = totalMuonLogLikelihoodNoBragg;

        std::cout << "[CALCLIKELIHOOD] true PDG: " << pdgCode << std::endl;
        std::cout << "[CALCLIKELIHOOD] proton likelihood: " << Bragg_p << std::endl;
        std::cout << "[CALCLIKELIHOOD] muon likelihood: " << Bragg_mu << std::endl;
        std::cout << "[CALCLIKELIHOOD] pion likelihood: " << Bragg_pi << std::endl;
        std::cout << "[CALCLIKELIHOOD] muon nobragg likelihood" << Bragg_MipNoBragg << std::endl;

        bool correctPDG = false;
        if (pdgCode == 2212){
          if (std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd) >= std::min(totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd)){
            correctPDG = false;
          }
          else{
            correctPDG = true;
          }
        }
        else if (pdgCode == 13 || pdgCode == -13 || pdgCode == 211 || pdgCode == -211)
        {
          if ((std::min(totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd) >= std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd))
              || (totalMuonLogLikelihoodNoBragg >= std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd))){
            correctPDG = false;
          }
          else{
            correctPDG = true;
          }
        }
        if (correctPDG) n_correct++;
        n_total++;

        if (pdgCode == 2212){

          double likelihoodProtonUnderMuonAssmp   = totalMuonLogLikelihoodFwd;
          double likelihoodProtonUnderProtonAssmp = totalProtonLogLikelihoodFwd;
          double likelihoodProtonUnderMuonBwdAssmp = totalMuonLogLikelihoodBwd;
          double likelihoodProtonUnderProtonBwdAssmp = totalProtonLogLikelihoodBwd;

          if (std::min(likelihoodProtonUnderProtonAssmp, likelihoodProtonUnderProtonBwdAssmp) >= std::min(likelihoodProtonUnderMuonAssmp, likelihoodProtonUnderMuonBwdAssmp)){

            fBadLikelihood->cd();
            TH2D* h = new TH2D(Form("h_ProtonMistakenForMuon_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodProtonUnderMuonAssmp < likelihoodProtonUnderMuonBwdAssmp)
                h->Fill(resrg.at(i), dedx.at(i));
              else
                h->Fill(resrg.at(resrg.size()-i-1), dedx.at(i));

            }
            h->Write();

          }
          else{

            fGoodLikelihood->cd();
            TH2D* h = new TH2D(Form("h_Proton_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodProtonUnderProtonAssmp < likelihoodProtonUnderProtonBwdAssmp)
                h->Fill(resrg.at(i), dedx.at(i));
              else
                h->Fill(resrg.at(resrg.size() - i -1), dedx.at(i));

            }
            h->Write();

          }

        }
        if (pdgCode == 13 || pdgCode == 211){

          double likelihoodMuonUnderMuonAssmp   = totalMuonLogLikelihoodFwd;
          double likelihoodMuonUnderProtonAssmp = totalProtonLogLikelihoodFwd;
          double likelihoodMuonUnderMuonBwdAssmp = totalMuonLogLikelihoodFwd;
          double likelihoodMuonUnderProtonBwdAssmp = totalMuonLogLikelihoodBwd;

          if (std::min(likelihoodMuonUnderProtonAssmp, likelihoodMuonUnderProtonBwdAssmp) <= std::min(likelihoodMuonUnderMuonAssmp, likelihoodMuonUnderMuonBwdAssmp)){

            fBadLikelihood->cd();
            TH2D* h = new TH2D(Form("h_MuonMistakenForProton_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodMuonUnderProtonAssmp < likelihoodMuonUnderProtonBwdAssmp)
                h->Fill(resrg.at(i), dedx.at(i));
              else
                h->Fill(resrg.at(resrg.size() - i -1), dedx.at(i));

            }
            h->Write();
          }
          else{

            fGoodLikelihood->cd();
            TH2D* h = new TH2D(Form("h_Muon_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodMuonUnderMuonAssmp < likelihoodMuonUnderMuonBwdAssmp)
                h->Fill(resrg.at(i), dedx.at(i));
              else
                h->Fill(resrg.at(resrg.size() - i -1), dedx.at(i));

            }
            h->Write();
          }

        }


        // Save plots for all tracks to root file (for debugging, for now)
        //        TH2D* h = new TH2D(Form("h_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

        //      for (size_t i = 0; i < dedx.size(); i++){

        // skip first and last hits in plots
        if (i == 0) continue;
        if (i == dedx.size()-1) continue;

        //          double min_neg2logl = std::min({totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd, totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd, totalPionLogLikelihoodFwd, totalPionLogLikelihoodBwd});
        /*
           double rr_shift;
           if (min_neg2logl == totalMuonLogLikelihoodFwd) rr_shift = rr_shift_mufwd;
           if (min_neg2logl == totalProtonLogLikelihoodFwd) rr_shift = rr_shift_pfwd;
           if (min_neg2logl == totalPionLogLikelihoodFwd) rr_shift = rr_shift_pifwd;
           if (min_neg2logl == totalMuonLogLikelihoodBwd) rr_shift = rr_shift_mubwd;
           if (min_neg2logl == totalProtonLogLikelihoodBwd) rr_shift = rr_shift_pbwd;
           if (min_neg2logl == totalPionLogLikelihoodBwd) rr_shift = rr_shift_pibwd;

           if ((min_neg2logl == totalMuonLogLikelihoodFwd) || (min_neg2logl == totalPionLogLikelihoodFwd) || (min_neg2logl == totalProtonLogLikelihoodFwd))
           h->Fill(resrg.at(i)+rr_shift, dedx.at(i));
           else
           h->Fill(resrg.at(resrg.size() - i -1)+rr_shift, dedx.at(i));

           }
           */
      }

    }
    eventCounter++;
  }

  std::cout << "Number of tracks correctly identified: " << n_correct << std::endl;
  std::cout << "Number of total tracks: " << n_total << std::endl;

  return 0;

}

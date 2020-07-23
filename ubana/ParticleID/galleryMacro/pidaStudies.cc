/*****************************************
 * Semi-toy MC for PID-related tests
 ****************************************/

// c++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

// root includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom3.h"

// art includes
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"

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


double getPida(std::vector<double> const resrg, std::vector<double> const dedx, double maxResRange, double minResRange, TString type, bool isDebug, bool isTrackFlipped, float offset, int pdg, int trackId, int run, int event, TFile* f_pidaVals, bool isSaveProtonPIDA){

  std::vector<float> pidaVals;

  for (size_t i = 0; i < resrg.size(); i++){

    if (resrg.at(i)+offset > maxResRange || resrg.at(i)+offset < minResRange || resrg.at(i)+offset < 0) continue;

    if (isDebug)
      std::cout << "true resrg " << resrg.at(i) << std::endl;

    float val;
    if (isTrackFlipped == 0){
      val = dedx.at(i)*std::pow(resrg.at(i)+offset, 0.42);
      if (isDebug) std::cout << "resrg: " << resrg.at(i)+offset << " dedx: " << dedx.at(i) << std::endl;
    }
    if (isTrackFlipped == 1){
      val = dedx.at(dedx.size()-i-1)*std::pow(resrg.at(i)+offset, 0.42);
      if (isDebug) std::cout << "resrg: " << resrg.at(i)+offset << " dedx: " << dedx.at(i) << std::endl;
      if (isDebug) std::cout << "resrg: " << resrg.at(i)+offset << " Fixed dedx: " << dedx.at(dedx.size()-i-1) << std::endl;
    }

    pidaVals.push_back(val);
  }

  if (pdg == 2212 && isSaveProtonPIDA == true && type == "mean"){


    f_pidaVals->cd();

    TString histName = Form("trackID%i_%i_%iproton", run, event, trackId);

    TH1D* h = new TH1D(histName, ";PIDA Vals;", 1000, 0, 50);

    for (size_t i = 0; i < pidaVals.size(); i++){

      h->Fill(pidaVals.at(i));

    }

    h->Write();
    h->Delete();
  }

  if (pidaVals.size() != 0 && type == "mean")
    return TMath::Mean(pidaVals.begin(), pidaVals.end());
  else if (pidaVals.size() != 0 && type == "median")
    return TMath::Median(pidaVals.size(), &pidaVals[0]);
  else{ 
    if (isDebug) std::cout << ">> resrg size: " << resrg.size() << std::endl;
    return -1;}

}

double getSegmentLength(double xprev, double yprev, double zprev, double x, double y, double z){

  return std::sqrt(std::pow((x-xprev),2) + std::pow((y-yprev),2) + std::pow((z-zprev),2));

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

double getTrackedLength(const simb::MCParticle* mcp){
  /* get tracked length */

  // get trajectory
  simb::MCTrajectory traj = mcp->Trajectory();
  double trackedLength=0.0;
  for (unsigned int i = 1; i < mcp->NumberTrajectoryPoints(); i++){

    double xprev = traj.X(i-1);
    double yprev = traj.Y(i-1);
    double zprev = traj.Z(i-1);
    double x = traj.X(i);
    double y = traj.Y(i);
    double z = traj.Z(i);

    if (isInTPC(x, y, z))
      trackedLength += getSegmentLength(xprev, yprev, zprev, x, y, z);
  }

  return trackedLength;

}

double getResolution(double mx, double my, double mz, double tx, double ty, double tz, double offset){

  double res = std::sqrt(std::pow(mx-tx,2) + std::pow(my-ty,2) + std::pow(mz-tz,2));
  return res+std::abs(offset);

}

int main(int argv, char** argc){

  bool isDebug = false;
  bool isSaveProtonPIDA = true;
  float gausMean = atof(argc[2]);
  float gausSigma = atof(argc[3]);
  float trackFlipCut = atof(argc[4]);

  std::vector<std::string> filenames;

  // extract input filename from input textfile
  std::string file_name;
  std::ifstream input_file(argc[1]);
  while (getline(input_file,file_name))
    filenames.push_back(file_name);

  art::InputTag trackTag { "pandoraNu" };
  art::InputTag hitTag { "gaushit" };
  art::InputTag hitMatcherTag { "gaushitTruthMatch" };
  art::InputTag caloTag { "pandoraNucalo" };

  TString saveString = Form("/uboone/data/users/alister1/ParticleId/output_gausMean_%0.2f_gausSigma%0.2f_trackFlipCut%0.2f.root", gausMean, gausSigma, trackFlipCut);
  TFile *f_pidaVals = new TFile("f_pidaVals.root", "UPDATE");
  TFile *fOutput = new TFile(saveString, "RECREATE");
  TH1D* meanPida = new TH1D("meanPida", ";PIDA value;", 101, -1, 30);
  TH1D* meanPidaMuon = new TH1D("meanPidaMuon", ";PIDA value;", 101, -1, 30);
  TH1D* meanPidaPion = new TH1D("meanPidaPion", ";PIDA value;", 101, -1, 30);
  TH1D* meanPidaProton = new TH1D("meanPidaProton", ";PIDA value;", 101, -1, 30);
  TH1D* meanPidaOther = new TH1D("meanPidaOther", ";PIDA value;", 101, -1, 30);
  TH1D* medianPida = new TH1D("medianPida", ";PIDA value;", 101, -1, 30);
  TH1D* medianPidaMuon = new TH1D("medianPidaMuon", ";PIDA value;", 101, -1, 30);
  TH1D* medianPidaPion = new TH1D("medianPidaPion", ";PIDA value;", 101, -1, 30);
  TH1D* medianPidaProton = new TH1D("medianPidaProton", ";PIDA value;", 101, -1, 30);
  TH1D* medianPidaOther = new TH1D("medianPidaOther", ";PIDA value;", 101, -1, 30);
  TH1D* startRes = new TH1D("startRes", ";Start Resolution (cm);", 100, 0, 20);
  TH1D* endRes = new TH1D("endRes", ";End Resolution (cm);", 100, 0, 20);
  TH1D* startResX = new TH1D("startResX", ";Start Resolution (cm);", 100, -10, 10);
  TH1D* endResX = new TH1D("endResX", ";End Resolution (cm);", 100, -10, 10);
  TH1D* startResY = new TH1D("startResY", ";Start Resolution (cm);", 100, -10, 10);
  TH1D* endResY = new TH1D("endResY", ";End Resolution (cm);", 100, -10, 10);
  TH1D* startResZ = new TH1D("startResZ", ";Start Resolution (cm);", 100, -10, 10);
  TH1D* endResZ = new TH1D("endResZ", ";End Resolution (cm);", 100, -10, 10);
  TH2D* endResPidaMuon = new TH2D("endResPidaMuon", ";End Resolution (cm); PIDA value", 50, 0, 10, 101, -1, 30);
  TH2D* endResPidaPion = new TH2D("endResPidaPion", ";End Resolution (cm); PIDA value", 50, 0, 10, 101, -1, 30);
  TH2D* endResPidaProton = new TH2D("endResPidaProton", ";End Resolution (cm); PIDA value", 50, 0, 10, 101, -1, 30);
  TH2D* endResPidaOther = new TH2D("endResPidaOther", ";End Resolution (cm); PIDA value", 50, 0, 10, 101, -1, 30);
  TH2D* trkThetaPidaMuon = new TH2D("trkThetaPidaMuon", ";Track Theta; PIDA value", 50, 0, 3.14, 101, -1, 30);
  TH2D* trkThetaPidaPion = new TH2D("trkThetaPidaPion", ";Track Theta; PIDA value", 50, 0, 3.14, 101, -1, 30);
  TH2D* trkThetaPidaProton = new TH2D("trkThetaPidaProton", ";Track Theta; PIDA value", 50, 0, 3.14, 101, -1, 30);
  TH2D* trkThetaPidaOther = new TH2D("trkThetaPidaOther", ";Track Theta; PIDA value", 50, 0, 3.14, 101, -1, 30);
  TH2D* trkPhiPidaMuon = new TH2D("trkPhiPidaMuon", ";Track Phi; PIDA value", 100, -3.14, 3.14, 101, -1, 30);
  TH2D* trkPhiPidaPion = new TH2D("trkPhiPidaPion", ";Track Phi; PIDA value", 100, -3.14, 3.14, 101, -1, 30);
  TH2D* trkPhiPidaProton = new TH2D("trkPhiPidaProton", ";Track Phi; PIDA value", 100, -3.14, 3.14, 101, -1, 30);
  TH2D* trkPhiPidaOther = new TH2D("trkPhiPidaOther", ";Track Phi; PIDA value", 100, -3.14, 3.14, 101, -1, 30);
  TH2D* trkLengthPidaMuon = new TH2D("trkLengthPidaMuon", ";Track Length; PIDA value", 100, 0, 100, 101, -1, 30);
  TH2D* trkLengthPidaPion = new TH2D("trkLengthPidaPion", ";Track Length; PIDA value", 100, 0, 100, 101, -1, 30);
  TH2D* trkLengthPidaProton = new TH2D("trkLengthPidaProton", ";Track Length; PIDA value", 100, 0, 100, 101, -1, 30);
  TH2D* trkLengthPidaOther = new TH2D("trkLengthPidaOther", ";Track Length; PIDA value", 100, 0, 100, 101, -1, 30);
  TH2D* mcpNoCaloYZStartPositions = new TH2D("mcpNoCaloYZStartPositions", ";Track Start Z (cm); Track Start Y (cm)", 50, 0, 1040, 25, -116.5, 116.5 );
  TH2D* mcpNoCaloYZEndPositions = new TH2D("mcpNoCaloYZEndPositions", ";Track End Z (cm); Track End Y (cm)", 50, 0, 1040, 25, -116.5, 116.5 );
  TH2D* trkNoCaloYZStartPositions = new TH2D("trkNoCaloYZStartPositions", ";Track Start Z (cm); Track Start Y (cm)", 50, 0, 1040, 25, -116.5, 116.5 );
  TH2D* trkNoCaloYZEndPositions = new TH2D("trkNoCaloYZEndPositions", ";Track End Z (cm); Track End Y (cm)", 50, 0, 1040, 25, -116.5, 116.5 );
  TH2D* trkNoCaloThetaPhi = new TH2D("trkNoCaloThetaPhi", ";Track theta; Track phi", 50, 0, 3.14, 100, -3.14, 3.14);
  TH1D* trkNoCaloLength = new TH1D("trkNoCaloLength", ";Track Length;", 100, 0, 100);
  TH2D* trkNoCaloLengthPhi = new TH2D("trkNoCaloLengthPhi", ";Phi; Length (cm)", 100, -3.14, 3.14, 100, 0, 100);
  TH2D* trkNoCaloLengthTheta = new TH2D("trkNoCaloLengthTheta", ";Theta; Length (cm)", 50, 0, 3.14, 100, 0, 100);
  TH2D* trkLowPidaYZStartPositions = new TH2D("trkLowPidaYZStartPositions", ";Track Start Z (cm); Track Start Y (cm)", 50, 0, 1040, 25, -116.5, 116.5 );
  TH2D* trkLowPidaYZEndPositions = new TH2D("trkLowPidaYZEndPositions", ";Track End Z (cm); Track End Y (cm)", 50, 0, 1040, 25, -116.5, 116.5 );
  TH2D* trkLowPidaThetaPhi = new TH2D("trkLowPidaThetaPhi", ";Track theta; Track phi", 50, 0, 3.14, 100, -3.14, 3.14);
  TH1D* trkLowPidaLength = new TH1D("trkLowPidaLength", ";Track Length;", 100, 0, 100);
  TH2D* trkLowPidaLengthPhi = new TH2D("trkLowPidaLengthPhi", ";Phi; Length (cm)", 100, -3.14, 3.14, 100, 0, 100);
  TH2D* trkLowPidaLengthTheta = new TH2D("trkLowPidaLengthTheta", ";Theta; Length (cm)", 50, 0, 3.14, 100, 0, 100);
  TH2D* trkLowPidaThetaxzThetayz = new TH2D("trkLowPidaThetaxzThetayz", ";Track theta xz; Track theta yz", 50, 0, 3.14, 100, 0, 3.14);


  // introduce randomness
  TRandom3 *r3 = new TRandom3();
  TF1* gaussian = new TF1("gaussian", "gaus", -100, 100);
  gaussian->SetParameters(1, gausMean, gausSigma);

  // begin event loop
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

//    if (isDebug)
      std::cout << "------ Processing "
        << "Run " << ev.eventAuxiliary().run() << ", " 
        << "Event " << ev.eventAuxiliary().event() << ", "
        << "Time " << ev.eventAuxiliary().time().timeHigh() << "------" << std::endl;

    const auto& trackHandle = ev.getValidHandle< std::vector<recob::Track> >(trackTag);
    const auto& trackVec(*trackHandle);
    const auto& hitHandle = ev.getValidHandle< std::vector<recob::Hit> >(hitTag);
    const auto& clusterHandle = ev.getValidHandle< std::vector<recob::Cluster> >(trackTag);
    const auto& clusterVec(*clusterHandle);

    art::FindManyP<recob::Hit> trackHitAssn(trackHandle, ev, trackTag);
    art::FindMany<anab::Calorimetry> trackCaloAssn(trackHandle, ev, caloTag);

    for (size_t it = 0; it < trackVec.size(); it++){

      if (isDebug) std::cout << "------ track " << it << " -------" << std::endl;

      // get tracks and mcparticles
      art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit(hitHandle, ev, hitMatcherTag);
      recob::Track track = trackVec.at(it);
      const simb::MCParticle* mcp = getMatchedParticle(trackHitAssn, particlesPerHit, hitMatcherTag, it);

      if (mcp == 0) continue;


      // get pdg code of matched mcp
      int mcpPdg = mcp->PdgCode();

      // here, check whether the mcp direction and the track direction matches.
      // if not, this indicates that the track has been flipped, so we flip it back
      double mcpStartX = mcp->Vx();
      double mcpStartY = mcp->Vy();
      double mcpStartZ = mcp->Vz();
      double mcpEndX   = mcp->EndX();
      double mcpEndY   = mcp->EndY();
      double mcpEndZ   = mcp->EndZ();
      double trkStartX = track.Trajectory().Start().X();
      double trkStartY = track.Trajectory().Start().Y();
      double trkStartZ = track.Trajectory().Start().Z();
      double trkEndX   = track.Trajectory().End().X();
      double trkEndY   = track.Trajectory().End().Y();
      double trkEndZ   = track.Trajectory().End().Z();
      double trkthetaxz = std::atan2(track.VertexDirection().X(), track.VertexDirection().Z());
      double trkthetayz = std::atan2(track.VertexDirection().Y(), track.VertexDirection().Z());

      if (mcp->EndE() > mcp->Mass() || mcp->EndProcess() != "FastScintillation" || !isInTPC(mcpEndX, mcpEndY, mcpEndZ) || !isInTPC(mcpEndX, mcpEndY, mcpEndZ)){
        std::cout << "------ MCP DID NOT PASS ------" << std::endl;
        std::cout << ">> End process:      " << mcp->EndProcess() << std::endl;
        std::cout << ">> particle with pdg:" << mcpPdg << std::endl;
        std::cout << ">> EndE:             " << mcp->EndE() << std::endl;
        std::cout << ">> Mass:             " << mcp->Mass() << std::endl;
        std::cout << ">> mcpStart X, Y, Z: " << mcpStartX << " " << mcpStartY << " " << mcpStartZ << std::endl;
        std::cout << ">> mcpEnd X, Y, Z:   " << mcpEndX << " " << mcpEndY << " " << mcpEndZ << std::endl;
        continue;
      }
      else {
        std::cout << "------ MCP PASSED ------" << std::endl;
        std::cout << ">> End process:      " << mcp->EndProcess() << std::endl;
        std::cout << ">> particle with pdg:" << mcpPdg << std::endl;
        std::cout << ">> EndE:             " << mcp->EndE() << std::endl;
        std::cout << ">> Mass:             " << mcp->Mass() << std::endl;
        std::cout << ">> mcpStart X, Y, Z: " << mcpStartX << " " << mcpStartY << " " << mcpStartZ << std::endl;
        std::cout << ">> mcpEnd X, Y, Z:   " << mcpEndX << " " << mcpEndY << " " << mcpEndZ << std::endl;
      }


      bool isTrackFlipped;

      if (trkStartZ > trkEndZ && mcpStartZ < mcpEndZ)
        isTrackFlipped = true;
      else if (trkStartZ < trkEndZ && mcpStartZ > mcpEndZ)
        isTrackFlipped = true;
      else isTrackFlipped = false;

      if (isDebug){ 
        std::cout << ">> MCP start Z: " << mcpStartZ << " MCP End Z: " << mcpEndZ << std::endl;
        std::cout << ">> Track start Z: " << trkStartZ << "Track end Z: " << trkEndZ << std::endl;
        std::cout << ">> Is reco track flipped with respect to MCP? " << isTrackFlipped << std::endl;
      }

      // get tracking offset, i.e. how much of the end of the track do we miss?
      // allows us to simulate the possibility that track end resolution is worse
      // in data than in simulation
      float offset = gaussian->GetRandom();
      double trackLength_reco = track.Length();
      double trackLength_true = getTrackedLength(mcp);

      if (isDebug){
        std::cout << ">> trackLength_reco: " << trackLength_reco << " trackLength_true: " << trackLength_true << std::endl;
        std::cout << ">> Offset: " << offset << std::endl;
      }

      trackLength_reco = trackLength_reco+offset;

      double startResolution = -1;
      double endResolution = -1;

      if (isTrackFlipped == 0){
        startResolution = getResolution(mcpStartX, mcpStartY, mcpStartZ, trkStartX, trkStartY, trkStartZ, 0);
        endResolution = getResolution(mcpEndX, mcpEndY, mcpEndZ, trkEndX, trkEndY, trkEndZ, offset);
        startResX->Fill(mcpStartX-trkStartX);
        startResY->Fill(mcpStartY-trkStartY);
        startResZ->Fill(mcpStartZ-trkStartZ);
        endResX->Fill(mcpEndX-trkEndX);
        endResY->Fill(mcpEndY-trkEndY);
        endResZ->Fill(mcpEndZ-trkEndZ);
      }
      else if (isTrackFlipped == 1){
        startResolution = getResolution(mcpStartX, mcpStartY, mcpStartZ, trkEndX, trkEndY, trkEndZ, 0);
        endResolution = getResolution(mcpEndX, mcpEndY, mcpEndZ, trkStartX, trkStartY, trkStartZ, offset);
        startResX->Fill(mcpStartX-trkEndX);
        startResY->Fill(mcpStartY-trkEndY);
        startResZ->Fill(mcpStartZ-trkEndZ);
        endResX->Fill(mcpEndX-trkStartX);
        endResY->Fill(mcpEndY-trkStartY);
        endResZ->Fill(mcpEndZ-trkStartZ);

      }
      startRes->Fill(startResolution);
      endRes->Fill(endResolution);

      if (startResolution > 5.0 || endResolution > 5.0) continue;

      // now throw a random number to determine whether to flip the track or not
      // allows us to simulate larger fraction of tracks being flipped in data 
      // than in simulation
      float rando = r3->Uniform(0,1);
      if (rando < trackFlipCut) isTrackFlipped = !isTrackFlipped;
      if (isDebug) std::cout << ">> is Track Flipped? " << isTrackFlipped << std::endl; 

      std::vector<const anab::Calorimetry*> calos;
      std::vector< art::Ptr<recob::Hit> > hits = trackHitAssn.at(it);
      trackCaloAssn.get(it, calos);

      int hitCounter = 0;
      for (size_t i = 0; i < hits.size(); i++){

        if (hits.at(i)->View() == 2)
          hitCounter++;

      }

      auto const& calo = calos.at(2); // ONLY TAKE COLLECTION PLANE
      const std::vector<double>& dedx = calo->dEdx();
      const std::vector<double>& resrg = calo->ResidualRange();

      std::cout << ">> number of calos:      " << dedx.size() << std::endl; 
      std::cout << ">> number of hits:      " << hitCounter << std::endl; 

      for (size_t i = 0; i < clusterVec.size(); i++){

        std::cout << ">> cluster " << i << " has " << clusterVec.at(i).NHits() << " hits" << std::endl;

      }

      // if the track is flipped then we want to take the start and start+30cm instead of
      // residual range

      double maxResRg = 30 + offset;
      double minResRg = 0 + offset;

      double pidaVal;
      pidaVal= getPida(resrg, dedx, maxResRg, minResRg, "mean", isDebug, isTrackFlipped, offset, mcpPdg, track.ID(), ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), f_pidaVals, isSaveProtonPIDA); 
      meanPida->Fill(pidaVal);
      if (mcpPdg == 13)   meanPidaMuon->Fill(pidaVal);
      else if (mcpPdg == 211)  meanPidaPion->Fill(pidaVal);
      else if (mcpPdg == 2212) meanPidaProton->Fill(pidaVal);
      else meanPidaOther->Fill(pidaVal);

      pidaVal = getPida(resrg, dedx, maxResRg, minResRg, "median", isDebug, isTrackFlipped, offset, mcpPdg, track.ID(), ev.eventAuxiliary().run(), ev.eventAuxiliary().run(), f_pidaVals, isSaveProtonPIDA); 
      medianPida->Fill(pidaVal);

      if (mcpPdg == 13){
        medianPidaMuon->Fill(pidaVal);
        endResPidaMuon->Fill(endResolution, pidaVal);
        trkThetaPidaMuon->Fill(track.Theta(), pidaVal);
        trkPhiPidaMuon->Fill(track.Phi(), pidaVal);
        trkLengthPidaMuon->Fill(track.Length(), pidaVal);
      }
      else if (mcpPdg == 211){
        medianPidaPion->Fill(pidaVal);
        endResPidaPion->Fill(endResolution, pidaVal);
        trkThetaPidaPion->Fill(track.Theta(), pidaVal);
        trkPhiPidaPion->Fill(track.Phi(), pidaVal);
        trkLengthPidaPion->Fill(track.Length(), pidaVal);
      }
      else if (mcpPdg == 2212){
        medianPidaProton->Fill(pidaVal);
        endResPidaProton->Fill(endResolution, pidaVal);
        trkThetaPidaProton->Fill(track.Theta(), pidaVal);
        trkPhiPidaProton->Fill(track.Phi(), pidaVal);
        trkLengthPidaProton->Fill(track.Length(), pidaVal);
      }
      else{
        medianPidaOther->Fill(pidaVal);
        endResPidaOther->Fill(endResolution, pidaVal);
        trkThetaPidaOther->Fill(track.Theta(), pidaVal);
        trkPhiPidaOther->Fill(track.Phi(), pidaVal);
        trkLengthPidaOther->Fill(track.Length(), pidaVal);
      }

      if (pidaVal == -1){
        mcpNoCaloYZStartPositions->Fill(mcpStartZ, mcpStartY);
        mcpNoCaloYZEndPositions->Fill(mcpEndZ, mcpEndY);
        trkNoCaloThetaPhi->Fill(track.Theta(), track.Phi());
        trkNoCaloLength->Fill(trackLength_reco);
        trkNoCaloLengthPhi->Fill(track.Phi(), trackLength_reco);
        trkNoCaloLengthTheta->Fill(track.Theta(), trackLength_reco);
        if (isTrackFlipped == 0){
          trkNoCaloYZStartPositions->Fill(trkStartZ, trkStartY);
          trkNoCaloYZEndPositions->Fill(trkEndZ, trkEndY);
        }
        if (isTrackFlipped == 1){
          trkNoCaloYZStartPositions->Fill(trkEndZ, trkEndY);
          trkNoCaloYZEndPositions->Fill(trkStartZ, trkStartY);
        }
      }
      if (pidaVal > 0 && pidaVal < 3){
        trkLowPidaThetaPhi->Fill(track.Theta(), track.Phi());
        trkLowPidaLength->Fill(trackLength_reco);
        trkLowPidaLengthPhi->Fill(track.Phi(), trackLength_reco);
        trkLowPidaLengthTheta->Fill(track.Theta(), trackLength_reco);
        trkLowPidaThetaxzThetayz->Fill(trkthetaxz, trkthetayz);
        if (isTrackFlipped == 0){
          trkLowPidaYZStartPositions->Fill(trkStartZ, trkStartY);
          trkLowPidaYZEndPositions->Fill(trkEndZ, trkEndY);
        }
        if (isTrackFlipped == 1){
          trkLowPidaYZStartPositions->Fill(trkEndZ, trkEndY);
          trkLowPidaYZEndPositions->Fill(trkStartZ, trkStartY);
        }
      }


    }
  }

  std::cout << ">>| Track starting resolution mean: " << startRes->GetMean() << " and sigma: " << startRes->GetStdDev() << std::endl;
  std::cout << ">>|>>| Track starting resolution X mean: " << startResX->GetMean() << " and sigma: " << startResX->GetStdDev() << std::endl;
  std::cout << ">>|>>| Track starting resolution Y mean: " << startResY->GetMean() << " and sigma: " << startResY->GetStdDev() << std::endl;
  std::cout << ">>|>>| Track starting resolution Z mean: " << startResZ->GetMean() << " and sigma: " << startResZ->GetStdDev() << std::endl;
  std::cout << ">>| Track ending resolution mean: " << endRes->GetMean() << " and sigma: " << endRes->GetStdDev() << std::endl;
  std::cout << ">>|>>| Track ending resolution X mean: " << endResX->GetMean() << " and sigma: " << endResX->GetStdDev() << std::endl;
  std::cout << ">>|>>| Track ending resolution Y mean: " << endResY->GetMean() << " and sigma: " << endResY->GetStdDev() << std::endl;
  std::cout << ">>|>>| Track ending resolution Z mean: " << endResZ->GetMean() << " and sigma: " << endResZ->GetStdDev() << std::endl;

  fOutput->cd();
  meanPida->Write();
  meanPidaMuon->Write();
  meanPidaPion->Write();
  meanPidaProton->Write();
  meanPidaOther->Write();
  medianPida->Write();
  medianPidaMuon->Write();
  medianPidaPion->Write();
  medianPidaProton->Write();
  medianPidaOther->Write();
  startRes->Write();
  startResX->Write();
  startResY->Write();
  startResZ->Write();
  endRes->Write();
  endResX->Write();
  endResY->Write();
  endResZ->Write();
  endResPidaMuon->Write();
  endResPidaPion->Write();
  endResPidaProton->Write();
  endResPidaOther->Write();
  trkThetaPidaMuon->Write();
  trkThetaPidaPion->Write();
  trkThetaPidaProton->Write();
  trkThetaPidaOther->Write();
  trkPhiPidaMuon->Write();
  trkPhiPidaPion->Write();
  trkPhiPidaProton->Write();
  trkPhiPidaOther->Write();
  trkLengthPidaMuon->Write();
  trkLengthPidaPion->Write();
  trkLengthPidaProton->Write();
  trkLengthPidaOther->Write();
  mcpNoCaloYZStartPositions->Write();
  mcpNoCaloYZEndPositions->Write();
  trkNoCaloYZStartPositions->Write();
  trkNoCaloYZEndPositions->Write();
  trkNoCaloThetaPhi->Write();
  trkNoCaloLength->Write();
  trkNoCaloLengthPhi->Write();
  trkNoCaloLengthTheta->Write();
  trkLowPidaYZStartPositions->Write();
  trkLowPidaYZEndPositions->Write();
  trkLowPidaThetaPhi->Write();
  trkLowPidaLength->Write();
  trkLowPidaLengthPhi->Write();
  trkLowPidaLengthTheta->Write();
  trkLowPidaThetaxzThetayz->Write();
  return 0;

}

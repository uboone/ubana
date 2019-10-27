// _________________________________________________________________________________________________________________________________________________________________________________________________

#ifndef PID_CXX
#define PID_CXX

#include "PID.h"
// _________________________________________________________________________________________________________________________________________________________________________________________________

void PID::Chi2(art::FindManyP<anab::ParticleID> PIDTotrackAsso, art::Ptr<recob::Track> track, TVector3 Trk_start_SCEcorr, TVector3 Trk_end_SCEcorr, int hits_dEdx_size_pl0, int hits_dEdx_size_pl1, int hits_dEdx_size_pl2){
  // Gain PID info of the track
  if(!PIDTotrackAsso.isValid()){
    throw cet::exception("[Numu0pi0p]") << "No matched PID - track information!" << std::endl;
  }
  auto Trk_pos = Trk_end_SCEcorr - Trk_start_SCEcorr;

  double theta_pl2 = std::atan2(Trk_pos.Z(), Trk_pos.Y()); // atan2(y,x)
  double theta_pl1 = theta_pl2 + M_PI/3; // If plan1 is -60 degree to Y, looking from outside to the TPC
  double theta_pl0 = theta_pl2 - M_PI/3; // If plan0 is +60 degree to Y, looking from outside to the TPC

  int w2 = 0; int w1 = 0; int w0 = 0;
  double sin2_pl2 = sin(theta_pl2) * sin(theta_pl2);
  double sin2_pl1 = sin(theta_pl1) * sin(theta_pl1);
  double sin2_pl0 = sin(theta_pl0) * sin(theta_pl0);
  
  // Get projected angle wrt to the wires (docdb 23008)
  if (sin2_pl2 >= 0.5) w2 = 1;
  if (sin2_pl1 >= 0.5) w1 = 1;
  if (sin2_pl0 >= 0.5) w0 = 1;

  // PID
  auto trkPID = PIDTotrackAsso.at(track.key());
  if (trkPID.size() == 0){
    std::cout << "No PID information for this selected track!" << std::endl;
  }
  std::vector<anab::sParticleIDAlgScores> vAlg_PID = trkPID.front()->ParticleIDAlgScores();
  double PIDChi2_mu[3] = {-999,-999,-999};
  double PIDChi2_p[3] = {-999,-999,-999};
  double PIDChi2_pi[3] = {-999,-999,-999};
  double PIDChi2_K[3] = {-999,-999,-999};
  for(int i_Alg_PID = 0; i_Alg_PID < (int) vAlg_PID.size(); i_Alg_PID++){
    anab::sParticleIDAlgScores Alg_PID = vAlg_PID.at(i_Alg_PID);
    for(int id_pl = 0; id_pl < 3; id_pl++){
      if (Alg_PID.fPlaneMask.test(id_pl) && Alg_PID.fAlgName == "Chi2"){
         if (abs(Alg_PID.fAssumedPdg) == 13){ // muon
           PIDChi2_mu[id_pl] = Alg_PID.fValue;
         }
         if (abs(Alg_PID.fAssumedPdg) == 2212){ // proton
           PIDChi2_p[id_pl] = Alg_PID.fValue;
         }
         if (abs(Alg_PID.fAssumedPdg) == 211){ // pion
           PIDChi2_pi[id_pl] = Alg_PID.fValue;
         }
         if (abs(Alg_PID.fAssumedPdg) == 321){ // kaon
           PIDChi2_K[id_pl] = Alg_PID.fValue;
         }
      }
    }
  }

   PID_Chi2Mu_pl0 = PIDChi2_mu[0];
   PID_Chi2Mu_pl1 = PIDChi2_mu[1];
   PID_Chi2Mu_pl2 = PIDChi2_mu[2];
   int ww0 = w0; int ww1 = w1; int ww2 = w2; // copy from the origin
   if (PID_Chi2Mu_pl0 < 0) ww0 = 0;
   if (PID_Chi2Mu_pl1 < 0) ww1 = 0;
   if (PID_Chi2Mu_pl2 < 0) ww2 = 0;
   if (ww0 + ww1 + ww2 == 0) PID_Chi2Mu_3pl = -999;
   else PID_Chi2Mu_3pl = (ww0 * PID_Chi2Mu_pl0 + ww1 * PID_Chi2Mu_pl1 + ww2 * PID_Chi2Mu_pl2) / (ww0 + ww1 + ww2);

   PID_Chi2P_pl0 = PIDChi2_p[0];
   PID_Chi2P_pl1 = PIDChi2_p[1];
   PID_Chi2P_pl2 = PIDChi2_p[2];
   ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
   if (PID_Chi2P_pl0 < 0) ww0 = 0;
   if (PID_Chi2P_pl1 < 0) ww1 = 0;
   if (PID_Chi2P_pl2 < 0) ww2 = 0;
   if (ww0 + ww1 + ww2 == 0) PID_Chi2P_3pl = -999;
   else PID_Chi2P_3pl = (ww0 * PID_Chi2P_pl0 + ww1 * PID_Chi2P_pl1 + ww2 * PID_Chi2P_pl2) / (ww0 + ww1 + ww2);

   PID_Chi2Pi_pl0 = PIDChi2_pi[0];
   PID_Chi2Pi_pl1 = PIDChi2_pi[1];
   PID_Chi2Pi_pl2 = PIDChi2_pi[2];
   ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
   if (PID_Chi2Pi_pl0 < 0) ww0 = 0;
   if (PID_Chi2Pi_pl1 < 0) ww1 = 0;
   if (PID_Chi2Pi_pl2 < 0) ww2 = 0;
   if (ww0 + ww1 + ww2 == 0) PID_Chi2Pi_3pl = -999;
   else PID_Chi2Pi_3pl = (ww0 * PID_Chi2Pi_pl0 + ww1 * PID_Chi2Pi_pl1 + ww2 * PID_Chi2Pi_pl2) / (ww0 + ww1 + ww2);

   PID_Chi2K_pl0 = PIDChi2_K[0];
   PID_Chi2K_pl1 = PIDChi2_K[1];
   PID_Chi2K_pl2 = PIDChi2_K[2];
   ww0 = w0; ww1 = w1; ww2 = w2; // copy from the origin
   if (PID_Chi2K_pl0 < 0) ww0 = 0;
   if (PID_Chi2K_pl1 < 0) ww1 = 0;
   if (PID_Chi2K_pl2 < 0) ww2 = 0;
   if (ww0 + ww1 + ww2 == 0) PID_Chi2K_3pl = -999;
   else PID_Chi2K_3pl = (ww0 * PID_Chi2K_pl0 + ww1 * PID_Chi2K_pl1 + ww2 * PID_Chi2K_pl2) / (ww0 + ww1 + ww2);

   std::vector<int> Nhit_3pl = {hits_dEdx_size_pl0, hits_dEdx_size_pl1, hits_dEdx_size_pl2};
   BestPlane_PID = std::max_element(Nhit_3pl.begin(), Nhit_3pl.end()) - Nhit_3pl.begin();
   if (BestPlane_PID == 0){
     Pl0_for_PID = true;
     Pl1_for_PID = false;
     Pl2_for_PID = false;
   }
   if (BestPlane_PID == 1){
     Pl0_for_PID = false;
     Pl1_for_PID = true;
     Pl2_for_PID = false;
   }
   if (BestPlane_PID == 2){
     Pl0_for_PID = false;
     Pl1_for_PID = false;
     Pl2_for_PID = true;
   }

   //-- Get minimum Chi2 and there corresponding particle type
   std::vector<double> PIDChi2_avg;// It follows the order of muon, proton, pion, kaon
   if (PID_Chi2Mu_3pl < 0) PIDChi2_avg.push_back(9999);
   else PIDChi2_avg.push_back(PID_Chi2Mu_3pl);
   if (PID_Chi2P_3pl < 0) PIDChi2_avg.push_back(9999);
   else PIDChi2_avg.push_back(PID_Chi2P_3pl);
   if (PID_Chi2Pi_3pl < 0) PIDChi2_avg.push_back(9999);
   else PIDChi2_avg.push_back(PID_Chi2Pi_3pl);
   if (PID_Chi2K_3pl < 0) PIDChi2_avg.push_back(9999);
   else PIDChi2_avg.push_back(PID_Chi2K_3pl);

   PID_avg_Chi2 = *std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end());
   int ID_PID = std::min_element(PIDChi2_avg.begin(), PIDChi2_avg.end()) - PIDChi2_avg.begin();
   if (ID_PID == 0) PID_Pdg_3pl = 13;
   if (ID_PID == 1) PID_Pdg_3pl = 2212;
   if (ID_PID == 2) PID_Pdg_3pl = 211;
   if (ID_PID == 3) PID_Pdg_3pl = 321;

   // Use Plane 2
   std::vector<double> PIDChi2_pl2;// It follows the order of muon, proton, pion, kaon
   if (PIDChi2_mu[2] < 0) PIDChi2_pl2.push_back(9999);
   else PIDChi2_pl2.push_back(PIDChi2_mu[2]);
   if (PIDChi2_p[2] < 0) PIDChi2_pl2.push_back(9999);
   else PIDChi2_pl2.push_back(PIDChi2_p[2]);
   if (PIDChi2_pi[2] < 0) PIDChi2_pl2.push_back(9999);
   else PIDChi2_pl2.push_back(PIDChi2_pi[2]);
   if (PIDChi2_K[2] < 0) PIDChi2_pl2.push_back(9999);
   else PIDChi2_pl2.push_back(PIDChi2_K[2]);

   PID_pl2_Chi2 = *std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end());
   int ID_PID_pl2 = std::min_element(PIDChi2_pl2.begin(), PIDChi2_pl2.end()) - PIDChi2_pl2.begin();
   if (ID_PID_pl2 == 0) PID_Pdg_pl2 = 13;
   if (ID_PID_pl2 == 1) PID_Pdg_pl2 = 2212;
   if (ID_PID_pl2 == 2) PID_Pdg_pl2 = 211;
   if (ID_PID_pl2 == 3) PID_Pdg_pl2 = 321;

   // Use plane 1
   std::vector<double> PIDChi2_pl1;// It follows the order of muon, proton, pion, kaon
   if (PIDChi2_mu[1] < 0) PIDChi2_pl1.push_back(9999);
   else PIDChi2_pl1.push_back(PIDChi2_mu[1]);
   if (PIDChi2_p[1] < 0) PIDChi2_pl1.push_back(9999);
   else PIDChi2_pl1.push_back(PIDChi2_p[1]);
   if (PIDChi2_pi[1] < 0) PIDChi2_pl1.push_back(9999);
   else PIDChi2_pl1.push_back(PIDChi2_pi[1]);
   if (PIDChi2_K[1] < 0) PIDChi2_pl1.push_back(9999);
   else PIDChi2_pl1.push_back(PIDChi2_K[1]);

   PID_pl1_Chi2 = *std::min_element(PIDChi2_pl1.begin(), PIDChi2_pl1.end());
   int ID_PID_pl1 = std::min_element(PIDChi2_pl1.begin(), PIDChi2_pl1.end()) - PIDChi2_pl1.begin();
   if (ID_PID_pl1 == 0) PID_Pdg_pl1 = 13;
   if (ID_PID_pl1 == 1) PID_Pdg_pl1 = 2212;
   if (ID_PID_pl1 == 2) PID_Pdg_pl1 = 211;
   if (ID_PID_pl1 == 3) PID_Pdg_pl1 = 321;

   // Use plane 0
   std::vector<double> PIDChi2_pl0;// It follows the order of muon, proton, pion, kaon
   if (PIDChi2_mu[0] < 0) PIDChi2_pl0.push_back(9999);
   else PIDChi2_pl0.push_back(PIDChi2_mu[0]);
   if (PIDChi2_p[0] < 0) PIDChi2_pl0.push_back(9999);
   else PIDChi2_pl0.push_back(PIDChi2_p[0]);
   if (PIDChi2_pi[0] < 0) PIDChi2_pl0.push_back(9999);
   else PIDChi2_pl0.push_back(PIDChi2_pi[0]);
   if (PIDChi2_K[0] < 0) PIDChi2_pl0.push_back(9999);
   else PIDChi2_pl0.push_back(PIDChi2_K[0]);

   PID_pl0_Chi2 = *std::min_element(PIDChi2_pl0.begin(), PIDChi2_pl0.end());
   int ID_PID_pl0 = std::min_element(PIDChi2_pl0.begin(), PIDChi2_pl0.end()) - PIDChi2_pl0.begin();
   if (ID_PID_pl0 == 0) PID_Pdg_pl0 = 13;
   if (ID_PID_pl0 == 1) PID_Pdg_pl0 = 2212;
   if (ID_PID_pl0 == 2) PID_Pdg_pl0 = 211;
   if (ID_PID_pl0 == 3) PID_Pdg_pl0 = 321;
}


#endif

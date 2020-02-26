// _________________________________________________________________________________________________________________________________________________________________________________________________

#ifndef PID_CXX
#define PID_CXX

#include "PID.h"
// _________________________________________________________________________________________________________________________________________________________________________________________________

int PID::Chi2(anab::Calorimetry){
  //Calorimetry Info
  // pandoracaliSCE has E-field and spatial correction
  if(assoCal.size()!=3){
    throw cet::exception("[Numu0pi0p]") << "Where are the three planes for the calorimetry!" << std::endl;
  }
  // induction 0 = 0, induction 1 = 1, collection = 2 (Check if the plane ID is correct)
  // assoCal[id_pl]->PlaneID().Plane == 2 (collection)
  // The vector is ordered by residual range from small to big (track end to track start)
  dEdx_pl0 = assoCal[0]->dEdx();
  dQdx_pl0 = assoCal[0]->dQdx();
  resRange_pl0 = assoCal[0]->ResidualRange();

  dEdx_pl1 = assoCal[1]->dEdx();
  dQdx_pl1 = assoCal[1]->dQdx();
  resRange_pl1 = assoCal[1]->ResidualRange();

  dEdx_pl2 = assoCal[2]->dEdx();
  dQdx_pl2 = assoCal[2]->dQdx();
  resRange_pl2 = assoCal[2]->ResidualRange();

  hits_dEdx_size_pl0 = dEdx_pl0.size();
  hits_dEdx_size_pl1 = dEdx_pl1.size();
  hits_dEdx_size_pl2 = dEdx_pl2.size();

  int half_size_pl0 = hits_dEdx_size_pl0 / 2;
  int half_size_pl1 = hits_dEdx_size_pl1 / 2;
  int half_size_pl2 = hits_dEdx_size_pl2 / 2;

  dEdx_pl0_start_half = std::accumulate(dEdx_pl0.end() - half_size_pl0, dEdx_pl0.end(), 0.) / half_size_pl0;
  dEdx_pl0_end_half = std::accumulate(dEdx_pl0.begin(), dEdx_pl0.begin() + half_size_pl0, 0. ) / half_size_pl0;
  dEdx_pl1_start_half = std::accumulate(dEdx_pl1.end() - half_size_pl1, dEdx_pl1.end(), 0.) / half_size_pl1;
  dEdx_pl1_end_half = std::accumulate(dEdx_pl1.begin(), dEdx_pl1.begin() + half_size_pl1, 0. ) / half_size_pl1;
  dEdx_pl2_start_half = std::accumulate(dEdx_pl2.end() - half_size_pl2, dEdx_pl2.end(), 0.) / half_size_pl2;
  dEdx_pl2_end_half = std::accumulate(dEdx_pl2.begin(), dEdx_pl2.begin() + half_size_pl2, 0. ) / half_size_pl2;
}


#endif

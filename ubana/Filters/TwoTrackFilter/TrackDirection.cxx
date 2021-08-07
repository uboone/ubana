//Class created by Yifan Chen (yifan.chen@lhep.unibe.ch)
//July 3rd, 2019

// _________________________________________________________________________________________________________________________________________________________________________________________________

#ifndef TRACKDIRECTION_CXX
#define TRACKDIRECTION_CXX

#include "TrackDirection.h"
// _________________________________________________________________________________________________________________________________________________________________________________________________

//1. If a pandora track contains multiple track: boolean (if dE/dx is not continuous)
//2. Where's the break point in terms the residual range: double (as the direction and vertex hasn't been moved)
//3. Where's the break point in terms the X, Y, Z: TVector3 
///// For non-merged track
//4. Directional judgement by dE/dx: boolean (if it's the same as pandora track direction)
//5. Directional judgement by MCS: boolean (if it's the same as pandora track direction)
//6. Directional judgement by combined: boolean (if it's the same as pandora track direction)
//     What if they are not consistent?
//7. Vertex: TVector3
//8. Angle theta: double
//9. Angle phi: double
//// If the vertex is in the TPC
//10. If there's a potential broken track which is prior to the vertex: boolean
//      The angle of the sub-tracks; The continuity of the sub-tracks; The direction of the sub-tracks
//11.  Is the potential new vertex in the TPC: boolean
//12.  Where's the potential new vertex: TVector3

void TrackDirection::TrackDir(art::Ptr<recob::Track>& ThisTrack, art::FindManyP<anab::Calorimetry>& trkToCalAsso){
  auto assoCal = trkToCalAsso.at(ThisTrack.key()); // take the only track
  std::cout<<"calo size: "<< assoCal.size()<<std::endl;
  if(assoCal.size()!=3){
    throw cet::exception("[Numu0pi0p]") << "Where are the three planes for the calorimetry!" << std::endl;
  }
  int id_pl_collect = -99; // Plane ID for collection plane, which is supposed to be 2
  for(int id_pl = 0; id_pl < 3; id_pl++){
    if(assoCal[id_pl]->PlaneID().Plane == 2){
      id_pl_collect = id_pl;
    }
  }
  std::cout<<"collection plane id: "<<id_pl_collect<<std::endl;

////////////////////
//  // Only use collection plane for continuity check for now 
//  std::vector<double> vdEdx = assoCal.at(id_pl_collect)->dEdx();
//  std::vector<double> vresRange = assoCal.at(id_pl_collect)->ResidualRange();
//  double Trk_Length = assoCal.at(id_pl_collect)->Range();
// 
//  // dE/dx continutiy
//  //
//  int size_dEdx = vdEdx.size();
//  int unit = size_dEdx / 4;
//  // Average 4 hits and compare 4 groups of average dEdx
//  std::vector<double> v_dEdx_seg;
//    if(size_dEdx > 16){
//      for(int i_seg = 0; i_seg < 4; i_seg++){
//
//        v_dEdx_seg.push_back(0.25 * (vdEdx[size_dEdx - unit * i_seg - 1]
//                                   + vdEdx[size_dEdx - unit * i_seg - 2]
//                                   + vdEdx[size_dEdx - unit * i_seg - 3]
//                                   + vdEdx[size_dEdx - unit * i_seg - 4]));
//      }
//      //////
//      if(vdEdx_seg[i_seg]<10 && vdEdx)
//      auto sorted_v_dEdx_seg = v_dEdx_seg;
//      sort(sorted_v_dEdx_seg.begin(), sorted_v_dEdx_seg.end());
//    }
}

#endif

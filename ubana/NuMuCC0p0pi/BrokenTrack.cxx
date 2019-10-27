//Class created by Yifan Chen (yifan.chen@lhep.unibe.ch)
//October 4th, 2019

// _________________________________________________________________________________________________________________________________________________________________________________________________

#ifndef BROKENTRACK_CXX
#define BROKENTRACK_CXX

#include "BrokenTrack.h"
#include "FiducialVolume.h"
#include <algorithm>
// _________________________________________________________________________________________________________________________________________________________________________________________________

//_fiducial_volume.Configure(pset.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
//                             geo->DetHalfHeight(),
//                             2.*geo->DetHalfWidth(),
//                             geo->DetLength());

void BrokenTrack::MatchTracks(art::Ptr<recob::Track>& ThisTrack, std::vector< art::Ptr<recob::Track>>& TrackCollection){
  
  // Add spatial correction to the track start and end
  TVector3 Trk_start = ThisTrack->Vertex<TVector3>();
  auto Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_start.X(), Trk_start.Y(), Trk_start.Z()));
  TVector3 Trk_start_SCEcorr;
  Trk_start_SCEcorr.SetX(Trk_start.X() - Trk_start_offset.X());
  Trk_start_SCEcorr.SetY(Trk_start.Y() + Trk_start_offset.Y());
  Trk_start_SCEcorr.SetZ(Trk_start.Z() + Trk_start_offset.Z());

  TVector3 Trk_end = ThisTrack->End<TVector3>();
  auto Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Trk_end.X(), Trk_end.Y(), Trk_end.Z()));
  TVector3 Trk_end_SCEcorr;
  Trk_end_SCEcorr.SetX(Trk_end.X() - Trk_end_offset.X());
  Trk_end_SCEcorr.SetY(Trk_end.Y() + Trk_end_offset.Y());
  Trk_end_SCEcorr.SetZ(Trk_end.Z() + Trk_end_offset.Z());  

  TVector3 Trk_vec = Trk_end_SCEcorr - Trk_start_SCEcorr;
 
  //Initialize
  newTrk = false;
  Nr_mergedTrk = 1;
  trk_temp_end1 = Trk_start_SCEcorr;
  trk_temp_end2 = Trk_end_SCEcorr;
  trk_temp_length = (Trk_start_SCEcorr - Trk_end_SCEcorr).Mag();
 
  //std::cout<<"Trk_start X: "<< trk_temp_end1.X()<<", Y: "<< trk_temp_end1.Y()<<", Z: "<<trk_temp_end1.Z()<<std::endl;
  //std::cout<<"Trk_end X: "<< trk_temp_end2.X()<<", Y: "<< trk_temp_end2.Y()<<", Z: "<<trk_temp_end2.Z()<<std::endl;  

  // For the moment, this part only supports two track merging
  for(int i_trk = 0; i_trk < (int) TrackCollection.size(); i_trk++){
    if (TrackCollection[i_trk] != ThisTrack){
      // Add spatial correction to the matching track start and end
      TVector3 Rolling_Trk_start = TrackCollection[i_trk]->Vertex<TVector3>();
      auto Rolling_Trk_start_offset = SCE->GetCalPosOffsets(geo::Point_t(Rolling_Trk_start.X(), Rolling_Trk_start.Y(), Rolling_Trk_start.Z()));
      TVector3 Rolling_Trk_start_SCEcorr;
      Rolling_Trk_start_SCEcorr.SetX(Rolling_Trk_start.X() - Rolling_Trk_start_offset.X());
      Rolling_Trk_start_SCEcorr.SetY(Rolling_Trk_start.Y() + Rolling_Trk_start_offset.Y());
      Rolling_Trk_start_SCEcorr.SetZ(Rolling_Trk_start.Z() + Rolling_Trk_start_offset.Z());

      TVector3 Rolling_Trk_end = TrackCollection[i_trk]->End<TVector3>();
      auto Rolling_Trk_end_offset = SCE->GetCalPosOffsets(geo::Point_t(Rolling_Trk_end.X(), Rolling_Trk_end.Y(), Rolling_Trk_end.Z()));
      TVector3 Rolling_Trk_end_SCEcorr;
      Rolling_Trk_end_SCEcorr.SetX(Rolling_Trk_end.X() - Rolling_Trk_end_offset.X());
      Rolling_Trk_end_SCEcorr.SetY(Rolling_Trk_end.Y() + Rolling_Trk_end_offset.Y());
      Rolling_Trk_end_SCEcorr.SetZ(Rolling_Trk_end.Z() + Rolling_Trk_end_offset.Z());

      TVector3 Rolling_Trk_vec = Rolling_Trk_end_SCEcorr - Rolling_Trk_start_SCEcorr;

      //------ The angle between two vector
      double cosA = (Trk_vec * Rolling_Trk_vec) / (Trk_vec.Mag() * Rolling_Trk_vec.Mag());

      //------ The distance of Rolling track start/end to the track vector
      // In case one end of the Rolling track is out of TPC, we don't correct spatial coord of Rolling track. (Spatial correction is 0 when the given point is out of the TPC). On the other hand we use spatial correction for the selected track
      TVector3 Rs_s = Rolling_Trk_start_SCEcorr - Trk_start_SCEcorr;
      double angle_Rs = acos(abs((Rs_s * Trk_vec)/(Rs_s.Mag() * Trk_vec.Mag()))); // Limit the angle to be 0 - pi/2
      double dis_Rs = Rs_s.Mag() * sin(angle_Rs);

      TVector3 Re_s = Rolling_Trk_end_SCEcorr - Trk_start_SCEcorr;
      double angle_Re = acos(abs((Re_s * Trk_vec)/(Re_s.Mag() * Trk_vec.Mag())));// Limit the angle to be 0 - pi/2
      double dis_Re = Re_s.Mag() * sin(angle_Re);

      //------ Check if two tracks overlap
      bool overlap = false;
      // Rolling track in the track
      if ((Rolling_Trk_start_SCEcorr.X() - Trk_start_SCEcorr.X()) * (Rolling_Trk_start_SCEcorr.X() - Trk_end_SCEcorr.X()) < 0 || (Rolling_Trk_end_SCEcorr.X() - Trk_start_SCEcorr.X()) * (Rolling_Trk_end_SCEcorr.X() - Trk_end_SCEcorr.X()) < 0){ overlap = true;}
      if ((Rolling_Trk_start_SCEcorr.Y() - Trk_start_SCEcorr.Y()) * (Rolling_Trk_start_SCEcorr.Y() - Trk_end_SCEcorr.Y()) < 0 || (Rolling_Trk_end_SCEcorr.Y() - Trk_start_SCEcorr.Y()) * (Rolling_Trk_end_SCEcorr.Y() - Trk_end_SCEcorr.Y()) < 0){ overlap = true;}
      if ((Rolling_Trk_start_SCEcorr.Z() - Trk_start_SCEcorr.Z()) * (Rolling_Trk_start_SCEcorr.Z() - Trk_end_SCEcorr.Z()) < 0 || (Rolling_Trk_end_SCEcorr.Z() - Trk_start_SCEcorr.Z()) * (Rolling_Trk_end_SCEcorr.Z() - Trk_end_SCEcorr.Z()) < 0){ overlap = true;}
      // the track in rolling track
      if((Trk_start_SCEcorr.X() - Rolling_Trk_start_SCEcorr.X()) * (Trk_start_SCEcorr.X() - Rolling_Trk_start_SCEcorr.X()) < 0 || (Trk_end_SCEcorr.X() - Rolling_Trk_start_SCEcorr.X()) * (Trk_end_SCEcorr.X() - Rolling_Trk_end_SCEcorr.X()) < 0){ overlap = true;}
      if((Trk_start_SCEcorr.Y() - Rolling_Trk_start_SCEcorr.Y()) * (Trk_start_SCEcorr.Y() - Rolling_Trk_start_SCEcorr.Y()) < 0 || (Trk_end_SCEcorr.Y() - Rolling_Trk_start_SCEcorr.Y()) * (Trk_end_SCEcorr.Y() - Rolling_Trk_end_SCEcorr.Y()) < 0){ overlap = true;}
      if((Trk_start_SCEcorr.Z() - Rolling_Trk_start_SCEcorr.Z()) * (Trk_start_SCEcorr.Z() - Rolling_Trk_start_SCEcorr.Z()) < 0 || (Trk_end_SCEcorr.Z() - Rolling_Trk_start_SCEcorr.Z()) * (Trk_end_SCEcorr.Z() - Rolling_Trk_end_SCEcorr.Z()) < 0){ overlap = true;}

      // 1. the distance of rolling track ends to the track vector is under 30cm
      // 2. the angle in between the vec of rolling track and the track is smaller than 30 degree
      /// If the two vector is within 30 degree (0.52359878) difference, cos(0.52359878) = 0.86602540
      // 3. the rolling track and the selected track should not be overlapped
      if(dis_Rs < 30 && dis_Re < 30 && abs(cosA) >= 0.86602540 && overlap == false){

        newTrk = true;
        Nr_mergedTrk++;

        std::cout<<"Number of merged track: "<< Nr_mergedTrk<<std::endl;
        std::vector<double> dis;

//        dis.push_back((Trk_start_SCEcorr - Rolling_Trk_start_SCEcorr).Mag()); // S_RS
//        dis.push_back((Trk_start_SCEcorr - Rolling_Trk_end_SCEcorr).Mag()); //S_RE
//        dis.push_back((Trk_end_SCEcorr - Rolling_Trk_start_SCEcorr).Mag()); //E_RS
//        dis.push_back((Trk_end_SCEcorr - Rolling_Trk_end_SCEcorr).Mag()); //E_RE
//
//        std::cout<<"start_Rstart: "<< dis[0]<<"; start_Rend: "<< dis[1]<<"; end_Rstart: "<< dis[2]<<"; end_Rend: "<< dis[3]<<std::endl;
//
//        std::cout<<"The max distance: "<< *std::max_element(dis.begin(), dis.end())<<std::endl;
//        int Max = std::distance(dis.begin(), std::max_element(dis.begin(), dis.end()));
//        std::cout<<"Max: "<<Max<<std::endl;
//        //int Max = std::max_element(dis.begin(), dis.end() - dis.begin());
//
//        if(Max == 0){
//          trk_temp_end1 = Trk_start_SCEcorr;
//          trk_temp_end2 = Rolling_Trk_start_SCEcorr;
//
//          trk_end1 = Trk_start_SCEcorr;
//          trk_end2 = Rolling_Trk_start_SCEcorr;
//        }
//        if(Max == 1){
//          trk_end1 = Trk_start_SCEcorr;
//          trk_end2 = Rolling_Trk_end_SCEcorr;
//        }
//        if(Max == 2){
//          trk_end1 = Trk_end_SCEcorr;
//          trk_end2 = Rolling_Trk_start_SCEcorr;
//        }
//        if(Max == 3){
//          trk_end1 = Trk_end_SCEcorr;
//          trk_end2 = Rolling_Trk_end_SCEcorr;
//        }



        //---------
        dis.push_back((trk_temp_end1 - Rolling_Trk_start_SCEcorr).Mag()); // S_RS
        dis.push_back((trk_temp_end1 - Rolling_Trk_end_SCEcorr).Mag()); //S_RE
        dis.push_back((trk_temp_end2 - Rolling_Trk_start_SCEcorr).Mag()); //E_RS
        dis.push_back((trk_temp_end2 - Rolling_Trk_end_SCEcorr).Mag()); //E_RE

        //std::cout<<"start_Rstart: "<< dis[0]<<"; start_Rend: "<< dis[1]<<"; end_Rstart: "<< dis[2]<<"; end_Rend: "<< dis[3]<<std::endl;

        trk_temp_length = *std::max_element(dis.begin(), dis.end());
        //std::cout<<"The max distance: "<< *std::max_element(dis.begin(), dis.end())<<std::endl;

        if (*std::max_element(dis.begin(), dis.end()) > trk_temp_length){

          trk_temp_length = *std::max_element(dis.begin(), dis.end());

          int Max = std::distance(dis.begin(), std::max_element(dis.begin(), dis.end()));
          //std::cout<<"Max: "<<Max<<std::endl;

          if(Max == 0){
            trk_temp_end1 = trk_temp_end1;
            trk_temp_end2 = Rolling_Trk_start_SCEcorr;
          }
          if(Max == 1){
            trk_temp_end1 = trk_temp_end1;
            trk_temp_end2 = Rolling_Trk_end_SCEcorr;
          }
          if(Max == 2){
            trk_temp_end1 = trk_temp_end2;
            trk_temp_end2 = Rolling_Trk_start_SCEcorr;
          }
          if(Max == 3){
            trk_temp_end1 = trk_temp_end2;
            trk_temp_end2 = Rolling_Trk_end_SCEcorr;
          }
          //std::cout<<"temp end1 X: "<< trk_temp_end1.X()<<", Y: "<< trk_temp_end1.Y()<<", Z: "<<trk_temp_end1.Z()<<std::endl;
          //std::cout<<"temp end2 X: "<< trk_temp_end2.X()<<", Y: "<< trk_temp_end2.Y()<<", Z: "<<trk_temp_end2.Z()<<std::endl;
        }

      } // matched tracks (assume it's only two...)
    } // not itself
  } // loop all the tracks in this event

  trk_end1 = trk_temp_end1;
  trk_end2 = trk_temp_end2;
  trk_length = trk_temp_length;

  //std::cout<<"trk end1 X: "<< trk_temp_end1.X()<<", Y: "<< trk_temp_end1.Y()<<", Z: "<<trk_temp_end1.Z()<<std::endl;
  //std::cout<<"trk end2 X: "<< trk_temp_end2.X()<<", Y: "<< trk_temp_end2.Y()<<", Z: "<<trk_temp_end2.Z()<<std::endl;

} // function

bool BrokenTrack::NewTrk(){
  return newTrk;
}

int BrokenTrack::NumberMergedTracks(){
  return Nr_mergedTrk;
}

double BrokenTrack::TrkLen(){
  return trk_length;
}

TVector3 BrokenTrack::TrkEnd1(){
  return trk_end1;
}

TVector3 BrokenTrack::TrkEnd2(){
  return trk_end2;
}

//bool BrokenTrack::NewTrkFV(){
//  return newTrk_FV;
//}
//
//bool BrokenTrack::NewTrkContained(){
//  return newTrk_contained;
//}
#endif

#ifndef FINDDEADREGIONS_CXX
#define FINDDEADREGIONS_CXX

// Art includes
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "TH2F.h"

#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "FindDeadRegions.h"

FindDeadRegions::FindDeadRegions()//fhicl::ParameterSet const & p, art::ActivityRegistry & areg)
{

  // Load wire geometry at initialization
  LoadWireGeometry();

  CSchannelVec.resize(8256);
  CSstatusVec.resize(8256);

}

void FindDeadRegions::Configure(fhicl::ParameterSet const& pset) {
  _use_file   = pset.get< bool   > ( "UseFile",   false );
  _tolerance  = pset.get< double > ( "Tolerance", 0.6   ); //cm
  _ch_thres   = pset.get< int >    ( "ChThres",   4     ); 
}


void FindDeadRegions::SetChannelStatus(unsigned int ch, int status) {

  if (status < _ch_thres) {
    CSstatusVec.at(ch) = false;
  } else {
    CSstatusVec.at(ch) = true;
  }

}

void FindDeadRegions::LoadWireGeometry() {

  std::cout << "[FindDeadRegions] Loading wires from " << (_use_file ? "files." : "database.") << std::endl;

  ::art::ServiceHandle<geo::Geometry> geo;

  std::cout << "[FindDeadRegions] Loading geometry." << std::endl;

  if (!_use_file) {

    // **********
    // From Geometry
    // **********

    // Loop over all the channels
    for (unsigned int channel = 0; channel < 8256; ++channel) {

      std::vector<geo::WireID> const wire_v = geo->ChannelToWire(channel);
      geo::WireGeo const& wire_g = geo->Wire(wire_v[0]);

      auto const start = wire_g.GetStart();
      auto const end = wire_g.GetEnd();

      channelVec.push_back(channel);
      planeVec.push_back(wire_v[0].Plane);
      wireVec.push_back(wire_v[0].Wire);
      sxVec.push_back(start.X());
      syVec.push_back(start.Y());
      szVec.push_back(start.Z());
      exVec.push_back(end.X());
      eyVec.push_back(end.Y());
      ezVec.push_back(end.Z());

    } // channel loop
  } else {

    // **********
    // From File
    // **********

    std::ifstream geofile;
    geofile.open("ChannelWireGeometry_v2.txt");
     if (!geofile.is_open())
      std::cerr << "Problem opening file ChannelWireGeometry_v2.txt."  << std::endl;

    std::string string_channel;
    std::string string_plane;
    std::string string_wire;
    std::string string_sx;
    std::string string_sy;
    std::string string_sz;
    std::string string_ex;
    std::string string_ey;
    std::string string_ez;

    std::string dummy;

    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);
    getline(geofile,dummy);

    while(geofile >> string_channel)
    {
      geofile >> string_plane;
      geofile >> string_wire;
      geofile >> string_sx;
      geofile >> string_sy;
      geofile >> string_sz;
      geofile >> string_ex;
      geofile >> string_ey;
      geofile >> string_ez;

      channelVec.push_back(atoi(string_channel.c_str()));
      planeVec.push_back(atoi(string_plane.c_str()));
      wireVec.push_back(atoi(string_wire.c_str()));
      sxVec.push_back(atof(string_sx.c_str()));
      syVec.push_back(atof(string_sy.c_str()));
      szVec.push_back(atof(string_sz.c_str()));
      exVec.push_back(atof(string_ex.c_str()));
      eyVec.push_back(atof(string_ey.c_str()));
      ezVec.push_back(atof(string_ez.c_str()));
    }

    geofile.close();
  }

}




void FindDeadRegions::LoadChannelStatus() {

  std::cout << "[FindDeadRegions] Loading channel statuses." << std::endl;

  if (!_use_file) {

    // **********
    // From Database
    // **********

    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

    CSchannelVec.resize(8256);
    CSstatusVec.resize(8256);

    for (unsigned int channel = 0; channel < 8256; channel++) {

      CSchannelVec.push_back(channel);

      // Channel statuses: 1=dead, 3=noisy, 4=good
      if (chanFilt.Status(channel) < _ch_thres) {
        CSstatusVec.at(channel) = false;
      } else {
        CSstatusVec.at(channel) = true;
      }
    }
  } else {

    // **********
    // From File
    // **********

    std::cout << "[FindDeadRegions] Reading channel status form file." << std::endl;
    std::ifstream chanstatfile;
    chanstatfile.open("ChanStatus.txt");
    if (!chanstatfile.is_open())
      std::cerr << "Problem opening file ChanStatus.txt." << std::endl;

    std::string string_CSchannel;
    std::string string_CSstatus;

    unsigned int CSchannel;
    bool CSstatus;

    CSstatusVec.clear();
    CSstatusVec.clear();

    while(chanstatfile >> string_CSchannel)
    {
      chanstatfile >> string_CSstatus;

      CSchannel = atoi(string_CSchannel.c_str());
      CSstatus = atoi(string_CSstatus.c_str());

      CSchannelVec.push_back(CSchannel);
      CSstatusVec.push_back(CSstatus);
    }

    chanstatfile.close();
  }

  std::cout << "[FindDeadRegions] Loading ended." << std::endl;

}


void FindDeadRegions::CreateBWires() {

  std::cout << "[FindDeadRegions] Now creating BWires" << std::endl;

  BWires_U.clear();
  BWires_V.clear();
  BWires_Y.clear();

  bool isGoodChannel;

  isGoodChannel = true;
  for(int i = 0; i < 2400; i++) {
    if((CSstatusVec.at(i) == false) && (isGoodChannel == true)) {
      isGoodChannel = false;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = true;

      BWires_U.push_back(BWire);
    }
    else if((CSstatusVec.at(i) == true) && (isGoodChannel == false)) {
      isGoodChannel = true;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i-1);
      BWire.y_start = syVec.at(i-1);
      BWire.z_start = szVec.at(i-1);
      BWire.y_end = eyVec.at(i-1);
      BWire.z_end = ezVec.at(i-1);
      BWire.isLowWire = false;

      BWires_U.push_back(BWire);
    }
    else if((i == 2399) && (CSstatusVec.at(i) == false) && (isGoodChannel == false)) {
      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = false;

      BWires_U.push_back(BWire);

    }
  }

  isGoodChannel = true;
  for(int i = 2400; i < 4800; i++) {
    if((CSstatusVec.at(i) == false) && (isGoodChannel == true)) {
      isGoodChannel = false;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = true;

      BWires_V.push_back(BWire);
    }
    else if((CSstatusVec.at(i) == true) && (isGoodChannel == false)) {
      isGoodChannel = true;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i-1);
      BWire.y_start = syVec.at(i-1);
      BWire.z_start = szVec.at(i-1);
      BWire.y_end = eyVec.at(i-1);
      BWire.z_end = ezVec.at(i-1);
      BWire.isLowWire = false;

      BWires_V.push_back(BWire);
    }
    else if((i == 4799) && (CSstatusVec.at(i) == false) && (isGoodChannel == false)) {
      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = false;

      BWires_V.push_back(BWire);
    }
  }

  isGoodChannel = true;
  for(int i = 4800; i < 8256; i++) {

    if((CSstatusVec.at(i) == false) && (isGoodChannel == true)) {
      isGoodChannel = false;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = true;

      BWires_Y.push_back(BWire);
    }
    else if((CSstatusVec.at(i) == true) && (isGoodChannel == false)) {
      isGoodChannel = true;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i-1);
      BWire.y_start = syVec.at(i-1);
      BWire.z_start = szVec.at(i-1);
      BWire.y_end = eyVec.at(i-1);
      BWire.z_end = ezVec.at(i-1);
      BWire.isLowWire = false;

      BWires_Y.push_back(BWire);
    }
    else if((i == 8255) && (CSstatusVec.at(i) == false) && (isGoodChannel == false)) {
      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = false;

      BWires_Y.push_back(BWire);
    }
  }

  std::cout << "[FindDeadRegions] Finishes creating BWires." << std::endl;

  std::cout << "[FindDeadRegions] Number of BWires in U " << BWires_U.size() << std::endl;
  std::cout << "[FindDeadRegions] Number of BWires in V " << BWires_V.size() << std::endl;
  std::cout << "[FindDeadRegions] Number of BWires in Y " << BWires_Y.size() << std::endl;

  return;
}


//___________________________________________________________________________________________________
bool FindDeadRegions::NearDeadReg2P(float yVal, float zVal, float tolerance) {

  float minDist_U = 100000.0;
  float minDist_V = 100000.0;
  float minDist_Y = 100000.0;

  for (unsigned int i = 0; i < BWires_U.size(); i++) {
    float m = (BWires_U[i].y_end-BWires_U[i].y_start)/(BWires_U[i].z_end-BWires_U[i].z_start);
    float b = BWires_U[i].y_start - m*BWires_U[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist_U) {
      minDist_U = dist;
    }

    if (BWires_U[i].isLowWire == true) {
      float m_next = (BWires_U[i+1].y_end-BWires_U[i+1].y_start)/(BWires_U[i+1].z_end-BWires_U[i+1].z_start);
      float b_next = BWires_U[i+1].y_start - m*BWires_U[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal <= y1) && (yVal >= y2)) {
        minDist_U = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_V.size(); i++) {
    float m = (BWires_V[i].y_end-BWires_V[i].y_start)/(BWires_V[i].z_end-BWires_V[i].z_start);
    float b = BWires_V[i].y_start - m*BWires_V[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist_V) {
      minDist_V = dist;
    }
 
    if (BWires_V[i].isLowWire == true) {
      float m_next = (BWires_V[i+1].y_end-BWires_V[i+1].y_start)/(BWires_V[i+1].z_end-BWires_V[i+1].z_start);
      float b_next = BWires_V[i+1].y_start - m*BWires_V[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal >= y1) && (yVal <= y2)) {
        minDist_V = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_Y.size(); i++) {
    float dist = fabs(zVal-BWires_Y[i].z_start);

    if (dist < minDist_Y) {
      minDist_Y = dist;
    }

    if (BWires_Y[i].isLowWire == true) {
      float z1 = BWires_Y[i].z_start;
      float z2 = BWires_Y[i+1].z_start;

      if ((zVal >= z1) && (zVal <= z2)) {
        minDist_Y = 0.0;
      }
    }
  }

  if ((minDist_U < tolerance) && (minDist_V < tolerance)) {
    return true;
  }
  else if ((minDist_U < tolerance) && (minDist_Y < tolerance)) {
    return true;
  }
  else if ((minDist_V < tolerance) && (minDist_Y < tolerance)) {
    return true;
  }
  else {
    return false;
  }
}


//____________________________________________________________________________________________________
bool FindDeadRegions::NearDeadRegCollection(float zVal, float tolerance) {

  float minDist_Y = 100000.0;

  for (unsigned int i = 0; i < BWires_Y.size(); i++) {
    float dist = fabs(zVal-BWires_Y[i].z_start);

    if (dist < minDist_Y) {
      minDist_Y = dist;
    }

    if (BWires_Y[i].isLowWire == true) {
      float z1 = BWires_Y[i].z_start;
      float z2 = BWires_Y[i+1].z_start;

      if ((zVal >= z1) && (zVal <= z2)) {
        minDist_Y = 0.0;
      }
    }
  }

  if (minDist_Y < tolerance) {
    return true;
  }
  
  return false;

}



//____________________________________________________________________________________________________
bool FindDeadRegions::NearDeadReg3P(float yVal, float zVal, float tolerance) {

  float minDist = 100000.0;

  for (unsigned int i = 0; i < BWires_U.size(); i++) {
    float m = (BWires_U[i].y_end-BWires_U[i].y_start)/(BWires_U[i].z_end-BWires_U[i].z_start);
    float b = BWires_U[i].y_start - m*BWires_U[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist) {
      minDist = dist;
    }

    if (BWires_U[i].isLowWire == true) {
      float m_next = (BWires_U[i+1].y_end-BWires_U[i+1].y_start)/(BWires_U[i+1].z_end-BWires_U[i+1].z_start);
      float b_next = BWires_U[i+1].y_start - m*BWires_U[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal <= y1) && (yVal >= y2)) {
        minDist = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_V.size(); i++) {
    float m = (BWires_V[i].y_end-BWires_V[i].y_start)/(BWires_V[i].z_end-BWires_V[i].z_start);
    float b = BWires_V[i].y_start - m*BWires_V[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist) {
      minDist = dist;
    }

    if (BWires_V[i].isLowWire == true) {
      float m_next = (BWires_V[i+1].y_end-BWires_V[i+1].y_start)/(BWires_V[i+1].z_end-BWires_V[i+1].z_start);
      float b_next = BWires_V[i+1].y_start - m*BWires_V[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal >= y1) && (yVal <= y2)) {
        minDist = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_Y.size(); i++) {
    float dist = fabs(zVal-BWires_Y[i].z_start);

    if (dist < minDist) {
      minDist = dist;
    }

    if (BWires_Y[i].isLowWire == true) {
      float z1 = BWires_Y[i].z_start;
      float z2 = BWires_Y[i+1].z_start;

      if ((zVal >= z1) && (zVal <= z2)) {
        minDist = 0.0;
      }
    }
  }

  if (minDist < tolerance) {
    return true;
  }
  else {
    return false;
  }
}



//__________________________________________________________________________________________________
TH2F* FindDeadRegions::GetDeadRegionHisto2P(){

  TH2F *deadReg2P = new TH2F("deadReg2P","",10350,0.0,1035.0,2300,-115.0,115.0);

  for (int i = 1; i <= deadReg2P->GetNbinsX(); i++) {
    for (int j = 1; j <= deadReg2P->GetNbinsY(); j++) {
      if (NearDeadReg2P(deadReg2P->GetYaxis()->GetBinCenter(j),deadReg2P->GetXaxis()->GetBinCenter(i),_tolerance) == true) {
        deadReg2P->SetBinContent(i,j,1.0);
      }
    }
  }

  return deadReg2P;
}

//__________________________________________________________________________________________________
void FindDeadRegions::GetDeadRegionHisto2P(TH2F* deadReg2P){

  for (int i = 1; i <= deadReg2P->GetNbinsX(); i++) {
    for (int j = 1; j <= deadReg2P->GetNbinsY(); j++) {
      if (NearDeadReg2P(deadReg2P->GetYaxis()->GetBinCenter(j),deadReg2P->GetXaxis()->GetBinCenter(i),_tolerance) == true) {
        deadReg2P->SetBinContent(i,j,1.0);
      }
    }
  }

  return;
}


//__________________________________________________________________________________________________
TH2F* FindDeadRegions::GetDeadRegionHisto3P(){

  TH2F *deadReg3P = new TH2F("deadReg3P","",10350,0.0,1035.0,2300,-115.0,115.0);

  for (int i = 1; i <= deadReg3P->GetNbinsX(); i++) {
    for (int j = 1; j <= deadReg3P->GetNbinsY(); j++) {
      if (NearDeadReg3P(deadReg3P->GetYaxis()->GetBinCenter(j),deadReg3P->GetXaxis()->GetBinCenter(i),_tolerance) == true) {
        deadReg3P->SetBinContent(i,j,1.0);
      }
    }
  }

  return deadReg3P;
}

//__________________________________________________________________________________________________
void FindDeadRegions::GetDeadRegionHisto3P(TH2F* deadReg3P){

  for (int i = 1; i <= deadReg3P->GetNbinsX(); i++) {
    for (int j = 1; j <= deadReg3P->GetNbinsY(); j++) {
      if (NearDeadReg3P(deadReg3P->GetYaxis()->GetBinCenter(j),deadReg3P->GetXaxis()->GetBinCenter(i),_tolerance) == true) {
        deadReg3P->SetBinContent(i,j,1.0);
      }
    }
  }

  return;
}

#endif

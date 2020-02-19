#ifndef FIDUCIALVOLUME_CXX
#define FIDUCIALVOLUME_CXX

#include "FiducialVolume.h"
#include <iostream>

namespace ubana {

  FiducialVolume::FiducialVolume()
  {
  }

  void FiducialVolume::Configure(fhicl::ParameterSet const& pset, double det_half_height, double det_width, double det_length)
  {
    _border_x_low    = pset.get<std::vector<double>>("BorderXLow");
    _border_x_high   = pset.get<std::vector<double>>("BorderXHigh");
    _border_y_low    = pset.get<std::vector<double>>("BorderYLow");
    _border_y_high   = pset.get<std::vector<double>>("BorderYHigh");
    _border_z_low    = pset.get<std::vector<double>>("BorderZLow");
    _border_z_high   = pset.get<std::vector<double>>("BorderZHigh");

    if ( (_border_x_low.size() != _border_x_high.size()) 
      || (_border_x_low.size() != _border_y_low.size())
      || (_border_x_low.size() != _border_y_high.size())
      || (_border_x_low.size() != _border_z_low.size())
      || (_border_x_low.size() != _border_z_high.size()) ) {

      std::cout << "[FiducialVolume] Input vectors size mismatch." << std::endl;
      throw std::exception();
    }

    // In this specific application, 
    // always assume the fix box is the containment requirement
    // and the following boxes are the forbidden regions for the vertex

    _n_fv = _border_x_low.size();
    if (_n_fv < 1){
      throw cet::exception("[Fiducial volume]")<< "The FV is not set correctly" << std::endl;
    }

    _det_half_height = det_half_height;
    _det_width       = det_width;
    _det_length      = det_length;

    _configured = true;
  }

  void FiducialVolume::PrintConfig() {

    std::cout << "--- FiducialVolume configuration:" << std::endl;
    std::cout << "---   Number of fiducial volumes = " << _n_fv << std::endl;
    std::cout << "---   _det_half_height = " << _det_half_height << std::endl;
    std::cout << "---   _det_width       = " << _det_width << std::endl;
    std::cout << "---   _det_length      = " << _det_length << std::endl;
    for (size_t i = 0; i < _n_fv; i++) {
      std::cout << "---   Fiducial volume number " << i << std::endl;
      std::cout << "---     _border_x_low    = " << _border_x_low.at(i) << std::endl;
      std::cout << "---     _border_x_high   = " << _border_x_high.at(i) << std::endl;
      std::cout << "---     _border_y_low    = " << _border_y_low.at(i) << std::endl;
      std::cout << "---     _border_y_high   = " << _border_y_high.at(i) << std::endl;
      std::cout << "---     _border_z_low    = " << _border_z_low.at(i) << std::endl;
      std::cout << "---     _border_z_high   = " << _border_z_high.at(i) << std::endl;
    }

  }

  //////////////////////////////
  //The 0th box is the containment requirement
  //The rest boxes are the "dead wire" regions which the vertices should be excluded in these regions
  ////////////////////////////

  bool FiducialVolume::VertexInActive(TVector3 x) {

    return this->VertexInActive(x.X(), x.Y(), x.Z());

  }

  bool FiducialVolume::VertexInFV(double* x) {

    return this->VertexInFV(x[0], x[1], x[2]);

  }

  bool FiducialVolume::VertexInFV(TVector3 x) {

    return this->VertexInFV(x.X(), x.Y(), x.Z());

  }

  bool FiducialVolume::PointContain(TVector3 x) {

    return this->PointContain(x.X(), x.Y(), x.Z());

  }

  // Assume x1 is track start and x2 is track end
  bool FiducialVolume::TrackContain(TVector3 x1, TVector3 x2) {

    return (this->VertexInFV(x1.X(), x1.Y(), x1.Z()) && this->PointContain(x2.X(), x2.Y(), x2.Z()));

  }

  /////////////////// For Dirt (if vertex in the active area or not)
  bool FiducialVolume::VertexInActive(double x, double y, double z) {

    // Construct a vector
    ::geoalgo::Vector the_point(x, y, z);

    // Vertex in the active volume
    // Construct the fiducial volume
    ::geoalgo::AABox fidvol(0., 
                            -1.*_det_half_height,
                            0.,
                            _det_width,
                            _det_half_height,
                            _det_length);

    // Check if the vector is in the FV
    if(fidvol.Contain(the_point))
      return true;
    else
      return false;
  }

  /////////////////// For Containment
  bool FiducialVolume::PointContain(double x, double y, double z) {

    // Construct a vector
    ::geoalgo::Vector the_point(x, y, z);

    // Containment only requires the point to be in the 0th box
    // Construct the fiducial volume
    ::geoalgo::AABox fidvol(_border_x_low.at(0), 
                            -1.*_det_half_height + _border_y_low.at(0), 
                            _border_z_low.at(0),
                            _det_width - _border_x_high.at(0), 
                            _det_half_height - _border_y_high.at(0), 
                            _det_length - _border_z_high.at(0));

    // Check if the vector is in the FV
    if(fidvol.Contain(the_point))
      return true;
    else
      return false;
  }

  ////////////////// For Vertex
  bool FiducialVolume::VertexInFV(double x, double y, double z) {

    // Construct a vector
    ::geoalgo::Vector the_point(x, y, z);

    if (_n_fv == 1){
      return this->PointContain(x, y, z);
    }

    if (_n_fv > 1){
      // Vertex shouldn't be in the other boxes than 0th
      for (size_t i = 1; i < _n_fv; i++) {
        ::geoalgo::AABox fidvol(_border_x_low.at(i),
                                -1.*_det_half_height + _border_y_low.at(i),
                                _border_z_low.at(i),
                                _det_width - _border_x_high.at(i),
                                _det_half_height - _border_y_high.at(i),
                                _det_length - _border_z_high.at(i));
        // Check if the vector in the sub-boxes
        if (fidvol.Contain(the_point)){
          return false;
        }  
      }
      return true;
    }
    
    // dumb return
    return true;
  }
}


#endif

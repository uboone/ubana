#ifndef FIDUCIALVOLUME_CXX
#define FIDUCIALVOLUME_CXX

#include "ubana/Utilities/FiducialVolume.h"
#include "fhiclcpp/ParameterSet.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoVector.h"

#include <iostream>

namespace ubana {

  FiducialVolume::FiducialVolume(fhicl::ParameterSet const& pset,
                                 double det_half_height,
                                 double det_width,
                                 double det_length)
  {
    _border_x_low = pset.get<std::vector<double>>("BorderXLow");
    _border_x_high = pset.get<std::vector<double>>("BorderXHigh");
    _border_y_low = pset.get<std::vector<double>>("BorderYLow");
    _border_y_high = pset.get<std::vector<double>>("BorderYHigh");
    _border_z_low = pset.get<std::vector<double>>("BorderZLow");
    _border_z_high = pset.get<std::vector<double>>("BorderZHigh");

    if ((_border_x_low.size() != _border_x_high.size()) ||
        (_border_x_low.size() != _border_y_low.size()) ||
        (_border_x_low.size() != _border_y_high.size()) ||
        (_border_x_low.size() != _border_z_low.size()) ||
        (_border_x_low.size() != _border_z_high.size())) {

      std::cout << "[FiducialVolume] Input vectors size mismatch." << std::endl;
      throw std::exception();
    }

    _n_fv = _border_x_low.size();

    _det_half_height = det_half_height;
    _det_width = det_width;
    _det_length = det_length;
  }

  void FiducialVolume::PrintConfig() const
  {

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

  bool FiducialVolume::InFV(double const* x) const { return InFV(x[0], x[1], x[2]); }
  bool FiducialVolume::InFV(TVector3 const& x) const { return InFV(x.X(), x.Y(), x.Z()); }
  bool FiducialVolume::InFV(geo::Point_t const& x) const { return InFV(x.X(), x.Y(), x.Z()); }

  bool FiducialVolume::InFV(TVector3 const& x1, TVector3 const& x2) const
  {
    return InFV(x1.X(), x1.Y(), x1.Z()) and InFV(x2.X(), x2.Y(), x2.Z());
  }

  bool FiducialVolume::InFV(double x, double y, double z) const
  {
    // Construct a vector
    ::geoalgo::Vector the_point(x, y, z);

    for (size_t i = 0; i < _n_fv; i++) {

      // Construct the fiducial volume
      ::geoalgo::AABox fidvol(_border_x_low.at(i),
                              -1. * _det_half_height + _border_y_low.at(i),
                              _border_z_low.at(i),
                              _det_width - _border_x_high.at(i),
                              _det_half_height - _border_y_high.at(i),
                              _det_length - _border_z_high.at(i));

      // Check if the vector is in the FV
      if (fidvol.Contain(the_point)) return true;
    }

    return false;
  }

}

#endif

#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include "fhiclcpp/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "TVector3.h"

#include <vector>

namespace ubana {

  class FiducialVolume {
  public:
    FiducialVolume(fhicl::ParameterSet const& p,
                   double det_half_height,
                   double det_width,
                   double det_length);

    /// Print the current configuration
    void PrintConfig() const;

    /// Returns true if the point is in the FV
    bool InFV(double x, double y, double z) const;
    bool InFV(double const* x) const;
    bool InFV(TVector3 const& x) const;
    bool InFV(geo::Point_t const& x) const;

    /// Returns true if BOTH points are in the FV
    bool InFV(TVector3 const& x1, TVector3 const& x2) const;

  protected:
    double _det_half_height;
    double _det_width;
    double _det_length;

    std::vector<double> _border_x_low;
    std::vector<double> _border_x_high;
    std::vector<double> _border_y_low;
    std::vector<double> _border_y_high;
    std::vector<double> _border_z_low;
    std::vector<double> _border_z_high;

    size_t _n_fv;
  };
}

#endif

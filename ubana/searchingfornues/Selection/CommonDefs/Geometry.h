#ifndef GEOMETRYFUNCS_H
#define GEOMETRYFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace searchingfornues
{
  float distance2d(const float& x1, const float& y1,
                  const float& x2, const float& y2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2));
  }

  float distance3d(const float& x1, const float& y1, const float& z1,
                  const float& x2, const float& y2, const float& z2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2) +
                (z1-z2)*(z1-z2));
  }

  float YZtoPlanecoordinate(const float y, const float z, const int plane)
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    double _wire2cm = geom->WirePitch(0, 0, 0);
    return geom->WireCoordinate(y, z, geo::PlaneID(0, 0, plane)) * _wire2cm;
  }
} // namespace searchingfornues

#endif

#include "spatial_discretization.h"

//######################################################################
std::string spatial_discretization_method_name(
    const SpatialDiscretizationMethod discretization_method)
{
  switch (spatial_discretization_method)
  {
    case SpatialDiscretizationMethod::FINITE_VOLUME:
      return "FINITE_VOLUME";
    case SpatialDiscretizationMethod::PIECEWISE_LINEAR_CONTINUOUS:
      return "PIECEWISE_LINEAR_CONTINUOUS";
    case SpatialDiscretizationMethod::PIECEWISE_LINEAR_DISCONTINIOUS:
      return "PIECEWISE_LINEAR_DISCONTINUOUS";
    case SpatialDiscretizationMethod::LAGRANGE_CONTINUOUS:
      return "LAGRANGE_CONTINUOUS";
    case SpatialDiscretizationMethod::LAGRANGE_DISCONTINUOUS:
      return "LAGRANGE_DISCONTINUOUS";
    default:
      return "UNDEFINED";
  }
}
#include "grid_structs.h"

std::string coordinate_system_name(const CoordinateSystem coordinate_system)
{
  switch (coordinate_system)
  {
    case CoordinateSystem::CARTESIAN:     return "CARTESIAN";
    case CoordinateSystem::CYLINDRICAL:   return "CYLINDRICAL";
    case CoordinateSystem::SPHERICAL:     return "SPHERICAL";
    default:                              return "UNDEFINED";
  }
}

#ifndef GRID_STRUCTS_H
#define GRID_STRUCTS_H

#include <string>
#include <vector>
#include <memory>


/// Coordinate system types.
enum class CoordinateSystem
{
  CARTESIAN   = 1,  ///< \f$ (x, y, z) \f$ coordinates
  CYLINDRICAL = 2,  ///< \f$ (r, z, \varphi) \f$ coordinates
  SPHERICAL   = 3   ///< \f$ (r, \varphi, \theta) \f$ coordinates
};

/// Get the coordinate system type as a string.
std::string coordinate_system_name(const CoordinateSystem coordinate_system);

//######################################################################

class Point;
typedef Point Vertex;
typedef Point Centroid;
typedef Point Normal;

class Face;
class Cell;
class Mesh;

//######################################################################

/**
 * \brief Create a 1D mesh from a list of vertices.
 * \param vertices A list of vertex locations.
 * \param coordinate_system The coordinate system type. The default is
 *                          Cartesian coordinates.
 */
std::shared_ptr<Mesh>
create_1d_mesh(const std::vector<double>& vertices,
               const CoordinateSystem coordinate_system =
                   CoordinateSystem::CARTESIAN);

/**
 * \brief Create a zoned 1D mesh.
 * \param zone_edges The edges of mesh zones. There should be one more
 *                   zone edge than number of zones.
 * \param zone_subdivisions The number of cells per zone.
 * \param material_ids The material ID per zone.
 * \param coordinate_system The coordinate system type. The default is
 *                          Cartesian coordinates.
 */
std::shared_ptr<Mesh>
create_1d_mesh(const std::vector<double>& zone_edges,
               const std::vector<size_t>& zone_subdivisions,
               const std::vector<int>& material_ids,
               const CoordinateSystem coordinate_system =
                   CoordinateSystem::CARTESIAN);


#include "point.h"

#endif //GRID_STRUCTS_H

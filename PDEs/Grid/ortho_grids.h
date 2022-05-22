#ifndef GRID_STRUCTS_H
#define GRID_STRUCTS_H

#include "mesh.h"
#include "point.h"

using namespace pdes::Grid;
using CoordSys = CoordinateSystem;


/**
 * Create a 1D mesh from a list of vertices.
 *
 * \param vertices A list of vertex locations.
 * \param coordinate_system The coordinate system type. The default is
 *                          Cartesian coordinates.
 */
Mesh
create_1d_mesh(const std::vector<double> vertices,
               const CoordSys coordinate_system = CoordSys::CARTESIAN,
               const bool verbose = false);


/**
 * Create a zoned 1D mesh.
 *
 * Zones are defined by edges, a number of subdivisions (cells), and a
 * material ID. This allows for non-uniform cells throughout the mesh and
 *
 * \param zone_edges The edges of mesh zones. There should be one more
 *                   zone edge than number of zones.
 * \param zone_subdivisions The number of cells per zone.
 * \param material_ids The material ID per zone.
 * \param coordinate_system The coordinate system type. The default is
 *                          Cartesian coordinates.
 */
Mesh
create_1d_mesh(const std::vector<double> zone_edges,
               const std::vector<size_t> zone_subdivisions,
               const std::vector<int> material_ids,
               const CoordSys coordinate_system = CoordSys::CARTESIAN,
               const bool verbose = false);

#endif //GRID_STRUCTS_H

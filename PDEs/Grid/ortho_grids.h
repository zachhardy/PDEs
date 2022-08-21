#ifndef GRID_STRUCTS_H
#define GRID_STRUCTS_H

#include "mesh.h"
#include "cartesian_vector.h"


namespace Grid
{

  /**
   * Create a 1D mesh from a list of vertices.
   *
   * \param vertices A list of vertex locations.
   * \param coordinate_system The coordinate system type. The default is
   *                          Cartesian coordinates.
   */
  std::shared_ptr<Mesh>
  create_1d_orthomesh(const std::vector<double>& vertices,
                      const CoordinateSystemType coordinate_system,
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
  std::shared_ptr<Mesh>
  create_1d_orthomesh(const std::vector<double>& zone_edges,
                      const std::vector<size_t>& zone_subdivisions,
                      const std::vector<int>& material_ids,
                      const CoordinateSystemType coordinate_system,
                      const bool verbose = false);


  /**
   * Create a 2D orthogonal mesh from a list of x and y vertices.
   *
   * This routine generates a mesh whose vertices are defined by the outer
   * product of the x and y vertices specified. For example, <tt> x_vertices =
   * [0.0, 1.0] </tt> and <tt> y_vertices = [0.0, 1.0] </tt>, this would form
   * a mesh with vertices <tt>(0.0, 0.0)</tt>, <tt>(1.0, 0.0)</tt>, <tt>(0.0,
   * 1.0)</tt>, and <tt>(1.0, 1.0)</tt>.
   */
  std::shared_ptr<Mesh>
  create_2d_orthomesh(const std::vector<double>& x_vertices,
                      const std::vector<double>& y_vertices,
                      const bool verbose = false);

}
#endif //GRID_STRUCTS_H

#ifndef MESH_H
#define MESH_H

#include "grid_structs.h"
#include "cell.h"

#include <vector>
#include <memory>

namespace grid
{

/**
 * \brief A class that represents a general computational mesh.
 *
 * A Mesh is characterized by a spatial dimension and coordinate system type
 * and is defined by a collection geometric objects up to its dimension. At the
 * lowest level, a Mesh is comprised of 0D vertices. Connections between
 * vertices define 1D edges. A closed collection of edges make up 2D surfaces.
 * Lastly, a bound collection of surfaces then make up 3D volumes. Of course,
 * a <tt>dim</tt>-dimensional mesh only contains objects up to dimension \p dim.
 */
class Mesh
{
public:
  const unsigned int dim;
  const CoordinateSystem coord_sys;

  std::vector<Vertex> vertices;
  std::vector<std::shared_ptr<Cell>> cells;
  std::vector<size_t> boundary_cell_ids;

public:
  /// Default constructor.
  Mesh(const unsigned int dimension, const CoordinateSystem coordinate_system)
      : coord_sys(coordinate_system), dim(dimension)
  {}

public:
  void compute_geometric_info();
  void establish_connectivity();
};

}

#endif //MESH_H

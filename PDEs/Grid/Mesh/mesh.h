#ifndef MESH_H
#define MESH_H

#include "../grid_structs.h"
#include "../Cell/cell.h"

#include <vector>
#include <memory>


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
public:   /*---------- Public Attributes ----------*/

  const unsigned int dim; ///< The spatial dimension of the mesh.
  const CoordinateSystem  coord_sys; ///< The coordinate system type.

  /// A collection of Point objects that define the Mesh.
  std::vector<Vertex> vertices;

  /// A collection of Cell objects that define the Mesh.
  std::vector<std::shared_ptr<Cell>> cells;

  /// IDs of the cells that lie on the boundary.
  std::vector<size_t> boundary_cell_ids;

public:   /*---------- Constructors, Destructors, and Assignments ----------*/

  /**
   * \brief Default constructor.
   * \param dimension The spatial dimension.
   * \param coordinate_system The coordinate system type.
   */
  Mesh(const unsigned int dimension, const CoordinateSystem coordinate_system)
    : coord_sys(coordinate_system), dim(dimension)
  {}

public:   /*---------- Routines ----------*/

  /// Compute the geometric information for the cells and faces.
  void compute_geometric_info();

  /// Establish the relationships between cells and their neighbors.
  void establish_connectivity();
};

#endif //MESH_H

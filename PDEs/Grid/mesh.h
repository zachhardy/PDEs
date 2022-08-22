#ifndef MESH_H
#define MESH_H

#include "cell.h"
#include "cartesian_vector.h"

#include <vector>
#include <memory>
#include <cstddef>


namespace Grid
{
  /**
   * Coordinate systems available for meshes.
   */
  enum class CoordinateSystemType
  {
    CARTESIAN = 0,  ///< \f$(x, y, z)\f$ coordinates.
    CYLINDRICAL = 1,  ///< \f$(r, z, \varphi)\f$ coordinates.
    SPHERICAL = 2   ///< \f$(r, \varphi, \theta)\f$ coordinates.
  };


  /**
   * Return the coordinate system type as a string.
   */
  std::string
  coordinate_system_str(const CoordinateSystemType coord_sys);


  /**
   * A class that represents a general computational mesh.
   *
   * A Mesh is characterized by a spatial dimension and coordinate system type
   * and is defined by a collection geometric objects up to its dimension. At
   * the lowest level, a Mesh comprises 0D vertices. Connections between
   * vertices define 1D edges. A closed collection of edges make up 2D surfaces.
   * Lastly, a bound collection of surfaces then make up 3D volumes. Of course,
   * a <tt>dim</tt>-dimensional mesh only contains objects up to dimension
   * \p dim.
   */
  class Mesh
  {
  public:
    /**
     * The spatial dimension of the mesh.
     */
    const unsigned int dim;

    /**
     * The coordinate system the mesh is defined in.
     */
    const CoordinateSystemType coord_sys;

    /**
     * A list of the vertices on the mesh.
     */
    std::vector<Vertex> vertices;

    /**
     * A list of the cells that make up the mesh. See \ref Cell.
     */
    std::vector<Cell> cells;

    /**
     * A list of the boundary cell IDs. This is used as a convenience for
     * access to boundary cells without having to iterate over interior cells.
     */
    std::vector<size_t> boundary_cell_ids;

    /**
     * A mapping from the cell ID to its row(i), column(j), level(k) index.
     * This is only used for orthogonal meshes.
     */
    std::vector<std::vector<size_t>> ijk_mapping;

    /**
     * Construct a mesh with the specified dimension and coordinate system.
     *
     * \param dimension The spatial dimension of the mesh.
     * \param coordinate_system The coordinate system of the mesh.
     */
    Mesh(const unsigned int dimension,
         const CoordinateSystemType coordinate_system);

    /**
     * Establish the relationships between cells and their neighbors.
     *
     * Connectivity is established by
     * 1. For each vertex, determine which cell objects the vertex belongs to.
     * 2. For each face of each cell, compare vertex ids to those of adjacent
     *    cells as defined by the step 1.
     * 3. If vertex IDs match, set the neighbor properties.
     *
     * \note This routine may be quite expensive. Only use this when
     *       connectivity cannot be established a-priori. Generally, this
     *       routine should only be utilized for unstructured meshes.
     */
    void
    compute_geometric_info();

    /** Compute the geometric properties of the cells and faces. */
    void
    establish_connectivity();
  };

}

#endif //MESH_H

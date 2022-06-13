#ifndef MESH_H
#define MESH_H

#include "cell.h"
#include "point.h"

#include <vector>
#include <memory>
#include <cstddef>


namespace Grid
{
  enum class CoordinateSystem
  {
    CARTESIAN = 0,  ///< \f$(x, y, z)\f$ coordinates.
    CYLINDRICAL = 1,  ///< \f$(r, z, \varphi)\f$ coordinates.
    SPHERICAL = 2   ///< \f$(r, \varphi, \theta)\f$ coordinates.
  };

  /** Return the coordinate system type as a string. */
  std::string coordinate_system_str(const CoordinateSystem coord_sys);

  //###########################################################################

  /**
   * A class that represents a general computational mesh.
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
    const size_t           dim;
    const CoordinateSystem coord_sys;

    std::vector<Vertex> vertices;
    std::vector<Cell>   cells;

    std::vector<size_t> boundary_cell_ids;

  public:
    /**
     * Construct a mesh with the specified dimension and coordinate system.
     *
     * \param dimension The spatial dimension of the mesh.
     * \param coordinate_system The coordinate system of the mesh.
     */
    explicit Mesh(const size_t dimension,
                  const CoordinateSystem coordinate_system);

  public:

    /**
     * Establish the relationships between cells and their neighbors.
     *
     * Connectivity is established by
     * 1. For each Vertex, determine which Cell objects the Vertex belongs to.
     * 2. For each Face of each Cell, compare vertex ids to those of adjacent cells
     *    as defined by the step 1.
     * 3. If vertex ids match, set the neighbor properties.
     *
     * \note This routine may be quite expensive. Only use this when connectivity
     *       cannot be established a-priori. Generally, this routine should only be
     *       utilized for unstructured meshes.
     */
    void compute_geometric_info();


    /** Compute the geometric properties of the cells and faces. */
    void establish_connectivity();
  };

}

#endif //MESH_H

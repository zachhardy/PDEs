#ifndef CELL_H
#define CELL_H

#include "cartesian_vector.h"
#include "face.h"

#include <vector>
#include <cstddef>


namespace Grid
{
  /**
   * Available cell geometries.
   */
  enum class CellType
  {
    SLAB = 0, ///< 1D Cartesian geometry.
    ANNULUS = 1, ///< 1D cylindrical geometry.
    SHELL = 2, ///< 1D spherical geometry.
    QUADRILATERAL = 3 ///< 2D quadrilateral.
  };


  /**
   * Return the cell type as a string.
   */
  std::string
  cell_type_str(const CellType cell_type);


  /**
   * A class representing a cell on a mesh.
   *
   * A cell is defined as a <tt>dim</tt>-dimensional object bound by
   * <tt>dim - 1</tt>-dimensional face objects. The cell type largely depends on
   * the dimension and the coordinate system type. Each cell is uniquely
   * identified by its \p id and can store a \p material_id to identify
   * material properties that live on the cell. Examples of cell per dimension
   * and coordinate system are:
   *  - 1D:
   *      - Cartesian:   Slab
   *      - Cylindrical: Annulus
   *      - Spherical:   Shell
   *  - 2D:
   *      - Cartesian:   Polygon
   *          - Traingle
   *          - Quadrilaterals
   *  - 3D:
   *      - Cartesian:   Polyhedron
   *          - Tetrahedra
   *          - Polyhedra
   */
  class Cell
  {
  public:
    /**
     * The cell geometry. This is used to ensure the correct cell volume and
     * face areas are computed upon initialization.
     */
    const CellType type;

    /**
     * A unique cell ID used to identify the cell. This is often used for
     * defining/mapping degrees of freedom.
     */
    size_t id;

    /**
     * A material ID used to ensure the correct material properties are used
     * when performing cell-wise computations.
     */
    unsigned int material_id = -1;

    /**
     * The coordinate of the center of the cell.
     */
    Centroid centroid;

    /**
     * The volume of the cell.
     */
    double volume = 0.0;

    /**
     * A list of the vertex IDs that live on the cell.
     */
    std::vector<size_t> vertex_ids;

    /**
     * A list of the faces that bound the cell. See \ref Face.
     */
    std::vector<Face> faces;

  public:
    /**
     * Construct cell with the specified geometry.
     */
    explicit Cell(const CellType cell_type);

    /**
     * Return the contents of the cell as a string.
     */
    std::string
    str() const;
  };


  /**
   * Insert the contents of a cell into an output stream.
   */
  std::ostream&
  operator<<(std::ostream& os, const Cell& cell);
}
#endif //CELL_H

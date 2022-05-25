#ifndef CELL_H
#define CELL_H

#include "point.h"
#include "face.h"

#include <vector>
#include <cstddef>


namespace pdes::Grid
{

/**
 * Available cell types.
 */
enum class CellType
{
  SLAB    = 0,  ///< 1D Cartesian geometry.
  ANNULUS = 1,  ///< 1D cylindrical geometry.
  SHELL   = 2   ///< 1D spherical geometry.
};


/**
 * Return the cell type as a string.
 */
std::string
cell_type_str(const CellType cell_type);


//######################################################################

/**
 * A class representing a cell on a Mesh.
 *
 * A Cell is defined as a <tt>dim</tt>-dimensional object bound by
 * <tt>dim - 1</tt>-dimensional Face objects. The Cell type largely depends on
 * the dimension and the coordinate system type. Each cell is uniquely
 * identified by its \p id and can store a \p material_id to identify
 * material properties that live on the Cell. Examples of cell per dimension
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
  const CellType type;
  size_t id;
  int material_id = -1;

  Centroid  centroid;
  double volume = 0.0;

  std::vector<size_t> vertex_ids;
  std::vector<Face> faces;

public:
  /**
   * Construct an empty cell of the specified type.
   */
  explicit
  Cell(const CellType cell_type);

  /**
   * Return the cell as a string.
   */
  std::string
  str() const;
};

/**
 * Insert a cell into an output stream.
 */
std::ostream&
operator<<(std::ostream& os, const Cell& cell);

}
#endif //CELL_H

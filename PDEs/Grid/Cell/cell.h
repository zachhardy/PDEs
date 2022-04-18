#ifndef CELL_H
#define CELL_H

#include "../grid_structs.h"
#include "../Face/face.h"

#include <vector>

/// The available Cell types.
enum class CellType
{
  SLAB = 1,     ///< 1D Cartesian geometry.
  ANNULUS = 2,  ///< 1D cylindrical geometry.
  SHELL = 3     ///< 1D spherical geometry.
};

std::string cell_type_name(const CellType cell_type);


//######################################################################
/**
 * \brief A class representing a cell on a Mesh.
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

  const CellType type; ///< The type of the cell.
  size_t id;           ///< A unique identifier tag.
  int material_id = 0; ///< A tag to map to material properties.

  // Geometric information
  Centroid  centroid;   ///< The centroid of the cell.
  double volume = 0.0; ///< The volume of the cell.

  std::vector<size_t> vertex_ids; ///< The vertex IDs that belong to the cell.
  std::vector<Face> faces;        ///< The faces that bound the cell.

public:

  /**
   * \brief Default constructor.
   * \param cell_type The Cell type.
   */
  explicit Cell(const CellType cell_type)
    : type(cell_type)
  {}

  Cell(const Cell& other);            ///< Copy constructor.
  Cell(Cell&& other);                 ///< Move constructor.
  Cell& operator=(const Cell& other); ///< Assignment operator.

public:
  /// Print the contents of the Cell to a string.
  std::string to_string() const;
};

#endif //CELL_H

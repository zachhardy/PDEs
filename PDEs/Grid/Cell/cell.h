#ifndef CELL_H
#define CELL_H

#include "../grid_structs.h"
#include "../Face/face.h"

#include <vector>

/// Different Cell types avaialable.
enum class CellType
{
  SLAB = 1,     ///< 1D Cartesian geometry.
  ANNULUS = 2,  ///< 1D cylindrical geometry.
  SHELL = 3     ///< 1D spherical geometry.
};


/// Get the Cell cell_type as a string.
std::string cell_type_name(const CellType cell_type);


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
public:   /*---------- Public Attributes ----------*/

  const CellType type; ///< The type of the cell.
  size_t id; ///< A unique identifier tag.
  int material_id = 0; ///< A tag to identify the correct material properties.

  // Geometric information
  Centroid  centroid;
  double volume = 0.0;

  /// A mapping to the Vertex objects contained within the Mesh object.
  std::vector<size_t> vertex_ids;
  std::vector<Face> faces;

public:   /*---------- Constructors, Destructors, and Assignments ----------*/

  /**
   * \brief Default constructor.
   * \param cell_type The Cell type.
   */
  explicit Cell(const CellType cell_type)
    : type(cell_type)
  {}

  /// Copy constructor.
  Cell(const Cell& other);

  /// Move constructor.
  Cell(Cell&& other);

  /// Assignment operator.
  Cell& operator=(const Cell& other);

public:   /*---------- Routines ----------*/

  /// Get the Cell information as a string.
  std::string to_string() const;
};

#endif //CELL_H

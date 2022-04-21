#ifndef CELL_H
#define CELL_H

#include "../grid_structs.h"
#include "../Face/face.h"

#include <vector>

namespace grid
{

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
  const CellType type;
  size_t id;
  int material_id = 0;

  Centroid  centroid;
  double volume = 0.0;

  std::vector<size_t> vertex_ids;
  std::vector<Face> faces;

public:
  explicit Cell(const CellType cell_type)
    : type(cell_type)
  {}

  Cell(const Cell& other);
  Cell(Cell&& other);
  Cell& operator=(const Cell& other);

public:
  std::string to_string() const;
};

}
#endif //CELL_H

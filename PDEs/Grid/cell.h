#ifndef CELL_H
#define CELL_H

#include "cartesian_vector.h"
#include "face.h"

#include <vector>
#include <cstddef>


namespace PDEs
{
  namespace Grid
  {
    /// Available cell geometries.
    enum class CellType
    {
      SLAB = 0, ///< 1D Cartesian geometry
      ANNULUS = 1, ///< 1D cylindrical geometry
      SHELL = 2, ///< 1D spherical geometry
      QUADRILATERAL = 3, ///< 2D quadrilateral
    };


    /// Return the cell type as a string.
    std::string
    cell_type_str(const CellType cell_type);

    //######################################################################

    /**
     * A class representing a cell on a mesh.
     *
     * A cell is defined as a <tt>dim</tt>-dimensional object bound by
     * (<tt>dim</tt> - 1)-dimensional face objects. The cell type largely
     * depends on the dimension and the coordinate system type. Each cell is
     * uniquely identified by its \p id and can store a \p material_id to
     * identify material properties that live on the cell.
     */
    class Cell
    {
    public:
      const CellType type;

      size_t id;
      unsigned int material_id = -1;
      std::vector<size_t> ijk;

      Centroid centroid;
      double volume = 0.0;

      std::vector<size_t> vertex_ids;
      std::vector<Face> faces;

    public:
      /// Construct cell with the specified geometry.
      explicit Cell(const CellType cell_type);

      /// Return the contents of the cell as a string.
      std::string str() const;

      /// Insert the contents of a cell into an output stream.
      friend std::ostream&
      operator<<(std::ostream& os, const Cell& cell);
    };


    /// Insert the contents of a cell into an output stream.
    std::ostream& operator<<(std::ostream& os, const Cell& cell);
  }
}
#endif //CELL_H

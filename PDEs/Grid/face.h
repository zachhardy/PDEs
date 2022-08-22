#ifndef FACE_H
#define FACE_H

#include "cartesian_vector.h"

#include <vector>
#include <cstddef>


namespace Grid
{

  /**
   * A class that represents a face on a cell.
   *
   * A face is defined as a <tt>dim - 1</tt>-dimensional object which, in
   * a collection, bounds a <tt>dim</tt>-dimensional cell. Faces in various
   * dimensions are:
   *  -   1D: Vertex, or 1 vertex
   *  -   2D: Edge, or connection of 2 vertices
   *  -   3D: Surface, or collection of several edges.
   *      -   Triangle
   *      -   Quadrilateral
   */
  class Face
  {
  public:
    /**
     * A list of the vertex IDs that live on this face.
     */
    std::vector<size_t> vertex_ids;

    /**
     * On interior faces, this stores the global ID of the cell opposite the
     * face. On boundary faces, this stores the boundary ID.
     */
    size_t neighbor_id = 0;

    /**
     * A flag for whether a face has a neighboring cell or is on boundary.
     */
    bool has_neighbor = false;

    /**
     * An outward normal vector from the perspective of the cell the face
     * belongs to.
     */
    Normal normal;

    /**
     * The coordinate of the center of the face.
     */
    Centroid centroid;

    /**
     * The area of the face.
     */
    double area = 0.0;

    /**
     * Return the contents of the face as a string.
     */
    std::string
    str() const;
  };


  /**
   * Insert the contents of the face into an output stream.
   */
  std::ostream&
  operator<<(std::ostream&, const Face& face);
}
#endif //FACE_H

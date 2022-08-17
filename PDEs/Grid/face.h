#ifndef FACE_H
#define FACE_H

#include "point.h"

#include <vector>
#include <cstddef>


namespace Grid
{

  /**
   * A class that represents a face on a Cell.
   *
   * A face is defined as a <tt>dim - 1</tt>-dimensional object which, in
   * a collection, bounds a <tt>dim</tt>-dimensional cell. face objects in
   * various dimensions are:
   *  -   1D: Vertex
   *  -   2D: Edge, or connection of 2 vertices
   *  -   3D: Surface, or collection of several edges.
   *      -   Triangle
   *      -   Quadrilateral
   */
  class Face
  {
  public:
    std::vector<size_t> vertex_ids;

    /**
     * The ID of the cell opposite this face in the direction of its normal if
     * an interior face. If a boundary face, the boundary ID.
     */
    size_t neighbor_id = 0;
    bool has_neighbor = false;

    Normal normal;
    Centroid centroid;
    double area = 0.0;

    std::string str() const;
  };


  std::ostream& operator<<(std::ostream&, const Face& face);
}
#endif //FACE_H

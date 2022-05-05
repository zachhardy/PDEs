#ifndef FACE_H
#define FACE_H

#include "grid_structs.h"

#include <vector>
#include <cinttypes>


namespace grid
{

/**
 * \brief A class that represents a face on a Cell.
 *
 * A Face is defined as a <tt>dim - 1</tt>-dimensional object which, in
 * a collection, bounds a <tt>dim</tt>-dimensional Cell. Face objects in various
 * dimensions are:
 *  -   1D: Vertex
 *  -   2D: Edge, or connection of 2 vertices
 *  -   3D: Surface, or collection of several edges.
 *      -   Triangle
 *      -   Quadrilateral
 */
class Face
{
public:
  std::vector<uint64_t> vertex_ids;

  bool has_neighbor = false;
  uint64_t neighbor_id = 0;  ///< The neighbor cell or boundary ID.

  Normal normal;
  Centroid centroid;
  double area = 0.0;

public:
  Face() = default;
  Face(const Face& other);
  Face(Face&& other);
  Face& operator=(const Face& other);

public:
  std::string to_string() const;
};

}
#endif //FACE_H

#ifndef FACE_H
#define FACE_H

#include "grid_structs.h"
#include <vector>

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
public:   /*---------- Public Members ----------*/

  /// A mapping to the Vertex objects contained within the Mesh object.
  std::vector<size_t> vertex_ids;

  // Neighbor information
  /// Flag for whether a neighbor exists. If false, this is a boundary.
  bool has_neighbor = false;
  /// The `id` of the neighbor. For boundary cells, this is a boundary ID.
  size_t neighbor_id = 0;

  // Geometric information
  Normal normal;   ///< The outward-pointing normal vector.
  Centroid centroid;
  double area = 0.0;

public:   /*---------- Constructors, Destructors, and Assignments ----------*/

  /// Default constructor.
  Face() = default;

  /// Copy constructor.
  Face(const Face& other);

  /// Move constructor.
  Face(Face&& other);

  /// Default destructor.
  ~Face() = default;

  /// Assignment operator.
  Face& operator=(const Face& other);

public:   /*---------- Routines ----------*/

  /// Get the Face information as a string.
  std::string to_string() const;
};


#endif //FACE_H

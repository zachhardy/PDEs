#include "mesh.h"
#include "cell.h"
#include "macros.h"

#include <set>
#include <cmath>


using namespace Grid;


std::string
coordinate_system_str(const CoordinateSystemType coord_sys)
{
  switch (coord_sys)
  {
    case CoordinateSystemType::CARTESIAN: return "CARTESIAN";
    case CoordinateSystemType::CYLINDRICAL: return "CYLINDRICAL";
    case CoordinateSystemType::SPHERICAL: return "SPHERICAL";
    default: return "UNDEFINED";
  }
}


//######################################################################


Mesh::Mesh(const size_t dimension,
           const CoordinateSystemType coordinate_system)
  : dim(dimension), coord_sys(coordinate_system)
{}


void
Mesh::establish_connectivity()
{
  std::cout << "Establishing cell connectivity.\n";

  /* Determine the cells which contain a specific vertex. This is done using a
   * list where each element represents a vertex and its value contains a set
   * which holds the unique cell ids in which the vertex is found. */
  const size_t n_vertices = vertices.size();
  std::vector<std::set<size_t>> vertex_cell_map(n_vertices);
  for (const auto& cell : cells)
    for (const auto& v_id : cell.vertex_ids)
      vertex_cell_map.at(v_id).insert(cell.id);

  /* Establish connectivity by going through each face of each cell and
   * finding faces on neighboring cells (as defined by vertex_cell_map) which
   * share vertex ids. */
  for (auto& cell : cells)
  {
    for (auto& face : cell.faces)
    {
      /* If there is a neighbor, this face has already been
       * encountered, continue to the next. */
      if (face.has_neighbor) continue;

      // Get the vertex ids for this face
      const std::set<size_t> v_ids(face.vertex_ids.begin(),
                                   face.vertex_ids.end());

      // Use the vertex_cell_map to find neighbor cells
      std::set<size_t> cells_to_search;
      for (const auto& v_id : face.vertex_ids)
        for (const auto& c_id : vertex_cell_map.at(v_id))
          if (c_id != cell.id)
            cells_to_search.insert(c_id);

      // Search the neighbor cells for matching faces
      for (const auto& adj_c_id : cells_to_search)
      {
        auto& adj_cell = cells.at(adj_c_id);

        // Go through neighbor cell faces
        for (auto& adj_face : adj_cell.faces)
        {
          // If neighbor has already been set, it is not a neighbor
          if (adj_face.has_neighbor) continue;

          // Get the adjacent face vertex ids
          const std::set<size_t> adj_v_ids(adj_face.vertex_ids.begin(),
                                           adj_face.vertex_ids.end());

          /* If this face and the adjacent face share vertex ids, they
           * are the same. Define neighbor properties. */
          if (v_ids == adj_v_ids)
          {
            face.neighbor_id = adj_c_id;
            face.has_neighbor = true;

            adj_face.neighbor_id = cell.id;
            adj_face.has_neighbor = true;

            goto found_neighbor;
          }
        }//for adjacent face
      }//for adjacent cells

      found_neighbor:;

    }//for face
  }//for cell
}


//######################################################################


void
Mesh::compute_geometric_info()
{
  std::cout << "Computing geometric information on cells and faces.\n";

  Assert(dim < 3, "Only 1D and 2D meshes are implemented.");

  // Loop over cells
  for (auto& cell : cells)
  {

    // Compute cell centroid
    cell.centroid *= 0.0;
    for (auto& v_id : cell.vertex_ids)
      cell.centroid += vertices[v_id];
    cell.centroid /= static_cast<double>(cell.vertex_ids.size());


    // Compute cell volume
    if (cell.vertex_ids.size() == 2)
    {
      const auto& v1 = vertices[cell.vertex_ids[1]].z;
      const auto& v0 = vertices[cell.vertex_ids[0]].z;

      if (cell.type == CellType::SLAB)
        cell.volume = v1 - v0;
      else if (cell.type == CellType::ANNULUS)
        cell.volume = M_PI*(v1*v1 - v0*v0);
      else if (cell.type == CellType::SHELL)
        cell.volume = 4.0/3.0*M_PI*(v1*v1*v1 - v0*v0*v0);
      else
        throw std::runtime_error("Unexpected 1D cell type.");
    }//if 1D
    else if (cell.vertex_ids.size() == 4)
    {
      assert(cell.type == CellType::QUADRILATERAL);

      const auto& vbl = vertices[cell.vertex_ids[0]];
      const auto& vtr = vertices[cell.vertex_ids[2]];
      cell.volume = (vtr.x - vbl.x) * (vtr.y - vbl.y);
    }//if 2D quad


    // Loop over faces
    for (auto& face : cell.faces)
    {
      // Compute face centroids
      face.centroid *= 0.0;
      for (auto& v_id : face.vertex_ids)
        face.centroid += vertices[v_id];
      face.centroid /= static_cast<double>(face.vertex_ids.size());

      // Compute face area
      if (face.vertex_ids.size() == 1)
      {
        const auto& v = vertices[face.vertex_ids[0]].z;

        if (cell.type == CellType::SLAB)
          face.area = 1.0;
        else if (cell.type == CellType::ANNULUS)
          face.area = 2.0*M_PI * v;
        else if (cell.type == CellType::SHELL)
          face.area = 4.0*M_PI * v*v;
      }// if 1D
      else if (face.vertex_ids.size() == 2)
      {
        assert(coord_sys == CoordinateSystemType::CARTESIAN);

        const auto& v0 = vertices[face.vertex_ids[0]];
        const auto& v1 = vertices[face.vertex_ids[1]];
        face.area = v1.distance(v0);
      }
    }//for faces
  }//for cell
}

#include "mesh.h"
#include "../Cell/cell.h"

#include <set>

/**
 * \brief Establish the relationships between cells and their neighbors.
 *
 * Connectivity is established by
 * 1. For each Vertex, determine which Cell objects the Vertex belongs to.
 * 2. For each Face of each Cell, compare vertex ids to those of adjacent cells
 *    as defined by the step 1.
 * 3. If vertex ids match, set the neighbor properties.
 *
 * \note This routine may be quite expensive. Only use this when connectivity
 *       cannot be established a-priori. Generally, this routine should only be
 *       utilized for unstructured meshes.
 */
void grid::Mesh::establish_connectivity()
{
  std::cout << "Establishing cell connectivity.\n";

  // Determine the cells which contain a specific vertex. This is done using a
  // list where each element represents a vertex and its value contains a set
  // which holds the unique cell ids in which the vertex is found.
  const size_t n_vertices = vertices.size();
  std::vector<std::set<size_t>> vertex_cell_map(n_vertices);
  for (const auto& cell : cells)
    for (const auto& v_id : cell->vertex_ids)
      vertex_cell_map.at(v_id).insert(cell->id);

  // Establish connectivity by going through each face of each cell and
  // finding faces on neighboring cells (as defined by vertex_cell_map) which
  // share vertex ids.

  // Go through cells
  for (auto& cell : cells)
  {
    // Go through faces on cell
    for (auto& face : cell->faces)
    {
      // If there is a neighbor, this face has already been encountered
      if (face.has_neighbor) continue;

      // Get the vertex ids for this face for comparison
      const std::set<size_t> v_ids(face.vertex_ids.begin(),
                                   face.vertex_ids.end());

      // Use the vertex_cell_map to find neighbor cells
      std::set<size_t> cells_to_search;
      for (const auto& v_id : face.vertex_ids)
        for (const auto& c_id : vertex_cell_map.at(v_id))
          if (c_id != cell->id)
            cells_to_search.insert(c_id);

      // Search the neighbor cells for matching faces
      for (const auto& adj_c_id : cells_to_search)
      {
        auto& adj_cell = cells.at(adj_c_id);

        // Go through neighbor cell faces
        for (auto& adj_face : adj_cell->faces)
        {
          // If neighbor has already been set, it is not a neighbor
          if (adj_face.has_neighbor) continue;

          // Get the adjacent face vertex ids
          const std::set<size_t> adj_v_ids(adj_face.vertex_ids.begin(),
                                           adj_face.vertex_ids.end());

          // If this face and the adjacent face share vertex ids, they
          // are the same. Define neighbor properties.
          if (v_ids == adj_v_ids)
          {
            face.neighbor_id = adj_c_id;
            face.has_neighbor = true;

            adj_face.neighbor_id = cell->id;
            adj_face.has_neighbor = true;

            goto found_neighbor;
          }
        }//for adjacent face
      }//for adjacent cells
      found_neighbor:;
    }//for face
  }//for cell

  std::cout << "Finished establishing cell connectivity.\n";
}

//######################################################################

/// Compute the geometric properties of the cells and faces.
void  grid::Mesh::compute_geometric_info()
{
  std::cout << "Computing geometric information on cells and faces.\n";

  // Go through each cell
  for (auto& cell : cells)
  {
    // Compute cell centroid
    cell->centroid *= 0.0;
    for (auto& v_id : cell->vertex_ids)
      cell->centroid += vertices[v_id];
    cell->centroid /= static_cast<double>(cell->vertex_ids.size());

    // Compute cell volume
    if (cell->vertex_ids.size() == 2)
    {
      auto cell_type = cell->type;
      auto& v1 = vertices[cell->vertex_ids[1]].z;
      auto& v0 = vertices[cell->vertex_ids[0]].z;

      // Compute vo
      switch (cell_type)
      {
        case CellType::SLAB:
        { cell->volume = v1 - v0; break; }
        case CellType::ANNULUS:
        { cell->volume = M_PI*(v1*v1 - v0*v0); break; }
        case CellType::SHELL:
        { cell->volume = 4.0/3.0*M_PI*(v1*v1*v1 - v0*v0*v0); break; }
      }
    }//if 1D
    else
    {
      std::stringstream err;
      err << "Mesh::" << __FUNCTION__ << ": "
          << "Only 1D cell volume routines are implemented.";
      throw std::runtime_error(err.str());
    }

    // Go through the cell's faces
    for (auto& face : cell->faces)
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

        switch (coord_sys)
        {
          case CoordinateSystem::CARTESIAN:
          { face.area = 1.0; break; }
          case CoordinateSystem::CYLINDRICAL:
          { face.area = 2.0*M_PI * v; break; }
          case CoordinateSystem::SPHERICAL:
          { face.area = 4.0*M_PI * v*v; break; }
        }
      }// if 1D
      else
      {
        std::stringstream err;
        err << "Mesh::" << __FUNCTION__ << ": "
            << "Only 1D face area routines are implemented.";
        throw std::runtime_error(err.str());
      }
    }//for faces
  }//for cell

  std::cout << "Finished computing geometric information on "
               "cells and faces.\n";
}

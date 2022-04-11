#include "mesh.h"
#include <cmath>


/**
 * For each Cell object, the Centroid and volume are computed. For each Face
 * object, the Centroid and area are computed.
 */
void Mesh::compute_geometric_info()
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
        {
          cell->volume = v1 - v0;
          break;
        }
        case CellType::ANNULUS:
        {
          cell->volume = M_PI*(v1*v1 - v0*v0);
          break;
        }
        case CellType::SHELL:
        {
          cell->volume = 4.0/3.0*M_PI*(v1*v1*v1 - v0*v0*v0);
          break;
        }
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
          {
            face.area = 1.0;
            break;
          }
          case CoordinateSystem::CYLINDRICAL:
          {
            face.area = 2.0*M_PI * v;
            break;
          }
          case CoordinateSystem::SPHERICAL:
          {
            face.area = 4.0*M_PI * v*v;
            break;
          }
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

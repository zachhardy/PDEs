#include "mesh.h"
#include "cell.h"

#include <set>
#include <cmath>

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>
#include <iomanip>


namespace PDEs::Grid
{


  std::string
  coordinate_system_str(const CoordinateSystemType coord_sys)
  {
    switch (coord_sys)
    {
      case CoordinateSystemType::CARTESIAN:   return "CARTESIAN";
      case CoordinateSystemType::CYLINDRICAL: return "CYLINDRICAL";
      case CoordinateSystemType::SPHERICAL:   return "SPHERICAL";
      default:                                return "UNDEFINED";
    }
  }

//######################################################################

  Mesh::Mesh(const unsigned int dimension,
             const CoordinateSystemType coordinate_system) :
      dimension(dimension), coordinate_system(coordinate_system)
  {}

//######################################################################

  void Mesh::establish_connectivity()
  {
    std::cout << "Establishing cell connectivity.\n";

    // Determine the cells which contain a specific vertex. This is done using a
    // list where each element represents a vertex and its value contains a set
    // which holds the unique cell ids in which the vertex is found.
    const size_t n_vertices = vertices.size();
    std::vector<std::set<size_t>> vertex_cell_map(n_vertices);
    for (const auto& cell: cells)
      for (const auto& v_id: cell.vertex_ids)
        vertex_cell_map.at(v_id).insert(cell.id);

    /* Establish connectivity by going through each face of each cell and
     * finding faces on neighboring cells (as defined by vertex_cell_map) which
     * share vertex ids. */
    for (auto& cell: cells)
    {
      for (auto& face: cell.faces)
      {
        /* if there is a neighbor, this face has already been
         * encountered, continue to the next. */
        if (face.has_neighbor) continue;

        // get the vertex ids for this face
        const std::set<size_t> v_ids(face.vertex_ids.begin(),
                                     face.vertex_ids.end());

        // use the vertex_cell_map to find neighbor cells
        std::set<size_t> cells_to_search;
        for (const auto& v_id: face.vertex_ids)
          for (const auto& c_id: vertex_cell_map.at(v_id))
            if (c_id != cell.id)
              cells_to_search.insert(c_id);

        // search the neighbor cells for matching faces
        for (const auto& adj_c_id: cells_to_search)
        {
          auto& adj_cell = cells.at(adj_c_id);

          // go through neighbor cell faces
          for (auto& adj_face: adj_cell.faces)
          {
            // if neighbor has already been set, it is not a neighbor
            if (adj_face.has_neighbor) continue;

            // get the adjacent face vertex ids
            const std::set<size_t> adj_v_ids(adj_face.vertex_ids.begin(),
                                             adj_face.vertex_ids.end());

            // set neighbor properties if the faces are the same
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


  void Mesh::compute_geometric_info()
  {
    assert(dimension < 3);
    std::cout << "Computing geometric information on cells and faces.\n";

    // loop over cells
    for (auto& cell: cells)
    {

      // compute cell centroid
      cell.centroid *= 0.0;
      for (auto& v_id: cell.vertex_ids)
        cell.centroid += vertices[v_id];
      cell.centroid /= static_cast<double>(cell.vertex_ids.size());


      // compute cell volume
      if (cell.vertex_ids.size() == 2)
      {
        const auto& v1 = vertices[cell.vertex_ids[1]].z();
        const auto& v0 = vertices[cell.vertex_ids[0]].z();

        if (cell.type == CellType::SLAB)
          cell.volume = v1 - v0;
        else if (cell.type == CellType::ANNULUS)
          cell.volume = M_PI * (v1 * v1 - v0 * v0);
        else if (cell.type == CellType::SHELL)
          cell.volume = 4.0 / 3.0 * M_PI * (v1 * v1 * v1 - v0 * v0 * v0);
        else
          throw std::runtime_error("Unexpected 1D cell type.");
      }//if 1D

      else if (cell.vertex_ids.size() == 4)
      {
        assert(cell.type == CellType::QUADRILATERAL);

        const auto& vbl = vertices[cell.vertex_ids[0]];
        const auto& vtr = vertices[cell.vertex_ids[2]];
        cell.volume = (vtr.x() - vbl.x()) * (vtr.y() - vbl.y());
      }//if 2D quad


      // loop over faces
      for (auto& face: cell.faces)
      {
        // compute face centroids
        face.centroid *= 0.0;
        for (auto& v_id: face.vertex_ids)
          face.centroid += vertices[v_id];
        face.centroid /= static_cast<double>(face.vertex_ids.size());

        // compute face area
        if (face.vertex_ids.size() == 1)
        {
          const auto& v = vertices[face.vertex_ids[0]].z();

          if (cell.type == CellType::SLAB)
            face.area = 1.0;
          else if (cell.type == CellType::ANNULUS)
            face.area = 2.0 * M_PI * v;
          else if (cell.type == CellType::SHELL)
            face.area = 4.0 * M_PI * v * v;
        }// if 1D

        else if (face.vertex_ids.size() == 2)
        {
          assert(coordinate_system == CoordinateSystemType::CARTESIAN);
          const auto& v0 = vertices[face.vertex_ids[0]];
          const auto& v1 = vertices[face.vertex_ids[1]];
          face.area = v1.distance(v0);
        }// if 2D Cartesian
      }//for faces
    }//for cell

    std::cout << "Done computing geometric information on cells and faces.\n";
  }

//######################################################################

  void Mesh::write_ascii(const std::string output_directory,
                         const std::string file_prefix) const
  {
    if (not std::filesystem::is_directory(output_directory))
      std::filesystem::create_directory(output_directory);
    assert(std::filesystem::is_directory(output_directory));

    std::string filepath = output_directory + "/" + file_prefix;
    if (filepath.find(".") != std::string::npos)
      filepath = file_prefix.substr(0, file_prefix.rfind("."));
    filepath = filepath + ".txt";

    std::ofstream file(filepath,
                       std::ofstream::out |
                       std::ofstream::trunc);
    assert(file.is_open());

    file << "##################################################\n"
            "# Mesh File\n"
            "# _________"
            "# For Each Vertex:\n"
            "#   VertexID  x  y  z\n"
            "# For Each Cell:\n"
            "#   CellID  MaterialID  <VertexIDs>\n"
         << "##################################################\n";

    file << "vertices\n";
    for (size_t i = 0; i < vertices.size(); ++i)
      file << std::left
           << std::setw(5) << i << "  "
           << std::setw(8) << vertices[i].x() << "  "
           << std::setw(8) << vertices[i].y() << "  "
           << std::setw(8) << vertices[i].z() << "\n";

    file << "cells\n";
    for (const auto& cell: cells)
    {
      file << std::left << std::setw(5) << cell.id;
      for (const auto& vid: cell.vertex_ids)
        file << "  " << std::setw(5) << vid;
      file << "\n";
    }

    file.close();
  }


  void Mesh::write_binary(const std::string output_directory,
                          const std::string file_prefix) const
  {

    if (not std::filesystem::is_directory(output_directory))
      std::filesystem::create_directory(output_directory);
    assert(std::filesystem::is_directory(output_directory));

    std::string filepath = output_directory + "/" + file_prefix;
    if (filepath.find(".") != std::string::npos)
      filepath = file_prefix.substr(0, file_prefix.rfind("."));
    filepath = filepath + ".txt";

    std::ofstream file(filepath,
                       std::ofstream::binary |
                       std::ofstream::out |
                       std::ofstream::trunc);
    assert(file.is_open());

    // write the header_info
    int size = 300;
    std::string header_info =
        "Mesh File\nHeader size: " + std::to_string(size) + "bytes\n";
    header_info +=
        "Structure(type - info)\n"
        "Vertices:\n"
        "  uint64_t  n_vertices\n"
        "  Each Verex:\n"
        "    uint64_t vertex_id\n"
        "    double   x_position\n"
        "    double   y_position\n"
        "    double   z_position\n"
        "Cells:\n"
        "  uint64_t  n_cells\n"
        "  Each Cell:\n"
        "    uint64_t      cell_id\n"
        "    unsigned int  material_id\n"
        "    Each Vertex ID\n"
        "      uint64_t  vertex_id\n";

    int header_size = (int) header_info.length();

    char header_bytes[size];
    memset(header_bytes, '-', size);
    strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
    header_bytes[size - 1] = '\0';

    const auto n_vertices = vertices.size();
    const auto n_cells = cells.size();

    file << header_bytes;
    file.write((char*) &n_vertices, sizeof(uint64_t));
    for (uint64_t i = 0; i < vertices.size(); ++i)
    {
      const auto x = vertices[i].x();
      const auto y = vertices[i].y();
      const auto z = vertices[i].z();

      file.write((char*) &i, sizeof(uint64_t));
      file.write((char*) &x, sizeof(double));
      file.write((char*) &y, sizeof(double));
      file.write((char*) &z, sizeof(double));
    }

    file.write((char*) &n_cells, sizeof(uint64_t));
    for (const auto& cell: cells)
    {
      file.write((char*) &cell.id, sizeof(uint64_t));
      file.write((char*) &cell.material_id, sizeof(unsigned int));
      for (const auto& vid: cell.vertex_ids)
        file.write((char*) &vid, sizeof(uint64_t));
    }
  }
}
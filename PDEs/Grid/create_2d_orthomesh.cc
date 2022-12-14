#include "ortho_grids.h"

#include <cassert>
#include <algorithm>


namespace PDEs::Grid
{

  std::shared_ptr<Mesh>
  create_2d_orthomesh(const std::vector<double>& x_vertices,
                      const std::vector<double>& y_vertices,
                      const bool verbose)
  {
    assert(!x_vertices.empty());
    assert(!y_vertices.empty());
    std::cout << "Creating 2D orthogonal mesh from vertices.\n";

    // create the mesh
    auto mesh = std::make_shared<Mesh>(2, CoordinateSystemType::CARTESIAN);

    // count the vertices in each dimension
    const size_t n_x = x_vertices.size();
    const size_t n_y = y_vertices.size();

    // initialize the vertices
    mesh->vertices.reserve(n_x * n_y);

    // define the vertex map
    std::vector<std::vector<size_t>> vmap;
    vmap.resize(n_y, std::vector<size_t>(n_x));
    for (size_t i = 0; i < n_y; ++i)
      for (size_t j = 0; j < n_x; ++j)
      {
        vmap[i][j] = mesh->vertices.size();
        mesh->vertices.emplace_back(x_vertices[j], y_vertices[i]);
      }

    // define mesh boundary points
    auto x_min = *std::min_element(x_vertices.begin(), x_vertices.end());
    auto x_max = *std::max_element(x_vertices.begin(), x_vertices.end());
    auto y_min = *std::min_element(y_vertices.begin(), y_vertices.end());
    auto y_max = *std::max_element(y_vertices.begin(), y_vertices.end());

    // reference normal vectors
    const Normal ihat = CartesianVector(1.0, 0.0, 0.0);
    const Normal jhat = CartesianVector(0.0, 1.0, 0.0);
    const Normal khat = CartesianVector(0.0, 0.0, 1.0);

    // create cells along rows, then columns
    for (size_t i = 0; i < n_y - 1; ++i)
      for (size_t j = 0; j < n_x - 1; ++j)
      {
        Cell cell(CellType::QUADRILATERAL);
        cell.id = i * (n_x - 1) + j;

        // vertices are numbered from lower-left, counter-clockwise
        cell.vertex_ids = {vmap[i][j], vmap[i][j + 1],
                           vmap[i + 1][j + 1], vmap[i + 1][j]};

        // initialize the faces in the same order of vertices
        for (size_t f = 0; f < 4; ++f)
        {
          Face face;

          // define the face vertices counter-clock-wise
          if (f < 3)
            face.vertex_ids = {cell.vertex_ids[f], cell.vertex_ids[f + 1]};
          else
            face.vertex_ids = {cell.vertex_ids[f], cell.vertex_ids[0]};

          // compute the face outward normal vector
          auto v0 = mesh->vertices[face.vertex_ids[0]];
          auto v1 = mesh->vertices[face.vertex_ids[1]];
          face.normal = khat.cross(v0 - v1).normalize();

          // define the face neighbors
          if (face.normal == -jhat)                 // bottom face
            face.neighbor_id = cell.id - (n_x - 1);
          else if (face.normal == ihat)             // right face
            face.neighbor_id = cell.id + 1;
          else if (face.normal == jhat)             // top face
            face.neighbor_id = cell.id + (n_x - 1);
          else if (face.normal == -ihat)            // left face
            face.neighbor_id = cell.id - 1;
          else
            throw std::runtime_error("Unexpected face normal encountered.");

          // define boundary faces
          if (v0.y() == y_min && v1.y() == y_min)      // bottom face
            face.neighbor_id = 0;
          else if (v0.x() == x_max && v1.x() == x_max) // right face
            face.neighbor_id = 1;
          else if (v0.y() == y_max && v1.y() == y_max) // top face
            face.neighbor_id = 2;
          else if (v0.x() == x_min && v1.x() == x_min) // left face
            face.neighbor_id = 3;
          else
            face.has_neighbor = true;

          // add the face to the cell
          cell.faces.push_back(face);
        }//for f

        // add the cell to the mesh
        mesh->ijk_mapping.push_back({j, i});
        mesh->cells.push_back(cell);
        for (const auto& face: cell.faces)
          if (not face.has_neighbor)
            mesh->boundary_cell_ids.push_back(cell.id);
      }

    // compute the cell and face geometric info
    mesh->compute_geometric_info();

    if (verbose)
    {
      double n_faces = 0;
      for (const auto& cell: mesh->cells)
        for (const auto& face: cell.faces)
          n_faces += (face.has_neighbor) ? 0.5 : 1.0;

      std::cout << "Mesh Details:\n"
                << "\t# of Vertices: " << mesh->vertices.size() << "\n"
                << "\t# of Cells:    " << mesh->cells.size() << "\n"
                << "\t# of Faces:    " << (size_t) n_faces << "\n";
    }
    return mesh;
  }

}
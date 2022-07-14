#include "ortho_grids.h"

#include <cassert>
#include <numeric>


std::shared_ptr<Grid::Mesh>
Grid::create_2d_orthomesh(const std::vector<double> x_vertices,
                          const std::vector<double> y_vertices,
                          const bool verbose)
{
  std::cout << "Creating 2D orthogonal mesh from vertices.\n";

  assert(!x_vertices.empty());
  assert(!y_vertices.empty());

  // Create the mesh
  auto mesh = std::make_shared<Mesh>(2, CoordinateSystemType::CARTESIAN);

  //========================================
  // Creat the vertices
  //========================================

  // Count the vertices in each dimension
  const size_t n_x = x_vertices.size();
  const size_t n_y = y_vertices.size();

  mesh->vertices.reserve(n_x * n_y);

  // Define the vertex map. This is a structure that defines the mapping of
  std::vector<std::vector<size_t>> vmap;
  vmap.resize(n_y, std::vector<size_t>(n_x));
  for (size_t i = 0; i < n_y; ++i)
    for (size_t j = 0; j < n_x; ++j)
    {
      vmap[i][j] = mesh->vertices.size();
      mesh->vertices.emplace_back(x_vertices[j], y_vertices[i]);
    }

  //========================================
  // Create the cells
  //========================================

  // Reference normal vectors
  Normal ihat = Point(1.0, 0.0, 0.0);
  Normal jhat = Point(0.0, 1.0, 0.0);
  Normal khat = Point(0.0, 0.0, 1.0);

  for (size_t i = 0; i < n_y; ++i)
    for (size_t j = 0; j < n_x; ++j)
    {
      // Initialize the cell
      Cell cell(CellType::QUADRILATERAL);
      cell.id = i*(n_x - 1) + j;

      // Vertices are numbered from lower-left, counter-clockwise
      cell.vertex_ids = {vmap[i][j], vmap[i][j + 1],
                         vmap[i + 1][j + 1], vmap[i + 1][j]};

      // Initialize the faces in the same order of vertices
      for (size_t f = 0; f < 4; ++f)
      {
        Face face;
        if (f < 3)
          face.vertex_ids = {cell.vertex_ids[f], cell.vertex_ids[f + 1]};
        else
          face.vertex_ids = {cell.vertex_ids[f], cell.vertex_ids[0]};


      }


    }



  return mesh;
}

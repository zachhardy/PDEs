#include "grid_structs.h"
#include "mesh.h"

#include <vector>
#include <numeric>


/**
 * \brief Create a 1D mesh from a list of vertices.
 * \param vertices A list of vertex locations.
 * \param coordinate_system The coordinate system type. The default is
 *                          Cartesian coordinates.
 */
std::shared_ptr<grid::Mesh>
grid::create_1d_mesh(const std::vector<double>& vertices,
                     const CoordinateSystem coordinate_system,
                     const bool verbose)
{
  std::cout << "Creating a 1D mesh from vertices...\n";

  // Check for empty input
  if (vertices.empty())
  {
    std::stringstream err;
    err << __FUNCTION__ << ": "
        << "No vertices provided for Mesh creation.";
    throw std::runtime_error(err.str());
  }

  // Create the Mesh
  std::shared_ptr<Mesh> mesh(new Mesh(1, coordinate_system));

  // Count the number of cells
  uint64_t n_cells = vertices.size() - 1;

  // Compute the cell widths
  std::vector<double> widths;
  widths.reserve(n_cells);
  for (uint64_t v = 0; v < n_cells; ++v)
    widths.push_back(vertices[v + 1] - vertices[v]);

  // Initialize the vertices
  mesh->vertices.reserve(vertices.size());
  mesh->vertices.emplace_back(0.0, 0.0, 0.0);

  // Compute the vertices, starting from 0.0
  double current_pos = 0.0;
  for (const auto& width : widths)
  {
    mesh->vertices.emplace_back(0.0, 0.0, current_pos + width);
    current_pos += width;
  }

  // Get the type of cell from the coordinate system
  CellType cell_type;
  switch (coordinate_system)
  {
    case CoordinateSystem::CARTESIAN:
    { cell_type = CellType::SLAB; break; }
    case CoordinateSystem::CYLINDRICAL:
    { cell_type = CellType::ANNULUS; break; }
    case CoordinateSystem::SPHERICAL:
    { cell_type = CellType::SHELL; break; }
  }

  // Create the cells
  for (uint64_t c = 0; c < n_cells; ++c)
  {
    // Initialize the cell and face objects
    std::shared_ptr<Cell> cell(new Cell(cell_type));
    Face left_face, right_face;

    // Define the cell info
    cell->id = c;
    cell->vertex_ids = {c, c + 1};

    // Define the left face info, add to cell
    left_face.vertex_ids = {c};
    left_face.has_neighbor = (c > 0);
    left_face.neighbor_id = (c > 0) ? c - 1 : 0;
    left_face.normal = Normal(0.0, 0.0, -1.0);
    cell->faces.push_back(left_face);

    // Define the right face info, add to cell
    right_face.vertex_ids = {c + 1};
    right_face.has_neighbor = {c < n_cells - 1};
    right_face.neighbor_id = (c < n_cells - 1) ? c + 1 : 1;
    right_face.normal = Normal(0.0, 0.0, 1.0);
    cell->faces.push_back(right_face);

    // Add cells to the mesh
    mesh->cells.push_back(cell);
  }

  // Define the boundary cells
  mesh->boundary_cell_ids = {0, n_cells - 1};

  // Compute the cell and face geometric info
  mesh->compute_geometric_info();

  std::cout << "Done creating mesh.\n";

  if (verbose)
    std::cout << "Mesh Details:\n"
              << "\t# of Vertices: " << mesh->vertices.size() << "\n"
              << "\t# of Cells:    " << mesh->cells.size() << "\n"
              << "\t# of Lines:    " << mesh->cells.size() << "\n";

  return mesh;
}


//######################################################################

/**
 * \brief Create a zoned 1D mesh.
 *
 * Zones are defined by edges, a number of subdivisions (cells), and a
 * material ID. This allows for non-uniform cells throughout the mesh and
 *
 * \param zone_edges The edges of mesh zones. There should be one more
 *                   zone edge than number of zones.
 * \param zone_subdivisions The number of cells per zone.
 * \param material_ids The material ID per zone.
 * \param coordinate_system The coordinate system type. The default is
 *                          Cartesian coordinates.
 */
std::shared_ptr<grid::Mesh>
grid::create_1d_mesh(const std::vector<double>& zone_edges,
                     const std::vector<uint64_t>& zone_subdivisions,
                     const std::vector<int>& material_ids,
                     const CoordinateSystem coordinate_system,
                     const bool verbose)
{
  std::cout << "Creating a 1D mesh from zones...\n";

  // Check for empty inputs
  if (zone_edges.empty() or
      zone_subdivisions.empty() or
      material_ids.empty())
  {
    std::stringstream err;
    err << __FUNCTION__ << ": Some inputs are empty.\n"
        << "zone_edges, zone_subdivisions, and material_ids "
        << "must all be provided.";
    throw std::runtime_error(err.str());
  }

  // Check compatibility
  if (zone_edges.size() != zone_subdivisions.size() + 1 or
      zone_subdivisions.size() != material_ids.size())
  {
    std::stringstream err;
    err << __FUNCTION__ << ": Incompatible inputs.\n"
        << "Ensure that there is one more zone edge than subdivision entry "
        << "and that there are as many subdivision entries ";
    throw std::runtime_error(err.str());
  }

  // Create the mesh
  auto mesh = std::make_shared<Mesh>(1, coordinate_system);

  // Count the number of cells
  uint64_t n_cells = std::accumulate(zone_subdivisions.begin(),
                                   zone_subdivisions.end(), 0);

  // Initialize the vertices
  mesh->vertices.reserve(n_cells);
  mesh->vertices.emplace_back(0.0, 0.0, zone_edges[0]);

  // Define the vertices, loop over each zone, then the cells per zone
  double current_pos = 0.0;
  for (uint64_t z = 0; z < zone_subdivisions.size(); ++z)
  {
    // Define the width of cells in this zone z
    double zone_width = zone_edges[z+1] - zone_edges[z];
    double n_zone_cells = static_cast<double>(zone_subdivisions[z]);
    double cell_width = zone_width / n_zone_cells;

    for (uint64_t c = 0; c < n_zone_cells; ++ c)
    {
      mesh->vertices.emplace_back(0.0, 0.0, current_pos+cell_width);
      current_pos += cell_width;
    }
  }

  // Get the type of cell from the coordinate system
  CellType cell_type;
  switch (coordinate_system)
  {
    case CoordinateSystem::CARTESIAN:
    { cell_type = CellType::SLAB; break; }
    case CoordinateSystem::CYLINDRICAL:
    { cell_type = CellType::ANNULUS; break; }
    case CoordinateSystem::SPHERICAL:
    { cell_type = CellType::SHELL; break; }
  }

  // Create the cells, loop over zones, then cells per zone
  uint64_t count = 0;
  for (uint64_t z = 0; z < zone_subdivisions.size(); ++z)
  {
    for (uint64_t c = 0; c < zone_subdivisions[z]; ++c)
    {
      // Initialize the cell and face objects
      auto cell = std::make_shared<Cell>(cell_type);
      Face left_face, right_face;

      // Define the cell info
      cell->id = count;
      cell->vertex_ids = {count, count + 1};
      cell->material_id = material_ids[z];

      // Define the left face info, add to cell
      left_face.vertex_ids = {count};
      left_face.has_neighbor = (count > 0);
      left_face.neighbor_id = (count > 0) ? count - 1 : 0;
      left_face.normal = Normal(0.0, 0.0, -1.0);
      cell->faces.push_back(left_face);

      // Define the right face info, add to cell
      right_face.vertex_ids = {count + 1};
      right_face.has_neighbor = {count < n_cells - 1};
      right_face.neighbor_id = (count < n_cells - 1) ? count + 1 : 1;
      right_face.normal = Normal(0.0, 0.0, 1.0);
      cell->faces.push_back(right_face);

      // Add the cell to the mesh
      mesh->cells.emplace_back(cell);
      ++count;
    }//for cell
  }//for zone

  // Define the boundary cells
  mesh->boundary_cell_ids = {0, n_cells - 1};

  // Compute the cell and face geometric info
  mesh->compute_geometric_info();


  std::cout << "Done creating mesh.\n";

  if (verbose)
    std::cout << "Mesh Details:\n"
              << "\t# of Vertices: " << mesh->vertices.size() << "\n"
              << "\t# of Cells:    " << mesh->cells.size() << "\n"
              << "\t# of Lines:    " << mesh->cells.size() << "\n";

  return mesh;
}

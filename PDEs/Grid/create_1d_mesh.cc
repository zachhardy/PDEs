#include "ortho_grids.h"
#include "macros.h"

#include <numeric>


std::shared_ptr<Grid::Mesh>
Grid::create_1d_mesh(const std::vector<double> vertices,
                     const CoordinateSystemType coordinate_system,
                     const bool verbose)
{
  std::cout << "Creating a 1D mesh from vertices.\n";

  Assert(!vertices.empty(), "No vertices provided.");

  // Create the Mesh
  auto mesh = std::make_shared<Mesh>(1, coordinate_system);

  // Count the number of cells
  size_t n_cells = vertices.size() - 1;

  // Compute the cell widths
  std::vector<double> widths;
  widths.reserve(n_cells);
  for (size_t v = 0; v < n_cells; ++v)
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
    case CoordinateSystemType::CARTESIAN:cell_type = CellType::SLAB;
      break;

    case CoordinateSystemType::CYLINDRICAL:cell_type = CellType::ANNULUS;
      break;

    case CoordinateSystemType::SPHERICAL:cell_type = CellType::SHELL;
      break;

  }

  // Create the cells
  for (size_t c = 0; c < n_cells; ++c)
  {
    // Initialize the cell and face objects
    Cell cell(cell_type);
    Face left_face, right_face;

    // Define the cell info
    cell.id = c;
    cell.vertex_ids = {c, c + 1};

    // Define the left face info, add to cell
    left_face.vertex_ids = {c};
    left_face.has_neighbor = (c > 0);
    left_face.neighbor_id = (c > 0)? c - 1 : 0;
    left_face.normal = Normal(0.0, 0.0, -1.0);
    cell.faces.push_back(left_face);

    // Define the right face info, add to cell
    right_face.vertex_ids = {c + 1};
    right_face.has_neighbor = {c < n_cells - 1};
    right_face.neighbor_id = (c < n_cells - 1)? c + 1 : 1;
    right_face.normal = Normal(0.0, 0.0, 1.0);
    cell.faces.push_back(right_face);

    // Add cells to the mesh
    mesh->cells.push_back(cell);
  }

  // Define the boundary cells
  mesh->boundary_cell_ids = {0, n_cells - 1};

  // Compute the cell and face geometric info
  mesh->compute_geometric_info();

  if (verbose)
    std::cout << "Mesh Details:\n"
              << "\t# of Vertices: " << mesh->vertices.size() << "\n"
              << "\t# of Cells:    " << mesh->cells.size() << "\n"
              << "\t# of Lines:    " << mesh->cells.size() << "\n";
  return mesh;
}

//######################################################################

std::shared_ptr<Grid::Mesh>
Grid::create_1d_mesh(const std::vector<double> zone_edges,
                     const std::vector<size_t> zone_subdivisions,
                     const std::vector<int> material_ids,
                     const CoordinateSystemType coordinate_system,
                     const bool verbose)
{
  std::cout << "Creating a 1D mesh from zones.\n";

  Assert(!zone_edges.empty(), "No zone edges provided.");
  Assert(!zone_subdivisions.empty(), "No zone subdivisions provided.");
  Assert(!material_ids.empty(), "No material IDs provided.");
  Assert(zone_edges.size() == zone_subdivisions.size() + 1,
         "Incompatible zone edges and zone subdivisions."
         "There must be one more edge than subdivisions entry.");
  Assert(zone_subdivisions.size() == material_ids.size(),
         "Incompatible zone subdivisions and material IDs."
         "There must the same number of zone subdivisions as material IDs.");

  // Create the mesh
  auto mesh = std::make_shared<Mesh>(1, coordinate_system);

  // Count the number of cells
  size_t n_cells = std::accumulate(zone_subdivisions.begin(),
                                   zone_subdivisions.end(), 0);

  // Initialize the vertices
  mesh->vertices.reserve(n_cells);
  mesh->vertices.emplace_back(0.0, 0.0, zone_edges[0]);

  // Define the vertices, loop over each zone, then the cells per zone
  double current_pos = 0.0;
  for (size_t z = 0; z < zone_subdivisions.size(); ++z)
  {
    // Define the width of cells in this zone z
    double zone_width = zone_edges[z + 1] - zone_edges[z];
    double n_zone_cells = static_cast<double>(zone_subdivisions[z]);
    double cell_width = zone_width/n_zone_cells;

    for (size_t c = 0; c < n_zone_cells; ++c)
    {
      mesh->vertices.emplace_back(0.0, 0.0, current_pos + cell_width);
      current_pos += cell_width;
    }
  }

  // Get the type of cell from the coordinate system
  CellType cell_type;
  switch (coordinate_system)
  {
    case CoordinateSystemType::CARTESIAN:cell_type = CellType::SLAB;
      break;

    case CoordinateSystemType::CYLINDRICAL:cell_type = CellType::ANNULUS;
      break;

    case CoordinateSystemType::SPHERICAL:cell_type = CellType::SHELL;
      break;

  }

  // Create the cells, loop over zones, then cells per zone
  size_t count = 0;
  for (size_t z = 0; z < zone_subdivisions.size(); ++z)
  {
    for (size_t c = 0; c < zone_subdivisions[z]; ++c)
    {
      // Initialize the cell and face objects
      Cell cell(cell_type);
      Face left_face, right_face;

      // Define the cell info
      cell.id = count;
      cell.vertex_ids = {count, count + 1};
      cell.material_id = material_ids[z];

      // Define the left face info, add to cell
      left_face.vertex_ids = {count};
      left_face.has_neighbor = (count > 0);
      left_face.neighbor_id = (count > 0)? count - 1 : 0;
      left_face.normal = Normal(0.0, 0.0, -1.0);
      cell.faces.push_back(left_face);

      // Define the right face info, add to cell
      right_face.vertex_ids = {count + 1};
      right_face.has_neighbor = {count < n_cells - 1};
      right_face.neighbor_id = (count < n_cells - 1)? count + 1 : 1;
      right_face.normal = Normal(0.0, 0.0, 1.0);
      cell.faces.push_back(right_face);

      // Add the cell to the mesh
      mesh->cells.emplace_back(cell);
      ++count;
    }//for cell
  }//for zone

  // Define the boundary cells
  mesh->boundary_cell_ids = {0, n_cells - 1};

  // Compute the cell and face geometric info
  mesh->compute_geometric_info();

  if (verbose)
    std::cout << "Mesh Details:\n"
              << "\t# of Vertices: " << mesh->vertices.size() << "\n"
              << "\t# of Cells:    " << mesh->cells.size() << "\n"
              << "\t# of Lines:    " << mesh->cells.size() << "\n";
  return mesh;
}

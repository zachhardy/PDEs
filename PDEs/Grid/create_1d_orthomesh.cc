#include "ortho_grids.h"

#include <cassert>
#include <numeric>


namespace PDEs::Grid
{

  std::shared_ptr<Mesh>
  create_1d_orthomesh(const std::vector<double>& vertices,
                      const CoordinateSystemType coordinate_system,
                      const bool verbose)
  {
    assert(!vertices.empty());
    std::cout << "Creating a 1D mesh from vertices.\n";

    // create the Mesh
    auto mesh = std::make_shared<Mesh>(1, coordinate_system);

    // count the number of cells
    size_t n_cells = vertices.size() - 1;

    // create the vertices
    mesh->vertices.reserve(vertices.size());
    for (const auto& vertex: vertices)
      mesh->vertices.emplace_back(0.0, 0.0, vertex);

    // get the type of cell from the coordinate system
    CellType cell_type;
    switch (coordinate_system)
    {
      case CoordinateSystemType::CARTESIAN:
        cell_type = CellType::SLAB;
        break;

      case CoordinateSystemType::CYLINDRICAL:
        cell_type = CellType::ANNULUS;
        break;

      case CoordinateSystemType::SPHERICAL:
        cell_type = CellType::SHELL;
        break;
    }

    // create the cells
    for (size_t c = 0; c < n_cells; ++c)
    {
      Cell cell(cell_type);
      Face left_face, right_face;

      cell.id = c;
      cell.vertex_ids = {c, c + 1};

      left_face.vertex_ids = {c};
      left_face.has_neighbor = (c > 0);
      left_face.neighbor_id = (c > 0) ? c - 1 : 0;
      left_face.normal = Normal(0.0, 0.0, -1.0);
      cell.faces.push_back(left_face);

      right_face.vertex_ids = {c + 1};
      right_face.has_neighbor = {c < n_cells - 1};
      right_face.neighbor_id = (c < n_cells - 1) ? c + 1 : 1;
      right_face.normal = Normal(0.0, 0.0, 1.0);
      cell.faces.push_back(right_face);

      mesh->ijk_mapping.push_back({c});
      mesh->cells.push_back(cell);
    }

    // define the boundary cells
    mesh->boundary_cell_ids = {0, n_cells - 1};

    // compute the cell and face geometric info
    mesh->compute_geometric_info();

    if (verbose)
      std::cout << "Mesh Details:\n"
                << "\t# of Vertices: " << mesh->vertices.size() << "\n"
                << "\t# of Cells:    " << mesh->cells.size() << "\n"
                << "\t# of Faces:    " << mesh->cells.size() + 1 << "\n";
    return mesh;
  }


  std::shared_ptr<Mesh>
  create_1d_orthomesh(const std::vector<double>& zone_edges,
                      const std::vector<size_t>& zone_subdivisions,
                      const std::vector<unsigned int>& material_ids,
                      const CoordinateSystemType coordinate_system,
                      const bool verbose)
  {
    assert(!zone_edges.empty());
    assert(!zone_subdivisions.empty());
    assert(!material_ids.empty());
    assert(zone_edges.size() == zone_subdivisions.size() + 1);
    assert(zone_subdivisions.size() == material_ids.size());

    std::cout << "Creating a 1D mesh from zones.\n";

    // create the mesh
    auto mesh = std::make_shared<Mesh>(1, coordinate_system);

    // count the number of cells
    auto n_cells = std::accumulate(
        zone_subdivisions.begin(),
        zone_subdivisions.end(),
        (size_t)0
    );

    // initialize the vertices
    mesh->vertices.reserve(n_cells + 1);
    mesh->vertices.emplace_back(0.0, 0.0, zone_edges[0]);

    // define the vertices, loop over each zone, then the cells per zone
    double current_pos = zone_edges[0];
    for (size_t z = 0; z < zone_subdivisions.size(); ++z)
    {
      // define the width of cells in this zone z
      double zone_width = zone_edges[z + 1] - zone_edges[z];
      size_t n_zone_cells = zone_subdivisions[z];
      double cell_width = zone_width / (double) n_zone_cells;

      for (size_t c = 0; c < n_zone_cells; ++c)
      {
        mesh->vertices.emplace_back(0.0, 0.0, current_pos + cell_width);
        current_pos += cell_width;
      }
    }

    // get the type of cell from the coordinate system
    CellType cell_type;
    switch (coordinate_system)
    {
      case CoordinateSystemType::CARTESIAN:
        cell_type = CellType::SLAB;
        break;

      case CoordinateSystemType::CYLINDRICAL:
        cell_type = CellType::ANNULUS;
        break;

      case CoordinateSystemType::SPHERICAL:
        cell_type = CellType::SHELL;
        break;
    }

    // create the cells, loop over zones, then cells per zone
    size_t count = 0;
    for (size_t z = 0; z < zone_subdivisions.size(); ++z)
      for (size_t c = 0; c < zone_subdivisions[z]; ++c)
      {
        Cell cell(cell_type);
        Face left_face, right_face;

        cell.id = count;
        cell.vertex_ids = {count, count + 1};
        cell.material_id = material_ids[z];

        left_face.vertex_ids = {count};
        left_face.has_neighbor = (count > 0);
        left_face.neighbor_id = (count > 0) ? count - 1 : 0;
        left_face.normal = Normal(0.0, 0.0, -1.0);
        cell.faces.push_back(left_face);

        right_face.vertex_ids = {count + 1};
        right_face.has_neighbor = {count < n_cells - 1};
        right_face.neighbor_id = (count < n_cells - 1) ? count + 1 : 1;
        right_face.normal = Normal(0.0, 0.0, 1.0);
        cell.faces.push_back(right_face);

        mesh->ijk_mapping.push_back({c});
        mesh->cells.emplace_back(cell);
        ++count;
      }//for cell

    // define the boundary cells
    mesh->boundary_cell_ids = {0, n_cells - 1};

    // compute the cell and face geometric info
    mesh->compute_geometric_info();

    if (verbose)
      std::cout << "Mesh Details:\n"
                << "\t# of Vertices: " << mesh->vertices.size() << "\n"
                << "\t# of Cells:    " << mesh->cells.size() << "\n"
                << "\t# of Faces:    " << mesh->cells.size() + 1 << "\n";
    return mesh;
  }

}

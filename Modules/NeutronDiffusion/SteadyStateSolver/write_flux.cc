#include "steadystate_solver.h"

#include <fstream>
#include <cstring>

void neutron_diffusion::SteadyStateSolver::write_solution(
    const std::string file_base, const math::Vector& output_flux) const
{

  std::string filename = file_base + ".data";

  std::cout << "Writing solution to file " << filename << "...\n";

  std::ofstream file(filename, std::ofstream::binary |
                               std::ofstream::out |
                               std::ofstream::trunc);

  if (not file.is_open())
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "Failed to open " << filename << ".";
    throw std::runtime_error(err.str());
  }

  // Write the header
  std::string header =
      "PDEs Steady State Neutron Diffusion:\n"
      "Header Size: 500 bytes\n"
      "Structure(type-info):\n"
      "uint64_t n_cells\n"
      "uint64_t n_nodes\n"
      "uint64_t n_groups\n"
      "Each Cell:\n"
      "  uint64_t cell_id\n"
      "  uint64_t material_id\n"
      "  uint64_t n_nodes\n"
      "  Each Node:\n"
      "    double x_position\n"
      "    double y_position\n"
      "    double z_position\n"
      "Each Scalar Flux:\n"
      "  uint64_t cell_id\n"
      "  uint64_t node_number\n"
      "  uint64_t group_number\n"
      "  double   scalar_flux_value\n";

  int header_size = (int)header.length();

  char header_bytes[500];
  memset(header_bytes, '-', 500);
  strncpy(header_bytes, header.c_str(), std::min(header_size,499));
  header_bytes[499]='\0';

  file << header_bytes;

  uint64_t num_nodes = discretization->n_nodes();
  uint64_t num_groups = n_groups;
  uint64_t num_cells = mesh->cells.size();

  file.write((char*)&num_cells, sizeof(uint64_t));
  file.write((char*)&num_nodes, sizeof(uint64_t));
  file.write((char*)&num_groups, sizeof(uint64_t));

  for (const auto& cell : mesh->cells)
  {
    uint64_t id = cell->id;
    file.write((char*)&id, sizeof(uint64_t));

    uint64_t material_id = cell->material_id;
    file.write((char*)&material_id, sizeof(uint64_t));

    uint64_t nodes_per_cell = discretization->nodes_per_cell();
    file.write((char*)&nodes_per_cell, sizeof(uint64_t));

    auto node_loc = cell->centroid;
    file.write((char*)&node_loc.x, sizeof(double));
    file.write((char*)&node_loc.y, sizeof(double));
    file.write((char*)&node_loc.z, sizeof(double));
  }

  for (const auto& cell : mesh->cells)
    for (size_t i = 0; i < discretization->nodes_per_cell(); ++i)
      for (size_t g = 0; g < num_groups; ++g)
      {
        uint64_t cell_id = cell->id;
        uint64_t dof_map = cell_id*num_groups + g;
        double value = output_flux[dof_map];

        file.write((char*)&cell_id, sizeof(uint64_t));
        file.write((char*)&i, sizeof(uint64_t));
        file.write((char*)&g, sizeof(uint64_t));
        file.write((char*)&value, sizeof(double));
      }

  file.close();
}
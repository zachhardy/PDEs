#include "transient_solver.h"

#include <fstream>
#include <cstring>
#include <cassert>

using namespace NeutronDiffusion;


void
TransientSolver::write(const size_t output_index) const
{
  std::string file_base(std::to_string(output_index));
  file_base.insert(0, 4 - file_base.size(), '0');
  std::string file_name = output_directory + "/" + file_base + ".data";

  // Open the file
  if (verbosity > 1)
    std::cout << "Writing snapshot to " << file_name << "\n";

  std::ofstream file(file_name,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header
  int size = 1400;
  std::string header_info =
    "NeutronDiffusion::TransientSolver: Snapshot file\n"
    "Header size: " + std::to_string(size) + " bytes\n";
  header_info +=
    "Structure(type-info):\n"
    "uint64_t      time_step\n"
    "double        time\n"
    "double        power\n"
    "double        peak_power_density\n"
    "double        average_power_density\n"
    "double        peak_fuel_temperature\n"
    "double        average_fuel_temperature\n"
    "unsigned_int  n_data_blocks\n"
    "unsigned int  n_cells\n"
    "unsigned int  n_nodes\n"
    "unsigned int  n_moments\n"
    "unsigned int  n_groups\n"
    "unsigned int  n_precursors\n"
    "unsigned int  max_precursors\n"
    "Each cell:\n"
    "  uint64_t     cell_id\n"
    "  uint64_t     material_id\n"
    "  unsigned int n_nodes\n"
    "  Centroid:\n"
    "    double  centroid_x_position\n"
    "    double  centroid_y_position\n"
    "    double  centroid_z_position\n"
    "  Each node:\n"
    "    double   x_position\n"
    "    double   y_position\n"
    "    double   z_position\n"
    "Scalar Flux Records:\n"
    "  unsigned int  record_type\n"
    "  uint64_t      n_data_blocks\n"
    "  Each Record:\n"
    "    uint64_t      cell_id\n"
    "    unsigned int  node\n"
    "    unsigned int  moment\n"
    "    unsigned int  group\n"
    "    double        value\n"
    "Precursor Records:\n"
    "  unsigned int record_type\n"
    "  uint64_t     n_data_blocks\n"
    "  Each Record:\n"
    "    uint64_t      cell_id\n"
    "    uint64_t      material_id\n"
    "    unsigned int  precursor\n"
    "    double        value\n"
    "Power Density Records:\n"
    "  unsigned int record_type\n"
    "  uint64_t     n_data_blocks\n"
    "  Each Record:\n"
    "    uint64_t cell_id\n"
    "    double   value\n";

  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  const unsigned int n_data_blocks = (use_precursors)? 4 : 3;
  uint64_t n_records;
  unsigned int record_type = 0;

  const uint64_t n_cells = mesh->cells.size();
  const uint64_t n_nodes = discretization->nodes_per_cell();

  const unsigned int n_moments = 1;


  // Write header_info and general information
  file << header_bytes;

  file.write((char*)&output_index, sizeof(uint64_t));
  file.write((char*)&time, sizeof(double));

  file.write((char*)&power, sizeof(double));
  file.write((char*)&peak_power_density, sizeof(double));
  file.write((char*)&average_power_density, sizeof(double));
  file.write((char*)&peak_fuel_temperature, sizeof(double));
  file.write((char*)&average_fuel_temperature, sizeof(double));

  file.write((char*)&n_data_blocks, sizeof(unsigned int));

  file.write((char*)&n_cells, sizeof(uint64_t));
  file.write((char*)&n_nodes, sizeof(uint64_t));

  file.write((char*)&n_moments, sizeof(unsigned int));
  file.write((char*)&n_groups, sizeof(unsigned int));
  file.write((char*)&n_precursors, sizeof(unsigned int));
  file.write((char*)&max_precursors, sizeof(unsigned int));

  // Write discretization information
  for (const auto& cell : mesh->cells)
  {
    const uint64_t cell_id = cell.id;
    const uint64_t material_id = cell.material_id;

    const auto& centroid = cell.centroid;

    const auto& nodes = discretization->nodes(cell);
    const unsigned int npc = nodes.size();

    file.write((char*)&cell_id, sizeof(uint64_t));
    file.write((char*)&material_id, sizeof(uint64_t));
    file.write((char*)&npc, sizeof(unsigned int));

    // Write the centroid
    file.write((char*)&centroid.x, sizeof(double));
    file.write((char*)&centroid.y, sizeof(double));
    file.write((char*)&centroid.z, sizeof(double));

    // Write the nodes
    for (const auto& node : nodes)
    {
      file.write((char*)&node.x, sizeof(double));
      file.write((char*)&node.y, sizeof(double));
      file.write((char*)&node.z, sizeof(double));
    }
  }//for cell

  // Write scalar flux data
  n_records = phi.size();

  file.write((char*)&record_type, sizeof(unsigned int));
  file.write((char*)&n_records, sizeof(uint64_t));

  for (const auto& cell : mesh->cells)
  {
    const uint64_t cell_id = cell.id;
    const unsigned int npc = discretization->nodes(cell).size();

    for (unsigned int i = 0; i < npc; ++i)
      for (unsigned int m = 0; m < n_moments; ++m)
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          const uint64_t dof = cell.id*npc*n_groups + i*n_groups + g;
          assert(dof < phi.size());

          file.write((char*)&cell_id, sizeof(uint64_t));
          file.write((char*)&i, sizeof(unsigned int));
          file.write((char*)&m, sizeof(unsigned int));
          file.write((char*)&g, sizeof(unsigned int));
          file.write((char*)&phi[dof], sizeof(double));
        }
  }//for cell
  ++record_type;

  // Write precursor data
  if (use_precursors)
  {
    n_records = precursors.size();

    file.write((char*)&record_type, sizeof(unsigned int));
    file.write((char*)&n_records, sizeof(uint64_t));

    for (const auto& cell : mesh->cells)
    {
      const uint64_t cell_id = cell.id;
      const uint64_t material_id = cell.material_id;

      const auto& xs = material_xs[matid_to_xs_map[material_id]];
      for (unsigned int j = 0; j < xs->n_precursors; ++j)
      {
        const uint64_t dof = cell.id*max_precursors + j;
        assert(dof < precursors.size());

        file.write((char*)&cell_id, sizeof(uint64_t));
        file.write((char*)&material_id, sizeof(uint64_t));
        file.write((char*)&j, sizeof(unsigned int));
        file.write((char*)&precursors[dof], sizeof(double));
      }
    }//for cell
  }
  ++record_type;

  // Write power density data
  n_records = fission_rate.size();

  file.write((char*)&record_type, sizeof(unsigned int));
  file.write((char*)&n_records, sizeof(uint64_t));

  for (const auto& cell : mesh->cells)
  {
    const uint64_t cell_id = cell.id;
    const double value = energy_per_fission * fission_rate[cell_id];

    file.write((char*)&cell_id, sizeof(uint64_t));
    file.write((char*)&value, sizeof(double ));
  }//for cell
  ++record_type;

  // Write temperature data
  n_records = temperature.size();

  file.write((char*)&record_type, sizeof(unsigned int));
  file.write((char*)&n_records, sizeof(uint64_t));

  for (const auto& cell : mesh->cells)
  {
    const uint64_t cell_id = cell.id;
    const double value = temperature[cell_id];

    file.write((char*)&cell_id, sizeof(uint64_t));
    file.write((char*)&value, sizeof(double));
  }//for cell
  ++record_type;

  file.close();
  std::cout << "HERE\n";
}

#include "steadystate_solver.h"

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>


using namespace NeutronDiffusion;


void
SteadyStateSolver::
write(const std::string output_directory,
      const std::string file_prefix) const
{
  if (not std::filesystem::is_directory(output_directory))
    std::filesystem::create_directory(output_directory);
  assert(std::filesystem::is_directory(output_directory));

  std::string filepath = output_directory + "/" + file_prefix;
  if (filepath.find(".") != std::string::npos)
    filepath = file_prefix.substr(0, file_prefix.rfind("."));
  filepath += ".data";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header_info
  int size = 1000;
  std::string header_info =
    "NeutronDiffusion::SteadyStateSolver: Snapshot File\n"
    "Header size: " + std::to_string(size) + " bytes\n";
  header_info +=
    "Structure(type-info):\n"
    "unsigned int  n_data_blocks\n"
    "uint64_t      n_cells\n"
    "uint64_t      n_nodes\n"
    "unsigned int  n_moments\n"
    "unsigned int  n_groups\n"
    "unsigned int  n_precursors\n"
    "unsigned int  max_precursors\n"
    "Each Cell:\n"
    "  uint64_t      cell_id\n"
    "  uint64_t      material_id\n"
    "  unsigned int  n_nodes\n"
    "  Centroid:\n"
    "    double  centroid_x_position\n"
    "    double  centroid_y_position\n"
    "    double  centroid_z_position\n"
    "  Each node:\n"
    "    double  x_position\n"
    "    double  y_position\n"
    "    double  z_position\n"
    "Scalar Flux Records:\n"
    "  unsigned int  record_type\n"
    "  uint64_t      n_records\n"
    "  Each Record:\n"
    "    uint64_t      cell_id\n"
    "    unsigned int  node\n"
    "    unsigned int  moment\n"
    "    unsigned int  group\n"
    "    double        value\n"
    "Precursor Records:\n"
    "  unsigned int  record_type\n"
    "  uint64_t      n_records\n"
    "  Each Record:\n"
    "    uint64_t      cell_id\n"
    "    uint64_t      material_id\n"
    "    unsigned int  precursor\n"
    "    double        value\n";

  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  uint64_t n_records;
  unsigned int record_type = 0;
  unsigned int n_data_blocks = 2;

  const uint64_t n_cells = mesh->cells.size();
  const uint64_t n_nodes = discretization->n_nodes();

  const unsigned int n_moments = 1;

  // Write header_info and general information
  file << header_bytes;
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
    double x = centroid.x(), y = centroid.y(), z = centroid.z();
    file.write((char*)&x, sizeof(double));
    file.write((char*)&y, sizeof(double));
    file.write((char*)&z, sizeof(double));

    // Write the nodes
    for (const auto& node : nodes)
    {
      x = node.x(), y = node.y(), z = node.z();
      file.write((char*)&x, sizeof(double));
      file.write((char*)&y, sizeof(double));
      file.write((char*)&z, sizeof(double));
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
      for (unsigned int j = 0; j < max_precursors; ++j)
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
  file.close();
}


void
SteadyStateSolver::
write_flux_moments(const std::string directory,
                   const std::string file_prefix) const
{
  if (not std::filesystem::is_directory(directory))
    std::filesystem::create_directory(directory);
  assert(std::filesystem::is_directory(directory));

  std::string filepath = directory + "/" + file_prefix;
  if (filepath.find(".") != std::string::npos)
    filepath = file_prefix.substr(0, file_prefix.rfind("."));
  filepath += ".data";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header
  int size = 200;
  std::stringstream ss;
  ss
      << "Multi-Group Flux Moment Snapshot\n"
         "Header size: " << size  << " bytes\n"
      << "Structure(type-info)\n"
         "uint64_t      n_nodes\n"
         "unsigned int  n_moments\n"
         "unsigned int  n_groups\n"
         "Each Record:\n"
         "  uint64_t      cell\n"
         "  unsigned int  node\n"
         "  unsigned int  moment\n"
         "  unsigned int  group\n"
         "  double        value\n";

  std::string header_info = ss.str();
  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  // Write the data
  file << header_bytes;

  const uint64_t n_nodes = discretization->n_nodes();
  file.write((char*)&n_nodes, sizeof(uint64_t));

  const unsigned int n_moments = 1;
  file.write((char*)&n_moments, sizeof(unsigned int));

  file.write((char*)&n_groups, sizeof(unsigned int));

  uint64_t uk_map = 0;
  for (const auto& cell : mesh->cells)
  {
    const auto nodes = discretization->nodes(cell);
    for (unsigned int i = 0; i < nodes.size(); ++i)
      for (unsigned int ell = 0; ell < 1; ++ell)
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          file.write((char*) &cell.id, sizeof(uint64_t));
          file.write((char*)&i, sizeof(unsigned int));
          file.write((char*)&ell, sizeof(unsigned int));
          file.write((char*)&g, sizeof(unsigned int));
          file.write((char*)&ell, sizeof(unsigned int));
          file.write((char*)&phi[uk_map++], sizeof(double));
        }//for g
  }//for cell
}

void
SteadyStateSolver::
write_precursors(const std::string directory,
                 const std::string file_prefix) const
{
  if (not std::filesystem::is_directory(directory))
    std::filesystem::create_directory(directory);
  assert(std::filesystem::is_directory(directory));

  std::string filepath = directory + "/" + file_prefix;
  if (filepath.find(".") != std::string::npos)
    filepath = file_prefix.substr(0, file_prefix.rfind("."));
  filepath += ".data";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header
  int size = 200;
  std::stringstream ss;
  ss
      << "Delayed Neutron Precursor Snapshot\n"
         "Header size: " << size  << " bytes\n"
      << "Structure(type-info)\n"
         "uint64_t      n_cells\n"
         "unsigned int  max_precursors\n"
         "Each Record:\n"
         "  uint64_t      cell\n"
         "  unsigned int  material_id\n"
         "  unsigned int  precursor\n"
         "  double        value\n";

  std::string header_info = ss.str();
  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  // Write the data
  file << header_bytes;

  const uint64_t n_cells = mesh->cells.size();
  file.write((char*)&n_cells, sizeof(uint64_t));

  file.write((char*)&max_precursors, sizeof(unsigned int));

  uint64_t uk_map = 0;
  for (const auto& cell : mesh->cells)
    for (unsigned int j = 0; j < max_precursors; ++j)
    {
      file.write((char*) &cell.id, sizeof(uint64_t));
      file.write((char*)&cell.material_id, sizeof(unsigned int));
      file.write((char*)&j, sizeof(unsigned int));
      file.write((char*)&precursors[uk_map++], sizeof(double));
    }//for j
}


void
SteadyStateSolver::
write_fission_rate(const std::string directory,
                   const std::string file_prefix) const
{
  if (not std::filesystem::is_directory(directory))
    std::filesystem::create_directory(directory);
  assert(std::filesystem::is_directory(directory));

  std::string filepath = directory + "/" + file_prefix;
  if (filepath.find(".") != std::string::npos)
    filepath = file_prefix.substr(0, file_prefix.rfind("."));
  filepath += ".data";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header
  int size = 200;
  std::stringstream ss;
  ss
      << "Fission Rate Snapshot\n"
         "Header size: " << size  << " bytes\n"
      << "Structure(type-info)\n"
         "uint64_t      n_cells\n"
         "Each Record:\n"
         "  uint64_t      cell\n"
         "  double        value\n";

  std::string header_info = ss.str();
  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  // Write the data
  file << header_bytes;

  const uint64_t n_cells = mesh->cells.size();
  file.write((char*)&n_cells, sizeof(uint64_t));

  for (const auto& cell : mesh->cells)
  {
    // Get cross-section data
    const auto xs_id = matid_to_xs_map[cell.material_id];
    const auto xs = material_xs[xs_id];

    // Compute fission rate cell-wise
    double f = 0.0;
    const double* sig_f = xs->sigma_f.data();
    for (unsigned int g = 0; g < n_groups; ++g)
      f += *sig_f++ * phi[cell.id * n_groups + g];

    // Write the fission rate
    file.write((char*) &cell.id, sizeof(uint64_t));
    file.write((char*)&f, sizeof(double));
  }//for cell
}

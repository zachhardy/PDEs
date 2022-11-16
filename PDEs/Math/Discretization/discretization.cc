#include "discretization.h"

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>


using namespace PDEs;
using namespace Math;


void
Math::Discretization::
write(const std::string directory,
      const std::string file_prefix) const
{
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

  int size = 400;
  std::string header_info =
      "PDEs Geometry File\n"
      "Header size: " + std::to_string(size) + " bytes\n";
  header_info +=
      "Structure(type-info):\n"
      "uint64_t   n_cells\n"
      "uint64_t   n_nodes\n"
      "Each Cell:\n"
      "   uint64_t      cell_id\n"
      "   unsigned int  material_id\n"
      "   unsigned int  n_nodes\n"
      "   Centroid:\n"
      "      double  x_position\n"
      "      double  y_position\n"
      "      double  z_position\n"
      "   Each Node:\n"
      "      double  x_position\n"
      "      double  y_position\n"
      "      double  z_position\n";

  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  file << header_bytes;

  const uint64_t n_cells = mesh->cells.size();
  file.write((char*)&n_cells, sizeof(uint64_t));

  const uint64_t n_nodes = this->n_nodes();
  file.write((char*)&n_nodes, sizeof(uint64_t));

  for (const auto& cell : mesh->cells)
  {
    file.write((char*)&cell.id, sizeof(uint64_t));
    file.write((char*)&cell.material_id, sizeof(unsigned int));

    const unsigned int npc = this->nodes_per_cell();
    file.write((char*)&npc, sizeof(unsigned int));

    file.write((char*)&cell.centroid.x(), sizeof(double));
    file.write((char*)&cell.centroid.y(), sizeof(double));
    file.write((char*)&cell.centroid.z(), sizeof(double));

    const auto& nodes = this->nodes(cell);
    for (const auto node : nodes)
    {
      file.write((char*)&node.x(), sizeof(double));
      file.write((char*)&node.y(), sizeof(double));
      file.write((char*)&node.z(), sizeof(double));
    }
  }
  file.close();
}

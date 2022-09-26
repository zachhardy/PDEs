#include "transient_solver.h"

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>

using namespace NeutronDiffusion;


void
TransientSolver::write(const unsigned int output_index) const
{
  std::string directory(std::to_string(output_index));
  directory.insert(0, 4 - directory.size(), '0');
  directory = output_directory + "/" + directory;
  if (!std::filesystem::is_directory(directory))
    std::filesystem::create_directories(directory);

  std::string filepath = output_directory + "/summary.txt";

  // Open the file
  std::ofstream file(filepath, std::ofstream::out | std::ofstream::app);
  assert(file.is_open());

  file.setf(std::ios::left);
  file.setf(std::ios::scientific, std::ios::floatfield);

  file << "  ";
  file << std::setw(18) << std::setprecision(10) << time
       << std::setw(18) << std::setprecision(10) << power
       << std::setw(18) << std::setprecision(10) << peak_power_density
       << std::setw(18) << std::setprecision(10) << average_power_density
       << std::setw(18) << std::setprecision(10) << peak_fuel_temperature
       << std::setw(18) << std::setprecision(10) << average_fuel_temperature
       << std::endl;

  file.close();

  write_flux_moments(directory);
  write_precursors(directory);
  write_fission_rate(directory);
  write_temperature(directory);
}


void
TransientSolver::
write_temperature(const std::string directory,
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

  // Write the header
  int size = 200;
  std::stringstream ss;
  ss
      << "Temperature Snapshot\n"
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
    file.write((char*) &cell.id, sizeof(uint64_t));
    file.write((char*)&temperature[cell.id], sizeof(double));
  }//for cell
}

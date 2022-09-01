#include "steadystate_solver.h"

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>


using namespace InfiniteMedium;


void
SteadyStateSolver::
write_angular_flux(const std::string directory,
                   const std::string file_prefix) const

{
  if (not std::filesystem::is_directory(directory))
    std::filesystem::create_directory(directory);
  assert(std::filesystem::is_directory(directory));

  std::string filepath = directory + "/" + file_prefix;
  if (filepath.find(".") != std::string::npos)
    filepath = file_prefix.substr(0, file_prefix.rfind("."));
  filepath = filepath + ".data";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header_info
  int size = 250;
  std::string header_info =
      "InfiniteMedium::SteadyStateSolver: Angular Flux File\n"
      "Header size: " + std::to_string(size) + " bytes\n";
  header_info +=
      "Structure(type-info):\n"
      "unsigned int   n_angles\n"
      "unsigned int   n_groups\n"
      "Angular Flux Records:\n"
      "  Each Record:\n"
      "    double value\n";

  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  file << header_bytes;
  file.write((char*)&n_angles, sizeof(unsigned int));
  file.write((char*)&n_groups, sizeof(unsigned int));

  for (const auto& value : psi)
    file.write((char*)&value, sizeof(double));

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
  filepath = filepath + ".data";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::binary |
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  // Write the header_info
  int size = 250;
  std::string header_info =
      "InfiniteMedium::SteadyStateSolver: Flux Moment File\n"
      "Header size: " + std::to_string(size) + " bytes\n";
  header_info +=
      "Structure(type-info):\n"
      "unsigned int   n_moments\n"
      "unsigned int   n_groups\n"
      "Flux Moment Records:\n"
      "  Each Record:\n"
      "    double value\n";

  int header_size = (int)header_info.length();

  char header_bytes[size];
  memset(header_bytes, '-', size);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, size - 1));
  header_bytes[size - 1] = '\0';

  file << header_bytes;
  file.write((char*)&n_moments, sizeof(unsigned int));
  file.write((char*)&n_groups, sizeof(unsigned int));

  for (const auto& value : phi)
    file.write((char*)&value, sizeof(double));

  file.close();
}

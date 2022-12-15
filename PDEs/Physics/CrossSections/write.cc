#include "cross_sections.h"

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>


using namespace PDEs;
using namespace Physics;


void
CrossSections::
write_group_structure(const std::string directory,
                      const std::string file_prefix) const
{
  if (not std::filesystem::is_directory(directory))
    std::filesystem::create_directory(directory);
  assert(std::filesystem::is_directory(directory));

  std::string filepath = directory + "/" + file_prefix;
  if (filepath.find(".") != std::string::npos)
    filepath = file_prefix.substr(0, file_prefix.rfind("."));
  filepath = filepath + ".txt";

  // Open the file
  std::ofstream file(filepath,
                     std::ofstream::out |
                     std::ofstream::trunc);
  assert(file.is_open());

  file << "########################################\n"
       << "# " << n_groups << " Group Structure (MeV) \n"
       << "########################################\n";
  for (unsigned int g = 0; g < n_groups + 1; ++g)
    file << g << "  " << e_bounds[g] << "\n";

  file.close();
}

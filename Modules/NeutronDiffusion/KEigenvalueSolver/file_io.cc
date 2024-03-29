#include "keigenvalue_solver.h"

#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>


using namespace NeutronDiffusion;


void
KEigenvalueSolver::
write(const std::string directory,
      const std::string file_prefix) const
{
  SteadyStateSolver::write(directory, file_prefix);

  std::string filepath = directory + "/k_eff.txt";
  std::ofstream file(filepath, std::ofstream::out | std::ofstream::trunc);
  assert(file.is_open());

  file << std::fixed << std::setprecision(12) << k_eff << "\0";
  file.close();

}
#include "transient_solver.h"

#include <filesystem>

using namespace InfiniteMedium;


void TransientSolver::initialize()
{
  SteadyStateSolver::initialize();

  // Set the initial condition
  for (const auto& ic : initial_conditions)
    for (unsigned int n = 0; n < n_angles; ++n)
    {
      const double x = ic.second;
      psi[n * n_groups + ic.first] = x;

      for (unsigned int ell = 0; ell < n_moments; ++ell)
        phi[ell * n_groups + ic.first] += discrete_to_moment[ell][n] * x;
    }
  psi_old = psi;
  phi_old = phi;



  // Initialize output directory
  if (write_outputs)
  {
    output_directory = std::filesystem::absolute(output_directory);

    std::cout << "Setting up output directories at "
              << output_directory << std::endl;

    if (not std::filesystem::is_directory(output_directory))
      std::filesystem::create_directory(output_directory);

    using DirectoryIterator = std::filesystem::directory_iterator;
    for (const auto& entry: DirectoryIterator(output_directory))
      std::filesystem::remove_all(entry.path());

    xs->write_group_structure(output_directory);
  }
}

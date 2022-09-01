#include "transient_solver.h"

#include <iomanip>
#include <fstream>
#include <cstring>
#include <cassert>
#include <filesystem>


using namespace InfiniteMedium;


void
TransientSolver::execute()
{
  std::cout
    << "\n**********************************************************\n"
    << "Executing the multi-group infinite medium transient solver"
    << "\n**********************************************************\n";

  const double eps = 1.0e-10;

  unsigned int count = 0;
  if (write_outputs)
    write_snapshot(count++);

  unsigned int step = 0;
  while (time < t_end - eps)
  {
    // Solve the time step
    execute_time_step();
    double balance = check_balance();

    if (write_outputs)
      write_snapshot(count++);

    // Post-process the time step
    time += dt;
    ++step;

    phi_old = phi_ell = phi;
    psi_old = psi;


    std::cout
      << "\n***** Time Step " << step << " *****\n"
      << "Simulation Time:        " << time << " s\n"
      << "Time Step Size :        " << dt << " s\n"
      << "Balance:                " << balance << "\n";
  }
}


void
TransientSolver::execute_time_step()
{
  unsigned int nit;
  double change;
  bool converged;

  // Start iterations
  phi_ell = phi_old;
  for (nit = 0; nit < max_inner_iterations; ++nit)
  {
    // Compute the RHS and solve
    q_moments = 0.0;
    set_source(APPLY_MATERIAL_SOURCE |
               APPLY_SCATTER_SOURCE |
               APPLY_FISSION_SOURCE);
    transient_sweep();

    change = l1_norm(phi - phi_ell);
    converged = change < inner_tolerance;
    phi_ell = phi;

    if (verbosity > 1)
      std::cout
        << std::left << (nit == 0? "\n" : "")
        << "inner::"
        << "Iteration  " << std::setw(5) << nit
        << "Change  " << std::setw(8) << change
        << (converged? "  CONVERGED" : "  ")
        << std::endl;

    if (converged) break;
  }//for nit
}


void
TransientSolver::transient_sweep()
{
  phi = 0.0;
  for (unsigned int n = 0; n < n_angles; ++n)
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      // Solve for the angular flux
      double x = xs->inv_velocity[g] / dt * psi_old[n * n_groups + g];
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        x += moment_to_discrete[ell][n] * q_moments[ell * n_groups + g];
      x /= (xs->inv_velocity[g] / dt + xs->sigma_t[g]);

      // Store the angular flux
      psi[n*n_groups + g] = x;

      // Accumulate flux moment
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        phi[ell*n_groups + g] += discrete_to_moment[ell][n] * x;
    }//for g
}


double
TransientSolver::balance_check()
{
  double balance = 0.0;
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    balance += src->values[g];
    balance += phi_old[g] / (xs->inv_velocity[g] * dt);
    balance -= phi[g] / (xs->inv_velocity[g] * dt);
    balance -= xs->sigma_a[g] * phi[g];
  }
  return balance;
}


void
TransientSolver::write_snapshot(const unsigned int index) const
{
  std::string directory(std::to_string(index));
  directory.insert(0, 4 - directory.size(), '0');
  directory = output_directory + "/" + directory;
  directory = std::filesystem::absolute(directory);

  if (not std::filesystem::is_directory(directory))
    std::filesystem::create_directory(directory);
  assert(std::filesystem::is_directory(directory));

  //================================================== Write snapshot summary

  std::string filepath = directory + "/" + "macro.txt";
  std::ofstream file(filepath, std::ofstream::out | std::ofstream::trunc);
  assert(file.is_open());

  file << "########################################\n"
       << "# Snapshot " << index << " Summary\n"
       << "########################################\n"
       << "Output= " << std::setprecision(5) << index << "\n"
       << "Time= " << std::setprecision(12) << time << "\n";

  file.close();

  //================================================== Write simulation data
  write_flux_moments(directory);
  write_angular_flux(directory);
}

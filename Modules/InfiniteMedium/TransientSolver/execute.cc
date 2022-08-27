#include "transient_solver.h"

#include <iomanip>


using namespace InfiniteMedium;


void
TransientSolver::execute()
{
  std::cout
    << "\n**********************************************************\n"
    << "Executing the multi-group infinite medium transient solver"
    << "\n**********************************************************\n";

  const double eps = 1.0e-10;

  unsigned int step = 0;
  while (time < t_end - eps)
  {
    // Solve the time step
    execute_time_step();

    // Post-process the time step
    time += dt;
    ++step;

    phi_old = phi_ell = phi;
    psi_old = psi;

    std::cout
      << "\n***** Time Step " << step << " *****\n"
      << "Simulation Time:        " << time << " s\n"
      << "Time Step Size :        " << dt << " s\n"
      << "Magnitude:              " << psi.l2_norm() << "\n";
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

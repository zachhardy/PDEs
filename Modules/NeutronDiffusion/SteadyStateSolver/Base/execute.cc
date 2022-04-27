#include "steadystate_solver.h"

#include <fstream>

/// Run the steady state multigroup diffusion simulation.
void neutron_diffusion::SteadyStateSolver::execute()
{
  std::cout << "\nExecuting solver...\n";

  // Initialize matrices
  for (auto& gs : groupsets)
  { assemble_matrix(gs); gs.linear_solver->setup(); }

  auto phi_ell = phi;
  double diff = 1.0;
  bool converged = false;

  //======================================== Loop over groupsets
  for (auto& groupset : groupsets)
    solve_groupset(groupset);

  std::cout << "\nDone executing solver.\n";
}

//######################################################################

/**
 * \brief Solve
 * \param groupset
 */
void neutron_diffusion::SteadyStateSolver::
solve_groupset(Groupset& groupset)
{
  std::cout << "***** Solving Groupset " << groupset.id << "\n";

  math::Vector phi_gs(groupset.rhs.size(), 0.0);
  auto phi_gs_ell = phi_gs;
  double diff = 1.0;
  bool converged = false;

  //======================================== Start iterations
  for (size_t nit = 0; nit < groupset.max_iterations; ++nit)
  {
    // Clear and reset the RHS
    groupset.rhs *= 0.0;
    set_source(groupset, groupset.rhs);

    // Solve the system
    phi_gs = groupset.linear_solver->solve(groupset.rhs);

    // Convergence check, finalize iteration
    diff = math::l2_norm(phi_gs - phi_gs_ell);
    scoped_transfer(groupset, phi_gs, phi);
    phi_gs_ell = phi_gs;

    if (diff < groupset.tolerance)
      converged = true;

    // Print iteration information
    std::stringstream iter_info;
    iter_info << "Iteration: " << std::setw(3) << nit << " "
              << "Difference: " << diff;
    if (converged) iter_info << " CONVERGED";
    std::cout << iter_info.str() << "\n";

    if (diff < groupset.tolerance)
      break;
  }
}

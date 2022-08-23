#include "steadystate_solver.h"

#include <iomanip>
#include <fstream>


void
NeutronDiffusion::SteadyStateSolver::execute()
{
  std::cout
      << "\n************************************************\n"
      <<   "Executing the multi-group diffusion steady-state solver"
      << "\n************************************************\n";

  // Initialize matrix and solve
  if (algorithm == Algorithm::DIRECT)
  {
    assemble_matrix(ASSEMBLE_SCATTER | ASSEMBLE_FISSION);
    linear_solver->set_matrix(A);

    set_source(APPLY_MATERIAL_SOURCE | APPLY_BOUNDARY_SOURCE);
    linear_solver->solve(phi, b);
  }
  else
  {
    assemble_matrix();
    linear_solver->set_matrix(A);

    iterative_solve(APPLY_MATERIAL_SOURCE | APPLY_BOUNDARY_SOURCE |
                    APPLY_SCATTER_SOURCE | APPLY_FISSION_SOURCE);
  }

  // Compute precursors
  if (use_precursors)
    compute_precursors();
}

//######################################################################

void
NeutronDiffusion::SteadyStateSolver::
iterative_solve(SourceFlags source_flags)
{
  const auto b_init = b;

  unsigned nit;
  double change;
  bool converged = false;

  // Start iterations
  auto x = phi;
  for (nit = 0; nit < max_inner_iterations; ++nit)
  {
    // Compute the RHS and solve
    b = b_init;
    set_source(source_flags);
    linear_solver->solve(phi, b);

    // Convergence check, finalize iteration
    change = l1_norm(phi - x);
    converged = change < inner_tolerance;

    // Print iteration information
    if (verbosity > 1)
    {
      std::stringstream iter_info;
      iter_info
        << std::left << "SourceIteration::"
        << "Step  " << std::setw(3) << nit << "   "
        << "Value  " << change;
      if (converged) iter_info << "  CONVERGED";
      std::cout << iter_info.str() << std::endl;
    }

    if (converged) break;
  }
}

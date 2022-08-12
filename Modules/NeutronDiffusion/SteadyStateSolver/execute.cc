#include "steadystate_solver.h"

#include "LinearSolvers/Direct/sparse_lu.h"
#include "LinearSolvers/Direct/sparse_cholesky.h"
#include "LinearSolvers/Iterative/jacobi.h"
#include "LinearSolvers/Iterative/gauss_seidel.h"
#include "LinearSolvers/Iterative/sor.h"
#include "LinearSolvers/Iterative/ssor.h"
#include "LinearSolvers/Iterative/cg.h"

#include "macros.h"

#include <iomanip>
#include <fstream>


using namespace Math;
using namespace Math::LinearSolver;


void
NeutronDiffusion::SteadyStateSolver::execute()
{
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

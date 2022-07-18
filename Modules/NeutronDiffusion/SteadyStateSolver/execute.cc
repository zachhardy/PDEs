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
  // Initialize matrices
  for (auto& groupset : groupsets)
  {
    if (solution_technique == SolutionTechnique::GROUPSET_WISE)
      assemble_matrix(groupset);
    else
      assemble_matrix(groupset, ASSEMBLE_SCATTER | ASSEMBLE_FISSION);

    linear_solver->set_matrix(groupset.A);
  }

  // Solve
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
    solve_full_system(APPLY_MATERIAL_SOURCE);
  else
    for (auto& groupset : groupsets)
      solve_groupset(groupset,
                     APPLY_MATERIAL_SOURCE |
                     APPLY_WGS_SCATTER_SOURCE | APPLY_AGS_SCATTER_SOURCE |
                     APPLY_WGS_FISSION_SOURCE | APPLY_AGS_FISSION_SOURCE);

  // Compute precursors
  if (use_precursors)
    compute_precursors();
}

//######################################################################

void
NeutronDiffusion::SteadyStateSolver::
solve_groupset(Groupset& groupset, SourceFlags source_flags)
{
  if (verbosity > 1)
    std::cout << "\n***** Solving Groupset "
              << groupset.id << std::endl << std::endl;

  Vector& b = groupset.b;
  const Vector b_init = b;

  size_t nit;
  double change;
  bool converged = false;

  // Start iterations
  for (nit = 0; nit < groupset.max_iterations; ++nit)
  {
    // Compute the RHS and solve
    b = b_init;
    set_source(groupset, source_flags);
    auto x = linear_solver->solve(b);

    // Convergence check, finalize iteration
    scoped_transfer(groupset, x, phi);
    change = compute_change(groupset);
    scoped_copy(groupset, phi, phi_ell);

    if (change < groupset.tolerance)
      converged = true;

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


//###########################################################################


void
NeutronDiffusion::SteadyStateSolver::
solve_full_system(SourceFlags source_flags)
{
  if (verbosity > 1)
    std::cout << "\n***** Solving Full System\n";

  set_source(groupsets.front(), source_flags);
  phi = linear_solver->solve(groupsets.front().b);
}

#include "steadystate_solver.h"

#include "LinearSolvers/Direct/sparse_lu.h"

#include "macros.h"

#include <iomanip>
#include <fstream>

using namespace pdes;
using namespace pdes::Math;
using namespace pdes::Math::LinearSolver;


void
NeutronDiffusion::SteadyStateSolver::execute()
{
  //======================================== Initialize matrices
  for (auto& gs : groupsets)
   assemble_matrix(gs);

  //======================================== Solve
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
    solve_full_system(APPLY_MATERIAL_SOURCE);
  else
    for (auto& groupset : groupsets)
      solve_groupset(groupset,
                     APPLY_MATERIAL_SOURCE |
                     APPLY_WGS_SCATTER_SOURCE | APPLY_AGS_SCATTER_SOURCE |
                     APPLY_WGS_FISSION_SOURCE | APPLY_AGS_FISSION_SOURCE);
}

//######################################################################

void
NeutronDiffusion::SteadyStateSolver::
solve_groupset(Groupset& groupset, SourceFlags source_flags)
{
  std::cout << "\n***** Solving Groupset " << groupset.id << "\n\n";

  SparseMatrix& A = groupset.matrix;
  Vector& b = groupset.rhs;

  auto solver = initialize_linear_solver(groupset);

  //======================================== Start iterations
  size_t nit; double diff; bool converged = false;
  for (nit = 0; nit < groupset.max_iterations; ++nit)
  {
    // Compute the RHS and solve
    b = 0.0;
    set_source(groupset, b, source_flags);
    auto x = solver->solve(b);

    // Convergence check, finalize iteration
    scoped_transfer(groupset, x, phi);
    diff = compute_change(groupset);
    scoped_copy(groupset, phi, phi_ell);

    if (diff < groupset.tolerance)
      converged = true;

    // Print iteration information
    std::stringstream iter_info;
    iter_info << "Iteration: " << std::setw(3) << nit << " "
              << "Change: " << diff;
    if (converged) iter_info << " CONVERGED\n";
    std::cout << iter_info.str() << "\n";

    if (converged) break;
  }
}

//###########################################################################

void
NeutronDiffusion::SteadyStateSolver::
solve_full_system(SourceFlags source_flags)
{
  std::cout << "\n***** Solving Full System\n";

  SparseMatrix& A = groupsets.front().matrix;
  Vector& b = groupsets.front().rhs;

  SparseLU lu(A);

  b = 0.0;
  set_source(groupsets.front(), b, source_flags);
  phi = lu.solve(b);
}

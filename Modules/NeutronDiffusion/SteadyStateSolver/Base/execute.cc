#include "steadystate_solver.h"

#include "LinearSolvers/Direct/lu.h"
#include "LinearSolvers/Direct/cholesky.h"
#include "LinearSolvers/Direct/sparse_lu.h"
#include "LinearSolvers/Direct/sparse_cholesky.h"
#include "LinearSolvers/Iterative/jacobi.h"
#include "LinearSolvers/Iterative/gauss_seidel.h"

#include <iomanip>
#include <fstream>

using namespace pdes;
using namespace Math;


void
NeutronDiffusion::SteadyStateSolver:: execute()
{
  std::cout << "Executing solver...\n";

  // Initialize matrices
  for (auto& gs : groupsets)
   assemble_matrix(gs);

  SourceFlags source_flags = APPLY_MATERIAL_SOURCE;
  if (solution_technique == SolutionTechnique::GROUPSET_WISE)
    source_flags = source_flags | APPLY_WGS_SCATTER_SOURCE |
                   APPLY_AGS_SCATTER_SOURCE | APPLY_WGS_FISSION_SOURCE |
                   APPLY_AGS_FISSION_SOURCE;

  //======================================== Solve
  if (solution_technique == SolutionTechnique::GROUPSET_WISE)
    for (auto& groupset : groupsets)
      solve_groupset(groupset, source_flags);
  else
    solve_full_system(source_flags);

  std::cout << "\nDone executing solver.\n";
}

//######################################################################


void
NeutronDiffusion::SteadyStateSolver::
solve_groupset(Groupset& groupset, SourceFlags source_flags)
{
  std::cout << "\n***** Solving Groupset " << groupset.id << "\n\n";

  double change;
  size_t k;
  bool converged = false;

  SparseMatrix& A = groupset.matrix;
  Vector& b = groupset.rhs;

  GaussSeidelSolver solver;
  Vector x(b.size());

  //======================================== Start iterations
  for (k = 0; k < groupset.max_iterations; ++k)
  {
    // Compute the RHS and solve
    b = 0.0;
    set_source(groupset, b, source_flags);
    solver.solve(A, b, x);

    // Convergence check, finalize iteration
    scoped_transfer(groupset, x, phi);
    change = compute_change(groupset);
    scoped_copy(groupset, phi, phi_ell);

    if (change < groupset.tolerance)
      converged = true;

    // Print iteration information
    std::stringstream iter_info;
    iter_info << "Iteration: " << std::setw(3) << k << " "
              << "Change: " << change;
    if (converged) iter_info << " CONVERGED";
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
  lu.factorize();

  b = 0.0;
  set_source(groupsets.front(), b, source_flags);
  lu.solve(b, phi);
}

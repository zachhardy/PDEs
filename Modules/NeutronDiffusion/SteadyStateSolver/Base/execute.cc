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
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
  {
    solve_full_system(source_flags);
    return;
  }

  for (auto& groupset : groupsets)
    solve_groupset(groupset, source_flags);


  std::cout << "\nDone executing solver.\n";
}


//######################################################################


void
NeutronDiffusion::SteadyStateSolver::
solve_groupset(Groupset& groupset, SourceFlags source_flags)
{
  std::cout << "\n***** Solving Groupset " << groupset.id << "\n\n";

  SparseMatrix& A = groupset.matrix;
  Vector& b = groupset.rhs;

  Vector init_b = b;
  set_source(groupset, init_b, APPLY_MATERIAL_SOURCE |
                               APPLY_AGS_SCATTER_SOURCE |
                               APPLY_AGS_FISSION_SOURCE);

  GaussSeidelSolver solver;
  Vector x(b.size(), 0.0);

  size_t nit;
  double diff;
  bool converged = false;

  //======================================== Start iterations
  for (nit = 0; nit < groupset.max_iterations; ++nit)
  {
    // Compute the RHS and solve
    b = init_b;
    set_source(groupset, b,
               APPLY_WGS_SCATTER_SOURCE |
               APPLY_WGS_SCATTER_SOURCE);
    solver.solve(A, b, x);

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

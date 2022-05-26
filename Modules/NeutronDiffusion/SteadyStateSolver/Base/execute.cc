#include "steadystate_solver.h"

#include "LinearSolvers/Direct/lu.h"
#include "LinearSolvers/Direct/cholesky.h"
#include "LinearSolvers/Direct/sparse_lu.h"
#include "LinearSolvers/Direct/sparse_cholesky.h"
#include "LinearSolvers/Iterative/jacobi.h"
#include "LinearSolvers/Iterative/gauss_seidel.h"

#include "macros.h"

#include <iomanip>
#include <fstream>

using namespace pdes;
using namespace Math;


void
NeutronDiffusion::SteadyStateSolver::execute()
{
  std::cout << "Executing solver...\n";

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

  std::shared_ptr<LinearSolverBase> solver;
  switch (linear_solver_type)
  {
    case LinearSolverType::LU:
      solver = std::make_shared<SparseLU>(A);
    case LinearSolverType::CHOLESKY:
      solver = std::make_shared<SparseCholesky>(A);
    case LinearSolverType::JACOBI:
      solver = std::make_shared<JacobiSolver>(A);
    case LinearSolverType::GAUSS_SEIDEL:
      solver = std::make_shared<GaussSeidelSolver>(A);
    default:
      Assert(true, "Linear solver not implemented.");
  }
  /* std::map<std::string, Varying> opts;
   * set_relaxation_factor()
   * solver.set_additional_options(AdditionalOptions) */

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

  b = 0.0;
  set_source(groupsets.front(), b, source_flags);
  phi = lu.solve(b);
}

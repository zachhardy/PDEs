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
  //======================================== Initialize matrices
  for (auto& groupset : groupsets)
  {
    if (solution_technique == SolutionTechnique::GROUPSET_WISE)
      assemble_matrix(groupset);
    else
      assemble_matrix(groupset, ASSEMBLE_SCATTER | ASSEMBLE_FISSION);
  }

  //======================================== Solve
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
    solve_full_system(APPLY_MATERIAL_SOURCE);
  else
    for (auto& groupset : groupsets)
      solve_groupset(groupset,
                     APPLY_MATERIAL_SOURCE |
                     APPLY_WGS_SCATTER_SOURCE | APPLY_AGS_SCATTER_SOURCE |
                     APPLY_WGS_FISSION_SOURCE | APPLY_AGS_FISSION_SOURCE);

  //======================================== Compute precursors
  if (use_precursors)
    compute_precursors();
}

//######################################################################

void
NeutronDiffusion::SteadyStateSolver::
solve_groupset(Groupset& groupset, SourceFlags source_flags)
{
  if (verbose)
    std::cout << "\n***** Solving Groupset " << groupset.id << "\n\n";

  SparseMatrix& A = groupset.A;
  Vector& b = groupset.b;
  const Vector b_init = b;

  size_t nit;
  double change;
  bool converged = false;

  //======================================== Start iterations
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
    if (verbose)
    {
      std::stringstream iter_info;
      iter_info << "Iteration: " << std::setw(3) << nit << " "
                << "Value: " << change;
      if (converged) iter_info << " CONVERGED\n";
      std::cout << iter_info.str() << "\n";
    }

    if (converged) break;
  }
}


//###########################################################################


void
NeutronDiffusion::SteadyStateSolver::
solve_full_system(SourceFlags source_flags)
{
  if (verbose)
    std::cout << "\n***** Solving Full System\n";

  set_source(groupsets.front(), source_flags);

  Vector& b = groupsets.front().b;
  phi = linear_solver->solve(b);
}


//###########################################################################


void
NeutronDiffusion::SteadyStateSolver::
create_linear_solver(Groupset& groupset)
{
  SparseMatrix& A = groupset.A;
  switch (linear_solver_type)
  {
    case LinearSolverType::LU:
    {
      linear_solver = std::make_shared<SparseLU>(A);
      break;
    }
    case LinearSolverType::CHOLESKY:
    {
      linear_solver = std::make_shared<SparseCholesky>(A);
      break;
    }
    case LinearSolverType::JACOBI:
    {
      linear_solver = std::make_shared<Jacobi>(A, linear_solver_opts);
      break;
    }
    case LinearSolverType::GAUSS_SEIDEL:
    {
      linear_solver = std::make_shared<GaussSeidel>(A, linear_solver_opts);
      break;
    }
    case LinearSolverType::SOR:
    {
      linear_solver = std::make_shared<SOR>(A, linear_solver_opts);
      break;
    }
    case LinearSolverType::SSOR:
    {
      linear_solver = std::make_shared<SSOR>(A, linear_solver_opts);
      break;
    }
    case LinearSolverType::CG:
    {
      linear_solver = std::make_shared<CG>(A, linear_solver_opts);
      break;
    }
    default:throw std::runtime_error("Invalid linear solver type encountered.");
  }
}

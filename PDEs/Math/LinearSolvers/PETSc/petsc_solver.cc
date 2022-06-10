#include "petsc_solver.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include "PETScUtils/petsc_utils.h"

#include <iomanip>
#include <algorithm>
#include <cassert>


using namespace Math;
using namespace LinearSolver;


PETScSolver::PETScSolver(const std::string solver_type,
                         const std::string preconditioner_type,
                         const Options& opts) :
  solver_type(solver_type),
  preconditioner_type(preconditioner_type),
  relative_residual_tolerance(opts.tolerance),
  max_iterations(opts.max_iterations),
  verbosity(opts.verbosity)
{}


void
PETScSolver::set_matrix(const SparseMatrix& matrix)
{
  PETScUtils::CreateMatrix(A, matrix);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetType(ksp, solver_type.c_str());

  KSPGetPC(ksp, &pc);
  PCSetType(pc, preconditioner_type.c_str());

  KSPSetTolerances(ksp, relative_residual_tolerance,
                   PETSC_DEFAULT, PETSC_DEFAULT, max_iterations);

  if (verbosity > 1)
    KSPMonitorSet(ksp, &KSPMonitor, NULL, NULL);

  KSPSetUp(ksp);
}


void
PETScSolver::solve(Vector& x, const Vector& b) const
{
  Vec rhs, solution;
  PETScUtils::CreateVector(rhs, b);
  PETScUtils::CreateVector(solution, x);

  KSPSolve(ksp, rhs, solution);

  // Check convergence
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  PetscInt it;
  KSPGetIterationNumber(ksp, &it);

  PetscReal res;
  KSPGetResidualNorm(ksp, &res);

  bool converged =
    (reason == KSP_CONVERGED_RTOL)? true : false;

  std::string solver_str = solver_type;
  std::transform(solver_str.begin(), solver_str.end(),
                 solver_str.begin(), ::toupper);

  std::string pc_str = preconditioner_type;
  std::transform(pc_str.begin(), pc_str.end(),
                 pc_str.begin(), ::toupper);

  if (converged && verbosity == 1)
    std::cout << "PETSc::" << solver_str << "::"
              << pc_str << "::  CONVERGED   "
              << "Iteration:  " << it << "   "
              << "Value:  " << res << std::endl;

  if (!converged)
  {
    std::stringstream err;
    err << "!!*!! " << "PETSc::"
        << solver_str << "::"
        << pc_str << " FAILURE !!*!!\n"
        << "# of Iterations:   " << it << std::endl
        << "Final Difference:  " << res << std::endl;
    throw std::runtime_error(err.str());
  }

  PETScUtils::CopyToVector(solution, x);

  VecDestroy(&rhs);
  VecDestroy(&solution);
}



PetscErrorCode
PETScSolver::KSPMonitor(KSP solver, PetscInt it,
                        PetscReal rnorm,
                        void* monitordestroy)
{
  Vec rhs;
  KSPGetRhs(solver, &rhs);

  double rhs_norm;
  VecNorm(rhs, NORM_2, &rhs_norm);
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  std::stringstream ss;
  ss << "PETSc::"
     << "Iteration:  " << std::setw(4) << it << "   "
     << "Value:   " << std::scientific
     << std::setprecision(7) << rnorm
     << std::endl;

  std::cout << ss.str();
  return 0;
}


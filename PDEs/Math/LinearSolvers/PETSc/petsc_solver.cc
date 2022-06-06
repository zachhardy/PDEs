#include "petsc_solver.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include "PETScUtils/petsc_utils.h"

#include <iomanip>


using namespace Math;


LinearSolver::PETScSolver::
PETScSolver(const std::string solver_type,
            const std::string preconditioner_type,
            const Options& opts) :
  solver_type(solver_type),
  preconditioner_type(preconditioner_type),
  relative_residual_tolerance(opts.tolerance),
  max_iterations(opts.max_iterations),
  verbosity(opts.verbosity)
{}


PetscErrorCode KSPMonitorAChiTech(
  KSP ksp, PetscInt n, PetscReal rnorm, void *monitordestroy)
{

  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  const auto ksp_name = "Diffusion";

  std::stringstream buff;
  buff
    << ksp_name
    << " iteration "
    << std::setw(4) << n
    << " - Residual "
    << std::scientific << std::setprecision(7) << rnorm // / rhs_norm
    << std::endl;

  std::cout << buff.str();
  return 0;
}


void
LinearSolver::PETScSolver::set_matrix(const SparseMatrix& matrix)
{
  PETScUtils::CreateMatrix(A, matrix);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetType(ksp, solver_type.c_str());

  KSPGetPC(ksp, &pc);
  PCSetType(pc, preconditioner_type.c_str());

  KSPSetTolerances(ksp, relative_residual_tolerance,
                   PETSC_DEFAULT, PETSC_DEFAULT, max_iterations);

  KSPMonitorSet(ksp, &KSPMonitorAChiTech, NULL, NULL);

  KSPSetUp(ksp);
}


void
LinearSolver::PETScSolver::solve(Vector& x, const Vector& b) const
{
  Vec rhs, solution;
  PETScUtils::CreateVector(rhs, b);
  PETScUtils::CreateVector(solution, x);

  KSPSolve(ksp, rhs, solution);

  if (verbosity > 0)
  {
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    PetscInt its;
    KSPGetIterationNumber(ksp, &its);


    std::stringstream ostr;
    switch (reason)
    {
      case KSP_CONVERGED_RTOL_NORMAL     :
        ostr << "KSP_CONVERGED_RTOL_NORMAL";
        break;
      case KSP_CONVERGED_ATOL_NORMAL     :
        ostr << "KSP_CONVERGED_ATOL_NORMAL";
        break;
      case KSP_CONVERGED_RTOL            :
        ostr << "KSP_CONVERGED_RTOL";
        break;
      case KSP_CONVERGED_ATOL            :
        ostr << "KSP_CONVERGED_ATOL";
        break;
      case KSP_CONVERGED_ITS             :
        ostr << "KSP_CONVERGED_ITS";
        break;
      default:
        ostr << "UNRECOGNIZED";
        break;
    }
    std::cout << ostr.str() << "   "
              << "Iterations:  " << its << std::endl;
  }

  PETScUtils::CopyToVector(solution, x);

  VecDestroy(&rhs);
  VecDestroy(&solution);
}

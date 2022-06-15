#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H

#include "LinearSolvers/linear_solver.h"
#include "PETScUtils/petsc_utils.h"

#include <petscksp.h>


namespace Math::LinearSolver
{

  /** Implementation of a PETSc solver. */
  class PETScSolver : public LinearSolverBase<SparseMatrix>
  {
  protected:
    //========== Solver parameters
    size_t verbosity;

    double relative_residual_tolerance;
    size_t max_iterations;

    std::string solver_type = KSPCG;
    std::string preconditioner_type = PCNONE;

    //========== PETSc solver objects
    Mat A;
    KSP ksp;
    PC pc;

  public:
    PETScSolver(const std::string solver_type = KSPCG,
                const std::string preconditioner_type = PCNONE,
                const Options& opts = Options());

    /** Attach a matrix to the solver. */
    void set_matrix(const SparseMatrix& matrix) override;

    /** Solve the system using PETSc. */
    void solve(Vector& x, const Vector& b) const override;


    using LinearSolverBase::solve;

  private:

    static PetscErrorCode
    KSPMonitor(KSP solver, PetscInt it,
               PetscReal rnorm,
               void* monitordestroy);

    static PetscErrorCode
    KSPConvergenceTest(KSP solver, PetscInt it,
                       PetscReal rnorm,
                       KSPConvergedReason* reason, void*);
  };
}

#endif //PETSC_SOLVER_H

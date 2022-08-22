#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H

#include "LinearSolvers/linear_solver.h"
#include "PETScUtils/petsc_utils.h"

#include <petscksp.h>


namespace Math::LinearSolver
{

  /**
   * Implementation of a PETSc solver.
   * */
  class PETScSolver : public LinearSolverBase<SparseMatrix>
  {
  public:
    using LinearSolverBase::solve;

    /**
     * Default constructor.
     *
     * \param solver_type A PETSc solver type macro.
     * \param preconditioner_type A PETSc preconditioner type macro.
     */
    PETScSolver(const std::string solver_type = KSPCG,
                const std::string preconditioner_type = PCNONE,
                const Options& opts = Options());

    /**
     * Attach a sparse matrix to the solver.
     */
    void
    set_matrix(const SparseMatrix& matrix) override;

    /**
     * Solve the system using PETSc.
     */
    void
    solve(Vector& x, const Vector& b) const override;

  protected:
    /**
     * The PETSc matrix.
     */
    Mat A;

    /**
     * The PETSc linear Krylov solver.
     */
    KSP ksp;

    /**
     * The PETSc preconditioner.
     */
    PC pc;

    /**
     * The PETSc macro that specifies the solver type.
     */
    std::string solver_type = KSPCG;

    /**
     * The PETSc macro that specifies the preconditioner type.
     */
    std::string preconditioner_type = PCNONE;

    /**
        * The iteration convergence tolerance.
        */
    double tolerance;

    /**
     * The maximum number of iterations to attempt to converge the solution.
     */
    unsigned int max_iterations;

    /**
     * The level of screen output from the linear solver.
     */
    unsigned int verbosity = 0;

  private:

    /**
     * A routine used to monitor the progress of the PETSc solver.
     */
    static PetscErrorCode
    KSPMonitor(KSP solver, PetscInt it, PetscReal rnorm,
               void* monitordestroy);

    /**
     * A routine used to check the convergence of the PETSc solver.
     */
    static PetscErrorCode
    KSPConvergenceTest(KSP solver, PetscInt it, PetscReal rnorm,
                       KSPConvergedReason* reason, void*);
  };
}

#endif //PETSC_SOLVER_H

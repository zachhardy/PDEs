#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H

#include "LinearSolvers/linear_solver.h"
#include "PETScUtils/petsc_utils.h"

#include <petscksp.h>


namespace PDEs
{
  namespace Math
  {
    namespace LinearSolvers
    {

      /**
       * Implementation of a PETSc solver.
       */
      class PETScSolver : public LinearSolverBase<SparseMatrix>
      {
      protected:
        Mat A;
        KSP ksp;
        PC pc;

        std::string solver_type = KSPCG;
        std::string preconditioner_type = PCNONE;

        double tolerance;
        unsigned int max_iterations;
        unsigned int verbosity = 0;

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

        /** Attach a sparse matrix to the solver. */
        void set_matrix(const SparseMatrix& matrix) override;

        /** Solve the system using PETSc. */
        void solve(Vector& x, const Vector& b) const override;

      private:

        /**A routine used to monitor the progress of the PETSc solver. */
        static PetscErrorCode
        KSPMonitor(KSP solver, PetscInt it,
                   PetscReal rnorm, void* monitordestroy);

        /** A routine used to check the convergence of the PETSc solver. */
        static PetscErrorCode
        KSPConvergenceTest(KSP solver, PetscInt it,
                           PetscReal rnorm, KSPConvergedReason* reason, void*);
      };
    }
  }
}

#endif //PETSC_SOLVER_H

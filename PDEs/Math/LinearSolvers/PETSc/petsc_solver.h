#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H

#include "LinearSolvers/linear_solver.h"
#include "PETScUtils/petsc_utils.h"

#include <petscksp.h>


namespace Math::LinearSolver
{

  /**
   * A class for a PETSc linear solver.
   */
  class PETScSolver : public LinearSolverBase<SparseMatrix>
  {

  protected:
    //========== Solver parameters
    size_t verbosity = 0;

    std::string solver_type = KSPCG;
    std::string preconditioner_type = PCNONE;

    double relative_residual_tolerance = 1.0e-6;
    size_t max_iterations = 200;

    //========== PETSc solver objects
    Mat   A;
    KSP   ksp;
    PC    pc;

  public:
    /**
     * Default constructor.
     */
    PETScSolver(const std::string solver_type = KSPCG,
                const std::string preconditioner_type = PCNONE,
                const Options& opts = Options());

    /**
     * Attach a matrix to the solver.
     */
    void
    set_matrix(const SparseMatrix& matrix) override;

    /**
     * Solve the system using PETSc.
     */
    void
    solve(Vector& x, const Vector& b) const override;



  };
}

#endif //PETSC_SOLVER_H

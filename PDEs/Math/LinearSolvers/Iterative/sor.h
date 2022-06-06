#ifndef SOR_H
#define SOR_H

#include "LinearSolvers/linear_solver.h"


namespace Math::LinearSolver
{

  /**
   * Implementation of a successive over relaxation (SOR) solver.
   */
  class SOR : public IterativeSolverBase
  {
  protected:
    double omega;

  public:
    /**
     * Default constructor.
     */
    SOR(const double omega = 1.5,
        const Options& opts = Options(),
        const std::string solver_name = "SOR");

    /**
     * Solve the system using the SOR iterative method.
     */
    virtual void
    solve(Vector& x, const Vector& b) const override;


    using LinearSolverBase::solve;
  };

}

#endif //SOR_H

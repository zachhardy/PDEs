#ifndef JACOBI_H
#define JACOBI_H

#include "LinearSolvers/linear_solver.h"


namespace Math::LinearSolver
{

  /** Implementation of the Jacobi iterative method. */
  class Jacobi : public IterativeSolverBase
  {
  public:
    Jacobi(const Options& opts = Options());

    /** Solve the system using the Jacobi iterative method. */
    void solve(Vector& x, const Vector& b) const override;


    using LinearSolverBase::solve;
  };
}

#endif //JACOBI_H

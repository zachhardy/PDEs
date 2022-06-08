#ifndef CG_H
#define CG_H

#include "LinearSolvers/linear_solver.h"


namespace Math::LinearSolver
{

  /** Implementation of the conjugate gradient (CG) method. */
  class CG : public IterativeSolverBase
  {
  public:
    CG(const Options& opts = Options());

    /** Solve the system using the CG method. */
    void solve(Vector& x, const Vector& b) const override;


    using LinearSolverBase::solve;
  };

}

#endif //CG_H

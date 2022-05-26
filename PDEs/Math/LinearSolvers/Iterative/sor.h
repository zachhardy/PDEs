#ifndef SOR_H
#define SOR_H

#include "LinearSolvers/linear_solver.h"

#include <cstddef>


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of a successive over relaxation (SOR) solver.
 */
class SOR : public IterativeSolver
{
protected:
  double omega;

public:

  /**
   * Default constructor.
   */
  SOR(const SparseMatrix& A,
      const double omega = 1.7,
      const double tolerance = 1.0e-8,
      const size_t max_iteration = 1000,
      const bool verbose = false);

  /**
   * Solve the system using the SOR iterative method.
   */
  void
  solve(Vector& x, const Vector& b) const override;


  using LinearSolverBase::solve;
};

}

#endif //SOR_H

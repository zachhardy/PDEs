#ifndef SOR_H
#define SOR_H

#include "LinearSolvers/linear_solver.h"

#include <cstddef>


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of a successive over relaxation (SOR) solver.
 */
class SOR : public LinearSolverBase
{
private:
  const SparseMatrix& A;
  double tolerance;
  size_t max_iterations;

  double omega;

public:
  /**
   * Default constructor.
   */
  SOR(const SparseMatrix& A,
      const double tolerance = 1.0e-8,
      const size_t max_iteration = 1000,
      const double omega = 0.5);

  /**
   * Solve the system using the SOR iterative method.
   */
  void
  solve(Vector& x, const Vector& b) const override;
};

}

#endif //SOR_H
#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "LinearSolvers/linear_solver.h"

#include <cstddef>

namespace pdes::Math::LinearSolver
{

/**
 * Implementation of the Gauss Seidel iterative method.
 */
class GaussSeidel : public LinearSolverBase
{
private:
  const SparseMatrix& A;
  double tol;
  size_t maxiter;

public:
  /**
   * Default constructor.
   */
  GaussSeidel(const SparseMatrix& A,
              const double tolerance = 1.0e-8,
              const size_t max_iterations = 1000);

  /**
   * Solve the system using the Gauss Seidel iterative method.
   */
  void
  solve(const Vector& b, Vector& x) const override;

  /**
   * Return the solution of the Gauss Seidel solve.
   */
  Vector
  solve(const Vector& b) const override;
};

}


#endif //GAUSS_SEIDEL_H

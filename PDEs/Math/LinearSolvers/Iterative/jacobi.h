#ifndef JACOBI_H
#define JACOBI_H

#include "LinearSolvers/linear_solver.h"

#include <cstddef>


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of the Jacobi iterative method.
 */
class Jacobi : public IterativeSolverBase
{
public:

  /**
   * Default constructor.
   */
  Jacobi(const SparseMatrix& A,
         const double tolerance = 1.0e-8,
         const size_t max_iterations = 1000,
         const bool verbose = false);

  /**
   * Solve the system using the Jacobi iterative method.
   */
  void
  solve(Vector& x, const Vector& b) const override;


  using LinearSolverBase::solve;
};
}

#endif //JACOBI_H

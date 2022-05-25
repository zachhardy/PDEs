#ifndef JACOBI_H
#define JACOBI_H

#include "sparse_matrix.h"
#include "vector.h"

#include <cstddef>

namespace pdes::Math
{

class JacobiSolver
{
private:
  double tol = 1.0e-8;
  size_t max_iter = 1000;

public:
  /**
   * Default constructor.
   */
  JacobiSolver() = default;

  /**
   * Constructor with specified iteration controls.
   */
  JacobiSolver(const double tolerance,
               const size_t max_iterations);

  /**
   * Solve the system using the Jacobi iterative method.
   */
  void
  solve(const SparseMatrix& A, const Vector& b, Vector& x);

};
}

#endif //JACOBI_H

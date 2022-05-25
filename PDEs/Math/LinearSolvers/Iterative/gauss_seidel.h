#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "sparse_matrix.h"
#include "vector.h"
#include "linear_solver.h"

namespace pdes::Math
{

class GaussSeidelSolver
{
public:
  static const LinearSolverType type = LinearSolverType::ITERATIVE;

private:
  double tol;
  size_t maxiter;

public:
  /**
   * Constructor with specified iteration controls.
   */
  GaussSeidelSolver(const double tolerance = 1.0e-8,
                    const size_t max_iterations = 1000);

  /**
   * Solve the system using the Gauss Seidel iterative method.
   */
  void
  solve(const SparseMatrix& A, Vector& x, const Vector& b) const;

  /**
   * Return the solution of the Gauss Seidel solve.
   */
  Vector
  solve(const SparseMatrix& A, const Vector& b) const;
};

}


#endif //GAUSS_SEIDEL_H

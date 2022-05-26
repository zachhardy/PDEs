#ifndef JACOBI_H
#define JACOBI_H

#include "sparse_matrix.h"
#include "vector.h"
#include "linear_solver.h"

#include <cstddef>

namespace pdes::Math
{

class JacobiSolver : public LinearSolverBase
{
public:
  static const LinearSolverType type = LinearSolverType::DIRECT;

private:
  const SparseMatrix& A;
  double tol;
  size_t maxiter;

public:
  /**
   * Default constructor.
   */
  JacobiSolver(const SparseMatrix& A,
               const double tolerance = 1.0e-8,
               const size_t max_iterations = 1000);

  /**
   * Solve the system using the Jacobi iterative method.
   */
  void
  solve(const Vector& b, Vector& x) const override;

  /**
   * Return the solution of the Jacobi solve.
   * \see Jacobi::solve
   */
  Vector
  solve(const Vector& b) const override;

};
}

#endif //JACOBI_H

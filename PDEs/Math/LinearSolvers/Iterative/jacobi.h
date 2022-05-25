#ifndef JACOBI_H
#define JACOBI_H

#include "sparse_matrix.h"
#include "vector.h"
#include "linear_solver.h"

#include <cstddef>

namespace pdes::Math
{

class JacobiSolver
{
public:
  static const LinearSolverType type = LinearSolverType::DIRECT;

private:
  double tol;
  size_t maxiter;

public:
  /**
   * Constructor with specified iteration controls.
   */
  JacobiSolver(const double tolerance = 1.0e-8,
               const size_t max_iterations = 1000);

  /**
   * Solve the system using the Jacobi iterative method.
   */
  void
  solve(const SparseMatrix& A, Vector& x, const Vector& b) const;

  /**
   * Return the solution of the Jacobi solve.
   * \see Jacobi::solve
   */
  Vector
  solve(const SparseMatrix& A, const Vector& b)const ;

};
}

#endif //JACOBI_H

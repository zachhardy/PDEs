#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "sparse_matrix.h"
#include "vector.h"
#include "linear_solver.h"

namespace pdes::Math
{

class GaussSeidelSolver : public LinearSolverBase
{
public:
  using value_type = SparseMatrix::value_type;

private:
  const SparseMatrix& A;
  value_type tol;
  size_t maxiter;

public:
  /**
   * Default constructor.
   */
  GaussSeidelSolver(const SparseMatrix& A,
                    const value_type tolerance = 1.0e-8,
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

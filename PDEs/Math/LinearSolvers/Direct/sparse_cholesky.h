#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#include "LinearSolvers/linear_solver.h"


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of a sparse Cholesky solver.
 */
class SparseCholesky : public LinearSolverBase
{
private:
  SparseMatrix& A;
  bool factorized = false;

public:

  /**
   * Default constructor.
   */
  SparseCholesky(SparseMatrix& A, const bool verbose = false);


  /**
   * Perform a Cholesky factorization on the matrix \f$ \boldsymbol{A} \f$.
   * \see Cholesky::solve
   */
  void
  factorize();

  /**
   * Solve the Cholesky factored linear system.
   * \see Cholesky::solve
   */
  void
  solve(Vector& x, const Vector& b) const override;


  using LinearSolverBase::solve;
};

}
#endif //SPARSE_CHOLESKY_H

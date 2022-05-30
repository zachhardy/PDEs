#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#include "LinearSolvers/linear_solver.h"


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of a sparse Cholesky solver.
 */
class SparseCholesky : public DirectSolverBase
{
public:
  /**
   * Default constructor.
   */
  SparseCholesky(SparseMatrix& A);


  /**
   * Perform a Cholesky factorization on the matrix \f$ \boldsymbol{A} \f$.
   * \see Cholesky::solve
   */
  void
  factorize() override;

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

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "../linear_solver.h"
#include "matrix.h"


namespace Math::LinearSolver
{

  /** Implementation of a Cholesky decomposition solver. */
  class Cholesky : public DirectSolverBase<Matrix>
  {
  public:
    Cholesky();

    /**
     * Perform a Cholesky factorization on the matrix \f$ \boldsymbol{A} \f$.
     *
     * Cholesky factorization is meant for symmetric positive definite matrices
     * into a lower triangular matrix and its transpose. This method is more
     * efficient than the LU decomposition when applicable.
     *
     * \note Checks are not performed to ensure symetric positive definiteness.
     *    The user is responsible for ensuring the matrix fits this criteria.
     */
    void factorize() override;

    /**
     * Solve the Cholesky factored linear system.
     *
     * The Cholesky solve is a specialization of the LU solve in that
     * \f$ \boldsymbol{U} = \boldsymbol{L}^T \f$.
     *
     * \see LU::solve
     */
    void solve(Vector& x, const Vector& b) const override;


    using LinearSolverBase::solve;
  };

}

#endif //CHOLESKY_H

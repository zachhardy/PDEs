#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "matrix.h"
#include "linear_solver.h"

namespace pdes::Math
{

/**
 * A class for a Choleky decomposition solver.
 */
class Cholesky : public Matrix
{
public:
  using value_type = Matrix::value_type;
  static const LinearSolverType type = LinearSolverType::DIRECT;

private:
  bool factorized = false;

public:
  /**
   * Copy construction from a matrix.
   */
  Cholesky(const Matrix& other);

  /**
   * Move construction from a matrix.
   */
  Cholesky(Matrix&& other);

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
  void
  factorize();

  /**
   * Solve the Cholesky factored linear system.
   *
   * The Cholesky solve is a specialization of the LU solve in that
   * \f$ \boldsymbol{U} = \boldsymbol{L}^T \f$. See \ref lu_solve for
   * implementation detail.
   *
   * \param b A vector of length \f$ n \f$.
   * \param x The destination vector.
   * \return The solution \f$ \vec{x} \f$ of
   *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
   */
  void
  solve(const Vector& b, Vector& x) const;

  /**
   * Return the solution of the Cholesky solve.
   * \see Cholesky::solve
   */
  Vector
  solve(const Vector& b) const;
};

}

#endif //CHOLESKY_H

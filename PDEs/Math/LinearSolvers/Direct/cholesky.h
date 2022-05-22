#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "matrix.h"

namespace pdes::Math
{

/** A class for a Choleky decomposition solver. */
class Cholesky : public Matrix
{
public:
  using value_type = Matrix::value_type;

private:
  bool factorized = false;

public:
  /**
   * Copy construction from a matrix.
   */
  Cholesky(const Matrix& other) : Matrix(other) {}

  /**
   * Move construction from a matrix.
   */
  Cholesky(Matrix&& other) : Matrix(std::move(other)) {}

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
  factorize()
  {
    size_t n = n_rows();
    for (size_t j = 0; j < n; ++j)
    {
      value_type sum = 0.0;
      value_type* a_j = data(j);
      for (size_t k = 0; k < j; ++k)
        sum += std::pow(a_j[k], 2.0);
      a_j[j] = std::sqrt(a_j[j] - sum);

      for (size_t i = j + 1; i < n; ++i)
      {
        sum = 0.0;
        value_type* a_i = data(i);
        for (size_t k = 0; k < j; ++k)
          sum += a_i[k] * a_j[k];
        a_i[j] = (a_i[j] - sum) / a_j[j];
      }
    }
    factorized = true;
  }

  /**
   * Solve the Cholesky factored linear system.
   *
   * The Cholesky solve is a specialization of the LU solve in that
   * \f$ \boldsymbol{U} = \boldsymbol{L}^T \f$. See \ref lu_solve for
   * implementation detail.
   *
   * \param b A vector of length \f$ n \f$.
   * \return The solution \f$ \vec{x} \f$ of
   *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
   */
  void
  solve(const Vector& b, Vector& x)
  {
    Assert(b.size() == n_rows(), "Dimension mismatch error.");
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    if (!factorized)
      factorize();

    //================================================== Forward solve
    size_t n = n_rows();
    for (size_t i = 0; i < n; ++i)
    {
      value_type value = b[i];
      const value_type* a_i = data(i);
      for (size_t j = 0; j < i; ++j)
        value -= a_i[j] * x[j];
      x[i] = value / diagonal(i);
    }

    //================================================== Backward solve
    for (size_t i = n - 1; i != -1; --i)
    {
      x[i] /= this->diagonal(i);
      const value_type* a_ij = data(i);
      for (size_t j = 0; j < i; ++j)
        x[j] -= *a_ij++ * x[i];
    }
  }

  Vector
  solve(const Vector& b)
  {
    Vector x(n_cols());
    solve(b, x);
    return x;
  }
};

}

#endif //CHOLESKY_H

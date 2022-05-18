#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "matrix.h"

namespace math
{

/** A class for a Choleky decomposition solver. */
template<typename number = double>
class Cholesky : public Matrix<number>
{
public:
  using value_type = typename Matrix<number>::value_type;
  using size_type = typename Matrix<number>::size_type;

private:
  bool factorized = false;

public:
  /// Default constructor.
  Cholesky() = default;

  /// Copy constructor.
  Cholesky(const Cholesky& other)
    : Matrix<value_type>(other.values)
  {}

  /// Move constructor.
  Cholesky(Cholesky&& other)
    : Matrix<value_type>(std::move(other.values))
  {}

  /// Copy from a Matrix.
  Cholesky(const Matrix<value_type>& other)
    : Matrix<value_type>(other)
  {}

  /// Move from a Matrix.
  Cholesky(Matrix<value_type>&& other)
    : Matrix<value_type>(std::move(other))
  {}

  /**
   * Perform a Cholesky factorization on the matrix \f$ \boldsymbol{A} \f$.
   *
   * Cholesky factorization is meant for symmetric positive definite matrices into
   * a lower triangular matrix and its transpose. This method is more efficient
   * than the LU decomposition when applicable.
   *
   * \note Checks are not performed to ensure symetric positive definiteness. The
   *       user is responsible for ensuring the matrix fits this criteria.
   */
  void
  factorize()
  {
    size_type n = this->n_rows();
    for (size_type j = 0; j < n; ++j)
    {
      // Set the diagonal element
      value_type sum = 0.0;
      for (size_type k = 0; k < j; ++k)
        sum += (*this)(j, k) * (*this)(j, k);
      (*this)(j, j) = std::sqrt((*this)(j, j) - sum);

      // Set the off-diagonals
      for (size_type i = j + 1; i < n; ++i)
      {
        sum = 0.0;
        for (size_type k = 0; k < j; ++k)
          sum += (*this)(i, k) * (*this)(j, k);
        (*this)(i, j) = ((*this)(i, j) - sum) / (*this)(j, j);
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
  solve(const Vector<value_type>& b, Vector<value_type>& x)
  {
    Assert(b.size() == this->n_rows(), "Dimension mismatch error.");
    Assert(x.size() == this->n_cols(), "Dimension mismatch error.");
    if (!factorized)
      factorize();

    size_type n = this->n_rows();
    value_type value = 0.0;

    //================================================== Forward solve
    for (size_type i = 0; i < n; ++i)
    {
      value = b[i];
      const value_type* a_ij = this->values[i].data();
      for (size_type j = 0; j < i; ++j)
        value -= *a_ij++ * x[j];
      x[i] = value / (*this)(i, i);
    }

    //================================================== Backward solve
    for (size_type i = n - 1; i != -1; --i)
    {
      value = x[i];
      const value_type* a_ij = this->values[i].data();
      for (size_type j = i + 1; j < n; ++j)
        value -= *a_ij++ * x[j];
      x[i] = value / (*this)(i, i);
    }
  }

  Vector<value_type>
  solve(const Vector<value_type>& b)
  {
    Vector<value_type> x(this->n_cols());
    solve(b, x);
    return x;
  }


};

}

#endif //CHOLESKY_H

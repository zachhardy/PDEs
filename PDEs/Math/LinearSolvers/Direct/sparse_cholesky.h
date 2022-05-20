#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#include "sparse_matrix.h"

namespace math
{

/** A class for a Choleky decomposition solver. */
template<typename number = double>
class SparseCholesky : public SparseMatrix<number>
{
public:
  using value_type = typename SparseMatrix<number>::value_type;
  using size_type = typename SparseMatrix<number>::size_type;

  using iterator = typename SparseMatrix<number>::iterator;
  using entry = typename SparseMatrix<number>::entry;

private:
  bool factorized = false;

public:
  /// Default constructor.
  SparseCholesky() = default;

  /// Copy constructor.
  SparseCholesky(const SparseCholesky& other)
      : SparseMatrix<value_type>(other.values)
  {}

  /// Move constructor.
  SparseCholesky(SparseCholesky&& other)
      : SparseMatrix<value_type>(std::move(other.values))
  {}

  /// Copy from a Matrix.
  SparseCholesky(const SparseMatrix<value_type>& other)
      : SparseMatrix<value_type>(other)
  {}

  /// Move from a Matrix.
  SparseCholesky(SparseMatrix<value_type>&& other)
      : SparseMatrix<value_type>(std::move(other))
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
      value_type* d = this->locate(j, j);
      Assert(d && *d != 0.0,
             "Singular matrix error.");

      // Set the diagonal
      value_type sum = 0.0;
      for (const auto elem : this->const_row_iterator(j))
        if (elem.column < j)
          sum += elem.value * elem.value;
      *d = std::sqrt(*d - sum);

      // Set the off-diagonals
      for (size_type i = j + 1; i < n; ++i)
      {
        sum = 0.0;
        for (const auto a_ik : this->const_row_iterator(i))
          if (a_ik.column < j)
            for (const auto a_jk : this->const_row_iterator(j))
              if (a_jk.column == a_ik.column)
                sum += a_ik.value * a_jk.value;

        value_type* a_ij = this->locate(i, j);
        value_type value = (a_ij) ? (*a_ij - sum) / *d : -sum/ *d;
        if (std::fabs(value) != 0.0)
          this->set(i, j, value);
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

    //================================================== Forward solve
    size_type n = this->n_rows();
    Vector<value_type> y(n);
    for (size_type i = 0; i < n; ++i)
    {
      value_type value = b[i];
      for (const auto elem : this->const_row_iterator(i))
        if (elem.column < i)
          value -= elem.value * x[elem.column];
      x[i] = value / (*this)(i, i);
    }

    //================================================== Backward solve
    for (size_type i = n - 1; i != -1; --i)
    {
      x[i] /= *this->diagonal(i);
      for (const auto a_ij : this->const_row_iterator(i))
        if (a_ij.column < i)
          x[a_ij.column] -= a_ij.value * x[i];
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
#endif //SPARSE_CHOLESKY_H

#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#include "sparse_matrix.h"

namespace pdes::Math
{

/** A class for a Choleky decomposition solver. */
template<typename number = double>
class SparseCholesky : public SparseMatrix
{
public:
  using value_type = typename SparseMatrix::value_type;

  using entry = typename SparseMatrix::iterator::entry;
  using const_entry = typename SparseMatrix::const_iterator::const_entry;

private:
  bool factorized = false;

public:
  /**
   * Copy construction from a sparse matrix.
   */
  SparseCholesky(const SparseMatrix& other) : SparseMatrix(other) {}

  /**
   * Move construction from a sparse matrix.
   */
  SparseCholesky(SparseMatrix&& other) : SparseMatrix(std::move(other)) {}

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
      value_type* d = locate(j, j);
      Assert(d && *d != 0.0, "Singular matrix error.");

      // Set the diagonal
      value_type sum = 0.0;
      for (const_entry elem : const_row_iterator(j))
        if (elem.column < j)
          sum += elem.value * elem.value;
      *d = std::sqrt(*d - sum);

      // Set the off-diagonals
      for (size_t i = j + 1; i < n; ++i)
      {
        sum = 0.0;
        for (const_entry a_ik : const_row_iterator(i))
          if (a_ik.column < j)
            for (const_entry a_jk : const_row_iterator(j))
              if (a_jk.column == a_ik.column)
                sum += a_ik.value * a_jk.value;

        value_type* a_ij = locate(i, j);
        value_type value = (a_ij) ? (*a_ij - sum) / *d : -sum/ *d;
        if (std::fabs(value) != 0.0)
          set(i, j, value);
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
  solve(const Vector& b, Vector& x) const
  {
    Assert(factorized, "Matrix must be factorized before solving.");
    Assert(b.size() == n_rows(), "Dimension mismatch error.");
    Assert(x.size() == n_cols(), "Dimension mismatch error.");

    //================================================== Forward solve
    size_t n = n_rows();
    Vector y(n);
    for (size_t i = 0; i < n; ++i)
    {
      value_type value = b[i];
      for (const_entry el : const_row_iterator(i))
        if (el.column < i)
          value -= el.value * x[el.column];
      x[i] = value / (*this)(i, i);
    }

    //================================================== Backward solve
    for (size_t i = n - 1; i != -1; --i)
    {
      x[i] /= *diagonal(i);
      for (const_entry a_ij : const_row_iterator(i))
        if (a_ij.column < i)
          x[a_ij.column] -= a_ij.value * x[i];
    }
  }

  Vector
  solve(const Vector& b) const
  {
    Vector x(n_cols());
    solve(b, x);
    return x;
  }
};

}
#endif //SPARSE_CHOLESKY_H

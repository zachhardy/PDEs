#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include "matrix.h"
#include "vector.h"

#include <cinttypes>


namespace math
{

/**
 * Solve a system using Gaussian elimination.
 *
 * This routine factors the matrix and right-hand side into row-echelon form,
 * then uses back substitution to directly obtain the solution. Row-echelon
 * form creates an upper triangular system. This is done by traversing through
 * each column and performing row-wise operations to eliminate all sub-diagonal
 * coefficients. For numerical stability, pivoting can be used to ensure that
 * for each subsequent column, division by the largest element is carried out.
 *
 * \param A An \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \param pivot A flag for whether pivoting is performed.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
template<typename value_type>
Vector<value_type>
gaussian_elimination(Matrix<value_type>& A,
                     Vector<value_type>& b,
                     const bool pivot)
{
  Assert(A.n_rows() == A.n_cols(),
         "Square matrices are required for Gaussian elimination.");
  Assert(b.size() == A.n_rows(), "Dimension mismatch error.");

  uint64_t n = b.size();

  //======================================== Row-echelon factorization
  // Go through the columns of the matrix
  for (uint64_t j = 0; j < n; ++j)
  {
    /* Find the row index for the largest magnitude entry in this column.
     * This is only done for sub-diagonal elements. */
    value_type max = 0.0;
    uint64_t argmax = j;
    for (uint64_t i = j; i < n; ++i)
    {
      if (std::fabs(A[i][j]) > max)
      {
        max = std::fabs(A[i][j]);
        argmax = i;
      }
    }

    // If the sub-diagonal is uniformly zero, throw error
    Assert(A[argmax][j] != 0.0,
           "Degenerate matrix. "
           "Specifically, all sub-diagonal elements are zero.");

    /* Swap the current row and the row containing the largest magnitude
     * entry corresponding for the current column. This is done to improve
     * the numerical stability of the algorithm. */
    if (pivot and argmax != j)
    {
      std::swap(b[j], b[argmax]);
      A.swap_row(j, argmax);
    }

    /* Perform row-wise operations such that all sub-diagonal values are zero.
     * This is done by subtracting the current row times the ratio of the
     * sub-diagonal and the current row's leading value. */
    for (uint64_t i = j + 1; i < n; ++i)
    {
      value_type factor = A[i][j] / A[j][j];
      for (uint64_t k = j; k < n; ++k)
        A[i][k] -= A[j][k] * factor;
      b[i] -= b[j] * factor;
    }
  }

  //======================================== Back substitution solve
  Vector<value_type> x(n, 0.0);
  for (int i = n - 1; i >= 0; --i)
  {
    value_type value = b[i];
    for (int j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}
}
#endif //GAUSSIAN_ELIMINATION_H

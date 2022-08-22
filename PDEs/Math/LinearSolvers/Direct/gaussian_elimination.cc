#include "gaussian_elimination.h"

#include "vector.h"
#include "matrix.h"
#include "Sparse/sparse_matrix.h"

#include <cmath>
#include <cassert>

using namespace Math;


template<>
Vector
Math::gaussian_elimination<Matrix>(Matrix& A,
                                   Vector& b,
                                   const bool pivot)
{
  assert(A.n_rows() == A.n_cols());
  assert(b.size() == A.n_rows());
  size_t n = b.size();

  //======================================== Row-echelon factorization
  // Go through the columns of the matrix
  for (size_t j = 0; j < n; ++j)
  {
    // Find the row index for the largest magnitude entry in this column.
    // This is only done for sub-diagonal elements.
    if (pivot)
    {
      size_t argmax = j;
      double max = A(j, j);
      for (size_t k = j; k < n; ++k)
      {
        const double a_kj = A(k, j);
        if (std::fabs(a_kj) > max)
        {
          argmax = k;
          max = std::fabs(a_kj);
        }
      }

      // If the sub-diagonal is uniformly zero, throw error
      assert(max != 0.0);

      // Swap the current row and the row containing the largest magnitude
      // entry corresponding for the current column. This is done to improve
      // the numerical stability of the algorithm.
      if (argmax != j)
      {
        std::swap(b[j], b[argmax]);
        A.swap_row(j, argmax);
      }
    }//if pivot

    const double* a_j = A.data(j); // accessor for row j
    const double a_jj = a_j[j]; // diagonal element for row j

    // Perform row-wise operations such that all sub-diagonal values are zero.
    // This is done by subtracting the current row times the ratio of the
    // sub-diagonal and the current row's leading value.
    for (size_t i = j + 1; i < n; ++i)
    {
      double* a_i = A.data(i); // accessor for row i

      const double factor = a_i[j] / a_jj;
      for (size_t k = j; k < n; ++k, ++a_i)
        *a_i -= a_j[k] * factor;
      b[i] -= b[j] * factor;
    }
  }

  //======================================== Back substitution solve
  Vector x(n, 0.0);
  for (size_t i = n - 1; i != -1; --i)
  {
    const double* a_i = A.data(i); // accessor for row i
    const double a_ii = a_i[i]; // diagonal element for row i
    a_i += i + 1; // increment to first element after diagonal

    double value = b[i];
    for (size_t j = i + 1; j < n; ++j, ++a_i)
      value -= *a_i * x[j];
    x[i] = value/a_ii;
  }
  return x;
}


template<>
Vector
Math::gaussian_elimination<SparseMatrix>(SparseMatrix& A,
                                         Vector& b,
                                         const bool pivot)
{
  assert(A.n_rows() == A.n_cols());
  assert(b.size() == A.n_rows());
  size_t n = b.size();

  //======================================== Row-echelon factorization
  // Go through the columns of the matrix
  for (size_t j = 0; j < n; ++j)
  {
    // Find the row index for the largest magnitude entry in this column.
    // This is only done for sub-diagonal elements.
    if (pivot)
    {
      const double a_jj = A.diag_el(j);

      size_t argmax = j;
      double max = std::fabs(a_jj);
      for (size_t k = j + 1; k < n; ++k)
      {
        const double a_kj = A.el(k, j);
        if (a_kj > max)
        {
          argmax = k;
          max = std::fabs(a_kj);
        }
      }

      // If the sub-diagonal is uniformly zero, throw an error
      assert(max != 0.0);

      // Swap the current row and the row containing the largest magnitude
      // entry corresponding for the current column. This is done to improve
      // the numerical stability of the algorithm.
      if (argmax != j)
      {
        std::swap(b[j], b[argmax]);
        A.swap_row(j, argmax);
      }
    }//if pivot

    const double a_jj = A.diag(j);

    // Perform row-wise operations such that all sub-diagonal values are zero.
    // This is done by subtracting the current row times the ratio of the
    // sub-diagonal and the current row's leading value.
    for (size_t i = j + 1; i < n; ++i)
    {
      const double a_ij = A.el(i, j);
      if (a_ij != 0.0)
      {
        const double factor = a_ij/a_jj;
        for (const auto el : A.row_iterator(j))
          if (el.column >= j)
            A.add(i, el.column, -el.value * factor);
        b[i] -= b[j]*factor;
      }
    }
  }

  //======================================== Back substitution solve
  Vector x(n, 0.0);
  for (size_t i = n - 1; i != -1; --i)
  {
    double value = b[i];
    for (const auto el : A.row_iterator(i))
      if (el.column > i)
        value -= el.value * x[el.column];
    x[i] = value / A.diag(i);
  }
  return x;
}

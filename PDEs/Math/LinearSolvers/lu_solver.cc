#include "lu.h"

//######################################################################

/**
 * \brief Factor the matrix in-place.
 *
 * The LU factorization froms a row-echelon form, or upper triangular, matrix
 * (\ref math::GaussianElimination::setup) and a lower triangular matrix which
 * contains the row-wise operations used to form the row-echelon form matrix.
 * This factorization is performed in place and is given by
 * \f$ \boldsymbol{A} = \boldsymbol{L} \boldsymbol{U} \f$.
 */
void math::LU::setup()
{
  if (not initailized)
  {
    std::stringstream err;
    err << "math::LU::" << __FUNCTION__ << ": "
        << "No matrix available to compute the row-echelon form.";
    throw std::runtime_error(err.str());
  }

  size_t n = A.n_rows();

  // Initialize permutation matrix such that each row maps to itself.
  P.clear();
  for (size_t i = 0; i < n; ++i)
    P.emplace_back(i);

  // Continue while more columns need to be traversed.
  size_t column = 0;
  while (column < n)
  {
    // Find the row with the largest column entry magnitude.
    double max = 0.0;
    size_t argmax = column;
    for (size_t i = column; i < n; ++i)
    {
      if (fabs(A[i][column]) > max)
      {
        max = fabs(A[i][column]);
        argmax = i;
      }
    }

    // If sub-diagonal is uniformly zero, throw error
    if (A[argmax][column] == 0.0)
    {
      std::stringstream err;
      err << "math::LU::" << __FUNCTION__ << ": "
          << "Matrix is ill-defined. Column " << column << " is uniformly zero";
      throw std::runtime_error(err.str());
    }

    /* Swap the current row and the row containing the largest magnitude
     * entry corresponding for the current column. This is done to improve
     * the numerical stability of the algorithm. */
    if (with_pivoting and argmax != column)
    {
      std::cout << "Swapping row " << column << " with row " << argmax << ".\n";
      std::swap(P[column], P[argmax]);
      std::swap(b[column], b[argmax]);
      A.swap_row(column, argmax);
    }

    /* Perform row-wise operations such that all sub-diagonal values are zero.
     * This is done by subtracting the current row times the ratio of the
     * sub-diagonal and the current row's leading value. */
    for (size_t i = column + 1; i < n; ++i)
    {
      double factor = A[i][column] / A[column][column];
      for (size_t j = column + 1; j < n; ++j)
        A[i][j] -= A[column][j] * factor;
      A[i][column] = factor;
    }
    ++column;
  }
}

//######################################################################

/**
 * \brief Solve a system using an LU factorized matrix.
 *
 * The LU factorization yields the system
 * \f$
 *      \boldsymbol{A} \vec{x} =
 *      \boldsymbol{L} \boldsymbol{U} \vec{x} =
 *      \vec{b}
 * \f$. This is solved in a two step process. First, the system
 * \f$ \boldsymbol{L} \vec{y} = \vec{b} \f$, where
 * \f$ \vec{y} = \boldsymbol{U} \vec{x} \f$,
 * is solved using forward substitution. Next, the solution is obtained
 * by solving \f$ \boldsymbol{U} \vec{x} = \vec{y} \f$ using backward
 * substitution.
 */
Vector math::LU::solve()
{
  size_t n = A.n_rows();
  Vector x(n, 0.0);

  // First solve using forward substitution
  for (int i = 0; i < n; ++i)
  {
    double value = b[i];
    for (int j = 0; j < i; ++j)
      value -= A[i][j] * x[j];
    x[i] = value;
  }

  // Now, solve using backward substitution
  for (int i = n - 1; i >= 0; --i)
  {
    double value = x[i];
    for (int j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }

  // Map a pivoted solution back to the original ordering.
  if (with_pivoting)
    for (size_t i = 0; i < n; ++i)
      x[P[i]] = x[i];

  return x;
}

#include "gaussian_elimination.h"

//######################################################################
/**
 * The setup routine factors the system into row-echelon form so that
 * back-substitution can be used to solve the system.
 * \see row_echelon
 */
void math::GaussianElimination::setup()
{
  if (not initailized)
  {
    std::stringstream err;
    err << "GaussianElimination::" << __FUNCTION__ << ": "
        << "No matrix available to compute the row-echelon form.";
    throw std::runtime_error(err.str());
  }
  row_echelon();
}

//######################################################################

void math::GaussianElimination::row_echelon()
{
  size_t n = A.n_rows();
  for (size_t i = 0; i < n; ++i)
    P.emplace_back(i);

  // While the number of pivots is less than the dimension of the system...
  size_t pivot = 0;
  while (pivot < n)
  {
    /* Find the row with the largest pivot column entry magnitude. */
    double max_val = 0.0;
    size_t argmax = pivot;
    for (size_t i = pivot; i < n; ++i)
    {
      if (A[i][pivot] > max_val)
      {
        max_val = A[i][pivot];
        argmax = i;
      }
    }

    // If sub-diagonal is uniformly zero, go to next column
    if (A[argmax][pivot] == 0.0) ++pivot;

    // Otherwise, continue with algorithm
    else
    {
      /* Swap the current row and the row containing the largest magnitude
       * entry corresponding to the pivot column. This is done to improve the
       * numerical stability of the algorithm. */
      if (with_pivoting)
      {
        std::swap(P[pivot], P[argmax]);
        A.swap_row(pivot, argmax);
      }

      /* Normalize the current pivot equation. */
      if (with_normalization)
      {
        b[pivot] /= A[pivot][pivot];
        for (size_t j = pivot + 1; j < n; ++j)
          A[pivot][j] /= A[pivot][pivot];
        A[pivot][pivot] = 1.0;
      }

      /* Perform row-wise operations such that all sub-diagonal values pivot
       * column are zero. This is done by subtracting the pivot row times the
       * ratio of the current row's pivot column value from the current row. */
      for (size_t i = pivot + 1; i < n; ++i)
      {
        double factor = A[i][pivot] / A[pivot][pivot];

        for (size_t j = pivot; j < n; ++j)
          A[i][j] -= A[pivot][j] * factor;
        b[i] -= b[pivot] * factor;
      }
      ++pivot;
    }
  }
}

//######################################################################

Vector math::GaussianElimination::back_substitution()
{
  size_t n = A.n_rows();
  Vector x(n, 0.0);

  /* The row-echelon system is upper triangular. The bottom equation is
   * not coupled to any other unknowns, so it is solved first. All subsequent
   * equations depend on those prior, and therefore, can be directly computed
   * with previous results. */
  for (int i = n - 1; i >= 0; --i)
  {
    double value = b[i];
    for (int j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }

  /* Map a pivoted solution back to the original ordering. */
  if (with_pivoting)
    for (size_t i = 0; i < n; ++i)
      x[P[i]] = x[i];

  return x;
}

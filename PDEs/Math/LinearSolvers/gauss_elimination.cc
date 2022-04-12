#include "gauss_elimination.h"


//######################################################################
void GaussElimination::row_echelon()
{
  size_t n = A.n_rows();
  size_t pivot = 0;

  // While the number of pivots is less than the dimension of the system...
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
      if (with_pivoting) A.swap_row(pivot, argmax);

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
Vector GaussElimination::backward_substitution_solve()
{
  size_t n = A.n_rows();
  Vector x(n, 0.0);

  for (int i = n - 1; i >= 0; --i)
  {
    double value = b[i];
    for (int j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}

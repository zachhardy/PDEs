#include "sparse_lu.h"
#include "macros.h"

#include <cmath>

using namespace pdes::Math;

//################################################## Constructors

SparseLU::SparseLU(const SparseMatrix& other, const bool pivot)
  : row_pivots(other.n_rows()), pivot_flag(pivot),
    SparseMatrix(other)
{}


SparseLU::SparseLU(SparseMatrix&& other, const bool pivot)
  : row_pivots(other.n_rows()), pivot_flag(pivot),
    SparseMatrix(other)
{}

//################################################## Properties

void
SparseLU::pivot(const bool flag)
{ pivot_flag = flag; }


bool
SparseLU::pivot() const
{ return pivot_flag; }

//################################################## Methods

void
SparseLU::factorize()
{
  size_t n = n_rows();

  // Initialize the pivot mappings such that each row maps to itself
  for (size_t i = 0; i < n; ++i)
    row_pivots[i] = i;

  //======================================== Apply Doolittle algorithm
  for (size_t j = 0; j < n; ++j)
  {
    /* Find the row index for the largest magnitude entry in this column.
     * This is only done for sub-diagonal elements. */
    if (pivot_flag)
    {
      const value_type* a_jj = diagonal(j);

      size_t argmax = j;
      value_type max = std::fabs((a_jj) ? *a_jj : 0.0);
      for (size_t k = j + 1; k < n; ++k)
      {
        const value_type* a_kj = locate(k, j);
        if (a_kj && *a_kj > max)
        {
          argmax = k;
          max = std::fabs(*a_kj);
        }
      }

      // If the sub-diagonal is uniformly zero, throw an error.
      Assert(max != 0.0, "Singular matrix error.");

      /* Swap the current row and the row containing the largest magnitude
       * entry corresponding for the current column. This is done to improve
       * the numerical stability of the algorithm. */
      if (argmax != j)
      {
        std::swap(row_pivots[j], row_pivots[argmax]);
        swap_row(j, argmax);
      }
    }//if pivot

    const value_type a_jj = *diagonal(j);

    // Compute the elements of the LU decomposition
    for (size_t i = j + 1; i < n; ++i)
    {
      value_type* a_ij = locate(i, j);
      if (a_ij && *a_ij != 0.0)
      {
        /* Lower triangular components. This represents the row operations
         * performed to attain the upper-triangular, row-echelon matrix. */
        *a_ij /= a_jj;

        /* Upper triangular components. This represents the row-echelon form
         * of the original matrix. Her*/
        for (const auto el : const_row_iterator(j))
          if (el.column > j)
            add(i, el.column, -(*a_ij) * el.value);
      }//if a_ij exists
    }//for rows > j
  }//for j
  factorized = true;
}


void
SparseLU::solve(const Vector& b, Vector& x) const
{
  Assert(factorized, "The matrix must be factorized before solve is called.");
  Assert(b.size() == n_rows(), "Dimension mismatch error.");
  Assert(x.size() == n_cols(), "Dimension mismatch error.");

  size_t n = n_rows();

  //================================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    value_type value = b[row_pivots[i]];
    for (const auto el: const_row_iterator(i))
      if (el.column < i)
        value -= el.value * x[el.column];
    x[i] = value;
  }

  //================================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    value_type value = x[i];
    for (const auto el : const_row_iterator(i))
      if (el.column > i)
        value -= el.value * x[el.column];
    x[i] = value / *diagonal(i);
  }
}


Vector
SparseLU::solve(const Vector& b) const
{
  Vector x(n_cols());
  solve(b, x);
  return x;
}

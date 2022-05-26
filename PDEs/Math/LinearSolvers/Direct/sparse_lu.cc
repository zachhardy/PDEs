#include "sparse_lu.h"

#include "vector.h"
#include "sparse_matrix.h"

#include "macros.h"

#include <cmath>


using namespace pdes::Math;

//################################################## Constructors

LinearSolver::SparseLU::SparseLU(SparseMatrix& other, const bool pivot)
  : A(other), row_pivots(other.n_rows()), pivot_flag(pivot)
{
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  factorize();
}


//################################################## Properties

void
LinearSolver::SparseLU::pivot(const bool flag)
{ pivot_flag = flag; }


bool
LinearSolver::SparseLU::pivot() const
{ return pivot_flag; }

//################################################## Methods

void
LinearSolver::SparseLU::factorize()
{
  size_t n = A.n_rows();

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
      const double* a_jj = A.diagonal(j);

      size_t argmax = j;
      double max = std::fabs((a_jj) ? *a_jj : 0.0);
      for (size_t k = j + 1; k < n; ++k)
      {
        const double* a_kj = A.locate(k, j);
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
        A.swap_row(j, argmax);
      }
    }//if pivot

    const double a_jj = *A.diagonal(j);

    // Compute the elements of the LU decomposition
    for (size_t i = j + 1; i < n; ++i)
    {
      double* a_ij = A.locate(i, j);
      if (a_ij && *a_ij != 0.0)
      {
        /* Lower triangular components. This represents the row operations
         * performed to attain the upper-triangular, row-echelon matrix. */
        *a_ij /= a_jj;

        /* Upper triangular components. This represents the row-echelon form
         * of the original matrix. Her*/
        for (const auto el : A.const_row_iterator(j))
          if (el.column > j)
            A.add(i, el.column, -(*a_ij) * el.value);
      }//if a_ij exists
    }//for rows > j
  }//for j
  factorized = true;
}


void
LinearSolver::SparseLU::solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  Assert(factorized, "The matrix must be factorized before solve is called.");
  Assert(b.size() == n, "Dimension mismatch error.");
  Assert(x.size() == n, "Dimension mismatch error.");

  //================================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    double value = b[row_pivots[i]];
    for (const auto el : A.const_row_iterator(i))
      if (el.column < i)
        value -= el.value * x[el.column];
    x[i] = value;
  }

  //================================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    double value = x[i];
    for (const auto el : A.const_row_iterator(i))
      if (el.column > i)
        value -= el.value * x[el.column];
    x[i] = value / *A.diagonal(i);
  }
}

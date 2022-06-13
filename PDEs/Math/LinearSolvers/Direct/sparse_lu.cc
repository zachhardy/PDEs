#include "sparse_lu.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include <cmath>
#include <cassert>


using namespace Math::LinearSolver;


//################################################## Setup


SparseLU::SparseLU(const bool pivot) : pivot_flag(pivot)
{}


void
SparseLU::set_matrix(const SparseMatrix& matrix)
{
  row_pivots.resize(matrix.n_rows());
  DirectSolverBase<SparseMatrix>::set_matrix(matrix);
}


//################################################## Methods


void
SparseLU::factorize()
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
      const double a_jj = A.diag_el(j);

      size_t argmax = j;
      double max = std::fabs((a_jj)? a_jj : 0.0);
      for (size_t k = j + 1; k < n; ++k)
      {
        const double a_kj = A.el(k, j);
        if (std::fabs(a_kj) > max)
        {
          argmax = k;
          max = std::fabs(a_kj);
        }
      }

      // If the sub-diagonal is uniformly zero, throw an error.
      assert(max != 0.0);

      /* Swap the current row and the row containing the largest magnitude
       * entry corresponding for the current column. This is done to improve
       * the numerical stability of the algorithm. */
      if (argmax != j)
      {
        std::cout << "Swapping row " << j
                  << " with row " << argmax << std::endl;
        std::cout << "Pre-Swap:\n";
        A.print_formatted();

        std::swap(row_pivots[j], row_pivots[argmax]);
        A.swap_row(j, argmax);

        std::cout << "Post-Swap:\n";
        A.print_formatted();
      }
    }//if pivot

    const double a_jj = A.diag(j);

    // Compute the elements of the LU decomposition
    for (size_t i = j + 1; i < n; ++i)
    {
      if (A.exists(i, j))
      {
        double& a_ij = A(i, j);

        // Get the
        /* Lower triangular components. This represents the row operations
         * performed to attain the upper-triangular, row-echelon matrix. */
        a_ij /= a_jj;

        /* Upper triangular components. This represents the row-echelon form
         * of the original matrix. */
        for (const auto el : A.row_iterator(j))
          if (el.column() > j)
            A.add(i, el.column(), -a_ij*el.value());
      }//if a_ij exists
    }//for rows > j
  }//for j
  factorized = true;
}


void
SparseLU::solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  assert(factorized);
  assert(b.size() == n);
  assert(x.size() == n);

  //======================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    double value = b[row_pivots[i]];
    for (const auto el : A.row_iterator(i))
      if (el.column() < i)
        value -= el.value()*x[el.column()];
    x[i] = value;
  }

  //======================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    double value = x[i];
    for (const auto el : A.row_iterator(i))
      if (el.column() > i)
        value -= el.value()*x[el.column()];
    x[i] = value/A.diag(i);
  }
}


//################################################## Properties


void
SparseLU::pivot(const bool flag)
{ pivot_flag = flag; }


bool
SparseLU::pivot() const
{ return pivot_flag; }

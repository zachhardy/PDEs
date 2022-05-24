#include "sparse_cholesky.h"
#include "macros.h"

#include <cmath>

using namespace pdes::Math;

//################################################## Constructors

SparseCholesky::SparseCholesky(const SparseMatrix& other)
  : SparseMatrix(other)
{}


SparseCholesky::SparseCholesky(SparseMatrix&& other)
  : SparseCholesky(other)
{}

//################################################## Methods

void
SparseCholesky::factorize()
{
  size_t n = n_rows();

  // Compute the factorization column by column
  for (size_t j = 0; j < n; ++j)
  {
    // Accessor for the diagonal element
    value_type* d = locate(j, j);
    Assert(d && *d != 0.0, "Singular matrix error.");

    // Compute the new diagonal term
    value_type sum = 0.0;
    for (const auto elem : const_row_iterator(j))
      if (elem.column < j)
        sum += elem.value * elem.value;
    *d = std::sqrt(*d - sum);

    // Set the lower-diagonal components
    for (size_t i = j + 1; i < n; ++i)
    {
      // Go through row i and j, add to sum when columns are equal
      sum = 0.0;
      for (const auto a_ik : const_row_iterator(i))
        if (a_ik.column < j)
          for (const auto a_jk : const_row_iterator(j))
            if (a_jk.column == a_ik.column)
              sum += a_ik.value * a_jk.value;

      // Set element i, j
      value_type* a_ij = locate(i, j);
      value_type value = (a_ij) ? (*a_ij - sum) / *d : -sum/ *d;
      if (std::fabs(value) != 0.0)
        set(i, j, value);
    }
  }
  factorized = true;
}


void
SparseCholesky::solve(const Vector& b, Vector& x) const
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
    for (const auto el : const_row_iterator(i))
      if (el.column < i)
        value -= el.value * x[el.column];
    x[i] = value / (*this)(i, i);
  }

  //================================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    x[i] /= *diagonal(i);
    for (const auto a_ij : const_row_iterator(i))
      if (a_ij.column < i)
        x[a_ij.column] -= a_ij.value * x[i];
  }
}


Vector
SparseCholesky::solve(const Vector& b) const
{
  Vector x(n_cols());
  solve(b, x);
  return x;
}

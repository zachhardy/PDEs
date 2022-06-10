#include "sparse_cholesky.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include <cmath>
#include <cassert>


using namespace Math::LinearSolver;


SparseCholesky::SparseCholesky() :
  DirectSolverBase<SparseMatrix>()
{}


void
SparseCholesky::factorize()
{
  size_t n = A.n_rows();

  // Compute the factorization column by column
  for (size_t j = 0; j < n; ++j)
  {
    // Accessor for the diagonal element
    double& d = A.diag(j);

    // Compute the new diagonal term
    double sum = 0.0;
    for (const auto el : A.row_iterator(j))
      if (el.column() < j)
        sum += el.value() * el.value();
    d = std::sqrt(d - sum);

    // Set the lower-diagonal components
    for (size_t i = j + 1; i < n; ++i)
    {
      // Go through row i and j, add to sum when columns are equal
      sum = 0.0;
      for (const auto a_ik : A.row_iterator(i))
        if (a_ik.column() < j)
          for (const auto a_jk : A.row_iterator(j))
            if (a_jk.column()== a_ik.column())
              sum += a_ik.value() * a_jk.value();

      // Set element i, j
      double a_ij = A.el(i, j);
      double value = (a_ij)? (a_ij - sum) / d : -sum / d;
      if (std::fabs(value) != 0.0)
        A.set(i, j, value);
    }
  }
  factorized = true;
}


void
SparseCholesky::solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  assert(factorized);
  assert(b.size() == n);
  assert(x.size() == n);

  //======================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    double value = b[i];
    for (const auto el : A.row_iterator(i))
      if (el.column() < i)
        value -= el.value() * x[el.column()];
    x[i] = value / A.diag(i);
  }

  //======================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    x[i] /= A.diag(i);
    for (const auto a_ij : A.row_iterator(i))
      if (a_ij.column() < i)
        x[a_ij.column()] -= a_ij.value() * x[i];
  }
}

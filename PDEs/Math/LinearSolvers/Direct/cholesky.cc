#include "cholesky.h"

#include "vector.h"
#include "matrix.h"

#include <cmath>
#include <cassert>


using namespace PDEs;
using namespace Math;
using namespace LinearSolvers;


Cholesky::Cholesky() :
    DirectSolverBase<Matrix>()
{}


SparseCholesky::SparseCholesky() :
    DirectSolverBase<SparseMatrix>()
{}


void
Cholesky::factorize()
{
  if (factorized) return;

  // Compute the factorization column by column
  size_t n = A.n_rows();
  for (size_t j = 0; j < n; ++j)
  {
    double* a_j = A.data(j); // accessor for row j

    // Compute the diagonal term
    double sum = 0.0;
    for (size_t k = 0; k < j; ++k)
      sum += a_j[k] * a_j[k];
    a_j[j] = std::sqrt(a_j[j] - sum);

    // Compute the lower diagonal terms
    for (size_t i = j + 1; i < n; ++i)
    {
      double* a_i = A.data(i); // accessor for row i

      sum = 0.0;
      for (size_t k = 0; k < j; ++k)
        sum += a_i[k] * a_j[k];
      a_i[j] = (a_i[j] - sum) / a_j[j];
    }
  }
  factorized = true;
}


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
    for (const auto el: A.row_iterator(j))
      if (el.column < j)
        sum += el.value * el.value;
    d = std::sqrt(d - sum);

    // Set the lower-diagonal components
    for (size_t i = j + 1; i < n; ++i)
    {
      // Go through row i and j, add to sum when columns are equal
      sum = 0.0;
      for (const auto a_ik: A.row_iterator(i))
        if (a_ik.column < j)
          for (const auto a_jk: A.row_iterator(j))
            if (a_jk.column == a_ik.column)
              sum += a_ik.value * a_jk.value;

      // Set element i, j
      double a_ij = A.el(i, j);
      double value = (a_ij) ? (a_ij - sum) / d : -sum / d;
      if (std::fabs(value) != 0.0)
        A.set(i, j, value);
    }
  }
  factorized = true;
}


void
Cholesky::solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  assert(factorized);
  assert(b.size() == n);
  assert(x.size() == n);

  //======================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    const double* a_i = A.data(i); // accessor for row i
    const double a_ii = a_i[i]; // diagonal element for row i

    double value = b[i];
    for (size_t j = 0; j < i; ++j)
      value -= *a_i++ * x[j];
    x[i] = value / a_ii;
  }

  //======================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    const double* a_i = A.data(i); // accessor for row i
    const double a_ii = a_i[i]; // diagonal element for row i

    double& x_i = x[i]; // accessor for element i

    x_i /= a_ii;
    for (size_t j = 0; j < i; ++j)
      x[j] -= *a_i++ * x_i;
  }
}


void
SparseCholesky::solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  assert(factorized);
  assert(b.size() == n);
  assert(x.size() == n);

  // Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    double value = b[i];
    for (const auto el: A.row_iterator(i))
      if (el.column < i)
        value -= el.value * x[el.column];
    x[i] = value / A.diag(i);
  }

  // Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    x[i] /= A.diag(i);
    for (const auto a_ij: A.row_iterator(i))
      if (a_ij.column < i)
        x[a_ij.column] -= a_ij.value * x[i];
  }
}


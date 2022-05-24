#include "cholesky.h"
#include <cmath>

using namespace pdes::Math;

//################################################## Constructors


Cholesky::Cholesky(const Matrix& other) : Matrix(other) {}


Cholesky::Cholesky(Matrix&& other) : Matrix(other) {}

//################################################## Methods

void
Cholesky::factorize()
{
  size_t n = n_rows();

  // Compute the factorization column by column
  for (size_t j = 0; j < n; ++j)
  {
    value_type* a_j = data(j); // accessor for row j

    // Compute the diagonal term
    value_type sum = 0.0;
    for (size_t k = 0; k < j; ++k)
      sum += a_j[k] * a_j[k];
    a_j[j] = std::sqrt(a_j[j] - sum);

    // Compute the lower diagonal terms
    for (size_t i = j + 1; i < n; ++i)
    {
      double* a_i = data(i); // accessor for row i

      sum = 0.0;
      for (size_t k = 0; k < j; ++k)
        sum += a_i[k] * a_j[k];
      a_i[j] = (a_i[j] - sum) / a_j[j];
    }
  }
  factorized = true;
}


void
Cholesky::solve(const Vector& b, Vector& x) const
{
  Assert(factorized, "Matrix must be factorized before solving.");
  Assert(b.size() == n_rows(), "Dimension mismatch error.");
  Assert(x.size() == n_cols(), "Dimension mismatch error.");

  size_t n = n_rows();

  //================================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    const value_type* a_i = data(i); // accessor for row i
    const value_type a_ii = a_i[i]; // diagonal element for row i

    value_type value = b[i];
    for (size_t j = 0; j < i; ++j)
      value -= *a_i++ * x[j];
    x[i] = value / a_ii;
  }

  //================================================== Backward solve
  for (size_t i = n - 1; i != -1; --i)
  {
    const value_type* a_i = data(i); // accessor for row i
    const value_type a_ii = a_i[i]; // diagonal element for row i

    value_type& x_i = x[i]; // accessor for element i

    x_i /= a_ii;
    for (size_t j = 0; j < i; ++j)
      x[j] -= *a_i++ * x_i;
  }
}


Vector
Cholesky::solve(const Vector& b) const
{
  Vector x(n_cols());
  solve(b, x);
  return x;
}

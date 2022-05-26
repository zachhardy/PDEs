#include "cholesky.h"

#include "vector.h"
#include "matrix.h"

#include "macros.h"

#include <cmath>


using namespace pdes::Math;

//################################################## Constructors


LinearSolver::Cholesky::
Cholesky(Matrix& A) : A(A)
{
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  factorize();
}


//################################################## Methods

void
LinearSolver::Cholesky::factorize()
{
  size_t n = A.n_rows();

  // Compute the factorization column by column
  for (size_t j = 0; j < n; ++j)
  {
    double * a_j = A.data(j); // accessor for row j

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
LinearSolver::Cholesky::
solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  Assert(factorized, "Matrix must be factorized before solving.");
  Assert(b.size() == n, "Dimension mismatch error.");
  Assert(x.size() == n, "Dimension mismatch error.");

  //================================================== Forward solve
  for (size_t i = 0; i < n; ++i)
  {
    const double* a_i = A.data(i); // accessor for row i
    const double a_ii = a_i[i]; // diagonal element for row i

    double value = b[i];
    for (size_t j = 0; j < i; ++j)
      value -= *a_i++ * x[j];
    x[i] = value / a_ii;
  }

  //================================================== Backward solve
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

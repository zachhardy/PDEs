#include "gauss_seidel.h"

#include <sstream>
#include <iomanip>
#include <cmath>

using namespace pdes::Math;

GaussSeidelSolver::
GaussSeidelSolver(const double tolerance,
                  const size_t max_iterations)
  : tol(tolerance), maxiter(max_iterations)
{}


void
GaussSeidelSolver::solve(const SparseMatrix& A,
                         Vector& x, const Vector& b) const
{
  Assert(A.n_rows() == A.n_cols(), "Only square matrices are allowed.")
  Assert(A.n_rows() == b.size(), "Dimension mismatch error.");
  Assert(A.n_cols() == x.size(), "Dimension mismatrch error.");

  size_t n = A.n_rows();

  //======================================== Iteration loop
  double diff; size_t nit; bool converged = false;
  for (nit = 0; nit < maxiter; ++nit)
  {
    diff = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
      double value = b[i];
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          value -= el.value * x[el.column];
      value /= *A.diagonal(i);

      diff += std::fabs(value - x[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    if (diff < tol)
    { converged = true; break;}
  }
  Assert(converged, "Gauss Seidel solver did not converge.");
}


Vector
GaussSeidelSolver::solve(const SparseMatrix& A, const Vector& b) const
{
  Vector x(A.n_cols(), 0.0);
  solve(A, x, b);
  return x;
}

#include "gauss_seidel.h"

#include <sstream>
#include <iomanip>
#include <cmath>

using namespace pdes::Math;

GaussSeidelSolver::GaussSeidelSolver(const double tolerance,
                                     const size_t max_iterations)
  : tol(tolerance), max_iter(max_iterations)
{}


void
GaussSeidelSolver::solve(const SparseMatrix& A,
                         const Vector& b, Vector& x)
{
  Assert(A.n_rows() == A.n_cols(), "Only square matrices are allowed.")
  Assert(A.n_rows() == b.size(), "Dimension mismatch error.");
  Assert(A.n_cols() == x.size(), "Dimension mismatrch error.");

  size_t n = A.n_rows();

  // Book-keeping parameters
  bool converged = false;

  // Iterate
  double diff; size_t k;
  for (k = 0; k < max_iter; ++k)
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
    {
      converged = true;
      break;
    }
  }
}

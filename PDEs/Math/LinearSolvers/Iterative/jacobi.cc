#include "jacobi.h"

#include <sstream>
#include <iomanip>
#include <cmath>

using namespace pdes::Math;


JacobiSolver::JacobiSolver(const double tolerance,
                           const size_t max_iterations)
  : tol(tolerance), max_iter(max_iterations)
{
  Assert(tol > 0.0, "Illegal negative tolerance specified.");
}


void JacobiSolver::solve(const SparseMatrix& A,
                         const Vector& b, Vector& x)
{
  Assert(A.n_rows() == A.n_cols(), "Only square matrices are allowed.")
  Assert(A.n_rows() == b.size(), "Dimension mismatch error.");
  Assert(A.n_cols() == x.size(), "Dimension mismatrch error.");

  size_t n = A.n_rows();

  // Book-keeping parameters
  Vector x_ell = x;
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
          value -= el.value * x_ell[el.column];
      value /= *A.diagonal(i);

      diff += std::fabs(value - x_ell[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    x_ell = x;
    if (diff < tol)
    {
      converged = true;
      break;
    }
  }

//  std::stringstream ss;
//  ss << "***** JacobiSolver ";
//  if (converged)
//    ss << "converged in " << k << " iterations";
//  else
//    ss << "did not converge. Final difference = " << diff;
//  std::cout << ss.str() << ".\n";

}


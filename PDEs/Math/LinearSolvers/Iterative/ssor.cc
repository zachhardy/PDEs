#include "ssor.h"

#include "vector.h"
#include "sparse_matrix.h"

#include "macros.h"

#include <cmath>
#include <iomanip>


using namespace pdes::Math;


LinearSolver::SSOR::
SSOR(const SparseMatrix& A, const Options& opts) :
  SOR(A, opts, "SSOR")
{
}


void LinearSolver::SSOR::
solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  Assert(b.size() == n, "Dimension mismatch error.");
  Assert(x.size() == n, "Dimension mismatrch error.");

  size_t nit;
  double change;
  Vector x_ell = x;

  //======================================== Iteration loop
  for (nit = 0; nit < max_iterations; ++nit)
  {
    change = 0.0;

    //==================== Compute forward sweep
    for (size_t i = 0; i < n; ++i)
    {
      double s = 0.0;
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          s += el.value * x[el.column];

      double a_ii = *A.diagonal(i);
      x[i] += omega * ((b[i] - s)/a_ii - x[i]);
    }

    //==================== Compute backward sweep
    for (size_t i = n - 1; i != -1; --i)
    {
      double s = 0.0;
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          s += el.value * x[el.column];

      double a_ii = *A.diagonal(i);
      x[i] +=  omega * ((b[i] - s)/a_ii - x[i]);
      change += std::fabs(x[i] - x_ell[i]) / std::fabs(x[i]);
    }

    //==================== Check convergence
    x_ell = x;
    if (check(nit, change))
      break;
  }
}

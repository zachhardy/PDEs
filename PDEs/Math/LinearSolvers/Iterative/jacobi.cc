#include "jacobi.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include "macros.h"

#include <cmath>


using namespace pdes::Math;


LinearSolver::Jacobi::
Jacobi(const SparseMatrix& A, const Options& opts) :
    IterativeSolverBase(A, opts, "Jacobi")
{}


void
LinearSolver::Jacobi::
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
    for (size_t i = 0; i < n; ++i)
    {
      //==================== Compute element-wise update
      double value = b[i];
      for (const auto el : A.const_row(i))
        if (el.column != i)
          value -= el.value * x_ell[el.column];
      value /= *A.diagonal(i);

      //==================== Increment difference
      change += std::fabs(value - x_ell[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    //==================== Check convergence
    x_ell = x;
    if (check(nit + 1, change))
      break;
  }
}



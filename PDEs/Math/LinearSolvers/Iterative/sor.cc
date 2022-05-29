#include "sor.h"

#include "vector.h"
#include "sparse_matrix.h"

#include "macros.h"

#include <cmath>


using namespace pdes::Math;

LinearSolver::SOR::
SOR(const SparseMatrix& A,
    const Options& opts,
    const std::string solver_name) :
    IterativeSolverBase(A, opts, solver_name),
    omega(opts.omega)
{
  Assert(omega > 0 && omega < 2, "Invalid relaxation parameter.");
}


void LinearSolver::SOR::
solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  Assert(b.size() == n, "Dimension mismatch error.");
  Assert(x.size() == n, "Dimension mismatrch error.");

  size_t nit;
  double change;

  //======================================== Iteration loop
  for (nit = 0; nit < max_iterations; ++nit)
  {
    change = 0.0;
    for (size_t i = 0; i < A.n_rows(); ++i)
    {
      //==================== Compute element-wise update
      double value = 0.0;
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          value += el.value * x[el.column];

      double a_ii = *A.diagonal(i);
      value = x[i] + omega * ((b[i] - value)/a_ii - x[i]);

      //==================== Increment difference
      change += std::fabs(value - x[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    //==================== Check convergence
    if (check(nit, change))
      break;
  }
}

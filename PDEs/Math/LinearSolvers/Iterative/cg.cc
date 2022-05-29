#include "cg.h"

#include "vector.h"
#include "Math/Sparse/sparse_matrix.h"

#include "macros.h"

#include <cmath>

using namespace pdes::Math;


LinearSolver::CG::
CG(const SparseMatrix& A, const Options& opts) :
    IterativeSolverBase(A, opts, "CG")
{}


void
LinearSolver::CG::
solve(Vector& x, const Vector& b) const
{
  size_t n = A.n_rows();
  Assert(b.size() == n, "Dimension mismatch error.");
  Assert(x.size() == n, "Dimension mismatrch error.");

  // Allocate data needed for the CG solver
  Vector r(x.size());
  Vector p(x.size());
  Vector q(x.size());

  double alpha;
  double res;
  double res_prev;

  /* Initialize residual, residual norms, and search directions.
   * If the residual norm is smaller than the tolerance, exit because
   * the initial guess is the solution. */
  r = (!x.all_zero())? b - A.vmult(x) : b;
  res = res_prev = r.dot(r);
  if (res < tolerance)
    return;
  p = r;

  //======================================== Iteration loop
  size_t nit;
  for (nit = 0; nit < max_iterations; ++nit)
  {
    // Precompute necessary matrix-vector product q = Ap
    A.vmult(p, q);

    // Recompute alpha factor
    alpha = res_prev / p.dot(q);

    // Update solution and residual vector
    x.add(p, alpha);
    r.add(q, -alpha);

    // Update residual norm
    res = r.dot(r);

    // Check convergence
    if(check(nit + 1, std::sqrt(res)))
      break;

    // If not converged, prep for next iteration
    p.sadd(res/res_prev, r);
    res_prev = res;
  }
}

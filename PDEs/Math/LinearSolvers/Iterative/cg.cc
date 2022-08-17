#include "cg.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include <cmath>
#include <cassert>


using namespace Math;


LinearSolver::CG::
CG(const Options& opts) : IterativeSolverBase(opts, "CG")
{}


void
LinearSolver::CG::
solve(Vector& x, const Vector& b) const
{
  size_t n = A->n_rows();
  assert(b.size() == n);
  assert(x.size() == n);

  double norm = b.l2_norm();

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
  r = (!x.all_zero())? b - A->vmult(x) : b;
  res = res_prev = r.dot(r);
  if (res < tolerance)
    return;
  p = r;

  //======================================== Iteration loop
  size_t nit;
  for (nit = 0; nit < max_iterations; ++nit)
  {
    // Precompute necessary matrix-vector product q = Ap
    A->vmult(p, q);

    // Recompute alpha factor
    alpha = res_prev/p.dot(q);

    // Update solution and residual vector
    x.add(p, alpha);

    r.add(q, -alpha);

    // Update residual norm
    res = r.dot(r);

    // Check convergence
    bool converged = check(nit + 1, std::sqrt(res)/norm);
    if (converged)
      break;

    // If not converged, prep for next iteration
    p.sadd(res/res_prev, r);
    res_prev = res;
  }
}

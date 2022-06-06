#include "cg.h"

#include "vector.h"
#include "Sparse/sparse_matrix.h"

#include "macros.h"
#include "Timer/timer.h"

#include <cmath>


using namespace Math;


LinearSolver::CG::
CG(const Options& opts) : IterativeSolverBase(opts, "CG")
{}


void
LinearSolver::CG::
solve(Vector& x, const Vector& b) const
{
  size_t n = A->n_rows();
  Assert(b.size() == n, "Dimension mismatch error.")
  Assert(x.size() == n, "Dimension mismatch error.")

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
  Timer timer;
  std::vector<double> times;
  for (nit = 0; nit < max_iterations; ++nit)
  {
    times.clear();

    // Precompute necessary matrix-vector product q = Ap
    timer.start();
    A->vmult(p, q);
    timer.stop();
    times.push_back(timer.get_time());

    // Recompute alpha factor
    timer.start();
    alpha = res_prev / p.dot(q);
    timer.stop();
    times.push_back(timer.get_time());

    // Update solution and residual vector
    timer.start();
    x.add(p, alpha);
    timer.stop();
    times.push_back(timer.get_time());

    timer.start();
    r.add(q, -alpha);
    timer.stop();
    times.push_back(timer.get_time());

    // Update residual norm
    timer.start();
    res = r.dot(r);
    timer.stop();
    times.push_back(timer.get_time());

    // Check convergence
    timer.start();
    bool converged = check(nit + 1, std::sqrt(res));
    timer.stop();
    times.push_back(timer.get_time());

    if (converged)
      break;

    // If not converged, prep for next iteration
    timer.start();
    p.sadd(res / res_prev, r);
    timer.stop();
    times.push_back(timer.get_time());

    timer.start();
    res_prev = res;
    timer.stop();
    times.push_back(timer.get_time());

    double total_time = 0.0;
    for (const auto& time : times)
      total_time += time;
  }
}

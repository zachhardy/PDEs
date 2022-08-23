#include "sor.h"

#include "vector.h"
#include "Math/sparse_matrix.h"

#include <cmath>
#include <cassert>


using namespace PDEs;
using namespace Math;
using namespace LinearSolvers;


SOR::SOR(const double omega,
         const Options& opts,
         const std::string solver_name) :
    IterativeSolverBase(opts, solver_name), omega(omega)
{ assert(omega > 0 && omega < 2); }


GaussSeidel::GaussSeidel(const Options& opts) :
    SOR(1.0, opts, "Gauss-Seidel")
{}


SSOR::SSOR(const double omega, const Options& opts) :
    SOR(omega, opts, "SSOR")
{}


void
SOR::solve(Vector& x, const Vector& b) const
{
  size_t n = A->n_rows();
  assert(b.size() == n);
  assert(x.size() == n);

  size_t nit;
  double change;

  // Iteration loop
  for (nit = 0; nit < max_iterations; ++nit)
  {
    change = 0.0;
    for (size_t i = 0; i < A->n_rows(); ++i)
    {
      // Compute element-wise update
      double value = 0.0;
      for (const auto el: A->row_iterator(i))
        if (el.column != i)
          value += el.value * x[el.column];

      double a_ii = A->diag(i);
      value = x[i] + omega * ((b[i] - value) / a_ii - x[i]);

      // Increment difference
      change += std::fabs(value - x[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    // Check convergence
    if (check(nit + 1, change))
      break;
  }
}


void
SSOR::solve(Vector& x, const Vector& b) const
{
  size_t n = A->n_rows();
  assert(b.size() == n);
  assert(x.size() == n);

  size_t nit;
  double change;
  Vector x_ell = x;

  // Iteration loop
  for (nit = 0; nit < max_iterations; ++nit)
  {
    change = 0.0;

    // Compute forward sweep
    for (size_t i = 0; i < n; ++i)
    {
      double s = 0.0;
      for (const auto el: A->row_iterator(i))
        if (el.column != i)
          s += el.value * x[el.column];

      double a_ii = A->diag(i);
      x[i] += omega * ((b[i] - s) / a_ii - x[i]);
    }

    // Compute backward sweep
    for (size_t i = n - 1; i != -1; --i)
    {
      double s = 0.0;
      for (const auto el: A->row_iterator(i))
        if (el.column != i)
          s += el.value * x[el.column];

      double a_ii = A->diag(i);
      x[i] += omega * ((b[i] - s) / a_ii - x[i]);
      change += std::fabs(x[i] - x_ell[i]) / std::fabs(b[i]);
    }

    // Check convergence
    x_ell = x;
    if (check(nit + 1, change))
      break;
  }
}

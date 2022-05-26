#include "jacobi.h"
#include "macros.h"

#include <sstream>
#include <iomanip>
#include <cmath>

using namespace pdes::Math;


JacobiSolver::
JacobiSolver(const SparseMatrix& A,
             const double tolerance,
             const size_t max_iterations)
  : A(A), tol(tolerance), maxiter(max_iterations)
{
  Assert(A.n_rows() == A.n_cols(), "Square matrix required.");
  Assert(tol > 0.0, "Illegal negative tolerance specified.");
}


void
JacobiSolver::solve(const Vector& b, Vector& x) const
{
  Assert(A.n_rows() == b.size(), "Dimension mismatch error.");
  Assert(A.n_cols() == x.size(), "Dimension mismatrch error.");

  size_t n = A.n_rows();

  //======================================== Iteration loop
  Vector x_ell = x;
  value_type diff; size_t nit; bool converged = false;
  for (nit = 0; nit < maxiter; ++nit)
  {
    diff = 0.0;

    //============================== Loop over equations
    for (size_t i = 0; i < n; ++i)
    {
      //==================== Compute solution update
      value_type value = b[i];
      for (const auto el : A.const_row_iterator(i))
        if (el.column != i)
          value -= el.value * x_ell[el.column];
      value /= *A.diagonal(i);

      //==================== Increment solution change
      diff += std::fabs(value - x_ell[i]) / std::fabs(b[i]);
      x[i] = value;
    }

    //============================== Prep for next iteration
    x_ell = x;
    if (diff < tol)
    { converged = true; break;}
  }
  Assert(converged, "Jacobi solver did not converge.");
}


Vector
JacobiSolver::solve(const Vector& b) const
{
  Vector x(A.n_cols(), 0.0);
  solve(b, x);
  return x;
}


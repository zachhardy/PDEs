#include "LinearSolvers/linear_solver.h"

#include "vector.h"

#include "macros.h"


using namespace pdes::Math;


Vector
LinearSolver::LinearSolverBase::
solve(const Vector& b) const
{
  Vector x(b.size(), 0.0);
  solve(x, b);
  return x;
}


LinearSolver::IterativeSolverBase::
IterativeSolverBase(const SparseMatrix& A,
                    const double tolerance,
                    const size_t max_iterations,
                    const bool verbose) :
  A(A),
  tolerance(tolerance),
  max_iterations(max_iterations),
  verbose(verbose)
{
  Assert(tolerance > 0.0, "Tolerance must be positive.")
}


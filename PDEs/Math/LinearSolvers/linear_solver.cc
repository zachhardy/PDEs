#include "LinearSolvers/linear_solver.h"

#include "vector.h"

using namespace pdes::Math;

Vector
LinearSolver::LinearSolverBase::
solve(const Vector& b) const
{
  Vector x(b.size(), 0.0);
  solve(x, b);
  return x;
}

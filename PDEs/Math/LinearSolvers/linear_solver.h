#ifndef LINEAR_SOLVER_BASE_H
#define LINEAR_SOLVER_BASE_H

#include "vector.h"

namespace pdes::Math
{

enum class LinearSolverType
{
  DIRECT = 0,
  ITERATIVE = 1
};


class LinearSolverBase
{
  virtual void
  solve(const Vector& b, Vector& x) const = 0;

  virtual Vector
  solve(const Vector& b) const = 0;
};

}

#endif //LINEAR_SOLVER_BASE_H

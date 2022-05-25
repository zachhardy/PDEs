#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

namespace pdes::Math
{

/**
 * Types of linear solvers.
 */
enum class LinearSolverType
{
  DIRECT = 0,
  ITERATIVE = 1
};

}
#endif //LINEAR_SOLVER_H

#ifndef LINEAR_SOLVER_BASE_H
#define LINEAR_SOLVER_BASE_H


//########## Forward declarations
namespace pdes::Math
{
  class Vector;
  class Matrix;
  class SparseMatrix;
}

namespace pdes::Math::LinearSolver
{

/**
 * Available types of linear solvers.
 */
enum class LinearSolverType
{
  LU            = 0,
  CHOLESKY      = 1,
  JACOBI        = 2,
  GAUSS_SEIDEL  = 3,
  SOR           = 4
};


/**
 * Base class from which all linear solvers must derive.
 */
class LinearSolverBase
{
public:
  virtual void
  solve(Vector& x, const Vector& b) const = 0;

  Vector
  solve(const Vector& b) const;
};

}

#endif //LINEAR_SOLVER_BASE_H

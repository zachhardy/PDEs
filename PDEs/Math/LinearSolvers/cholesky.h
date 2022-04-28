#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "Math//math.h"
#include "linear_solver.h"

namespace math
{

/**
 * A Cholesky decomposition solver.
 * See \ref math::cholesky_factorization and \ref math::cholesky_solve for
 * implementation details.
 */
class Cholesky : public LinearSolver
{
public:
  Cholesky(Matrix<double>& matrix) : LinearSolver(matrix) {}

public:
  /**
   * Perform the Cholesky factorization.
   * See \ref math::cholesky_factorization for implementation details.
   */
  void setup() override { cholesky_factorization(A); initialized = true; }

  /**
   * Solve the Cholesky factored linear system.
   * See \ref math::cholesky_solve for implementation details.
   * \param The right-hand side vector of the linear system.
   * \return The solution to the linear system.
   */
  Vector<double> solve(const Vector<double>& b) override
  {
    if (not initialized) setup();
    return cholesky_solve(A, b);
  }
};

}

#endif //CHOLESKY_H

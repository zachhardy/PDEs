#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "Math//math.h"
#include "linear_solver.h"

namespace math
{

class Cholesky : public LinearSolver
{
private:
  bool initialized = false;

public:
  Cholesky(Matrix& matrix) : LinearSolver(matrix) {}

public:
  /**
   * \brief Perform the Cholesky factorization.
   *
   * See \ref math::cholesky_factorization for implementation details.
   */
  void setup() override { cholesky_factorization(A); initialized = true; }

  /**
   * \brief Solve the Cholesky factored linear system
   *
   * See \ref math::cholesky_solve for implementation details.
   *
   * \param The right-hand side vector of the linear system.
   * \return The solution to the linear system.
   */
  Vector solve(const Vector& b) override
  {
    if (not initialized) setup();
    return cholesky_solve(A, b);
  }
};

}

#endif //CHOLESKY_H

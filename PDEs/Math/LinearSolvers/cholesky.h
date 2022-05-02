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
template<typename value_type>
class Cholesky : public LinearSolver<value_type>
{
public:
  Cholesky(Matrix<value_type>& matrix) : LinearSolver<value_type>(matrix) {}

public:
  /**
   * Perform the Cholesky factorization.
   * See \ref math::cholesky_factorization for implementation details.
   */
  void setup() override;
//  { cholesky_factorization(A); initialized = true; }

  /**
   * Solve the Cholesky factored linear system.
   * See \ref math::cholesky_solve for implementation details.
   * \param The right-hand side vector of the linear system.
   * \return The solution to the linear system.
   */
  Vector<value_type> solve(const Vector<value_type>& b) override;
//  {
//    if (not initialized) setup();
//    return cholesky_solve(A, b);
//  }
};

}

#endif //CHOLESKY_H

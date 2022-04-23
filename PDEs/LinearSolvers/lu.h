#ifndef LU_H
#define LU_H

#include "Math//math.h"
#include "linear_solver.h"

namespace linear_solvers
{
using namespace math;

/**
 * \property linear_solvers::LU::row_pivots
 * \brief The pivot mapping vector.
 *
 * The index corresponds to the initial row number and the value to the pivoted
 * row number. This is used to map the right-hand side vector to the correct row
 * when solving.
 */

/**
 * \brief An LU decomposition solver.
 *
 * See \ref math::lu_factorization and \ref lu_solve for implementation details.
 */
class LU : public LinearSolver
{
private:
 bool initialized = false;
 bool pivot = true;

 std::vector<size_t> row_pivots;

public:
 LU(Matrix& matrix, const bool pivot_option = true)
   : LinearSolver(matrix), pivot(pivot_option)
 {}

public:
 void set_pivot_option(const bool pivot_option) { pivot = pivot_option; }
 bool get_pivot_option() const { return pivot; }

public:
 /**
  * \brief Perform the LU factorization.
  *
  * See \ref math::lu_factorization for implementation details.
  */
  void setup() override
  {
    row_pivots = lu_factorization(A, pivot);
    initialized = true;
  }

  /**
   * \brief Solve the LU factored linear system.
   *
   * See \ref math::lu_solve for implementation details.
   *
   * \param b The right-hand side vector of the linear system.
   * \return The solution to the linear system.
   */
  Vector solve(const Vector& b) override
  {
    if (not initialized) setup();
    return lu_solve(A, b, row_pivots);
  }
};

}
#endif //LU_H

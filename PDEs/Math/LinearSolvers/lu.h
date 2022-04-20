#ifndef LU_H
#define LU_H

#include "Math//math.h"
#include "linear_solver.h"

namespace math
{

class LU : public LinearSolver
{
private:
  /// A flag to determine if the solver has been initialized.
  bool initialized = false;

  /// An option for enabling row pivoting. Default is on.
  bool pivot = true;

  /** The pivot mapping vector. The index corresponds to the initial row number
   *  and the value to the pivoted row number. This is used to map the
   *  right-hand side vector to the correct row when solving. */
  std::vector<size_t> row_pivots;

public:
  LU(Matrix& matrix, const bool pivot_option = true)
    : LinearSolver(matrix), pivot(pivot_option)
  {}

public:
  /// Set the pivoting option.
  void set_pivot_option(const bool pivot_option) { pivot = pivot_option; }
  /// Get the pivoting option.
  bool get_pivot_option() const { return pivot; }

public:
  /**
   * \brief Perform the LU factorization and set initialized.
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
   * \param b The right-hand side vector of the linear system.
   * \return The solution to the linear system.
   *
   * See \ref math::lu_solve for implementation details.
   */
  Vector solve(const Vector& b) override
  {
    if (not initialized) setup();
    return lu_solve(A, b, row_pivots);
  }

public:

};

}


#endif //LU_H

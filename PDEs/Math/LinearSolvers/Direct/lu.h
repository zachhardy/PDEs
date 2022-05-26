#ifndef LU_H
#define LU_H

#include "LinearSolvers/linear_solver.h"

#include <vector>
#include <cstddef>


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of an LU decomposition solver.
 */
class LU : public LinearSolverBase
{
private:
  Matrix& A;
  bool factorized = false;
  bool pivot_flag;

  /**
   * The pivot mapping vector.
   * The index corresponds to the initial row number and the value to the
   * pivoted row number. This is used to map the right-hand side vector to the
   * correct row when solving.
   */
  std::vector<size_t> row_pivots;

public:

  /**
   * Default constructor.
   */
  LU(Matrix& other, const bool pivot = true);

  /**
   * Set the pivot option.
   */
  void
  pivot(const bool flag);

  /**
   * Get the pivot option.
   */
  bool
  pivot() const;

  /**
   * Factor the matrix \f$ \boldsymbol{A} \f$ into an upper and lower
   * triangular form in-place.
   *
   * An LU factorization defines the relationship
   * \f$ \boldsymbol{A} = \boldsymbol{L} \boldsymbol{U} \f$ where
   * \f$ \boldsymbol{L} \f$ is a lower triangular matrix and
   * \f$ \boldsymbol{U} \f$ is an upper triangular matrix. The factoization is
   * performed in place rather than creating an additional Matrix object.
   *
   * The algorithm used to do perform this factorization is an extension of the
   * formation of a row-echelon form matrix in that the upper triangular matrix
   * is identical to the row-echelon form. The lower triangular matrix then
   * contains the row operations used to form upper triangular system.
   */
  void
  factorize();

  /**
  * Solve an LU factored linear system.
  *
  * The LU factored linear system is define by \f$ \boldsymbol{A} \vec{x} =
  * \boldsymbol{L} \boldsymbol{U} \vec{x} = \vec{b} \f$. The solution
  * \f$ \vec{x} \f$ is obtained in a two step process. First, define \f$ \vec{y}
  * = \boldsymbol{U} \vec{x} \f$ and plug this in to obtain \f$ \boldsymbol{L}
  * \vec{y} = \vec{b} \f$. The vector \f$ \vec{y} \f$ can be obtained using
  * forward substitution after reordering the right-hand side vector
  * \f$ \vec{b} \f$ according to the pivot mapping vector. Next, the solution
  * \f$ \vec{x} \f$ is computed using the previous definition
  * \f$ \boldsymbol{U} \vec{x} = \vec{y} \f$ where \f$ \vec{y} \f$ is now the
  * source term. This system can be solved using back substitution.
  */
  void
  solve(Vector& x, const Vector& b) const override;
};

}
#endif //LU_H

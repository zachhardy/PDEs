#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include "linear_solver.h"

namespace math
{

//######################################################################
/**
 * \brief Gauss elimination direct linear solver.
 *
 * The Gauss elimination solver factors a matrix and right-hand side vector
 * into row-echelon form and then uses back-substitution to solve the system.
 * If pivoting is used, the operations are stored within a permutation matrix
 * which is used to map solutions to the original ordering after a solve.
 */
class GaussianElimination : public LinearSolver
{
protected:
  /// The permutation matrix to store pivoting operations.
  std::vector<size_t> P;

public:
  /// A flag for using pivoting in the system factorization.
  bool with_pivoting = false;
  /// A flag for normalizing the leading coefficients to one.
  bool with_normalization = false;

public:
  /// Default constructor with a matrix and right-hand side vector.
  GaussianElimination(Matrix& matrix, Vector& rhs)
    : LinearSolver(matrix, rhs)
  {}

  /// Default constructor plus options.
  GaussianElimination(Matrix& matrix, Vector& rhs, bool pivoting, bool normalize)
    : LinearSolver(matrix, rhs),
      with_pivoting(pivoting),
      with_normalization(normalize)
  {}

  void setup() override;
  Vector solve() override { return back_substitution(); }

private:
  /// Factor the matrix and right-hand side vector to row-echelon form.
  void row_echelon();

  /// Solve the row-echelon system by using back-substitution.
  Vector back_substitution();

};

}

#endif //GAUSSIAN_ELIMINATION_H

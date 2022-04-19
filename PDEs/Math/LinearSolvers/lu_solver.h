#ifndef LU_H
#define LU_H

#include "linear_solver.h"

namespace math
{

/**
 * \brief LU decomposition direct solver.
 *
 * The LU decomposition factors a matrix into a lower triangular and upper
 * triangular matrix. The upper triangular matrix is identical to that from
 * Gaussian elimination and the lower triangular matrix stores the
 * row-operations used to form the upper triangular system.
 */
class LU : public LinearSolver
{
protected:
  /// The permutation matrix to store pivoting operations.
  std::vector<size_t> P;

public:
  /// A flag for using pivoting in the system factorization.
  bool with_pivoting = false;

public:
  /// Default constructor.
  LU(Matrix& matrix, Vector& rhs) : LinearSolver(matrix, rhs) {}

  /// Default constructor plus options.
  LU(Matrix& matrix, Vector& rhs, bool pivoting)
    : LinearSolver(matrix, rhs), with_pivoting(pivoting)
  {}

  void setup() override;
  Vector solve() override;
};

}

#endif //LU_H

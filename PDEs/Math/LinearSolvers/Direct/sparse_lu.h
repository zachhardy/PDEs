#ifndef SPARSE_LU_H
#define SPARSE_LU_H

#include "LinearSolvers/linear_solver.h"

#include <vector>
#include <cstddef>


namespace pdes::Math::LinearSolver
{

/**
 * Implementation of a sparse LU solver.
 */
class SparseLU : public LinearSolverBase
{
private:
  SparseMatrix& A;
  bool factorized = false;
  bool pivot_flag = true;

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
  SparseLU(SparseMatrix& other, const bool pivot = true);

  /**
   * Set the pivot option.
   */
  void
  pivot(const bool flag);

  /**
   * Return the pivot option.
   */
  bool
  pivot() const;

  /**
   * Factor the matrix \f$ \boldsymbol{A} \f$ in-place.
   *
   * \see LU::factorize
   */
  SparseLU&
  factorize();

  /**
   * Solve the LU factored linear system.
   *
   * \param b A vector of length \f$ n \f$.
   * \param x The solution \f$ \vec{x} \f$ of
   *          \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
   *
   * \see LU::solve
   */
  void
  solve(const Vector& b, Vector& x) const override;

  /**
   * Return the solution of the LU solve.
   * \see SparseLU::solve LU::solve
   */
  Vector
  solve(const Vector& b) const override;
};

}

#endif //SPARSE_LU_H

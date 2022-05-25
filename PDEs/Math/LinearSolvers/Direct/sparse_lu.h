#ifndef SPARSE_LU_H
#define SPARSE_LU_H

#include "sparse_matrix.h"


namespace pdes::Math
{

/**
 * A class for a sparse LU decomposition solver.
 */
class SparseLU : public SparseMatrix
{
public:
  using value_type = SparseMatrix::value_type;


private:
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
  /** Copy construction from a sparse matrix. */
  SparseLU(const SparseMatrix& other, const bool pivot = true);

  /** Move construction from a sparse matrix. */
  SparseLU(SparseMatrix&& other, const bool pivot = true);

  /** Set the pivot option. */
  void
  pivot(const bool flag);

  /** Return the pivot option. */
  bool
  pivot() const;

  /**
   * Factor the matrix \f$ \boldsymbol{A} \f$ into an upper and lower triangular
   * form in-place.
   *
   * \see LU::factorize
   */
  void
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
  solve(const Vector& b, Vector& x) const;

  /**
   * Return the solution of the LU solve.
   * \see SparseLU::solve LU::solve
   */
  Vector
  solve(const Vector& b) const;
};

}

#endif //SPARSE_LU_H
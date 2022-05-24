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
  /**
   * Copy construction from a sparse matrix.
   */
  SparseLU(const SparseMatrix& other, const bool pivot = true);

  /**
   * Move construction from a sparse matrix.
   */
  SparseLU(SparseMatrix&& other, const bool pivot = true);

  /**
   * Factor the matrix \f$ \boldsymbol{A} \f$ into an upper and lower triangular
   * form in-place.
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

  void
  solve(const Vector& b, Vector& x) const;

  Vector
  solve(const Vector& b) const;
};

}

#endif //SPARSE_LU_H

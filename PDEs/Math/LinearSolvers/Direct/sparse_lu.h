#ifndef SPARSE_LU_H
#define SPARSE_LU_H

#include "sparse_matrix.h"


namespace pdes::Math
{

/**
 * A class for a sparse LU decomposition solver.
 */
template<typename number>
class SparseLU : public SparseMatrix
{
public:
  using value_type = typename SparseMatrix::value_type;

private:
  bool factorized = false;
  bool pivoting = true;

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
  SparseLU(const bool pivot = true) :
      row_pivots(0), pivoting(pivot), SparseMatrix()
  {}

  /**
   * Copy construction from a sparse matrix.
   */
  SparseLU(const SparseMatrix& other, const bool pivot = true) :
      row_pivots(other.n_rows()), pivoting(pivot), SparseMatrix(other)
  {}

  /**
   * Move construction from a sparse matrix.
   */
  SparseLU(SparseMatrix&& other, const bool pivot = true) :
      row_pivots(other.n_rows()), pivoting(pivot), SparseMatrix(other)
  {}

  /**
   * Set the pivoting flag.
   */
  void
  pivot(const bool pivot)
  { pivoting = pivot; }

  /**
   * Return the pivoting flag.
   */
  bool
  pivot() const
  { return pivoting; }

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
  factorize()
  {
    size_t n = n_rows();

    // Initialize the pivot mappings such that each row maps to itself
    for (size_t i = 0; i < n; ++i)
      row_pivots[i] = i;

    //======================================== Apply Doolittle algorithm
    for (size_t j = 0; j < n; ++j)
    {
      if (pivoting)
      {
        const value_type* a_jj = diagonal(j);

        /* Find the row containing the largest magnitude entry for column i.
         * This is only done for the sub-diagonal elements. */
        size_t argmax = j;
        value_type max = std::fabs((a_jj) ? *a_jj : 0.0);
        for (size_t k = j + 1; k < n; ++k)
        {
          const value_type* a_kj = locate(k, j);
          if (a_kj && *a_kj > max)
          {
            argmax = k;
            max = std::fabs(*a_kj);
          }
        }

        // If the sub-diagonal is uniformly zero throw error.
        Assert(max != 0.0, "Singular matrix error.");

        /* Swap the current row and the row containing the largest magnitude
         * entry corresponding for the current column. This is done to improve
         * the numerical stability of the algorithm. */
        if (argmax != j)
        {
          std::swap(row_pivots[j], row_pivots[argmax]);
          swap_row(j, argmax);
        }
      }//if pivot

      // Compute the elements of the LU decomposition
      for (size_t i = j + 1; i < n; ++i)
      {
        value_type* a_ij = locate(i, j);

        if (a_ij && *a_ij != 0.0)
        {
          /* Lower triangular components. This represents the row operations
           * performed to attain the upper-triangular, row-echelon matrix. */
          *a_ij /= *diagonal(j);


          /* Upper triangular components. This represents the row-echelon form
           * of the original matrix. Her*/
          for (const auto entry : const_row_iterator(j))
            if (entry.column > j)
              add(i, entry.column, -(*a_ij) * entry.value);
        }//if a_ij exists
      }//for rows > j
    }//for j
    factorized = true;
  }

  void
  solve(const Vector& b, Vector& x)
  {
    Assert(b.size() == n_rows(), "Dimension mismatch error.");
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    if (!factorized)
      factorize();

    //================================================== Forward solve
    size_t n = n_rows();
    b.print();
    for (size_t i = 0; i < n; ++i)
    {
      value_type value = b[row_pivots[i]];
      for (const auto elem : const_row_iterator(i))
        if (elem.column < i)
          value -= elem.value * x[elem.column];
      x[i] = value;
    }

    //================================================== Backward solve
    for (size_t i = n - 1; i != -1; --i)
    {
      value_type value = x[i];
      for (const auto elem : const_row_iterator(i))
        if (elem.column > i)
          value -= elem.value * x[elem.column];
      x[i] = value / *diagonal(i);
    }
  }

  Vector
  solve(const Vector& b)
  {
    Vector x(n_cols());
    solve(b, x);
    return x;
  }

};

}

#endif //SPARSE_LU_H

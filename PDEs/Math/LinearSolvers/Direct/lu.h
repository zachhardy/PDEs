#ifndef LU_H
#define LU_H

#include "matrix.h"


namespace math
{

/** A class for an LU decomposition solver. */
template<typename number = double>
class LU : public Matrix<number>
{
public:
  using value_type = typename Matrix<number>::value_type;
  using size_type = typename Matrix<number>::size_type;

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
  /// Default constructor.
  LU(const bool pivot = true)
    : row_pivots(0) , pivoting(pivot),
      Matrix<value_type>()
  {}

  /// Copy constructor.
  LU(const LU& other)
    : row_pivots(other.n_rows()), pivoting(other.pivoting),
      Matrix<value_type>(other.values)
  {}

  /// Move constructor.
  LU(LU&& other)
    : row_pivots(other.n_rows()), pivoting(other.pivoting),
      Matrix<value_type>(std::move(other.values))
  {}

  /// Copy from a Matrix.
  LU(const Matrix<value_type>& other, const bool pivot = true)
    : row_pivots(other.n_rows()), pivoting(pivot),
      Matrix<value_type>(other)
  {}

  /// Move from a Matrix.
  LU(Matrix<value_type>&& other, const bool pivot = true)
    : row_pivots(other.n_rows()), pivoting(pivot),
      Matrix<value_type>(std::move(other))
  {}

  /// Set the pivoting flag.
  void
  pivot(const bool pivot)
  { pivoting = pivot; }

  /// Get the pivoting flag.
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
    size_type n = this->n_rows();

    // Initialize the pivot mappings such that each row maps to itself
    for (size_type i = 0; i < n; ++i)
      row_pivots.emplace_back(i);

    //======================================== Apply Doolittle algorithm
    for (size_type i = 0; i < n; ++i)
    {
      if (pivoting)
      {
        /* Find the row containing the largest magnitude entry for column i.
         * This is only done for the sub-diagonal elements. */
        size_type argmax = i;
        value_type max = std::fabs((*this)(i, i));
        for (size_type k = i + 1; k < n; ++k)
        {
          if (std::fabs((*this)(k, i)) > max)
          {
            argmax = k;
            max = std::fabs((*this)(k, i));
          }
        }

        // If the sub-diagonal is uniformly zero, throw error
        Assert((*this)[argmax][i] != 0.0, "Singular matrix error.");

        /* Swap the current row and the row containing the largest magnitude
         * entry corresponding for the current column. This is done to improve
         * the numerical stability of the algorithm. */
        if (argmax != i)
        {
          std::swap(row_pivots[i], row_pivots[argmax]);
          this->swap_row(i, argmax);
        }
      }//if pivotong

      // Compute the elements of the LU decomposition
      for (size_type j = i + 1; j < n; ++j)
      {
        /* Lower triangular components. This represents the row operations
         * performed to attain the upper-triangular, row-echelon matrix. */
        (*this)(j, i) /= (*this)(i, i);

        /* Upper triangular components. This represents the row-echelon form of
         * the original matrix. */
        for (size_type k = i + 1; k < n; ++k)
          (*this)(j, k) -= (*this)(j, i) * (*this)(i, k);
      }
    }
    factorized = true;
  }

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
  *
  * \param b A vector of length \f$ n \f$.
  * \param x The solution \f$ \vec{x} \f$ of
  *          \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
  */
  void
  solve(const Vector<number>& b, Vector<number>& x)
  {
    Assert(b.size() == this->n_rows(), "Dimension mismatch error.");
    Assert(x.size() == this->n_cols(), "Dimension mismatch error.");
    if (!factorized)
      factorize();

    size_type n = this->n_rows();
    value_type value = 0.0;

    //================================================== Forward solve
    for (size_type i = 0; i < n; ++i)
    {
      value = b[row_pivots[i]];
      const value_type* a_ij = this->values[i].data();
      for (size_type j = 0; j < i; ++j)
        value -= *a_ij++ * x[j];
      x[i] = value;
    }

    //================================================== Backward solve
    for (size_type i = n - 1; i != -1; --i)
    {
      value = x[i];
      const value_type* a_ij = this->values[i].data();
      for (size_type j = i + 1; j < n; ++j)
        value -= *a_ij++ * x[j];
      x[i] = value / (*this)(i, i);
    }
  }

  Vector<number>
  solve(const Vector<number>& b)
  {
    Vector<value_type> x(this->n_rows(), 0.0);
    solve(b, x);
    return x;
  }
};

}
#endif //LU_H

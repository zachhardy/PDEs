#ifndef LU_H
#define LU_H

#include "matrix.h"


namespace pdes::Math
{

/** A class for an LU decomposition solver. */
class LU : public Matrix
{
public:
  using value_type = typename Matrix::value_type;

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
  LU(const bool pivot = true) :
      row_pivots(0) , pivoting(pivot), Matrix()
  {}

  /**
   * Copy construction from a matrix.
   */
  LU(const Matrix& other, const bool pivot = true) :
      row_pivots(other.n_rows()), pivoting(pivot), Matrix(other)
  {}

  /**
   * Move construction from a matrix.
   */
  LU(Matrix&& other, const bool pivot = true) :
      row_pivots(other.n_rows()), pivoting(pivot), Matrix(other)
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
      std::cout << "Starting column " << j << ":\n";

      if (pivoting)
      {
        const value_type a_jj = (*this)(j, j);

        /* Find the row containing the largest magnitude entry for column j.
         * This is only done for the sub-diagonal elements. */
        size_t argmax = j;
        value_type max = std::fabs(a_jj);
        for (size_t k = j + 1; k < n; ++k)
        {
          const value_type a_kj = (*this)(k, j);
          if (std::fabs(a_kj) > max)
          {
            argmax = k;
            max = std::fabs(a_kj);
          }
        }

        // If the sub-diagonal is uniformly zero, throw error
        Assert(max != 0.0, "Singular matrix error.");

        /* Swap the current row and the row containing the largest magnitude
         * entry corresponding for the current column. This is done to improve
         * the numerical stability of the algorithm. */
        if (argmax != j)
        {
          std::cout << "Swapping row " << j
                    << " with row " << argmax << ".\n";

          std::swap(row_pivots[j], row_pivots[argmax]);
          swap_row(j, argmax);
        }
      }//if pivoting

      // Compute the elements of the LU decomposition
      for (size_t i = j + 1; i < n; ++i)
      {
        value_type& a_ij = (*this)(i, j);

        /* Lower triangular components. This represents the row operations
         * performed to attain the upper-triangular, row-echelon matrix. */
        a_ij /= (*this)(j, j);

        /* Upper triangular components. This represents the row-echelon form of
         * the original matrix. */
        const value_type* a_jk = &coeffs[j][j + 1];
        for (size_t k = j + 1; k < n; ++k)
          (*this)(i, k) -= a_ij * *a_jk++;
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
  solve(const Vector& b, Vector& x)
  {
    Assert(b.size() == n_rows(), "Dimension mismatch error.");
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    if (!factorized)
      factorize();

    //================================================== Forward solve
    size_t n = n_rows();
    for (size_t i = 0; i < n; ++i)
    {
      value_type value = b[row_pivots[i]];
      const value_type* a_ij = data(i);
      for (size_t j = 0; j < i; ++j)
        value -= *a_ij++ * x[j];
      x[i] = value;
    }

    //================================================== Backward solve
    for (size_t i = n - 1; i != -1; --i)
    {
      value_type value = x[i];
      const value_type* a_ij = &coeffs[i][i + 1];
      for (size_t j = i + 1; j < n; ++j)
        value -= *a_ij++ * x[j];
      x[i] = value / (*this)(i, i);
    }
  }

  Vector
  solve(const Vector& b)
  {
    Vector x(n_rows(), 0.0);
    solve(b, x);
    return x;
  }
};

}
#endif //LU_H

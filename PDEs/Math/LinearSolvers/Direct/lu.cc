#include "lu.h"

#include <cinttypes>

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
template<typename value_type>
void math::LU<value_type>::setup()
{
  auto& A = this->A;
  uint64_t n = A.n_rows();

  // Initialize the pivot mappings such that each row maps to itself
  for (uint64_t i = 0; i < n; ++i)
    row_pivots.emplace_back(i);

  //======================================== Apply Doolittle algorithm
  for (uint64_t i = 0; i < n; ++i)
  {
    if (pivot)
    {
      /* Find the row containing the largest magnitude entry for column i.
       * This is only done for the sub-diagonal elements. */
      uint64_t argmax = i;
      value_type max = std::fabs(A[i][i]);
      for (uint64_t k = i + 1; k < n; ++k)
      {
        if (std::fabs(A[k][i]) > max)
        {
          argmax = k;
          max = std::fabs(A[k][i]);
        }
      }

      // If the sub-diagonal is uniformly zero, throw error
      Assert(A[argmax][i] != 0.0, "Singular matrix error.");

      /* Swap the current row and the row containing the largest magnitude
       * entry corresponding for the current column. This is done to improve
       * the numerical stability of the algorithm. */
      if (argmax != i)
      {
        std::swap(row_pivots[i], row_pivots[argmax]);
        A.swap_row(i, argmax);
      }
    }//if pivot

    // Compute the elements of the LU decomposition
    for (uint64_t j = i + 1; j < n; ++j)
    {
      /* Lower triangular components. This represents the row operations
       * performed to attain the upper-triangular, row-echelon matrix. */
      A[j][i] /= A[i][i];

      /* Upper triangular components. This represents the row-echelon form of
       * the original matrix. */
      for (uint64_t k = i + 1; k < n; ++k)
        A[j][k] -= A[j][i] * A[i][k];
    }
  }
  this->initialized = true;
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
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
template<typename value_type>
math::Vector<value_type>
math::LU<value_type>::solve(const Vector<value_type>& b)
{
  auto& A = this->A;
  Assert(b.size() == A.n_rows(), "Mismatched size error.");
  if (not this->initialized)
    this->setup();

  uint64_t n = A.n_rows();
  value_type value = 0.0;

  // Forward solve
  Vector x(n, 0.0);
  for (uint64_t i = 0; i < n; ++i)
  {
    value = b[row_pivots[i]];
    for (uint64_t j = 0; j < i; ++j)
      value -= A[i][j] * x[j];
    x[i] = value;
  }

  // Backward solve
  for (int i = n - 1; i >= 0; --i)
  {
    value = x[i];
    for (uint64_t j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}


template class math::LU<double>;
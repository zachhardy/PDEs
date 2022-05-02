#include "lu.h"


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
  size_t n = A.n_rows();

  // Initialize the pivot mappings such that each row maps to itself
  for (size_t i = 0; i < n; ++i)
    row_pivots.emplace_back(i);

  // Go through each column in the matrix.
  for (size_t j = 0; j < n; ++j)
  {
    /* Find the row index for the largest magnitude entry in this column.
     * This is only done for sub-diagonal elements. */
    value_type max = 0.0;
    size_t argmax = j;
    for (size_t i = j; i < n; ++i)
    {
      if (std::fabs(A[i][j]) > max)
      {
        max = std::fabs(A[i][j]);
        argmax = i;
      }
    }

    // If the sub-diagonal is uniformly zero, throw error
    Assert(A[argmax][j] != 0.0,
           "Degenerate matrix encountered. "
           "Specifically, all elements on the sub-diagonal were zero.");

    /* Swap the current row and the row containing the largest magnitude
     * entry corresponding for the current column. This is done to improve
     * the numerical stability of the algorithm. */
    if (pivot and argmax != j)
    {
      std::swap(row_pivots[j], row_pivots[argmax]);
      A.swap_row(j, argmax);
    }

    /* Perform row-wise operations such that all sub-diagonal values are zero.
     * This is done by subtracting the current row times the ratio of the
     * sub-diagonal and the current row's leading value. */
    for (size_t i = j + 1; i < n; ++i)
    {
      value_type factor = A[i][j] / A[j][j];
      for (size_t k = j + 1; k < n; ++k)
        A[i][k] -= A[j][k] * factor;
      A[i][j] = factor;
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

  size_t n = A.n_rows();
  value_type value = 0.0;

  // Forward solve
  Vector x(n, 0.0);
  for (size_t i = 0; i < n; ++i)
  {
    value = b[row_pivots[i]];
    for (size_t j = 0; j < i; ++j)
      value -= A[i][j] * x[j];
    x[i] = value;
  }

  // Backward solve
  for (int i = n - 1; i >= 0; --i)
  {
    value = x[i];
    for (size_t j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}


template class math::LU<double>;

#include "cholesky.h"

#include <cinttypes>

/**
 * Perform a Cholesky factorization on the matrix \f$ \boldsymbol{A} \f$.
 *
 * Cholesky factorization is meant for symmetric positive definite matrices into
 * a lower triangular matrix and its transpose. This method is more efficient
 * than the LU decomposition when applicable.
 *
 * \note Checks are not performed to ensure symetric positive definiteness. The
 *       user is responsible for ensuring the matrix fits this criteria.
 */
template<typename number>
void math::Cholesky<number>::setup()
{
  auto& A = this->A;
  uint64_t n = A.n_rows();

  // Factor the matrix column-wise
  for (uint64_t j = 0; j < n; ++j)
  {
    // Set the diagonal element
    number sum = 0.0;
    for (uint64_t k = 0; k < j; ++k)
      sum += A[j][k] * A[j][k];
    A[j][j] = std::sqrt(A[j][j] - sum);

    // Set the off-diagonals
    for (uint64_t i = j + 1; i < n; ++i)
    {
      sum = 0.0;
      for (uint64_t k = 0; k < j; ++k)
        sum += A[i][k] * A[j][k];
      A[i][j] = (A[i][j] - sum) / A[j][j];
    }
  }
  this->initialized = true;
}


/**
 * Solve the Cholesky factored linear system.
 *
 * The Cholesky solve is a specialization of the LU solve in that
 * \f$ \boldsymbol{U} = \boldsymbol{L}^T \f$. See \ref lu_solve for
 * implementation detail.
 *
 * \param b A vector of length \f$ n \f$.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
template<typename number>
math::Vector<number>
math::Cholesky<number>::solve(const Vector<number>& b)
{
  auto& A = this->A;
  Assert(b.size() == A.n_rows(), "Mismatched size error.");
  if (not this->initialized)
    this->setup();

  uint64_t n = A.n_rows();
  number value = 0.0;

  // Forward solve
  Vector x(n, 0.0);
  for (uint64_t i = 0; i < n; ++i)
  {
    value = b[i];
    for (uint64_t j = 0; j < i; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }

  // Backward solve
  for (int i = n - 1; i >= 0; --i)
  {
    value = x[i];
    for (uint64_t j = i + 1; j < n; ++j)
      value -= A[j][i] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}

template class math::Cholesky<double>;

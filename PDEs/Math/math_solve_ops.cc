#include "math.h"
#include "matrix.h"
#include "vector.h"

#include <sstream>

/**
 * \brief Solve a linear system using back substitution.
 * \param A An upper triangular \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 *
 * Back substitution involves solving each equation in the system from last to
 * first. This requires an upper triangular matrix, where each equation is
 * decoupled from those above it. This allows for the direct computation of all
 * unknowns when solving from the last equation to the first.
 */
math::Vector<double>
math::back_substitution(const Matrix<double>& A,
                        const Vector<double>& b)
{
  bool is_upper = true;
  for (size_t i = 0; i < b.size(); ++i)
  {
    for (size_t j = 0; j < i; ++j)
      if (A[i][j] != 0.0) { is_upper = false; break; }
    if (not is_upper) break;
  }
  if (not is_upper)
  {
    std::stringstream err;
    err << __FUNCTION__ << ": "
        << "The matrix must be upper triangular to use back substitution.";
    throw std::runtime_error(err.str());
  }

  Vector<double> x(b.size(), 0.0);
  for (int i = b.size() - 1; i >= 0; --i)
  {
    double value = b[i];
    for (int j = i + 1; j < b.size(); ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}

/**
 * \brief Solve a linear system using forward substitution.
 * \param A An lower triangular \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 *
 * Forward substitution involves solving each equation in the system from first
 * to last. This requires a lower triangular matrix, where each equation is
 * decoupled from those below it. This allows for the direct computation of all
 * unknowns when solving from the first equation to the last.
 */
math::Vector<double>
math::forward_substitution(const Matrix<double>& A,
                           const Vector<double>& b)
{
  bool is_lower = true;
  for (size_t i = 0; i < b.size(); ++i)
  {
    for (size_t j = i + 1; j < b.size(); ++j)
      if (A[i][j] != 0.0) { is_lower = false; break; }
    if (not is_lower) break;
  }
  if (not is_lower)
  {
    std::stringstream err;
    err << __FUNCTION__ << ": "
        << "The matrix must be lower triangular to use forward substitution.";
    throw std::runtime_error(err.str());
  }

  Vector<double> x(b.size(), 0.0);
  for (size_t i = 0; i < b.size(); ++i)
  {
    double value = b[i];
    for (size_t j = 0; j < i; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}

/**
 * \brief Solve a linear system using Gaussian elimination.
 * \param A An \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
math::Vector<double>
math::gaussian_elimination(Matrix<double>& A,
                           Vector<double>& b,
                           const bool pivot)
{
  row_echelon_form(A, b, pivot);
  return back_substitution(A, b);
}

/**
 * \brief Solve an LU factored linear system.
 * \param A An LU factored \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \param P The pivot mapping vector to reorder to right-hand side.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
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
 */
math::Vector<double> math::lu_solve(const Matrix<double>& A,
                            const Vector<double>& b,
                            const std::vector<size_t> P)
{
  size_t n = b.size();

  // Forward solve
  Vector<double> y(b.size(), 0.0);
  for (size_t i = 0; i < b.size(); ++i)
  {
    double value = b[i];
    for (size_t j = 0; j < i; ++j)
      value -= A[i][j] * y[j];
    y[i] = value;
  }

  // Backward solve
  Vector<double> x(b.size(), 0.0);
  for (int i = n - 1; i >= 0; --i)
  {
    double value = y[i];
    for (size_t j = i + 1; j < n; ++j)
      value -= A[i][j] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}

/**
 * \brief Solve the Cholesky factored linear system.
 * \param A A Cholesky factored \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \return The solution \f$ \vec{x} \f$ of
 *         \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 *
 * The Cholesky solve is a specialization of the LU solve in that
 * \f$ \boldsymbol{U} = \boldsymbol{L}^T \f$. See \ref lu_solve for
 * implementation detail.
 */
math::Vector<double>
math::cholesky_solve(const Matrix<double>& A,
                     const Vector<double>& b)
{
  size_t n = b.size();

  // Forward solve
  Vector<double> y(n, 0.0);
  for (size_t i = 0; i < n; ++i)
  {
    double value = b[i];
    for (size_t j = 0; j < i; ++j)
      value -= A[i][j] * y[j];
    y[i] = value / A[i][i];
  }

  // Backward solve
  Vector<double> x(n, 0.0);
  for (int i = n - 1; i >= 0; --i)
  {
    double value = y[i];
    for (size_t j = i + 1; j < n; ++j)
      value -= A[j][i] * x[j];
    x[i] = value / A[i][i];
  }
  return x;
}
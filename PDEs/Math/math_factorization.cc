#include "math.h"
#include "matrix.h"
#include "vector.h"

/**
 * \brief Factor a matrix and vector into row-echelon form.
 * \param A An \f$ n \times n \f$ matrix.
 * \param b A vector of length \f$ n \f$.
 * \param pivot A flag for whether pivoting is performed.
 *
 * Row-echelon form creates an upper triangular system. This is done by
 * traversing through each column and performing row-wise operations to
 * eliminate all sub-diagonal coefficients. For numerical stability, pivoting
 * can be used to ensure that for each subsequent column, division by the
 * largest element is carried out. When this is performed, the pivoting is
 * stored in a vector which maps the initial row to the final, pivoted position.
 * This is used to reorder solutions into the initial ordering. Additionally,
 * normalization can be performed to ensure that the diagonal elements are
 * unity.
 *
 */
void math::row_echelon_form(Matrix& A, Vector& b, const bool pivot)
{
  if (A.n_rows() != A.n_cols() or A.n_rows() != b.size())
  {
    std::stringstream err;
    err << __FUNCTION__ << ": Invalid inputs. "
        << "The matrix must be square and agree with the size of the vector.";
    throw std::runtime_error(err.str());
  }

  size_t n = A.n_rows();

  // Go through each column in the matrix.
  for (size_t column = 0; column < n; ++column)
  {
    /* Find the row index for the largest magnitude entry in this column.
     * This is only done for sub-diagonal elements. */
    double max = 0.0;
    size_t argmax = column;
    for (size_t i = column; i < n; ++i)
    {
      if (std::fabs(A[i][column]) > max)
      {
        max = std::fabs(A[i][column]);
        argmax = i;
      }
    }

    // If the sub-diagonal is uniformly zero, throw error
    if (A[argmax][column] == 0.0)
    {
      std::stringstream err;
      err << __FUNCTION__ << ": Degenerate matrix."
          << "\nThis was encountered on column " << column << "."
          << "\nSpecifically, no non-zero entries were found on either the "
          << "diagonal or the sub-diagonal elements.";
      throw std::runtime_error(err.str());
    }

    /* Swap the current row and the row containing the largest magnitude
     * entry corresponding for the current column. This is done to improve
     * the numerical stability of the algorithm. */
    if (pivot and argmax != column)
    {
      std::swap(b[column], b[argmax]);
      A.swap_row(column, argmax);
    }

    /* Perform row-wise operations such that all sub-diagonal values are zero.
     * This is done by subtracting the current row times the ratio of the
     * sub-diagonal and the current row's leading value. */
    for (size_t i = column + 1; i < n; ++i)
    {
      double factor = A[i][column] / A[column][column];
      for (size_t j = column; j < n; ++j)
        A[i][j] -= A[column][j] * factor;
      b[i] -= b[column] * factor;
    }
  }
}

/**
 * \brief Factor the matrix \f$ \boldsymbol{A} \f$ into an upper and lower
 *        triangular form in place.
 * \param A A \f$ n \times n \f$ matrix.
 * \param pivot A flag for whether pivoting is performed.
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
 *
 */
std::vector<size_t> math::lu_factorization(Matrix& A, const bool pivot)
{
  if (A.n_rows() != A.n_cols())
  {
    std::stringstream err;
    err << __FUNCTION__ << ": "
        << "LU factorization can only be performed on square matrices.";
    throw std::runtime_error(err.str());
  }

  size_t n = A.n_rows();

  // Initialize the pivot mappings to map each row to itself
  std::vector<size_t> P;
  for (size_t i = 0; i < n; ++i)
    P.emplace_back(i);

  // Go through each column in the matrix.
  for (size_t column = 0; column < n; ++column)
  {
    /* Find the row index for the largest magnitude entry in this column.
     * This is only done for sub-diagonal elements. */
    double max = 0.0;
    size_t argmax = column;
    for (size_t i = column; i < n; ++i)
    {
      if (std::fabs(A[i][column]) > max)
      {
        max = std::fabs(A[i][column]);
        argmax = i;
      }
    }

    // If the sub-diagonal is uniformly zero, throw error
    if (A[argmax][column] == 0.0)
    {
      std::stringstream err;
      err << __FUNCTION__ << ": Degenerate matrix."
          << "\nThis was encountered on column " << column << "."
          << "\nSpecifically, no non-zero entries were found on either the "
          << "diagonal or the sub-diagonal elements.";
      throw std::runtime_error(err.str());
    }

    /* Swap the current row and the row containing the largest magnitude
     * entry corresponding for the current column. This is done to improve
     * the numerical stability of the algorithm. */
    if (pivot and argmax != column)
    {
      std::swap(P[column], P[argmax]);
      A.swap_row(column, argmax);
    }

    /* Perform row-wise operations such that all sub-diagonal values are zero.
     * This is done by subtracting the current row times the ratio of the
     * sub-diagonal and the current row's leading value. */
    for (size_t i = column + 1; i < n; ++i)
    {
      double factor = A[i][column] / A[column][column];
      for (size_t j = column + 1; j < n; ++j)
        A[i][j] -= A[column][j] * factor;
      A[i][column] = factor;
    }
  }
  return P;
}

/**
 * \brief Perform a Cholesky factorization on the matrix \f$ \boldsymbol{A} \f$.
 * \param A A \f$ n \times n \f$ symmetric positive defined matrix.
 *
 * Cholesky factorization is meant for symmetric positive definite matrices into
 * a lower triangular matrix and its transpose. This method is more efficient
 * than the LU decomposition when applicable.
 *
 * \note Checks are not performed to ensure symetric positive definiteness. The
 *       user is responsible for ensuring the matrix fits this criteria.
 */
void math::cholesky_factorization(Matrix& A)
{
  if (A.n_rows() != A.n_cols())
  {
    std::stringstream err;
    err << __FUNCTION__ << ": "
        << "LU factorization can only be performed on square matrices.";
    throw std::runtime_error(err.str());
  }

  // Factor the matrix column-wise
  for (size_t j = 0; j < A.n_cols(); ++j)
  {
    // Set the diagonal element
    double sum = 0.0;
    for (size_t k = 0; k < j; ++k)
      sum += A[j][k] * A[j][k];
    A[j][j] = std::sqrt(A[j][j] - sum);

    // Set the off-diagonals
    for (size_t i = j + 1; i < A.n_rows(); ++i)
    {
      sum = 0.0;
      for (size_t k = 0; k < j; ++k)
        sum += A[i][k] * A[j][k];
      A[i][j] = (A[i][j] - sum) / A[j][j];
      A[j][i] = 0.0;
    }
  }
}


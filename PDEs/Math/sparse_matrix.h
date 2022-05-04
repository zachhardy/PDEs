#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "vector.h"
#include "exceptions.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

namespace math
{

template<typename value_type>
class SparseMatrix
{
private:
  size_t rows;
  size_t cols;

  std::vector<std::vector<value_type>> m_data;
  std::vector<std::vector<size_t>> m_indices;

public:
  /** Default constructor. */
  SparseMatrix() = default;

  /** Construct a square sparse matrix with dimension \p n. */
  explicit SparseMatrix(const size_t n)
      : rows(n), cols(n), m_data(n), m_indices(n)
  {}

  /** Construct a sparse matrix with \p n_rows and \p n_cols. */
  explicit SparseMatrix(const size_t n_rows, const size_t n_cols)
    : rows(n_rows), cols(n_cols), m_data(n_rows), m_indices(n_rows)
  {}

  /** Construct a sparsity pattern. */
  SparseMatrix(const std::vector<std::vector<size_t>>& pattern)
    : rows(pattern.size()), m_data(pattern.size()), m_indices(pattern)
  {
    cols = 0;
    for (size_t i = 0; i < rows; ++i)
    {
      // Sort column indices for the row
      std::stable_sort(m_indices[i].begin(), m_indices[i].end());

      // Resize the data vector for the row
      m_data[i].resize(m_indices[i].size(), 0.0);

      // Check the largest column index
      for (const auto& j : m_indices[i])
        if (j > cols) cols = j;
    }
  }

  /** Copy constructor. */
  SparseMatrix(const SparseMatrix& other)
    : rows(other.rows), cols(other.cols),
      m_data(other.m_data), m_indices(other.m_indices)
  {}

  /** Move constructor. */
  SparseMatrix(SparseMatrix&& other)
    : rows(other.rows), cols(other.cols),
      m_data(std::move(other.m_data)),
      m_indices(std::move(other.m_indices))
  {}

  /** Copy assignment operator. */
  SparseMatrix& operator=(const SparseMatrix& other)
  {
    rows = other.rows;
    cols = other.cols;
    m_data = other.m_data;
    m_indices = other.m_indices;
    return *this;
  }

  /** Move assignment operator. */
  SparseMatrix& operator=(SparseMatrix&& other)
  {
    rows = other.rows;
    cols = other.cols;
    m_data = std::move(other.m_data);
    m_indices = std::move(other.m_indices);
    return *this;
  }

public:
  /** \name Accessors */
  /** @{ */

  /** Read access for the element at row \p i and column \p j. */
  value_type operator()(const size_t i, const size_t j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // If row is uninitialized, return zero
    if (m_indices[i].empty())
      return 0.0;

    // Otherwise, look for the specified column in the row
    const auto& colnums = m_indices[i];
    auto rel_loc = std::lower_bound(colnums.begin(), colnums.end(), j);

    // If the column is not in the row, return zero
    if (rel_loc == colnums.end())
      return 0.0;

    // Otherwise, return the element
    size_t jr = rel_loc - colnums.begin();
    return m_data[i][jr];
  }

  /** Read/write access for the element at row \p i and column \p j. */
  value_type& operator()(const size_t i, const size_t j)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // Check whether the row exists
    const auto& colnums = m_indices[i];
    Assert(not colnums.empty(),
           "Invalid access attempt. Element not initialized.");

    // Check whether the column exists on the row
    auto rel_loc = std::lower_bound(colnums.begin(), colnums.end(), j);
    Assert(rel_loc != colnums.end(),
           "Invalid access attempt. Element not initialized.");

    size_t jr = rel_loc - colnums.begin();
    return m_data[i][jr];
  }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /** Clear all data from the sparse matrix. */
  void clear()
  {
    rows = cols = 0;
    m_data.clear();
    m_indices.clear();
  }

  /** Set the element at row \p i and column \p j to \p value. */
  void set(const size_t i, const size_t j, const value_type value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

   /* If the row is empty or the column number is larger than all current
    * entries on the row, add to the back of the row. */
   if (m_indices[i].size() == 0 or m_indices[i].back() < j)
   {
     m_indices[i].push_back(j);
     m_data[i].push_back(value);
     return;
   }

   // Find the index to insert which maintains sorting
   auto& colnums = m_indices[i];
   auto rel_loc = std::lower_bound(colnums.begin(), colnums.end(), j);

   // If this points to an existing column, add to it
   if (*rel_loc == j)
   {
     size_t jr = rel_loc - colnums.begin();
     m_data[i][jr] += value;
     return;
   }

   // Otherwise, insert into its sorted location
   m_indices[i].insert(rel_loc, j);
   m_data[i].insert(rel_loc, value);
  }

  /** Set a list of elements. See \ref set.*/
  void set(const std::vector<size_t>& row_indices,
           const std::vector<size_t>& col_indices,
           const std::vector<value_type>& values)
  {
    Assert(row_indices.size() == col_indices.size() &&
           row_indices.size() == values.size(),
           "All inputs must be of the same length.");

    for (size_t i = 0; i < row_indices.size(); ++i)
      this->set(row_indices[i], col_indices[i], values[i]);
  }

  /** @} */
  /** \name Information */
  /** @{ */

  /** Return the number of rows the sparse matrix represents. */
  size_t n_rows() const
  {
    return rows;
  }

  /** Return the number of columns the sparse matrix represents. */
  size_t n_cols() const
  {
    return cols;
  }

  /** Return the number of nonzero elements in the sparse matrix. */
  size_t nnz() const
  {
    size_t nnz = 0;
    for (const auto& colnums : m_indices)
      nnz += colnums.size();
    return nnz;
  }

  /** @} */
  /** \name Scalar Operations */
  /** @{ */

  /** Element-wise multiplication by a scalar. */
  SparseMatrix operator*(const value_type value) const
  {
    SparseMatrix A(*this);
    for (auto& row_data : A)
      for (auto& elem : row_data)
        elem *= value;
    return A;
  }

  /** Element-wise multiplication by a scalar in-place. */
  SparseMatrix& operator*=(const value_type value)
  {
    for (auto& row_data : m_data)
      for (auto& elem : row_data)
        elem *= value;
    return *this;
  }

  /** Element-wise division by a scalar. */
  SparseMatrix operator/(const value_type value) const
  {
    Assert(value != 0.0, "Zero division error.");

    SparseMatrix A(*this);
    for (auto& row_data : A)
      for (auto& elem : row_data)
        elem /= value;
    return A;
  }

  /** Element-wise division by a scalar in-place. */
  SparseMatrix& operator/=(const value_type value)
  {
    Assert(value != 0.0, "Zero division error.");

    for (auto& row_data : m_data)
      for (auto& elem : row_data)
        elem /= value;
    return *this;
  }

  /** @} */
  /** \name Linear Algebra */
  /** @{ */

  /** Element-wise addition of two sparse matrices. */
  SparseMatrix operator+(const SparseMatrix& other) const
  {
    Assert(rows == other.n_rows() &&
           cols == other.n_cols(),
           "Dimension mismatch error.");

    SparseMatrix A(*this);
    for (size_t i = 0; i < rows; ++i)
      for (size_t jr = 0; jr < other.m_indices[i].size(); ++jr)
        A.set(i, other.m_indices[i][jr], other.m_data[i][jr]);
    return A;
  }

  /** Element-wise addition of two sparse matrices in-place. */
  SparseMatrix& operator+=(const SparseMatrix& other)
  {
    Assert(rows == other.n_rows() &&
           cols == other.n_cols(),
           "Dimension mismatch error.");

    for (size_t i = 0; i < rows; ++i)
      for (size_t jr = 0; jr < other.m_indices[i].size(); ++jr)
        this->set(i, other.m_indices[i][jr], other.m_data[i][jr]);
    return *this;
  }

  /** Element-wise subtraction of two sparse matrices. */
  SparseMatrix operator-(const SparseMatrix& other) const
  {
    Assert(rows == other.n_rows() &&
           cols = other.n_cols(),
           "Dimension mismatch error.");

    SparseMatrix A(*this);
    for (size_t i = 0; i < rows; ++i)
      for (size_t jr = 0; jr < other.m_indices[i].size(); ++jr)
        A.set(i, other.m_indices[i][jr], -other.m_data[i][jr]);
    return A;
  }

  /** Element-wise subtraction of two sparse matrices in-place. */
  SparseMatrix& operator-=(const SparseMatrix& other)
  {
    Assert(rows == other.n_rows() &&
           cols == other.n_cols(),
           "Dimension mismatch error.");

    for (size_t i = 0; i < rows; ++i)
      for (size_t jr = 0; jr < other.m_indices[i].size(); ++jr)
        this->set(i, other.m_indices[i][jr], -other.m_data[i][jr]);
    return *this;
  }

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector<value_type> operator*(const Vector<value_type>& x) const
  {
    Assert(cols == x.size(), "Dimension mismatch error.");

    Vector<value_type> Ax(cols);
    for (size_t i = 0; i < rows; ++i)
      for (size_t jr = 0; jr < m_indices[i].size(); ++jr)
        Ax[i] += m_data[i][jr] * x[m_indices[i][jr]];
    return Ax;
  }

  /**
   * Compute a matrix-matrix product.
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}
   *     c_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  SparseMatrix operator*(const SparseMatrix& other) const
  {
    Assert(cols == other.n_rows(), "Dimension mismatch error.");

    SparseMatrix A(rows, other.n_cols());
    for (size_t i = 0; i < rows; ++i)
    {
      // Loop over relative column indices for row i of this
      for (size_t jr = 0; jr < m_indices[i].size(); ++jr)
      {
        size_t j = m_indices[i][jr];
        value_type a_ij = m_data[i][jr];

        // Loop over relative column indices of row j of other
        for (size_t kr = 0; kr < other.m_indices[j].size(); ++kr)
        {
          // Compute c_{ik} += a_{ij} b_{jk}
          size_t k = other.m_indices[j][kr];
          value_type c_ik = a_ij * other.m_data[j][kr];
          A.set(i, k, c_ik);
        }//for columns in row j of other matrix
      }//for columns in row i of this matrix
    }//for row
  }

  void vmult(const Vector<value_type>& x,
             Vector<value_type>        dst)
  {
    Assert(rows == x.size(), "Dimension mismatch error.");
    Assert(cols == dst.size(), "Dimension mismatch error.");

    for (size_t i = 0; i < rows; ++i)
    {
      value_type value = 0.0;
      for (size_t jr = 0; jr < m_indices[i].size(); ++jr)
        value += m_data[i][jr] * x[m_indices[i][jr]];
      dst[i] = value;
    }
  }

  /** @} */
};


/*-------------------- Inline Implementations --------------------*/


/** Element-wise multiplication by a scalar value. */
template<typename value_type>
inline SparseMatrix<value_type>
operator*(const value_type value, const SparseMatrix<value_type>& A)
{
  return A * value;
}

}
#endif //SPARSE_MATRIX_H

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
  SparseMatrix(const std::vector<std::vector<size_t>>& pattern);

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
  SparseMatrix& operator=(const SparseMatrix& other);

  /** Move assignment operator. */
  SparseMatrix& operator=(SparseMatrix&& other);

public:
  /** \name Accessors */
  /** @{ */

  /** Read access for the element at row \p i and column \p j. */
  value_type operator()(const size_t i, const size_t j) const;

  /** Read/write access for the element at row \p i and column \p j. */
  value_type& operator()(const size_t i, const size_t j);

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /** Clear all data from the sparse matrix. */
  void clear();

  /** Set the element at row \p i and column \p j to \p value. */
  void set(const size_t i, const size_t j, const value_type value);

  /** Set a list of elements. See \ref set.*/
  void set(const std::vector<size_t>& row_indices,
           const std::vector<size_t>& col_indices,
           const std::vector<value_type>& values);


  /** @} */
  /** \name Information */
  /** @{ */

  /** Return the number of rows the sparse matrix represents. */
  size_t n_rows() const { return rows; }

  /** Return the number of columns the sparse matrix represents. */
  size_t n_cols() const { return cols; }

  /** Return the total number of elements in the sparse matrix. */
  size_t size() const { return rows * cols; }

  /** Return the number of nonzero elements in the sparse matrix. */
  size_t n_nonzero_elements() const;

  /** @} */
  /** \name Scalar Operations */
  /** @{ */

  /** Element-wise multiplication by a scalar. */
  SparseMatrix operator*(const value_type value) const;

  /** Element-wise multiplication by a scalar in-place. */
  SparseMatrix& operator*=(const value_type value);

  /** Element-wise division by a scalar. */
  SparseMatrix operator/(const value_type value) const;

  /** Element-wise division by a scalar in-place. */
  SparseMatrix& operator/=(const value_type value);

  /** @} */
  /** \name Linear Algebra */
  /** @{ */

  /** Element-wise addition of two sparse matrices. */
  SparseMatrix operator+(const SparseMatrix& other) const;

  /** Element-wise addition of two sparse matrices in-place. */
  SparseMatrix& operator+=(const SparseMatrix& other);

  /** Element-wise subtraction of two sparse matrices. */
  SparseMatrix operator-(const SparseMatrix& other) const;

  /** Element-wise subtraction of two sparse matrices in-place. */
  SparseMatrix& operator-=(const SparseMatrix& other);

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector<value_type> operator*(const Vector<value_type>& x) const;

  /**
   * Compute a matrix-matrix product.
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}
   *     c_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
   * \f]
   */
   SparseMatrix operator*(const SparseMatrix& other) const;

  /** @} */
};

/*-------------------- Inline Implementations --------------------*/

template<typename value_type>
inline SparseMatrix<value_type>::
SparseMatrix<value_type>(const std::vector<std::vector<size_t>>& pattern)
  : rows(pattern.size()), m_data(pattern.size()), m_indices(pattern)
{
  cols = 0;
  for (size_t i = 0; i < rows; ++i)
  {
    // Sort column indices for the row
    std::stable_sort(m_indices[i].begin(), m_indices[i].end());

    // Resize the data vector for the row
    m_data[i].resize(m_indices[i].size(), 0.0);

    // Check for largest column index
    for (const auto& j : m_indices[i])
      if (j > cols) cols = j;
  }
}



template<typename value_type>
inline SparseMatrix<value_type>&
SparseMatrix<value_type>::operator=(const SparseMatrix<value_type>& other)
{
  rows = other.rows;
  cols = other.cols;
  m_data = other.m_data;
  m_indices = other.m_indices;
  return *this;
}



template<typename value_type>
inline SparseMatrix<value_type>&
SparseMatrix<value_type>::operator=(SparseMatrix<value_type>&& other)
{
  rows = other.rows;
  cols = other.cols;
  m_data = std::move(other.m_data);
  m_indices = std::move(other.m_indices);
  return *this;
}



template<typename value_type>
inline value_type
SparseMatrix<value_type>::operator()(const size_t i, const size_t j) const
{
  Assert(i < rows and j < cols, "Indices exceed the matrix dimensions.");

  // If no entries on the row, return the default.
  if (m_indices[i].empty())
    return 0.0;

  // Otherwise, look for the specified column in the row
  const auto& colnums = m_indices[i];
  auto rel_loc = std::find(colnums.begin(), colnums.end(), j);

  // If the column is not in the row, return the default
  if (rel_loc == colnums.end())
    return 0.0;

  // Otherwise, return the matrix element
  size_t jr = rel_loc - colnums.begin();
  return m_data[i][jr];
}



template<typename value_type>
inline value_type&
SparseMatrix<value_type>::operator()(const size_t i, const size_t j)
{
  Assert(i < rows and j < cols, "Indices exceed the matrix dimensions.");

  // Check whether the row exists
  const auto& colnums = m_indices[i];
  Assert(not colnums.empty(), "Invalid access attempt.");

  // Check whether the column exists on the row
  auto rel_loc = std::lower_bound(colnums.begin(), colnums.end(), j);
  Assert(rel_loc != colnums.end(), "Invalid access attempt.");

  size_t jr = rel_loc - colnums.begin();
  return m_data[i][jr];
}



template<typename value_type>
inline void SparseMatrix<value_type>::clear()
{
  rows = 0; cols = 0;
  m_data.clear(); m_indices.clear();
}



template<typename value_type>
inline void
SparseMatrix<value_type>::set(const size_t i, const size_t j,
                              const value_type value)
{
  Assert(i < rows and j < cols, "Indices exceed the matrix dimensions.");

  /* If the row is empty or the column number is larger than all current
   * entries on the row, add to the back of the row. */
  if (m_indices[i].size() == 0 or m_indices[i].back() < j)
  {
    m_indices[i].push_back(j);
    m_data[i].push_back(value);
    return;
  }

  // Find the index which maintains sorting.
  auto it = std::lower_bound(m_indices[i].begin(),
                             m_indices[i].end(), j);

  // If this points to an existing column, add to it
  if (*it == j)
  {
    size_t jr = it - m_indices[i].begin();
    m_data[i][jr] += value;
    return;
  }

  // Otherwise, insert into its sorted location
  m_indices[i].insert(it, j);
  m_data[i].insert(it, j);
}



template<typename value_type>
inline void
SparseMatrix<value_type>::set(const std::vector<size_t>& row_indices,
                              const std::vector<size_t>& col_indices,
                              const std::vector<value_type>& values)
{
  Assert(row_indices.size() == col_indices.size() and
         row_indices.size() == values.size(),
         "Invalid inputs.");
  for (size_t i = 0; i < row_indices.size(); ++i)
    this->set(row_indices[i], col_indices[i], values[i]);
}



template<typename value_type>
inline size_t
SparseMatrix<value_type>::n_nonzero_elements() const
{
  size_t count = 0;
  for (const auto& row : m_indices)
    count += row.size();
  return count;
}



template<typename value_type>
inline SparseMatrix<value_type>
SparseMatrix<value_type>::operator*(const value_type value) const
{
  SparseMatrix m(*this);
  for (auto& row : m)
    for (auto& elem : row)
      elem *= value;
  return m;
}



template<typename value_type>
inline SparseMatrix<value_type>&
SparseMatrix<value_type>::operator*=(const value_type value)
{
  for (auto& row : m_data)
    for (auto& elem : row)
      elem *= value;
  return *this;
}



template<typename value_type>
inline SparseMatrix<value_type>
SparseMatrix<value_type>::operator/(const value_type value) const
{
  Assert(value != 0.0, "Zero division error.");

  SparseMatrix m(*this);
  for (auto& row : m)
    for (auto& elem : row)
      elem /= value;
  return m;
}



template<typename value_type>
inline SparseMatrix<value_type>&
SparseMatrix<value_type>::operator/=(const value_type value)
{
  Assert(value != 0.0, "Zero division error.");

  for (auto& row : m_data)
    for (auto& elem : row)
      elem /= value;
  return *this;
}



template<typename value_type>
inline SparseMatrix<value_type>
SparseMatrix<value_type>::
operator+(const SparseMatrix<value_type>& other) const
{
  Assert(rows == other.n_rows() and
         cols == other.n_cols(),
         "Mismatched size error.");

  SparseMatrix m(*this);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < other.m_indices[i].size(); ++j)
      m.set(i, other.m_indices[i][j], other.m_data[i][j]);
  return m;
}



template<typename value_type>
inline SparseMatrix<value_type>&
SparseMatrix<value_type>::operator+=(const SparseMatrix<value_type>& other)
{
  Assert(rows == other.n_rows() and
         cols == other.n_cols(),
         "Mismatched size error.");

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < other.m_indices[i].size(); ++j)
      this->set(i, other.m_indices[i][j], other.m_data[i][j]);
  return *this;
}



template<typename value_type>
inline SparseMatrix<value_type>
SparseMatrix<value_type>::
operator-(const SparseMatrix<value_type>& other) const
{
  Assert(rows == other.n_rows() and
         cols == other.n_cols(),
         "Mismatched size error.");

  SparseMatrix m(*this);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < other.m_indices[i].size(); ++j)
      m.set(i, other.m_indices[i][j], -other.m_data[i][j]);
  return m;
}



template<typename value_type>
inline SparseMatrix<value_type>&
SparseMatrix<value_type>::operator-=(const SparseMatrix<value_type>& other)
{
  Assert(rows == other.n_rows() and
         cols == other.n_cols(),
         "Mismatched size error.");

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < other.m_indices[i].size(); ++j)
      this->set(i, other.m_indices[i][j], -other.m_data[i][j]);
  return *this;
}



template<typename value_type>
inline Vector<value_type>
SparseMatrix<value_type>::operator*(const Vector<value_type>& x) const
{
  Assert(this->n_cols() == x.size(), "Mismatched size error.");

  Vector<value_type> y(x.size());
  for (size_t i = 0; i < m_indices.size(); ++i)
    for (size_t j = 0; j < m_indices[i].size(); ++j)
      y[i] += m_data[i][j] * x[m_indices[i][j]];
}



template<typename value_type>
inline SparseMatrix<value_type>
SparseMatrix<value_type>::operator*(const SparseMatrix<value_type>& other) const
{
  Assert(cols == other.n_rows(), "Mismatched size error.");

  SparseMatrix m(rows, other.n_cols());
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
        // Compute c_{ik} += a_{ij} * b_{jk}
        size_t k = other.m_indices[j][kr];
        value_type c_ik += a_ij * other.m_data[j][kr];
        m.set(i, k, c_ik);
      }
    }
  }
}

}
#endif //SPARSE_MATRIX_H

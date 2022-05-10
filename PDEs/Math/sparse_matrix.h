#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "vector.h"
#include "matrix.h"
#include "exceptions.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cinttypes>


namespace math
{

template<typename value_type>
class SparseMatrix
{
private:
  typedef std::vector<std::vector<uint64_t>> SparsityPattern;

private:
  uint64_t rows;
  uint64_t cols;

  SparsityPattern m_colnums;
  std::vector<std::vector<value_type>> m_data;


public:
  /** Default constructor. */
  SparseMatrix() : rows(0), cols(0) {}

  /** Construct a square sparse matrix with dimension \p n. */
  explicit SparseMatrix(const uint64_t n)
      : rows(n), cols(n), m_data(n), m_colnums(n)
  {}

  /** Construct a sparse matrix with \p n_rows and \p n_cols. */
  explicit SparseMatrix(const uint64_t n_rows, const uint64_t n_cols)
    : rows(n_rows), cols(n_cols), m_data(n_rows),
      m_colnums(n_rows)
  {}

  /** Construct from a sparsity pattern. */
  SparseMatrix(const SparsityPattern& sparsity_pattern)
    : rows(sparsity_pattern.size()), m_data(sparsity_pattern.size()),
      m_colnums(sparsity_pattern)
  {
    cols = 0;
    for (uint64_t i = 0; i < rows; ++i)
    {
      // Sort column indices for the row
      auto& colnums = m_colnums[i];
      std::stable_sort(colnums.begin(), colnums.end());

      // Resize the data vector for the row
      m_data[i].resize(m_colnums[i].size(), 0.0);

      // Check the largest column index
      for (const auto& j : m_colnums[i])
        if (j > cols) cols = j;
    }
  }

  /** Copy constructor. */
  SparseMatrix(const SparseMatrix& other)
    : rows(other.rows), cols(other.cols), m_data(other.m_data),
      m_colnums(other.m_colnums)
  {}

  /** Move constructor. */
  SparseMatrix(SparseMatrix&& other)
    : rows(other.rows), cols(other.cols),
      m_data(std::move(other.m_data)),
      m_colnums(std::move(other.m_colnums))
  {}

  /** Initialize with a dense matrix. */
  SparseMatrix(const Matrix<value_type>& other)
    : rows(other.n_rows()), cols(other.n_cols()),
      m_data(other.n_rows()),
      m_colnums(other.n_rows())
  {
    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t j = 0; j < cols; ++j)
        if (other[i][j] != 0.0)
        {
          m_data[i].push_back(other[i][j]);
          m_colnums[i].push_back(j);
        }
  }

  /** Copy assignment operator. */
  SparseMatrix& operator=(const SparseMatrix& other)
  {
    rows = other.rows;
    cols = other.cols;
    m_data = other.m_data;
    m_colnums = other.m_colnums;
    return *this;
  }

  /** Move assignment operator. */
  SparseMatrix& operator=(SparseMatrix&& other)
  {
    rows = other.rows;
    cols = other.cols;
    m_data = std::move(other.m_data);
    m_colnums = std::move(other.m_colnums);
    return *this;
  }

public:
  /** \name Accessors */
  /** @{ */

  /** Read access for the element at row \p i and column \p j. */
  value_type operator()(const uint64_t i, const uint64_t j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // If row is uninitialized, return zero
    if (m_colnums[i].empty())
      return 0.0;

    // Otherwise, look for the specified column in the row
    const auto& colnums = m_colnums[i];
    auto rel_loc = std::lower_bound(colnums.begin(), colnums.end(), j);

    // If the column is not in the row, return zero
    if (rel_loc == colnums.end())
      return 0.0;

    // Otherwise, return the element
    uint64_t jr = rel_loc - colnums.begin();
    return m_data[i][jr];
  }

  /** Read/write access for the element at row \p i and column \p j. */
  value_type& operator()(const uint64_t i, const uint64_t j)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // Check whether the row exists
    const auto& colnums = m_colnums[i];
    Assert(not colnums.empty(),
           "Invalid access attempt. Element not initialized.");

    // Check whether the column exists on the row
    auto rel_loc = std::lower_bound(colnums.begin(), colnums.end(), j);
    Assert(rel_loc != colnums.end(),
           "Invalid access attempt. Element not initialized.");

    uint64_t jr = rel_loc - colnums.begin();
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
    m_colnums.clear();
  }

  /** Reinitialize the sparse matrix with dimension \p n. */
  void reinit(const uint64_t n)
  {
    rows = cols = n;
    m_colnums.clear();
    m_data.clear();

    m_colnums.resize(n);
    m_data.resize(n);
  }

  /** Reinitialize the sparse matrix with \p n_rows and \p n_cols. */
  void reinit(const uint64_t n_rows, const uint64_t n_cols)
  {
    rows = n_rows; cols = n_cols;
    m_colnums.clear();
    m_data.clear();

    m_colnums.resize(n_rows);
    m_data.resize(n_rows);
  }

  /** Set the element at row \p i and column \p j to \p value. */
  void insert(const uint64_t i, const uint64_t j,
              const value_type value,
              const bool adding = true)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    /* If the row is empty or the column number is larger than all current
     * entries on the row, add to the back of the row. */
    if (m_colnums[i].size() == 0 or m_colnums[i].back() < j)
    {
      m_colnums[i].push_back(j);
      m_data[i].push_back(value);
      return;
    }

    // Find the index to insert which maintains sorting
    auto& colnums = m_colnums[i];
    auto jloc = std::lower_bound(colnums.begin(), colnums.end(), j);
    uint64_t jr = jloc - colnums.begin();

    // If this points to an existing column, add to it
    if (*jloc == j)
    {
      m_data[i][jr] = (adding)? m_data[i][jr] + value : value;
      return;
    }

    // Insert the entries into the data structures maintaining sorting
    m_colnums[i].insert(jloc, j);
    m_data[i].insert(m_data[i].begin() + jr, value);
  }

  /** Set a list of elements. See \ref insert.*/
  void insert(const std::vector<uint64_t>& row_indices,
              const std::vector<uint64_t>& col_indices,
              const std::vector<value_type>& values,
              const bool adding = true)
  {
    Assert(row_indices.size() == col_indices.size() &&
           row_indices.size() == values.size(),
           "All inputs must be of the same length.");

    for (uint64_t i = 0; i < row_indices.size(); ++i)
      insert(row_indices[i], col_indices[i], values[i], adding);
  }

  /** Swap the elements of two rows. */
  void swap_row(const uint64_t i0, const uint64_t i1)
  {
    Assert(i0 < rows && i1 < rows, "Invalid row indices provided.");
    m_colnums[i0].swap(m_colnums[i1]);
    m_data[i0].swap(m_data[i1]);
  }

  /** Swap the elements of this sparse matrix with another. */
  void swap(SparseMatrix& other)
  {
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
    m_colnums.swap(other.m_colnums);
    m_data.swap(other.m_data);
  }

  /** @} */
  /** \name Information */
  /** @{ */

  /** Return the number of rows the sparse matrix represents. */
  uint64_t n_rows() const
  {
    return rows;
  }

  /** Return the number of columns the sparse matrix represents. */
  uint64_t n_cols() const
  {
    return cols;
  }

  /** Return the number of nonzero elements in the sparse matrix. */
  uint64_t nnz() const
  {
    uint64_t nnz = 0;
    for (const auto& colnums : m_colnums)
      nnz += colnums.size();
    return nnz;
  }

  /** Returns whether an element exists or not. */
  bool exists(const uint64_t i, const uint64_t j)
  {
    return std::binary_search(m_colnums[i].begin(), m_colnums[i].end(), j);
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
    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t jr = 0; jr < other.m_colnums[i].size(); ++jr)
        A.insert(i, other.m_colnums[i][jr], other.m_data[i][jr]);
    return A;
  }

  /** Element-wise addition of two sparse matrices in-place. */
  SparseMatrix& operator+=(const SparseMatrix& other)
  {
    Assert(rows == other.n_rows() &&
           cols == other.n_cols(),
           "Dimension mismatch error.");

    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t jr = 0; jr < other.m_colnums[i].size(); ++jr)
        this->insert(i, other.m_colnums[i][jr], other.m_data[i][jr]);
    return *this;
  }

  /** Element-wise subtraction of two sparse matrices. */
  SparseMatrix operator-(const SparseMatrix& other) const
  {
    Assert(rows == other.n_rows() &&
           cols = other.n_cols(),
           "Dimension mismatch error.");

    SparseMatrix A(*this);
    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t jr = 0; jr < other.m_colnums[i].size(); ++jr)
        A.insert(i, other.m_colnums[i][jr], -other.m_data[i][jr]);
    return A;
  }

  /** Element-wise subtraction of two sparse matrices in-place. */
  SparseMatrix& operator-=(const SparseMatrix& other)
  {
    Assert(rows == other.n_rows() &&
           cols == other.n_cols(),
           "Dimension mismatch error.");

    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t jr = 0; jr < other.m_colnums[i].size(); ++jr)
        this->insert(i, other.m_colnums[i][jr], -other.m_data[i][jr]);
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
    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t jr = 0; jr < m_colnums[i].size(); ++jr)
        Ax[i] += m_data[i][jr] * x[m_colnums[i][jr]];
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
    for (uint64_t i = 0; i < rows; ++i)
    {
      // Loop over relative column indices for row i of this
      for (uint64_t jr = 0; jr < m_colnums[i].size(); ++jr)
      {
        uint64_t j = m_colnums[i][jr];
        value_type a_ij = m_data[i][jr];

        // Loop over relative column indices of row j of other
        for (uint64_t kr = 0; kr < other.m_colnums[j].size(); ++kr)
        {
          // Compute c_{ik} += a_{ij} b_{jk}
          uint64_t k = other.m_colnums[j][kr];
          value_type c_ik = a_ij * other.m_data[j][kr];
          A.insert(i, k, c_ik);
        }//for columns in row j of other matrix
      }//for columns in row i of this matrix
    }//for row
  }

  std::string to_string() const
  {
    std::stringstream ss;
    auto fill = std::cout.fill(' ');
    ss << std::left << std::setw(8)
       << "Row"    << fill << std::setw(8)
       << "Column" << fill << std::setw(10)
       << std::right << "Value" << fill << "\n";
    ss << std::left << std::setw(8)
       << "---" << fill << std::setw(8)
       << "------" << fill << std::setw(10)
       << std::right << "-----" << fill << "\n";

    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t jr = 0; jr < m_colnums[i].size(); ++jr)
        ss << std::left << std::setw(8)
           << i << fill << std::setw(8)
           << m_colnums[i][jr] << fill << std::setw(10)
           << std::right << m_data[i][jr] << fill << "\n";
    return ss.str();
  }

  void print() const
  {
    std::cout << to_string();
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

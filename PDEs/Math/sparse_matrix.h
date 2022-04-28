#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

class SparseMatrix
{
private:
  size_t rows;
  size_t cols;

  std::vector<std::vector<double>> m_data;
  std::vector<std::vector<size_t>> m_indices;

public:
  /// Default constructor.
  SparseMatrix() = default;

  /// Construct a sparse matrix with dimension \p n.
  explicit SparseMatrix(const size_t n)
    : rows(n), cols(n), m_data(n), m_indices(n)
  {}
  /// Construct a sparse matrix with \p n_rows rows and \p n_cols columns.
  explicit SparseMatrix(const size_t n_rows, const size_t n_cols)
    : rows(n_rows), cols(n_cols), m_data(n_rows), m_indices(n_rows)
  {}

  /// Construct a sparsity pattern.
  SparseMatrix(const std::vector<std::vector<size_t>>& pattern)
    : rows(pattern.size()), m_data(pattern.size()), m_indices(pattern)
  {
    for (size_t i = 0; i < rows; ++i)
    {
      cols = 0;
      m_data[i].resize(m_indices[i].size(), 0.0);
      for (const auto& j : m_indices[i])
        if (j > cols) cols = j;
    }
  }

  /// Copy constructor.
  SparseMatrix(const SparseMatrix& other)
    : rows(other.rows), cols(other.cols),
      m_data(other.m_data), m_indices(other.m_indices)
  {}
  /// Move constructor.
  SparseMatrix(SparseMatrix&& other)
    : rows(other.rows), cols(other.cols),
      m_data(std::move(other.m_data)),
      m_indices(std::move(other.m_indices))
  {}

  /// Copy assignment operator.
  SparseMatrix& operator=(const SparseMatrix& other)
  {
    rows = other.rows;
    cols = other.cols;
    m_data = other.m_data;
    m_indices = other.m_indices;
    return *this;
  }
  /// Move assignment operator.
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

  inline double operator()(const size_t i, const size_t j) const;
  inline double& operator()(const size_t i, const size_t j);

  /** @} */
  /** \name Modifiers */
  /** @{ */

  inline void clear();

  inline void set(const size_t i, const size_t j, const double value);
  inline void set(const std::vector<size_t>& row_indices,
                  const std::vector<size_t>& col_indices,
                  const std::vector<double>& values);
  inline void set(const size_t i, const std::vector<size_t> col_indices,
                  const std::vector<double> values);

  inline reinit(std::vector<size_t> prealloc);
  inline void compress();

  /** @} */
  /** \name Information */
  /** @{ */

  size_t n_rows() const { return rows; }
  size_t n_cols() const { return cols; }
  size_t size() const
  {
    size_t count = 0;
    for (const auto& row : m_data)
      count += row.size();
    return count;
  }

  /** @} */

private:
  inline void validate_indices(const size_t i, const size_t j,
                               const std::string func_name) const
  {
    if (i > rows or j > cols)
    {
      std::stringstream err;
      err << "SparseMatrix::" << func_name << ": "
          << "The provided indices exceed the dimensions of the matrix.";
      throw std::length_error(err.str());
    }
  }
};

/*-------------------- Inline Implementations --------------------*/

/// Clear all data within the sparse matrix.
inline void SparseMatrix::clear()
{
  rows = 0; cols = 0;
  m_data.clear();
  m_indices.clear();
}

/// Read access to element <tt>(i, j)</tt>.
inline double
SparseMatrix::operator()(const size_t i, const size_t j) const
{
  this->validate_indices(i, j, __FUNCTION__);

  // Check if index exists in structure
  double value = 0.0;
  if (not m_data[i].empty())
  {
    auto rel_loc = std::find(m_data[i].begin(), m_data[i].end(), j);
    bool elem_exists = (rel_loc != m_data[i].end());
    if (elem_exists)
    {
      size_t jr = rel_loc - m_data[i].begin();
      value = m_data[i][jr];
    }
  }
  return value;
}

/// Read/write access to element <tt>(i, j)</tt>.
inline double&
SparseMatrix::operator()(const size_t i, const size_t j)
{
  // Are indices valid?
  this->validate_indices(i, j, __FUNCTION__);

  // Does the row have nonzero entries?
  if (m_data[i].empty())
  {
    std::stringstream err;
    err << "SparseMatrix::" << __FUNCTION__ << ": "
        << "Cannot access an unset element. Row " << i << " is empty.";
    throw std::length_error(err.str());
  }

  // Has element (i, j) been set?
  auto rel_loc = std::find(m_data[i].begin(), m_data[i].end(), j);
  bool elem_exists = (rel_loc != m_data[i].end());
  if (not elem_exists)
  {
    std::stringstream err;
    err << "SparseMatrix::" << __FUNCTION__ << ": "
        << "Cannot access an unset element. "
        << "Element (" << i << "," << j << ") has not been set.";
    throw std::length_error(err.str());
  }

  // Return the element
  size_t jr = rel_loc - m_data[i].begin();
  return m_data[i][jr];
}

/// Set the element <tt>(i, j)</tt>.
inline void
SparseMatrix::set(const size_t i, const size_t j, const double value)
{
  this->validate_indices(i, j, __FUNCTION__);

  auto rel_loc = std::find(m_data[i].begin(), m_data[i].end(), j);
  bool elem_exists = (rel_loc != m_data[i].end());
  if (elem_exists)
  {
    size_t jr = rel_loc - m_data[i].begin();
    m_data[i][jr] += value;
  }
  else
  {
    m_indices[i].push_back(j);
    m_data[i].push_back(value);
  }

}

/// Set a number of elements based on a set of row and column indices.
inline void
SparseMatrix::set(const std::vector<size_t>& row_indices,
                  const std::vector<size_t>& col_indices,
                  const std::vector<double>& values)
{
  if (row_indices.size() != col_indices.size() or
      row_indices.size() != values.size())
  {
    std::stringstream err;
    err << "SparseMatrix::" << __FUNCTION__ << ": "
        << "There must be the same number of row indices, column indices, "
        << "and values.";
    throw std::runtime_error(err.str());
  }

  for (size_t e = 0; e < values.size(); ++e)
    this->set(row_indices[e], col_indices[e], values[e]);
}

/// Set a number of column elements on a single row.
inline void
SparseMatrix::set(const size_t i,
                  const std::vector<size_t> col_indices,
                  const std::vector<double> values)
{
  if (col_indices.size() != values.size())
  {
    std::stringstream err;
    err << "SparseMatrix::" << __FUNCTION__ << ": "
        << "There must be the same number of column indices, "
        << "and values.";
    throw std::runtime_error(err.str());
  }

  for (size_t j = 0; j < values.size(); ++j)
    this->set(i, col_indices[j], values[j]);
}

/// Reinitialize the matrix with a preallocation vector.
inline void
SparseMatrix::reinit(std::vector<size_t> prealloc)
{
  this->clear();
  rows = prealloc.size();
  m_data.resize(rows);
  m_indices.resize(rows);

  for (size_t i = 0; i < rows; ++i)
  {
    m_data[i].resize(prealloc[i]);
    m_indices[i].resize(prealloc[i]);
    for (const auto& j : prealloc[i])
      if (j > cols) cols = j;
  }
}


#endif //SPARSE_MATRIX_H

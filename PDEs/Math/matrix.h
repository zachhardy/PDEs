#ifndef MATRIX_H
#define MATRIX_H

#include "Math/vector.h"

#include <cmath>
#include <vector>
#include <iomanip>

namespace math
{

template<typename value_type>
class Matrix
{
private:
  std::vector<std::vector<value_type>> m_data;

public:
  /** Default constructor. */
  Matrix() = default;

  /** Construct a square matrix of dimension \p n set to \p value. */
  explicit Matrix(const size_t n, const value_type value = 0.0)
      : m_data(n, std::vector<value_type>(n, value))
  {}

  /** Construct a matrix with \p n_rows and \p n_cols set to \p value. */
  explicit Matrix(const size_t n_rows,
                  const size_t n_cols,
                  const value_type value = 0.0)
      : m_data(n_rows, std::vector<value_type>(n_cols, value))
  {}

  /** Copy constructor. */
  Matrix(const Matrix& other)
      : m_data(other.m_data)
  {}

  /** Move constructor. */
  Matrix(Matrix&& other)
      : m_data(std::move(other.m_data))
  {}

  /** Copy construction from an STL vector. */
  Matrix(const std::vector<std::vector<value_type>>& other)
  {
    check_stl(other, __PRETTY_FUNCTION__);
    m_data = other;
  }

  /** Move construction from an STL vector. */
  Matrix(std::vector<std::vector<value_type>>&& other)
  {
    check_stl(other, __PRETTY_FUNCTION__);
    m_data = std::move(other);
  }

  /** Construction from a nested initializer list. */
  Matrix(std::initializer_list<std::initializer_list<value_type>>& list)
  {
    m_data.clear();
    for (auto& row : list)
      m_data.emplace_back(row);
    check_stl(m_data, __PRETTY_FUNCTION__);
  }

  /** Copy assignment operator. */
  Matrix& operator=(const Matrix& other)
  {
    m_data = other.m_data;
    return *this;
  }

  /** Move assignment operator. */
  Matrix& operator=(Matrix&& other)
  {
    m_data = std::move(other.m_data);
    return *this;
  }

  /** Copy assignment from an STL vector. */
  Matrix& operator=(const std::vector<std::vector<value_type>>& other)
  {
    check_stl(other, __PRETTY_FUNCTION__);
    m_data = other;
    return *this;
  }

  /** Move assignment from an STL vector. */
  Matrix& operator=(std::vector<std::vector<value_type>>&& other)
  {
    check_stl(other, __PRETTY_FUNCTION__);
    m_data = std::move(other);
    return *this;
  }

public:
  /** \name Accessors */
  /** @{ */

  /** Read/write access for row \p i. */
  std::vector<value_type>& operator[](const size_t i)
  {
    return m_data[i];
  }

  /** Read only access for row \p i. */
  std::vector<value_type> operator[](const size_t i) const
  {
    return m_data[i];
  }

  /** Read/write access for row \p i. */
  std::vector<value_type>& operator()(const size_t i)
  {
    return m_data[i];
  }

  /** Read only access for row \p i. */
  std::vector<value_type> operator()(const size_t i) const
  {
    return m_data[i];
  }

  /** Read/write access for row \p i with bounds checking. */
  std::vector<value_type>& at(const size_t i)
  {
    return m_data.at(i);
  }

  /** Read only access for row \p i with bounds checking. */
  std::vector<value_type> at(const size_t i) const
  {
    return m_data.at(i);
  }

  /** Read/write access for row \p i and column \p j. */
  value_type& operator()(const size_t i, const size_t j)
  {
    return m_data[i][j];
  }

  /** Read only access for row \p i and column \p j. */
  value_type operator()(const size_t i, const size_t j) const
  {
    return m_data[i][j];
  }

  /** Read/write access for row \p i and column \p j with bounds checking. */
  value_type& at(const size_t i, const size_t j)
  {
    return m_data.at(i).at(j);
  }

  /** Read only access for row \p i and column \p j with bounds checking. */
  value_type at(const size_t i, const size_t j) const
  {
    return m_data.at(i).at(j);
  }

  /** Read/write access to the <tt>i</tt>'th diagonal element. */
  value_type& diagonal(const size_t i)
  {
    return m_data[i][i];
  }

  /** Read access to the <tt>i</tt>'th diagonal element. */
  value_type diagonal(const size_t i) const
  {
    return m_data[i][i];
  }

  /** Return the diagonal of the matrix. */
  std::vector<value_type> diagonal() const
  {
    if (m_data.empty())
      return std::vector<value_type>();

    // Compute diagonal dimension = min(n_rows, n_cols)
    size_t min_dim = std::min(n_rows(), n_cols());

    // Populate the diagonal vector
    std::vector<value_type> diag;
    for (size_t i = 0; i < min_dim; ++i)
      diag.push_back(m_data[i][i]);
    return diag;
  }

  /** Access the underlying matrix data. */
  std::vector<value_type>* data()
  {
    return m_data.data();
  }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /** Clear the matrix. */
  void clear()
  {
    m_data.clear();
  }

  /** Remove the last row of the matrix. */
  void pop_back()
  {
    m_data.pop_back();
  }

  /** Add a new row from \p values to the back of the matrix. */
  void push_back(const std::vector<value_type>& values)
  {
    m_data.push_back(values);
  }

  /** Add a new row from \p values in place to the back. */
  void emplace_back(const std::vector<value_type>& values)
  {
    m_data.emplace_back(values);
  }

  /** Resize to dimension \p n, setting new elements to \p value. */
  void resize(const size_t n, const value_type value = 0.0)
  {
    m_data.resize(n, std::vector<value_type>(n, value));
  }

  /** Resize to \p n_rows and \p n_cols, setting new elements to \p value. */
  void resize(const size_t n_rows, const size_t n_cols,
              const value_type value = 0.0)
  {
    m_data.resize(n_rows, std::vector<value_type>(n_cols, value));
  }

  /** Swap the elements of two rows. */
  void swap_row(const size_t i0, const size_t i1)
  {
    Assert(i0 < n_rows() && i1 < n_rows(),
           "Invalid row indices provided.");
    m_data[i0].swap(m_data[i1]);
  }

  /** Swap the elements of two columns. */
  void swap_column(const size_t j0, const size_t j1)
  {
    Assert(j0 < n_cols() && j1 < n_cols(),
           "Invalid column indices provided.");
    for (size_t i = 0; i < n_rows(); ++i)
      std::swap(m_data[i][j0], m_data[i][j1]);
  }

  /** Swap the elements of this matrix with another. */
  void swap(Matrix& other)
  {
    m_data.swap(other.m_data);
  }

  /** Set the diagonal of the matrix. */
  void set_diagonal(const Vector<value_type>& diag)
  {
    // Set the matrix from a diagonal
    if (m_data.empty())
    {
      m_data.resize(diag.size(), std::vector<value_type>(diag.size()));
      for (size_t i = 0; i < diag.size(); ++i)
        m_data[i][i] = diag[i];
    }

      // Override the diagonal of an initialized matrix
    else
    {
      size_t min_dim = std::min(n_rows(), n_cols());
      Assert(diag.size() == min_dim, "Dimension mismatch error.");
      for (size_t i = 0; i < min_dim; ++i)
        m_data[i][i] = diag[i];
    }
  }

  /** Set the diagonal of the matrix with an STL vector. */
  void set_diagonal(const std::vector<value_type>& diag)
  {
    set_diagonal(Vector(diag));
  }

  /** Set the diagonal of the matrix with a fixed scalar value. */
  void set_diagonal(const value_type value)
  {
    if (m_data.empty())
      m_data.resize(1, std::vector<value_type>(1, value));

    size_t min_dim = std::min(n_rows(), n_cols());
    for (size_t i = 0; i < min_dim; ++i)
      m_data[i][i] = value;
  }

  /** @} */
  /** \name Information */
  /** @{ */

  /** Return the number of rows. */
  size_t n_rows() const
  {
    return m_data.size();
  }

  /** Return the number of columns. */
  size_t n_cols() const
  {
    return m_data.front().size();
  }

  /** Return the number of elements <tt> n_rows * n_cols </tt>. */
  size_t size() const
  {
    return m_data.size() * m_data.front().size();
  }

  /** Return the number of non-zeros. */
  size_t nnz() const
  {
    size_t nnz = 0;
    for (const auto& row : m_data)
      for (const auto& elem : row)
        if (elem != 0.0) ++nnz;
    return nnz;
  }

  /** Return whether the matrix is empty. */
  bool empty() const noexcept
  {
    return m_data.empty();
  }

  /** @} */
  /** \name Iterators */
  /** @{ */

  /** Mutable iterator over rows at the start of the matrix. */
  typename std::vector<std::vector<value_type>>::iterator
  begin()
  {
    return m_data.begin();
  }

  /** Mutable iterator one past the last row of the matrix. */
  typename std::vector<std::vector<value_type>>::iterator
  end()
  {
    return m_data.end();
  }

  /** Mutable iterator at the start of row \p i. */
  typename std::vector<value_type>::iterator
  begin(const size_t i)
  {
    return m_data[i].begin();
  }

  /** Mutable iterator one past the end of row \p i. */
  typename std::vector<value_type>::iterator
  end(const size_t i)
  {
    return m_data[i].end();
  }

  /** Constant iterator over rows at the start of the matrix. */
  typename std::vector<std::vector<value_type>>::const_iterator
  cbegin() const
  {
    return m_data.cbegin();
  }

  /** Constant iterator one past the last row of the matrix. */
  typename std::vector<std::vector<value_type>>::const_iterator
  cend() const
  {
    return m_data.cend();
  }

  /** Constant iterator at the start of the row \p i. */
  typename std::vector<value_type>::const_iterator
  cbegin(const size_t i) const
  {
    return m_data[i].cbegin();
  }

  /** Constant iterator one past the end of row \p i. */
  typename std::vector<value_type>::const_iterator
  cend(const size_t i) const
  {
    return m_data[i].cend();
  }

  /** @} */
  /** \name Scalar Operations */
  /** @{ */

  /** Element-wise negation. */
  Matrix operator-() const
  {
    Matrix dst(m_data);
    for (auto& row : dst)
      for (auto& elem : row)
        elem = -elem;
    return dst;
  }

  /** Element-wise negation in-place. */
  Matrix& operator-()
  {
    for (auto& row : m_data)
      for (auto& elem : row)
        elem = -elem;
    return *this;
  }

  /** Element-wise multiplication by a scalar. */
  Matrix operator*(const value_type value) const
  {
    Matrix dst(m_data);
    for (auto& row : dst)
      for (auto& elem : row)
        elem *= value;
    return dst;
  }

  /** Element-wise multiplication by a scalar in-place. */
  Matrix& operator*=(const value_type value)
  {
    for (auto& row : m_data)
      for (auto& elem : row)
        elem *= value;
    return *this;
  }

  /** Element-wise division by a scalar. */
  Matrix operator/(const value_type value) const
  {
    Assert(value != 0.0, "Division by zero error.");

    Matrix dst(m_data);
    for (auto& row : dst)
      for (auto& elem : row)
        elem /= value;
    return dst;
  }

  /** Element-wise division by a scalar in-place. */
  Matrix& operator/=(const value_type value)
  {
    Assert(value != 0.0, "Division by zero error.");

    for (auto& row : m_data)
      for (auto& elem : row)
        elem /= value;
    return *this;
  }

  /** @} */
  /** \name Linear Algebra Operations */
  /** @{ */

  /** Element-wise addition of two matrices. */
  Matrix operator+(const Matrix& other) const
  {
    Assert(n_rows() == other.n_rows() &&
           n_cols() == other.n_cols(),
           "Dimension mismatch error.");

    Matrix dst(m_data);
    for (size_t i = 0; i < dst.n_rows(); ++i)
      for (size_t j = 0; j < dst.n_cols(); ++j)
        dst[i][j] += other[i][j];
    return dst;
  }

  /** Element-wise addition of two matrices in-place. */
  Matrix& operator+=(const Matrix& other)
  {
    Assert(n_rows() == other.n_rows() &&
           n_cols() == other.n_cols(),
           "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
      for (size_t j = 0; j < n_cols(); ++j)
        m_data[i][j] += other[i][j];
    return *this;
  }

  /** Element-wise subtraction of two matrices. */
  Matrix operator-(const Matrix& other) const
  {
    Assert(n_rows() == other.n_rows() and
           n_rows() == other.n_cols(),
           "Dimension mismatch error.");

    Matrix dst(m_data);
    for (size_t i = 0; i < dst.n_rows(); ++i)
      for (size_t j = 0; j < dst.n_cols(); ++j)
        dst[i][j] -= other[i][j];
    return dst;
  }

  /** Element-wise subtraction of two matrices in-place. */
  Matrix& operator-=(const Matrix& other)
  {
    Assert(n_rows() == other.n_rows() and
           n_rows() == other.n_cols(),
           "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
      for (size_t j = 0; j < n_cols(); ++j)
        m_data[i][j] -= other[i][j];
    return *this;
  }

  /**
   * Return the matrix-matrix product.
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B} \\
   *     c_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  Matrix operator*(const Matrix& other) const
  {
    Assert(n_cols() == other.n_rows(), "Dimension mismatch error.");

    Matrix dst(n_rows(), other.n_cols());
    for (size_t i = 0; i < dst.n_rows(); ++i)
    {
      for (size_t j = 0; j < dst.n_cols(); ++j)
      {
        value_type c_ij = 0.0;
        for (size_t k = 0; k < dst.n_rows(); ++k)
          c_ij += m_data[i][k] * other[k][j];
        dst[i][j] += c_ij;
      }
    }
    return dst;
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
    Assert(x.size() == n_cols(), "Dimension mismatch error.");

    Vector<value_type> dst(n_rows());
    for (size_t i = 0; i < n_rows(); ++i)
      for (size_t j = 0; j < n_cols(); ++j)
        dst[i] += m_data[i][j] * x[j];
    return dst;
  }

  /**
   * Return the transpose of a matrix.
   * This is computed via
   * \f[ \boldsymbol{B} = \boldsymbol{A}^T \\
   *     b_{ij} = a_{ji}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  Matrix transpose() const
  {
    Matrix dst(n_cols(), n_rows());
    for (size_t i = 0; i < dst.n_rows(); ++i)
      for (size_t j = 0; j < dst.n_cols(); ++j)
        dst[i][j] = m_data[j][i];
    return dst;
  }

  void vmult(const Vector<value_type>& x,
             Vector<value_type>& dst)
  {
    Assert(n_cols() == x.size(), "Dimension mismatch error.");
    Assert(n_rows() == dst.size(), "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
    {
      value_type value = 0.0;
      for (size_t j = 0; j < n_cols(); ++j)
        value += m_data[i][j] * x[j];
      dst[i] = value;
    }

  }

  /** @} */
  /** \name Print Utilities */
  /** @{ */

  /** Return the matrix as a string. */
  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < m_data.size(); ++i)
    {
      ss << ((i == 0)? "[" : " [");
      for (size_t j = 0; j < m_data[i].size() - 1; ++j)
        ss << std::setw(10) << std::setprecision(6) << m_data[i][j] << " ";
      ss << std::setw(10) << std::setprecision(6) << m_data[i].back() << "]"
         << ((i == m_data.size()-1)? "]\n" : "\n");
    }
    return ss.str();
  }

  /** Print the vector to `std::cout`. */
  void print() const
  {
    std::cout << to_string();
  }

  /** @} */
private:
  static void check_stl(const std::vector<std::vector<value_type>>& A,
                        const std::string func_name)
  {
    bool is_valid = true;
    auto m = A.front().size();
    for (size_t i = 0; i < A.size(); ++i)
      if (A[i].size() != m) { is_valid = false; break; }

    Assert(is_valid,
           "Dimension mismatch error. "
           "All rows must be the same size.");
  }
};


/*-------------------- Inline Implementations --------------------*/


/** Element-wise multiplication by a scalar value. */
template<typename value_type>
inline Matrix<value_type>
operator*(const value_type value, const Matrix<value_type>& A)
{
  return A * value;
}

}
#endif //MATRIX_H

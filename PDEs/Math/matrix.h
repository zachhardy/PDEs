#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"
#include "exceptions.h"

#include <cmath>
#include <vector>
#include <iomanip>
#include <cinttypes>


namespace math
{


template<typename number>
class Matrix
{
public:
  using row_iterator = typename Vector<number>::iterator ;
  using const_row_iterator = typename Vector<number>::const_iterator;

  using iterator = typename std::vector<Vector<number>>::iterator;
  using const_iterator = typename std::vector<Vector<number>>::const_iterator;
  
private:
  std::vector<Vector<number>> m_data;

public:
  /** Default constructor. */
  Matrix() = default;

  /** Construct a matrix with \p n_rows and \p n_cols set to \p value. */
  explicit Matrix(const uint64_t n_rows,
                  const uint64_t n_cols,
                  const number value = 0.0)
    : m_data(n_rows, Vector<number>(n_cols, value))
  {}

  /** Construct a square matrix of dimension \p n set to \p value. */
  explicit Matrix(const uint64_t n, const number value)
    : Matrix(n, n, value)
  {}

  /** Copy constructor. */
  Matrix(const Matrix& other) : m_data(other.m_data) {}

  /** Move constructor. */
  Matrix(Matrix&& other) : m_data(std::move(other.m_data)) {}

  /** Copy construction from an STL vector. */
  Matrix(const std::vector<std::vector<number>>& other)
  {
    Assert(valid_stl_matrix(other),
           "Invalid input. Ensure all rows are the same length.");
    m_data = other;
  }

  /** Move construction from an STL vector. */
  Matrix(std::vector<std::vector<number>>&& other)
  {
    Assert(valid_stl_matrix(other),
           "Invalid input. Ensure all rows are the same length.");
    m_data = std::move(other);
  }

  /** Construction from a nested initializer list. */
  Matrix(std::initializer_list<std::initializer_list<number>>& init_list)
  {
    m_data.clear();
    for (auto& row : init_list)
      m_data.push_back(row);

    uint64_t m = m_data.front().size();
    for (const auto& row : m_data)
      Assert(row.size() == m,
             "Invalid input. Ensure all rows are the same length.");
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
  Matrix& operator=(const std::vector<std::vector<number>>& other)
  {
    Assert(valid_stl_matrix(other),
           "Invalid STL input. Ensure all rows are the same length.");
    m_data.clear();
    for (const auto& row : other)
      m_data.push_back(row);
  }

  /** Move assignment from an STL vector. */
  Matrix& operator=(std::vector<std::vector<number>>&& other)
  {
    Assert(valid_stl_matrix(other),
           "Invalid STL input. Ensure all rows are the same length.");
    m_data.clear();
    for (auto& row : other)
      m_data.push_back(std::move(row));
  }

  //###########################################################################
  /** \name Information */
  // @{

  /** Return the number of rows. */
  uint64_t n_rows() const { return m_data.size(); }

  /** Return the number of columns. */
  uint64_t n_cols() const { return m_data.front().size(); }

  /** Return the number of elements in the matrix. */
  uint64_t size() const { return n_rows() * n_cols(); }

  /** Return the number of nonzero elements in the matrix. */
  uint64_t nnz() const
  {
    uint64_t count = 0;
    for (const auto& row : m_data)
      for (const auto& el : row)
        if (el != 0.0) ++count;
    return count;
  }

  /** Return whether the matrix is empty. */
  bool empty() const { return m_data.empty(); }


  //###########################################################################
  /** \name Iterators */
  // @{

  /** Mutable iterator to the start of row \p i. */
  row_iterator begin(const uint64_t i)  { return m_data[i].begin(); }

  /** Mutable iterator to the end of row \p i. */
  row_iterator end(const uint64_t i) { return m_data[i].end(); }

  /** Mutable iterator to the first row of the matrix. */
  iterator begin() { return m_data.begin(); }

  /** Mutable iterator to the end of the matrix. */
  iterator end() { return m_data.end(); }

  /** Constant iterator to the start of row \p i. */
  const_row_iterator begin(const uint64_t i) const { return m_data[i].cbegin(); }

  /** Constant iterator to the end of row \p i. */
  const_row_iterator end(const uint64_t i) const { return m_data[i].cend(); }

  /** Constant iterator to the first row of the matrix. */
  const_iterator begin() const { return m_data.cbegin(); }

  /** Constant iterator to the end of the matrix. */
  const_iterator end() const { return m_data.cend(); }

  // @}

  //###########################################################################
  /** \name Accessors */
  // @{

  /** Read access for row \p i. */
  Vector<number> operator[](const uint64_t i) const { return m_data[i]; }

  /** Read/write access for row \p i. */
  Vector<number>& operator[](const uint64_t i) { return m_data[i]; }

  /** Read access for row \p i. */
  Vector<number> operator()(const uint64_t i) const { return m_data[i]; }

  /** Read/write access for row \p i. */
  Vector<number>& operator()(const uint64_t i) { return m_data[i]; }

  /** Read access for row \p i with bounds checking. */
  Vector<number> at(const uint64_t i) const { return m_data.at(i); }

  /** Read/write access for row \p i with bounds checking. */
  Vector<number>& at(const uint64_t i) { return m_data.at(i); }

  /** Read access for element \p(i, \p j). */
  number operator()(const uint64_t i, const uint64_t j) const { return m_data[i][j]; }

  /** Read/write access for element \p(i, \p j). */
  number& operator()(const uint64_t i, const uint64_t j) { return m_data[i][j]; }

  /** Read access for element \p(i, \p j) with bounds checking. */
  number at(const uint64_t i, const uint64_t j) const { return m_data.at(i).at(j); }

  /** Read/write access for element \p(i, \p j) with bounds checking. */
  number& at(const uint64_t i, const uint64_t j) { return m_data.at(i).at(j); }

  /** Read access to the <tt>i</tt>'th diagonal. */
  number diagonal(const uint64_t i) const { return m_data.at(i).at(i); }

  /** Read/write access to the <tt>i</tt>'th diagonal. */
  number& diagonal(const uint64_t i) { return m_data.at(i).at(i); }

  /** Return the diagonal of the matrix. */
  Vector<number> diagonal() const
  {
    std::vector<number> diag;
    uint64_t min_dim = std::min(n_rows(), n_cols());
    for (uint64_t i = 0; i < min_dim; ++i)
      diag.push_back(m_data[i][i]);
    return diag;
  }

  /** Return the underlying matrix m_data. */
  Vector<number>* data() { return m_data.data(); }

  /** Return the underlying m_data for row \p i. */
  number* row_data(const uint64_t i) { return m_data[i].data(); }

  // @}

  //###########################################################################
  /** \name Modifiers */
  // @{

  /** Clear the matrix m_data. */
  void clear() { m_data.clear(); }

  /** Remove the last row of the matrix. */
  void pop_back() { m_data.pop_back(); }

  /** Add a new row to the back of the matrix. */
  void push_back(const Vector<number>& row) { m_data.push_back(row); }

  void push_back(const std::vector<number>& row) { m_data.push_back(row); }

  /** Add a new row in-place to the back of the matrix. */
  void emplace_back(Vector<number>&& row) { m_data.emplace_back(row); }

  void emplace_back(std::vector<number>&& row) { m_data.emplace_back(row); }


  /** Resize to \p n_rows and \p n_cols, setting new elements to \p value. */
  void resize(const uint64_t n_rows, const uint64_t n_cols,
              const number value = 0.0)
  { m_data.resize(n_rows, Vector<number>(n_cols, value)); }

  /** Resize to dimension \p n, setting new elements to \p value. */
  void resize(const uint64_t n, const number value = 0.0)
  { resize(n, n, value); }

  /** Swap the elements of two rows. */
  void swap_row(const uint64_t i0, const uint64_t i1)
  {
    Assert(i0 < n_rows() && i1 < n_rows(), "Out of range error.");
    m_data[i0].swap(m_data[i1]);
  }

  /** Swap the elements of two columns. */
  void swap_column(const uint64_t j0, const uint64_t j1)
  {
    Assert(j0 < n_cols() && j1 < n_cols(), "Out of range error.");
    for (uint64_t i = 0; i < n_rows(); ++i)
      std::swap(m_data[i][j0], m_data[i][j1]);
  }

  /** Swap the elements of this matrix with another. */
  void swap(Matrix& other) { m_data.swap(other.m_data); }

  /** Set the diagonal of a matrix. */
  void set_diagonal(const Vector<number>& diag)
  {
    if (m_data.empty())
    {
      resize(diag.size());
      for (uint64_t i = 0; i < diag.size(); ++i)
        m_data[i][i] = diag[i];
    }
    else
    {
      uint64_t min_dim = std::min(n_rows(), n_cols());
      Assert(diag.size() == min_dim, "Dimension mismatch error.");
      for (uint64_t i = 0; i < min_dim; ++i)
        m_data[i][i] = diag[i];
    }
  }

  void set_diagonal(const std::vector<number>& diag)
  { set_diagonal(Vector<number>(diag)); }

  /** Set the diagonal with a fixed scalar value. */
  void set_diagonal(const number value)
  {
    if (m_data.empty()) {
      resize(1, value);
      return;
    }
    else
    {
      uint64_t min_dim = std::min(n_rows(), n_cols());
      for (uint64_t i = 0; i < min_dim; ++i)
        m_data[i][i] = value;
    }
  }

  // @}

  //###########################################################################
  /** \name Scalar Operations */
  // @{

  /** Element-wise negation in-place. */
  Matrix& operator-()
  {
    for (auto& row : m_data)
      for (auto& entry : row)
        entry = -entry;
    return *this;
  }

  /** Element-wise negation */
  Matrix operator-() const
  {
    Matrix A = *this;
    return -A;
  }

  /** Element-wise multiplication by a scalar in-place. */
  Matrix& operator*=(const number value)
  {
    for (auto& row : *this)
      for (auto& entry : row)
        entry *= value;
    return *this;
  }

  /** Element-wise multiplication by a scalar. */
  Matrix operator*(const number value) const
  {
    Matrix A = *this;
    A *= value;
    return A;
  }

  /** Element-wise division by a scalar in-place. */
  Matrix& operator/=(const number value)
  {
    Assert(value != 0.0, "Division by zero error.");
    for (auto& row : *this)
      for (auto& entry : row)
        entry /= value;
    return *this;
  }

  /** Element-wise division by a scalar. */
  Matrix operator/(const number value) const
  {
    Matrix A = *this;
    A /= value;
    return A;
  }

  //@}

  //###########################################################################
  /** \name Linear Algebra */
  // @{

  /** Element-wise addition of two matrices in-place. */
  Matrix& operator+=(const Matrix& other)
  {
    Assert(n_rows() == other.n_rows() &&
           n_cols() == other.n_cols(),
           "Dimension mismatch error.");
    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < n_cols(); ++j)
        m_data[i][j] += other[i][j];
    return *this;
  }

  /** Element-wise addition of two matrices. */
  Matrix operator+(const Matrix& other) const
  {
    Matrix A = *this;
    A += other;
    return A;
  }

  /** Element-wise subtraction of two matrices in-place. */
  Matrix& operator-=(const Matrix& other)
  {
    Assert(n_rows() == other.n_rows() &&
           n_cols() == other.n_cols(),
           "Dimension mismatch error.");
    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < n_cols(); ++j)
        m_data[i][j] -= other[i][j];
    return *this;
  }

  /** Element-wise subtraction of two matrices. */
  Matrix operator-(const Matrix& other) const
  {
    Matrix A = *this;
    A -= other;
    return A;
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
    Matrix C(n_rows(), other.n_cols());
    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < other.n_cols(); ++j)
      {
        number c_ij = 0.0;
        for (uint64_t k = 0; k < n_cols(); ++k)
          c_ij += m_data[i][k] * other[k][j];
        C[i][j] += c_ij;
      }
    return C;
  }

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector<number> operator*(const Vector<number>& x) const
  {
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    Vector<number> y(n_rows());
    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < n_cols(); ++j)
        y[i] += m_data[i][j] * x[j];
    return y;
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
    Matrix At(n_cols(), n_rows());
    for (uint64_t i = 0; i < At.n_rows(); ++i)
      for (uint64_t j = 0; j < At.n_cols(); ++j)
        At[i][j] = m_data[j][i];
    return At;
  }

  // @}

  //###########################################################################
  /** Print the matrix. */
  void print(std::ostream& os,
             const bool scientific = false,
             const unsigned int precision = 3,
             const unsigned int width = 0) const
  {
    unsigned int w                   = width;
    std::ios::fmtflags old_flags     = os.flags();
    unsigned int       old_precision = os.precision(precision);

    if (scientific)
    {
      os.setf(std::ios::scientific, std::ios::floatfield);
      w = (!width)? precision + 7 : w;
    }
    else
    {
      os.setf(std::ios::fixed, std::ios::floatfield);
      w = (!width)? precision + 4 : w;
    }

    for (uint64_t i = 0; i < n_rows(); ++i)
    {
      for (uint64_t j = 0; j < n_cols(); ++j)
        os << std::setw(w) << m_data[i][j];
      os << std::endl;
    }
    os << std::endl;
  }

private:
  static bool valid_stl_matrix(const std::vector<std::vector<number>>& matix)
  {
    uint64_t m = matix.front().size();
    for (const auto& row : matix)
      if (row.size() != m) return false;
    return true;
  }
};

}

#endif //MATRIX_H

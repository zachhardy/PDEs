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
  using value_type          = number;
  using size_type           = uint64_t;
  using pointer             = Vector<value_type>*;
  using const_pointer       = const Vector<value_type>*;
  using reference           = Vector<value_type>&;
  using const_reference     = const Vector<value_type>&;
  using iterator            = Vector<value_type>*;
  using const_iterator      = const Vector<value_type>*;

private:
  std::vector<Vector<value_type>> values;

public:
  /// Default constructor.
  Matrix()
    : values(0) {}

  /// Copy constructor.
  Matrix(const Matrix& other)
    : values(other.values) {}

  /// Move constructor.
  Matrix(Matrix&& other)
    : values(std::move(other.values)) {}

  /// Construct a Matrix with \p n_rows and \p n_cols.
  explicit
  Matrix(const size_type n_rows, const size_type n_cols)
    : values(n_rows, Vector<number>(n_cols)) {}

  /// Construct a Matrix with \p n_rows and \p n_cols set to \p value.
  explicit
  Matrix(const size_t n_rows, const size_type n_cols, const value_type value)
    : values(n_rows, Vector<number>(n_cols, value)) {}

  /// Construct a square Matrix with dimension \p n.
  explicit
  Matrix(const size_type n)
    : Matrix<number>(n, n) {}

  /// Construct a square Matrix with dimension \p n set to \p values.
  explicit
  Matrix(const size_type n, const value_type value)
    : Matrix<number>(n, n, value) {}

  /// Copy construction from an STL vector.
  Matrix(const std::vector<std::vector<value_type>>& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");
    values = other;
  }

  /// Move construction from an STL vector.
  Matrix(std::vector<std::vector<value_type>>&& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");
    values = other;
  }

  /// Construction from a nested initializer list.
  Matrix(std::initializer_list<std::initializer_list<value_type>> list)
  {
    std::vector<std::vector<value_type>> tmp;
    for (auto& row : list)
      tmp.push_back(row);

    Assert(valid_dimensions(tmp), "All rows must be the same length.");
    values = tmp;
  }

  /// Copy assignment.
  Matrix&
  operator=(const Matrix& other)
  {
    values = other.values;
    return *this;
  }

  /// Move assignment.
  Matrix&
  operator=(Matrix&& other)
  {
    values = std::move(other.values);
    return *this;
  }

  /// Copy assignment from an STL vector.
  Matrix&
  operator=(const std::vector<std::vector<value_type>>& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");

    values.clear();
    for (const auto& row : other)
      values.push_back(row);
  }

  /// Move assignment from an STL vector.
  Matrix&
  operator=(std::vector<std::vector<value_type>>&& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");

    values.clear();
    for (auto& row : other)
      values.push_back(std::move(row));
  }

  /// Equality comparison operator.
  bool
  operator==(const Matrix& other)
  { return values == other.values; }

  /// Inequality comparison operator.
  bool
  operator!=(const Matrix& other)
  { return values != other.values; }

  /** \name Information */
  // @{

  /// Return the number of rows in the Matrix.
  size_type
  n_rows() const
  { return values.size(); }

  /// Return the number of columns in the Matrix.
  size_type
  n_cols() const
  { return values.front().size(); }

  /// Return the number of elements in the Matrix.
  size_type
  size() const
  { return n_rows() * n_cols(); }

  /// Return the number of nonzero elements in the Matrix.
  size_type
  nnz() const
  {
    size_type count = 0;
    for (const auto& row : values)
      for (const auto& elem : row)
        if (elem != 0.0) ++count;
    return count;
  }

  /// Return whether the Matrix is empty.
  bool
  empty() const
  { return values.empty(); }

  /// Return whether the Matrix is uniformly zero.
  bool
  all_zero() const
  {return (nnz() == 0)? true : false; }

  // @}
  /** \name Iterators */
  // @{

  /// Mutable iterator to the start of the Matrix.
  iterator
  begin()
  { return values.begin(); }

  /// Constant iterator to the start of the Matrix.
  const_iterator
  begin() const
  { return values.begin(); }

  /// Mutable iterator to the end of the Matrix.
  iterator
  end()
  { return values.end(); }

  /// Constant iterator to the end of the Matrix.
  const_iterator
  end() const
  { return values.end(); }

  /// Mutable iterator to the start of row \p i.
  typename Vector<value_type>::iterator
  begin(const size_type i)
  { return values[i].begin(); }

  /// Constant iterator to the end of row \p i.
  typename Vector<value_type>::const_iterator
  begin(const size_type i) const
  { return values[i].begin(); }

  /// Mutable iterator to the end of row \p i.
  typename Vector<value_type>::iterator
  end(const size_type i)
  { return values[i].end(); }

  /// Constant iterator to the end of row \p i.
  typename Vector<value_type>::const_iterator
  end(const size_type i) const
  { return values[i].end(); }

  // @}
  /** \name Accessors */
  // @{

  /// Read/write access for row \p i.
  reference
  operator[](const size_type i)
  { return values[i]; }

  /// Read only access for row \p i.
  const_reference
  operator[](const size_type i) const
  { return values[i]; }

  /// Read/write access for row \p i. */
  reference
  operator()(const size_type i)
  { return values[i]; }

  /// Read access for row \p i. */
  const_reference
  operator()(const size_type i) const
  { return values[i]; }

  /// Read/write access for row \p i with bounds checking.
  reference
  at(const size_type i)
  { return values.at(i); }

  /// Read only access for row \p i with bounds checking.
  const_reference
  at(const size_type i) const
  { return values.at(i); }

  /// Read/write access for element <tt>(i, j)</tt>.
  typename Vector<value_type>::reference
  operator()(const size_type i, const size_type j)
  { return values[i][j]; }

  /// Read only access for element <tt>(i, j)</tt>.
  typename Vector<value_type>::const_reference
  operator()(const size_type i, const size_type j) const
  { return values[i][j]; }

  /// Read/write access for element <tt>(i, j)</tt> with bounds checking.
  typename Vector<value_type>::reference
  at(const size_type i, const size_type j)
  { return values.at(i).at(j); }

  /// Read only access for element <tt>(i, j)</tt> with bounds checking.
  typename Vector<value_type>::const_reference
  at(const size_type i, const size_type j) const
  { return values.at(i).at(j); }

  /// Read/write access to the <tt>i</tt>'th diagonal element.
  typename Vector<value_type>::reference
  diagonal(const size_type i)
  { return values.at(i).at(i); }

  /// Read access to the <tt>i</tt>'th diagonal element.
  typename Vector<value_type>::const_reference
  diagonal(const size_type i) const
  { return values.at(i).at(i); }

  /// Return the diagonal of the Matrix.
  Vector<value_type>
  diagonal() const
  {
    Vector<number> diag;
    uint64_t min_dim = std::min(n_rows(), n_cols());
    for (uint64_t i = 0; i < min_dim; ++i)
      diag.push_back(values[i][i]);
    return diag;
  }

  /// Return the underlying Matrix data.
  pointer
  data()
  { return values.data(); }

  /// Return the underlying Matrix data.
  const_pointer
  data() const
  { return values.data(); }

  /// Return the underlying data for row \p i.
  typename Vector<value_type>::pointer
  data(const size_type i)
  { return values[i].data(); }

  /// Return the underlying data for row \p i.
  typename Vector<value_type>::const_pointer
  data(const size_type i) const
  { return values[i].data(); }

  // @}
  /** \name Modifiers */
  // @{

  /// Return the Matrix to an uninitialized state.
  void
  clear()
  { values.clear(); }

  /// Remove the last row of the Matrix.
  void
  pop_back()
  { values.pop_back(); }

  /// Add a row to the back of the Matrix.
  void
  push_back(const Vector<value_type>& row)
  {
    Assert(row.size() == n_cols(), "Dimension mismatch error.");
    values.push_back(row);
  }

  /// Move a row to the back of the Matrix.
  void
  push_back(Vector<number>&& row)
  {
    Assert(row.size() == n_cols(), "Dimension mismatch error.")
    values.emplace_back(row);
  }

  /// Resize to \p n_rows and \p n_cols.
  void
  resize(const size_type n_rows, const size_type n_cols)
  { values.resize(n_rows, Vector<value_type>(n_cols)); }

  /// Resize to dimension \p n.
  void
  resize(const size_type n)
  { resize(n, n); }

  /// Resize to \p n_rows and \p n_cols, setting new elements to \p value. */
  void
  resize(const size_type n_rows,
         const size_type n_cols,
         const value_type value)
  { values.resize(n_rows, Vector<value_type>(n_cols, value)); }

  /// Resize to dimension \p n, setting new elements to \p value. */
  void
  resize(const size_type n, const value_type value)
  { resize(n, n, value); }

  /// Swap the elements of two rows.
  void
  swap_row(const size_type i0, const size_type i1)
  {
    Assert(i0 < n_rows() && i1 < n_rows(), "Out of range error.");
    values[i0].swap(values[i1]);
  }

  /// Swap the elements of two columns.
  void
  swap_column(const size_type j0, const size_type j1)
  {
    Assert(j0 < n_cols() && j1 < n_cols(), "Out of range error.");
    for (uint64_t i = 0; i < n_rows(); ++i)
      std::swap(values[i][j0], values[i][j1]);
  }

  /// Swap the elements of this matrix with another.
  void
  swap(Matrix& other)
  { values.swap(other.values); }

  /// Set the diagonal of a matrix.
  void
  set_diagonal(const Vector<number>& diag)
  {
    if (values.empty())
    {
      resize(diag.size());
      for (uint64_t i = 0; i < diag.size(); ++i)
        values[i][i] = diag[i];
    }
    else
    {
      uint64_t min_dim = std::min(n_rows(), n_cols());
      Assert(diag.size() == min_dim, "Dimension mismatch error.");
      for (uint64_t i = 0; i < min_dim; ++i)
        values[i][i] = diag[i];
    }
  }

  /// Set the diagonal with a fixed scalar value.
  void
  set_diagonal(const value_type value)
  {
    if (values.empty()) {
      resize(1, value);
      return;
    }
    else
    {
      uint64_t min_dim = std::min(n_rows(), n_cols());
      for (uint64_t i = 0; i < min_dim; ++i)
        values[i][i] = value;
    }
  }

  // @}
  /** \name Scalar Operations */
  // @{

  /// Element-wise negation in-place.
  Matrix&
  operator-()
  {
    for (auto& row : values)
      for (auto& entry : row)
        entry = -entry;
    return *this;
  }

  /// Element-wise negation.
  Matrix
  operator-() const
  {
    Matrix A = *this;
    return -A;
  }

  /// Element-wise multiplication by a scalar in-place.
  Matrix&
  operator*=(const value_type value)
  {
    for (auto& row : *this)
      for (auto& entry : row)
        entry *= value;
    return *this;
  }

  /// Element-wise multiplication by a scalar.
  Matrix
  operator*(const value_type value) const
  {
    Matrix A = *this;
    A *= value;
    return A;
  }

  /// Element-wise division by a scalar in-place.
  Matrix&
  operator/=(const value_type value)
  {
    Assert(value != 0.0, "Division by zero error.");

    for (auto& row : *this)
      for (auto& entry : row)
        entry /= value;
    return *this;
  }

  /// Element-wise division by a scalar.
  Matrix
  operator/(const number value) const
  {
    Matrix A = *this;
    A /= value;
    return A;
  }

  //@}
  /** \name Linear Algebra */
  // @{

  /// Element-wise addition of two matrices in-place.
  Matrix&
  operator+=(const Matrix& other)
  {
    Assert(n_rows() == other.n_rows() &&
           n_cols() == other.n_cols(),
           "Dimension mismatch error.");

    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < n_cols(); ++j)
        values[i][j] += other[i][j];
    return *this;
  }

  /// Element-wise addition of two matrices.
  Matrix
  operator+(const Matrix& other) const
  {
    Matrix A = *this;
    A += other;
    return A;
  }

  /// Element-wise subtraction of two matrices in-place.
  Matrix&
  operator-=(const Matrix& other)
  {
    Assert(n_rows() == other.n_rows() &&
           n_cols() == other.n_cols(),
           "Dimension mismatch error.");

    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < n_cols(); ++j)
        values[i][j] -= other[i][j];
    return *this;
  }

  /// Element-wise subtraction of two matrices.
  Matrix
  operator-(const Matrix& other) const
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
  Matrix
  operator*(const Matrix& other) const
  {
    Assert(n_cols() == other.n_rows(), "Dimension mismatch error.");

    Matrix C(n_rows(), other.n_cols());
    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < other.n_cols(); ++j)
      {
        number c_ij = 0.0;
        for (uint64_t k = 0; k < n_cols(); ++k)
          c_ij += values[i][k] * other[k][j];
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
  Vector<value_type>
  operator*(const Vector<value_type>& x) const
  {
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    Vector<number> y(n_rows());
    for (uint64_t i = 0; i < n_rows(); ++i)
      for (uint64_t j = 0; j < n_cols(); ++j)
        y[i] += values[i][j] * x[j];
    return y;
  }

  /**
   * Return the transpose of a matrix.
   * This is computed via
   * \f[ \boldsymbol{B} = \boldsymbol{A}^T \\
   *     b_{ij} = a_{ji}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  Matrix
  transpose() const
  {
    Matrix At(n_cols(), n_rows());
    for (uint64_t i = 0; i < At.n_rows(); ++i)
      for (uint64_t j = 0; j < At.n_cols(); ++j)
        At[i][j] = values[j][i];
    return At;
  }

  // @}

  //###########################################################################
  /** Print the matrix. */
  void
  print(std::ostream& os,
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
        os << std::setw(w) << values[i][j];
      os << std::endl;
    }
    os << std::endl;
  }

private:
  static bool
  valid_dimensions(const std::vector<std::vector<value_type>>& matix)
  {
    uint64_t m = matix.front().size();
    for (const auto& row : matix)
      if (row.size() != m) return false;
    return true;
  }
};

}

#endif //MATRIX_H

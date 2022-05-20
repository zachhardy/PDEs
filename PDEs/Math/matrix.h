#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"
#include "macros.h"

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
  using value_type = number;
  using size_type = uint64_t;
  using pointer = Vector<value_type>*;
  using const_pointer = const Vector<value_type>*;
  using reference = Vector<value_type>&;
  using const_reference = const Vector <value_type>&;
  using iterator = Vector <value_type>*;
  using const_iterator = const Vector<value_type>*;

protected:
  std::vector<Vector<value_type>> values;

public:
  /// Default constructor.
  Matrix()
      : values(0)
  {}

  /// Copy constructor.
  Matrix(const Matrix& other)
      : values(other.values)
  {}

  /// Move constructor.
  Matrix(Matrix&& other)
      : values(std::move(other.values))
  {}

  /// Construct a Matrix with \p n_rows and \p n_cols.
  explicit
  Matrix(const size_type n_rows, const size_type n_cols)
      : values(n_rows, Vector<number>(n_cols))
  {}

  /// Construct a Matrix with \p n_rows and \p n_cols set to \p value.
  explicit
  Matrix(const size_t n_rows, const size_type n_cols, const value_type value)
      : values(n_rows, Vector<number>(n_cols, value))
  {}

  /// Copy construction from an STL vector.
  Matrix(const std::vector<std::vector<value_type>>& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");

    values.clear();
    for (const auto& row : other)
      values.push_back(row);
  }

  /// Move construction from an STL vector.
  Matrix(std::vector<std::vector<value_type>>&& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");

    values.clear();
    for (const auto& row : other)
      values.push_back(std::move(row));
  }

  /// Construction from a nested initializer list.
  Matrix(std::initializer_list<std::initializer_list<value_type>> list)
  {
    std::vector<std::vector<value_type>> tmp;
    for (auto& row : list)
      tmp.push_back(row);
    Assert(valid_dimensions(tmp), "All rows must be the same length.");

    values.clear();
    for (const auto& row : tmp)
      values.push_back(row);
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

  /// Assignment to a scalar value.
  Matrix&
  operator=(const value_type value)
  {
    Assert(!empty(), "Empty matrix error.");
    for (size_type i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = values[i].data();
      for (size_type j = 0; j < n_cols(); ++j)
        *a_ij++ = value;
    }
    return *this;
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
  { return (nnz() == 0) ? true : false; }

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
  {
    Assert(i < n_rows(), "Out of range error.");
    return values[i].begin();
  }

  /// Constant iterator to the end of row \p i.
  typename Vector<value_type>::const_iterator
  begin(const size_type i) const
  {
    Assert(i < n_rows(), "Out of range error.")
    return values[i].begin();
  }

  /// Mutable iterator to the end of row \p i.
  typename Vector<value_type>::iterator
  end(const size_type i)
  {
    Assert(i < n_rows(), "Out of range error.")
    return values[i].end();
  }

  /// Constant iterator to the end of row \p i.
  typename Vector<value_type>::const_iterator
  end(const size_type i) const
  {
    Assert(i < n_rows(), "Out of range error.")
    return values[i].end();
  }

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
  Vector <value_type>
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
  push_back(const Vector <value_type>& row)
  {
    Assert(row.size() == n_cols(), "Dimension mismatch error.");
    values.push_back(row);
  }

  /// Move a row to the back of the Matrix.
  void
  push_back(Vector <number>&& row)
  {
    Assert(row.size() == n_cols(), "Dimension mismatch error.")
    values.emplace_back(row);
  }

  /// Resize to \p n_rows and \p n_cols.
  void
  resize(const size_type n_rows, const size_type n_cols)
  { values.resize(n_rows, Vector<value_type>(n_cols)); }

  /// Resize to \p n_rows and \p n_cols, setting new elements to \p value. */
  void
  resize(const size_type n_rows,
         const size_type n_cols,
         const value_type value)
  { values.resize(n_rows, Vector<value_type>(n_cols, value)); }

  /// See \ref resize
  void
  reinit(const size_type n_rows, const size_type n_cols)
  { resize(n_rows, n_cols); }

  /// See \ref resize
  void
  reinit(const size_type n_rows,
         const size_type n_cols,
         const value_type value)
  { resize(n_rows, n_cols, value); }

  /// Swap the elements of two rows.
  void
  swap_row(const size_type i, const size_type k)
  {
    Assert(i < n_rows() && k < n_rows(), "Out of range error.");
    values[i].swap(values[k]);
  }

  /// Swap the elements of two columns.
  void
  swap_column(const size_type j, const size_type k)
  {
    Assert(j < n_cols() && k < n_cols(), "Out of range error.");
    for (uint64_t i = 0; i < n_rows(); ++i)
      std::swap(values[i][j], values[i][k]);
  }

  /// Swap the elements of this matrix with another.
  void
  swap(Matrix& other)
  { values.swap(other.values); }

  /// Set the diagonal of a matrix.
  void
  set_diagonal(const Vector <number>& diag)
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
    Assert(!empty(), "Empty matrix error.");

    uint64_t min_dim = std::min(n_rows(), n_cols());
    for (uint64_t i = 0; i < min_dim; ++i)
      values[i][i] = value;
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

  /**
   * Addition of a scaled Matrix.
   *
   * \f[ \boldsymbol{A} += \alpha \boldsymbol{B} \\
   *     a_{ij} +=  a_{ij} + \alpha b_{ij}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  void add(const Matrix& B, const value_type factor = 1.0)
  {
    Assert(!empty(), "Empty matrix error.")
    Assert(n_rows() == B.n_rows(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_cols(), "Dimension mismatch error.");

    for (size_type i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = values[i].data();
      const value_type* b_ij = B.data(i);
      for (size_type j = 0; j < n_cols(); ++j)
        *a_ij++ += factor * *b_ij++;
    }
  }

  /**
   * Addition of a scaled transpose Matrix.
   *
   * \f[ \boldsymbol{A} += \alpha \boldsymbol{B}^T \\
   *     a_{ij} += a_{ij} + \alpha b_{ji}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  void Tadd(const Matrix& B, const value_type factor = 1.0)
  {
    Assert(n_rows() == B.n_cols(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

    for (size_type i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = values[i].data();
      for (size_type j = 0; j < n_cols(); ++j)
        *a_ij++ += factor * B[j][i];
    }
  }

  /// Element-wise addition of two matrices in-place.
  Matrix&
  operator+=(const Matrix& B)
  {
    Assert(n_rows() == B.n_rows() &&
           n_cols() == B.n_cols(),
           "Dimension mismatch error.");

    for (uint64_t i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = values[i].data();
      const value_type* b_ij = B.data(i);
      for (uint64_t j = 0; j < n_cols(); ++j)
        *a_ij++ += *b_ij++;
    }
    return *this;
  }

  /// Element-wise addition of two matrices.
  Matrix
  operator+(const Matrix& B) const
  {
    Matrix A = *this;
    A += B;
    return A;
  }

  /// Element-wise subtraction of two matrices in-place.
  Matrix&
  operator-=(const Matrix& B)
  {
    Assert(n_rows() == B.n_rows() &&
           n_cols() == B.n_cols(),
           "Dimension mismatch error.");

    for (uint64_t i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = values[i].data();
      const value_type* b_ij = B.data(i);
      for (uint64_t j = 0; j < n_cols(); ++j)
        *a_ij++ -= *b_ij++;
    }
    return *this;
  }

  /// Element-wise subtraction of two matrices.
  Matrix
  operator-(const Matrix& B) const
  {
    Matrix A = *this;
    A -= B;
    return A;
  }

  /**
   * Compute a matrix-matrix multiplication.
   *
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B} \\
   *     c_{ij} = \sum_{k=0}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  void
  mmult(const Matrix& B, Matrix& C,
        const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_rows(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

    for (size_type i = 0; i < C.n_rows(); ++i)
    {
      const value_type* a_i = values[i].data();
      for (size_type j = 0; j < C.n_cols(); ++j)
      {
        value_type c_ij = adding ? C(i, j) : 0.0;
        for (size_type k = 0; k < n_cols(); ++k)
          c_ij += a_i[k] * B(k, j);
        C(i, j) = c_ij;
      }
    }
  }

  /// See \ref mmult.
  Matrix
  mmult(const Matrix& B)
  {
    Matrix C(n_rows(), B.n_cols());
    mmult(B, C);
    return C;
  }

  /**
   * Compute a transpose matrix-matrix multiplication.
   *
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B} \\
   *     c_{ij} = \sum_{k=0}^{n} a_{ki} b_{kj}, \hspace{0.25cm}, \forall i, j
   * \f]
   */
  void
  Tmmult(const Matrix& B, Matrix& C,
         const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_cols(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.");
    Assert(n_rows() == B.n_rows(), "Dimension mismatch error.");

    for (size_type i = 0; i < C.n_rows(); ++i)
      for (size_type j = 0; j < C.n_cols(); ++j)
      {
        value_type c_ij = adding ? C(i, j) : 0.0;
        for (size_type k = 0; k < n_rows(); ++k)
          c_ij += values[k][i] * B(k, j);
        C(i, j) = c_ij;
      }
  }

  /// See \ref Tmmult.
  Matrix
  Tmmult(const Matrix& B)
  {
    Matrix C(n_cols(), B.n_cols());
    Tmmult(B, C);
    return C;
  }

  /**
   * Compute a matrix-transpose matrix multiplication.
   *
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}^T \\
   *     c_{ij} = \sum_{k=1}^{n} a_{ik} b_{jk}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  void
  mTmult(const Matrix& B, Matrix& C,
         const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_rows(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_cols(), "Dimension mismatch error.");

    for (size_type i = 0; i < C.n_rows(); ++i)
    {
      const value_type* a_i = values[i].data();
      for (size_type j = 0; j < C.n_cols(); ++j)
      {
        const value_type* b_j = B.data(j);
        value_type c_ij = adding ? C(i, j) : 0.0;
        for (size_type k = 0; k < n_cols(); ++k)
          c_ij += a_i[k] * *b_j++;
        C(i, j) = c_ij;
      }
    }
  }

  /// See \ref mTmult.
  Matrix
  mTmult(const Matrix& B)
  {
    Matrix C(n_rows(), B.n_rows());
    mTmult(B, C);
    return C;
  }

  /**
   * Compute a transpose matrix-transpose matrix multiplication.
   *
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T \\
   *     c_{ij} = \sum_{k=1}^{n} a_{ki} b_{jk}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  void
  TTmult(const Matrix& B, Matrix& C,
         const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_cols(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.");
    Assert(n_rows() == B.n_cols(), "Dimension mismatch error.");

    for (size_type i = 0; i < C.n_rows(); ++i)
      for (size_type j = 0; j < C.n_cols(); ++j)
      {
        const value_type* b_j = B.data(j);
        value_type c_ij = adding ? C(i, j) : 0.0;
        for (size_type k = 0; k < n_rows(); ++k)
          c_ij += values[k][i] * *b_j++;
        C(i, j) = c_ij;
      }
  }

  /// See \ref TTmult.
  Matrix
  TTmult(const Matrix& B)
  {
    Matrix C(n_cols(), B.n_rows());
    TTmult(B, C);
    return C;
  }

  /**
   * Return the matrix-matrix product.
   * \see mmult
   */
  Matrix
  operator*(const Matrix& B) const
  {
    Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

    Matrix C(n_rows(), B.n_cols());
    for (uint64_t i = 0; i < C.n_rows(); ++i)
    {
      const value_type* a_i = values[i].data();
      for (uint64_t j = 0; j < C.n_cols(); ++j)
      {
        value_type c_ij = 0.0;
        for (uint64_t k = 0; k < n_cols(); ++k)
          c_ij += a_i[k] * B(k, j);
        C(i, j) = c_ij;
      }
    }
    return C;
  }

  /**
   * Compute a matrix-vector product.
   *
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  void
  vmult(const Vector<value_type>& x,
        Vector<value_type>& y,
        const bool adding = false)
  {
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    Assert(y.size() == n_rows(), "Dimension mismatch error.");

    for (size_type i = 0; i < n_rows(); ++i)
    {
      value_type v = adding ? y[i] : 0.0;
      const value_type* a_ij = values[i].data();
      for (size_type j = 0; j < n_cols(); ++j)
        v += *a_ij++ * x[j];
      y[i] = v;
    }
  }

  Vector<value_type>
  vmult(const Vector<value_type>& x)
  {
    Vector<value_type> y(n_rows());
    vmult(x, y);
    return y;
  }

  /**
   * Add a matrix-vector product to the destination vector.
   * \see vmult
   */
  void
  vmult_add(const Vector<value_type>& x,
            Vector<value_type>& dst)
  { vmult(x, dst, true); }

  /**
   * Compute a transpose matrix-vector product.
   *
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A}^T \vec{x} \\
   *     y_i = \sum_{i=1}^{n} a_{ji} x_i, \hspace{0.25cm} \forall i
   * \f]
   */
  void
  Tvmult(const Vector<value_type>& x,
         Vector<value_type>& y,
         const bool adding = false)
  {
    Assert(x.size() == n_rows(), "Dimension mismatch error.");
    Assert(y.size() == n_cols(), "Dimension mismatch error.");

    y = adding ? y : 0.0;
    for (size_type i = 0; i < n_rows(); ++i)
    {
      const value_type x_i = x[i];
      const value_type* a_ij = values[i].data();
      for (size_type j = 0; j < n_cols(); ++j)
        y[j] += *a_ij++ * x_i;
    }
  }

  /// See \ref Tvmult.
  Vector<value_type>
  Tvmult(const Vector<value_type>& x)
  {
    Vector<value_type> y(n_cols());
    Tvmult(x, y);
    return y; }

  /**
   * Add a transpose matrix-vector product to the destination vector.
   * \see Tvmult
   */
  void
  Tvmult_add(const Vector<value_type>& x,
         Vector<value_type>& dst)
  { Tvmult(x, dst, true); }

  /**
   * Compute a matrix-vector product.
   * \see vmult
   */
  Vector<value_type>
  operator*(const Vector<value_type>& x) const
  {
    Assert(x.size() == n_cols(), "Dimension mismatch error.");

    Vector<number> y(n_rows());
    for (uint64_t i = 0; i < n_rows(); ++i)
    {
      value_type y_i = 0.0;
      const value_type* a_ij = values[i].data();
      for (uint64_t j = 0; j < n_cols(); ++j)
        y_i += *a_ij++ * x[j];
      y[i] = y_i;
    }
    return y;
  }

  /// Return the transpose of the Matrix.
  Matrix
  transpose() const
  {
    Matrix A_T(n_cols(), n_rows());
    for (size_type i = 0; i < n_rows(); ++i)
    {
      const value_type* a_ij = values[i].data();
      for (size_type j = 0; j < n_cols(); ++j)
        A_T[j][i] = *a_ij++;
    }
    return A_T;
  }

  // @}
  /** \name Printing Utilities */
  // @{

  /// Print the matrix.
  void
  print(std::ostream& os = std::cout,
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
      w = (!width)? precision + 10 : w;
    }
    else
    {
      os.setf(std::ios::fixed, std::ios::floatfield);
      w = (!width)? precision + 5 : w;
    }

    for (uint64_t i = 0; i < n_rows(); ++i)
    {
      const value_type* a_ij = values[i].data();
      for (uint64_t j = 0; j < n_cols(); ++j)
        os << std::setw(w) << *a_ij++;
      os << std::endl;
    }
    os << std::endl;
  }

  // @}

private:
  static bool
  valid_dimensions(const std::vector<std::vector<value_type>>& A)
  {
    uint64_t m = A.front().size();
    for (const auto& row : A)
      if (row.size() != m) return false;
    return true;
  }
};

}

#endif //MATRIX_H

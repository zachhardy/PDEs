#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"
#include "macros.h"

#include <cmath>
#include <vector>
#include <iomanip>
#include <cinttypes>


namespace pdes::Math
{

class Matrix
{
public:
  using value_type = double;

protected:
  std::vector<Vector> coeffs;

public:
  /**
   * Default constructor.
   */
  Matrix() = default;

  /**
   * Construct a matrix with \p n_rows and \p n_cols.
   */
  explicit
  Matrix(const size_t n_rows, const size_t n_cols) :
    coeffs(n_rows, Vector(n_cols))
  {}

  /**
   * Construct a Matrix with \p n_rows and \p n_cols set to \p value.
   */
  explicit
  Matrix(const size_t n_rows,
         const size_t n_cols,
         const value_type value) :
    coeffs(n_rows, Vector(n_cols, value))
  {}

  /**
   * Copy construction with underlying data
   */
  Matrix(const std::vector<Vector>& other) : coeffs(other)
  {
    Assert(valid_dimensions(coeffs), "All rows must be the same length.")
  }

  /**
   * Copy construction with nested STL vectors.
   */
  Matrix(const std::vector<std::vector<value_type>>& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");
    for (const auto& row : other)
      coeffs.push_back(row);
  }

  /**
   * Move construction from nested STL vectors.
   */
  Matrix(std::vector<std::vector<value_type>>&& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");
    for (auto& row : other)
      coeffs.push_back(std::move(row));
  }

  /**
   * Construction from nested initializer lists.
   */
  Matrix(const std::initializer_list<std::initializer_list<value_type>> list)
  {
    std::vector<std::vector<value_type>> tmp;
    for (const auto& row : list)
      tmp.push_back(row);
    Assert(valid_dimensions(tmp), "All rows must be the same length.");

    for (const auto& row : tmp)
      coeffs.push_back(row);
  }

  /**
   * Copy assignment from nested STL vectors.
   */
  Matrix&
  operator=(const std::vector<std::vector<value_type>>& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");
    for (const auto& row : other)
      coeffs.push_back(row);
    return *this;
  }

  /**
   * Move assignment from nested STL vectors.
   */
  Matrix&
  operator=(std::vector<std::vector<value_type>>&& other)
  {
    Assert(valid_dimensions(other), "All rows must be the same length.");
    for (auto row : other)
      coeffs.push_back(std::move(row));
    return *this;
  }

  /**
   * Assignment to a scalar value.
   */
  Matrix&
  operator=(const value_type value)
  {
    Assert(!empty(), "Empty matrix error.");

    for (auto& row : coeffs)
      for (auto& el : row)
        el = value;
    return *this;
  }

  /**
   * Test the equality of two matrices.
   */
  bool
  operator==(const Matrix& other)
  { return (coeffs == other.coeffs); }

  /**
   * Test the inequality of two matrics.
   */
  bool
  operator!=(const Matrix& other)
  { return (coeffs != other.coeffs); }

  /** \name Information */
  // @{

  /**
   * Return the number of rows in the matrix.
   */
  size_t
  n_rows() const
  { return coeffs.size(); }

  /**
   * Return the number of columns in the matrix.
   */
  size_t 
  n_cols() const
  { return coeffs.front().size(); }

  /**
   * Return the number of elements in the matrix.
   */
  size_t
  size() const
  { return n_rows() * n_cols(); }

  /**
   * Return the number of non-zero elements in the matrix.
   */
  size_t
  nnz() const
  {
    size_t count = 0;
    for (const auto& row : coeffs)
      for (const auto& elem : row)
        if (elem != 0.0) ++count;
    return count;
  }

  /**
   * Return whether the matrix is empty.
   */
  bool
  empty() const
  { return coeffs.empty(); }

  /**
   * Return whether the matrix is uniformly zero.
   */
  bool
  all_zero() const
  { return (nnz() == 0) ? true : false; }

  // @}
  /** \name Iterators */
  // @{

  /**
   * Mutable iterator to the start of the matrix.
   */
  std::vector<Vector>::iterator
  begin()
  { return coeffs.begin(); }

  /**
   * Mutable iterator to the end of the matrix.
   */
   std::vector<Vector>::iterator
   end()
  { return coeffs.end(); }

  /**
   * Constant iterator to the start of the matrix.
   */
  std::vector<Vector>::const_iterator
  begin() const
  { return coeffs.begin(); }

  /**
   * Constant iterator to the end of the matrix.
   */
  std::vector<Vector>::const_iterator
  end() const
  { return coeffs.end(); }

  /**
   * Mutable iterator to the start of row \p i of the matrix.
   */
  std::vector<double>::iterator
  begin(const size_t i)
  {
    Assert(i < n_rows(), "Out of range error.");
    return coeffs[i].begin();
  }

  /**
   * Mutable iterator to the end of row \p i of the matrix.
   */
  std::vector<double>::iterator
  end(const size_t i)
  {
    Assert(i < n_rows(), "Out of range error.");
    return coeffs[i].end();
  }

  /**
   * Constant iterator to the start of row \p i of the matrix.
   */
  std::vector<double>::const_iterator
  begin(const size_t i) const
  {
    Assert(i < n_rows(), "Out of range error.");
    return coeffs[i].begin();
  }

  /**
   * Constant iterator to the end of row \p i of the matrix.
   */
  std::vector<double>::const_iterator
  end(const size_t i) const
  {
    Assert(i < n_rows(), "Out of range error.");
    return coeffs[i].end();
  }

  // @}
  /** \name Accessors */
  // @{

  /**
   * Read and write access for row \p i of the matrix.
   */
  Vector&
  operator[](const size_t i)
  { return coeffs[i]; }

  /**
   * Read access for row \p i of the matrix.
   */
  const Vector&
  operator[](const size_t i) const
  { return coeffs[i]; }

  /**
   * Read and write access for row \p i of the matrix.
   */
  Vector&
  operator()(const size_t i)
  { return coeffs[i]; }

  /**
   * Read access for row \p i of the matrix.
   */
  const Vector&
  operator()(const size_t i) const
  { return coeffs[i]; }

  /**
   * Read and write access for row \p i of the matrix with bounds checking.
   */
  Vector&
  at(const size_t i)
  { return coeffs.at(i); }

  /**
   * Read access for row \p i of the matrix with bounds checking.
   */
  const Vector&
  at(const size_t i) const
  { return coeffs.at(i); }

  /**
   * Read and write access for element <tt>(i, j)</tt>.
   */
  value_type&
  operator()(const size_t i, const size_t j)
  { return coeffs[i][j]; }

  /**
   * Read and write access for element <tt>(i, j)</tt>.
   */
  const value_type&
  operator()(const size_t i, const size_t j) const
  { return coeffs[i][j]; }

  /**
   * Read and write access for element <tt>(i, j)</tt> with bounds checking.
   */
  value_type&
  at(const size_t i, const size_t j)
  { return coeffs.at(i).at(j); }

  /**
   * Read and write access for element <tt>(i, j)</tt> with bounds checking.
   */
  const value_type&
  at(const size_t i, const size_t j) const
  { return coeffs.at(i).at(j); }

  /**
   * Read and write access to the <tt>i</tt>'th diagonal element.
   */
   value_type&
   diagonal(const size_t i)
  { return coeffs.at(i).at(i); }

  /**
   * Read access to the <tt>i</tt>'th diagonal element.
   */
  const value_type&
  diagonal(const size_t i) const
  { return coeffs.at(i).at(i); }

  /**
   * Return the diagonal of the matrix.
   */
  Vector
  diagonal() const
  {
    Vector diag;
    size_t min_dim = std::min(n_rows(), n_cols());
    for (size_t i = 0; i < min_dim; ++i)
      diag.push_back(coeffs[i][i]);
    return diag;
  }

  /**
   * Return a mutable pointer to the underlying matrix data.
   */
  Vector*
  data()
  { return coeffs.data(); }

  /**
   * Return a constant pointer to the underlying matrix data.
   */
  const Vector*
  data() const
  { return coeffs.data(); }

  /**
   * Return a pointer to the underlying data for row \p i of the matrix.
   */
  value_type*
  data(const size_t i)
  {
    Assert(i < n_rows(), "Out of range error.");
    return coeffs[i].data();
  }

  /**
   * Return a pointer to the underlying data for row \p i of the matrix.
   */
  const value_type*
  data(const size_t i) const
  {
    Assert(i < n_rows(), "Out of range error.");
    return coeffs[i].data();
  }

  // @}
  /** \name Modifiers */
  // @{

  /**
   * Return the matrix to its uninitialized state.
   */
  void
  clear()
  { coeffs.clear(); }

  /**
   * Remove the last row from the matrix.
   */
  void
  pop_back()
  { coeffs.pop_back(); }

  /**
   * Add a row to the back of the matrix.
   */
  void
  push_back(const Vector& row)
  {
    Assert(row.size() == n_cols(), "Dimension mismatch error.");
    coeffs.push_back(row);
  }

  /**
   * Move a row to the back of the matrix.
   */
  void
  push_back(Vector&& row)
  {
    Assert(row.size() == n_cols(), "Dimension mismatch error.")
    coeffs.push_back(row);
  }

  /**
   * Resize the matrix to \p n_rows and \p n_cols.
   */
  void
  resize(const size_t n_rows, const size_t n_cols)
  { coeffs.resize(n_rows, Vector(n_cols)); }

  /**
   * Resize to \p n_rows and \p n_cols, setting new elements to \p value.
   */
  void
  resize(const size_t n_rows,
         const size_t n_cols,
         const value_type value)
  { coeffs.resize(n_rows, Vector(n_cols, value)); }

  /**
   * Alias to \ref resize.
   */
  void
  reinit(const size_t n_rows, const size_t n_cols)
  { resize(n_rows, n_cols); }

  /**
   * Alias to \ref resize.
   */
  void
  reinit(const size_t n_rows,
         const size_t n_cols,
         const value_type value)
  { resize(n_rows, n_cols, value); }

  /**
   * Swap two rows of the matrix.
   */
  void
  swap_row(const size_t i, const size_t k)
  {
    Assert(i < n_rows() && k < n_rows(), "Out of range error.");
    coeffs[i].swap(coeffs[k]);
  }

  /**
   * Swap two columns of the matrix.
   */
  void
  swap_column(const size_t j, const size_t k)
  {
    Assert(j < n_cols() && k < n_cols(), "Out of range error.");
    for (uint64_t i = 0; i < n_rows(); ++i)
      std::swap(coeffs[i][j], coeffs[i][k]);
  }

  /**
   * Swap the contents of this matrix with another.
   */
  void
  swap(Matrix& other)
  { coeffs.swap(other.coeffs); }

  /**
   * Set the diagonal of the matrix with a vector.
   */
  void
  set_diagonal(const Vector& diag)
  {
    if (coeffs.empty())
    {
      resize(diag.size(), diag.size());
      for (size_t i = 0; i < diag.size(); ++i)
        coeffs[i][i] = diag[i];
    }
    else
    {
      size_t min_dim = std::min(n_rows(), n_cols());
      Assert(diag.size() == min_dim, "Dimension mismatch error.");

      for (size_t i = 0; i < min_dim; ++i)
        coeffs[i][i] = diag[i];
    }
  }

  /**
   * Set the diagonal of the matrix with a fixed scalar value.
   */
  void
  set_diagonal(const value_type value)
  {
    if (empty())
      resize(1, 1, value);
    else
    {
      size_t min_dim = std::min(n_rows(), n_cols());
      for (size_t i = 0; i < min_dim; ++i)
        coeffs[i][i] = value;
    }
  }

  // @}
  /** \name Scalar Operations */
  // @{

  /**
   * Negate the elements of the matrix. This is compute via
   * \f$ \boldsymbol{A} = -\boldsymbol{A} = -a_{ij}, ~ \forall i, j \f$.
   */
  Matrix&
  operator-()
  {
    for (auto& row : coeffs)
      row = -row;
    return *this;
  }

  /**
   * Return a matrix with the negated elements.
   * \see Matrix::operator-()
   */
  Matrix
  operator-() const
  { return -Matrix(*this); }

  /**
   * Multiply the elements of the matrix by a scalar value.
   * This is computed via
   * \f$ \boldsymbol{A} = \alpha \boldsymbol{A}
   *                    = \alpha a_{ij} ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator*=(const value_type value)
  {
    for (auto& row : *this)
      row *= value;
    return *this;
  }

  /**
   * Divide the elements of the matrix by a scalar value.
   * This is computed via
   * \f$ \boldsymbol{A} = \frac{\boldsymbol{A}}{\alpha}
   *                    = \frac{a_{ij}}{\alpha}, ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator/=(const value_type value)
  {
    Assert(value != 0.0, "Division by zero error.");
    for (auto& row : *this)
      row /= value;
    return *this;
  }

  //@}
  /** \name Linear Algebra */
  // @{

  /**
   * Addition of a scaled matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} + \alpha \boldsymbol{B}
   *                    = a_{ij} + \alpha b_{ij}, ~ \forall i, j
   * \f$.
   */
  void add(const Matrix& B, const value_type factor = 1.0)
  {
    Assert(!empty(), "Empty matrix error.")
    Assert(n_rows() == B.n_rows(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_cols(), "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = data(i);
      const value_type* b_ij = B.data(i);
      for (size_t j = 0; j < n_cols(); ++j)
        *a_ij++ += factor * *b_ij++;
    }
  }

  /**
   * Addition of a scaled transpose Matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} + \alpha \boldsymbol{B}^T
   *                    = a_{ij} + \alpha b{ji}, ~ \forall i, j
   * \f$.
   */
  void Tadd(const Matrix& B, const value_type factor = 1.0)
  {
    Assert(n_rows() == B.n_cols(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = data(i);
      for (size_t j = 0; j < n_cols(); ++j)
        *a_ij++ += factor * B[j][i];
    }
  }

  /**
   * Add another matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} + \boldsymbol{B}
   *                    = a_{ij} + b_{ij}, ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator+=(const Matrix& B)
  {
    Assert(n_rows() == B.n_rows() &&
           n_cols() == B.n_cols(),
           "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = data(i);
      const value_type* b_ij = B.data(i);
      for (size_t j = 0; j < n_cols(); ++j)
        *a_ij++ += *b_ij++;
    }
    return *this;
  }

  /**
   * Subtract another matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} - \boldsymbol{B}
   *                    = a_{ij} - b_{ij}, ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator-=(const Matrix& B)
  {
    Assert(n_rows() == B.n_rows() &&
           n_cols() == B.n_cols(),
           "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
    {
      value_type* a_ij = data(i);
      const value_type* b_ij = B.data(i);
      for (size_t j = 0; j < n_cols(); ++j)
        *a_ij++ -= *b_ij++;
    }
    return *this;
  }

  /**
   * Compute a matrix-matrix multiplication. This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}
   *                    = \sum_{k=0}^{n} a_{ik} b_{kj}, ~ \forall i, j
   * \f$.
   */
  void
  mmult(const Matrix& B, Matrix& C,
        const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_rows(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

    for (size_t i = 0; i < C.n_rows(); ++i)
    {
      const value_type* a_i = data(i);
      for (size_t j = 0; j < C.n_cols(); ++j)
      {
        value_type c_ij = adding? C(i, j) : 0.0;
        for (size_t k = 0; k < n_cols(); ++k)
          c_ij += a_i[k] * B(k, j);
        C(i, j) = c_ij;
      }
    }
  }

  /**
   * Return a matrix-matrix product.
   * \see Matrix::mmult
   */
  Matrix
  mmult(const Matrix& B) const
  {
    Matrix C(n_rows(), B.n_cols());
    mmult(B, C);
    return C;
  }

  /**
   * Compute a transpose matrix-matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B} \\
   *                    = \sum_{k=0}^{n} a_{ki} b_{kj}, ~ \forall i, j
   * \f$.
   */
  void
  Tmmult(const Matrix& B, Matrix& C,
         const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_cols(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.");
    Assert(n_rows() == B.n_rows(), "Dimension mismatch error.");

    for (size_t i = 0; i < C.n_rows(); ++i)
      for (size_t j = 0; j < C.n_cols(); ++j)
      {
        value_type c_ij = adding? C(i, j) : 0.0;
        for (size_t k = 0; k < n_rows(); ++k)
          c_ij += coeffs[k][i] * B(k, j);
        C(i, j) = c_ij;
      }
  }

  /**
   * Return a transpose matrix-matrix product.
   * \see Matrix::Tmmult
   */
  Matrix
  Tmmult(const Matrix& B) const
  {
    Matrix C(n_cols(), B.n_cols());
    Tmmult(B, C);
    return C;
  }

  /**
   * Compute a matrix-transpose matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}^T \\
   *                    = \sum_{k=1}^{n} a_{ik} b_{jk}, ~ \forall i, j
   * \f$.
   */
  void
  mTmult(const Matrix& B, Matrix& C,
         const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_rows(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.");
    Assert(n_cols() == B.n_cols(), "Dimension mismatch error.");

    for (size_t i = 0; i < C.n_rows(); ++i)
    {
      const value_type* a_i = data(i);
      for (size_t j = 0; j < C.n_cols(); ++j)
      {
        const value_type* b_j = B.data(j);
        value_type c_ij = adding ? C(i, j) : 0.0;
        for (size_t k = 0; k < n_cols(); ++k)
          c_ij += a_i[k] * *b_j++;
        C(i, j) = c_ij;
      }
    }
  }

  /**
   * Return a matrix-transpose matrix multiplication.
   * \see Matrix::mTmult
   */
  Matrix
  mTmult(const Matrix& B) const
  {
    Matrix C(n_rows(), B.n_rows());
    mTmult(B, C);
    return C;
  }

  /**
   * Compute a transpose matrix-transpose matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T \\
   *                    = \sum_{k=1}^{n} a_{ki} b_{jk}, \forall i, j
   * \f$.
   */
  void
  TTmult(const Matrix& B, Matrix& C,
         const bool adding = false) const
  {
    Assert(!empty(), "Empty matrix error.");
    Assert(C.n_rows() == n_cols(), "Dimension mismatch error.");
    Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.");
    Assert(n_rows() == B.n_cols(), "Dimension mismatch error.");

    for (size_t i = 0; i < C.n_rows(); ++i)
      for (size_t j = 0; j < C.n_cols(); ++j)
      {
        const value_type* b_j = B.data(j);
        value_type c_ij = adding ? C(i, j) : 0.0;
        for (size_t k = 0; k < n_rows(); ++k)
          c_ij += coeffs[k][i] * *b_j++;
        C(i, j) = c_ij;
      }
  }

  /**
   * Return a transpose matrix-transpose matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T
   *                    = \sum_{k=0}^{n} a_{ki} b_{jk}, ~ \forall i, j
   * \f$.
   */
  Matrix
  TTmult(const Matrix& B) const
  {
    Matrix C(n_cols(), B.n_rows());
    TTmult(B, C);
    return C;
  }

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f$ \vec{y} = \boldsymbol{A} \vec{x} \\
   *             = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i
   * \f$.
   */
  void
  vmult(const Vector& x, Vector& y,
        const bool adding = false) const
  {
    Assert(x.size() == n_cols(), "Dimension mismatch error.");
    Assert(y.size() == n_rows(), "Dimension mismatch error.");

    for (size_t i = 0; i < n_rows(); ++i)
    {
      value_type v = adding ? y[i] : 0.0;
      const value_type* a_ij = data(i);
      for (size_t j = 0; j < n_cols(); ++j)
        v += *a_ij++ * x[j];
      y[i] = v;
    }
  }

  /**
   * Return a matrix-vector product.
   * \see Matrix::vmult
   */
  Vector
  vmult(const Vector& x) const
  {
    Vector y(n_rows());
    vmult(x, y);
    return y;
  }

  /**
   * Add a matrix-vector product to the destination vector.
   * \see Matrix::vmult
   */
  void
  vmult_add(const Vector& x, Vector& dst) const
  { vmult(x, dst, true); }

  /**
   * Compute a transpose matrix-vector product.
   * This is computed via
   * \f$ \vec{y} = \boldsymbol{A}^T \vec{x}
   *             = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall i
   * \f$.
   */
  void
  Tvmult(const Vector& x, Vector& y,
         const bool adding = false) const
  {
    Assert(x.size() == n_rows(), "Dimension mismatch error.");
    Assert(y.size() == n_cols(), "Dimension mismatch error.");

    if (!adding) y = 0.0;
    for (size_t i = 0; i < n_rows(); ++i)
    {
      const value_type x_i = x[i];
      const value_type* a_ij = data(i);
      for (size_t j = 0; j < n_cols(); ++j)
        y[j] += *a_ij++ * x_i;
    }
  }

  /**
   * Return a transpose matrix-vector product.
   * \see Matrix::Tvmult
   */
  Vector
  Tvmult(const Vector& x) const
  {
    Vector y(n_cols());
    Tvmult(x, y);
    return y;
  }

  /**
   * Add a transpose matrix-vector product to the destination vector.
   * \see Matrix::Tvmult
   */
  void
  Tvmult_add(const Vector& x, Vector& dst)
  { Tvmult(x, dst, true); }

  /**
   * Compute a matrix-vector product.
   * \see Matrix::vmult
   */
  Vector
  operator*(const Vector& x) const
  { return vmult(x); }

  /**
   * Return the transpose of the matrix.
   */
  Matrix
  transpose() const
  {
    Matrix A_T(n_cols(), n_rows());
    for (size_t i = 0; i < n_rows(); ++i)
    {
      const value_type* a_ij = coeffs[i].data();
      for (size_t j = 0; j < n_cols(); ++j)
        A_T[j][i] = *a_ij++;
    }
    return A_T;
  }

  // @}
  /** \name Printing Utilities */
  // @{

  /**
   * Return the matrix as a string.
   */
  std::string
  str(const bool scientific = false,
      const unsigned int precision = 3,
      const unsigned int width = 0) const
  {
    std::stringstream ss;
    print(ss, scientific, precision, width);
    return ss.str();
  }

  /**
   * Return the matrix as a string.
   *
   * \param os The output stream to print the matrix to.
   * \param scientific A flag for scientific notation.
   * \param precision The precision of the digits to display.
   * \param width The spacing between entries.
   */
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
      const value_type* a_ij = coeffs[i].data();
      for (uint64_t j = 0; j < n_cols(); ++j)
        os << std::setw(w) << *a_ij++;
      os << std::endl;
    }
    os << std::endl;
    os.flags(old_flags);
    os.precision(old_precision);
  }

  // @}

private:
  static bool
  valid_dimensions(const std::vector<std::vector<value_type>>& A)
  {
    size_t m = A.front().size();
    for (const auto& row : A)
      if (row.size() != m) return false;
    return true;
  }

  static bool
  valid_dimensions(const std::vector<Vector>& A)
  {
    size_t m = A.front().size();
    for (const auto& row : A)
      if (row.size() != m) return false;
    return true;
  }

};

/*-------------------- Methods -------------------- */

/**
 * Add two matrices together.
 * \see Matrix::operator+=
 */
inline Matrix
operator+(const Matrix& A, const Matrix& B)
{ return Matrix(A) += B; }


/**
 * Subtract two matrices.
 * \see Matrix::operator-=
 */
inline Matrix
operator-(const Matrix& A, const Matrix& B)
{ return Matrix(A) -= B; }


/**
 * Return a matrix-matrix multiplication.
 * \see Matrix::operator*(const Matrix&) Matrix::mmult
 */
inline Matrix
operator*(const Matrix& A, const Matrix& B)
{ return A.mmult(B); }


/**
 * Insert a matrix into an output stream.
 */
inline std::ostream&
operator<<(std::ostream& os, const Matrix& A)
{ return os << A.str(); }

}

#endif //MATRIX_H

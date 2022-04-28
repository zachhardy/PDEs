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

  /** Construct a square matrix of dimension \p n. */
  explicit Matrix(const size_t n) : m_data(n, std::vector<value_type>(n)) {}

  /** Construct a square matrix of dimension \p n set to \p value. */
  explicit Matrix(const size_t n, const value_type value)
    : m_data(n, std::vector<value_type>(n, value))
  {}

  /** Construct a matrix with \p n_rows and \p n_cols. */
  explicit Matrix(const size_t n_rows, const size_t n_cols)
    : m_data(n_rows, std::vector<value_type>(n_cols))
  {}

  /** Construct a matrix with \p n_rows and \p n_cols set to \p value. */
  explicit Matrix(const size_t n_rows,
                  const size_t n_cols,
                  const value_type value)
    : m_data(n_rows, std::vector<value_type>(n_cols, value))
  {}

  /** Copy constructor. */
  Matrix(const Matrix& other) : m_data(other.m_data) {}

  /** Move constructor. */
  Matrix(Matrix&& other) :  m_data(std::move(other.m_data)) {}

  /** Copy construction from an STL vector. */
  Matrix(const std::vector<std::vector<value_type>>& other);

  /** Move construction from an STL vector. */
  Matrix(std::vector<std::vector<value_type>>&& other);

  /** Construction from a nested initializer list. */
  Matrix(std::initializer_list<std::initializer_list<value_type>> list);

  /** Copy assignment operator. */
  Matrix& operator=(const Matrix& other);

  /** Move assignment operator. */
  Matrix& operator=(Matrix&& other);

  /** Copy assignment from an STL vector. */
  Matrix& operator=(const std::vector<std::vector<value_type>>& other);

  /** Move assignment from an STL vector. */
  Matrix& operator=(std::vector<std::vector<value_type>>&& other);

public:
  /** \name Accessors */
  /** @{ */

  /** Read/write access for row \p i. */
  std::vector<value_type>& operator[](const size_t i) { return m_data[i]; }

  /** Read only access for row \p i. */
  std::vector<value_type> operator[](const size_t i) const { return m_data[i]; }

  /** Read/write access for row \p i. */
  std::vector<value_type>& operator()(const size_t i) { return m_data[i]; }

  /** Read only access for row \p i. */
  std::vector<value_type> operator()(const size_t i) const { return m_data[i]; }

  /** Read/write access for row \p i with bounds checking. */
  std::vector<value_type>& at(const size_t i) { return m_data.at(i); }

  /** Read only access for row \p i with bounds checking. */
  std::vector<value_type> at(const size_t i) const { return m_data.at(i); }

  /** Read/write access for row \p i and column \p j. */
  value_type& operator()(const size_t i, const size_t j)
  { return m_data[i][j]; }

  /** Read only access for row \p i and column \p j. */
  value_type operator()(const size_t i, const size_t j) const
  { return m_data[i][j]; }

  /** Read/write access for row \p i and column \p j with bounds checking. */
  value_type& at(const size_t i, const size_t j) { return m_data.at(i).at(j); }

  /** Read only access for row \p i and column \p j with bounds checking. */
  value_type at(const size_t i, const size_t j) const
  { return m_data.at(i).at(j); }

  /** Read/write access to the <tt>i</tt>'th diagonal element. */
  value_type& diagonal(const size_t i) { return m_data[i][i]; }

  /** Read access to the <tt>i</tt>'th diagonal element. */
  value_type diagonal(const size_t i) const { return m_data[i][i]; }

  /** Return the diagonal of the matrix. */
  std::vector<value_type> diagonal() const;

  /** Access the underlying matrix data. */
  std::vector<value_type>* data() { return m_data.data(); }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /** Clear the matrix. */
  void clear() { m_data.clear(); }

  /** Remove the last row of the matrix. */
  void pop_back() { m_data.pop_back(); }

  /** Add a new row from \p values to the back of the matrix. */
  void push_back(const std::vector<value_type>& values)
  { m_data.push_back(values); }

  /** Add a new row from \p values in place to the back. */
  void emplace_back(const std::vector<value_type>& values)
  { m_data.emplace_back(values); }

  /** Resize to dimension \p n with default new values. */
  void resize(const size_t n)
  { m_data.resize(n, std::vector<value_type>(n)); }

  /** Resize to dimension \p n, setting new elements to \p value. */
  void resize(const size_t n, const value_type value)
  { m_data.resize(n, std::vector<value_type>(n, value)); }

  /** Resize to \p n_rows and \p n_cols with default new values. */
  void resize(const size_t n_rows, const size_t n_cols)
  { m_data.resize(n_rows, std::vector<value_type>(n_cols)); }

  /** Resize to \p n_rows and \p n_cols, setting new elements to \p value. */
  void resize(const size_t n_rows, const size_t n_cols, const value_type value)
  { m_data.resize(n_rows, std::vector<value_type>(n_cols, value)); }

  /** Swap the elements of two rows. */
  void swap_row(const size_t i1, const size_t i2)
  { m_data[i1].swap(m_data[i2]); }

  /** Swap the elements of a row and another vector. */
  void swap_row(const size_t i, Vector<value_type>& other);

  /** Swap the elements of a row with an STL vector. */
  void swap_row(const size_t i, std::vector<value_type>& other);

  /** Swap the elements of two columns. */
  void swap_column(const size_t j1, const size_t j2);

  /** Swap the elements of a column and another vector. */
  void swap_column(const size_t j, Vector<value_type>& other);

  /** Swap the elements of a column with an STL vector. */
  void swap_column(const size_t j, std::vector<value_type>& other);

  /** Swap the elements of this matrix with another. */
  void swap(Matrix& other) { m_data.swap(other.m_data); }

  /** Swap the elements of this matrix and an STL vector. */
  void swap(std::vector<std::vector<value_type>>& other) { m_data.swap(other); }

  /** Set the diagonal of the matrix. */
  void set_diagonal(const Vector<value_type>& diagonal);

  /** Set the diagonal of the matrix with an STL vector. */
  void set_diagonal(const std::vector<value_type>& diagonal)
  { return this->set_diagonal(Vector(diagonal)); }

  /** Set the diagonal of the matrix with a fixed scalar value. */
  void set_diagonal(const value_type value);

  /** @} */
  /** \name Information */
  /** @{ */

  /** Return the number of rows. */
  size_t n_rows() const { return m_data.size(); }

  /** Return the number of columns. */
  size_t n_cols() const { return m_data.front().size(); }

  /** Return the number of elements <tt> n_rows * n_cols </tt>. */
  size_t size() const { return this->n_rows() * this->n_cols(); }

  /** Return the number of non-zeros. */
  size_t n_nonzero_elements() const;

  /** Return whether the matrix is empty. */
  bool empty() const noexcept { return m_data.empty(); }

  /** @} */
  /** \name Iterators */
  /** @{ */

  /** Mutable iterator over rows at the start of the matrix. */
  typename std::vector<std::vector<value_type>>::iterator
  begin() { return m_data.begin(); }

  /** Mutable iterator one past the last row of the matrix. */
  typename std::vector<std::vector<value_type>>::iterator
  end() { return m_data.end(); }

  /** Mutable iterator at the start of row \p i. */
  typename std::vector<value_type>::iterator
  begin(const size_t i) { return m_data[i].begin(); }

  /** Mutable iterator one past the end of row \p i. */
  typename std::vector<value_type>::iterator
  end(const size_t i) { return m_data[i].end(); }

  /** Constant iterator over rows at the start of the matrix. */
  typename std::vector<std::vector<value_type>>::const_iterator
  cbegin() const { return m_data.cbegin(); }

  /** Constant iterator one past the last row of the matrix. */
  typename std::vector<std::vector<value_type>>::const_iterator
  cend() const { return m_data.cend(); }

  /** Constant iterator at the start of the row \p i. */
  typename std::vector<value_type>::const_iterator
  cbegin(const size_t i) const { return m_data[i].cbegin(); }

  /** Constant iterator one past the end of row \p i. */
  typename std::vector<value_type>::const_iterator
  cend(const size_t i) const { return m_data[i].cend(); }


  /** @} */
  /** \name Scalar Operations */
  /** @{ */

  /** Element-wise negation. */
  Matrix operator-() const;

  /** Element-wise negation in-place. */
  Matrix& operator-();

  /** Element-wise multiplication by a scalar. */
  Matrix operator*(const value_type value) const;

  /** Element-wise multiplication by a scalar in-place. */
  Matrix& operator*=(const value_type value);

  /** Element-wise division by a scalar. */
  Matrix operator/(const value_type value) const;

  /** Element-wise division by a scalar in-place. */
  Matrix& operator/=(const value_type value);

  /** @} */
  /** \name Linear Algebra Operations */
  /** @{ */

   /** Element-wise addition of two matrices. */
  Matrix operator+(const Matrix& other) const;

  /** Element-wise addition of two matrices in-place. */
  Matrix& operator+=(const Matrix& other);

  /** Element-wise subtraction of two matrices. */
  Matrix operator-(const Matrix& other) const;

  /** Element-wise subtraction of two matrices in-place. */
  Matrix& operator-=(const Matrix& other);

  /**
   * Return the matrix-matrix product.
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B} \\
   *     c_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  Matrix operator*(const Matrix& other) const;

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector<value_type> operator*(const Vector<value_type>& x) const;

  /**
   * Return the transpose of a matrix.
   * This is computed via
   * \f[ \boldsymbol{B} = \boldsymbol{A}^T \\
   *     b_{ij} = a_{ji}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  Matrix transpose() const;

  /** @} */
  /** \name Print Utilities */
  /** @{ */

  /** Return the matrix as a string. */
  std::string to_string() const;

  /** Print the vector to `std::cout`. */
  void print() const { std::cout << this->to_string(); }

  /** @} */
private:

  /** Check STL matrix inputs. */
  static void
  validate_stl_input(const std::vector<std::vector<value_type>>& other,
                     const std::string func_name);

};

/*-------------------- Inline Implementations --------------------*/

template<typename value_type>
inline Matrix<value_type>::
Matrix(const std::vector<std::vector<value_type>>& other)
{
  this->validate_stl_input(other, __FUNCTION__);
  m_data = other;
}


template<typename value_type>
inline Matrix<value_type>::
Matrix(std::vector<std::vector<value_type>>&& other)
{
  this->validate_stl_input(other, __FUNCTION__);
  m_data = std::move(other);
}


template<typename value_type>
inline Matrix<value_type>::
Matrix(std::initializer_list<std::initializer_list<value_type>> list)
{
  std::vector<std::vector<value_type>> other;
  for (auto& row : list)
    other.push_back(row);

  this->validate_stl_input(other, __FUNCTION__);
  m_data = std::move(other);
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator=(const Matrix<value_type>& other)
{
  m_data = other.m_data;
  return *this;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator=(Matrix<value_type>&& other)
{
  m_data = std::move(other.m_data);
  return *this;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator=(const std::vector<std::vector<value_type>>& other)
{
  this->validate_stl_input(other, __FUNCTION__);
  m_data = other;
  return *this;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator=(std::vector<std::vector<value_type>>&& other)
{
  this->validate_stl_input(other, __FUNCTION__);
  m_data = std::move(other);
  return *this;
}


template<typename value_type>
inline std::vector<value_type>
Matrix<value_type>::diagonal() const
{
  if (this->empty())
    return std::vector<value_type>();

  // Compute minimum dimension = diagonal size
  size_t min_dim = std::min(this->n_rows(), this->n_cols());

  // Populate the diagonal vector
  std::vector<value_type> v(min_dim);
  for (size_t i = 0; i < min_dim; ++i)
    v.push_back(m_data[i][i]);
  return v;
}


template<typename value_type>
inline void
Matrix<value_type>::swap_row(const size_t i, Vector<value_type>& other)
{
  if (other.size() != this->n_cols())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  for (size_t j = 0; j < this->n_cols(); ++j)
    std::swap(m_data[i][j], other[j]);
}


template<typename value_type>
inline void
Matrix<value_type>::swap_row(const size_t i, std::vector<value_type>& other)
{
  if (other.size() != this->n_cols())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }
  m_data[i].swap(other);
}


template<typename value_type>
inline void
Matrix<value_type>::swap_column(const size_t j1, const size_t j2)
{
  for (size_t i = 0; i < this->n_rows(); ++i)
    std::swap(m_data[i][j1], m_data[i][j1]);
}


template<typename value_type>
inline void
Matrix<value_type>::swap_column(const size_t j, Vector<value_type>& other)
{
  if (other.size() != this->n_rows())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  for (size_t i = 0; i < this->n_rows(); ++i)
    std::swap(m_data[i][j], other[i]);
}


template<typename value_type>
inline void
Matrix<value_type>::swap_column(const size_t j, std::vector<value_type>& other)
{
  if (other.size() != this->n_rows())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  for (size_t i = 0; i < this->n_rows(); ++i)
    std::swap(m_data[i][j], other[i]);
}


template<typename value_type>
inline void
Matrix<value_type>::set_diagonal(const Vector <value_type>& diagonal)
{
  // If empty, define a square matrix with diagonal
  if (this->empty())
  {
    this->resize(diagonal.size());
    for (size_t i = 0; i < diagonal.size(); ++i)
      m_data[i][i] = diagonal[i];
  }

  // Handle initialized matrices
  else
  {
    size_t min_dim = std::min(this->n_rows(), this->n_cols());
    if (diagonal.size() != min_dim)
    {
      std::stringstream err;
      err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
      throw std::length_error(err.str());
    }

    for (size_t i = 0; i < min_dim; ++i)
      m_data[i][i] = diagonal[i];
  }
}


template<typename value_type>
inline void Matrix<value_type>::set_diagonal(const value_type value)
{
  if (this->empty())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": "
        << "Cannot set an empty matrix with a scalar value.";
    throw std::runtime_error(err.str());
  }

  size_t min_dim = std::min(this->n_rows(), this->n_cols());
  for (size_t i = 0; i < min_dim; ++i)
    m_data[i][i] = value;
}


template<typename value_type>
inline size_t Matrix<value_type>::n_nonzero_elements() const
{
  size_t nnz = 0;
  for (const auto& row : m_data)
    for (const auto& entry : row)
      if (entry != 0.0) ++nnz;
  return nnz;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::operator-() const
{
  Matrix m(m_data);
  for (auto& row : m)
    for (auto& elem : row)
      elem = -elem;
  return m;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator-()
{
  for (auto& row : m_data)
    for (auto& elem : row)
      elem = -elem;
  return *this;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::operator*(const value_type value) const
{
  Matrix m(m_data);
  for (auto& row : m)
    for (auto& elem : row)
      elem *= value;
  return m;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator*=(const value_type value)
{
  for (auto& row : m_data)
    for (auto& elem : row)
      elem *= value;
  return *this;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::operator/(const value_type value) const
{
  if (value == 0.0)
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Zero division encountered.";
    throw std::runtime_error(err.str());
  }

  Matrix m(m_data);
  for (auto& row : m)
    for (auto& elem : row)
      elem /= value;
  return m;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator/=(const value_type value)
{
  if (value == 0.0)
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Zero division encountered.";
    throw std::runtime_error(err.str());
  }

  for (auto& row : m_data)
    for (auto& elem : row)
      elem /= value;
  return *this;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::operator+(const Matrix<value_type>& other) const
{
  if (this->n_rows() != other.n_rows() or
      this->n_cols() != other.n_cols())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  Matrix m(this->n_rows(), this->n_cols());
  for (size_t i = 0; i < m.n_rows(); ++i)
    for (size_t j = 0; j < m.n_cols(); ++j)
      m[i][j] = m_data[i][j] + other[i][j];
  return m;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator+=(const Matrix<value_type>& other)
{
  if (this->n_rows() != other.n_rows() or
      this->n_cols() != other.n_cols())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  for (size_t i = 0; i < this->n_rows(); ++i)
    for (size_t j = 0; j < this->n_cols(); ++j)
      m_data[i][j] += other[i][j];
  return *this;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::operator-(const Matrix<value_type>& other) const
{
  if (this->n_rows() != other.n_rows() or
      this->n_cols() != other.n_cols())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  Matrix m(this->n_rows(), this->n_cols());
  for (size_t i = 0; i < m.n_rows(); ++i)
    for (size_t j = 0; j < m.n_cols(); ++j)
      m[i][j] = m_data[i][j] - other[i][j];
  return m;
}


template<typename value_type>
inline Matrix<value_type>&
Matrix<value_type>::operator-=(const Matrix<value_type>& other)
{
  if (this->n_rows() != other.n_rows() or
      this->n_cols() != other.n_cols())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  for (size_t i = 0; i < this->n_rows(); ++i)
    for (size_t j = 0; j < this->n_cols(); ++j)
      m_data[i][j] -= other[i][j];
  return *this;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::operator*(const Matrix<value_type>& other) const
{
  if (this->n_cols() != other.n_rows())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  Matrix m(this->n_rows(), other.n_cols(), 0.0);
  for (size_t i = 0; i < m.n_rows(); ++i)
  {
    for (size_t j = 0; j < m.n_cols(); ++j)
    {
      value_type value = 0.0;
      for (size_t k = 0; k < other.n_rows(); ++k)
        value += m_data[i][k] * other[k][j];
      m[i][j] += value;
    }
  }
  return m;
}


template<typename value_type>
inline Vector<value_type>
Matrix<value_type>::operator*(const Vector<value_type>& x) const
{
  if (this->n_cols() != x.size())
  {
    std::stringstream err;
    err << "Matrix::" << __FUNCTION__ << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

  Vector v(this->n_rows(), 0.0);
  for (size_t i = 0; i < this->n_rows(); ++i)
    for (size_t j = 0; j < this->n_cols(); ++j)
      v[i] += m_data[i][j] * x[j];
  return v;
}


template<typename value_type>
inline Matrix<value_type>
Matrix<value_type>::transpose() const
{
  Matrix m(this->n_cols(), this->n_rows(), 0.0);
  for (size_t i = 0; i < m.n_rows(); ++i)
    for (size_t j = 0; j < m.n_cols(); ++j)
      m[i][j] = m_data[j][i];
  return m;
}


template<typename value_type>
inline std::string Matrix<value_type>::to_string() const
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


template<typename value_type>
inline void
Matrix<value_type>::
validate_stl_input(const std::vector<std::vector<value_type>>& other,
                   const std::string func_name)
{
  bool valid = true;
  size_t ref_val = other.front().size();
  for (auto& row : other)
  {
    if (row.size() != ref_val )
    { valid = false; break; }
  }

  // Throw error if not a valid matrix
  if (not valid)
  {
    std::stringstream err;
    err << "Matrix::" << func_name << ": "
        << "Invalid STL input encountered. All columns (inner STL vectors) "
        << "must have the same size.";
    throw std::runtime_error(err.str());
  }
}


/** Element-wise multiplication by a scalar value. */
template<typename value_type>
inline Matrix<value_type>
operator*(const value_type value, const Matrix<value_type>& A)
{ return A * value; }


/** Multiply a matrix by a scalar value. */
template<typename value_type>
inline Matrix<value_type>
multiply(const Matrix<value_type>& A, const value_type value)
{ return value * A; }


/** Multiply a matrix by a vector. */
template<typename value_type>
inline Vector<value_type>
multiply(const Matrix<value_type> A, const Vector<value_type> x)
{ return A * x; }


/** Multiply a matrix by a matrix. */
template<typename value_type>
inline Vector<value_type>
multiply(const Matrix<value_type> A, const Matrix<value_type> B)
{ return A * B; }

}
#endif //MATRIX_H

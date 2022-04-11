#ifndef MATRIX_H
#define MATRIX_H

#include "Math/vector.h"

#include <cmath>
#include <vector>
#include <iomanip>

class Matrix
{
private:
  /// The underlying matrix data.
  std::vector<std::vector<double>> m_data;

public:
  using iterator = std::vector<std::vector<double>>::iterator;
  using const_iterator = std::vector<std::vector<double>>::const_iterator;
  using reverse_iterator = std::vector<std::vector<double>>::reverse_iterator;
  using const_reverse_iterator =
      std::vector<std::vector<double>>::const_reverse_iterator;

public:
  Matrix() = default;  ///< Default constructor.

  /// Construct a matrix of dimension \p n set to default.
  explicit Matrix(size_t n) : m_data(n, std::vector<double>(n)) {}
  /// Construct a matrix of dimension \p n set to \p value.
  explicit Matrix(size_t n, const double value)
    : m_data(n, std::vector<double>(n, value))
  {}
  /// Construct a matrix with \p n_rows and \p n_cols set to default.
  explicit Matrix(size_t n_rows, size_t n_cols)
    : m_data(n_rows, std::vector<double>(n_cols))
  {}
  /// Construct a matrix with \p n_rows and \p n_cols set to \p value.
  explicit Matrix(size_t n_rows, size_t n_cols, const double value)
    : m_data(n_rows, std::vector<double>(n_cols, value))
  {}
  /// Copy constructor.
  Matrix(const Matrix& other) : m_data(other.m_data) {}
  /// Move constructor
  Matrix(Matrix&& other) :  m_data(std::move(other.m_data)) {}
  /// Construct a matrix from an STL vector.
  Matrix(const std::vector<std::vector<double>>& other)
  {
    if (not this->valid_stl(other))
      this->invalid_stl_error(__FUNCTION__);
    m_data = other;
  }
  /// Construct a matrix from an initialized list.
  Matrix(std::initializer_list<std::initializer_list<double>> list)
  {
    std::vector<std::vector<double>> other;
    for (auto& row : list)
      other.push_back(row);

    if (not this->valid_stl(other))
      this->invalid_stl_error(__FUNCTION__);

    m_data = std::move(other);
  }

  /// Copy assignment operator.
  Matrix& operator=(const Matrix& other)
  {
    m_data = other.m_data;
    return *this;
  }
  /// Move assignment operator.
  Matrix& operator=(Matrix&& other)
  {
    m_data = std::move(other.m_data);
    return *this;
  }
  /// Copy assignment from an STL vector
  Matrix& operator=(const std::vector<std::vector<double>>& other)
  {
    if (not this->valid_stl(other))
      this->invalid_stl_error(__FUNCTION__);

    m_data = other;
    return *this;
  }
  /// Move assignment from an STL vector
  Matrix& operator=(std::vector<std::vector<double>>&& other)
  {
    if (not this->valid_stl(other))
      this->invalid_stl_error(__FUNCTION__);
    m_data = std::move(other);
    return *this;
  }

  ~Matrix() = default; ///< Default destructor.

public:
  /** \name Access Operators */
  /** @{ */

  /// Read/write access for row \p row.
  std::vector<double>& operator[](const size_t row) { return m_data[row]; }
  /// Read only access for row \p row.
  std::vector<double> operator[](const size_t row) const { return m_data[row]; }

  /// Read/write access for row \p row with bounds checking.
  std::vector<double>& at(const size_t row) { return m_data.at(row); }
  /// Read only access for row \p row with bounds checking.
  std::vector<double> at(const size_t row) const { return m_data.at(row); }

  /// Read/write access for row \p row and column \p col with bounds checking.
  double& at(const size_t row, const size_t col)
  {
    return m_data.at(row).at(col);
  }
  /// Read only access for row \p row and column \p col with bounds checking.
  double at(const size_t row, const size_t col) const
  {
    return m_data.at(row).at(col);
  }

  /// Access the underlying matrix data.
  std::vector<double>* data() { return m_data.data(); }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /// Clear the matrix
  void clear() { m_data.clear(); }

  /// Resize the matrix to dimension \p n with default new values.
  void resize(const size_t n)
  {
    m_data.resize(n);
    for (auto& row : m_data)
      row.resize(n);
  }
  /// Resize the matrix to dimension \p n, setting new elements to \p value.
  void resize(const size_t n, const double value)
  {
    m_data.resize(n);
    for (auto& row : m_data)
      row.resize(n, value);
  }
  /// Resize the matrix to \p n_rows and \p n_cols with default new values.
  void resize(const size_t n_rows, const size_t n_cols)
  {
    m_data.resize(n_rows);
    for (auto& row : m_data)
      row.resize(n_cols);
  }
  /// Resize the matrix to \p n_rows and \p n_cols, setting new elements to \p value.
  void resize(const size_t n_rows, const size_t n_cols, const double value)
  {
    m_data.resize(n_rows);
    for (auto& row : m_data)
      row.resize(n_cols, value);
  }

  /** @} */
  /** \name Memory */
  /** @{ */

  /// Reserve storage for a matrix with dimension \p n.
  void reserve(const size_t n) { return m_data.reserve(n); }

  /// Return the number of rows.
  size_t n_rows() const { return m_data.size(); }
  /// Return the number of columns.
  size_t n_cols() const { return m_data.front().size(); }
  /// Return the number of elements in the matirx.
  size_t size() const { return this->n_rows() * this->n_cols(); }

  /// Get the number of rows that can be held in allocated memory.
  size_t capacity() const { return m_data.capacity(); }

  /// Return whether the matrix is empty.
  bool empty() const { return m_data.empty(); }

  /** @} */
  /** \name Iterators */
  /** @{ */

  /// Iterator to the first element.
  iterator begin() { return m_data.begin(); }
  /// Iterator to the last element.
  iterator end() { return m_data.end(); }

  /// Constant iterator to the first element.
  const_iterator cbegin() const { return m_data.cbegin(); }
  /// Constant iterator to the last element.
  const_iterator  cend() const { return m_data.cend(); }

  /// Reverse iterator to the first element.
  reverse_iterator rbegin() { return m_data.rbegin(); }
  /// Reverse iterator to the last element.
  reverse_iterator rend() { return m_data.rend(); }

  /// Constant reverse iterator to the first element.
  const_reverse_iterator crbegin() { return m_data.crbegin(); }
  /// Constant reverse iterator to the last element.
  const_reverse_iterator crend() { return m_data.crend(); }

  /** @} */
public:
  /** \name Scalar Operations */
  /** @{ */

  /**
   * \brief Element-wise negation
   * \f[ \boldsymbol{B} = -\boldsymbol{A} \f]
   */
  Matrix operator-() const { return -Matrix(m_data); }
  /// See \ref operator-() const
  Matrix& operator-()
  {
    for (auto& row : m_data)
      for (auto& v : row)
        v = -v;
    return *this;
  }

  /**
   * \brief Element-wise multiplication by a scalar value.
   * \f[ \boldsymbol{B} = \boldsymbol{A} \alpha \f]
   */
  Matrix operator*(const double value) const
  {
    Matrix m(m_data);
    for (auto& row : m)
      for (auto& entry : row)
        entry *= value;
    return m;
  }
  /// See \ref operator*(const double value) const
  Matrix& operator*=(const double value)
  {
    for (auto& row : m_data)
      for (auto& entry : row)
        entry *= value;
    return *this;
  }

  /**
   * \brief Element-wise division by a scalar value.
   * \f[ \boldsymbol{B} = \frac{\boldsymbol{A}}{\alpha} \f]
   */
  Matrix operator/(const double value) const
  {
    if (value == 0.0)
      this->zero_division_error(__FUNCTION__);

    Matrix m(m_data);
    for (auto& row : m)
      for (auto& entry : row)
        entry /= value;
    return m;
  }
  /// See \ref operator/(const double value) const
  Matrix& operator/=(const double value)
  {
    if (value == 0.0)
      this->zero_division_error(__FUNCTION__);

    for (auto& row : m_data)
      for (auto& entry : row)
        entry /= value;
    return *this;
  }

  /** @} */
  /** \name Linear Algebra Operations */
  /** @{ */

  /**
   * \brief Element-wise addition of two matrices.
   * \f[ \boldsymbol{C} = \boldsymbol{A} + \boldsymbol{B} \f]
   */
  Matrix operator+(const Matrix& other) const
  {
    if (this->n_rows() != other.n_rows() or
        this->n_cols() != other.n_cols())
      this->mismatched_size_error(__FUNCTION__);

    Matrix m(this->n_rows(), this->n_cols());
    for (size_t i = 0; i < m.n_rows(); ++i)
      for (size_t j = 0; j < m.n_cols(); ++j)
        m[i][j] = m_data[i][j] + other[i][j];
    return m;
  }
  /// See \ref operator+(const Matrix& other) const
  Matrix& operator+=(const Matrix& other)
  {
    if (this->n_rows() != other.n_rows() or
        this->n_cols() != other.n_cols())
      this->mismatched_size_error(__FUNCTION__);

    for (size_t i = 0; i < this->n_rows(); ++i)
      for (size_t j = 0; j < this->n_cols(); ++j)
        m_data[i][j] += other[i][j];
    return *this;
  }

  /**
   * \brief Element-wise subtraction of two matrices.
   * \f[ \boldsymbol{C} = \boldsymbol{A} - \boldsymbol{B} \f]
   */
  Matrix operator-(const Matrix& other) const
  {
    if (this->n_rows() != other.n_rows() or
        this->n_cols() != other.n_cols())
      this->mismatched_size_error(__FUNCTION__);

    Matrix m(this->n_rows(), this->n_cols());
    for (size_t i = 0; i < m.n_rows(); ++i)
      for (size_t j = 0; j < m.n_cols(); ++j)
        m[i][j] = m_data[i][j] - other[i][j];
    return m;
  }
  /// See \ref operator-(const Matrix& other) const
  Matrix& operator-=(const Matrix& other)
  {
    if (this->n_rows() != other.n_rows() or
        this->n_cols() != other.n_cols())
      this->mismatched_size_error(__FUNCTION__);

    for (size_t i = 0; i < this->n_rows(); ++i)
      for (size_t j = 0; j < this->n_cols(); ++j)
        m_data[i][j] -= other[i][j];
    return *this;
  }

  /**
   * \brief Compute a matrix-matrix product.
   * \f[
   *    \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B} \\
   *    c_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  Matrix operator*(const Matrix& other) const
  {
    if (this->n_cols() != other.n_rows())
      this->mismatched_size_error(__FUNCTION__);

    Matrix m(this->n_rows(), other.n_cols(), 0.0);
    for (size_t i = 0; i < m.n_rows(); ++i)
    {
      for (size_t j = 0; m.n_cols(); ++j)
      {
        double value = 0.0;
        for (size_t k = 0; k < other.n_rows(); ++k)
          value += m_data[i][k] * other[k][j];
        m[i][j] += value;
      }
    }
    return m;
  }
  /// See \ref operator*(const Matrix& other) const
  Vector operator*(const Vector& vector) const
  {
    if (this->n_cols() != vector.size())
      this->mismatched_size_error(__FUNCTION__);

    Vector v(this->n_rows(), 0.0);
    for (size_t i = 0; i < this->n_rows(); ++i)
      for (size_t j = 0; j < this->n_cols(); ++j)
        v[i] += m_data[i][j] * vector[j];
    return v;
  }

  /** @} */
public:
  /** \name Miscellaneous Routines */
  /** @{ */

  /// Return the diagonal of the matrix.
  Vector diagonal() const
  {
    if (this->empty())
      return Vector();

    // Compute minimum dimension = diagonal size
    size_t min_dim = std::min(this->n_rows(), this->n_cols());

    // Populate the diagonal vector
    Vector v(min_dim);
    for (size_t i = 0; i < min_dim; ++i)
      v.push_back(m_data[i][i]);
    return v;
  }

  /// Set the diagonal of the matrix.
  void set_diagonal(const Vector& diagonal)
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
        this->mismatched_size_error(__FUNCTION__);

      for (size_t i = 0; i < min_dim; ++i)
        m_data[i][i] = diagonal[i];
    }
  }
  /// Set the diagonal of the matrix with a fixed scalar.
  void set_diagonal(const double value)
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

  /// Return the number of non-zeros.
  size_t nnz() const
  {
    if (this->empty())
      return 0;

    size_t count = 0;
    for (const auto& row : m_data)
      for (const auto& entry : row)
        if (entry != 0.0) ++count;
    return count;
  }

  /** @} */
  /** \name Print Utilities */
  /** @{ */

  /// Get the matrix as a string.
  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < m_data.size(); ++i)
    {
      ss << ((i == 0)? "[" : " [");
      for (size_t j = 0; j < m_data[i].size() - 1; ++j)
        ss << std::setw(8) << m_data[i][j] << " ";
      ss << m_data[i].back() << "]"
         << ((i == m_data.size()-1)? "]\n" : "\n");
    }
    return ss.str();
  }

  /// Print the vector to `std::cout`.
  void print() const { std::cout << this->to_string(); }


  /** @} */

private:

  /// Check STL matrix inputs.
  static bool valid_stl(const std::vector<std::vector<double>>& other)
  {
    bool valid = true;
    size_t ref_val = other.front().size();
    for (auto& row : other)
    {
      if (row.size() != ref_val )
      {
        valid = false;
        break;
      }
    }
    return valid;
  }

  /// Determine whether zero elements exist.
  bool has_zero_elements() const
  {
    bool has_zeros = false;
    for (const auto& row : m_data)
    {
      for (const auto& entry : row)
      {
        if (entry == 0.0) { has_zeros = true; break; }
      }
    }
    return has_zeros;
  }

  /// Invalid STL matrix error.
  static void invalid_stl_error(const std::string func_name)
  {
    std::stringstream err;
    err << "Matrix::" << func_name << ": "
        << "All inner STL vectors must be of the same size.";
    throw std::length_error(err.str());
  }

  /// Throw an error for division by zero.
  static void zero_division_error(const std::string func_name)
  {
    std::stringstream err;
    err << "Vector::" << func_name << ": Zero division encountered.";
    throw std::runtime_error(err.str());
  }

  /// Throw an error for mismatched sizes.
  static void mismatched_size_error(const std::string func_name)
  {
    std::stringstream err;
    err << "Matrix::" << func_name << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }
};


#endif //MATRIX_H

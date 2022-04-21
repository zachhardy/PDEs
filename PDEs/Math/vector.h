#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace math
{

/// A class representing a general vector.
struct Vector
{
public:
  typedef std::vector<double> vector;

private:
  /// The underlying vector data.
  vector m_data;

public:
  /// Default constructor.
  Vector() = default;
  /// Construct \p size elements set to default values.
  explicit Vector(const size_t size) : m_data(size) {}
  /// Construct \p size elements set to \p value.
  explicit Vector(const size_t size, const double value) : m_data(size, value) {}

  /// Copy constructor.
  Vector(const Vector& other) : m_data(other.m_data) {}
  /// Move constructor.
  Vector(Vector&& other) : m_data(std::move(other.m_data)) {}

  /// Copy construction using an STL vector.
  Vector(const vector& other) : m_data(other) {}
  /// Move construction using an STL vector.
  Vector(vector&& other) : m_data(std::move(other)) {}

  /// Construct using an initializer list.
  Vector(std::initializer_list<double> list) : m_data(list) {}

  /// Copy assignment operator.
  Vector& operator=(const Vector& other)
  {
    m_data = other.m_data;
    return *this;
  }
  /// Move assignment operator.
  Vector& operator=(Vector&& other)
  {
    m_data = std::move(other.m_data);
    return *this;
  }

  /// Copy assignment from an STL vector.
  Vector& operator=(const vector & other)
  {
    m_data = other;
    return *this;
  }
  /// Move assignment from an STL vector
  Vector& operator=(vector&& other)
  {
    m_data = std::move(other);
    return *this;
  }

  /// Assigment using an initializer list.
  Vector& operator=(std::initializer_list<double> list)
  {
    m_data = list;
    return *this;
  }

  /// Comparison operator.
  bool operator==(const Vector& other)
  {
    if (other.size() != this->size())
      this->mismatched_size_error(__FUNCTION__);

    bool is_equal = true;
    for (size_t i = 0; i < this->size(); ++i)
      if (other[i] != m_data[i]) { is_equal = false; break; }
    return is_equal;
  }

public:
  /** \name Access Operators */
  /** @{ */

  /// Read/write access for element \p index.
  double& operator[](const size_t index) { return m_data[index]; }
  /// Read only access for element \p index.
  double operator[](const size_t index) const { return m_data[index]; }

  /// Read/write access for element \p index with bounds checking.
  double& at(const size_t index) { return m_data.at(index); }
  /// Read only access for element \p index with bounds checking.
  double at(const size_t index) const { return m_data.at(index); }

  /// Read/write access for the first element.
  double& front() { return m_data.front(); }
  /// Read only access for the first element.
  double front() const { return m_data.front(); }

  /// Read/write access for the last element.
  double& back() { return m_data.back(); }
  /// Read only access for the last element.
  double back() const { return m_data.back(); }

  /// Access the underlying data.
  double* data() { return m_data.data(); }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /// Clear the elements.
  void clear() { m_data.clear(); }

  /// Add a new element set to \p value to the back.
  void push_back(const double value) { m_data.push_back(value); }

  /// Add a new element set to \p value in place to the back.
  void emplace_back(const double value) { m_data.emplace_back(value); }

  /// Remove the last element.
  void pop_back() { m_data.pop_back(); }

  /// Resize to \p new_size elements, setting new elements to default.
  void resize(const size_t new_size) { m_data.resize(new_size); }
  /// Resize to \p new_size elements, setting new elements to default.
  void resize(const size_t new_size, const double value)
  { m_data.resize(new_size, value); }

  /// Swap the elements of this vector and another.
  void swap(Vector& other) { m_data.swap(other.m_data); }
  /// Swap the elements of this vector and an STL vector.
  void swap(vector& other) { m_data.swap(other); }

  /** @} */
  /** \name Memory */
  /** @{ */

  /// Allocate memory for \p new_size elements.
  void reserve(const size_t new_size) { m_data.reserve(new_size); }

  /// Return the number of elements.
  size_t size() const { return m_data.size(); }

  /// Return the number of allocated elements.
  size_t capacity() const { return m_data.capacity(); }

  /// Return whether the vector is empty.
  bool empty() const { return m_data.empty(); }

  /** @} */
  /** \name Iterators */
  /** @{ */

  vector::iterator begin() { return m_data.begin(); }
  vector::iterator end() { return m_data.end(); }

  vector::const_iterator cbegin() const { return m_data.cbegin(); }
  vector::const_iterator cend() const { return m_data.cend(); }

  vector::reverse_iterator rbegin() { return m_data.rbegin(); }
  vector::reverse_iterator rend() { return m_data.rend(); }

  vector::const_reverse_iterator crbegin() const { return m_data.crbegin(); }
  vector::const_reverse_iterator crend() const { return m_data.crend(); }

  /** @} */
public:
  /** \name Scalar Operations */
  /** @{ */

  /**
   * \brief Element-wise negation.
   * \f[
   *    \vec{y} = -\vec{x} \\
   *    y_i = -x_i, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator-() const { return -Vector(m_data); }
  /// See \ref operator-() const
  Vector& operator-()
  {
    for (auto& entry : m_data)
      entry = -entry;
    return *this;
  }

  /**
   * \brief Element-wise multiplication by a scalar a.
   * \f[
   *    \vec{y} = \vec{x} c \\
   *    y_i = x_i c, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator*(const double c) const
  {
    Vector v(m_data);
    for (auto& entry : v)
      entry *= c;
    return v;
  }
  /// See \ref operator*(const double c) const
  Vector& operator*=(const double c)
  {
    for (auto& entry : m_data)
      entry *= c;
    return *this;
  }

  /**
   * \brief Element-wise division by a scalar value.
   * \f[
   *    \vec{y} = \frac{\vec{x}}{c} \\
   *    y_i = \frac{x_i}{c}, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator/(const double c) const
  {
    if (c == 0.0)
      this->zero_division_error(__FUNCTION__);

    Vector v(m_data);
    for (auto& entry : v)
      entry /= c;
    return v;
  }
  /// See \ref  operator/(const double c) const
  Vector& operator/=(const double c)
  {
    if (c == 0.0)
      this->zero_division_error(__FUNCTION__);

    for (auto& entry : m_data)
      entry /= c;
    return *this;
  }

  /** @} */
  /** \name Vector-Vector Operations */
  /** @{ */

  /**
   * \brief Element-wise addition of two vectors.
   * \f[
   *    \vec{z} = \vec{x} + \vec{y}  \\
   *    z_i = x_i + y_i, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator+(const Vector& other) const
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    Vector v(m_data.size());
    for (size_t i = 0; i < v.size(); ++i)
      v[i] = m_data[i] + other[i];
    return v;
  }
  /// See \ref operator+(const Vector& other) const
  Vector& operator+=(const Vector& other)
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] += other[i];
    return *this;
  }

  /**
   * \brief Element-wise subtraction of two vectors.
   * \f[
   *    \vec{z} = \vec{x} - \vec{y} \\
   *    z_i = x_i - y_i, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator-(const Vector& other) const
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    Vector v(m_data.size());
    for (size_t i = 0; i < v.size(); ++i)
      v[i] = m_data[i] - other[i];
    return v;
  }
  /// See \ref operator-(const Vector& other) const
  Vector& operator-=(const Vector& other)
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] -= other[i];
    return *this;
  }

  /**
   * \brief Element-wise multiplication of two vectors.
   * \f[
   *    \vec{z} = \vec{x} \vec{y} \\
   *    z_i = x_i y_i, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator*(const Vector& other) const
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    Vector v(m_data.size());
    for (size_t i = 0; i < v.size(); ++i)
      v[i] = m_data[i] * other[i];
    return v;
  }
  /// See \ref operator*(const Vector& other) const
  Vector& operator*=(const Vector& other)
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] *= other[i];
    return *this;
  }

  /**
   * \brief Element-wise division of two vectors.
   * \f[
   *    \vec{z} = \frac{\vec{x}}{\vec{y}} \\
   *    z_i = \frac{x_i}{y_i} \hspace{0.25cm} \forall i
   * \f]
   */
  Vector operator/(const Vector& other) const
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);
    if (other.has_zero_elements())
      this->zero_division_error(__FUNCTION__);

    Vector v(m_data.size());
    for (size_t i = 0; i < v.size(); ++i)
      v[i] = m_data[i] / other[i];
    return v;
  }
  /// See \ref operator/(const Vector& other) const
  Vector& operator/(const Vector& other)
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);
    if (other.has_zero_elements())
      this->zero_division_error(__FUNCTION__);

    for (size_t i = 0; m_data.size(); ++i)
      m_data[i] /= other[i];
    return *this;
  }

  /**
   * \brief Return the dot product between this and another vector.
   * \f[ c = \vec{x} \cdot \vec{y} = \sum_i x_i y_i .\f]
   */
  double dot(const Vector& other) const
  {
    if (m_data.size() != other.size())
      this->mismatched_size_error(__FUNCTION__);

    double c = 0.0;
    for (size_t i = 0; i < m_data.size(); ++i)
      c += m_data[i] * other[i];
    return c;
  }

  /** @} */
  /** \name  Norms */
  /** @{ */

  /**
   * \brief Compute the \f$ \ell_\infty \f$-norm.
   * \f[ ||\vec{v}||_{\ell_\infty} = \max_i |v_i| \f]
   */
  double linf_norm() const
  {
    double norm = 0.0;
    for (const auto& v : m_data) if (std::fabs(v) > norm) norm = std::fabs(v);
    return norm;
  }

  /**
   * \brief Compute the \f$ \ell_1 \f$-norm.
   * \f[ ||\vec{v}||_{\ell_1} = \sum_i |v_i| \f]
   */
  double l1_norm() const
  {
    double norm = 0.0;
    for (const auto& v : m_data) norm += std::fabs(v);
    return norm;
  }

  /**
   * \brief Compute the \f$ \ell_2 \f$-norm.
   * \f[ ||\vec{v}||_{\ell_2} = \sqrt{ \sum_i |v_i|^2 } \f]
   */
  double l2_norm() const { return lp_norm(2.0); }

  /**
   * \brief Compute the \f$ \ell_{\ell_p} \f$-norm.
   * \f[ ||\vec{v}||_{\ell_p} = \left( \sum_i |v_i|^p \right)^{1/p} \f]
   */
  double lp_norm(const double p) const
  {
    double norm = 0.0;
    for (const auto& v : m_data) norm += std::pow(std::fabs(v), p);
    return std::pow(norm, 1.0/p);
  }

  /** @} */
  /** \name Vector Operations */
  /** @{ */

  /**
   * \brief Normalize this vector to unit length in place
   * \f[ \hat{v} = \frac{\vec{v}}{||\vec{v}||_{\ell_2}} \f]
   */
  Vector& normalize()
  {
    double length = this->l2_norm();
    if (length == 0.0) return *this;
    else return this->operator/=(length);
  }

  /// Element-wise absolute value in place.
  Vector& fabs()
  {
    for (auto& v : m_data) v = std::fabs(v);
    return *this;
  }

  /** @} */
public:
  /** /name Print Utilities */
  /** @{ */

  /// Get the vector as a string.
  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < m_data.size() - 1; ++i)
      ss << std::setprecision(6) << m_data[i] << " ";
    ss << std::setprecision(6) << m_data.back() << "]" << std::endl;
    return ss.str();
  }

  /// Print the vector to `std::cout`.
  void print() const { std::cout << this->to_string(); }

  /** @} */
private:

  /// Determine whether zero elements exist.
  bool has_zero_elements() const
  {
    bool has_zeros = false;
    for (const auto& v : m_data)
    {
      if (v == 0.0) { has_zeros = true; break; }
    }
    return has_zeros;
  }

  /// Throw an error for division by zero.
  static void zero_division_error(const std::string func_name)
  {
    std::stringstream err;
    err << "Vector::" << func_name << ": Zero division encountered.";
    throw std::runtime_error(err.str());
  }

  /// Throw an error for of mismatched sizes.
  static bool mismatched_size_error(const std::string func_name)
  {
    std::stringstream err;
    err << "Vector::" << func_name << ": Mismatched sizes encountered.";
    throw std::length_error(err.str());
  }

};

/*-------------------- Inline Implementations --------------------*/

/**
 * \brief Return the element-wise product of a vector and a scalar value.
 * \f[ \vec{y} = c \vec{x} \f]
 */
inline Vector operator*(const double c, const Vector& x) { return x * c; }

/**
 * \brief Return the dot product between two vectors
 * \f[ c = \vec{x} \cdot \vec{y} = \sum_i x_i y_i \f]
 */
inline double dot(const Vector& x, const Vector& y) { return x.dot(y); }

/**
 * \brief Return the \f$ \ell_\infty \f$-norm of a vector.
 * \f[ ||\vec{x}||_{\ell_\infty} = \max_i |x_i| \f]
 */
inline double linf_norm(const Vector& x) { return x.linf_norm(); }

/**
 * \brief Return the \f$ \ell_1 \f$ of a vector.
 * \f[ ||\vec{x}||_{\ell_1} = \sum_i |x_i| \f]
 */
inline double l1_norm(const Vector& x) { return x.l1_norm(); }

/**
 * \brief Return the \f$ \ell_2 \f$ of a vector.
 * \f[ ||\vec{x}||_{\ell_2} = \sqrt{\sum_i |x_i|^2} \f]
 */
inline double l2_norm(const Vector& x) { return x.l2_norm(); }

/**
 * \brief Return the \f$ \ell_{p} \f$ of a vector.
 * \f[ ||\vec{v}||_{\ell_p} = \left( \sum_i |v_i|^p \right)^{1/p} \f]
 */
inline double lp_norm(const Vector& x, const double p) { return x.lp_norm(p); }

/**
 * \brief Retrun the vector \p x normalized to its length. If the length of
 *        the vector is zero, a copy of the vector is returned.
 * \f[ \hat{y} = \vec{x} / ||\vec{x}||_{\ell_2} \f]
 */
inline Vector normalize(const Vector& x) { return Vector(x).normalize(); }

/**
 * \brief Return the absolute value of a vector.
 * \f[ \vec{y} = |\vec{x}| \f]
 */
inline Vector fabs(const Vector& x) { return Vector(x).fabs(); }

}
#endif //VECTOR_H

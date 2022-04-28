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

/** A class representing a general vector. */
template<typename value_type>
struct Vector
{
private:
  std::vector<value_type> m_data;

public:
  /** Default constructor. */
  Vector() = default;

  /** Construct a vector with \p n elements. */
  explicit Vector(const size_t n) : m_data(n) {}

  /** Construct a vector with \p n elements set to \p value */
  explicit Vector(const size_t n, const value_type value) : m_data(n, value) {}

  /** Copy constructor. */
  Vector(const Vector& other) : m_data(other.m_data) {}

  /** Move constructor. */
  Vector(Vector&& other) : m_data(std::move(other.m_data)) {}

  /** Copy construction using an STL vector. */
  Vector(const std::vector<value_type>& other) : m_data(other) {}

  /** Move construction using an STL vector. */
  Vector(std::vector<value_type>&& other) : m_data(std::move(other)) {}

  /** Construct using an initializer list. */
  Vector(std::initializer_list<value_type> list) : m_data(list) {}

  /** Copy assignment operator. */
  Vector& operator=(const Vector& other);

  /** Move assignment operator. */
  Vector& operator=(Vector&& other);

  /** Copy assignment from an STL vector. */
  Vector& operator=(const std::vector<value_type>& other);

  /** Move assignment from an STL vector. */
  Vector& operator=(std::vector<value_type>&& other);

  /** Assigment using an initializer list. */
  Vector& operator=(std::initializer_list<value_type>& list);

public:
  /** \name Access Operators */
  /** @{ */

  /** Read/write access for element \p i. */
  double& operator[](const size_t i) { return m_data[i]; }

  /** Read only access for element \p i. */
  double operator[](const size_t i) const { return m_data[i]; }

  /** Read/write access for element \p i with bounds checking. */
  double& at(const size_t i) { return m_data.at(i); }

  /** Read only access for element \p i with bounds checking. */
  double at(const size_t i) const { return m_data.at(i); }

  /** Read/write access for the first element. */
  double& front() { return m_data.front(); }

  /** Read only access for the first element. */
  double front() const { return m_data.front(); }

  /** Read/write access for the last element. */
  double& back() { return m_data.back(); }

  /** Read only access for the last element. */
  double back() const { return m_data.back(); }

  /** Access the underlying data. */
  double* data() { return m_data.data(); }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /** Clear the elements. */
  void clear() { m_data.clear(); }

  /** Add a new element set to \p value to the back. */
  void push_back(const value_type value) { m_data.push_back(value); }

  /** Add a new element set to \p value in place to the back. */
  void emplace_back(const value_type value) { m_data.emplace_back(value); }

  /** Remove the last element. */
  void pop_back() { m_data.pop_back(); }

  /** Resize to \p new_size elements, setting new elements to default. */
  void resize(const size_t new_size) { m_data.resize(new_size); }

  /** Resize to \p new_size elements, setting new elements to default. */
  void resize(const size_t new_size, const value_type value)
  { m_data.resize(new_size, value); }

  /** Swap the elements of this vector and another. */
  void swap(Vector& other) { m_data.swap(other.m_data); }

  /** Swap the elements of this vector and an STL vector. */
  void swap(std::vector<value_type>& other) { m_data.swap(other); }

  /** @} */
  /** \name Memory */
  /** @{ */

  /** Allocate memory for \p new_size elements. */
  void reserve(const size_t new_size) { m_data.reserve(new_size); }

  /** Return the number of elements. */
  size_t size() const { return m_data.size(); }

  /** Return whether the vector is empty. */
  bool empty() const { return m_data.empty(); }

  /** @} */
  /** \name Iterators */
  /** @{ */

  /** Mutable iterator at the start of the vector. */
  typename std::vector<value_type>::iterator
  begin() { return m_data.begin(); }

  /** Mutable iterator one past the end of the vector. */
  typename std::vector<value_type>::iterator
  end() { return m_data.end(); }

  /** Constant iterator at the start of the vector. */
  typename std::vector<value_type>::const_iterator
  cbegin() const { return m_data.cbegin(); }

  /** Constant iterator at the end of the vector. */
  typename std::vector<value_type>::const_iterator
  cend() const { return m_data.cend(); }

  /** @} */
public:
  /** \name Scalar Operations */
  /** @{ */

  /** Element-wise negation. */
  Vector operator-() const;

  /** Element-wise negation in-place. */
  Vector& operator-();

  /** Element-wise multiplication by a scalar. */
  Vector operator*(const value_type value) const;

  /** Element-wise multiplication by a scalar in-place. */
  Vector& operator*=(const value_type value);

  /** Element-wise division by a scalar. */
  Vector operator/(const value_type value) const;

  /** Element-wise division by a scalar in-place. */
  Vector& operator/=(const value_type value);

  /** @} */
  /** \name Vector-Vector Operations */
  /** @{ */

  /** Element-wise addition of two vectors. */
  Vector operator+(const Vector& other) const;

  /** Element-wise addition of two vectors in-place. */
  Vector& operator+=(const Vector& other);

  /** Element-wise subtraction of two vectors. */
  Vector operator-(const Vector& other) const;

  /** Element-wise subtraction of two vectors in-place. */
  Vector& operator-=(const Vector& other);

  /** Element-wise multiplication of two vectors. */
  Vector operator*(const Vector& other) const;

  /** Element-wise multiplication of two vectors in-place. */
  Vector& operator*=(const Vector& other);

  /** Element-wise division of two vectors. */
  Vector operator/(const Vector& other) const;

  /** Element-wise division of two vectors in-place. */
  Vector& operator/=(const Vector& other);

  /**
   * Return the dot product between this and another vector.
   * \f$ c = \vec{x} \cdot \vec{y} = \sum_i x_i y_i .\f$
   */
  value_type dot(const Vector& other) const;

  /** @} */
  /** \name  Norms */
  /** @{ */

  /** Compute the \f$ \ell_\infty \f$-norm.
   *  \f$ ||\vec{v}||_{\ell_\infty} = \max_i |v_i| \f$ */
  value_type linf_norm() const;

  /** Compute the \f$ \ell_1 \f$-norm.
   *  \f$ ||\vec{v}||_{\ell_1} = \sum_i |v_i| \f$ */
  value_type l1_norm() const;

  /* Compute the \f$ \ell_2 \f$-norm.
   * \f$ ||\vec{v}||_{\ell_2} = \sqrt{ \sum_i |v_i|^2 } \f$ */
  value_type l2_norm() const;

  /** Compute the \f$ \ell_{\ell_p} \f$-norm.
   *  \f$ ||\vec{v}||_{\ell_p} = \left( \sum_i |v_i|^p \right)^{1/p} \f$ */
  value_type lp_norm(const value_type p) const;

  /** @} */
  /** \name Vector Operations */
  /** @{ */

  /** Normalize this vector to unit length in place.
   *  \f$ \hat{v} = \frac{\vec{v}}{||\vec{v}||_{\ell_2}} \f$
   *  \note If the vector is uniformly zero, nothing is done.
   */
  Vector& normalize();

  /** Element-wise absolute value. */
  Vector& fabs();

  /** @} */
public:
  /** /name Print Utilities */
  /** @{ */

  /** Return the vector as a string. */
  std::string to_string() const;

  /** Print the vector to `std::cout`. */
  void print() const { std::cout << this->to_string(); }

  /** @} */
private:

  /** Return whether or not the vector has zero elements. */
  bool has_zero_elements() const;

  /** Throw an error for division by zero. */
  static void zero_division_error(const std::string func_name);

  /** Throw an error for of mismatched sizes. */
  static bool mismatched_size_error(const std::string func_name);

};

/*-------------------- Inline Implementations --------------------*/

template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator=(const Vector<value_type>& other)
{
  m_data = other.m_data;
  return *this;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator=(Vector<value_type>&& other)
{
  m_data = std::move(other.m_data);
  return *this;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator=(const std::vector<value_type>& other)
{
  m_data = other;
  return *this;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator=(std::vector<value_type>&& other)
{
  m_data = std::move(other);
  return *this;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator=(std::initializer_list<value_type>& list)
{
  m_data = list;
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator-() const
{
  Vector<value_type> v(m_data);
  for (auto& elem : v) elem = -elem;
  return v;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator-()
{
  for (auto& elem : m_data) elem = -elem;
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator*(const value_type value) const
{
  Vector<value_type> v(m_data);
  for (auto& entry : v) entry *= value;
  return v;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator*=(const value_type value)
{
  for (auto& elem : m_data) elem *= value;
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator/(const value_type value) const
{
  if (value == 0.0) this->zero_division_error(__FUNCTION__);
  Vector<value_type> v(m_data);
  for (auto& elem : v) elem /= value;
  return v;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator/=(const value_type value)
{
  if (value == 0.0) this->zero_division_error(__FUNCTION__);
  for (auto& elem : m_data) elem /= value;
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator+(const Vector<value_type>& other) const
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  Vector<value_type> v(m_data.size());
  for (size_t i = 0; i < v.size(); ++i)
    v[i] = m_data[i] + other[i];
  return v;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator+=(const Vector<value_type>& other)
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  for (size_t i = 0; i < m_data.size(); ++i)
    m_data[i] += other[i];
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator-(const Vector<value_type>& other) const
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  Vector v(m_data.size());
  for (size_t i = 0; i < v.size(); ++i)
    v[i] = m_data[i] - other[i];
  return v;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator-=(const Vector<value_type>& other)
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  for (size_t i = 0; i < m_data.size(); ++i)
    m_data[i] -= other[i];
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator*(const Vector<value_type>& other) const
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  Vector v(m_data.size());
  for (size_t i = 0; i < v.size(); ++i)
    v[i] = m_data[i] * other[i];
  return v;
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator*=(const Vector<value_type>& other)
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  for (size_t i = 0; i < m_data.size(); ++i)
    m_data[i] *= other[i];
  return *this;
}


template<typename value_type>
inline Vector<value_type>
Vector<value_type>::operator/(const Vector<value_type>& other) const
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


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::operator/=(const Vector<value_type>& other)
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);
  if (other.has_zero_elements())
    this->zero_division_error(__FUNCTION__);

  for (size_t i = 0; i < m_data.size(); ++i)
    m_data[i] /= other[i];
  return *this;
}


template<typename value_type>
inline value_type
Vector<value_type>::dot(const Vector<value_type>& other) const
{
  if (m_data.size() != other.size())
    this->mismatched_size_error(__FUNCTION__);

  value_type value = 0.0;
  for (size_t i = 0; i < m_data.size(); ++i)
    value += m_data[i] * other[i];
  return *this;
}


template<typename value_type>
inline value_type Vector<value_type>::linf_norm() const
{
  value_type norm = 0.0;
  for (const auto& elem : m_data)
    if (std::fabs(elem) > norm)
      norm = std::fabs(elem);
  return norm;
}


template<typename value_type>
inline value_type Vector<value_type>::l1_norm() const
{
  value_type norm = 0.0;
  for (const auto& elem : m_data)
    norm += std::fabs(elem);
  return norm;
}


template<typename value_type>
inline value_type Vector<value_type>::l2_norm() const
{
  value_type norm = 0.0;
  for (const auto& elem : m_data)
    norm += elem * elem;
  return std::sqrt(norm);
}


template<typename value_type>
inline value_type Vector<value_type>::lp_norm(const value_type p) const
{
  value_type norm = 0.0;
  for (const auto& elem : m_data)
    norm += std::pow(std::fabs(elem), p);
  return std::pow(norm, 1.0/p);
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::normalize()
{
  value_type length = this->l2_norm();
  if (length == 0.0) return *this;
  else return this->operator/=(length);
}


template<typename value_type>
inline Vector<value_type>&
Vector<value_type>::fabs()
{
  for (auto& elem : m_data) elem = std::fabs(elem);
  return *this;
}


template<typename value_type>
inline std::string Vector<value_type>::to_string() const
{
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < m_data.size() - 1; ++i)
    ss << std::setprecision(6) << m_data[i] << " ";
  ss << std::setprecision(6) << m_data.back() << "]\n";
  return ss.str();
}


template<typename value_type>
inline bool Vector<value_type>::has_zero_elements() const
{
  bool has_zeros = false;
  for (const auto& v : m_data)
    if (v == 0.0) { has_zeros = true; break; }
  return has_zeros;
}


template<typename value_type>
inline void
Vector<value_type>::zero_division_error(const std::string func_name)
{
  std::stringstream err;
  err << "Vector::" << func_name << ": Zero division encountered.";
  throw std::runtime_error(err.str());
}


template<typename value_type>
inline bool
Vector<value_type>::mismatched_size_error(const std::string func_name)
{
  std::stringstream err;
  err << "Vector::" << func_name << ": Mismatched sizes encountered.";
  throw std::length_error(err.str());
}


/** Element-wise multiplication by a scalar. */
template<typename value_type>
inline Vector<value_type> operator*(const value_type value,
                                    const Vector<value_type>& x)
{ return x * value; }


/** Return the dot product between two vectors. */
template<typename value_type>
inline value_type dot(const Vector<value_type>& x,
                      const Vector<value_type>& y)
{ x.dot(y); }


/** Return the \f$ \ell_\infty \f$-norm of a vector. */
template<typename value_type>
inline value_type linf_norm(const Vector<value_type>& x) { x.linf_norm(); }


/** Return the \f$ \ell_1 \f$-norm of a vector. */
template<typename value_type>
inline value_type l1_norm(const Vector<value_type>& x) { x.l1_norm(); }


/** Return the \f$ \ell_2 \f$-norm of a vector. */
template<typename value_type>
inline value_type l2_norm(const Vector<value_type>& x) { x.l2_norm(); }


/** Return the \f$ \ell_p \f$-norm of a vector. */
template<typename value_type>
inline value_type lp_norm(const Vector<value_type>& x,
                          const value_type p) { x.lp_norm(p); }


/** Return the direction (unit) vector of a vector. */
template<typename value_type>
inline Vector<value_type> normalize(const Vector<value_type>& x)
{ return Vector<value_type>(x).normalize(); }


/** Return the element-wise absolute value of a vector. */
template<typename value_type>
inline Vector<value_type> fabs(const Vector<value_type>& x)
{ return Vector<value_type>(x).fabs(); }

}
#endif //VECTOR_H

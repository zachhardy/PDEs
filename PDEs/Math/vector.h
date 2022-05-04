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

#include "exceptions.h"

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

  /** Construct a vector with \p n elements insert to \p value */
  explicit Vector(const size_t n, const value_type value = 0.0)
    : m_data(n, value)
  {}

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
  Vector& operator=(const Vector& other)
  {
    m_data = other.m_data;
    return *this;
  }

  /** Move assignment operator. */
  Vector& operator=(Vector&& other)
  {
    m_data = std::move(other.m_data);
    return *this;
  }

  /** Copy assignment from an STL vector. */
  Vector& operator=(const std::vector<value_type>& other)
  {
    m_data = other;
    return *this;
  }

  /** Move assignment from an STL vector. */
  Vector& operator=(std::vector<value_type>&& other)
  {
    m_data = std::move(other);
    return *this;
  }

  /** Assigment using an initializer list. */
  Vector& operator=(std::initializer_list<value_type>& list)
  {
    m_data = list;
    return *this;
  }

public:
  /** \name Access Operators */
  /** @{ */

  /** Read/write access for element \p i. */
  double& operator[](const size_t i)
  {
    return m_data[i];
  }

  /** Read only access for element \p i. */
  double operator[](const size_t i) const
  {
    return m_data[i];
  }

  /** Read/write access for element \p i with bounds checking. */
  double& at(const size_t i)
  {
    return m_data.at(i);
  }

  /** Read only access for element \p i with bounds checking. */
  double at(const size_t i) const
  {
    return m_data.at(i);
  }

  /** Read/write access for the first element. */
  double& front()
  {
    return m_data.front();
  }

  /** Read only access for the first element. */
  double front() const
  {
    return m_data.front();
  }

  /** Read/write access for the last element. */
  double& back()
  {
    return m_data.back();
  }

  /** Read only access for the last element. */
  double back() const
  {
    return m_data.back();
  }

  /** Access the underlying data. */
  double* data()
  {
    return m_data.data();
  }

  /** @} */
  /** \name Modifiers */
  /** @{ */

  /** Clear the elements. */
  void clear()
  {
    m_data.clear();
  }

  /** Add a new element insert to \p value to the back. */
  void push_back(const value_type value)
  {
    m_data.push_back(value);
  }

  /** Add a new element insert to \p value in place to the back. */
  void emplace_back(const value_type value)
  {
    m_data.emplace_back(value);
  }

  /** Remove the last element. */
  void pop_back()
  {
    m_data.pop_back();
  }

  /** Resize to \p new_size elements, setting new elements to default. */
  void resize(const size_t new_size, const value_type value = 0.0)
  {
    m_data.resize(new_size, value);
  }

  /** Swap the elements of this vector and another. */
  void swap(Vector& other)
  {
    m_data.swap(other.m_data);
  }

  /** Swap the elements of this vector and an STL vector. */
  void swap(std::vector<value_type>& other)
  {
    m_data.swap(other);
  }

  /** @} */
  /** \name Memory */
  /** @{ */

  /** Allocate memory for \p new_size elements. */
  void reserve(const size_t new_size)
  {
    m_data.reserve(new_size);
  }

  /** Return the number of elements. */
  size_t size() const
  {
    return m_data.size();
  }

  /** Return whether the vector is empty. */
  bool empty() const noexcept
  {
    return m_data.empty();
  }

  /** @} */
  /** \name Iterators */
  /** @{ */

  /** Mutable iterator at the start of the vector. */
  typename std::vector<value_type>::iterator begin()
  {
    return m_data.begin();
  }

  /** Mutable iterator one past the end of the vector. */
  typename std::vector<value_type>::iterator end()
  {
    return m_data.end();
  }

  /** Constant iterator at the start of the vector. */
  typename std::vector<value_type>::const_iterator cbegin() const
  {
    return m_data.cbegin();
  }

  /** Constant iterator at the end of the vector. */
  typename std::vector<value_type>::const_iterator cend() const
  {
    return m_data.cend();
  }

  /** @} */
public:
  /** \name Scalar Operations */
  /** @{ */

  /** Element-wise negation. */
  Vector operator-() const
  {
    Vector x(m_data);
    for (auto& elem : x)
      elem -= elem;
    return x;
  }

  /** Element-wise negation in-place. */
  Vector& operator-()
  {
    for (auto& elem : m_data)
      elem -= elem;
    return *this;
  }

  /** Element-wise multiplication by a scalar. */
  Vector operator*(const value_type value) const
  {
    Vector x(m_data);
    for (auto& elem : x)
      elem *= value;
    return x;
  }

  /** Element-wise multiplication by a scalar in-place. */
  Vector& operator*=(const value_type value)
  {
    for (auto& elem : m_data)
      elem *= value;
    return *this;
  }

  /** Element-wise division by a scalar. */
  Vector operator/(const value_type value) const
  {
    Assert(value != 0.0, "Zero division error.");

    Vector x(m_data);
    for (auto& elem : x)
      elem /= value;
    return x;
  }

  /** Element-wise division by a scalar in-place. */
  Vector& operator/=(const value_type value)
  {
    Assert(value != 0.0, "Zero division error.");

    for (auto& elem : m_data)
      elem /= value;
    return *this;
  }

  /** @} */
  /** \name Vector-Vector Operations */
  /** @{ */

  /** Element-wise addition of two vectors. */
  Vector operator+(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    Vector x(m_data);
    for (size_t i = 0; i < x.size(); ++i)
      x[i] += other[i];
    return x;
  }

  /** Element-wise addition of two vectors in-place. */
  Vector& operator+=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] += other[i];
    return *this;
  }

  /** Element-wise subtraction of two vectors. */
  Vector operator-(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    Vector x(m_data);
    for (size_t i = 0; i < x.size(); ++i)
      x[i] -= other[i];
    return x;
  }

  /** Element-wise subtraction of two vectors in-place. */
  Vector& operator-=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] -= other[i];
    return *this;
  }

  /** Element-wise multiplication of two vectors. */
  Vector operator*(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    Vector x(m_data);
    for (size_t i = 0; i < x.size(); ++i)
      x[i] *= other[i];
    return x;
  }

  /** Element-wise multiplication of two vectors in-place. */
  Vector& operator*=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] *= other[i];
    return *this;
  }

  /** Element-wise division of two vectors. */
  Vector operator/(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Assert(not other.has_zero_elements(), "Zero division error.");

    Vector x(m_data);
    for (size_t i = 0; i < x.size(); ++i)
      x[i] /= other[i];
    return x;
  }

  /** Element-wise division of two vectors in-place. */
  Vector& operator/=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Assert(not other.has_zero_elements(), "Zero division error.");

    for (size_t i = 0; i < m_data.size(); ++i)
      m_data[i] /= other[i];
    return *this;
  }

  /**
   * Return the dot product between this and another vector.
   * \f$ c = \vec{x} \cdot \vec{y} = \sum_i x_i y_i .\f$
   */
  value_type dot(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");

    value_type c = 0.0;
    for (size_t i = 0; i < m_data.size(); ++i)
      c += m_data[i] * other[i];
    return c;
  }

  /** @} */
  /** \name  Norms */
  /** @{ */

  /** Compute the \f$ \ell_\infty \f$-norm.
   *  \f$ ||\vec{x}||_{\ell_\infty} = \max_i |x_i| \f$ */
  value_type linf_norm() const
  {
    value_type norm = 0.0;
    for (const auto& elem : m_data)
      if (std::fabs(elem) > norm)
        norm = std::fabs(elem);
    return norm;
  }

  /** Compute the \f$ \ell_1 \f$-norm.
   *  \f$ ||\vec{x}||_{\ell_1} = \sum_i |x_i| \f$ */
  value_type l1_norm() const
  {
    value_type norm = 0.0;
    for (const auto& elem : m_data)
      norm += std::fabs(elem);
    return norm;
  }

  /** Compute the \f$ \ell_2 \f$-norm.
   * \f$ ||\vec{x}||_{\ell_2} = \sqrt{ \sum_i |x_i|^2 } \f$ */
  value_type l2_norm() const
  {
    value_type norm = 0.0;
    for (const auto& elem : m_data)
      norm += std::fabs(elem * elem);
    return std::sqrt(norm);
  }

  /** Compute the \f$ \ell_{\ell_p} \f$-norm.
   *  \f$ ||\vec{x}||_{\ell_p} = \left( \sum_i |x_i|^p \right)^{1/p} \f$ */
  value_type lp_norm(const value_type p) const
  {
    value_type norm = 0.0;
    for (const auto& elem : m_data)
      norm += std::pow(std::fabs(elem), p);
    return std::pow(norm, 1.0/p);
  }

  /** @} */
  /** \name Vector Operations */
  /** @{ */

  /** Return the unit-length direction vector.
   *  \f$ \hat{x} = \frac{\vec{x}}{||\vec{x}||_{\ell_2} \f$
   *  \note If the vector is uniformly zero, nothing is done. */
  Vector direction() const
  {
    Vector x(m_data);
    value_type norm = x.l2_norm();
    return (norm == 0.0)? x : x/norm;
  }

  /** Normalize this vector to unit length in-place.
   *  \f$ \hat{x} = \frac{\vec{x}}{||\vec{x}||_{\ell_2}} \f$
   *  \note If the vector is uniformly zero, nothing is done. */
  Vector& normalize()
  {
    value_type norm = l2_norm();
    return (norm == 0.0)? *this : operator/=(norm);
  }

  /** Element-wise absolute value. */
  Vector fabs() const
  {
    Vector x(m_data);
    for (auto& elem : x)
      elem = std::fabs(elem);
    return x;
  }

  /** Element-wise absolute value in-place. */
  Vector& fabs()
  {
    for (auto& elem : m_data)
      elem = std::fabs(elem);
    return *this;
  }

  /** @} */
public:
  /** /name Print Utilities */
  /** @{ */

  /** Return the vector as a string. */
  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < m_data.size() - 1; ++i)
      ss << std::setprecision(6) << m_data[i] << " ";
    ss << std::setprecision(6) << m_data.back() << "]\n";
    return ss.str();
  }

  /** Print the vector to `std::cout`. */
  void print() const
  {
    std::cout << to_string();
  }

  /** @} */
private:

  /** Return whether the vector has zero elements. */
  bool has_zero_elements() const
  {
    for (const auto& elem : m_data)
      if (elem == 0.0) return true;
    return false;
  }
};


/*-------------------- Inline Implementations --------------------*/


/** Element-wise multiplication by a scalar. */
template<typename value_type>
inline Vector<value_type> operator*(const value_type value,
                                    const Vector<value_type>& x)
{
  return x * value;
}


}
#endif //VECTOR_H

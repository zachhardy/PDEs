#ifndef VECTOR_H
#define VECTOR_H

#include "exceptions.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cinttypes>

#include <iostream>
#include <iomanip>
#include <sstream>


namespace math
{

/** A class representing a general vector. */
template<typename number>
struct Vector
{
public:
  using iterator = typename std::vector<number>::iterator;
  using const_iterator = typename std::vector<number>::const_iterator;

private:
  std::vector<number> m_data;

public:
  /** Default constructor. */
  Vector() = default;

  /** Construct a vector with \p n elements insert to \p value */
  explicit Vector(const uint64_t n, const number value = 0.0)
    : m_data(n, value)
  {}

  /** Copy constructor. */
  Vector(const Vector& other) : m_data(other.m_data) {}

  /** Move constructor. */
  Vector(Vector&& other) : m_data(std::move(other.m_data)) {}

  /** Copy construction using an STL vector. */
  Vector(const std::vector<number>& other) : m_data(other) {}

  /** Move construction using an STL vector. */
  Vector(std::vector<number>&& other) : m_data(std::move(other)) {}

  /** Construct using an initializer list. */
  Vector(std::initializer_list<number> list) : m_data(list) {}

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
  Vector& operator=(const std::vector<number>& other)
  {
    m_data = other;
    return *this;
  }

  /** Move assignment from an STL vector. */
  Vector& operator=(std::vector<number>&& other)
  {
    m_data = std::move(other);
    return *this;
  }

  /** Assigment using an initializer init_list. */
  Vector& operator=(std::initializer_list<number>& init_list)
  {
    m_data = init_list;
    return *this;
  }

public:
  //###########################################################################
  /** \name Access Operators */
  // @{

  /** Read/write access for element \p i. */
  double& operator[](const uint64_t i) { return m_data[i]; }

  /** Read only access for element \p i. */
  double operator[](const uint64_t i) const { return m_data[i]; }

  /** Read/write access for element \p i with bounds checking. */
  double& at(const uint64_t i) { return m_data.at(i); }

  /** Read only access for element \p i with bounds checking. */
  double at(const uint64_t i) const { return m_data.at(i); }

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

  // @}

  //###########################################################################
  /** \name Modifiers */
  // @{

  /** Clear the elements. */
  void clear() { m_data.clear(); }

  /** Insert a new element at the back of the vector. */
  void push_back(const number value) { m_data.push_back(value); }

  /** Insert a new element in place at the back of the vector. */
  void emplace_back(const number value) { m_data.emplace_back(value); }

  /** Remove the last element. */
  void pop_back() { m_data.pop_back(); }

  /** Resize to \p n elements, setting new elements to default. */
  void resize(const uint64_t n, const number value = 0.0)
  { m_data.resize(n, value); }

  /** Swap the elements of this vector and another. */
  void swap(Vector& other) { m_data.swap(other.m_data); }

  /** Swap the elements of this vector and an STL vector. */
  void swap(std::vector<number>& other) { m_data.swap(other); }

  // @}

  //###########################################################################
  /** \name Memory */
  // @{

  /** Allocate memory for \p n elements. */
  void reserve(const uint64_t n) { m_data.reserve(n); }

  /** Return the number of elements. */
  uint64_t size() const { return m_data.size(); }

  /** Return whether the vector is empty. */
  bool empty() const { return m_data.empty(); }

  // @}

  //###########################################################################
  /** \name Iterators */
  // @{

  /** Mutable iterator at the start of the vector. */
  iterator begin() { return m_data.begin(); }

  /** Mutable iterator one past the end of the vector. */
  iterator end() { return m_data.end(); }

  /** Constant iterator at the start of the vector. */
  const_iterator begin() const { return m_data.begin(); }

  /** Constant iterator at the end of the vector. */
  const_iterator end() const { return m_data.end(); }

  // @}

  //###########################################################################
  /** \name Scalar Operations */
  // @{

  /** Element-wise negation. */
  Vector operator-() const
  {
    Vector x(m_data);
    for (auto& elem : x) elem -= elem;
    return x;
  }

  /** Element-wise negation in-place. */
  Vector& operator-()
  {
    for (auto& elem : m_data) elem -= elem;
    return *this;
  }

  /** Element-wise multiplication by a scalar. */
  Vector operator*(const number value) const
  {
    Vector x(m_data);
    for (auto& elem : x) elem *= value;
    return x;
  }

  /** Element-wise multiplication by a scalar in-place. */
  Vector& operator*=(const number value)
  {
    for (auto& elem : m_data) elem *= value;
    return *this;
  }

  /** Element-wise division by a scalar. */
  Vector operator/(const number value) const
  {
    Assert(value != 0.0, "Zero division error.");
    Vector x(m_data);
    for (auto& elem : x) elem /= value;
    return x;
  }

  /** Element-wise division by a scalar in-place. */
  Vector& operator/=(const number value)
  {
    Assert(value != 0.0, "Zero division error.");
    for (auto& elem : m_data) elem /= value;
    return *this;
  }

  // @}

  //###########################################################################
  /** \name Vector-Vector Operations */
  // @{

  /** Element-wise addition of two vectors. */
  Vector operator+(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Vector x(m_data);
    for (uint64_t i = 0; i < x.size(); ++i)
      x[i] += other[i];
    return x;
  }

  /** Element-wise addition of two vectors in-place. */
  Vector& operator+=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    for (uint64_t i = 0; i < m_data.size(); ++i)
      m_data[i] += other[i];
    return *this;
  }

  /** Element-wise subtraction of two vectors. */
  Vector operator-(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Vector x(m_data);
    for (uint64_t i = 0; i < x.size(); ++i)
      x[i] -= other[i];
    return x;
  }

  /** Element-wise subtraction of two vectors in-place. */
  Vector& operator-=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    for (uint64_t i = 0; i < m_data.size(); ++i)
      m_data[i] -= other[i];
    return *this;
  }

  /** Element-wise multiplication of two vectors. */
  Vector operator*(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Vector x(m_data);
    for (uint64_t i = 0; i < x.size(); ++i)
      x[i] *= other[i];
    return x;
  }

  /** Element-wise multiplication of two vectors in-place. */
  Vector& operator*=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    for (uint64_t i = 0; i < m_data.size(); ++i)
      m_data[i] *= other[i];
    return *this;
  }

  /** Element-wise division of two vectors. */
  Vector operator/(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Assert(not other.has_zero_elements(), "Zero division error.");
    Vector x(m_data);
    for (uint64_t i = 0; i < x.size(); ++i)
      x[i] /= other[i];
    return x;
  }

  /** Element-wise division of two vectors in-place. */
  Vector& operator/=(const Vector& other)
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    Assert(not other.has_zero_elements(), "Zero division error.");
    for (uint64_t i = 0; i < m_data.size(); ++i)
      m_data[i] /= other[i];
    return *this;
  }

  /**
   * Return the dot product between this and another vector.
   * \f$ c = \vec{x} \cdot \vec{y} = \sum_i x_i y_i .\f$
   */
  number dot(const Vector& other) const
  {
    Assert(other.size() == m_data.size(), "Dimension mismatch error.");
    number c = 0.0;
    for (uint64_t i = 0; i < m_data.size(); ++i)
      c += m_data[i] * other[i];
    return c;
  }

  // @}

  //###########################################################################
  /** \name  Norms */
  // @{

  /** Compute the \f$ \ell_\infty \f$-norm.
   *  \f$ ||\vec{x}||_{\ell_\infty} = \max_i |x_i| \f$ */
  number linf_norm() const
  {
    number norm = 0.0;
    for (const auto& elem : m_data)
      if (std::fabs(elem) > norm)
        norm = std::fabs(elem);
    return norm;
  }

  /** Compute the \f$ \ell_1 \f$-norm.
   *  \f$ ||\vec{x}||_{\ell_1} = \sum_i |x_i| \f$ */
  number l1_norm() const
  {
    number norm = 0.0;
    for (const auto& elem : m_data)
      norm += std::fabs(elem);
    return norm;
  }

  /** Compute the \f$ \ell_2 \f$-norm.
   * \f$ ||\vec{x}||_{\ell_2} = \sqrt{ \sum_i |x_i|^2 } \f$ */
  number l2_norm() const
  {
    number norm = 0.0;
    for (const auto& elem : m_data)
      norm += std::fabs(elem * elem);
    return std::sqrt(norm);
  }

  /** Compute the \f$ \ell_{\ell_p} \f$-norm.
   *  \f$ ||\vec{x}||_{\ell_p} = \left( \sum_i |x_i|^p \right)^{1/p} \f$ */
  number lp_norm(const number p) const
  {
    number norm = 0.0;
    for (const auto& elem : m_data)
      norm += std::pow(std::fabs(elem), p);
    return std::pow(norm, 1.0/p);
  }

  // @}

  //###########################################################################
  /** \name Vector Operations */
  // @{

  /** Return the unit-length direction vector.
   *  \f$ \hat{x} = \frac{\vec{x}}{||\vec{x}||_{\ell_2} \f$
   *  \note If the vector is uniformly zero, nothing is done. */
  Vector direction() const
  {
    Vector x(m_data);
    number norm = x.l2_norm();
    return (norm == 0.0)? x : x/norm;
  }

  /** Normalize this vector to unit length in-place.
   *  \f$ \hat{x} = \frac{\vec{x}}{||\vec{x}||_{\ell_2}} \f$
   *  \note If the vector is uniformly zero, nothing is done. */
  Vector& normalize()
  {
    number norm = l2_norm();
    return (norm == 0.0)? *this : operator/=(norm);
  }

  /** Element-wise absolute value. */
  Vector fabs() const
  {
    Vector x(m_data);
    for (auto& elem : x) elem = std::fabs(elem);
    return x;
  }

  /** Element-wise absolute value in-place. */
  Vector& fabs()
  {
    for (auto& elem : m_data) elem = std::fabs(elem);
    return *this;
  }

  // @}

  //###########################################################################

  /** Print the vector. */
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

    os << "[";
    for (uint64_t i = 0; i < size() - 1; ++i)
      os << std::setw(w) << m_data[i];
    os << std::setw(w) << m_data.back() << "]" << std::endl;
  }

  // @}

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

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
  using value_type      = number;
  using size_type       = uint64_t;

  using pointer         = value_type*;
  using const_pointer   = const value_type*;
  using reference       = value_type&;
  using const_reference = const value_type&;

  using iterator        = value_type*;
  using const_iterator  = const value_type*;

private:
  std::vector<value_type> values;

public:
  /// Default constructor.
  Vector()
    : values(0)
  {}

  /// Copy constructor.
  Vector(const Vector& other)
    : values(other.values)
  {}

  /// Move constructor.
  Vector(Vector&& other)
    : values(std::move(other.values))
  {}

  /// Construct a Vector with \p n uninitialized elements.
  explicit Vector(const uint64_t n)
    : values(n)
  {}

  /// Construct a Vector with \p n elements set to \p value.
  explicit
  Vector(const uint64_t n, const value_type value)
    : values(n, value) {}

  /// Copy constructor using an STL vector.
  Vector(const std::vector<value_type>& other)
    : values(other) {}

  /// Move construction using an STL vector.
  Vector(std::vector<value_type>&& other)
    : values(other) {}

  /// Construct using an initializer list.
  Vector(std::initializer_list<value_type> list)
    : values(list) {}

  /// Copy assignment.
  Vector&
  operator=(const Vector& other)
  {
    values = other.values;
    return *this;
  }

  /// Move assignment.
  Vector&
  operator=(Vector&& other)
  {
    values = std::move(other.values);
    return *this;
  }

  /** Copy assignment from an STL vector. */
  Vector&
  operator=(const std::vector<value_type>& other)
  {
    values = other;
    return *this;
  }

  /// Move assignment from an STL vector.
  Vector&
  operator=(std::vector<value_type>&& other)
  {
    values = other;
    return *this;
  }

  /// Copy assigment using an initializer list.
  Vector&
  operator=(std::initializer_list<value_type> list)
  {
    values = list;
    return *this;
  }

  /// Assign a value to the vector.
  Vector&
  operator=(const value_type value)
  {
    for (auto& elem : values)
      elem = value;
    return *this;
  }

  /// Equality comparison operator.
  bool
  operator==(const Vector& other)
  { return values == other.values; }

  /// Inequality comparison operator.
  bool
  operator!=(const Vector& other)
  { return values != other.values; }

  /** \name Information */
  // @{

  /// Return the number of elements in the Vector.
  size_type
  size() const
  { return values.size(); }

  /// Return whether the Vector is empty.
  bool
  empty() const
  { return values.empty(); }

  /// Return whether the Vector is uniformly zero.
  bool
  all_zero() const
  {
    for (const auto& entry : values)
      if (entry != 0.0) return false;
    return true;
  }

  // @}
  /** \name Iterators */
  // @{

  /// Mutable iterator to the start of the Vector.
  iterator
  begin()
  { return values.begin(); }

  /// Constant iterator to the start of the Vector.
  const_iterator
  begin() const
  { return values.begin(); }

  /// Mutable iterator to the end of the Vector.
  iterator
  end()
  { return values.end(); }

  /// Constant iterator to the end of the Vector.
  const_iterator
  end() const
  { return values.end(); }

  // @}
  /** \name Modifiers */
  // @{

  /// Return the Vector to an uninitialized state.
  void
  clear()
  { values.clear(); }

  /// Insert a new element to the back of the Vector.
  void
  push_back(const value_type value)
  { values.push_back(value); }

  /// Move a new element to the back of the Vector.
  void
  push_back(value_type&& value)
  { values.emplace_back(value); }

  /// Remove the last element from the Vector.
  void
  pop_back()
  { values.pop_back(); }

  /// Resize the Vector to \p n elements. New elements remain uninitialized.
  void
  resize(const size_type n)
  { values.resize(n); }


  /// Resize the Vector to \p n elements, setting new elements to \p value.
  void
  resize(const size_type n, const value_type value)
  { values.resize(n, value); }

  /// Swap the elements of this Vector and another.
  void
  swap(Vector& other)
  { values.swap(other.values); }

  // @}
  /** \name Access Operators */
  // @{

  /// Read/write access for element \p i.
  reference
  operator[](const uint64_t i)
  { return values[i]; }

  /// Read only access for element \p i.
  const_reference
  operator[](const uint64_t i) const
  { return values[i]; }

  /// Read/write access for element \p i with bounds checking.
  reference
  at(const uint64_t i)
  { return values.at(i); }

  /// Read only access for element \p i with bounds checking.
  const_reference
  at(const uint64_t i) const
  { return values.at(i); }

  /// Read/write access for the first element.
  reference
  front()
  { return values.front(); }

  /// Read only access for the first element.
  const_reference
  front() const
  { return values.front(); }

  /// Read/write access for the last element.
  reference
  back()
  { return values.back(); }

  /// Read only access for the last element.
  const_reference
  back() const
  { return values.back(); }

  /// Read/write access to the underlying data.
  pointer
  data()
  { return values.data(); }

  /// Read only access to the underlying data.
  const_pointer
  data() const
  { return values.data(); }

  // @}
  /** \name Scalar Operations */
  // @{

  /// Element-wise negation in-place.
  Vector&
  operator-()
  {
    for (auto& elem : values)
      elem -= elem;
    return *this;
  }

  /// Element-wise negation.
  Vector
  operator-() const
  {
    Vector y = values;
    return -y;
  }

  /// Element-wise multiplication by a scalar in-place.
  Vector&
  operator*=(const value_type value)
  {
    for (auto& elem : values)
      elem *= value;
    return *this;
  }

  /// Element-wise multiplication by a scalar.
  Vector
  operator*(const value_type value) const
  {
    Vector y = values;
    y *= value;
    return y;
  }

  /// Element-wise division by a scalar in-place.
  Vector&
  operator/=(const value_type value)
  {
    Assert(value != 0.0, "Zero division error.");
    for (auto& elem : values)
      elem /= value;
    return *this;
  }

  /// Element-wise division by a scalar.
  Vector
  operator/(const value_type value) const
  {
    Vector y = values;
    y /= value;
    return y;
  }

  // @}
  /** \name Vector-Vector Operations */
  // @{

  /// Element-wise addition of two vectors in-place.
  Vector&
  operator+=(const Vector& y)
  {
    Assert(y.size() == values.size(), "Dimension mismatch error.");
    for (uint64_t i = 0; i < values.size(); ++i)
      values[i] += y[i];
    return *this;
  }

  /// Element-wise addition of two vectors.
  Vector
  operator+(const Vector& y) const
  {
    Vector z = values;
    z += y;
    return z;
  }

  /// Element-wise subtraction of two vectors in-place.
  Vector&
  operator-=(const Vector& y)
  {
    Assert(y.size() == values.size(), "Dimension mismatch error.");
    for (uint64_t i = 0; i < values.size(); ++i)
      values[i] -= y[i];
    return *this;
  }

  /// Element-wise subtraction of two vectors.
  Vector
  operator-(const Vector& y) const
  {
    Vector z = values;
    z -= y;
    return z;
  }

  /// Element-wise multiplication of two vectors in-place.
  Vector&
  operator*=(const Vector& y)
  {
    Assert(y.size() == values.size(), "Dimension mismatch error.");
    for (uint64_t i = 0; i < values.size(); ++i)
      values[i] *= y[i];
    return *this;
  }

  /// Element-wise multiplication of two vectors.
  Vector operator*(const Vector& y) const
  {
    Vector z = values;
    z *= y;
    return z;
  }

  /// Element-wise division of two vectors in-place.
  Vector&
  operator/=(const Vector& y)
  {
    Assert(y.size() == values.size(), "Dimension mismatch error.");
    Assert(not y.all_zero(), "Zero division error.");
    for (uint64_t i = 0; i < values.size(); ++i)
      values[i] /= y[i];
    return *this;
  }

  /// Element-wise division of two vectors. */
  Vector
  operator/(const Vector& y) const
  {
    Vector z = values;
    z /= y;
    return z;
  }

  /**
   * Return the dot product between this and another vector.
   * \f$ c = \vec{x} \cdot \vec{y} = \sum_i x_i y_i .\f$
   */
  value_type
  dot(const Vector& y) const
  {
    Assert(y.size() == values.size(), "Dimension mismatch error.");
    value_type c = 0.0;
    for (uint64_t i = 0; i < values.size(); ++i)
      c += values[i] * y[i];
    return c;
  }

  // @}
  /** \name  Norms */
  // @{

  /**
   * Compute the \f$ \ell_\infty \f$-norm.
   * \f$ ||\vec{x}||_{\ell_\infty} = \max_i |x_i| \f$
   */
  value_type
  linf_norm() const
  {
    value_type norm = 0.0;
    for (const auto& elem : values)
      if (std::fabs(elem) > norm)
        norm = std::fabs(elem);
    return norm;
  }

  /**
   * Compute the \f$ \ell_1 \f$-norm.
   * \f$ ||\vec{x}||_{\ell_1} = \sum_i |x_i| \f$
   */
  value_type
  l1_norm() const
  {
    value_type norm = 0.0;
    for (const auto& elem : values)
      norm += std::fabs(elem);
    return norm;
  }

  /**
   * Compute the \f$ \ell_2 \f$-norm.
   * \f$ ||\vec{x}||_{\ell_2} = \sqrt{ \sum_i |x_i|^2 } \f$
   */
  value_type
  l2_norm() const
  {
    value_type norm = 0.0;
    for (const auto& elem : values)
      norm += std::fabs(elem * elem);
    return std::sqrt(norm);
  }

  /**
   * Compute the \f$ \ell_{\ell_p} \f$-norm.
   * \f$ ||\vec{x}||_{\ell_p} = \left( \sum_i |x_i|^p \right)^{1/p} \f$
   */
  value_type
  lp_norm(const value_type p) const
  {
    value_type norm = 0.0;
    for (const auto& elem : values)
      norm += std::pow(std::fabs(elem), p);
    return std::pow(norm, 1.0/p);
  }

  // @}
  /** \name Vector Operations */
  // @{

  /**
   * Normalize this vector to unit length in-place.
   * \f$ \hat{x} = \frac{\vec{x}}{||\vec{x}||_{\ell_2}} \f$
   *
   * \note If the vector is uniformly zero, nothing is done.
   */
  Vector&
  normalize()
  {
    value_type norm = l2_norm();
    return (norm == 0.0)? *this : this->operator/=(norm);
  }

  /**
   * Return the unit-length direction vector.
   * \f$ \hat{x} = \frac{\vec{x}}{||\vec{x}||_{\ell_2} \f$
   *
   * \note If the vector is uniformly zero, nothing is done.
   */
  Vector
  direction() const
  { return Vector(values).normalize(); }

  /// Element-wise absolute value in-place.
  Vector&
  fabs()
  {
    for (auto& elem : values)
      elem = std::fabs(elem);
    return *this;
  }

  /// Element-wise absolute value.
  Vector
  fabs() const
  { return Vector<number>(values).fabs(); }

  // @}

  /// Print the vector.
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
      w = (!width)? precision + 7 : w;
    }
    else
    {
      os.setf(std::ios::fixed, std::ios::floatfield);
      w = (!width)? precision + 4 : w;
    }

    os << "[";
    for (uint64_t i = 0; i < size() - 1; ++i)
      os << std::setw(w) << values[i];
    os << std::setw(w) << values.back() << "]" << std::endl;
  }
};


/*-------------------- Inline Implementations --------------------*/


/** Element-wise multiplication by a scalar. */
template<typename number>
inline Vector<number>
operator*(const number value,
          const Vector<number>& x)
{
  return x * value;
}


}
#endif //VECTOR_H

#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <cinttypes>
#include <cmath>

#include <iostream>
#include <sstream>
#include <iomanip>

#include "macros.h"


namespace pdes::Math
{

class Vector
{
public:
  using value_type = double;

private:
  std::vector<value_type> elements;

public:
  /**
   * Default constructor.
   */
  Vector() = default;

  /**
   * Construct a vector with \p n uninitialized elements.
   */
  explicit
  Vector(const size_t n) : elements(n) {}

  /**
   * Construct a vector with \p n elements set to \p value.
   */
  explicit
  Vector(const size_t n, const value_type value) :
    elements(n, value)
  {}

  /**
   * Copy constructor from an STL vector.
   */
  Vector(const std::vector<value_type>& other) : elements(other) {}

  /**
   * Move constructor from an STL vector.
   */
  Vector(std::vector<value_type>&& other) :
    elements(other)
  {}

  /**
   * Construct from an initializer list.
   */
  Vector(const std::initializer_list<value_type> list) :
    elements(list)
  {}

  /**
   * Copy assignment from an STL vector.
   */
  Vector&
  operator=(const std::vector<value_type>& other)
  {
    elements = other;
    return *this;
  }

  /**
   * Move assignment from an STL vector.
   */
  Vector&
  operator=(std::vector<value_type>&& other)
  {
    elements = other;
    return *this;
  }

  /**
   * Copy assignment from an initializer list.
   */
  Vector&
  operator=(const std::initializer_list<value_type> list)
  {
    elements = list;
    return *this;
  }

  /**
   * Assign a value to all elements of the vector.
   */
  Vector&
  operator=(const value_type value)
  {
    for (auto& el : elements)
      el = value;
    return *this;
  }

  /**
   * Test the equality of two vectors.
   */
  bool
  operator==(const Vector& y) const
  { return (elements == y.elements); }

  /**
   * Test the inequality of two vectors.
   */
  bool
  operator!=(const Vector& y) const
  { return (elements != y.elements); }

  /** \name Characteristics */
  // @{

  /**
   * Return the number of elements in the vector.
   */
  size_t
  size() const
  { return elements.size(); }

  /**
   * Return whether the vector is empty.
   */
  bool
  empty() const
  { return elements.empty(); }

  /**
   * Return whether the vector is uniformly zero.
   */
  bool
  all_zero() const
  {
    for (const auto& el : elements)
      if (el != 0.0) return false;
    return true;
  }

  // @}
  /** \name Iterators */
  // @{

  /**
   * Mutable iterator to the first element of the vector.
   */
  std::vector<value_type>::iterator
  begin()
  { return elements.begin(); }

  /**
   * Mutable iterator to the end of the vector.
   */
  std::vector<value_type>::iterator
  end()
  { return elements.end(); }

  /**
   * Constant iterator to the start of the vector.
   */
  std::vector<value_type>::const_iterator
  begin() const
  { return elements.begin(); }

  /**
   * Constant iterator to the end of the vector.
   */
  std::vector<value_type>::const_iterator
  end() const
  { return elements.end(); }

  // @}
  /** \name Accessors */
  // @{

  /**
   * Read and write access for element \p i.
   */
  value_type&
  operator[](const size_t i)
  { return elements[i]; }

  /**
   * Read access for element \p i.
   */
  const value_type&
  operator[](const size_t i) const
  { return elements[i]; }

  /**
   * Read and write access for element \p i.
   */
  value_type&
  operator()(const size_t i)
  { return elements[i]; }

  /**
   * Read access for element \p i.
   */
  const value_type&
  operator()(const size_t i) const
  { return elements[i]; }

  /**
   * Read and write access for element \p i with bounds checking.
   */
  value_type&
  at(const size_t i)
  { return elements.at(i); }

  /**
   * Read access for element \p i with bounds checking
   */
  const value_type&
  at(const size_t i) const
  { return elements.at(i); }

  /**
   * Read and write access for the first element of the vector.
   */
  value_type&
  front()
  { return elements.front(); }

  /**
   * Read access for the first element of the vector.
   */
  const value_type&
  front() const
  { return elements.front(); }

  /**
   * Read and write access for the last element of the vector.
   */
  value_type&
  back()
  { return elements.back(); }

  /**
   * Read access for the last element of the vector.
   */
  const value_type&
  back() const
  { return elements.back(); }

  /**
   * Mutable pointer to the underlying data.
   */
  value_type*
  data()
  { return elements.data(); }

  /**
   * Constant pointer to the underlying data.
   */
  const value_type*
  data() const
  { return elements.data(); }

  // @}
  /** \name Modifiers */
  // @{

  /**
   * Return the vector to an uninitialized state.
   */
  void
  clear()
  { elements.clear(); }

  /**
   * Insert a new element at the back of the vector.
   */
  void
  push_back(const value_type value)
  { elements.push_back(value); }

  /**
   * Remove the last element from the vector.
   */
  void
  pop_back()
  { elements.pop_back(); }

  /**
   * Resize the vector to \p n elements. New elements remain uninitialized.
   */
  void
  resize(const size_t n)
  { elements.resize(n); }

  /**
   * Resize the vector to \p n elements with new elements set to \p value.
   */
  void
  resize(const size_t n, const value_type value)
  { elements.resize(n, value); }

  /**
   * Swap the elements with another vector.
   */
  void
  swap(Vector& y)
  { elements.swap(y.elements); }

  // @}
  /** \name Scalar Operations */
  // @{

  /**
   * Negate all elements of the vector. This is computed via
   * \f$ \vec{x} = -\vec{x} = -x_i, ~ \forall i \f$.
   */
  Vector&
  operator-()
  {
    for (auto& el : elements)
      el = -el;
    return *this;
  }

  /**
   * Return the negated vector.
   * \see Vector::operator-()
   */
  Vector
  operator-() const
  { return -Vector(elements); }

  /**
   * Multiply the elements of the vector by a scalar.
   * This is computed via
   * \f$ \vec{x} = \alpha \vec{x} = \alpha x_i, ~ \forall i \f$.
   */
  Vector&
  operator*=(const value_type factor)
  {
    for (auto& el : elements)
      el *= factor;
    return *this;
  }

  /**
   * Divide the elements of the vector by a scalar.
   * This is computed via
   * \f$ \vec{x} = \frac{\vec{x}}{\alpha}
   *             = \frac{x_i}{\alpha}, ~ \forall i
   * \f$.
   */
  Vector&
  operator/=(const value_type factor)
  {
    Assert(factor != 0.0, "Zero division error.");
    for (auto& el : elements)
      el /= factor;
    return *this;
  }

  // @}
  /** \name Linear Algebra Operations */
  // @{

  /**
   * Add another vector. This is computed via
   * \f$ \vec{x} = \vec{x} + \vec{y} = x_i + y_i, ~ \forall i \f$.
   */
  Vector&
  operator+=(const Vector& y)
  {
    Assert(size() == y.size(), "Dimension mismatch error.");
    for (size_t i = 0; i < size(); ++i)
      elements[i] += y.elements[i];
    return *this;
  }

  /**
   * Subtract another vector. This is computed via
   * \f$ \vec{x} = \vec{x} - \vec{y} = x_i - y_i, ~ \forall i \f$.
   */
  Vector&
  operator-=(const Vector& y)
  {
    Assert(size() == y.size(), "Dimension mismatch error.");
    for (size_t i = 0; i < size(); ++i)
      elements[i] -= y.elements[i];
    return *this;
  }

  /**
   * Return the dot product with another vector. This is computed via
   * \f$ \vec{x} \cdot \vec{y} = \sum_i x_i y_i \f$.
   */
  value_type
  dot(const Vector& y) const
  {
    Assert(size() == y.size(), "Dimension mismatch error.");
    double c = 0.0;
    for (size_t i = 0; i < size(); ++i)
      c += elements[i] * y.elements[i];
    return c;
  }

  // @}
  /** \name Vector-Norms */
  // @{

  /**
   * Compute the \f$ \ell_\infty \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_\infty} = \max_i |x_i| \f$.
   */
  value_type
  linf_norm() const
  {
    double norm = 0.0;
    for (const auto& el : elements)
      if (std::fabs(el) > norm)
        norm = std::fabs(el);
    return norm;
  }

  /**
   * Compute the \f$ \ell_1 \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_1} = \sum_i |x_i| \f$.
   */
  value_type
  l1_norm() const
  {
    double norm = 0.0;
    for (const auto& el : elements)
      norm += std::fabs(el);
    return norm;
  }

  /**
   * Compute the \f$ \ell_2 \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_2} = \sqrt{ \sum_i |x_i|^2 } \f$.
   */
  value_type
  l2_norm() const
  {
    double norm = 0.0;
    for (const auto& el : elements)
      norm += std::fabs(el) * std::fabs(el);
    return std::sqrt(norm);
  }

  /**
   * Compute the \f$ \ell_p \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_p} = \left( \sum_i |x_i|^p \right)^{1/p} \f$.
   */
  value_type
  lp_norm(const value_type p) const
  {
    double norm = 0.0;
    for (const auto& el : elements)
      norm += std::pow(std::fabs(el), p);
    return std::pow(norm, 1.0/p);
  }

  // @}
  /** \name Miscellaneous Operations */
  // @{

  /**
   * Normalize the vector to unit-length. This is computed via
   * \f$ \vec{x} = \hat{x} = \frac{\vec{x}}{|| \vec{x} ||_{\ell_2}} \f$.
   */
  Vector&
  normalize()
  {
    double norm = l2_norm();
    return *this /= (norm != 0.0) ? norm : 1.0;
  }

  /**
   * Return the unit-length vector.
   * \see Vector::normalize
   */
  Vector
  unit() const
  { return Vector(elements).normalize(); }

  /**
   * Take the absolute value of all the elements. This is computed via
   * \f$ \vec{x} = | \vec{x} | = |x_i|, ~ \forall i \f$.
   */
  Vector&
  fabs()
  {
    for (auto& el : elements)
      el = std::fabs(el);
    return *this;
  }

  /**
   * Return the absolute value of the vector.
   * \see Vector::fabs
   */
  Vector
  fabs() const
  { return Vector(elements).fabs(); }

  // @}
  /** \name Print Utilities */
  // @{

  /**
   * Return the vector as a string.
   *
   * \param scientific A flag for using scientific notation.
   * \param precision The precision to use when printing elements.
   * \param width The width between elements.
   *
   * \see Vector::print
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
   * Print the vector to an output stream.
   *
   * \param os The output stream to print the vector in.
   * \param scientific A flag for using scientific notation.
   * \param precision The precision to use when printing elements.
   * \param width The width between elements.
   *
   * \see Vector::str
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
      w = (!width) ? precision + 7 : w;
    }
    else
    {
      os.setf(std::ios::fixed, std::ios::floatfield);
      w = (!width) ? precision + 4 : w;
    }

    os << "[";
    for (size_t i = 0; i < size() - 1; ++i)
      os << std::setw(w) << elements[i];
    os << std::setw(w) << elements.back() << "]\n";

    os.flags(old_flags);
    os.precision(old_precision);
  }

  // @}

};

/*-------------------- Method Declarations --------------------*/

/**
 * Multiply each element of the vector by a scalar value.
 * \see Vector::operator*=
 */
inline Vector
operator*(const Vector& x, const double factor)
{ return Vector(x) *= factor; }

/**
 * Multiply each element of the vector by a scalar value.
 * \see Vector::operator*=
 */
inline Vector
operator*(const double factor, const Vector& x)
{ return Vector(x) *= factor; }

/**
 * Divide each element of the vector by a scalar value.
 * \see Vector::operator/=
 */
inline Vector
operator/(const Vector& x, const double factor)
{ return Vector(x) /= factor; }

/**
 * Add two vectors together.
 * \see Vector::operator+=
 */
inline Vector
operator+(const Vector& x, const Vector& y)
{ return Vector(x) += y; }

/**
 * Subtract two vectors.
 * \see Vector::operator-=
 */
inline Vector
operator-(const Vector& x, const Vector& y)
{ return Vector(x) -= y; }

/**
 * Return the dot product of two vectors.
 * \see Vector::dot
 */
inline double
dot(const Vector& x, const Vector& y)
{ return x.dot(y); }

/**
 * Return the absolute value of a vector.
 * \see Vector::fabs
 */
inline Vector
fabs(const Vector& x)
{ return x.fabs(); }

/**
 * Return the unit-length vector.
 * \see Vector::unit
 */
inline Vector
unit(const Vector& x)
{ return x.unit(); }

/**
 * Insert the vector into an output stream.
 * \see Vector::str Vector::print
 */
inline std::ostream&
operator<<(std::ostream& os, const Vector& x)
{ return os << x.str(); }

}
#endif //VECTOR_H

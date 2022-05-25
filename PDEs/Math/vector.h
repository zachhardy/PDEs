#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace pdes::Math
{

class Vector
{
public:
  using value_type = double;

protected:
  std::vector<value_type> elements;

public:

  //================================================== Constructors

  /**
   * Default constructor.
   */
  Vector() = default;

  /**
   * Construct a vector with \p n uninitialized elements.
   */
  explicit
  Vector(const size_t n);

  /**
   * Construct a vector with \p n elements set to \p value.
   */
  explicit
  Vector(const size_t n, const value_type value);

  /**
   * Copy constructor from an STL vector.
   */
  Vector(const std::vector<value_type>& other);

  /**
   * Move constructor from an STL vector.
   */
  Vector(std::vector<value_type>&& other);

  /**
   * Construct from an initializer list.
   */
  Vector(const std::initializer_list<value_type> list);

  //================================================== Assignment

  /**
   * Copy assignment from an STL vector.
   */
  Vector&
  operator=(const std::vector<value_type>& other);

  /**
   * Move assignment from an STL vector.
   */
  Vector&
  operator=(std::vector<value_type>&& other);

  /**
   * Copy assignment from an initializer list.
   */
  Vector&
  operator=(const std::initializer_list<value_type> list);

  /**
   * Assign a value to all elements of the vector.
   */
  Vector&
  operator=(const value_type value);

  //================================================== Comparison

  /**
   * Test the equality of two vectors.
   */
  bool
  operator==(const Vector& y) const;

  /**
   * Test the inequality of two vectors.
   */
  bool
  operator!=(const Vector& y) const;

  //================================================== Characteristics

  /** \name Characteristics */
  // @{

  /**
   * Return the number of elements in the vector.
   */
  size_t
  size() const;

  /**
   * Return whether the vector is empty.
   */
  bool
  empty() const;

  /**
   * Return whether the vector is uniformly zero.
   */
  bool
  all_zero() const;

  // @}

  //================================================== Iterators

  /** \name Iterators */
  // @{

  /**
   * Mutable iterator to the first element of the vector.
   */
  std::vector<value_type>::iterator
  begin();

  /**
   * Mutable iterator to the end of the vector.
   */
  std::vector<value_type>::iterator
  end();

  /**
   * Constant iterator to the start of the vector.
   */
  std::vector<value_type>::const_iterator
  begin() const;

  /**
   * Constant iterator to the end of the vector.
   */
  std::vector<value_type>::const_iterator
  end() const;

  // @}

  //================================================== Accessors

  /** \name Accessors */
  // @{

  /**
   * Read and write access for element \p i.
   */
  value_type&
  operator[](const size_t i);

  /**
   * Read access for element \p i.
   */
  const value_type&
  operator[](const size_t i) const;

  /**
   * Read and write access for element \p i.
   */
  value_type&
  operator()(const size_t i);

  /**
   * Read access for element \p i.
   */
  const value_type&
  operator()(const size_t i) const;

  /**
   * Read and write access for element \p i with bounds checking.
   */
  value_type&
  at(const size_t i);

  /**
   * Read access for element \p i with bounds checking.
   */
  const value_type&
  at(const size_t i) const;

  /**
   * Read and write access for the first element of the vector.
   */
  value_type&
  front();

  /**
   * Read access for the first element of the vector.
   */
  const value_type&
  front() const;

  /**
   * Read and write access for the last element of the vector.
   */
  value_type&
  back();

  /**
   * Read access for the last element of the vector.
   */
  const value_type&
  back() const;

  /**
   * Mutable pointer to the underlying data.
   */
  value_type*
  data();

  /**
   * Constant pointer to the underlying data.
   */
  const value_type*
  data() const;

  // @}

  //================================================== Modifiers

  /** \name Modifiers */
  // @{

  /**
   * Return the vector to an uninitialized state.
   */
  void
  clear();

  /**
   * Insert a new element at the back of the vector.
   */
  void
  push_back(const value_type value);

  /**
   * Remove the last element from the vector.
   */
  void
  pop_back();

  /**
   * Resize the vector to \p n elements. New elements remain uninitialized.
   */
  void
  resize(const size_t n);

  /**
   * Resize the vector to \p n elements with new elements set to \p value.
   */
  void
  resize(const size_t n, const value_type value);

  /**
   * Swap the elements with another vector.
   */
  void
  swap(Vector& y);

  /**
   * Normalize the vector to unit-length. This is computed via
   * \f$ \vec{x} = \hat{x} = \frac{\vec{x}}{|| \vec{x} ||_{\ell_2}} \f$.
   */
  Vector&
  normalize();

  /**
   * Return the unit-length vector.
   * \see Vector::normalize
   */
  Vector
  unit() const;

  /**
   * Take the absolute value of all the elements. This is computed via
   * \f$ \vec{x} = | \vec{x} | = |x_i|, ~ \forall i \f$.
   */
  Vector&
  fabs();

  /**
   * Return the absolute value of the vector.
   * \see Vector::fabs
   */
  Vector
  fabs() const;

  // @}

  //================================================== Scalar Operations

  /** \name Scalar Operations */
  // @{

  /**
   * Negate all elements of the vector. This is computed via
   * \f$ \vec{x} = -\vec{x} = -x_i, ~ \forall i \f$.
   */
  Vector&
  operator-();

  /**
   * Return the negated vector.
   * \see Vector::operator-()
   */
  Vector
  operator-() const;

  /**
   * Multiply the elements of the vector by a scalar.
   * This is computed via
   * \f$ \vec{x} = \alpha \vec{x} = \alpha x_i, ~ \forall i \f$.
   */
  Vector&
  operator*=(const value_type factor);

  /**
   * Divide the elements of the vector by a scalar.
   * This is computed via
   * \f$ \vec{x} = \frac{\vec{x}}{\alpha}
   *             = \frac{x_i}{\alpha}, ~ \forall i
   * \f$.
   */
  Vector&
  operator/=(const value_type factor);

  // @}

  //================================================== Linear Algebra

  /** \name Linear Algebra */
  // @{

  /**
   * Add another vector. This is computed via
   * \f$ \vec{x} = \vec{x} + \vec{y} = x_i + y_i, ~ \forall i \f$.
   */
  Vector&
  operator+=(const Vector& y);

  /**
   * Subtract another vector. This is computed via
   * \f$ \vec{x} = \vec{x} - \vec{y} = x_i - y_i, ~ \forall i \f$.
   */
  Vector&
  operator-=(const Vector& y);

  /**
   * Return the dot product with another vector. This is computed via
   * \f$ \vec{x} \cdot \vec{y} = \sum_i x_i y_i \f$.
   */
  value_type
  dot(const Vector& y) const;

  // @}

  //================================================== Vector Norms

  /** \name Vector Norms */
  // @{

  /**
   * Compute the \f$ \ell_\infty \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_\infty} = \max_i |x_i| \f$.
   */
  value_type
  linf_norm() const;

  /**
   * Compute the \f$ \ell_1 \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_1} = \sum_i |x_i| \f$.
   */
  value_type
  l1_norm() const;

  /**
   * Compute the \f$ \ell_2 \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_2} = \sqrt{ \sum_i |x_i|^2 } \f$.
   */
  value_type
  l2_norm() const;

  /**
   * Compute the \f$ \ell_p \f$-norm. This is computed via
   * \f$ || \vec{x} ||_{\ell_p} = \left( \sum_i |x_i|^p \right)^{1/p} \f$.
   */
  value_type
  lp_norm(const value_type p) const;

  // @}

  //================================================== Print Utilities

  /** \name Print Utilities */
  // @{

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
        const unsigned int width = 0) const;

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
      const unsigned int width = 0) const;

  // @}

};

//================================================== Methods

/**
 * Multiply each element of the vector by a scalar value.
 * \see Vector::operator*=
 */
Vector
operator*(const Vector& x, const double factor);

/**
 * Multiply each element of the vector by a scalar value.
 * \see Vector::operator*=
 */
Vector
operator*(const double factor, const Vector& x);

/**
 * Divide each element of the vector by a scalar value.
 * \see Vector::operator/=
 */
Vector
operator/(const Vector& x, const double factor);

/**
 * Add two vectors together.
 * \see Vector::operator+=
 */
Vector
operator+(const Vector& x, const Vector& y);

/**
 * Subtract two vectors.
 * \see Vector::operator-=
 */
Vector
operator-(const Vector& x, const Vector& y);

/**
 * Return the dot product of two vectors.
 * \see Vector::dot
 */
double
dot(const Vector& x, const Vector& y);

/**
 * Return the absolute value of a vector.
 * \see Vector::fabs
 */
Vector
fabs(const Vector& x);

/**
 * Return the unit-length vector.
 * \see Vector::unit
 */
Vector
unit(const Vector& x);

/**
 * Insert the vector into an output stream.
 * \see Vector::str Vector::print
 */
std::ostream&
operator<<(std::ostream& os, const Vector& x);

}
#endif //VECTOR_H

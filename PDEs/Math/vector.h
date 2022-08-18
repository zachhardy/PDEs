#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace Math
{
  /** Implementation of a general linear algebra vector. */
  class Vector
  {
  public:
    using iterator = std::vector<double>::iterator;
    using const_iterator = std::vector<double>::const_iterator;

  protected:
    std::vector<double> vals;

  public:

    //################################################## Initialization

    /** \name Construction and Initialization */
    // @{

    Vector() = default;

    /** Construct a vector with \p n uninitialized elements. */
    explicit Vector(const size_t n);

    /** Construct a vector with \p n elements set to \p value. */
    explicit Vector(const size_t n, const double value);

    Vector(const std::vector<double>& other);
    Vector(std::vector<double>&& other);

    Vector(const std::initializer_list<double> list);

    Vector& operator=(const std::vector<double>& other);
    Vector& operator=(std::vector<double>&& other);

    Vector& operator=(const std::initializer_list<double> list);

    /** Element-wise assignment to a scalar. */
    Vector& operator=(const double value);

    // @}

    //################################################## Information

    /** \name Information */
    // @{

    size_t size() const;

    /** Return the number of nonzero elements in the Vector. */
    size_t nnz() const;

    bool empty() const;
    bool all_zero() const;

    bool operator==(const Vector& other) const;
    bool operator!=(const Vector& other) const;

    // @}

    //################################################## Accessors

    /** \name Accessors */
    // @{

    double& operator[](const size_t i);
    const double& operator[](const size_t i) const;

    double& operator()(const size_t i);
    const double& operator()(const size_t i) const;

    double& at(const size_t i);
    const double& at(const size_t i) const;

    double& front();
    const double& front() const;

    double& back();
    const double& back() const;

    double* data();
    const double* data() const;

    // @}

    //################################################## Iterators

    /** \name Iterators */
    // @{

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    // @}

    //################################################## Modifiers

    /** \name Modifiers */
    // @{

    void clear();

    void push_back(const double value);

    void pop_back();

    /**
     * Resize the Vector to \p n elements. If \p n is less than the current
     * size, elements are deleted from the back. If \p n is greater than the
     * current size, all new elements are uninitialized.
     */
    void resize(const size_t n);

    /**
     * Resize the Vector to \p n elements. If new elements are created, they
     * are set to \p value.
     */
    void resize(const size_t n, const double value);

    void swap(Vector& y);

    // @}

    //################################################## Scalar Operations

    /** \name Scalar Operations and Norms */
    // @{

    /**
     * Take the dot product with another Vector via \f$ c = \vec{x} \cdot
     * \vec{y} = \sum_{i=0}^{N} x_i y_i \f$.
     */
    double dot(const Vector& y) const;

    /**
     * Return the \f$ \ell_\infty \f$-norm via \f$ || \vec{x} ||_{\infty} =
     * \max_i |x_i| \f$.
     */
    double linfty_norm() const;

    /**
     * Return the \f$ \ell_1 \f$-norm via \f$ || \vec{x} ||_{\ell_1} =
     * \sum_i |x_i| \f$.
     */
    double l1_norm() const;

    /**
     * Return the \f$ \ell_2 \f$-norm via \f$ || \vec{x} ||_{\ell_2} =
     * \sqrt{ \sum_i |x_i|^2 } \f$.
     */
    double l2_norm() const;

    /**
     * Return the \f$ \ell_p \f$-norm via \f$ || \vec{x} ||_{\ell_p} =
     * \left( \sum_i |x_i|^p \right)^{1/p} \f$.
     */
    double lp_norm(const double p) const;

    // @}

    //################################################## Linear Algebra

    /** \name Linear Algebra */
    // @{

    /**
     * Element-wise multiplication by a scalar in place.
     */
    Vector& scale(const double factor);

    /** Element-wise multiplication by a list of scalars in place. */
    Vector& scale(const Vector scaling_factors);

    /** Element-wise addition of a scalar value in place. */
    Vector& add(const double value);

    /** Element-wise addition with another Vector scaled by \p a in place. */
    Vector& add(const Vector& y, const double a);

    /**
     * Element-wise multiplication by a scalar followed by element-wise
     * addition with a Vector in place.
     */
    Vector& sadd(const double a, const Vector& y);

    /**
     * Element-wise multiplication by a scalar \p a followed by element-wise
     * addition with a Vector scaled by \p b in place.
     */
    Vector& sadd(const double a, const double b, const Vector& y);

    /** Element-wise assignment to a scaled Vector. */
    Vector& equal(const Vector& y, const double factor = 1.0);

    /** Element-wise absolute value in place. */
    Vector& fabs();

    /** Return a Vector containing the absolute value of the elements. */
    Vector fabs() const;

    /** Element-wise negation in place. */
    Vector& operator-();

    /** Return a Vector containing the negated elements. */
    Vector operator-() const;

    /** Element-wise multiplication by a scalar in place. */
    Vector& operator*=(const double factor);

    /** Element-wise division by a scalar in place. */
    Vector& operator/=(const double factor);

    /** Element-wise addition with another Vector. */
    Vector& operator+=(const Vector& y);

    /** Element-wise subtraction with another Vector. */
    Vector& operator-=(const Vector& y);

    // @}

    //################################################## Print Utilities

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
    void print(std::ostream& os = std::cout,
               const bool scientific = true,
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
    std::string str(const bool scientific = true,
                    const unsigned int precision = 3,
                    const unsigned int width = 0) const;

    // @}
  };

  /** Element-wise multiplication by a scalar. */
  Vector operator*(const Vector& x, const double factor);

  /** Element-wise multiplication by a scalar. */
  Vector operator*(const double factor, const Vector& x);

  /** Element-wise division by a scalar. */
  Vector operator/(const Vector& x, const double factor);

  /** Element-wise addition. */
  Vector operator+(const Vector& x, const Vector& y);

  /** Element-wise subtraction. */
  Vector operator-(const Vector& x, const Vector& y);

  /** Compute a dot product. \see Vector::dot */
  double dot(const Vector& x, const Vector& y);

  /** Return the absolute value the elements of a Vector. */
  Vector fabs(const Vector& x);

  /**
   * Return the \f$ \ell_{\infty} \f$-norm of a vector.
   * \see Vector::linfty_norm
   */
  double linfty_norm(const Vector& x);

  /** Return the \f$\ell_1\f$-norm of a vector. \see Vector::l1_norm */
  double l1_norm(const Vector& x);

  /**
   * Return the \f$\ell_2\f$-norm of a vector. \see Vector::l2_norm */
  double l2_norm(const Vector& x);

  /** Return the \f$\ell_p\f$-norm of a vector. \see Vector::lp_norm */
  double lp_norm(const Vector& x, const double p);

  std::ostream& operator<<(std::ostream& os, const Vector& x);
}
#endif //VECTOR_H

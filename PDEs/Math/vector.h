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
    using value_type = double;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

  protected:
    std::vector<value_type> vals;

  public:

    //################################################## Initialization

    /** \name Construction and Initialization */
    // @{

    Vector() = default;

    /* Construct a vector with \p n uninitialized elements. */
    explicit Vector(const size_t n);

    /** Construct a vector with \p n elements set to \p value. */
    explicit Vector(const size_t n, const value_type value);

    /** Copy constructor from an STL vector. */
    Vector(const std::vector<value_type>& other);

    /** Move constructor from an STL vector. */
    Vector(std::vector<value_type>&& other);

    /** Construct from an initializer list. */
    Vector(const std::initializer_list<value_type> list);

    /** Copy assignment from an STL vector. */
    Vector& operator=(const std::vector<value_type>& other);

    /** Move assignment from an STL vector. */
    Vector& operator=(std::vector<value_type>&& other);

    /** Copy assignment from an initializer list. */
    Vector& operator=(const std::initializer_list<value_type> list);

    /** Assign a value to all elements of the vector. */
    Vector& operator=(const value_type value);

    // @}

    //################################################## Information

    /** \name Information */
    // @{

    size_t size() const;

    /** Return the number of nonzero elements in the vector. */
    size_t nnz() const;

    bool empty() const;
    bool all_zero() const;

    bool operator==(const Vector& other) const;
    bool operator!=(const Vector& other) const;

    // @}

    //################################################## Accessors

    /** \name Accessors */
    // @{

    value_type& operator[](const size_t i);
    const value_type& operator[](const size_t i) const;

    value_type& operator()(const size_t i);
    const value_type& operator()(const size_t i) const;

    value_type& at(const size_t i);
    const value_type& at(const size_t i) const;

    value_type& front();
    const value_type& front() const;

    value_type& back();
    const value_type& back() const;

    value_type* data();
    const value_type* data() const;

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

    void push_back(const value_type value);

    void pop_back();

    void resize(const size_t n);
    void resize(const size_t n, const value_type value);

    void swap(Vector& y);

    // @}

    //################################################## Scalar Operations

    /** \name Scalar Operations and Norms */
    // @{

    /**
     * Return the dot product with another vector. The sum of element-wise
     * products, given by \f$ \vec{x} \cdot \vec{y} = \sum_i x_i y_i \f$.
     */
    value_type dot(const Vector& y) const;

    /**
     * Return the \f$ \ell_\infty \f$-norm. The maximum absolute value, given by
     * \f$ || \vec{x} ||_{\ell_\infty} = \max_i |x_i| \f$.
     */
    value_type linfty_norm() const;

    /**
     * Compute the \f$ \ell_1 \f$-norm. The sum of absolute values, given by
     * \f$ || \vec{x|| ||_{\ell_1} = \sum_i |x_i| \f$.
     */
    value_type l1_norm() const;

    /**
     * Compute the \f$ \ell_2 \f$-norm. The square-root of the sum of squares,
     * given by \f$ || \vec{x} ||_{\ell_2} = \sqrt{ \sum_i |x_i|^2 } \f$.
     */
    value_type l2_norm() const;

    /**
     * Compute the \f$ \ell_p \f$-norm. The <tt>p</tt>'th root of the sum
     * of the <tt>p</tt>'th power of the absolute value of the elements, given by
     * \f$ || \vec{x} ||_{\ell_p} = \left( \sum_i |x_i|^p \right)^{1/p} \f$.
     */
    value_type lp_norm(const value_type p) const;

    // @}

    //################################################## Linear Algebra

    /** \name Linear Algebra */
    // @{

    /**
     * Scale the vector by a scalar factor. This is computed via \f$ \vec{x}
     * = \alpha \vec{x} = \alpha x_i, ~ \forall i \f$.
     */
    Vector& scale(const value_type factor);

    /**
     * Scale by the specified scaling factors. This is computed via \f$ \vec{x}
     * = f_i x_i, \forall i \f$.
     */
    Vector& scale(const Vector scaling_factors);

    /** Add a scalar value to each element of the vector. */
    Vector& add(const value_type value);

    /**
     * Add a multiple of a vector. This is computed via \f$ \vec{x} = \vec{x} +
     * \alpha \vec{y} = x_i + \alpha y_i, ~ \forall i \f$. The default behavior
     * is \f$ \alpha = 1.0 \f$.
     */
    Vector& add(const Vector& y, const value_type a);

    /**
     * Multiply by a scalar value and add another vector. This is computed via
     * \f$ \vec{x} = \alpha \vec{x} + \vec{y} = \alpha x_i + y_i, ~ \forall i \f$.
     */
    Vector& sadd(const value_type a, const Vector& y);

    /**
     * Scale the vector and add another scaled vector to it. This is computed via
     * \f$ \vec{x} = \alpha \vec{x} + \beta \vec{y} = \alpha x_i + y_i, ~
     * \forall i \f$.
     */
    Vector& sadd(const value_type a, const value_type b, const Vector& y);

    /**
     * Set the vector to a multiple of another. This is computed via \f$ \vec{x} =
     * \alpha \vec{y} = \alpha y_i \f$.
     */
    Vector& equal(const Vector& y, const value_type factor = 1.0);

    /**
     * Take the absolute value of each element of the vector. This is computed
     * via \f$ \vec{x} = | \vec{x} | = |x_i|, ~ \forall i \f$.
     */
    Vector& fabs();

    /**
     * Negate all elements of the vector. This is computed via \f$ \vec{x} =
     * -\vec{x} = -x_i, ~ \forall i \f$.
     */
    Vector& operator-();

    /**
     * Return a vector with negated elements.
     *
     * \see Vector::operator-()
     * */
    Vector operator-() const;

    /**
     * Multiply the elements of the vector by a scalar. This is computed via
     * \f$ \vec{x} = \alpha \vec{x} = \alpha x_i, ~ \forall i \f$.
     */
    Vector& operator*=(const value_type factor);

    /**
     * Divide the elements of the vector by a scalar. This is computed via \f$
     * \vec{x} = \frac{\vec{x}}{\alpha} = \frac{x_i}{\alpha}, ~ \forall i \f$.
     */
    Vector& operator/=(const value_type factor);

    /**
     * Add another vector. This is computed via \f$ \vec{x} = \vec{x} + \vec{y}
     * = x_i + y_i, ~ \forall i \f$.
     */
    Vector& operator+=(const Vector& y);

    /**
     * Subtract another vector. This is computed via \f$ \vec{x} = \vec{x} -
     * \vec{y} = x_i - y_i, ~ \forall i \f$.
     */
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

  /** Multiply each element of the vector by a scalar value. */
  Vector operator*(const Vector& x, const double factor);

  /** Multiply each element of the vector by a scalar value. */
  Vector operator*(const double factor, const Vector& x);

  /** Divide each element of the vector by a scalar value. */
  Vector operator/(const Vector& x, const double factor);

  /** Add two vectors together. */
  Vector operator+(const Vector& x, const Vector& y);

  /** Subtract two vectors. */
  Vector operator-(const Vector& x, const Vector& y);

  /** Return the dot product of two vectors. */
  double dot(const Vector& x, const Vector& y);

  /** Return the absolute value of a vector. */
  Vector fabs(const Vector& x);

  /** Return the \f$\ell_{\infty}\f$-norm of a vector. */
  double linfty_norm(const Vector& x);

  /** Return the \f$\ell_1\f$-norm of a vector. */
  double l1_norm(const Vector& x);

  /** Return the \f$\ell_2\f$-norm of a vector. */
  double l2_norm(const Vector& x);

  /** Return the \f$\ell_p\f$-norm of a vector. */
  double lp_norm(const Vector& x, const double p);

  /** Insert the vector into an output stream. */
  std::ostream& operator<<(std::ostream& os, const Vector& x);
}
#endif //VECTOR_H

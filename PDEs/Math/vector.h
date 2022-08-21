#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace Math
{
  /**
   * Implementation of a general linear algebra vector.
   */
  class Vector
  {
  public:
    /**
     * Alias to an iterator over an STL vector of doubles.
     */
    using iterator = std::vector<double>::iterator;

    /*
     * Alias to a constant iterator over an STL vector of doubles.
     */
    using const_iterator = std::vector<double>::const_iterator;

  public:

    /**
     * \name Constructors and assignment
     */
    /* @{ */

    /**
     * Default constructor.
     */
    Vector() = default;

    /**
     * Construct a vector with \p n elements.
     */
    explicit Vector(const size_t n);

    /**
     * Construct a vector with \p n elements set to \p value.
     */
    explicit Vector(const size_t n, const double value);

    /**
     * Copy construction from an STL vector.
     */
    Vector(const std::vector<double>& other);

    /**
     * Move construction from an STL vector.
     * \param other
     */
    Vector(std::vector<double>&& other);

    /**
     * Copy construction from an initializer list.
     */
    Vector(const std::initializer_list<double> list);

    /**
     * Copy assignment with an STL vector.
     */
    Vector&
    operator=(const std::vector<double>& other);

    /**
     * Move assignment with an STL vector.
     */
    Vector&
    operator=(std::vector<double>&& other);

    /**
     * Copy assignment with an initializer list.
     */
    Vector&
    operator=(const std::initializer_list<double> list);

    /**
     * Element-wise assignment to a scalar.
     */
    Vector& operator=(const double value);

    /* @} */
    /**
     * \name Characteristics
     */
    // @{

    /**
     * Return the number of elements in the vector.
     */
    size_t
    size() const;

    /**
     * Return he number of non-zero elements in the vector.
     */
    size_t
    n_nonzero_elements() const;

    /**
     * Return whether the vector is empty (no allocated elements) or not.
     */
    bool
    empty() const;

    /* @} */
    /**
     * \name Comparison
     */
    /* @{ */

    /**
     * Return whether all elements of two vectors are equivalent.
     */
    bool
    operator==(const Vector& other) const;

    /**
     * Return whether any elements of two vectors are different.
     */
    bool
    operator!=(const Vector& other) const;

    /* @} */
    /**
     * \name Accessors
     */
    /* @{ */

    /**
     * Read and write access for element \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    double&
    operator[](const size_t i);

    /**
     * Read access for element \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    const double&
    operator[](const size_t i) const;

    /**
     * Read and write access for element \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    double&
    operator()(const size_t i);

    /**
     * Read access for element \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    const double&
    operator()(const size_t i) const;

    /**
     * Read and write access for element \p i with bounds checking.
     */
    double&
    at(const size_t i);

    /**
     * Read access for element \p i with bounds checking.
     */
    const double&
    at(const size_t i) const;

    /**
     * Read and write access to the first element.
     */
    double&
    front();

    /**
     * Read access to the first element.
     */
    const double&
    front() const;

    /**
     * Read and write access to the last element.
     */
    double&
    back();

    /**
     * Read access to the last elements.
     */
    const double&
    back() const;

    /**
     * Return a pointer to the underlying vector data.
     */
    double*
    data();

    /**
     * Return a constant pointer to the underlying vector data.
     */
    const double*
    data() const;

    /* @} */
    /**
     * \name Iterators
     */
    /* @{ */

    /**
     * Return an iterator to the start of the vector.
     */
    iterator
    begin();

    /**
     * Return an iterator that designates the end of the vector.
     */
    iterator
    end();

    /**
     * Return a constant iterator to the start of the vector.
     */
    const_iterator
    begin() const;

    /**
     * Return a constant iterator that designates the end of the vector.
     */
    const_iterator
    end() const;

    /* @} */
    /**
     * \name Modifiers
     */
    /* @{ */

    /**
     * Delete the contents of the vector.
     */
    void
    clear();

    /**
     * Add an element set to \p value to the back of the vector. This
     * increases the size of the vector.
     */
    void
    push_back(const double value);

    /**
     * Remove the last element from the vector. This decreases the size of
     * the vector.
     */
    void
    pop_back();

    /**
     * Resize the vector to \p n elements. If \p n is less than the current
     * number of elements, elements are deleted from the back. If \p n is
     * greater than the current size, all new elements are uninitialized.
     */
    void
    resize(const size_t n);

    /**
     * Resize the vector to \p n elements. If new elements are created, they
     * are set to \p value.
     */
    void
    resize(const size_t n, const double value);

    /**
     * Swap the contents of two vectors.
     */
    void
    swap(Vector& other);

    /**
     * Set this vector to \f$ \vec{y} \f$ scaled by \p a.
     */
    Vector&
    equal(const Vector& y, const double a = 1.0);

    /**
     * Element-wise absolute value in place.
     */
    Vector&
    fabs();

    /**
     * Return the absolute value of this vector.
     */
    Vector
    fabs() const;


    /* @} */
    /**
     * \name Scalar operations and norms
     */
    /* @{ */

    /**
     * Take the dot product with another Vector via \f$ c = \vec{x} \cdot
     * \vec{y} = \sum_i x_i y_i \f$.
     */
    double
    dot(const Vector& y) const;

    /**
     * Return the \f$ \ell_\infty \f$-norm via \f$ ||\vec{x}||_{\infty} =
     * \max_i |x_i| \f$.
     */
    double
    linfty_norm() const;

    /**
     * Return the \f$ \ell_1 \f$-norm via \f$ ||\vec{x}||_{\ell_1} =
     * \sum_i |x_i| \f$.
     */
    double
    l1_norm() const;

    /**
     * Return the \f$ \ell_2 \f$-norm via \f$ ||\vec{x}||_{\ell_2} =
     * \sqrt{ \sum_i |x_i|^2 } \f$.
     */
    double
    l2_norm() const;

    /**
     * Return the \f$ \ell_p \f$-norm via \f$ ||\vec{x}||_{\ell_p} =
     * \left( \sum_i |x_i|^p \right)^{1/p} \f$.
     */
    double
    lp_norm(const double p) const;

    /* @} */
    /**
     * \name Scaling operations
     */
    /* @{ */

    /**
     * Element-wise multiplication by a scalar in place.
     */
    Vector&
    scale(const double factor);

    /**
     * Element-wise multiplication by a list of scalars in place.
     */
    Vector&
    scale(const Vector scaling_factors);

    /**
     * Element-wise negation in place.
     */
    Vector&
    operator-();

    /**
     * Return the negative of a vector.
     */
    Vector
    operator-() const;

    /**
     * Element-wise multiplication by a scalar in place.
     */
    Vector&
    operator*=(const double factor);

    /**
     * Return a vector multiplied by a scalar.
     */
    Vector
    operator*(const double factor) const;

    /**
     * Element-wise division by a non-zero scalar in place.
     */
    Vector&
    operator/=(const double factor);

    /**
     * Return a vector divided by a non-zero scalar.
     */
    Vector
    operator/(const double factor) const;

    /* @} */
    /**
     * \name Addition and subtraction operations
     */
    /* @{ */

    /**
     * Shift each element in a vector by \p such that \f$ x_i = x_i + a \f$.
     */
    Vector&
    shift(const double value);

    /**
     * Element-wise addition by vector \f$ a \vec{y} \f$ in place.
     * This is computed via \f$ \vec{x} = \vec{x} + a \vec{y} = \sum_i x_i +
     * a y_i \f$.
     */
    Vector&
    add(const Vector& y, const double a = 1.0);

    /**
     * Element-wise addition with another vector in place.
     */
    Vector&
    operator+=(const Vector& y);

    /**
     * Return the sum of two vector.
     */
    Vector
    operator+(const Vector& y) const;

    /**
     * Element-wise subtraction by another vector in place.
     */
    Vector&
    operator-=(const Vector& y);

    /**
     * Return the difference between two vectors.
     */
    Vector
    operator-(const Vector& y) const;

    /**
     * Scale this vector by \p a and add vector \f$ \vec{y} \f$ in place.
     */
    Vector&
    sadd(const double a, const Vector& y);

    /**
     * Scale this vector by \p a and add vector \f$ \vec{y} \f$ scaled by \p b
     * in place.
     */
    Vector&
    sadd(const double a, const double b, const Vector& y);

    /* @} */
    /**
     * \name Print utilities
     */
    /* @{ */

    /**
     * Print the vector to an output stream.
     *
     * \param os The output stream to print the vector in.
     * \param scientific A flag for using scientific notation.
     * \param precision The precision to use when printing elements.
     * \param width The width between elements.
     */
    void
    print(std::ostream& os = std::cout,
          const bool scientific = true,
          const unsigned int precision = 3,
          const unsigned int width = 0) const;

    /**
     * Return the vector as a string. See \ref print.
     *
     * \param scientific A flag for using scientific notation.
     * \param precision The precision to use when printing elements.
     * \param width The width between elements.
     */
    std::string
    str(const bool scientific = true,
        const unsigned int precision = 3,
        const unsigned int width = 0) const;

    /* @} */

  protected:
    /**
     * The underlying vector data.
     */
    std::vector<double> values;
  };


  /**
   * Multiply a vector by a scalar value.
   */
  Vector
  operator*(const double factor, const Vector& x);

  /**
   * Take the dot product between two vectors. See \ref Vector::dot.
   */
  double
  dot(const Vector& x, const Vector& y);

  /**
   * Return the absolute value of a vector.
   */
  Vector
  fabs(const Vector& x);

  /**
   * Return the \f$ \ell_{\infty} \f$-norm of a vector. See \ref
   * Vector::linfty_norm.
   */
  double
  linfty_norm(const Vector& x);

  /**
   * Return the \f$ \ell_1 \f$-norm of a vector. See \ref Vector::l1_norm.
   */
  double
  l1_norm(const Vector& x);

  /**
   * Return the \f$ \ell_2 \f$-norm of a vector. See \ref Vector::l2_norm.
   */
  double
  l2_norm(const Vector& x);

  /**
   * Return the \f$\ell_p\f$-norm of a vector. See \ref Vector::lp_norm.
   */
  double
  lp_norm(const Vector& x, const double p);

  /**
   * Insert a vector into an output stream.
   */
  std::ostream&
  operator<<(std::ostream& os, const Vector& x);
}
#endif //VECTOR_H

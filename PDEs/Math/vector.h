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
     * Default constructor. Create an empty vector.
     */
    Vector() = default;

    /**
     * Copy constructor. Copy the internal data from another vector.
     */
    Vector(const Vector& other) = default;

    /**
     * Move constructor. Steal the internal data from another vector.
     */
    Vector(Vector&& other) = default;

    /**
     * Construct a vector with \p n elements, optionally set to \p value.
     */
    Vector(const size_t n, const double value = 0.0);

    /**
     * Copy construction from an initializer list.
     */
    Vector(const std::initializer_list<double>& list);

    /**
     * Constructor from iterators.
     */
    template<typename InputIterator>
    Vector(const InputIterator first, const InputIterator last);

    /**
     * Construct a vector from \p n contiguously stored elements.
     */
    Vector(const size_t n, const double* value_ptr);

    /**
     * Copy assignment from another vector.
     */
    Vector&
    operator=(const Vector& other);

    /**
     * Move assignment from another vector.
     */
    Vector&
    operator=(Vector&& other);

    /**
     * Assign each element of the vector to the specified \p value. If the
     * vector is empty, this will initialize a single element set to \p value.
     */
    Vector&
    operator=(const double value);

    /**
     * Reinitialize the vector to have \p n elements set to \p value. This
     * first clears the vector, then resizes it.
     */
    void
    reinit(const size_t n, const double value = 0.0);

    /* @} */
    /**
     * \name Information about the vector
     */
    // @{

    /**
     * Return the number of elements.
     */
    size_t
    size() const;

    /**
     * Return he number of non-zero elements.
     */
    size_t
    n_nonzero_elements() const;

    /**
     * Return whether the vector is empty (no allocated elements) or not.
     */
    bool
    empty() const;

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
     * \name Accessors and iterators
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
     * \name Modifying the vector
     */
    /* @{ */

    /**
     * Delete the contents of the vector.
     */
    void
    clear();

    /**
     * Add an element set to \p value to the back of the vector.
     */
    void
    push_back(const double value);

    /**
     * Remove the last element from the vector.
     */
    void
    pop_back();

    /**
     * Resize the vector to \p n elements. If \p n is less than the current
     * number of elements, elements are deleted from the back. If \p n is
     * greater than the current size, all new elements are uninitialized.
     */
    void
    resize(const size_t n, const double value = 0.0);

    /**
     * Swap the contents of two vectors.
     */
    void
    swap(Vector& other);

    /**
     * Set the vector to \f$ \vec{x} = a \vec{y} \f$
     */
    Vector&
    equal(const Vector& y, const double factor = 1.0);

    /**
     * Take the absolute value of each element of the vector in place.
     */
    Vector&
    fabs();

    /**
     * Return th absolute value of the vector.
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
     * Multiply each element of the vector by a scalar such that \f$ \vec{x} =
     * a \vec{x} \f$.
     */
    Vector&
    scale(const double factor);

    /**
     * Multiply each element of the vector by the corresponding element in the
     * argument such that \f$ x_i = a_i x_i \f$ where \f$ \vec{a} = (a_i, ...,
     * a_n) \f$ is of the same length as \f$ \vec{x} \f$.
     */
    Vector&
    scale(const Vector& scaling_factors);

    /**
     * Negate each element of the vector. This is equivalent to scaling by -1.0.
     * See \ref scale.
     */
    Vector&
    operator-();

    /**
     * Return a vector containing the negated elements of this vector. See
     * \ref scale.
     */
    Vector
    operator-() const;

    /**
     * Multiply each element of the vector by the argument. See \ref scale.
     */
    Vector&
    operator*=(const double factor);

    /**
     * Return a vector that contains the elements of this vector multiplied by
     * a scalar. See \ref scale.
     */
    Vector
    operator*(const double factor) const;

    /**
     * Divide the elements of the vector by a non-zero scalar. See \ref scale.
     */
    Vector&
    operator/=(const double factor);

    /**
     * Return a vector that contains the elements of this vector divided by a
     * non-zero scalar. See \ref scale.
     */
    Vector
    operator/(const double factor) const;

    /* @} */
    /**
     * \name Addition and subtraction operations
     */
    /* @{ */

    /**
     * Shift each element in a vector by \p such that \f$ x_i = x_i + a, ~
     * \forall i \f$.
     */
    Vector&
    shift(const double value);

    /**
     * Add a scaled vector to this one such that \f$ \vec{x} = \vec{x} + a
     * \vec{y} \f$.
     */
    Vector&
    add(const double a, const Vector& y);

    /**
     * Add another vector to this one. This is equivalent to adding a vector
     * scaled by 1.0. See \ref add
     */
    Vector&
    operator+=(const Vector& y);

    /**
     * Return the sum of this vector and another. See \ref add.
     */
    Vector
    operator+(const Vector& y) const;

    /**
     * Subtract another vector from this one. This is equivalent to adding a
     * vector scaled by -1.0. See \ref add.
     */
    Vector&
    operator-=(const Vector& y);

    /**
     * Return the difference between this vector and another. See \ref add.
     */
    Vector
    operator-(const Vector& y) const;

    /**
     * Multiply this vector by a scalar and add another scaled vector to it such
     * that \f$ \vec{x} = a \vec{x} + b \vec{y} \f$.
     */
    Vector&
    sadd(const double a, const double b, const Vector& y);

    /**
     * Scale this vector by a scalar and add another vector in place such
     * that \f$ \vec{x} = a \vec{x} + \vec{y} \f$. This is equivalent to
     * scaling the other vector by 1.0 using the more general \ref sadd routine.
     * See \ref sadd.
     */
    Vector&
    sadd(const double a, const Vector& y);

    /* @} */
    /**
     * \name Print utilities
     */
    /* @{ */

    /**
     * Return the vector as a string.
     *
     * \param scientific A flag for using scientific notation.
     * \param precision The precision to use when printing elements.
     * \param width The width between elements.
     */
    std::string
    str(const bool scientific = true,
        const unsigned int precision = 3,
        const unsigned int width = 0) const;

    /**
     * Print the vector to an output stream. See \ref str.
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

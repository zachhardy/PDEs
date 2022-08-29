#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace PDEs
{
  namespace Math
  {
    /**
     * Implementation of a general linear algebra vector.
     */
    class Vector
    {
    public:
      using iterator = std::vector<double>::iterator;
      using const_iterator = std::vector<double>::const_iterator;

    public:
      //################################################## Constructors

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
       * Construct a vector with \p n elements.
       */
      explicit Vector(const size_t n);

      /**
       * Construct a vector with \p n elements set to \p value.
       */
      Vector(const size_t n, const double value);

      /**
       * Construct a vector from an initializer list.
       */
      Vector(const std::initializer_list<double>& list);

      /**
       * Construct with the values pointed to by iterators in the range
       * <tt>[first, last)</tt>.
       */
      template<typename InputIterator>
      Vector(const InputIterator first, const InputIterator last);

      /**
       * Construct a vector from \p n contiguously stored entries.
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
       * Assign each element to the specified \p value.
       */
      Vector&
      operator=(const double value);

      //################################################## Capacity

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
       * Return whether the vector is empty (no allocated entries) or not.
       */
      bool
      empty() const;

      //################################################## Data Access

      /**
       * Read and write access for element \p i/
       */
      double&
      operator[](const size_t i);

      /**
       * Read access for element \p i.
       */
      const double&
      operator[](const size_t i) const;

      /**
       * Read and write access for element \p i.
       */
      double&
      operator()(const size_t i);

      /**
       * Read access for element \p i.
       */
      const double&
      operator()(const size_t i) const;

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

      //################################################## Modifiers

      /**
       * Delete the contents of the vector.
       */
      void
      clear();

      /**
       * Resize the vector to \p n elements.
       *
       * If \p n is less than the current number of elements, elements are
       * deleted from the back. If \p n is greater than the current size, new
       * elements are allocated and left in an unspecified state.
       */
      void
      resize(const size_t n);

      /**
       * Resize the vector to \p n elements and set any new elements to the
       * specified \p value.
       *
       * If \p n is less than the current number of elements, elements are
       * deleted from the back. If \p n is greater than the current size, new
       * elements are allocated and set to \p value.
       */
      void
      resize(const size_t n, const double value);

      /**
       * Swap the contents of two vectors.
       */
      void
      swap(Vector& other);

      /**
       * Set the vector to \f$ x = a y \f$
       */
      Vector&
      equal(const Vector& y, const double factor = 1.0);

      /**
       * Element-wise absolute value of the vector.
       */
      Vector&
      fabs();

      /**
       * Return a vector with the absolute value of the vector.
       */
      Vector
      fabs() const;

      //################################################## Scalar Operations
      //                                                   and Norms

      /**
       * Return the dot product with another vector.
       *
       * This is computed via
       * \f[ c = x \cdot y = \sum_i x_i y_i \forall i. \f]
       *
       * This operation is only successful if the vectors are of the same size.
       */
      double
      dot(const Vector& y) const;

      /**
       * Return the \f$ \ell_\infty \f$-norm.
       */
      double
      linfty_norm() const;

      /**
       * Return the \f$ \ell_1 \f$-norm.
       */
      double
      l1_norm() const;

      /**
       * Return the \f$ \ell_2 \f$-norm.
       */
      double
      l2_norm() const;

      /**
       * Return the \f$ \ell_p \f$-norm.
       */
      double
      lp_norm(const double p) const;

      //################################################## Linear Algebra

      /**
       * Element-wise multiplication by a scalar.
       */
      Vector&
      scale(const double factor);

      /**
       * Multiply each entry of the vector by the corresponding entry in the
       * argument.
       *
       * This operation is only successful if the two vectors are of the
       * same size.
       */
      Vector&
      scale(const Vector& scaling_factors);

      /**
       * Element-wise negation.
       */
      Vector&
      operator-();

      /**
       * Return a vector with the negated elements.
       */
      Vector
      operator-() const;

      /**
       * Element-wise multiplication by a scalar.
       */
      Vector&
      operator*=(const double factor);

      /**
       * Return a vector with the elements multiplied by a scalar.
       */
      Vector
      operator*(const double factor) const;

      /**
       * Element-wise division by a non-zero scalar.
       */
      Vector&
      operator/=(const double factor);

      /**
       * Return a vector with the elements divided by a non-zero scalar.
       */
      Vector
      operator/(const double factor) const;

      //################################################## Vector Operations

      /**
       * Element-wise addition by a scalar.
       */
      Vector&
      shift(const double value);

      /**
       * Element-wise multiplication by a scalar and addition of another scaled
       * vector, i.e. \f$ x = a x + b y \f$.
       */
      Vector&
      sadd(const double a, const double b, const Vector& y);

      /**
       * Element-wise multiplication by a scalar and addition of another, i.e.
       * \f$ x = a x + y \f$.
       */
      Vector&
      sadd(const double a, const Vector& y);

      /**
       * Element-wise addition of a scaled vector, i.e. \f$ x = x + b y \f$.
       */
      Vector&
      add(const double b, const Vector& y);

      /**
       * Element-wise addition of another vector, i.e. \f$ x = x + y \f$.
       */
      Vector&
      operator+=(const Vector& y);

      /**
       * Return the sum of two vectors.
       */
      Vector
      operator+(const Vector& y) const;

      /**
       * Element-wise subtraction by another vector, i.e. \f$ x = x - y \f$.
       */
      Vector&
      operator-=(const Vector& y);

      /**
       * Return the difference between two vectors.
       */
      Vector
      operator-(const Vector& y) const;

      //################################################## Print Utilities

      /**
       * Return the vector as a string with the specified formatting.
       */
      std::string
      str(const bool scientific = true,
          const unsigned int precision = 3,
          const unsigned int width = 0) const;

      /**
       * Print the vector to an output stream with the specified formatting.
       */
      void
      print(std::ostream& os = std::cout,
            const bool scientific = true,
            const unsigned int precision = 3,
            const unsigned int width = 0) const;

      //################################################## Comparison

      /**
       * Return whether all entries of two vectors are equivalent.
       */
      bool
      operator==(const Vector& other) const;

      /**
       * Return whether any entries of two vectors are different.
       */
      bool
      operator!=(const Vector& other) const;

      friend Vector
      operator*(const double factor, const Vector& x);

      friend std::ostream&
      operator<<(std::ostream& os, const Vector& x);

    protected:
      std::vector<double> values;
    };


    /**
     * Multiply a vector by a scalar value.
     */
    Vector
    operator*(const double factor, const Vector& x);

    /**
     * Take the dot product between two vectors.
     */
    double
    dot(const Vector& x, const Vector& y);

    /**
     * Return the absolute value of a vector.
     */
    Vector
    fabs(const Vector& x);

    /**
     * Return the \f$ \ell_{\infty} \f$-norm of a vector.
     */
    double
    linfty_norm(const Vector& x);

    /**
     * Return the \f$ \ell_1 \f$-norm of a vector.
     */
    double
    l1_norm(const Vector& x);

    /**
     * Return the \f$ \ell_2 \f$-norm of a vector.
     */
    double
    l2_norm(const Vector& x);

    /**
     * Return the \f$ \ell_p \f$-norm of a vector.
     */
    double
    lp_norm(const Vector& x, const double p);

    /**
     * Insert a vector into an output stream.
     */
    std::ostream&
    operator<<(std::ostream& os, const Vector& x);
  }
}

#endif //VECTOR_H

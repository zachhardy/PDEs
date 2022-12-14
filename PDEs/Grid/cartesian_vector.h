#ifndef CARTESIAN_VECTOR_H
#define CARTESIAN_VECTOR_H

#include <iostream>
#include <array>

namespace PDEs
{
  namespace Grid
  {
    // useful aliases
    class CartesianVector;

    using Point = CartesianVector;
    using Vertex = CartesianVector;
    using Node = CartesianVector;
    using Centroid = CartesianVector;
    using Normal = CartesianVector;
    using Gradient = CartesianVector;

    /**
     * A class representing a Cartesian vector in space.
     *
     * This class is meant to be generic and all-encompassing for representing
     * vertices on a mesh, nodes on a discretization, centroids on cells,
     * and normal vectors on faces, and gradient vectors.
     */
    class CartesianVector
    {
    protected:
      std::array<double, 3> xyz;

    public:
      /// \name Constructors and assignment
      /// @{

      /// Construct a zero vector <tt>(0, 0, 0)</tt>.
      CartesianVector();

      /// Construct the vector <tt>(a, 0, 0)</tt>
      explicit CartesianVector(const double a);

      /// Construct the vector <tt>(a, b, 0)</tt>.
      CartesianVector(const double a, const double b);

      /// Construct the vector <tt>(a, b, c)</tt>.
      CartesianVector(const double a, const double b, const double c);

      /// Copy constructor.
      CartesianVector(const CartesianVector&) = default;

      /// Move constructor.
      CartesianVector(CartesianVector&&) = default;

      /// Assignment to a scalar value.
      CartesianVector& operator=(const double value);

      /// Copy assignment.
      CartesianVector& operator=(const CartesianVector&) = default;

      /// Move assignment.
      CartesianVector& operator=(CartesianVector&&) = default;

      /**
       * Construct a unit vector in the specified dimension.
       *
       * The specified \p axis must be 0, 1, or 2 for the x, y, or
       * z directions.
       */
      static CartesianVector unit_vector(const unsigned int axis);

      /// @}
      /// \name Comparison operator
      /// @{

      /// Return whether two vectors are the same.
      bool operator==(const CartesianVector& other) const;

      /// Return whether two vectors are different.
      bool operator!=(const CartesianVector& other) const;

      /// @}
      /// Accessors
      /// @{

      /// Read access to element \p i.
      const double& operator[](const unsigned int i) const;

      /// Read/write access to element \p i.
      double& operator[](const unsigned int i);

      /// Read access to element \p i.
      const double& operator()(const unsigned int i) const;

      /// Read/write access to element \p i.
      double& operator()(const unsigned int i);

      /// Read access to the x-coordinate.
      const double& x() const;

      /// Read/write access to the x-coordinate.
      double& x();

      /// Read access to the y-coordinate
      const double& y() const;

      /// Read/write access to the y-coordinate.
      double& y();

      /// Read access to the z-coordinate
      const double& z() const;

      /// Read/write access to the z-coordinate.
      double& z();

      /// Return a pointer to the underlying data.
      double* data();

      /// Return a constant pointer to the underlying data.
      const double* data() const;

      /// @}
      /// \name Addition and subtractions
      /// @{

      /// In-place addition by another vector.
      CartesianVector& operator+=(const CartesianVector& other);

      /// Return the sum of two vectors.
      CartesianVector operator+(const CartesianVector& other) const;

      /// In-place subtraction by another vector.
      CartesianVector& operator-=(const CartesianVector& other);

      /// Return the difference between two vectors.
      CartesianVector operator-(const CartesianVector& other) const;

      /// @}
      /// \name Multiplication and division
      /// @{

      /// Scale a vector by a scalar value.
      CartesianVector& scale(const double factor);

      /// Return a vector scaled by a scalar value.
      CartesianVector scale(const double factor) const;

      /// Negate the elements of the vector.
      CartesianVector& operator-();

      /// Return the negative of the vector.
      CartesianVector operator-() const;

      /// In-place scalar multiplication.
      CartesianVector& operator*=(const double factor);

      /// Return the vector multiplied by a scalar value.
      CartesianVector operator*(const double factor) const;

      /// In-place scalar division.
      CartesianVector& operator/=(const double factor);

      /// Return the vector divided by a scalar value.
      CartesianVector operator/(const double factor) const;

      /// @}
      /// \name Other operations
      /// @{

      /// In-place element-wise absolute value.
      CartesianVector& fabs();

      /// Return the absolute value of a vector.
      CartesianVector fabs() const;

      /// Return the dot product between two vectors.
      double dot(const CartesianVector& other) const;

      /// Return the cross product between two vectors.
      CartesianVector cross(const CartesianVector& other) const;

      /// Return the length of a vector.
      double length() const;

      /// Return the squared length of a vector.
      double length_squared() const;

      /// Normalize this vector to its length.
      CartesianVector& normalize();

      /// Return the direction of a vector.
      CartesianVector direction() const;

      /// Return the distance between two vectors.
      double distance(const CartesianVector& other) const;

      /// Return the squared distance between two vectors.
      double distance_squared(const CartesianVector& other) const;

      /// @}
      /// \name Print utilities
      /// @{

      /// Return a vector as a string.
      std::string str() const;

      /// Print a vector to the standard output.
      void print() const;

      /// @}

      /// Left-multiply a vector by a scalar.
      friend CartesianVector
      operator*(const double factor, const CartesianVector& p);

      /// Insert a vector into an output stream.
      friend std::ostream&
      operator<<(std::ostream& os, const CartesianVector& p);
    };


    /// Return the dot product between two vectors.
    double dot(const CartesianVector& p,
               const CartesianVector& q);

    /// Return the cross product between two vectors.
    CartesianVector cross(const CartesianVector& p,
                          const CartesianVector& q);

    /// Return the length of a vector.
    double length(const CartesianVector& p);

    /// Return the squared length of a vector.
    double length_squared(const CartesianVector& p);

    /// Return the Euclidean distance between two vectors.
    double distance(const CartesianVector& p,
                    const CartesianVector& q);

    /// Return the squared Euclidean distance between two vectors.
    double distance_squared(const CartesianVector& p,
                            const CartesianVector& q);

    /// Return the absolute value of a vector.
    CartesianVector fabs(const CartesianVector& p);
  }
}

#endif //CARTESIAN_VECTOR_H

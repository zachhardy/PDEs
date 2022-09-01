#ifndef CARTESIAN_VECTOR_H
#define CARTESIAN_VECTOR_H

#include <iostream>
#include <array>


namespace PDEs
{
  namespace Grid
  {
    // Useful Aliases
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
      //######################################## Constructors

      /** Construct the point <tt>(0, 0, 0)</tt>. */
      CartesianVector();

      /** Construct the point <tt>(a, 0, 0)</tt>. */
      explicit CartesianVector(const double a);

      /** Construct the point <tt>(a, b, 0)</tt>. */
      CartesianVector(const double a, const double b);

      /** Construct the point <tt>(a, b, c)</tt>. */
      CartesianVector(const double a, const double b, const double c);

      /**
       * Construct a unit vector in the specified dimension.
       *
       * The specified \p axis must be 0, 1, or 2 for the x, y, or
       * z directions.
       */
      static CartesianVector unit_vector(const unsigned int axis);

      /** Entry-wise assignment to a scalar value. */
      CartesianVector& operator=(const double value);

      //################################################## Data Access

      /** Read and write access to entry \p i. */
      double& operator[](const unsigned int i);

      /** Read access to entry \p i. */
      const double& operator[](const unsigned int i) const;

      /** Read and write access to entry \p i. */
      double& operator()(const unsigned int i);

      /** Read access to entry \p i. */
      const double& operator()(const unsigned int i) const;

      /** Read and write access to the x-coordinate. */
      double& x();

      /** Read access to the x-coordinate. */
      const double& x() const;

      /** Read and write access to the y-coordinate. */
      double& y();

      /** Read access to the y-coordinate. */
      const double& y() const;

      /** Read and write access to the z-coordinate. */
      double& z();

      /** Read access to the z-coordinate. */
      const double& z() const;

      //################################################## Scalar Operations

      /** Return the Euclidean distance to the origin. */
      double length() const;

      /** Return the Euclidean distance to the origin squared. */
      double length_squared() const;

      /** Return the Euclidean distance between two Cartesian vectors. */
      double distance(const CartesianVector& other) const;

      /**
       * Return the squared Euclidean distance between two points Cartesian
       * vectors.
       */
      double distance_squared(const CartesianVector& other) const;

      /**
       * Take the dot product between two Cartesian vectors.
       *
       * This is computed via
       * \f[ c = p \cdot q = p_x q_x + p_y q_y + p_z q_z. \f]
       */
      double dot(const CartesianVector& other) const;

      //################################################## Cartesian Vector
      //                                                   Operations

      /** Entry-wise multiplication by a scalar. */
      CartesianVector& operator*=(const double factor);

      /**
       * Return a Cartesian vector with the elements multiplied by a scalar.
       */
      CartesianVector operator*(const double factor) const;

      /** Entry-wise division by a non-zero scalar. */
      CartesianVector&
      operator/=(const double factor);

      /**
       * Return a Cartesian vector with the elements divided by a non-zero
       * scalar.
       */
      CartesianVector operator/(const double factor) const;

      /** Entry-wise negation. */
      CartesianVector& operator-();

      /** Return a Cartesian vector with the negated elements. */
      CartesianVector operator-() const;

      /** Entry-wise addition of another Cartesian vector. */
      CartesianVector& operator+=(const CartesianVector& other);

      /** Return the element-wise sum of two Cartesian vectors. */
      CartesianVector operator+(const CartesianVector& other) const;

      /** Entry-wise subtraction by another Cartesian vector. */
      CartesianVector& operator-=(const CartesianVector& other);

      /** Return the element-wise difference of two Cartesian vectors. */
      CartesianVector operator-(const CartesianVector& other) const;

      /** Entry-wise absolute value. */
      CartesianVector& fabs();

      /** Return a Cartesian vector with the element-wise absolute value. */
      CartesianVector fabs() const;

      /**
       * Normalize the Cartesian vector to its length.
       *
       * This performs and element-wise division by the computed length of the
       * Cartesian vector. This results in a transformation to a unit, or
       * direction vector.
       */
      CartesianVector& normalize();

      /** Return a copy of the normalized Cartesian vector. */
      CartesianVector direction() const;

      /**
       * Take the cross product between two Cartesian vectors.
       *
       * This is computed via
       * \f[ r = p \times q
       *       = (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x).
       * \f]
       */
      CartesianVector cross(const CartesianVector& other) const;

      //################################################## Print Utilities

      /** Return a Cartesian vector as a string. */
      std::string str() const;

      /** Print a Cartesian vector into an output stream. */
      void print(std::ostream& os = std::cout) const;

      //################################################## Comparisons

      /**
       * Return whether all elements of two Cartesian vectors are equivalent.
       */
      bool operator==(const CartesianVector& other) const;

      /**
       * Return whether any elements of two Cartesian vectors are different.
       */
      bool operator!=(const CartesianVector& other) const;

      //################################################## Friends

      friend CartesianVector
      operator*(const double factor, const CartesianVector& p);

      friend std::ostream&
      operator<<(std::ostream& os, const CartesianVector& p);
    };

    /** Multiply a Cartesian vector by a scalar. */
    CartesianVector
    operator*(const double factor, const CartesianVector& p);

    /** Take the dot product between two Cartesian vectors. */
    double dot(const CartesianVector& p, const CartesianVector& q);

    /** Take the cross product between two Cartesian vectors. */
    CartesianVector cross(const CartesianVector& p, const CartesianVector& q);

    /** Return the Euclidean distance between two Cartesian vectors. */
    double distance(const CartesianVector& p, const CartesianVector& q);

    /** Return the absolute value of a Cartesian vector. */
    CartesianVector fabs(const CartesianVector& p);

    /** Return the direction of a Cartesian vector. */
    CartesianVector direction(const CartesianVector& p);

    /** Insert a Cartesian vector into an output stream. */
    std::ostream& operator<<(std::ostream& os, const CartesianVector& p);
  }
}
#endif //CARTESIAN_VECTOR_H

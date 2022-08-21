#ifndef CARTESIAN_VECTOR_H
#define CARTESIAN_VECTOR_H

#include <iostream>
#include <array>


namespace Grid
{

  /**
   * A class representing a Cartesian vector in space.
   *
   * This class is meant to be generic and all-encompassing for representing
   * vertices on a mesh, nodes on a discretization, centroids on cells,
   * and normal vectors on faces, and gradient vectors.
   */
  class CartesianVector
  {
  public:
    /**
     * \name Constructors and assignment
     */
    /* @{ */

    /**
     * Construct the point <tt>(0, 0, 0)</tt>.
     */
    CartesianVector();

    /**
     * Construct the point <tt>(a, 0, 0)</tt>.
     */
    explicit CartesianVector(const double a);

    /**
     * Construct the point <tt>(b, 0, 0)</tt>.
     */
    explicit CartesianVector(const double a, const double b);

    /**
     * Construct the point <tt>(a, b, c)</tt>.
     */
    explicit CartesianVector(const double a, const double b, const double c);

    /**
     * Construct a unit vector in the specified dimension.
     *
     * \param axis The axis to construct the unit vector for. 0 corresponds to
     *    the x-direction, 1 to the y-direction, and 2 to the z-direction.
     */
    static CartesianVector
    unit_vector(const unsigned int axis);

    /**
     * Element-wise assignment to a scalar value.
     */
    CartesianVector&
    operator=(const double value);

    /* @} */
    /**
     * \name Accessors
     */
    /* @} */

    /**
     * Read and write access to element \p i with bounds checking.
     */
    double&
    operator[](const unsigned int i);

    /**
     * Read access to element \p i with bounds checking.
     */
    const double&
    operator[](const unsigned int i) const;

    /**
     * Read and write access to element \p i with bounds checking.
     */
     double&
     operator()(const unsigned int i);

    /**
    * Read access to element \p i with bounds checking.
    */
    const double&
    operator()(const unsigned int i) const;

    /**
     * Return the x-coordinate.
     */
    const double&
    x() const;

    /**
     * Return the y-coordinate.
     */
    const double&
    y() const;

    /**
     * Return the z-coordinate.
     */
    const double&
    z() const;


    /* @} */
    /**
     * \name Comparisons
     */
    /* @{ */

    /**
     * Return whether all elements of two points are equivalent.
     */
    bool
    operator==(const CartesianVector& other) const;

    /**
     * Return whether any elements of two points are different.
     */
    bool
    operator!=(const CartesianVector& other) const;

    /**
     * Return the Euclidean distance between two points. This is given by \f$
     * d = \sqrt{(p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2} \f$.
     */
    double
    distance(const CartesianVector& other) const;

    /**
     * Return the Euclidean distance between two points squared. This is given
     * by \f$ d^2 = (p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2 \f$.
     */
    double
    distance_squared(const CartesianVector& other) const;

    /* @} */
    /**
     * \name Characteristics
     */
    /* @{ */

    /**
     * Return the Euclidean distance to the origin via \f$ \ell = \sqrt{ x^2 +
     * y^2 + z^2} \f$. See \ref length_squared.
     */
    double
    length() const;

    /**
     * Return the Euclidean distance to the origin squared. via \f$ \ell^2 =
     * x^2 + y^2 + z^2 \f$. See \ref length.
     */
    double
    length_squared() const;

    /* @} */
    /**
     * \name Arithmetic with points.
     */
    /* @{ */

    /**
     * Element-wise multiplication by a scalar.
     */
    CartesianVector&
    operator*=(const double factor);

    /**
     * Return a point multiplied by a scalar.
     */
    CartesianVector
    operator*(const double factor) const;

    /**
     * Element-wise division by a non-zero scalar.
     */
    CartesianVector&
    operator/=(const double factor);

    /**
     * Return a point divided by a non-zero scalar.
     */
     CartesianVector
     operator/(const double factor) const;

    /**
      * Element-wise negation in place.
      */
    CartesianVector&
    operator-();

    /**
     * Return the negative of a point.
     */
    CartesianVector
    operator-() const;

    /**
     * Element-wise addition in place.
     */
    CartesianVector&
    operator+=(const CartesianVector& other);

    /**
     * Return the sum of two points.
     */
    CartesianVector
    operator+(const CartesianVector& other) const;

    /**
     * Element-wise subtraction in place.
     */
    CartesianVector&
    operator-=(const CartesianVector& other);

    /**
     * Return the difference of two points.
     */
    CartesianVector
    operator-(const CartesianVector& other) const;

    /**
     * Element-wise absolute value in place.
     */
    CartesianVector&
    fabs();

    /**
     * Return the absolute value of a point.
     */
    CartesianVector
    fabs() const;

    /**
     * Normalize the point to its length. If the length is zero (e.g. the point
     * is the origin), return a copy. See \ref length.
     */
    CartesianVector&
    normalize();

    /**
     * Return the direction of a point. See \ref normalize.
     */
    CartesianVector
    direction() const;

    /**
     * Take the dot product between two points via \f$ c = p \cdot q = p_x q_x
     * + p_y q_y + p_z q_z \f$.
     */
    double
    dot(const CartesianVector& other) const;

    /**
     * Take the cross product between two points via \f$ r = p \times q =
     * (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x) \f$.
     */
    CartesianVector
    cross(const CartesianVector& other) const;

    /* @} */
    /**
     * \name Print utilities
     */
    /* @{ */

    /**
     * Return a point as a string.
     */
    std::string
    str() const;

    /**
     * Print a point into an output stream.
     */
    void
    print(std::ostream& os = std::cout) const;

    /**
     * Insert a point into an output stream.
     */

    /* @} */

  protected:
    /**
     * The x, y, z coordinates.
     */
    std::array<double, 3> values;
  };

  /**
   * Multiply a point by a scalar.
   */
  CartesianVector
  operator*(const double factor, const CartesianVector& p);

  /**
   * Take the dot product between two points. See \ref CartesianVector::dot.
   */
  double
  dot(const CartesianVector& p, const CartesianVector& q);

  /**
   * Take the cross product between two points. See \ref CartesianVector::cross
   */
  CartesianVector
  cross(const CartesianVector& p, const CartesianVector& q);

  /**
   * Return the Euclidean distance between two points. See \ref CartesianVector::distance.
   */
  double
  distance(const CartesianVector& p, const CartesianVector& q);

  /**
   * Return the absolute value of a point. See \ref CartesianVector::fabs.
   */
  CartesianVector
  fabs(const CartesianVector& p);

  /**
   * Return the direction of a vector pointing to the specified point.
   * See \ref CartesianVector::directions
   */
  CartesianVector
  direction(const CartesianVector& p);

  /**
   * Insert a point into an output stream.
   */
  std::ostream&
  operator<<(std::ostream& os, const CartesianVector& p);

  //================================================== Useful Aliases

  using Point = CartesianVector;
  using Vertex = CartesianVector;
  using Node = CartesianVector;
  using Centroid = CartesianVector;
  using Normal = CartesianVector;
  using Gradient = CartesianVector;

}
#endif //CARTESIAN_VECTOR_H

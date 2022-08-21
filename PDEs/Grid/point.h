#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <array>


namespace Grid
{

  /**
   * A class representing a Point in Cartesian space.
   *
   * This is intended to be a generic three-vector which can be used to
   * represent a Point, a Vertex on a Mesh, a Node in a Discretization,
   * a Normal on a Face, or a gradient.
   */
  class Point
  {
  public:
    /**
     * \name Constructors and assignment
     */
    /* @{ */

    /**
     * Construct the point <tt>(0, 0, 0)</tt>.
     */
    Point();

    /**
     * Construct the point <tt>(a, 0, 0)</tt>.
     */
    explicit Point(const double a);

    /**
     * Construct the point <tt>(b, 0, 0)</tt>.
     */
    explicit Point(const double a, const double b);

    /**
     * Construct the point <tt>(a, b, c)</tt>.
     */
    explicit Point(const double a, const double b, const double c);

    /**
     * Construct a unit vector in the specified dimension.
     *
     * \param axis The axis to construct the unit vector for. 0 corresponds to
     *    the x-direction, 1 to the y-direction, and 2 to the z-direction.
     */
    static Point
    unit_vector(const unsigned int axis);

    /**
     * Element-wise assignment to a scalar value.
     */
    Point&
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
    operator==(const Point& other) const;

    /**
     * Return whether any elements of two points are different.
     */
    bool
    operator!=(const Point& other) const;

    /**
     * Return the Euclidean distance between two points. This is given by \f$
     * d = \sqrt{(p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2} \f$.
     */
    double
    distance(const Point& other) const;

    /**
     * Return the Euclidean distance between two points squared. This is given
     * by \f$ d^2 = (p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2 \f$.
     */
    double
    distance_squared(const Point& other) const;

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
    Point&
    operator*=(const double factor);

    /**
     * Return a point multiplied by a scalar.
     */
    Point
    operator*(const double factor) const;

    /**
     * Element-wise division by a non-zero scalar.
     */
    Point&
    operator/=(const double factor);

    /**
     * Return a point divided by a non-zero scalar.
     */
     Point
     operator/(const double factor) const;

    /**
      * Element-wise negation in place.
      */
    Point&
    operator-();

    /**
     * Return the negative of a point.
     */
    Point
    operator-() const;

    /**
     * Element-wise addition in place.
     */
    Point&
    operator+=(const Point& other);

    /**
     * Return the sum of two points.
     */
    Point
    operator+(const Point& other) const;

    /**
     * Element-wise subtraction in place.
     */
    Point&
    operator-=(const Point& other);

    /**
     * Return the difference of two points.
     */
    Point
    operator-(const Point& other) const;

    /**
     * Element-wise absolute value in place.
     */
    Point&
    fabs();

    /**
     * Return the absolute value of a point.
     */
    Point
    fabs() const;

    /**
     * Normalize the point to its length. If the length is zero (e.g. the point
     * is the origin), return a copy. See \ref length.
     */
    Point&
    normalize();

    /**
     * Return the direction of a point. See \ref normalize.
     */
    Point
    direction() const;

    /**
     * Take the dot product between two points via \f$ c = p \cdot q = p_x q_x
     * + p_y q_y + p_z q_z \f$.
     */
    double
    dot(const Point& other) const;

    /**
     * Take the cross product between two points via \f$ r = p \times q =
     * (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x) \f$.
     */
    Point
    cross(const Point& other) const;

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
  Point
  operator*(const double factor, const Point& p);

  /**
   * Take the dot product between two points. See \ref Point::dot.
   */
  double
  dot(const Point& p, const Point& q);

  /**
   * Take the cross product between two points. See \ref Point::cross
   */
  Point
  cross(const Point& p, const Point& q);

  /**
   * Return the Euclidean distance between two points. See \ref Point::distance.
   */
  double
  distance(const Point& p, const Point& q);

  /**
   * Return the absolute value of a point. See \ref Point::fabs.
   */
  Point
  fabs(const Point& p);

  /**
   * Return the direction of a vector pointing to the specified point.
   * See \ref Point::directions
   */
  Point
  direction(const Point& p);

  /**
   * Insert a point into an output stream.
   */
  std::ostream&
  operator<<(std::ostream& os, const Point& p);

  //================================================== Useful Aliases

  using Vertex = Point;
  using Node = Point;
  using Centroid = Point;
  using Normal = Point;
  using Gradient = Point;

}
#endif //POINT_H

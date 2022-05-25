#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cstddef>

namespace pdes::Grid
{

/**
 * A class representing a point vector in space.
 *
 * In this context, a point is defined as a generic three-vector. This carries
 * a number of different interpretations, such as a scalar point or a vector.
 * Example uses for this class are representing vertices on a mesh, centroids
 * on cells or faces, nodes in a spatial discretization, outward normal vectors
 * on faces, or gradient vectors.
 */
class Point
{
public:
  /**
   * Access to the underlying value type.
  */
  using value_type = double;

public:
  value_type x; ///< The x-coordinate
  value_type y; ///< The y-coordinate
  value_type z; ///< The z-coordinate

public:

  Point();

  /** Construct a 1D point <tt>(a, 0, 0)</tt>. */
  explicit
  Point(const value_type a);

  /** Construct a 2D point <tt>(a, b, 0)</tt>. */
  explicit
  Point(const value_type a,
        const value_type b);

  /** Construct a 3D point <tt>(a, b, c)</tt>. */
  explicit
  Point(const value_type a,
        const value_type b,
        const value_type c);

  /** Set all elements of the point to a scalar value. */
  Point&
  operator=(const double value);

  /**
   * Static method to construct a unit vector in the specified dimension.
   *
   * \param axis The axis to construct the unit vector for. 0 corresponds to
   *    the x-direction, 1 to the y-direction, and 2 to the z-direction.
   * \throw Out of range when <tt>i > 2</tt>.
   */
  static Point
  unit_vector(const size_t axis);

  /** Test the equality of two points. */
  bool
  operator==(const Point& q) const;

  /** Test the inequality of two points. */
  bool
  operator!=(const Point& q) const;

  /** \name Accessors */
  // @{

  /**
   * Read and write access for element \p i of the point.
   *
   * \param i The coordinate to access. 0 corresponds to the x-coordinate, 1 to
   *    the y-coordinate, and 2 to the z-coordinate.
   * \throw Out of range when <tt>i > 2</tt>.
   */
  value_type&
  operator[](const size_t i);

  /**
   * Read access to element \p i of the point.
   * \see Point::operator[](const size_t i)
   */
  const value_type&
  operator[](const size_t i) const;


  /**
   * Read and write access to element \p i of the point.
   * \see Point::operator[](const size_t i)
   */
  value_type&
  operator()(const size_t i);

  /**
   * Read access to element \p i of the point.
   * \see Point::operator[](const size_t i) const
   */
  const value_type&
  operator()(const size_t i) const;

  // @}
  /** \name Characteristics */
  // @{

  /**
   * Return the length of the point vector. This is equivalent to returning
   * the Euclidean distance to the origin, computed via
   * \f$ d = || p ||_{\ell_2} = \sqrt{ x^2 + y^2 + z^2 | \f$.
   *
   * \see Point::length_squared
   */
  value_type
  length() const;

  /**
   * Return the length squared of the point vector. This is equivalent to
   * returning the squared Euclidean distance to the origin, given by
   * \f$ d^2 = || p ||^2_{\ell_2} = x^2 + y^2 + z^2 \f$.
   *
   * \see Point::length
   */
  value_type
  length_squared() const;

  // @}
  /** \name Scalar Operations */
  //@{

  /**
   * Negate each element in the point. This is computed via
   * \f$ p = -p = (-x, -y, -z) \f$.
   */
  Point&
  operator-();

  /**
   * Return the negated point.
   * \see operator-()
   */
  Point
  operator-() const;


  /**
   * Multiply each element in the point by \p factor. This computes
   * \f$ p = \alpha p = (\alpha x, \alpha y, \alpha z) \f$.
   */
  Point&
  operator*=(const value_type factor);

  /**
   * Divide each element in the point by \p factor. This computes
   * \f$ p = \frac{p}{\alpha}
   *       = (\frac{x}{\alpha}, \frac{y}{\alpha}, \frac{z}{\alpha})
   * \f$.
   */
  Point&
  operator/=(const value_type factor);

  // @}
  /** \name Point Operations */
  // @{

  /**
   * Add another point. This computes
   * \f$ p = p + q = (p_x + q_x, p_y + q_y, p_z + q_z) \f$.
   */
  Point&
  operator+=(const Point& q);

  /**
   * Subtract another point. This computes
   * \f$ p = p - q = (p_x - q_x, p_y - q_y, p_z - q_z) \f$.
   */
  Point&
  operator-=(const Point& q);

  /**
   * Return the dot product between this and another point. This computes
   * \f$ c = p_x q_x + p_y q_y + p_z q_z \f$.
   */
  value_type
  dot(const Point& q) const;

  /**
   * Return the cross product between this and another point.
   * This computes
   * \f$ r = p \times q
   *       = (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x)
   * \f$.
   */
  Point
  cross(const Point& q) const;

  /**
   * Return the Eulicdean distance between this point and another.
   * This computes
   * \f$ d = || p - q ||_{\ell_2}
   *       = \sqrt{ (p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2 }
   * \f$.
   *
   * \see Point::distance_squared
   */
  value_type
  distance(const Point& q) const;

  /**
   * Return the squared Euclidean distance between this point and another.
   * This computes
   * \f$ d = || p - q ||^2_{\ell_2}
   *       = (p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2
   * \f$.
   *
   * \see Point::distance
   */
  double
  distance_squared(const Point& q) const;

  /**
   * Take the absolute value of each element of the point. This computes
   * \f$ p = | p | = (|x|, |y|, |z|) \f$.
   */
  Point&
  fabs();

  /**
   * Return a point containing the absolute value of each element of the point.
   * \see Point::fabs()
   */
  Point
  fabs() const;

  /**
   * Normalize the point to its length. If the point is uniformly zero, nothing
   * is done and the original point is returned. This computes
   * \f$ p = \frac{p}{||p||_{\ell_2}} \f$.
   *
   * \see Point::length Point::direction
   */
  Point&
  normalize();

  /**
   * Return the direction, or unit-length, vector for this point. This computes
   * \f$ \hat{p} = \frac{p}{|| p ||_{\ell_2}} \f$.
   *
   * \see Point::length Point::normalize
   */
  Point
  direction() const;

  // @}
  /** \name Print Utilities */
  // @{

  /** Return the point represented as a string. */
  std::string
  str() const;

  /** Print the point to the standard output. */
  void
  print(std::ostream& os = std::cout) const;

  // @}

};

/*-------------------- Method Declarations --------------------*/

/**
 * Multiply each element of a point by a scalar value.
 * \see Point::operator*=
 */
Point
operator*(const Point& p, const double factor);

/**
 * Multiply each element of a point by a scalar value.
 * \see Point::operator*=
 */
Point
operator*(const double factor, const Point& p);

/**
 * Divide each element of a point by a scalar value.
 * \see Point::operator/=
 */
Point
operator/(const Point& p, const double factor);


/**
 * Add two points together.
 * \see Point::operator+=
 */
Point
operator+(const Point& p, const Point& q);

/**
 * Subtract two points.
 * \see Point::operator-=
 */
Point
operator-(const Point& p, const Point& q);

/**
 * Return the dot product of two points.
 * \see Point::dot
 */
double
dot(const Point& p, const Point& q);

/**
 * Return the cross product of two points.
 * \see Point::cross
 */
Point
cross(const Point& p, const Point& q);

/**
 * Return the Euclidean distance between two points.
 * \see Point::distance
 */
double
distance(const Point& p, const Point& q);

/**
 * Return the absolute value of a point.
 * \see Point::fabs
 */
Point
fabs(const Point& p);

/**
 * Return the direction, or unit length, point vector.
 * \see Point::direction
 */
Point
direction(const Point& p);

/**
 * Insert the string into an output stream.
 * \see Point::print Print::str
 */
std::ostream&
operator<<(std::ostream& os, const Point& p);

/*-------------------- Useful Aliases --------------------*/

using Vertex = Point;
using Node = Point;
using Centroid = Point;
using Normal = Point;
using Gradient = Point;

}
#endif //POINT_H

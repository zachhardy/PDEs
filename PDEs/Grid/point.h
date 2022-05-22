#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <sstream>

#include <cmath>
#include <cinttypes>

#include "macros.h"

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
  /**
   * Construct a point at the origin <tt>(0, 0, 0)</tt>.
   */
  Point() :
      x(0.0), y(0.0), z(0.0)
  {}

  /**
   * Construct a 1D point <tt>(a, 0, 0)</tt>.
   * \param a The x-coordinate.
   */
  explicit
  Point(const value_type a) :
      x(a), y(0.0), z(0.0)
  {}

  /**
   * Construct a 2D point <tt>(a, b, 0)</tt>.
   * \param a The x-coordinate.
   * \param b The y-coordinate.
   */
  explicit
  Point(const value_type a,
        const value_type b) :
      x(a), y(b), z(0.0)
  {}

  /**
   * Construct a 3D point <tt>(a, b, c)</tt>.
   * \param a The x-coordinate.
   * \param b The y-coordinate.
   * \param c The z-coordinate.
   */
  explicit
  Point(const value_type a,
        const value_type b,
        const value_type c) :
      x(a), y(b), z(c)
  {}

  /**
   * Set all elements of the point to a single scalar value.
   * \param value The value to set each element to.
   */
  Point&
  operator=(const double value)
  {
    x = value;
    y = value;
    z = value;
    return *this;
  }

  /**
   * Static method to construct a unit vector in the specified dimension.
   * \param axis The axis to construct the unit vector for. 0 corresponds to
   *    the x-direction, 1 to the y-direction, and 2 to the z-direction.
   * \throw Out of range when <tt>i > 2</tt>.
   */
  static Point
  unit_vector(const size_t axis)
  {
    Assert(axis < 3, "Invalid dimension provided.");
    if (axis == 0) return Point(1.0, 0.0, 0.0);
    else if (axis == 1) return Point(0.0, 1.0, 0.0);
    else return Point(0.0, 0.0, 1.0);
  }

  /**
   * Test the equality of two points.
   */
  bool
  operator==(const Point& q) const
  { return (x == q.x && y == q.y && z == q.z); }

  /**
   * Test the inequality of two points.
   */
  bool
  operator!=(const Point& q) const
  { return (x != q.x || y != q.y || z != q.z); }

  /** \name Accessors */
  // @{

  /**
   * Read and write access for element \p i of the point.
   * \param i The coordinate to access. 0 corresponds to the x-coordinate, 1 to
   *    the y-coordinate, and 2 to the z-coordinate.
   * \throw Out of range when <tt>i > 2</tt>.
   */
  value_type&
  operator[](const size_t i)
  {
    Assert(i < 3, "Invalid dimension provided.");
    if (i == 0) return x;
    else if (i == 1) return y;
    else return z;
  }

  /**
   * Read access to element \p i of the point.
   * \see Point::operator[](const size_t i)
   */
  const value_type&
  operator[](const size_t i) const
  {
    Assert(i < 3, "Invalid dimension provided.");
    if (i == 0) return x;
    else if (i == 1) return y;
    else return z;
  }


  /**
   * Read and write access to element \p i of the point.
   * \see Point::operator[](const size_t i)
   */
  value_type&
  operator()(const size_t i)
  { return (*this)[i]; }

  /**
   * Read access to element \p i of the point.
   * \see Point::operator[](const size_t i) const
   */
  const value_type&
  operator()(const size_t i) const
  { return (*this)[i]; }


  // @}
  /** \name Characteristics */
  // @{

  /**
   * Return the length of the point vector. This is equivalent to returning the
   * Euclidean distance to the origin, computed via
   * \f$ d = || p ||_{\ell_2} = \sqrt{ x^2 + y^2 + z^2 | \f$.
   * \see Point::length_squared
   */
  value_type
  length() const
  { return std::sqrt(length_squared()); }

  /**
   * Return the length squared of the point vector. This is equivalent to
   * returning the squared Euclidean distance to the origin, computed via
   * \f$ d^2 = || p ||^2_{\ell_2} = x^2 + y^2 + z^2 \f$.
   * \see Point::length
   */
  value_type
  length_squared() const
  { return x*x + y*y + z*z; }

  // @}
  /** \name Scalar Operations */
  //@{

  /**
   * Negate each element in the point. This computes
   * \f$ p = -p = (-x, -y, -z) \f$.
   */
  Point&
  operator-()
  {
    x = -x;
    y = -y;
    z = -z;
    return *this;
  }

  /**
   * Return the negated point.
   * \see operator-()
   */
  Point
  operator-() const
  { return Point(-x, -y, -z); }


  /**
   * Multiply each element in the point by \p factor. This computes
   * \f$ p = \alpha p = (\alpha x, \alpha y, \alpha z) \f$.
   */
  Point&
  operator*=(const value_type factor)
  {
    x *= factor;
    y *= factor;
    z *= factor;
    return *this;
  }

  /**
   * Divide each element in the point by \p factor. This computes
   * \f$ p = \frac{p}{\alpha}
   *       = (\frac{x}{\alpha}, \frac{y}{\alpha}, \frac{z}{\alpha})
   * \f$.
   */
  Point&
  operator/=(const value_type factor)
  {
    Assert(factor != 0.0, "Zero division error.");
    x /= factor;
    y /= factor;
    z /= factor;
    return *this;
  }

  // @}
  /** \name Point Operations */
  // @{

  /**
   * Add another point. This computes
   * \f$ p = p + q = (p_x + q_x, p_y + q_y, p_z + q_z) \f$.
   */
  Point&
  operator+=(const Point& q)
  {
    x += q.x;
    y += q.y;
    z += q.z;
    return *this;
  }

  /**
   * Subtract another point. This computes
   * \f$ p = p - q = (p_x - q_x, p_y - q_y, p_z - q_z) \f$.
   */
  Point&
  operator-=(const Point& q)
  {
    x -= q.x;
    y -= q.y;
    z -= q.z;
    return *this;
  }

  /**
   * Return the dot product between this and another point. This computes
   * \f$ c = p_x q_x + p_y q_y + p_z q_z \f$.
   */
  value_type
  dot(const Point& q) const
  { return x * q.x + y * q.y * z * q.z; }

  /**
   * Return the cross product between this and another point. This computes
   * \f$ r = p \times q
   *       = (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x)
   * \f$.
   */
  inline Point
  cross(const Point& q) const
  {
    return Point(y * q.z - z * q.y,
                 z * q.x - x * q.z,
                 x * q.y - y * q.x);
  }

  /**
   * Return the Eulicdean distance between this point and another.
   * This computes
   * \f$ d = || p - q ||_{\ell_2}
   *       = \sqrt{ (p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2 }
   * \f$.
   * \see Point::distance_squared
   */
  value_type
  distance(const Point& q) const
  { return std::sqrt(distance_squared(q)); }

  /**
   * Return the squared Eulicdean distance between this point and another.
   * This computes
   * \f$ d = || p - q ||^2_{\ell_2}
   *       = (p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2
   * \f$.
   * \see Point::distance
   */
  double
  distance_squared(const Point& q) const
  {
    double dx = x - q.x;
    double dy = y - q.y;
    double dz = z - q.z;
    return dx*dx + dy*dy + dz*dz;
  }

  /**
   * Take the absolute value of each element of the point. This computes
   * \f$ p = | p | = (|x|, |y|, |z|) \f$.
   */
  Point&
  fabs()
  {
    x = std::fabs(x);
    y = std::fabs(y);
    z = std::fabs(z);
    return *this;
  }

  /**
   * Return a point containing the absolute value of each element of the point.
   * \see Point::fabs()
   */
  inline Point
  fabs() const
  { return Point(x, y, z).fabs(); }

  /**
   * Normalize the point to its length. If the point is uniformly zero, nothing
   * is done and the original point is returned. This computs
   * \f$ p = \frac{p}{||p||_{\ell_2}} \f$.
   * \see Point::length Point::direction
   */
  Point&
  normalize()
  {
    double len = length();
    *this /= (len != 0.0)? len : 1.0;
    return *this;
  }

  /**
   * Return the direction, or unit-length, vector for this point. This computes
   * \f$ \hat{p} = \frac{p}{|| p ||_{\ell_2}} \f$.
   * \see Point::length Point::normalize
   */
  Point
  direction() const
  { return Point(x, y, z).normalize(); }

  // @}
  /** \name Print Utilities */
  // @{

  /**
   * Return the point represented as a string.
   */
  std::string
  str() const
  {
    std::stringstream ss;
    ss << "Point(" << x << " " << y << " " << z << ")\n";
    return ss.str();
  }


  /**
   * Print the point to the standard output.
   */
  void
  print(std::ostream& os = std::cout) const
  { os << str(); }

  // @}

};

/*-------------------- Method Declarations --------------------*/

/**
 * Multiply each element of a point by a scalar value.
 * \see Point::operator*=
 */
inline Point
operator*(const Point& p, const double factor)
{ return Point(p) *= factor; }

/**
 * Multiply each element of a point by a scalar value.
 * \see Point::operator*=
 */
inline Point
operator*(const double factor, const Point& p)
{ return Point(p) *= factor; }

/**
 * Divide each element of a point by a scalar value.
 * \see Point::operator/=
 */
inline Point
operator/(const Point& p, const double factor)
{ return Point(p) /= factor; }


/**
 * Add two points together.
 * \see Point::operator+=
 */
inline Point
operator+(const Point& p, const Point& q)
{ return Point(p) += q; }

/**
 * Subtract two points.
 * \see Point::operator-=
 */
inline Point
operator-(const Point& p, const Point& q)
{ return Point(p) -= q; }

/**
 * Return the dot product of two points.
 * \see Point::dot
 */
inline double
dot(const Point& p, const Point& q)
{ return p.dot(q); }

/**
 * Return the cross product of two points.
 * \see Point::cross
 */
inline Point
cross(const Point& p, const Point& q)
{ return p.cross(q); }

/**
 * Return the Euclidean distance between two points.
 * \see Point::distance
 */
inline double
distance(const Point& p, const Point& q)
{ return p.distance(q); }

/**
 * Return the absolute value of a point.
 * \see Point::fabs
 */
inline Point
fabs(const Point& p)
{ return p.fabs(); }

/**
 * Return the direction, or unit length, point vector.
 * \see Point::direction
 */
inline Point
direction(const Point& p)
{ return p.direction(); }

/**
 * Insert the string into an output stream.
 */
inline std::ostream&
operator<<(std::ostream& os, const Point& p)
{ return os << p.str(); }

/*-------------------- Useful Aliases --------------------*/

using Vertex = Point;
using Node = Point;
using Centroid = Point;
using Normal = Point;
using Gradient = Point;

}
#endif //POINT_H

#ifndef POINT_H
#define POINT_H

#include <array>
#include <cmath>

#include <iostream>
#include <sstream>
#include <cinttypes>

namespace grid
{

/**
 * A struct for representint a point vector in Cartesian space.
 *
 * In this context, a point is defined as a generic three-vector. This implies
 * that a point represents a vector from the origin to a specified (\p x, \p y,
 * \p z) coordinate. Example use cases for this struct are vertices or
 * centroids for cells and faces of a spatial mesh, nodes for spatial
 * discretizations, outward normal vectors on faces, gradient vectors, etc.
 */
class Point
{
public:   /*---------- Public Attributes ----------*/

  double x; ///< The x-coordinate.
  double y; ///< The y-coordinate.
  double z; ///< The z-coordinate.

public:

  /// Define a point at the origin `(0.0, 0.0, 0.0)`.
  Point() : x(0.0), y(0.0), z(0.0) {}
  /// Define the 1D point <tt>(a, 0.0, 0.0)</tt>.
  explicit Point(const double a) : x(a), y(0.0), z(0.0) {}
  /// Define the 2D point <tt>(a, b, 0.0)</tt>.
  explicit Point(const double a, const double b) : x(a), y(b), z(0.0) {}
  /// Define the 3D point <tt>(a, b, c)</tt>.
  explicit Point(const double a, const double b, const double c)
    : x(a), y(b), z(c)
  {}

  /// Copy constructor.
  Point(const Point& other)
    : x(other.x), y(other.y), z(other.z)
  {}

  /// Assignment operator.
  Point& operator=(const Point& other)
  {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

public:
  /** \name Access Operators */
  /** @{ */

  /// Read/write access for element \p index.
  double& operator[](const uint64_t index)
  {
    if (index > 2)
      this->indexing_error(__FUNCTION__);

    if      (index == 0)  return x;
    else if (index == 1)  return y;
    else              return z;
  }
  /// Read only access for element \p index.
  double operator[](const uint64_t index) const
  {
    if (index > 2)
      this->indexing_error(__FUNCTION__);

    if      (index == 0)  return x;
    else if (index == 1)  return y;
    else              return z;
  }

  /** @} */
public:
  /** \name Scalar Operations */
  /** @{ */

  /**
   * \brief Element-wise negation.
   * \f[ q = -p = (-p_x, -p_y, -p_z) \f]
   */
  Point operator-() const { return Point(-x, -y, -z); }
  /// See \ref operator-() const
  Point& operator-() { x = -x; y = -y; z = -z; return *this;}

  /**
   * \brief Element-wise multiplication by a scalar value.
   * \f[ q = p \alpha = (p_x \alpha, p_y \alpha, p_z \alpha) \f]
   */
  Point operator*(const double value) const
  {
    Point p;
    p.x = value * x;
    p.y = value * y;
    p.z = value * z;
    return p;
  }
  /// See \ref operator*(const double value) const
  Point& operator*=(const double value)
  {
    x *= value;
    y *= value;
    z *= value;
    return *this;
  }

  /**
   * \brief Element-wise division by a scalar value.
   * \f[ q = \frac{p}{\alpha}
   *       = \left(
   *            \frac{p_x}{\alpha},
   *            \frac{p_y}{\alpha},
   *            \frac{p_z}{\alpha}
   *         \right)
   * \f]
   */
  Point operator/(const double value) const
  {
    if (value == 0.0)
      this->zero_division_error(__FUNCTION__);

    Point p;
    p.x = x / value;
    p.y = y / value;
    p.z = z / value;
    return p;
  }
   /// See \ref operator/(const double value) const
  Point& operator/=(const double value)
  {
    if (value == 0.0)
      this->zero_division_error(__FUNCTION__);

    x /= value;
    y /= value;
    z /= value;
    return *this;
  }

  /** @} */
  /** \name Point Operations */
  /** @{ */

  /**
   * \brief Element-wise addition of two points.
   * \f[ s = p + q = (p_x + q_x, p_y + q_y, p_z + q_z) \f]
   */
  Point operator+(const Point& other) const
  {
    Point p;
    p.x = x + other.x;
    p.y = y + other.y;
    p.z = z + other.z;
    return p;
  }
  /// See \ref operator+(const Point& other) const
  Point& operator+=(const Point& other)
  {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  /**
   * \brief Element-wise subtraction of two points.
   * \f[ s = p - q = (p_x - q_x, p_y - q_y, p_z - q_z) \f]
   */
  Point operator-(const Point& other) const
  {
    Point p;
    p.x = x - other.x;
    p.y = y - other.y;
    p.z = z - other.z;
    return p;
  }
  /// See \ref operator-(const Point& other) const
  Point& operator-=(const Point& other)
  {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  /**
   * \brief Element-wise multiplication of two points.
   * \f[ s = p q = (p_x q_x, p_y q_y, p_z q_z) \f]
   */
  Point operator*(const Point& other) const
  {
    Point p;
    p.x = x * other.x;
    p.y = y * other.y;
    p.z = z * other.z;
    return p;
  }
  /// See \ref operator*(const Point& other) const
  Point& operator*=(const Point& other)
  {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
  }

  /**
   * \brief Element-wise division of two points.
   * \f[
   *    s = \frac{p}{q}
   *      = \left( \frac{p_x}{q_x}, \frac{p_y}{q_y}, \frac{p_z}{q_z} \right)
   * \f]
   */
  Point operator/(const Point& other) const
  {
    if (other.has_zero_elements())
      this->zero_division_error(__FUNCTION__);

    Point p;
    p.x = x / other.x;
    p.y = y / other.y;
    p.z = z / other.z;
    return p;
  }
  /// See \ref operator/(const Point& other) const
  Point& operator/=(const Point& other)
  {
    if (other.has_zero_elements())
      this->zero_division_error(__FUNCTION__);

    x /= other.x;
    y /= other.y;
    z /= other.z;
    return *this;
  }

  /**
   * \brief Return the length, or \f$\ell_2\f$-norm, of this point.
   * \f[ l = ||p||_{\ell_2} = \sqrt{ p_x^2 + p_y^2 + p_z^2 } \f]
   */
  double length() const { return sqrt(x*x + y*y + z*z); }

  /// Return the square of \ref length.
  double length_squared() const { return x*x + y*y + z*z; }

  /**
   * \brief Normalize the point vector to its length.
   * \f[ \hat{p} = \frac{p}{\sqrt{ p_x^2 + p_y^2 + p_z^2 }} \f]
   */
  Point& normalize()
  {
    double length = this->length();
    if (length == 0.0)
      return *this;

    x /= length;
    y /= length;
    z /= length;
    return *this;
  }

  /// Element-wise absolute value in place.
  Point& abs() { x = fabs(x); y = fabs(y); z = fabs(z); return *this; }

  /**
   * \brief Return the dot product between two points.
   * \f[
   * c = \vec{p} \cdot \vec{q}
   *   = p_x q_x + p_y q_y + p_z q_z
   * \f]
   */
  double dot(const Point& other) const
  { return x*other.x + y*other.y + z*other.z; }

  /**
   * \brief Return the cross product between two points.
   * \f[
   * \vec{s} = \vec{p} \times \vec{q}
   *         = (p_y q_z - p_z q_y,
   *            p_z q_x - p_x q_z,
   *            p_x q_y - p_y q_x)
   * \f]
   */
  Point cross(const Point& other) const
  {
    Point p;
    p.x = y*other.z - z*other.y;
    p.y = z*other.x - x*other.z;
    p.z = x*other.y - y*other.x;
    return p;
  }

  /**
   * \brief Return the Euclidean distance between two points.
   * \f[
   * d = || \vec{p} - \vec{q} ||_{\ell_2}
   *   = \sqrt{(p_x - q_x)^2 + (p_y - q_y)^2 + (p_z - q_z)^2}.
   * \f]
   */
  double distance(const Point& other) const
  {return this->operator-(other).length();}

  /// Return the square of \ref distance.
  double distance_squared(const Point& other) const
  {return this->operator-(other).length_squared();}

  /** @} */
public:
  /** \name Print Utilities */
  /** @{ */

  /// Return the point as a string.
  std::string to_string() const
  {
    std::stringstream ss;
    ss  << "[" << x << " " << y << " " << z << "]";
    return ss.str();
  }

  /// Print this point to `std::cout`.
  void print() const { std::cout << this->to_string(); }

  /// Insert this point into an output stream.
  friend std::ostream& operator<<(std::ostream& os, const Point& p)
  {
    os << "[" << p.x << " " << p.y << " " << p.z << "]";
    return os;
  }

  /** @} */
private:

  /// Determine whether zero elements exist.///
  bool has_zero_elements() const { return (x==0.0 or y==0.0 or z==0.0); }

  /// Throw an error for division by zero.///
  void zero_division_error(const std::string func_name) const
  {
    std::stringstream err;
    err << "Point::" << func_name << ": Zero division encountered.";
    throw std::runtime_error(err.str());
  }

  /// Throw an error for out of range indexing.///
  void indexing_error(const std::string func_name) const
  {
    std::stringstream err;
    err << "Point::" << func_name << ": Invalid index encountered.";
    throw std::out_of_range(err.str());
  }

};

/*-------------------- Inline Implementations --------------------*/
/**
 * \brief Element-wise multiplication of a point and a scalar value.
 * \f[ q = \alpha p = (\alpha p_x, \alpha p_y, \alpha p_z) \f]
 */
inline Point operator*(const double value, const Point& p) { return p * value; }

/**
 * \brief Return the dot product between two points.
 * \f[ v = p \cdot q = p_x q_x + p_y q_y + p_z q_z \f]
 */
inline double dot(const Point& p, const Point& q) { return p.dot(q); }

/**
 * \brief Return the cross product between two points.
 * \f[ s = p \times q = (p_y q_z - p_z q_y,
 *                       p_z q_x - p_x q_z,
 *                       p_x q_y - p_y q_x)
 * \f]
 */
inline Point cross(const Point& p, const Point& q) { return p.cross(q); }

/**
 * \brief Return the Euclidean distance between two points.
 * \f[
 * d = || \vec{v} - \vec{q} ||_{\ell_2}
 *   = \sqrt{(v_x - q_x)^2 + (v_y - q_y)^2 + (p_z - q_z)^2}.
 * \f]
 */
inline double distance(const Point& p, const Point& q) { return p.distance(q); }

/**
 * \brief Return the length, or \f$ \ell_2 \f$-norm of a point.
 * \f[ l = ||p||_{\ell_2} = \sqrt{ p_x^2 + p_y^2 + p_z^2 } \f]
 */
inline double length(const Point& p) { return p.length(); }

/**
 * \brief Return the direction of a point vector.
 * \f[ \hat{q} = \frac{p}{||p||_{\ell_2}
 *             = \frac{p}{\sqrt{ p_x^2 + p_y^2 + p_z^2 }}
 * \f]
 */
inline Point direction(const Point& p) { return Point(p)/p.length(); }

/// Return the absolute value of a point vector.
inline Point abs(const Point& p) { return Point(fabs(p.x), fabs(p.y), fabs(p.z)); }

}
#endif //POINT_H

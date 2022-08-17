#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cstddef>


namespace Grid
{

  /**
   * A class representing a point in space or a vector.
   *
   * A point is defined as a three-vector. This is intended to be used to
   * define various quantities such as vertices, centroids, nodes, normal
   * vectors, or gradient vectors.
   */
  class Point
  {
  public:
    double x;
    double y;
    double z;

  public:

    //################################################## Constructors

    /** \name Constructors and Initialization */
    // @{

    /** Construct a point at the origin <tt>(0, 0, 0)</tt>. */
    Point();

    /** Construct a 1D point <tt>(a, 0, 0)</tt>. */
    explicit Point(const double a);

    /** Construct a 2D point <tt>(a, b, 0)</tt>. */
    explicit Point(const double a, const double b);

    /** Construct a 3D point <tt>(a, b, c)</tt>. */
    explicit Point(const double a,
                   const double b,
                   const double c);

    /**
     * Construct a unit vector in the specified dimension.
     *
     * \param axis The axis to construct the unit vector for. 0 corresponds to
     *    the x-direction, 1 to the y-direction, and 2 to the z-direction.
     */
    static Point unit_vector(const unsigned int axis);

    /** Element-wise assignment to a scalar value. */
    Point& operator=(const double value);

    // @}

    //================================================== Accessors

    /** \name Accessors */
    // @{

    double& operator[](const unsigned int i);
    const double& operator[](const unsigned int i) const;

    double& operator()(const unsigned int i);
    const double& operator()(const unsigned int i) const;

    // @}

    //################################################## Information

    /** \name Information */
    // @{

    bool operator==(const Point& q) const;
    bool operator!=(const Point& q) const;

    /** Return the Euclidean distance to the origin. */
    double length() const;

    /** Return the Euclidean distance to the origin squared. */
    double length_squared() const;

    // @}

    //================================================== Scalar Operations

    /** \name Scalar Operations */
    //@{

    /** Element-wise negation in place. */
    Point& operator-();

    /** Return a Point containing the negated elements. */
    Point operator-() const;

    /** Element-wise multiplication by a scalar in place. */
    Point& operator*=(const double factor);

    /** Element-wise division by a scalar in place. */
    Point& operator/=(const double factor);

    // @}

    //################################################## Point Operations

    /** \name Point Operations */
    // @{

    /** Element-wise addition with another Point in place. */
    Point& operator+=(const Point& q);

    /** Element-wise subtraction with another Point in place. */
    Point& operator-=(const Point& q);

    /**
     * Take the dot product with another Point via \f$ c = p_x q_x + p_y q_y +
     * p_z q_z \f$.
     */
    double dot(const Point& q) const;

    /**
     * Take the cross product with another Point via  \f$ r = p \times q =
     * (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x) \f$.
     */
    Point cross(const Point& q) const;

    /** Return the Euclidean distance to another Point. */
    double distance(const Point& q) const;

    /** Return the Euclidean distance to another Point squared. */
    double distance_squared(const Point& q) const;

    /** Element-wise absolute value in place. */
    Point& fabs();

    /** Return a Point containing the absolute value of the elements. */
    Point fabs() const;

    /** Normalize the Point to its length. */
    Point& normalize();

    /**
     * Return a Point containing the normalized elements. This is the direction
     * of a vector pointing from the origin to the coordinates defined by the
     * Point.
     */
    Point direction() const;

    // @}

    //================================================== Print Utilities

    /** \name Print Utilities */
    // @{

    std::string str() const;
    void print(std::ostream& os = std::cout) const;

    // @}

  };

  //================================================== Methods

  /** Element-wise multiplication by a scalar. */
  Point operator*(const Point& p, const double factor);

  /** Element-wise multiplication by a scalar. */
  Point operator*(const double factor, const Point& p);

  /** Element-wise division by a scalar. */
  Point operator/(const Point& p, const double factor);

  /** Element-wise addition. */
  Point operator+(const Point& p, const Point& q);

  /** Element-wise subtraction. */
  Point operator-(const Point& p, const Point& q);

  /** Compute a dot product. \see Point::dot */
  double dot(const Point& p, const Point& q);

  /** Compute a cross product. \see Point::cross */
  Point cross(const Point& p, const Point& q);

  /** Compute the distance between Point \p p and \p q. */
  double distance(const Point& p, const Point& q);

  /** Return the absolute value of the elements of a Point. */
  Point fabs(const Point& p);

  /** Return the direction of the vector pointing from the origin to a Point. */
  Point direction(const Point& p);

  std::ostream& operator<<(std::ostream& os, const Point& p);

  //================================================== Useful Aliases

  using Vertex = Point;
  using Node = Point;
  using Centroid = Point;
  using Normal = Point;
  using Gradient = Point;

}
#endif //POINT_H

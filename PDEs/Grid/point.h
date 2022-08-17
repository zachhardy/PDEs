#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cstddef>


namespace Grid
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
     * Static method to construct a unit vector in the specified dimension.
     *
     * \param axis The axis to construct the unit vector for. 0 corresponds to
     *    the x-direction, 1 to the y-direction, and 2 to the z-direction.
     */
    static Point unit_vector(const unsigned int axis);

    /** Set all elements of the point to a scalar value. */
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

    /** Return the distance of the point from the origin. */
    double length() const;

    /** Return the squared distance of the point from the origin. */
    double length_squared() const;

    // @}

    //================================================== Scalar Operations

    /** \name Scalar Operations */
    //@{

    /** Element-wise negation in-place. */
    Point& operator-();

    /** Return an element-wise negated Point. */
    Point operator-() const;

    /** Element-wise multiplication by a scalar in-place. */
    Point& operator*=(const double factor);

    /** Element-wise division by a scalar in-place. */
    Point& operator/=(const double factor);

    // @}

    //################################################## Point Operations

    /** \name Point Operations */
    // @{

    /** Element-wise addition in-place. */
    Point& operator+=(const Point& q);

    /** Element-wise subtraction in-place. */
    Point& operator-=(const Point& q);

    /**
     * Return the dot product between this and another point via
     * \f$ c = p_x q_x + p_y q_y + p_z q_z \f$.
     */
    double dot(const Point& q) const;

    /**
     * Return the cross product between this and another point via
     * \f$ r = p \times q
     *       = (p_y q_z - p_z q_y, p_z q_x - p_x q_z, p_x q_y - p_y q_x) \f$.
     */
    Point cross(const Point& q) const;

    /** Return the Euclidean distance between two points. */
    double distance(const Point& q) const;

    /** Return the squared Euclidean distance between two points. */
    double distance_squared(const Point& q) const;

    /** Element-wise absolute value in-place. */
    Point& fabs();

    /** Return the absolute value of this Point. */
    Point fabs() const;

    /** Normalize the Point to its length. */
    Point& normalize();

    /**
     * Return the direction of the vector from the origin to this Point.
     * This simply copies Point and then normalizes it.
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

  /** Element-wise multiplication by a scalar value. */
  Point operator*(const Point& p, const double factor);

  /** Element-wise multiplication by a scalar value. */
  Point operator*(const double factor, const Point& p);

  /** Element-wise division by a scalar value. */
  Point operator/(const Point& p, const double factor);

  /** Element-wise addition. */
  Point operator+(const Point& p, const Point& q);

  /** Element-wise subtraction. */
  Point operator-(const Point& p, const Point& q);

  /** Return the dot product of two points. \see Point::dot */
  double dot(const Point& p, const Point& q);

  /** Return the cross product of two points. \see Point::cross */
  Point cross(const Point& p, const Point& q);

  /** Return the Euclidean distance between two points. */
  double distance(const Point& p, const Point& q);

  /** Return the element-wise absolute value of a point. */
  Point fabs(const Point& p);

  /** Return the direction, or unit length vector. */
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

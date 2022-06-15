#include "point.h"
#include "macros.h"

#include <cmath>


using namespace Grid;


//################################################## Constructors


Point::Point()
  : x(0.0), y(0.0), z(0.0)
{}


Point::Point(const double a)
  : x(a), y(0.0), z(0.0)
{}


Point::Point(const double a,
             const double b)
  : x(a), y(b), z(0.0)
{}


Point::Point(const double a,
             const double b,
             const double c)
  : x(a), y(b), z(c)
{}


//################################################## Assignment


Point&
Point::operator=(const double value)
{
  x = value;
  y = value;
  z = value;
  return *this;
}


//################################################## Static Methods


Point
Point::unit_vector(const std::size_t axis)
{
  Assert(axis < 3, "Invalid dimension provided.");
  if (axis == 0) return Point(1.0, 0.0, 0.0);
  else if (axis == 1) return Point(0.0, 1.0, 0.0);
  else return Point(0.0, 0.0, 1.0);
}


//################################################## Comparison Operators


bool
Point::operator==(const Point& q) const
{ return (x == q.x && y == q.y && z == q.z); }


bool
Point::operator!=(const Point& q) const
{ return !(*this == q); }


//################################################## Accessors


double&
Point::operator[](const size_t i)
{
  Assert(i < 3, "Invalid dimension provided.");
  if (i == 0) return x;
  else if (i == 1) return y;
  else return z;
}


const double&
Point::operator[](const size_t i) const
{
  Assert(i < 3, "Invalid dimension provided.");
  if (i == 0) return x;
  else if (i == 1) return y;
  else return z;
}


double&
Point::operator()(const size_t i)
{ return (*this)[i]; }


const double&
Point::operator()(const size_t i) const
{ return (*this)[i]; }


//################################################## Information


double
Point::length() const
{ return std::sqrt(length_squared()); }


double
Point::length_squared() const
{ return x*x + y*y + z*z; }


//################################################## Scalar Operations


Point&
Point::operator-()
{
  x = -x;
  y = -y;
  z = -z;
  return *this;
}


Point
Point::operator-() const
{ return Point(-x, -y, -z); }


Point&
Point::operator*=(const double factor)
{
  x *= factor;
  y *= factor;
  z *= factor;
  return *this;
}


Point&
Point::operator/=(const double factor)
{
  Assert(factor != 0.0, "Zero division error.");
  x /= factor;
  y /= factor;
  z /= factor;
  return *this;
}


//################################################## Point Operations


Point&
Point::operator+=(const Point& q)
{
  x += q.y;
  y += q.y;
  z += q.z;
  return *this;
}


Point&
Point::operator-=(const Point& q)
{
  x -= q.y;
  y -= q.y;
  z -= q.z;
  return *this;
}


double
Point::dot(const Point& q) const
{ return x*q.x + y*q.y*z*q.z; }


Point
Point::cross(const Point& q) const
{
  return Point(y*q.z - z*q.y,
               z*q.x - x*q.z,
               x*q.y - y*q.x);
}


double
Point::distance(const Point& q) const
{ return std::sqrt(distance_squared(q)); }


double
Point::distance_squared(const Point& q) const
{
  double dx = x - q.x;
  double dy = y - q.y;
  double dz = z - q.z;
  return dx*dx + dy*dy + dz*dz;
}


Point&
Point::fabs()
{
  x = std::fabs(x);
  y = std::fabs(y);
  z = std::fabs(z);
  return *this;
}


Point
Point::fabs() const
{ return Point(x, y, z).fabs(); }


Point&
Point::normalize()
{
  double len = length();
  return *this /= (len != 0.0)? len : 1.0;
}


Point
Point::direction() const
{ return Point(x, y, z).normalize(); }


//################################################## Print Utilities


std::string
Point::str() const
{
  std::stringstream ss;
  ss << "Point(" << x << " " << y << " " << z << ")\n";
  return ss.str();
}


void
Point::print(std::ostream& os) const
{ os << str(); }


//################################################## Methods


Point
Grid::operator*(const Point& p, const double factor)
{ return Point(p) *= factor; }


Point
Grid::operator*(const double factor, const Point& p)
{ return Point(p) *= factor; }


Point
Grid::operator/(const Point& p, const double factor)
{ return Point(p) /= factor; }


Point
Grid::operator+(const Point& p, const Point& q)
{ return Point(p) += q; }


Point
Grid::operator-(const Point& p, const Point& q)
{ return Point(p) -= q; }


double
Grid::dot(const Point& p, const Point& q)
{ return p.dot(q); }


Point
Grid::cross(const Point& p, const Point& q)
{ return p.cross(q); }


double
Grid::distance(const Point& p, const Point& q)
{ return p.distance(q); }


Point
Grid::fabs(const Point& p)
{ return p.fabs(); }


Point
Grid::direction(const Point& p)
{ return p.direction(); }


std::ostream&
Grid::operator<<(std::ostream& os, const Point& p)
{ return os << p.str(); }

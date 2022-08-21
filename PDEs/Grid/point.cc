#include "point.h"

#include <sstream>
#include <cmath>
#include <cassert>


using namespace Grid;


Point::Point() :
  values({0.0, 0.0, 0.0})
{}


Point::Point(const double a) :
  values({a, 0.0, 0.0})
{}


Point::Point(const double a, const double b) :
  values({a, b, 0.0})
{}


Point::Point(const double a, const double b, const double c) :
  values({a, b, c})
{}


Point
Point::unit_vector(const unsigned int axis)
{
  assert(axis <= 2);
  if (axis == 0) return Point(1.0, 0.0, 0.0);
  else if (axis == 1) return Point(0.0, 1.0, 0.0);
  else return Point(0.0, 0.0, 1.0);
}


Point&
Point::operator=(const double value)
{
  values = {value, value, value};
  return *this;
}


double&
Point::operator[](const unsigned int i)
{
  return values.at(i);
}


const double&
Point::operator[](const unsigned int i) const
{
  return values.at(i);
}


double&
Point::operator()(const unsigned int i)
{
  return values.at(i);
}


const double&
Point::operator()(const unsigned int i) const
{
  return values.at(i);
}


const double&
Point::x() const
{
  return values[0];
}


const double&
Point::y() const
{
  return values[1];
}


const double&
Point::z() const
{
  return values[2];
}


bool
Point::operator==(const Point& other) const
{
  return (values == other.values);
}


bool
Point::operator!=(const Point& other) const
{
  return !(values == other.values);
}


double
Point::distance(const Point& other) const
{
  return std::sqrt(this->distance_squared(other));
}


double
Point::distance_squared(const Point& other) const
{
  double retval = 0.0;
  const auto difference = *this - other;
  for (unsigned int i = 0; i <= 2; ++i)
    retval = difference[i] * difference[i];
  return retval;
}


double
Point::length() const
{
  return std::sqrt(this->length_squared());
}


double
Point::length_squared() const
{
  double retval = 0.0;
  for (const auto& v : values)
    retval += v * v;
  return retval;
}


Point&
Point::operator*=(const double factor)
{
  for (unsigned int i = 0; i <= 2; ++i)
    values[i] *= factor;
  return *this;
}


Point
Point::operator*(const double factor) const
{
  return Point(*this) *= factor;
}


Point&
Point::operator/=(const double factor)
{
  assert(factor != 0.0);
  return *this *= 1.0/factor;
}


Point
Point::operator/(const double factor) const
{
  return Point(*this) /= factor;
}


Point&
Point::operator-()
{
  return *this *= -1.0;
}


Point
Point::operator-() const
{
  return -Point(*this);
}


Point&
Point::operator+=(const Point& other)
{
  for (unsigned int i = 0; i <= 2; ++i)
    values[i] += other.values[i];
  return *this;
}


Point
Point::operator+(const Point& other) const
{
  return Point(*this) += other;
}


Point&
Point::operator-=(const Point& other)
{
  for (unsigned int i = 0; i <= 2; ++i)
    values[i] -= other.values[i];
  return *this;
}


Point
Point::operator-(const Point& other) const
{
  return Point(*this) -= other;
}


Point&
Point::fabs()
{
  for (auto& v : values)
    v = std::fabs(v);
  return *this;
}


Point
Point::fabs() const
{
  return Point(*this).fabs();
}


Point&
Point::normalize()
{
  double len = this->length();
  return (len > 0.0)? *this /= len : *this;
}


Point
Point::direction() const
{
  return Point(*this).normalize();
}


double
Point::dot(const Point& other) const
{
  double retval = 0.0;
  for (unsigned int i = 0; i <= 2; ++i)
    retval += values[i] * other.values[i];
  return retval;
}


Point
Point::cross(const Point& other) const
{
  double x = this->y()*other.z() - this->z()*other.y();
  double y = this->z()*other.x() - this->x()*other.z();
  double z = this->x()*other.y() - this->y()*other.x();
  return Point(x, y, z);
}


std::string
Point::str() const
{
  std::stringstream ss;
  ss << "Point("
     << this->x() << ", "
     << this->y() << ", "
     << this->z() << ")" << std::endl;
  return ss.str();
}


void
Point::print(std::ostream& os) const
{
  os << this->str();
}


Point
Grid::operator*(const double factor, const Point& p)
{
  return p * factor;
}


double
Grid::dot(const Point& p, const Point& q)
{
  return p.dot(q);
}


Point
Grid::cross(const Point& p, const Point& q)
{
  return p.cross(q);
}


double
Grid::distance(const Point& p, const Point& q)
{
  return p.distance(q);
}


Point
Grid::fabs(const Point& p)
{
  return p.fabs();
}


Point
Grid::direction(const Point& p)
{
  return p.direction();
}


std::ostream&
Grid::operator<<(std::ostream& os, const Point& p)
{
  return os << p.str();
}

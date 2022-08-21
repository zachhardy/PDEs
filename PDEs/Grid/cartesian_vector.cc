#include "cartesian_vector.h"

#include <sstream>
#include <cmath>
#include <cassert>


using namespace Grid;


CartesianVector::CartesianVector() :
  values({0.0, 0.0, 0.0})
{}


CartesianVector::CartesianVector(const double a) :
  values({a, 0.0, 0.0})
{}


CartesianVector::CartesianVector(const double a, const double b) :
  values({a, b, 0.0})
{}


CartesianVector::CartesianVector(const double a, const double b, const double c) :
  values({a, b, c})
{}


CartesianVector
CartesianVector::unit_vector(const unsigned int axis)
{
  assert(axis <= 2);
  if (axis == 0) return CartesianVector(1.0, 0.0, 0.0);
  else if (axis == 1) return CartesianVector(0.0, 1.0, 0.0);
  else return CartesianVector(0.0, 0.0, 1.0);
}


CartesianVector&
CartesianVector::operator=(const double value)
{
  values = {value, value, value};
  return *this;
}


double&
CartesianVector::operator[](const unsigned int i)
{
  return values.at(i);
}


const double&
CartesianVector::operator[](const unsigned int i) const
{
  return values.at(i);
}


double&
CartesianVector::operator()(const unsigned int i)
{
  return values.at(i);
}


const double&
CartesianVector::operator()(const unsigned int i) const
{
  return values.at(i);
}


const double&
CartesianVector::x() const
{
  return values[0];
}


const double&
CartesianVector::y() const
{
  return values[1];
}


const double&
CartesianVector::z() const
{
  return values[2];
}


bool
CartesianVector::operator==(const CartesianVector& other) const
{
  return (values == other.values);
}


bool
CartesianVector::operator!=(const CartesianVector& other) const
{
  return !(values == other.values);
}


double
CartesianVector::distance(const CartesianVector& other) const
{
  return std::sqrt(this->distance_squared(other));
}


double
CartesianVector::distance_squared(const CartesianVector& other) const
{
  double dx = this->x() - other.x();
  double dy = this->y() - other.y();
  double dz = this->z() - other.z();
  return dx*dx + dy*dy + dz*dz;
}


double
CartesianVector::length() const
{
  return std::sqrt(this->length_squared());
}


double
CartesianVector::length_squared() const
{
  double retval = 0.0;
  for (const auto& v : values)
    retval += v * v;
  return retval;
}


CartesianVector&
CartesianVector::operator*=(const double factor)
{
  for (unsigned int i = 0; i <= 2; ++i)
    values[i] *= factor;
  return *this;
}


CartesianVector
CartesianVector::operator*(const double factor) const
{
  return CartesianVector(*this) *= factor;
}


CartesianVector&
CartesianVector::operator/=(const double factor)
{
  assert(factor != 0.0);
  return *this *= 1.0/factor;
}


CartesianVector
CartesianVector::operator/(const double factor) const
{
  return CartesianVector(*this) /= factor;
}


CartesianVector&
CartesianVector::operator-()
{
  return *this *= -1.0;
}


CartesianVector
CartesianVector::operator-() const
{
  return -CartesianVector(*this);
}


CartesianVector&
CartesianVector::operator+=(const CartesianVector& other)
{
  for (unsigned int i = 0; i <= 2; ++i)
    values[i] += other.values[i];
  return *this;
}


CartesianVector
CartesianVector::operator+(const CartesianVector& other) const
{
  return CartesianVector(*this) += other;
}


CartesianVector&
CartesianVector::operator-=(const CartesianVector& other)
{
  for (unsigned int i = 0; i <= 2; ++i)
    values[i] -= other.values[i];
  return *this;
}


CartesianVector
CartesianVector::operator-(const CartesianVector& other) const
{
  return CartesianVector(*this) -= other;
}


CartesianVector&
CartesianVector::fabs()
{
  for (auto& v : values)
    v = std::fabs(v);
  return *this;
}


CartesianVector
CartesianVector::fabs() const
{
  return CartesianVector(*this).fabs();
}


CartesianVector&
CartesianVector::normalize()
{
  double len = this->length();
  return (len > 0.0)? *this /= len : *this;
}


CartesianVector
CartesianVector::direction() const
{
  return CartesianVector(*this).normalize();
}


double
CartesianVector::dot(const CartesianVector& other) const
{
  double retval = 0.0;
  for (unsigned int i = 0; i <= 2; ++i)
    retval += values[i] * other.values[i];
  return retval;
}


CartesianVector
CartesianVector::cross(const CartesianVector& other) const
{
  double x = this->y()*other.z() - this->z()*other.y();
  double y = this->z()*other.x() - this->x()*other.z();
  double z = this->x()*other.y() - this->y()*other.x();
  return CartesianVector(x, y, z);
}


std::string
CartesianVector::str() const
{
  std::stringstream ss;
  ss << "CartesianVector("
     << this->x() << ", "
     << this->y() << ", "
     << this->z() << ")" << std::endl;
  return ss.str();
}


void
CartesianVector::print(std::ostream& os) const
{
  os << this->str();
}


CartesianVector
Grid::operator*(const double factor, const CartesianVector& p)
{
  return p * factor;
}


double
Grid::dot(const CartesianVector& p, const CartesianVector& q)
{
  return p.dot(q);
}


CartesianVector
Grid::cross(const CartesianVector& p, const CartesianVector& q)
{
  return p.cross(q);
}


double
Grid::distance(const CartesianVector& p, const CartesianVector& q)
{
  return p.distance(q);
}


CartesianVector
Grid::fabs(const CartesianVector& p)
{
  return p.fabs();
}


CartesianVector
Grid::direction(const CartesianVector& p)
{
  return p.direction();
}


std::ostream&
Grid::operator<<(std::ostream& os, const CartesianVector& p)
{
  return os << p.str();
}

#include "cartesian_vector.h"

#include <sstream>
#include <cmath>
#include <cassert>

using namespace PDEs;
namespace PDEs::Grid
{

  //------------------------------------------------------------
  // Constructors and assignment
  //------------------------------------------------------------

  CartesianVector::CartesianVector() : xyz({0.0, 0.0, 0.0}) {}


  CartesianVector::CartesianVector(const double a) :
      xyz({a, 0.0, 0.0})
  {}


  CartesianVector::
  CartesianVector(const double a, const double b) :
      xyz({a, b, 0.0})
  {}


  CartesianVector::
  CartesianVector(const double a, const double b, const double c) :
      xyz({a, b, c})
  {}


  CartesianVector
  CartesianVector::unit_vector(const unsigned int axis)
  {
    assert(axis < 3);
    return axis == 0? CartesianVector(1.0, 0.0, 0.0) :
           axis == 1? CartesianVector(0.0, 1.0, 0.0) :
                      CartesianVector(0.0, 0.0, 1.0);
  }

  //------------------------------------------------------------
  // Comparison
  //------------------------------------------------------------

  bool CartesianVector::
  operator==(const CartesianVector& other) const
  {
    return this->xyz == other.xyz;
  }


  bool CartesianVector::
  operator!=(const CartesianVector& other) const
  {
    return this->xyz != other.xyz;
  }

  //------------------------------------------------------------
  // Accessors
  //------------------------------------------------------------

  const double& CartesianVector::
  operator[](const unsigned int i) const
  {
    return xyz.at(i);
  }


  double& CartesianVector::operator[](const unsigned int i)
  {
    return xyz.at(i);
  }


  const double& CartesianVector::
  operator()(const unsigned int i) const
  {
    return xyz.at(i);
  }


  double& CartesianVector::
  operator()(const unsigned int i)
  {
    return xyz.at(i);
  }


  const double& CartesianVector::x() const { return xyz[0]; }
  double& CartesianVector::x() { return xyz[0]; }

  const double& CartesianVector::y() const { return xyz[1]; }
  double& CartesianVector::y() { return xyz[1]; }

  const double& CartesianVector::z() const { return xyz[2]; }
  double& CartesianVector::z() { return xyz[2]; }

  double* CartesianVector::data() { return xyz.data(); }
  const double* CartesianVector::data() const { return xyz.data(); }

  //------------------------------------------------------------
  // Addition and Subtraction
  //------------------------------------------------------------

  CartesianVector& CartesianVector::
  operator+=(const CartesianVector& other)
  {
    this->x() += other.x();
    this->y() += other.y();
    this->z() += other.z();
    return *this;
  }


  CartesianVector CartesianVector::
  operator+(const CartesianVector& other) const
  {
    return CartesianVector(*this) += other;
  }


  CartesianVector& CartesianVector::
  operator-=(const CartesianVector& other)
  {
    this->x() -= other.x();
    this->y() -= other.y();
    this->z() -= other.z();
    return *this;
  }


  CartesianVector CartesianVector::
  operator-(const CartesianVector& other) const
  {
    return CartesianVector(*this) -= other;
  }

  //------------------------------------------------------------
  // Multiplication and division
  //------------------------------------------------------------

  CartesianVector& CartesianVector::
  scale(const double factor)
  {
    this->x() *= factor;
    this->y() *= factor;
    this->z() *= factor;
    return *this;
  }

  CartesianVector CartesianVector::
  scale(const double factor) const
  {
    return CartesianVector(*this).scale(factor);
  }


  CartesianVector& CartesianVector::operator-()
  {
    return this->scale(-1.0);
  }

  CartesianVector CartesianVector::operator-() const
  {
    return -CartesianVector(*this);
  }


  CartesianVector& CartesianVector::
  operator*=(const double factor)
  {
    return this->scale(factor);
  }


  CartesianVector CartesianVector::
  operator*(const double factor) const
  {
    return CartesianVector(*this) *= factor;
  }


  CartesianVector& CartesianVector::
  operator/=(const double factor)
  {
    assert(factor != 0.0);
    return this->scale(1.0 / factor);
  }


  CartesianVector CartesianVector::
  operator/(const double factor) const
  {
    return CartesianVector(*this) /= factor;
  }

  //------------------------------------------------------------
  // Other operations
  //------------------------------------------------------------

  CartesianVector& CartesianVector::fabs()
  {
    this->x() = std::fabs(this->x());
    this->y() = std::fabs(this->y());
    this->z() = std::fabs(this->z());
    return *this;
  }


  CartesianVector CartesianVector::fabs() const
  {
    return CartesianVector(*this).fabs();
  }


  double CartesianVector::dot(const CartesianVector& other) const
  {
    return this->x() * other.x() +
           this->y() * other.y() +
           this->z() * other.z();
  }

  CartesianVector CartesianVector::
  cross(const CartesianVector& other) const
  {
    double x = this->y() * other.z() - this->z() * other.y();
    double y = this->z() * other.x() - this->x() * other.z();
    double z = this->x() * other.y() - this->y() * other.x();
    return CartesianVector(x, y, z);
  }


  double CartesianVector::length() const
  {
    return std::sqrt(this->dot(*this));
  }


  double CartesianVector::length_squared() const
  {
    return this->dot(*this);
  }


  CartesianVector& CartesianVector::normalize()
  {
    double len = this->length();
    if (len == 0.0)
      return *this;

    this->x() /= len;
    this->y() /= len;
    this->z() /= len;
    return *this;
  }


  CartesianVector CartesianVector::direction() const
  {
    return CartesianVector(*this).normalize();
  }


  double CartesianVector::distance(const CartesianVector& other) const
  {
    return std::sqrt(this->distance_squared(other));
  }


  double CartesianVector::
  distance_squared(const CartesianVector& other) const
  {
    double dx = this->x() - other.x();
    double dy = this->y() - other.y();
    double dz = this->z() - other.z();
    return dx * dx + dy * dy + dz * dz;
  }

  //------------------------------------------------------------
  // Print utilities
  //------------------------------------------------------------

  std::string CartesianVector::str() const
  {
    std::stringstream ss;
    ss << "("
       << this->x() << ", " << this->y() << ", " << this->z()
       << ")";
    return ss.str();
  }


  void CartesianVector::print() const
  {
    std::cout << this->str();
  }

  //------------------------------------------------------------
  // Friends
  //------------------------------------------------------------

  CartesianVector operator*(const double factor,
                            const CartesianVector& p)
  {
    return p * factor;
  }


  std::ostream& operator<<(std::ostream& os,
                           const CartesianVector& p)
  {
    return os << p.str();
  }

  double dot(const CartesianVector& p,
             const CartesianVector& q)
  {
    return p.dot(q);
  }


  CartesianVector cross(const CartesianVector& p,
                        const CartesianVector& q)
  {
    return p.cross(q);
  }


  double length(const CartesianVector& p)
  {
    return p.length();
  }


  double length_squared(const CartesianVector& p)
  {
    return p.length_squared();
  }


  double distance(const CartesianVector&p,
                  const CartesianVector& q)
  {
    return p.distance(q);
  }


  double distance_squared(const CartesianVector& p,
                          const CartesianVector& q)
  {
    return p.distance_squared(q);
  }


  CartesianVector fabs(const CartesianVector& p)
  {
    return p.fabs();
  }
}

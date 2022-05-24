#include "vector.h"

#include <cmath>

#include <iostream>
#include <sstream>
#include <iomanip>

#include "macros.h"


using namespace pdes::Math;

//################################################## Constructors

Vector::Vector(const size_t n) :
  elements(n)
{}


Vector::Vector(const size_t n, const double value) :
  elements(n, value)
{}


Vector::Vector(const std::vector<double>& other) :
  elements(other)
{}


Vector::Vector(std::vector<double>&& other) :
  elements(other)
{}


Vector::Vector(const std::initializer_list<double> list) :
  elements(list)
{}

//################################################## Assignment

Vector&
Vector::operator=(const std::vector<double>& other)
{
  elements = other;
  return *this;
}


Vector&
Vector::operator=(std::vector<double>&& other)
{
  elements = other;
  return *this;
}


Vector&
Vector::operator=(const std::initializer_list<double> list)
{
  elements = list;
  return *this;
}


Vector&
Vector::operator=(const double value)
{
  for (auto& el : elements)
    el = value;
  return *this;
}

//################################################## Comparison Operators

bool
Vector::operator==(const Vector& y) const
{ return (elements == y.elements); }


bool
Vector::operator!=(const Vector& y) const
{ return (elements != y.elements); }


//################################################## Characteristics

size_t
Vector::size() const
{ return elements.size(); }


bool
Vector::empty() const
{ return elements.empty(); }


bool
Vector::all_zero() const
{
  for (const auto& el : elements)
    if (el != 0.0) return false;
  return true;
}

//################################################## Iterators

std::vector<double>::iterator
Vector::begin()
{ return elements.begin(); }


std::vector<double>::iterator
Vector::end()
{ return elements.end(); }


std::vector<double>::const_iterator
Vector::begin() const
{ return elements.begin(); }


std::vector<double>::const_iterator
Vector::end() const
{ return elements.end(); }


//################################################## Accessors

double&
Vector::operator[](const size_t i)
{ return elements[i]; }


const double&
Vector::operator[](const size_t i) const
{ return elements[i]; }


double&
Vector::operator()(const size_t i)
{ return elements[i]; }


const double&
Vector::operator()(const size_t i) const
{ return elements[i]; }


double&
Vector::at(const size_t i)
{ return elements.at(i); }


const double&
Vector::at(const size_t i) const
{ return elements.at(i); }


double&
Vector::front()
{ return elements.front(); }


const double&
Vector::front() const
{ return elements.front(); }


double&
Vector::back()
{ return elements.back(); }


const double&
Vector::back() const
{ return elements.back(); }


double*
Vector::data()
{ return elements.data(); }


const double*
Vector::data() const
{ return elements.data(); }

//################################################## Modifiers

void
Vector::clear()
{ elements.clear(); }


void
Vector::push_back(const double value)
{ elements.push_back(value); }


void
Vector::pop_back()
{ elements.pop_back(); }


void
Vector::resize(const size_t n)
{ elements.resize(n); }


void
Vector::resize(const size_t n, const double value)
{ elements.resize(n, value); }


void
Vector::swap(Vector& y)
{ elements.swap(y.elements); }


Vector&
Vector::normalize()
{
  double norm = l2_norm();
  return *this /= (norm != 0.0) ? norm : 1.0;
}


Vector
Vector::unit() const
{ return Vector(elements).normalize(); }


Vector&
Vector::fabs()
{
  for (auto& el : elements)
    el = std::fabs(el);
  return *this;
}


Vector
Vector::fabs() const
{ return Vector(elements).fabs(); }

//################################################## Scalar Operations

Vector&
Vector::operator-()
{
  for (auto& el : elements)
    el = -el;
  return *this;
}


Vector
Vector::operator-() const
{ return -Vector(elements); }


Vector&
Vector::operator*=(const double factor)
{
  for (auto& el : elements)
    el *= factor;
  return *this;
}


Vector&
Vector::operator/=(const double factor)
{
  Assert(factor != 0.0, "Zero division error.");
  for (auto& el : elements)
    el /= factor;
  return *this;
}

//################################################## Linear Algebra Operations

Vector&
Vector::operator+=(const Vector& y)
{
  Assert(size() == y.size(), "Dimension mismatch error.");
  for (size_t i = 0; i < size(); ++i)
    elements[i] += y.elements[i];
  return *this;
}


Vector&
Vector::operator-=(const Vector& y)
{
  Assert(size() == y.size(), "Dimension mismatch error.");
  for (size_t i = 0; i < size(); ++i)
    elements[i] -= y.elements[i];
  return *this;
}


double
Vector::dot(const Vector& y) const
{
  Assert(size() == y.size(), "Dimension mismatch error.");
  double c = 0.0;
  for (size_t i = 0; i < size(); ++i)
    c += elements[i] * y.elements[i];
  return c;
}

//################################################## Vector Norms

double
Vector::linf_norm() const
{
  double norm = 0.0;
  for (const auto& el : elements)
    if (std::fabs(el) > norm)
      norm = std::fabs(el);
  return norm;
}


double
Vector::l1_norm() const
{
  double norm = 0.0;
  for (const auto& el : elements)
    norm += std::fabs(el);
  return norm;
}


double
Vector::l2_norm() const
{
  double norm = 0.0;
  for (const auto& el : elements)
    norm += std::fabs(el) * std::fabs(el);
  return std::sqrt(norm);
}


double
Vector::lp_norm(const double p) const
{
  double norm = 0.0;
  for (const auto& el : elements)
    norm += std::pow(std::fabs(el), p);
  return std::pow(norm, 1.0/p);
}

//################################################## Print Utilities

void
Vector::print(std::ostream& os,
              const bool scientific,
              const unsigned int precision,
              const unsigned int width) const
{
  unsigned int w                   = width;
  std::ios::fmtflags old_flags     = os.flags();
  unsigned int       old_precision = os.precision(precision);

  if (scientific)
  {
    os.setf(std::ios::scientific, std::ios::floatfield);
    w = (!width) ? precision + 7 : w;
  }
  else
  {
    os.setf(std::ios::fixed, std::ios::floatfield);
    w = (!width) ? precision + 4 : w;
  }

  os << "[";
  for (size_t i = 0; i < size() - 1; ++i)
    os << std::setw(w) << elements[i];
  os << std::setw(w) << elements.back() << "]\n";

  os.flags(old_flags);
  os.precision(old_precision);
}


std::string
Vector::str(const bool scientific,
            const unsigned int precision,
            const unsigned int width) const
{
  std::stringstream ss;
  print(ss, scientific, precision, width);
  return ss.str();
}


//################################################## Methods

Vector
pdes::Math::operator*(const Vector& x, const double factor)
{ return Vector(x) *= factor; }


Vector
pdes::Math::operator*(const double factor, const Vector& x)
{ return Vector(x) *= factor; }


Vector
pdes::Math::operator/(const Vector& x, const double factor)
{ return Vector(x) /= factor; }


Vector
pdes::Math::operator+(const Vector& x, const Vector& y)
{ return Vector(x) += y; }


Vector
pdes::Math::operator-(const Vector& x, const Vector& y)
{ return Vector(x) -= y; }


double
pdes::Math::dot(const Vector& x, const Vector& y)
{ return x.dot(y); }


Vector
pdes::Math::fabs(const Vector& x)
{ return x.fabs(); }


Vector
pdes::Math::unit(const Vector& x)
{ return x.unit(); }


std::ostream&
pdes::Math::operator<<(std::ostream& os, const Vector& x)
{ return os << x.str(); }

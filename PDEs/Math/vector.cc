#include "vector.h"
#include "macros.h"

#include <cmath>
#include <iomanip>
#include <algorithm>


using namespace Math;

//################################################## Constructors

Vector::Vector(const size_t n) :
  vals(n)
{}


Vector::Vector(const size_t n, const double value) :
  vals(n, value)
{}


Vector::Vector(const std::vector<double>& other) :
  vals(other)
{}


Vector::Vector(std::vector<double>&& other) :
  vals(other)
{}


Vector::Vector(const std::initializer_list<double> list) :
  vals(list)
{}


//################################################## Assignment


Vector&
Vector::operator=(const std::vector<double>& other)
{
  vals = other;
  return *this;
}


Vector&
Vector::operator=(std::vector<double>&& other)
{
  vals = other;
  return *this;
}


Vector&
Vector::operator=(const std::initializer_list<double> list)
{
  vals = list;
  return *this;
}


Vector&
Vector::operator=(const double value)
{
  for (auto& el : vals)
    el = value;
  return *this;
}


//################################################## Characteristics


size_t
Vector::size() const
{ return vals.size(); }


size_t
Vector::nnz() const
{
  return std::count_if(vals.begin(), vals.end(),
                       [](const value_type v)
                       { return v != 0.0; });
}


bool
Vector::empty() const
{ return vals.empty(); }


bool
Vector::all_zero() const
{
  for (const auto& el : vals)
    if (el != 0.0) return false;
  return true;
}


bool
Vector::operator==(const Vector& y) const
{ return (vals == y.vals); }


bool
Vector::operator!=(const Vector& y) const
{ return (vals != y.vals); }


//################################################## Accessors


double&
Vector::operator[](const size_t i)
{ return vals[i]; }


const double&
Vector::operator[](const size_t i) const
{ return vals[i]; }


double&
Vector::operator()(const size_t i)
{ return vals[i]; }


const double&
Vector::operator()(const size_t i) const
{ return vals[i]; }


double&
Vector::at(const size_t i)
{ return vals.at(i); }


const double&
Vector::at(const size_t i) const
{ return vals.at(i); }


double&
Vector::front()
{ return vals.front(); }


const double&
Vector::front() const
{ return vals.front(); }


double&
Vector::back()
{ return vals.back(); }


const double&
Vector::back() const
{ return vals.back(); }


double*
Vector::data()
{ return vals.data(); }


const double*
Vector::data() const
{ return vals.data(); }


Vector::iterator
Vector::begin()
{ return vals.begin(); }


Vector::iterator
Vector::end()
{ return vals.end(); }


Vector::const_iterator
Vector::begin() const
{ return vals.begin(); }


Vector::const_iterator
Vector::end() const
{ return vals.end(); }


//################################################## Modifiers


void
Vector::clear()
{ vals.clear(); }


void
Vector::push_back(const double value)
{ vals.push_back(value); }


void
Vector::pop_back()
{ vals.pop_back(); }


void
Vector::resize(const size_t n)
{ vals.resize(n); }


void
Vector::resize(const size_t n, const double value)
{ vals.resize(n, value); }


void
Vector::swap(Vector& y)
{ vals.swap(y.vals); }


//################################################## Scalar Operations and
//                                                   Vector Norms

double
Vector::dot(const Vector& y) const
{
  Assert(size() == y.size(), "Dimension mismatch error.");
  double c = 0.0;
  for (size_t i = 0; i < size(); ++i)
    c += vals[i]*y.vals[i];
  return c;
}


double
Vector::linfty_norm() const
{
  double norm = 0.0;
  for (const auto& el : vals)
    if (std::fabs(el) > norm)
      norm = std::fabs(el);
  return norm;
}


double
Vector::l1_norm() const
{
  double norm = 0.0;
  for (const auto& el : vals)
    norm += std::fabs(el);
  return norm;
}


double
Vector::l2_norm() const
{
  double norm = 0.0;
  for (const auto& el : vals)
    norm += std::fabs(el)*std::fabs(el);
  return std::sqrt(norm);
}


double
Vector::lp_norm(const double p) const
{
  double norm = 0.0;
  for (const auto& el : vals)
    norm += std::pow(std::fabs(el), p);
  return std::pow(norm, 1.0/p);
}


//################################################## Linear Algebra Operations


Vector&
Vector::scale(const value_type factor)
{
  for (auto& el : vals)
    el *= factor;
  return *this;
}


Vector&
Vector::scale(const Vector scaling_factors)
{
  Assert(scaling_factors.size() == size(), "Dimension mismatch error.");

  // Get pointers for faster access
  value_type* el_ptr = data();
  const value_type* f_ptr = scaling_factors.data();
  value_type* end_ptr = data() + vals.size();

  // Perform the vector scaling
  for (; el_ptr != end_ptr; ++el_ptr, ++f_ptr)
    *el_ptr *= *f_ptr;
  return *this;
}


Vector&
Vector::add(const value_type value)
{
  for (auto& el : vals)
    el += value;
  return *this;
}


Vector&
Vector::add(const Vector& y, const value_type a)
{
  Assert(y.size() == size(), "Dimension mismatch error.");

  // Get pointers for fast access
  value_type* x_ptr = data();
  const value_type* y_ptr = y.data();
  value_type* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr += a**y_ptr;
  return *this;
}


Vector&
Vector::sadd(const value_type a, const Vector& y)
{
  Assert(y.size() == size(), "Dimension mismatch error.");

  // Get pointers for fast access
  value_type* x_ptr = data();
  const value_type* y_ptr = y.data();
  value_type* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr = a**x_ptr + *y_ptr;
  return *this;
}


Vector&
Vector::sadd(const value_type a, const value_type b, const Vector& y)
{
  Assert(y.size() == size(), "Dimension mismatch error.");

  // Get pointers for fast access
  value_type* x_ptr = data();
  const value_type* y_ptr = y.data();
  value_type* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr = a**x_ptr + b**y_ptr;
  return *this;
}


Vector&
Vector::equal(const Vector& y, const value_type factor)
{
  assert(y.size() == size());

  // Get pointers for fast access
  value_type* x_ptr = data();
  const value_type* y_ptr = y.data();
  value_type* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr;)
    *x_ptr++ = factor**y_ptr++;
  return *this;
}


Vector&
Vector::fabs()
{
  for (auto& el : vals)
    el = std::fabs(el);
  return *this;
}


Vector&
Vector::operator-()
{ return scale(-1.0); }


Vector
Vector::operator-() const
{ return -Vector(vals); }


Vector&
Vector::operator*=(const double factor)
{ return scale(factor); }


Vector&
Vector::operator/=(const double factor)
{
  Assert(factor != 0.0, "Zero division error.");
  return scale(1.0/factor);
}


Vector&
Vector::operator+=(const Vector& y)
{
  Assert(y.size() == size(), "Dimension mismatch error.");

  // Get pointers for fast access
  value_type* x_ptr = data();
  const value_type* y_ptr = y.data();
  value_type* end_ptr = data() + size();

  // Add the other vector
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr += *y_ptr;
  return *this;
}


Vector&
Vector::operator-=(const Vector& y)
{
  Assert(y.size() == size(), "Dimension mismatch error.");

  // Get pointers for fast access
  value_type* x_ptr = data();
  const value_type* y_ptr = y.data();
  value_type* end_ptr = data() + size();

  // Subtract the other vector
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr -= *y_ptr;
  return *this;
}


//################################################## Print Utilities


void
Vector::print(std::ostream& os,
              const bool scientific,
              const unsigned int precision,
              const unsigned int width) const
{
  unsigned int w = width;
  std::ios::fmtflags old_flags = os.flags();
  unsigned int old_precision = os.precision(precision);

  if (scientific)
  {
    os.setf(std::ios::scientific, std::ios::floatfield);
    w = (!width)? precision + 10 : w;
  }
  else
  {
    os.setf(std::ios::fixed, std::ios::floatfield);
    w = (!width)? precision + 5 : w;
  }

  os << "[";
  for (size_t i = 0; i < size() - 1; ++i)
    os << std::setw(w) << vals[i];
  os << std::setw(w) << vals.back() << "]\n";

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
Math::operator*(const Vector& x, const double factor)
{ return Vector(x) *= factor; }


Vector
Math::operator*(const double factor, const Vector& x)
{ return Vector(x) *= factor; }


Vector
Math::operator/(const Vector& x, const double factor)
{ return Vector(x) /= factor; }


Vector
Math::operator+(const Vector& x, const Vector& y)
{ return Vector(x) += y; }


Vector
Math::operator-(const Vector& x, const Vector& y)
{ return Vector(x) -= y; }


double
Math::dot(const Vector& x, const Vector& y)
{ return x.dot(y); }


Vector
Math::fabs(const Vector& x)
{ return Vector(x).fabs(); }


double
Math::linfty_norm(const Vector& x)
{ return x.linfty_norm(); }


double
Math::l1_norm(const Vector& x)
{ return x.l1_norm(); }


double
Math::l2_norm(const Vector& x)
{ return x.l2_norm(); }


double
Math::lp_norm(const Vector& x, const double p)
{ return x.lp_norm(p); }


std::ostream&
Math::operator<<(std::ostream& os, const Vector& x)
{ return os << x.str(); }

#include "vector.h"
#include "macros.h"

#include <cmath>
#include <iomanip>
#include <algorithm>


using namespace Math;


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


size_t
Vector::size() const
{
  return vals.size();
}


size_t
Vector::n_nonzero_elements() const
{
  return std::count_if(vals.begin(), vals.end(),
                       [](const double v)
                       { return v != 0.0; });
}


bool
Vector::empty() const
{
  return vals.empty();
}


bool
Vector::operator==(const Vector& y) const
{
  return (vals == y.vals);
}


bool
Vector::operator!=(const Vector& y) const
{
  return (vals != y.vals);
}


double&
Vector::operator[](const size_t i)
{
  return vals[i];
}


const double&
Vector::operator[](const size_t i) const
{
  return vals[i];
}


double&
Vector::operator()(const size_t i)
{
  return vals[i];
}


const double&
Vector::operator()(const size_t i) const
{
  return vals[i];
}


double&
Vector::at(const size_t i)
{
  return vals.at(i);
}


const double&
Vector::at(const size_t i) const
{
  return vals.at(i);
}


double&
Vector::front()
{
  return vals.front();
}


const double&
Vector::front() const
{
  return vals.front();
}


double&
Vector::back()
{
  return vals.back();
}


const double&
Vector::back() const
{
  return vals.back();
}


double*
Vector::data()
{
  return vals.data();
}


const double*
Vector::data() const
{
  return vals.data();
}


Vector::iterator
Vector::begin()
{
  return vals.begin();
}


Vector::iterator
Vector::end()
{
  return vals.end();
}


Vector::const_iterator
Vector::begin() const
{
  return vals.begin();
}


Vector::const_iterator
Vector::end() const
{
  return vals.end();
}


void
Vector::clear()
{
  vals.clear();
}


void
Vector::push_back(const double value)
{
  vals.push_back(value);
}


void
Vector::pop_back()
{
  vals.pop_back();
}


void
Vector::resize(const size_t n)
{
  vals.resize(n);
}


void
Vector::resize(const size_t n, const double value)
{
  vals.resize(n, value);
}


void
Vector::swap(Vector& other)
{
  vals.swap(other.vals);
}


double
Vector::dot(const Vector& y) const
{
  assert(this->size() == y.size());
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


Vector&
Vector::scale(const double factor)
{
  for (auto& el : vals)
    el *= factor;
  return *this;
}


Vector&
Vector::scale(const Vector scaling_factors)
{
  assert(scaling_factors.size() == this->size());

  // Get pointers for faster access
  double* el_ptr = data();
  const double* f_ptr = scaling_factors.data();
  double* end_ptr = data() + vals.size();

  // Perform the vector scaling
  for (; el_ptr != end_ptr; ++el_ptr, ++f_ptr)
    *el_ptr *= *f_ptr;
  return *this;
}


Vector&
Vector::shift(const double value)
{
  for (auto& el : vals)
    el += value;
  return *this;
}


Vector&
Vector::add(const Vector& y, const double a)
{
  assert(y.size() == this->size());

  // Get pointers for fast access
  double* x_ptr = data();
  const double* y_ptr = y.data();
  double* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr += a * *y_ptr;
  return *this;
}


Vector&
Vector::sadd(const double a, const Vector& y)
{
  assert(y.size() == this->size());

  // Get pointers for fast access
  double* x_ptr = data();
  const double* y_ptr = y.data();
  double* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr = a * *x_ptr + *y_ptr;
  return *this;
}


Vector&
Vector::sadd(const double a, const double b, const Vector& y)
{
  assert(y.size() == this->size());

  // Get pointers for fast access
  double* x_ptr = data();
  const double* y_ptr = y.data();
  double* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr = a * *x_ptr + b * *y_ptr;
  return *this;
}


Vector&
Vector::equal(const Vector& y, const double a)
{
  assert(y.size() == this->size());

  // Get pointers for fast access
  double* x_ptr = data();
  const double* y_ptr = y.data();
  double* end_ptr = data() + size();

  // Perform the add operation
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr = a * *y_ptr;
  return *this;
}


Vector&
Vector::fabs()
{
  for (auto& el : vals)
    el = std::fabs(el);
  return *this;
}


Vector
Vector::fabs() const
{
  return Vector(*this).fabs();
}


Vector&
Vector::operator-()
{
  return this->scale(-1.0);
}


Vector
Vector::operator-() const
{ return -Vector(vals); }


Vector&
Vector::operator*=(const double factor)
{
  return this->scale(factor);
}


Vector
Vector::operator*(const double factor) const
{
  return Vector(*this).scale(factor);
}


Vector&
Vector::operator/=(const double factor)
{
  assert(factor != 0.0);
  return this->scale(1.0/factor);
}


Vector
Vector::operator/(const double factor) const
{
  return Vector(*this) /= factor;
}


Vector&
Vector::operator+=(const Vector& y)
{
  assert(y.size() == this->size());

  // Get pointers for fast access
  double* x_ptr = data();
  const double* y_ptr = y.data();
  double* end_ptr = data() + size();

  // Add the other vector
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr += *y_ptr;
  return *this;
}


Vector
Vector::operator+(const Vector& y) const
{
  return Vector(*this) += y;
}


Vector&
Vector::operator-=(const Vector& y)
{
  assert(y.size() == this->size());

  // Get pointers for fast access
  double* x_ptr = data();
  const double* y_ptr = y.data();
  double* end_ptr = data() + size();

  // Subtract the other vector
  for (; x_ptr != end_ptr; ++x_ptr, ++y_ptr)
    *x_ptr -= *y_ptr;
  return *this;
}


Vector
Vector::operator-(const Vector& y) const
{
  return Vector(*this) -= y;
}


std::string
Vector::str(const bool scientific,
            const unsigned int precision,
            const unsigned int width) const
{
  std::stringstream ss;
  unsigned int w = width;

  if (scientific)
  {
    ss.setf(std::ios::scientific, std::ios::floatfield);
    w = (!width)? precision + 10 : w;
  }
  else
  {
    ss.setf(std::ios::fixed, std::ios::floatfield);
    w = (!width)? precision + 5 : w;
  }

  ss << "[";
  for (size_t i = 0; i < size() - 1; ++i)
    ss << std::setw(w) << vals[i];
  ss << std::setw(w) << vals.back() << "]\n";

  return ss.str();
}


void
Vector::print(std::ostream& os,
              const bool scientific,
              const unsigned int precision,
              const unsigned int width) const
{
  os << this->str(scientific, precision, width);
}


Vector
Math::operator*(const double factor, const Vector& x)
{
  return Vector(x) *= factor;
}


double
Math::dot(const Vector& x, const Vector& y)
{ return x.dot(y); }


Vector
Math::fabs(const Vector& x)
{
  return Vector(x).fabs();
}


double
Math::linfty_norm(const Vector& x)
{
  return x.linfty_norm();
}


double
Math::l1_norm(const Vector& x)
{
  return x.l1_norm();
}


double
Math::l2_norm(const Vector& x)
{
  return x.l2_norm();
}


double
Math::lp_norm(const Vector& x, const double p)
{
  return x.lp_norm(p);
}


std::ostream&
Math::operator<<(std::ostream& os, const Vector& x)
{
  return os << x.str();
}

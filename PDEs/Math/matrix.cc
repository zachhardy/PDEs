#include "matrix.h"
#include "macros.h"

#include <cmath>
#include <iomanip>

using namespace pdes::Math;

//################################################## Constructors

Matrix::Matrix(const size_t n_rows,
               const size_t n_cols) :
  coeffs(n_rows, Vector(n_cols))
{}


Matrix::Matrix(const size_t n_rows,
               const size_t n_cols,
               const double value) :
  coeffs(n_rows, Vector(n_cols, value))
{}


Matrix::Matrix(const STLMatrix& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.");

  for (const auto& row : other)
    coeffs.push_back(row);
}


Matrix::Matrix(STLMatrix&& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.");

  for (auto& row : other)
    coeffs.push_back(row);
}


Matrix::Matrix(const InitializerMatrix list)
{
  for (auto& row : list)
  {
    if (!coeffs.empty())
      Assert(row.size() == coeffs.front().size(),
             "Invalid input. All rows must be the same length.");

    coeffs.push_back(row);
  }
}

//################################################## Assignment Operators

Matrix&
Matrix::operator=(const STLMatrix& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.");

  for (const auto& row : other)
    coeffs.push_back(row);
  return *this;
}


Matrix&
Matrix::operator=(STLMatrix&& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.");

  for (auto& row : other)
    coeffs.push_back(row);
  return *this;
}


Matrix&
Matrix::operator=(const double value)
{
  Assert(!coeffs.empty(), "Cannot set an empty matrix to a scalar.");

  for (auto& row : coeffs)
    row = value;
  return *this;
}

//################################################## Comparison Operators

bool
Matrix::operator==(const Matrix& other) const
{ return (coeffs == other.coeffs); }


bool
Matrix::operator!=(const Matrix& other) const
{ return (coeffs != other.coeffs); }

//################################################## Characteristics

size_t
Matrix::n_rows() const
{ return coeffs.size(); }


size_t
Matrix::n_cols() const
{ return coeffs.front().size(); }


size_t
Matrix::size() const
{ return n_rows() * n_cols(); }


size_t
Matrix::nnz() const
{
  size_t count = 0;
  for (const auto& row : coeffs)
    for (const auto& el : row)
      if (el != 0.0) ++count;
  return count;
}


bool
Matrix::empty() const
{ return coeffs.empty(); }


bool
Matrix::all_zero() const
{ return (nnz() == 0)? true : false; }

//################################################## Iterators

std::vector<Vector>::iterator
Matrix::begin()
{ return coeffs.begin(); }


std::vector<Vector>::iterator
Matrix::end()
{ return coeffs.end(); }


std::vector<Vector>::const_iterator
Matrix::begin() const
{ return coeffs.begin(); }


std::vector<Vector>::const_iterator
Matrix::end() const
{ return coeffs.end(); }


std::vector<double>::iterator
Matrix::begin(const size_t i)
{ return coeffs.at(i).begin(); }


std::vector<double>::iterator
Matrix::end(const size_t i)
{ return coeffs.at(i).end(); }


std::vector<double>::const_iterator
Matrix::begin(const size_t i) const
{ return coeffs.at(i).begin(); }


std::vector<double>::const_iterator
Matrix::end(const size_t i) const
{ return coeffs.at(i).end(); }

//################################################## Accessors

Vector&
Matrix::operator[](const size_t i)
{ return coeffs[i]; }


const Vector&
Matrix::operator[](const size_t i) const
{ return coeffs[i]; }


Vector&
Matrix::operator()(const size_t i)
{ return coeffs[i]; }


const Vector&
Matrix::operator()(const size_t i) const
{ return coeffs[i]; }


Vector&
Matrix::at(const size_t i)
{ return coeffs.at(i); }


const Vector&
Matrix::at(const size_t i) const
{ return coeffs.at(i); }


double&
Matrix::operator()(const size_t i, const size_t j)
{ return coeffs[i][j]; }


const double&
Matrix::operator()(const size_t i, const size_t j) const
{ return coeffs[i][j]; }


double&
Matrix::at(const size_t i, const size_t j)
{ return coeffs.at(i).at(j); }


const double&
Matrix::at(const size_t i, const size_t j) const
{ return coeffs.at(i).at(j); }


double&
Matrix::diagonal(const size_t i)
{ return coeffs.at(i).at(i); }


const double&
Matrix::diagonal(const size_t i) const
{ return coeffs.at(i).at(i); }


Vector
Matrix::diagonal() const
{
  Vector diag;
  size_t min_dim = std::min(n_rows(), n_cols());
  for (size_t i = 0; i < min_dim; ++i)
    diag.push_back(coeffs[i][i]);
  return diag;
}


Vector*
Matrix::data()
{ return coeffs.data(); }


const Vector*
Matrix::data() const
{ return coeffs.data(); }


double*
Matrix::data(const size_t i)
{ return coeffs.at(i).data(); }


const double*
Matrix::data(const size_t i) const
{ return coeffs.at(i).data(); }

//################################################## Modifiers

void
Matrix::clear()
{ coeffs.clear(); }


void
Matrix::pop_back()
{ coeffs.pop_back(); }


void
Matrix::push_back(const Vector& row)
{
  Assert(row.size() == n_cols(), "Dimension mismatch error.");
  coeffs.push_back(row);
}


void
Matrix::push_back(Vector&& row)
{
  Assert(row.size() == n_cols(), "Dimension mismatch error.");
  coeffs.push_back(row);
}


void
Matrix::resize(const size_t n_rows,
               const size_t n_cols)
{ coeffs.resize(n_rows, Vector(n_cols)); }


void
Matrix::resize(const size_t n_rows,
               const size_t n_cols,
               const value_type value)
{ coeffs.resize(n_rows, Vector(n_cols, value)); }


void
Matrix::reinit(const size_t n_rows,
               const size_t n_cols)
{ resize(n_rows, n_cols); }


void
Matrix::reinit(const size_t n_rows,
               const size_t n_cols,
               const value_type value)
{ resize(n_rows, n_cols, value); }


void
Matrix::swap_row(const size_t i, const size_t k)
{ coeffs.at(i).swap(coeffs.at(k)); }


void
Matrix::swap_column(const size_t j, const size_t k)
{
  Assert(j < n_cols() && k < n_cols(), "Out of range error.");
  for (size_t i = 0; i < n_rows(); ++i)
    std::swap(coeffs[i][j], coeffs[i][k]);
}


void
Matrix::swap(Matrix& other)
{ coeffs.swap(other.coeffs); }


void
Matrix::set_diagonal(const Vector& diag)
{
  if (coeffs.empty())
  {
    resize(diag.size(), diag.size());
    for (size_t i = 0; i < diag.size(); ++i)
      coeffs[i][i] = diag[i];
  }
  else
  {
    size_t min_dim = std::min(n_rows(), n_cols());
    Assert(diag.size() == min_dim, "Dimension mismatch error.");

    for (size_t i = 0; i < min_dim; ++i)
      coeffs[i][i] = diag[i];
  }
}


void
Matrix::set_diagonal(const value_type value)
{
  Assert(!coeffs.empty(), "Cannot set an empty matrix with a scalar.");

  size_t min_dim = std::min(n_rows(), n_cols());
  for (size_t i = 0; i < min_dim; ++i)
    coeffs[i][i] = value;
}

//################################################## Scalar Operations

Matrix&
Matrix::operator-()
{
  for (auto& row : coeffs)
    row = -row;
  return *this;
}


Matrix
Matrix::operator-() const
{ return -Matrix(*this); }


Matrix&
Matrix::operator*=(const double value)
{
  for (auto& row : coeffs)
    row *= value;
  return *this;
}


Matrix&
Matrix::operator/=(const value_type value)
{
  {
    Assert(value != 0.0, "Division by zero error.");

    for (auto& row : *this)
      row /= value;
    return *this;
  }
}

//################################################## Linear Algebra

void
Matrix::add(const Matrix& B, const value_type factor)
{
  Assert(n_rows() == B.n_rows(), "Dimension mismatch error.");
  Assert(n_cols() == B.n_cols(), "Dimension mismatch error.");

  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    const double* b_ij = B.data(i);
    for (size_t j = 0; j < n_cols(); ++j)
      *a_ij++ += factor * *b_ij++;
  }
}


void
Matrix::Tadd(const Matrix& B, const value_type factor)
{
  Assert(n_rows() == B.n_cols(), "Dimension mismatch error.");
  Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

  for (size_t i = 0; i < n_rows(); ++i)
  {
    value_type* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j)
      *a_ij++ += factor * B[j][i];
  }
}


Matrix&
Matrix::operator+=(const Matrix& B)
{ add(B); return *this; }


Matrix&
Matrix::operator-=(const Matrix& B)
{ add(B, -1.0); return *this; }


void
Matrix::mmult(const Matrix& B, Matrix& C,
              const bool adding) const
{
  Assert(C.n_rows() == n_rows(), "Dimension mismatch error.");
  Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.");
  Assert(n_cols() == B.n_rows(), "Dimension mismatch error.");

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    value_type* c_ij = C.data(i);
    const value_type* a_i = data(i);
    for (size_t j = 0; j < C.n_cols(); ++j)
    {
      value_type value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_cols(); ++k)
        value += a_i[k] * B(k, j);
     *c_ij++ = value;
    }
  }
}


Matrix
Matrix::mmult(const Matrix& B) const
{
  Matrix C(n_rows(), B.n_cols());
  mmult(B, C);
  return C;
}


void
Matrix::Tmmult(const Matrix& B, Matrix& C,
               const bool adding) const

{
  Assert(C.n_rows() == n_cols(), "Dimension mismatch error.");
  Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.");
  Assert(n_rows() == B.n_rows(), "Dimension mismatch error.");

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    value_type* c_ij = C.data(i);
    for (size_t j = 0; j < C.n_cols(); ++j)
    {
      value_type value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_rows(); ++k)
        value += coeffs[k][i] * B(k, j);
      *c_ij++ = value;
    }
  }
}


Matrix
Matrix::Tmmult(const Matrix& B) const
{
  Matrix C(n_cols(), B.n_cols());
  Tmmult(B, C);
  return C;
}


void
Matrix::mTmult(const Matrix& B, Matrix& C,
               const bool adding) const

{
  Assert(C.n_rows() == n_rows(), "Dimension mismatch error.");
  Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.");
  Assert(n_cols() == B.n_cols(), "Dimension mismatch error.");

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    value_type* c_ij = C.data(i);
    const value_type* a_i = data(i);
    for (size_t j = 0; j < C.n_cols(); ++j)
    {
      const value_type* b_j = B.data(j);
      value_type value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_cols(); ++k)
        value += a_i[k] * *b_j++;
      *c_ij++ = value;
    }
  }
}


Matrix
Matrix::mTmult(const Matrix& B) const
{
  Matrix C(n_rows(), B.n_rows());
  mTmult(B, C);
  return C;
}


void
Matrix::TTmult(const Matrix& B, Matrix& C,
               const bool adding) const

{
  Assert(C.n_rows() == n_cols(), "Dimension mismatch error.");
  Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.");
  Assert(n_rows() == B.n_cols(), "Dimension mismatch error.");

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    value_type* c_ij = C.data(i);
    for (size_t j = 0; j < C.n_cols(); ++j)
    {
      const value_type* b_j = B.data(j);
      value_type value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_rows(); ++k)
        value += coeffs[k][i] * *b_j++;
     *c_ij++ = value;
    }
  }
}


Matrix
Matrix::TTmult(const Matrix& B) const
{
  Matrix C(n_cols(), B.n_rows());
  TTmult(B, C);
  return C;
}


void
Matrix::vmult(const Vector& x, Vector& y,
              const bool adding) const

{
  Assert(x.size() == n_cols(), "Dimension mismatch error.");
  Assert(y.size() == n_rows(), "Dimension mismatch error.");

  value_type* y_i = y.data();
  for (size_t i = 0; i < n_rows(); ++i)
  {
    value_type v = adding? *y_i : 0.0;
    const value_type* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j)
      v += *a_ij++ * x[j];
    *y_i++ = v;
  }
}


Vector
Matrix::vmult(const Vector& x) const
{
  Vector y(n_rows());
  vmult(x, y);
  return y;
}


void
Matrix::vmult_add(const Vector& x, Vector& y) const
{ vmult(x, y, true); }


void
Matrix::Tvmult(const Vector& x, Vector& y,
               const bool adding) const
{
  Assert(x.size() == n_rows(), "Dimension mismatch error.");
  Assert(y.size() == n_cols(), "Dimension mismatch error.");

  if (!adding) y = 0.0;
  for (size_t i = 0; i < n_rows(); ++i)
  {
    const value_type x_i = x[i];
    const value_type* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j)
      y[j] += *a_ij++ * x_i;
  }
}


Vector
Matrix::Tvmult(const Vector& x) const
{
  Vector y(n_cols());
  Tvmult(x, y);
  return y;
}


void
Matrix::Tvmult_add(const Vector& x, Vector& y)
{ Tvmult(x, y, true); }


Vector
Matrix::operator*(const Vector& x) const
{ return vmult(x); }


Matrix
Matrix::transpose() const
{
  Matrix A_T(n_cols(), n_rows());
  for (size_t i = 0; i < n_rows(); ++i)
  {
    const value_type* a_ij = coeffs[i].data();
    for (size_t j = 0; j < n_cols(); ++j)
      A_T[j][i] = *a_ij++;
  }
  return A_T;
}

//################################################## Print Utilities

void
Matrix::print(std::ostream& os,
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
    w = (!width)? precision + 10 : w;
  }
  else
  {
    os.setf(std::ios::fixed, std::ios::floatfield);
    w = (!width)? precision + 5 : w;
  }

  for (uint64_t i = 0; i < n_rows(); ++i)
  {
    const value_type* a_ij = coeffs[i].data();
    for (uint64_t j = 0; j < n_cols(); ++j)
      os << std::setw(w) << *a_ij++;
    os << std::endl;
  }
  os << std::endl;
  os.flags(old_flags);
  os.precision(old_precision);
}


std::string
Matrix::str(const bool scientific,
            const unsigned int precision,
            const unsigned int width) const
{
  std::stringstream ss;
  print(ss, scientific, precision, width);
  return ss.str();
}

//################################################## Validations

bool
Matrix::valid_dimensions(const STLMatrix& A)
{
  size_t m = A.front().size();
  for (const auto& row : A)
    if (row.size() != m) return false;
  return true;
}

//################################################## Methods

Matrix
pdes::Math::operator+(const Matrix& A, const Matrix& B)
{ return Matrix(A) += B; }


Matrix
pdes::Math::operator-(const Matrix& A, const Matrix& B)
{ return Matrix(A) -= B; }


Matrix
pdes::Math::operator*(const Matrix& A, const Matrix& B)
{ return A.mmult(B); }


void
pdes::Math::mmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.mmult(B, C); }


Matrix
pdes::Math::mmult(const Matrix& A, const Matrix& B)
{ return A.mmult(B); }


void
pdes::Math::Tmmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.Tmmult(B, C); }


Matrix
pdes::Math::Tmmult(const Matrix& A, const Matrix& B)
{ return A.Tmmult(B); }


void
pdes::Math::mTmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.mTmult(B, C); }


Matrix
pdes::Math::mTmult(const Matrix& A, const Matrix& B)
{ return A.mTmult(B); }


void
pdes::Math::TTmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.TTmult(B, C); }


Matrix
pdes::Math::TTmult(const Matrix& A, const Matrix& B)
{ return A.TTmult(B); }


void
pdes::Math::vmult(const Matrix& A, const Vector& x, Vector& y)
{ A.vmult(x, y); }


Vector
pdes::Math::vmult(const Matrix& A, const Vector& x)
{ return A.vmult(x); }


void
pdes::Math::Tvmult(const Matrix& A, const Vector& x, Vector& y)
{ A.Tvmult(x, y); }


Vector
pdes::Math::Tvmult(const Matrix& A, const Vector& x)
{ return A.Tvmult(x); }


std::ostream&
pdes::Math::operator<<(std::ostream& os, const Matrix& A)
{ return os << A.str(); }

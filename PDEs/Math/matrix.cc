#include "matrix.h"
#include "vector.h"
#include "macros.h"

#include <cmath>
#include <iomanip>


using namespace Math;


//################################################## Constructors


Matrix::Matrix(const size_t n_rows, const size_t n_cols) :
  vals(n_rows, Vector(n_cols))
{}


Matrix::Matrix(const size_t n_rows,
               const size_t n_cols,
               const double value) :
  vals(n_rows, Vector(n_cols, value))
{}


Matrix::Matrix(const STLMatrix& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.")

  for (const auto& row : other)
    vals.push_back(row);
}


Matrix::Matrix(STLMatrix&& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.")

  for (auto& row : other)
    vals.push_back(row);
}


Matrix::Matrix(const InitializerMatrix list)
{
  for (auto& row : list)
  {
    if (!vals.empty())
    Assert(row.size() == vals.front().size(),
           "Invalid input. All rows must be the same length.")

    vals.push_back(row);
  }
}


//################################################## Assignment Operators


Matrix&
Matrix::operator=(const STLMatrix& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.")

  for (const auto& row : other)
    vals.push_back(row);
  return *this;
}


Matrix&
Matrix::operator=(STLMatrix&& other)
{
  Assert(valid_dimensions(other),
         "Invalid input. All rows must be the same length.")

  for (auto& row : other)
    vals.push_back(row);
  return *this;
}


Matrix&
Matrix::operator=(const double value)
{
  Assert(!vals.empty(), "Cannot set an empty matrix to a scalar.")

  for (auto& row : vals)
    row = value;
  return *this;
}


//################################################## Comparison Operators


bool
Matrix::operator==(const Matrix& other) const
{ return (vals == other.vals); }


bool
Matrix::operator!=(const Matrix& other) const
{ return (vals != other.vals); }


//################################################## Characteristics


size_t
Matrix::n_rows() const
{ return vals.size(); }


size_t
Matrix::n_cols() const
{ return vals.front().size(); }


size_t
Matrix::size() const
{ return n_rows() * n_cols(); }


size_t
Matrix::nnz() const
{
  size_t count = 0;
  for (const auto& row : vals)
    for (const auto& el : row)
      if (el != 0.0) ++count;
  return count;
}


bool
Matrix::empty() const
{ return vals.empty(); }


bool
Matrix::all_zero() const
{ return (nnz() == 0)? true : false; }


//################################################## Accessors


Vector&
Matrix::operator[](const size_t i)
{ return vals[i]; }


const Vector&
Matrix::operator[](const size_t i) const
{ return vals[i]; }


Vector&
Matrix::operator()(const size_t i)
{ return vals[i]; }


const Vector&
Matrix::operator()(const size_t i) const
{ return vals[i]; }


Vector&
Matrix::at(const size_t i)
{ return vals.at(i); }


const Vector&
Matrix::at(const size_t i) const
{ return vals.at(i); }


double&
Matrix::operator()(const size_t i, const size_t j)
{ return vals[i][j]; }


const double&
Matrix::operator()(const size_t i, const size_t j) const
{ return vals[i][j]; }


double&
Matrix::at(const size_t i, const size_t j)
{ return vals.at(i).at(j); }


const double&
Matrix::at(const size_t i, const size_t j) const
{ return vals.at(i).at(j); }


double&
Matrix::diagonal(const size_t i)
{ return vals.at(i).at(i); }


const double&
Matrix::diagonal(const size_t i) const
{ return vals.at(i).at(i); }


Vector
Matrix::diagonal() const
{
  Vector diag;
  size_t min_dim = std::min(n_rows(), n_cols());
  for (size_t i = 0; i < min_dim; ++i)
    diag.push_back(vals[i][i]);
  return diag;
}


Vector*
Matrix::data()
{ return vals.data(); }


const Vector*
Matrix::data() const
{ return vals.data(); }


double*
Matrix::data(const size_t i)
{ return vals.at(i).data(); }


const double*
Matrix::data(const size_t i) const
{ return vals.at(i).data(); }


Matrix::iterator
Matrix::begin()
{ return vals.begin(); }


Matrix::iterator
Matrix::end()
{ return vals.end(); }


Matrix::const_iterator
Matrix::begin() const
{ return vals.begin(); }


Matrix::const_iterator
Matrix::end() const
{ return vals.end(); }


Vector::iterator
Matrix::begin(const size_t i)
{ return vals.at(i).begin(); }


Vector::iterator
Matrix::end(const size_t i)
{ return vals.at(i).end(); }


Vector::const_iterator
Matrix::begin(const size_t i) const
{ return vals.at(i).begin(); }


Vector::const_iterator
Matrix::end(const size_t i) const
{ return vals.at(i).end(); }


//################################################## Modifiers


void
Matrix::clear()
{ vals.clear(); }


void
Matrix::pop_back()
{ vals.pop_back(); }


void
Matrix::push_back(const Vector& row)
{
  Assert(row.size() == n_cols(), "Dimension mismatch error.")
  vals.push_back(row);
}


void
Matrix::push_back(Vector&& row)
{
  Assert(row.size() == n_cols(), "Dimension mismatch error.")
  vals.push_back(row);
}


void
Matrix::resize(const size_t n_rows,
               const size_t n_cols)
{ vals.resize(n_rows, Vector(n_cols)); }


void
Matrix::resize(const size_t n_rows,
               const size_t n_cols,
               const double value)
{ vals.resize(n_rows, Vector(n_cols, value)); }


void
Matrix::reinit(const size_t n_rows,
               const size_t n_cols)
{
  clear();
  resize(n_rows, n_cols);
}


void
Matrix::reinit(const size_t n_rows,
               const size_t n_cols,
               const double value)
{
  clear();
  resize(n_rows, n_cols, value);
}


void
Matrix::swap_row(const size_t i, const size_t k)
{ vals.at(i).swap(vals.at(k)); }


void
Matrix::swap_column(const size_t j, const size_t k)
{
  Assert(j < n_cols() && k < n_cols(), "Out of range error.")
  for (size_t i = 0; i < n_rows(); ++i)
    std::swap(vals[i][j], vals[i][k]);
}


void
Matrix::swap(Matrix& other)
{ vals.swap(other.vals); }


void
Matrix::set_diagonal(const Vector& diag)
{
  if (vals.empty())
  {
    resize(diag.size(), diag.size());
    for (size_t i = 0; i < diag.size(); ++i)
      vals[i][i] = diag[i];
  }
  else
  {
    size_t min_dim = std::min(n_rows(), n_cols());
    Assert(diag.size() == min_dim, "Dimension mismatch error.")

    for (size_t i = 0; i < min_dim; ++i)
      vals[i][i] = diag[i];
  }
}


void
Matrix::set_diagonal(const double value)
{
  Assert(!vals.empty(), "Cannot set an empty matrix with a scalar.")

  size_t min_dim = std::min(n_rows(), n_cols());
  for (size_t i = 0; i < min_dim; ++i)
    vals[i][i] = value;
}


//################################################## Linear Algebra


Matrix&
Matrix::scale(const double factor)
{
  for (auto& row : *this)
    row *= factor;
  return *this;
}


Matrix&
Matrix::add(const Matrix& B, const double a)
{
  Assert(n_rows() == B.n_rows(), "Dimension mismatch error.")
  Assert(n_cols() == B.n_cols(), "Dimension mismatch error.")


  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    const double* b_ij = B.data(i);

    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++b_ij)
      *a_ij += a * *b_ij;
  }
  return *this;
}


Matrix&
Matrix::Tadd(const Matrix& B, const double a)
{
  Assert(n_rows() == B.n_cols(), "Dimension mismatch error.")
  Assert(n_cols() == B.n_rows(), "Dimension mismatch error.")

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij)
      *a_ij += a * B(j, i);
  }
  return *this;
}


Matrix&
Matrix::sadd(const double a, const Matrix& B)
{
  Assert(B.n_rows() == n_rows(), "Dimension mismatch error.")
  Assert(B.n_cols() == n_cols(), "Dimension mismatch error.")

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    const double* b_ij = B.data(i);

    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++b_ij)
      *a_ij = a * *a_ij + *b_ij;
  }
  return *this;
}


Matrix&
Matrix::sadd(const double a, const double b, const Matrix& B)
{
  Assert(B.n_rows() == n_rows(), "Dimension mismatch error.")
  Assert(B.n_cols() == n_cols(), "Dimension mismatch error.")

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    const double* b_ij = B.data(i);

    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++b_ij)
      *a_ij = a * *a_ij + b * *b_ij;
  }
  return *this;
}


Matrix&
Matrix::sTadd(const double a, const Matrix& B)
{
  Assert(n_rows() == B.n_cols(), "Dimension mismatch error.")
  Assert(n_cols() == B.n_rows(), "Dimension mismatch error.")

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij)
      *a_ij = a * *a_ij * B(j, i);
  }
  return *this;
}


Matrix&
Matrix::sTadd(const double a, const double b, const Matrix& B)
{
  Assert(n_rows() == B.n_cols(), "Dimension mismatch error.")
  Assert(n_cols() == B.n_rows(), "Dimension mismatch error.")

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij)
      *a_ij = a * *a_ij * b * B(j, i);
  }
  return *this;
}


void
Matrix::mmult(const Matrix& B, Matrix& C,
              const bool adding) const
{
  Assert(C.n_rows() == n_rows(), "Dimension mismatch error.")
  Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.")
  Assert(n_cols() == B.n_rows(), "Dimension mismatch error.")

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);
    const double* a_i = data(i);

    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_cols(); ++k)
        value += a_i[k] * B(k, j);
      *c_ij = value;
    }
  }
}


void
Matrix::Tmmult(const Matrix& B, Matrix& C,
               const bool adding) const
{
  Assert(C.n_rows() == n_cols(), "Dimension mismatch error.")
  Assert(C.n_cols() == B.n_cols(), "Dimension mismatch error.")
  Assert(n_rows() == B.n_rows(), "Dimension mismatch error.")


  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);
    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_rows(); ++k)
        value += vals[k][i] * B(k, j);
      *c_ij = value;
    }
  }
}


void
Matrix::mTmult(const Matrix& B, Matrix& C,
               const bool adding) const
{
  Assert(C.n_rows() == n_rows(), "Dimension mismatch error.")
  Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.")
  Assert(n_cols() == B.n_cols(), "Dimension mismatch error.")

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);
    const double* a_i = data(i);

    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      const double* b_jk = B.data(j);

      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_cols(); ++k, ++b_jk)
        value += a_i[k] * *b_jk;
      *c_ij = value;
    }
  }
}


void
Matrix::TTmult(const Matrix& B, Matrix& C,
               const bool adding) const
{
  Assert(C.n_rows() == n_cols(), "Dimension mismatch error.")
  Assert(C.n_cols() == B.n_rows(), "Dimension mismatch error.")
  Assert(n_rows() == B.n_cols(), "Dimension mismatch error.")

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);

    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      const double* b_jk = B.data(j);

      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_rows(); ++k, ++b_jk)
        value += vals[k][i] * *b_jk;
      *c_ij = value;
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


Matrix
Matrix::Tmmult(const Matrix& B) const
{
  Matrix C(n_cols(), B.n_cols());
  Tmmult(B, C);
  return C;
}


Matrix
Matrix::mTmult(const Matrix& B) const
{
  Matrix C(n_rows(), B.n_rows());
  mTmult(B, C);
  return C;
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
  Assert(x.size() == n_cols(), "Dimension mismatch error.")
  Assert(y.size() == n_rows(), "Dimension mismatch error.")

  double* y_i = y.data();

  for (size_t i = 0; i < n_rows(); ++i, ++y_i)
  {
    const double* a_ij = data(i);
    const double* x_j = x.data();

    double v = adding? *y_i : 0.0;
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++x_j)
      v += *a_ij * *x_j;
    *y_i = v;
  }
}


void
Matrix::Tvmult(const Vector& x, Vector& y,
               const bool adding) const
{
  Assert(x.size() == n_rows(), "Dimension mismatch error.")
  Assert(y.size() == n_cols(), "Dimension mismatch error.")

  if (!adding) y = 0.0;
  for (size_t i = 0; i < n_rows(); ++i)
  {
    const double x_i = x[i];
    const double* a_ij = data(i);
    double* y_j = y.data();

    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++y_j)
      *y_j += *a_ij * x_i;
  }
}


Vector
Matrix::vmult(const Vector& x) const
{
  Vector y(n_rows());
  vmult(x, y);
  return y;
}


Vector
Matrix::Tvmult(const Vector& x) const
{
  Vector y(n_cols());
  Tvmult(x, y);
  return y;
}


void
Matrix::vmult_add(const Vector& x, Vector& y) const
{ vmult(x, y, true); }


void
Matrix::Tvmult_add(const Vector& x, Vector& y)
{ Tvmult(x, y, true); }


Matrix&
Matrix::operator-()
{ return scale(-1.0); }


Matrix
Matrix::operator-() const
{ return -Matrix(*this); }


Matrix&
Matrix::operator*=(const double factor)
{ return scale(factor); }


Matrix&
Matrix::operator/=(const double factor)
{
  Assert(factor != 0.0, "Zero division error!")
  return scale(1.0 / factor);
}


Matrix&
Matrix::operator+=(const Matrix& B)
{ return add(B); }


Matrix&
Matrix::operator-=(const Matrix& B)
{ return add(B, -1.0); }


Vector
Matrix::operator*(const Vector& x) const
{ return vmult(x); }


Matrix
Matrix::transpose() const
{
  Matrix A_T(n_cols(), n_rows());
  for (size_t i = 0; i < n_rows(); ++i)
  {
    const double* a_ij = vals[i].data();
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

  for (uint64_t i = 0; i < n_rows(); ++i)
  {
    const double* a_ij = vals[i].data();
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
Math::operator+(const Matrix& A, const Matrix& B)
{ return Matrix(A) += B; }


Matrix
Math::operator-(const Matrix& A, const Matrix& B)
{ return Matrix(A) -= B; }


Matrix
Math::operator*(const Matrix& A, const Matrix& B)
{ return A.mmult(B); }


void
Math::mmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.mmult(B, C); }


void
Math::Tmmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.Tmmult(B, C); }


void
Math::mTmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.mTmult(B, C); }


void
Math::TTmult(const Matrix& A, const Matrix& B, Matrix& C)
{ A.TTmult(B, C); }


Matrix
Math::mmult(const Matrix& A, const Matrix& B)
{ return A.mmult(B); }


Matrix
Math::Tmmult(const Matrix& A, const Matrix& B)
{ return A.Tmmult(B); }


Matrix
Math::mTmult(const Matrix& A, const Matrix& B)
{ return A.mTmult(B); }


Matrix
Math::TTmult(const Matrix& A, const Matrix& B)
{ return A.TTmult(B); }


void
Math::vmult(const Matrix& A, const Vector& x, Vector& y)
{ A.vmult(x, y); }


void
Math::Tvmult(const Matrix& A, const Vector& x, Vector& y)
{ A.Tvmult(x, y); }


Vector
Math::vmult(const Matrix& A, const Vector& x)
{ return A.vmult(x); }


Vector
Math::Tvmult(const Matrix& A, const Vector& x)
{ return A.Tvmult(x); }


std::ostream&
Math::operator<<(std::ostream& os, const Matrix& A)
{ return os << A.str(); }

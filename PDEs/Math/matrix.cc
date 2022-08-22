#include "matrix.h"
#include "vector.h"

#include <cmath>
#include <cassert>
#include <iomanip>


using namespace Math;


Matrix::
Matrix(const std::initializer_list<std::initializer_list<double>>& list) :
  values(list.begin(), list.end())
{}


template<typename InputIterator>
Matrix::Matrix(const InputIterator first, const InputIterator last)
{
  reinit(std::distance(first, last), (*first).size());
  std::copy(first, last, this->begin());
  for (const auto& row : values)
    assert(row.size() == this->n_cols());
}
template Matrix::Matrix(const Vector*, const Vector*);


Matrix::Matrix(const size_t n_rows,
               const size_t n_cols,
               const double value) :
  values(n_rows, Vector(n_cols, value))
{}


Matrix::Matrix(const size_t n_rows,
               const size_t n_cols,
               const double* value_ptr)
{
  reinit(n_rows, n_cols);
  for (size_t i = 0; i < n_rows; ++i)
    for (size_t j = 0; j < n_cols; ++j, ++value_ptr)
      values[i][j] = *value_ptr;
}


Matrix::Matrix(const Vector& diagonal)
{
  size_t n = diagonal.size();
  reinit(n, n);
  for (size_t i = 0; i < n; ++i)
    values[i][i] = diagonal[i];
}


Matrix::Matrix(Vector&& diagonal)
{
  size_t n = diagonal.size();
  reinit(n, n);
  for (size_t i = 0; i < n; ++i)
    values[i][i] = std::move(diagonal[i]);
}


Matrix::Matrix(const std::initializer_list<double>& diagonal)
{
  *this = Matrix(Vector(diagonal));
}


void
Matrix::reinit(const size_t n_rows,
               const size_t n_cols,
               const double value)
{
  values.clear();
  values.resize(n_rows, Vector(n_cols, value));
}


Matrix&
Matrix::operator=(const Matrix& other)
{
  values = other.values;
  return *this;
}


Matrix&
Matrix::operator=(Matrix&& other)
{
  values = std::move(other.values);
  return *this;
}


Matrix&
Matrix::operator=(const double value)
{
  if (empty())
    reinit(1, 1, value);
  else
    for (auto& row : values)
      row = value; // full row assignment to a scalar
  return *this;
}


size_t
Matrix::n_rows() const
{
  return values.size();
}


size_t
Matrix::n_cols() const
{
  return values.front().size();
}


size_t
Matrix::size() const
{
  return n_rows() * n_cols();
}


size_t
Matrix::n_nonzero_elements() const
{
  size_t count = 0;
  for (const auto& row : values)
    count += row.n_nonzero_elements();
  return count;
}


bool
Matrix::empty() const
{
  return values.empty();
}


Matrix
Matrix::transpose() const
{
  Matrix A_T(n_cols(), n_rows());
  for (size_t i = 0; i < n_rows(); ++i)
  {
    const double* a_ij = values[i].data();
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij)
      A_T[j][i] = *a_ij;
  }
  return A_T;
}


bool
Matrix::operator==(const Matrix& other) const
{
  return (values == other.values);
}


bool
Matrix::operator!=(const Matrix& other) const
{
  return (values != other.values);
}



Vector&
Matrix::operator[](const size_t i)
{
  return values[i];
}


const Vector&
Matrix::operator[](const size_t i) const
{
  return values[i];
}


Vector&
Matrix::operator()(const size_t i)
{
  return values[i];
}


const Vector&
Matrix::operator()(const size_t i) const
{
  return values[i];
}


Vector&
Matrix::at(const size_t i)
{
  return values.at(i);
}


const Vector&
Matrix::at(const size_t i) const
{
  return values.at(i);
}


double&
Matrix::operator()(const size_t i, const size_t j)
{
  return values[i][j];
}


const double&
Matrix::operator()(const size_t i, const size_t j) const
{
  return values[i][j];
}


double&
Matrix::at(const size_t i, const size_t j)
{
  return values.at(i).at(j);
}


const double&
Matrix::at(const size_t i, const size_t j) const
{
  return values.at(i).at(j);
}


double&
Matrix::diag(const size_t i)
{
  return values.at(i).at(i);
}


const double&
Matrix::diag(const size_t i) const
{
  return values.at(i).at(i);
}


Vector*
Matrix::data()
{
  return values.data();
}


const Vector*
Matrix::data() const
{
  return values.data();
}


double*
Matrix::data(const size_t i)
{
  return values.at(i).data();
}


const double*
Matrix::data(const size_t i) const
{
  return values.at(i).data();
}


Matrix::iterator
Matrix::begin()
{
  return values.begin();
}


Matrix::iterator
Matrix::end()
{
  return values.end();
}


Matrix::const_iterator
Matrix::begin() const
{
  return values.begin();
}


Matrix::const_iterator
Matrix::end() const
{
  return values.end();
}


Vector::iterator
Matrix::begin(const size_t i)
{
  return values.at(i).begin();
}


Vector::iterator
Matrix::end(const size_t i)
{
  return values.at(i).end();
}


Vector::const_iterator
Matrix::begin(const size_t i) const
{
  return values.at(i).begin();
}


Vector::const_iterator
Matrix::end(const size_t i) const
{
  return values.at(i).end();
}


void
Matrix::clear()
{
  values.clear();
}


void
Matrix::push_back(const Vector& row)
{
  assert(row.size() == n_cols());
  values.push_back(row);
}


void
Matrix::push_back(Vector&& row)
{
  assert(row.size() == n_cols());
  values.push_back(row);
}


void
Matrix::pop_back()
{
  values.pop_back();
}


void
Matrix::resize(const size_t n_rows,
               const size_t n_cols,
               const double value)
{
  values.resize(n_rows, Vector(n_cols, value));
}


void
Matrix::swap_row(const size_t i, const size_t k)
{
  values.at(i).swap(values.at(k));
}


void
Matrix::swap_column(const size_t j, const size_t k)
{
  assert(j < n_cols() && k < n_cols());
  for (size_t i = 0; i < n_rows(); ++i)
    std::swap(values[i][j], values[i][k]);
}


void
Matrix::swap(Matrix& other)
{
  values.swap(other.values);
}


void
Matrix::set_diag(const Vector& diagonal)
{
  if (empty())
  {
    this->reinit(diagonal.size(), diagonal.size());
    for (size_t i = 0; i < diagonal.size(); ++i)
      values[i][i] = diagonal[i];
  }
  else
  {
    size_t min_dim = std::min(n_rows(), n_cols());
    assert(diagonal.size() == min_dim);
    for (size_t i = 0; i < min_dim; ++i)
      values[i][i] = diagonal[i];
  }
}


void
Matrix::set_diag(Vector&& diagonal)
{
  if (empty())
  {
    size_t n = diagonal.size();
    reinit(n, n);
    for (size_t i = 0; i < n; ++i)
      values[i][i] = std::move(diagonal[i]);
  }
  else
  {
    size_t min_dim = std::min(n_rows(), n_cols());
    assert(diagonal.size() == min_dim);
    for (size_t i = 0; i < min_dim; ++i)
      values[i][i] = std::move(diagonal[i]);
  }
}


void
Matrix::set_diag(const std::initializer_list<double>& diagonal)
{
  this->set_diag(Vector(diagonal));
}


void
Matrix::set_diag(const double value)
{
  if (empty())
    reinit(1, 1, value);
  else
  {
    size_t min_dim = std::min(n_rows(), n_cols());
    for (size_t i = 0; i < min_dim; ++i)
      values[i][i] = value;
  }
}


Matrix&
Matrix::scale(const double factor)
{
  for (auto& row : values)
    row *= factor;
  return *this;
}


Matrix&
Matrix::operator-()
{
  return scale(-1.0);
}


Matrix
Matrix::operator-() const
{
  return Matrix(*this).scale(-1.0);
}


Matrix&
Matrix::operator*=(const double factor)
{
  return scale(factor);
}


Matrix
Matrix::operator*(const double factor) const
{
  return Matrix(*this).scale(factor);
}


Matrix&
Matrix::operator/=(const double factor)
{
  assert(factor != 0.0);
  return scale(1.0/factor);
}


Matrix
Matrix::operator/(const double factor) const
{
  assert(factor != 0.0);
  return Matrix(*this).scale(1.0/factor);
}


Matrix&
Matrix::add(const double b, const Matrix& B)
{
  assert(B.n_rows() == n_rows());
  assert(B.n_cols() == n_cols());

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    const double* b_ij = B.data(i);
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++b_ij)
      *a_ij += b * *b_ij;
  }
  return *this;
}


Matrix&
Matrix::operator+=(const Matrix& B)
{
  return add(1.0, B);
}


Matrix
Matrix::operator+(const Matrix& B) const
{
  return Matrix(*this).add(1.0, B);
}


Matrix&
Matrix::operator-=(const Matrix& B)
{
  return add(-1.0, B);
}


Matrix
Matrix::operator-(const Matrix& B) const
{
  return Matrix(*this).add(-1.0, B);
}


Matrix&
Matrix::Tadd(const double b, const Matrix& B)
{
  assert(B.n_rows() == n_cols());
  assert(B.n_cols() == n_rows());

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij)
      *a_ij += b * B(j, i);
  }
  return *this;
}


Matrix&
Matrix::sadd(const double a, const double b, const Matrix& B)
{
  assert(B.n_rows() == n_rows());
  assert(B.n_cols() == n_cols());

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
Matrix::sadd(const double a, const Matrix& B)
{
  return sadd(a, 1.0, B);
}


Matrix&
Matrix::sTadd(const double a, const double b, const Matrix& B)
{
  assert(B.n_rows() == n_cols());
  assert(B.n_cols() == n_rows());

  // Perform the operation
  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = data(i);

    for (size_t j = 0; j < n_cols(); ++j, ++a_ij)
      *a_ij = a**a_ij*b*B(j, i);
  }
  return *this;
}



Matrix&
Matrix::sTadd(const double a, const Matrix& B)
{
  return sTadd(a, 1.0, B);
}


void
Matrix::mmult(const Matrix& B, Matrix& C, const bool adding) const
{
  assert(C.n_rows() == n_rows());
  assert(C.n_cols() == B.n_cols());
  assert(B.n_rows() == n_cols());

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


Matrix
Matrix::mmult(const Matrix& B) const
{
  Matrix C(n_rows(), B.n_cols());
  mmult(B, C);
  return C;
}


void
Matrix::Tmmult(const Matrix& B, Matrix& C, const bool adding) const
{
  assert(C.n_rows() == n_cols());
  assert(C.n_cols() == B.n_cols());
  assert(B.n_rows() == n_rows());


  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);

    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_rows(); ++k)
        value += values[k][i]*B(k, j);
      *c_ij = value;
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
Matrix::mTmult(const Matrix& B, Matrix& C, const bool adding) const
{
  assert(C.n_rows() == n_rows());
  assert(C.n_cols() == B.n_rows());
  assert(B.n_cols() == n_cols());

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);
    const double* a_i = data(i);

    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      const double* b_jk = B.data(j);

      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_cols(); ++k, ++b_jk)
        value += a_i[k]**b_jk;
      *c_ij = value;
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
  assert(C.n_rows() == n_cols());
  assert(C.n_cols() == B.n_rows());
  assert(B.n_cols() == n_rows());

  for (size_t i = 0; i < C.n_rows(); ++i)
  {
    double* c_ij = C.data(i);

    for (size_t j = 0; j < C.n_cols(); ++j, ++c_ij)
    {
      const double* b_jk = B.data(j);

      double value = adding? *c_ij : 0.0;
      for (size_t k = 0; k < n_rows(); ++k, ++b_jk)
        value += values[k][i]**b_jk;
      *c_ij = value;
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
  assert(x.size() == n_cols());
  assert(y.size() == n_rows());

  double* y_i = y.data();
  for (size_t i = 0; i < n_rows(); ++i, ++y_i)
  {
    const double* a_ij = data(i);
    const double* x_j = x.data();

    double v = adding? *y_i : 0.0;
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++x_j)
      v += *a_ij**x_j;
    *y_i = v;
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
{
  vmult(x, y, true);
}


Vector
Matrix::operator*(const Vector& x) const
{
  return vmult(x);
}



void
Matrix::Tvmult(const Vector& x, Vector& y,
               const bool adding) const
{
  assert(x.size() == n_rows());
  assert(y.size() == n_cols());

  if (!adding) y = 0.0;
  for (size_t i = 0; i < n_rows(); ++i)
  {
    const double x_i = x[i];
    const double* a_ij = data(i);

    double* y_j = y.data();
    for (size_t j = 0; j < n_cols(); ++j, ++a_ij, ++y_j)
      *y_j += *a_ij*x_i;
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
{
  Tvmult(x, y, true);
}


std::string
Matrix::str(const bool scientific,
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

  for (size_t i = 0; i < n_rows(); ++i)
  {
    const double* a_ij = values[i].data();
    for (uint64_t j = 0; j < n_cols(); ++j)
      ss << std::setw(w) << *a_ij++;
    ss << std::endl;
  }
  ss << std::endl;
  return ss.str();
}



void
Matrix::print(std::ostream& os,
              const bool scientific,
              const unsigned int precision,
              const unsigned int width) const
{
  os << str(scientific, precision, width);
}


void
Math::mmult(const Matrix& A, const Matrix& B, Matrix& C)
{
  A.mmult(B, C);
}


Matrix
Math::mmult(const Matrix& A, const Matrix& B)
{
  return A.mmult(B);
}


void
Math::Tmmult(const Matrix& A, const Matrix& B, Matrix& C)
{
  A.Tmmult(B, C);
}


Matrix
Math::Tmmult(const Matrix& A, const Matrix& B)
{
  return A.Tmmult(B);
}


void
Math::mTmult(const Matrix& A, const Matrix& B, Matrix& C)
{
  A.mTmult(B, C);
}


Matrix
Math::mTmult(const Matrix& A, const Matrix& B)
{
  return A.mTmult(B);
}


void
Math::TTmult(const Matrix& A, const Matrix& B, Matrix& C)
{
  A.TTmult(B, C);
}



Matrix
Math::TTmult(const Matrix& A, const Matrix& B)
{
  return A.TTmult(B);
}


void
Math::vmult(const Matrix& A, const Vector& x, Vector& y)
{
  A.vmult(x, y);
}


Vector
Math::vmult(const Matrix& A, const Vector& x)
{
  return A.vmult(x);
}


void
Math::Tvmult(const Matrix& A, const Vector& x, Vector& y)
{
  A.Tvmult(x, y);
}


Vector
Math::Tvmult(const Matrix& A, const Vector& x)
{
  return A.Tvmult(x);
}


std::ostream&
Math::operator<<(std::ostream& os, const Matrix& A)
{
  return os << A.str();
}

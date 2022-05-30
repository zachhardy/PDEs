#include "sparse_matrix.h"
#include "matrix.h"
#include "vector.h"
#include "macros.h"

#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace pdes::Math;

//################################################## Constructors


SparseMatrix::SparseMatrix() :
    rows(0), cols(0), colnums(), coeffs()
{}


SparseMatrix::SparseMatrix(const size_t n_rows,
                           const size_t n_cols,
                           const size_t default_row_length) :
    rows(n_rows), cols(n_cols),
    colnums(n_rows, std::vector<size_t>(default_row_length)),
    coeffs(n_rows, std::vector<double>(default_row_length))
{}


SparseMatrix::SparseMatrix(const size_t n,
                           const size_t default_row_length) :
    SparseMatrix(n, n, default_row_length)
{}


SparseMatrix::SparseMatrix(SparsityPattern sparsity_pattern) :
    rows(sparsity_pattern.size()), colnums(sparsity_pattern),
    coeffs(sparsity_pattern.size())
{
  cols = 0;
  {
    cols = 0;
    for (size_t i = 0; i < rows; ++i)
    {
      // Sort column indices and allocate data
      std::sort(colnums[i].begin(), colnums[i].end());
      coeffs[i].resize(colnums[i].size());

      // Set number of columns
      for (const auto& col : colnums[i])
        cols = (col > cols)? col : cols;
    }
  }
}


SparseMatrix::SparseMatrix(const Matrix& other) :
    rows(other.n_rows()), cols(other.n_cols()),
    colnums(other.n_rows()), coeffs(other.n_rows())
{
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      if (other[i][j] != 0.0)
      {
        colnums[i].push_back(j);
        coeffs[i].push_back(other[i][j]);
      }
}


//################################################## Assignment


SparseMatrix&
SparseMatrix::operator=(const Matrix& other)

{
  colnums.clear();
  coeffs.clear();

  rows = other.n_rows();
  cols = other.n_cols();
  colnums.resize(rows);
  coeffs.resize(rows);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      if (other[i][j] != 0.0)
      {
        colnums[i].push_back(j);
        coeffs[i].push_back(other[i][j]);
      }
  return *this;
}


SparseMatrix&
SparseMatrix::operator=(const double value)
{
  Assert(!empty(), "Cannot set an empty matrix to a scalar.");
  for (entry el : *this)
    el.value = value;
  return *this;
}


//################################################## Comparison


bool
SparseMatrix::operator==(const SparseMatrix& other) const
{
  return (rows    == other.rows &&
          cols    == other.cols &&
          colnums == other.colnums &&
          coeffs == other.coeffs);
}


bool
SparseMatrix::operator!=(const SparseMatrix& other) const
{ return !(*this == other); }


//################################################## Characteristics


size_t
SparseMatrix::n_rows() const
{ return rows; }


size_t
SparseMatrix::n_cols() const
{ return cols; }


size_t
SparseMatrix::nnz() const
{
  size_t count = 0;
  for (size_t i = 0; i < rows; ++i)
    count += row_length(i);
  return count;
}


size_t
SparseMatrix::row_length(const size_t i) const
{ return colnums.at(i).size(); }


bool
SparseMatrix::empty() const
{
  return (rows == 0 && cols == 0 &&
          colnums.empty() && coeffs.empty());
}


//################################################## Accessors


const size_t&
SparseMatrix::column(const size_t i, const size_t jr) const
{ return colnums.at(i).at(jr); }


double&
SparseMatrix::value(const size_t i, const size_t jr)
{ return coeffs.at(i).at(jr); }


const double&
SparseMatrix::value(const size_t i, const size_t jr) const
{ return coeffs.at(i).at(jr); }


double*
SparseMatrix::locate(const size_t i, const size_t j)
{
  Assert(i < rows && j < cols, "Out of range error.");
  for (size_t jr = 0; jr < row_length(i); ++jr)
    if (colnums[i][jr] == j)
      return &coeffs[i][jr];
  return nullptr;
}


const double*
SparseMatrix::locate(const size_t i, const size_t j) const
{
  Assert(i < rows && j < cols, "Out of range error.");
  for (uint64_t jr = 0; jr < row_length(i); ++jr)
    if (colnums[i][jr] == j)
      return &coeffs[i][jr];
  return nullptr;
}


double&
SparseMatrix::operator()(const size_t i, const size_t j)
{
  double* value_ptr = locate(i, j);
  Assert(value_ptr != nullptr,
         "Cannot access an uninitialized element.");

  return *value_ptr;
}


const double&
SparseMatrix::operator()(const size_t i, const size_t j) const
{
  const double* value_ptr = locate(i, j);
  Assert(value_ptr != nullptr,
         "Cannot access an uninitialized element.");

  return *value_ptr;
}


double*
SparseMatrix::diagonal(const size_t i)
{ return locate(i, i); }


const double*
SparseMatrix::diagonal(const size_t i) const
{ return locate(i, i); }


std::vector<double*>
SparseMatrix::diagonal()
{
  size_t min_dim = std::min(rows, cols);
  std::vector<double*> diag(min_dim);
  for (size_t i = 0; i < min_dim; ++i)
    diag[i] = diagonal(i);
  return diag;
}


std::vector<const double*>
SparseMatrix::diagonal() const
{
  size_t min_dim = std::min(rows, cols);
  std::vector<const double*> diag(min_dim);
  for (size_t i = 0; i < min_dim; ++i)
    diag[i] = diagonal(i);
  return diag;
}


SparseMatrix::iterator
SparseMatrix::begin()
{ return (!empty()) ? begin(0) : end(); }


SparseMatrix::iterator
SparseMatrix::end()
{ return {this}; }


SparseMatrix::iterator
SparseMatrix::begin(const size_t i)
{
  Assert(i < rows, "Out of range error.");
  return {this, i};
}


SparseMatrix::iterator
SparseMatrix::end(const size_t i)
{
  Assert(i < rows, "Out of range error.");
  if (i + 1 == rows) return end();
  else return begin(i + 1);
}


SparseMatrix::row_accessor
SparseMatrix::row(const size_t i)
{ return {*this, i}; }


SparseMatrix::const_iterator
SparseMatrix::begin() const
{ return (!empty()) ? begin(0) : end(); }


SparseMatrix::const_iterator
SparseMatrix::end() const
{ return {this}; }


SparseMatrix::const_iterator
SparseMatrix::begin(const size_t i) const
{
  Assert(i < rows, "Out of range error.");
  return {this, i};
}


SparseMatrix::const_iterator
SparseMatrix::end(const size_t i) const
{
  Assert(i < rows, "Out of range error.");
  if (i + 1 == rows) return end();
  else return begin(i + 1);
}


SparseMatrix::const_row_accessor
SparseMatrix::const_row(const size_t i) const
{ return {*this, i}; }


//################################################## Modifiers


void
SparseMatrix::clear()
{
  rows = cols = 0;
  colnums.clear();
  coeffs.clear();
}


void
SparseMatrix::reinit(const size_t n_rows,
                     const size_t n_cols,
                     const size_t default_row_length)
{
  clear();
  rows = n_rows;
  cols = n_cols;
  colnums.resize(n_rows, std::vector<size_t>(default_row_length));
  coeffs.resize(n_rows, std::vector<double>(default_row_length));
}


void
SparseMatrix::reinit(const size_t n,
                     const size_t default_row_length)
{ reinit(n, n, default_row_length); }


void
SparseMatrix::set(const size_t i, const size_t j, const double value)
{
  Assert(i < rows && j < cols, "Out of range error.");
  if (value == 0.0)
    return;

  /* If the row is empty or the column number is larger than all current
   * entries on the row, add to the back of the row. */
  if (colnums[i].size() == 0 or colnums[i].back() < j)
  {
    colnums[i].push_back(j);
    coeffs[i].push_back(value);
    return;
  }

  // Find the index to insert the entry which maintains sorting.
  auto it = std::lower_bound(colnums[i].begin(), colnums[i].end(), j);
  uint64_t jr = it - colnums[i].begin();

  // If this points to an existing column, override the value
  if (*it == j)
  {
    coeffs[i][jr] = value;
    return;
  }

  // Insert the entries into the data structures maintaining sorting
  colnums[i].insert(it, j);
  coeffs[i].insert(coeffs[i].begin() + jr, value);
}


void
SparseMatrix::add(const size_t i, const size_t j, const double value)
{
  Assert(i < rows && j < cols, "Out of range error.");

  // Locate the element
  double* value_ptr = locate(i, j);

  // Set if uninitialized, otherwise, add.
  if (locate(i, j) == nullptr) set(i, j, value);
  else *value_ptr += value;
}


void
SparseMatrix::swap_row(const size_t i, const size_t k)
{
  colnums.at(i).swap(colnums.at(k));
  coeffs.at(i).swap(coeffs.at(k));
}


void
SparseMatrix::swap(SparseMatrix& other)
{
  std::swap(rows, other.rows);
  std::swap(cols, other.cols);
  colnums.swap(other.colnums);
  coeffs.swap(other.coeffs);
}


void
SparseMatrix::eliminate_zeros()
{

  for (size_t i = 0; i < rows; ++i)
  {
    std::vector<std::vector<size_t>::iterator> zero_colnums;
    std::vector<std::vector<double>::iterator> zero_coeffs;

    std::vector<size_t>::iterator col_ptr = colnums[i].begin();
    std::vector<double>::iterator coeff_ptr = coeffs[i].begin();
    std::vector<size_t>::iterator end_ptr = colnums[i].end();

    while (col_ptr != end_ptr)
    {
      if (std::fabs(*coeff_ptr) < 1.0e-14)
      {
        zero_colnums.push_back(col_ptr);
        zero_coeffs.push_back(coeff_ptr);
      }
      ++col_ptr; ++coeff_ptr;
    }

    for (size_t jr = 0; jr < zero_coeffs.size(); ++jr)
    {
      colnums[i].erase(zero_colnums[jr]);
      coeffs[i].erase(zero_coeffs[jr]);
    }
  }
}


//################################################## Linear Algebra


SparseMatrix&
SparseMatrix::scale(const value_type factor)
{
  for (auto el : *this)
    el.value *= factor;
  return *this;
}


SparseMatrix&
SparseMatrix::add(const SparseMatrix& B, const double a)
{
  Assert(colnums == B.colnums,
         "Adding two different sparse matrices with different "
         "sparsity patterns is not allowed.");

  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = coeffs[i].data();
    const double* b_ij = B.coeffs[i].data();

    for (size_t jr = 0; jr < row_length(i); ++jr, ++a_ij, ++b_ij)
      *a_ij += a * *b_ij;
  }
  return *this;
}


SparseMatrix&
SparseMatrix::sadd(const value_type a, const SparseMatrix& B)
{
  Assert(colnums == B.colnums,
         "Adding two different sparse matrices with different "
         "sparsity patterns is not allowed.");

  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = coeffs[i].data();
    const double* b_ij = B.coeffs[i].data();

    for (size_t jr = 0; jr < row_length(i); ++jr, ++a_ij, ++b_ij)
      *a_ij = a * *a_ij + *b_ij;
  }
  return *this;
}


SparseMatrix&
SparseMatrix::
sadd(const value_type a, const value_type b, const SparseMatrix& B)
{
  Assert(colnums == B.colnums,
         "Adding two different sparse matrices with different "
         "sparsity patterns is not allowed.");

  for (size_t i = 0; i < n_rows(); ++i)
  {
    double* a_ij = coeffs[i].data();
    const double* b_ij = B.coeffs[i].data();

    for (size_t jr = 0; jr < row_length(i); ++jr, ++a_ij, ++b_ij)
      *a_ij = a * *a_ij + b * *b_ij;
  }
  return *this;
}


void
SparseMatrix::vmult(const Vector& x, Vector& y, const bool adding) const
{
  Assert(x.size() == cols, "Dimension mismatch error.");
  Assert(y.size() == rows, "Dimension mismatch error.");

  if (!adding) y = 0.0;
  for (const const_entry el : *this)
    y[el.row] += el.value * x[el.column];
}


Vector
SparseMatrix::vmult(const Vector& x) const
{
  Vector y(rows);
  vmult(x, y);
  return y;
}


void
SparseMatrix::vmult_add(const Vector& x, Vector& y) const
{ vmult(x, y, true); }


void
SparseMatrix::Tvmult(const Vector& x, Vector& y, const bool adding) const
{
  Assert(x.size() == rows, "Dimension mismatch error.");
  Assert(y.size() == cols, "Dimension mismatch error.");

  if (!adding) y = 0.0;
  for (const const_entry el : *this)
    y[el.column] += el.value * x[el.row];
}


Vector
SparseMatrix::Tvmult(const Vector& x) const
{
  Vector y(cols);
  Tvmult(x, y);
  return y;
}


void
SparseMatrix::Tvmult_add(const Vector& x, Vector& y) const
{ Tvmult(x, y, true); }


SparseMatrix&
SparseMatrix::operator-()
{ return scale(-1.0); }


SparseMatrix
SparseMatrix::operator-() const
{ return -SparseMatrix(*this); }


SparseMatrix&
SparseMatrix::operator*=(const double factor)
{ return scale(factor); }


SparseMatrix&
SparseMatrix::operator/=(const double factor)
{
  Assert(factor != 0.0, "Zero division error.");
  return scale(1.0/factor);
}


SparseMatrix&
SparseMatrix::operator+=(const SparseMatrix& B)
{ return add(B); }


SparseMatrix&
SparseMatrix::operator-=(const SparseMatrix& B)
{ return add(B, -1.0); }


Vector
SparseMatrix::operator*(const Vector& x) const
{ return vmult(x); }


//################################################## Printing Utilities


void
SparseMatrix::print(std::ostream& os, 
                    const bool scientific, 
                    const unsigned int precision) const

{
  unsigned int w                   = precision + 7;
  std::ios::fmtflags old_flags     = os.flags();
  unsigned int       old_precision = os.precision(precision);

  if (scientific)
    os.setf(std::ios::scientific, std::ios::floatfield);
  else
    os.setf(std::ios::fixed, std::ios::floatfield);

  os.setf(std::ios::left, std::ios::adjustfield);

  os << std::setw(w) << "Row"
     << std::setw(w) << "Column"
     << std::setw(w) << "Value" << std::endl;
  os << std::setw(w) << "---"
      << std::setw(w) << "------"
      << std::setw(w) << "-----" << std::endl;

  for (const auto elem : *this)
    os << std::setw(w) << elem.row
       << std::setw(w) << elem.column
       << std::setw(w)<< elem.value << std::endl;

  os.flags(old_flags);
  os.precision(old_precision);
}


void 
SparseMatrix::print_formatted(std::ostream& os,
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
    w = (!width)? precision + 7 : w;
  }
  else
  {
    os.setf(std::ios::fixed, std::ios::floatfield);
    w = (!width)? precision + 4 : w;
  }

  // Loop over rows and columns and print element-wise
  for (size_t i = 0; i < rows; ++i)
  {
    for (size_t j = 0; j < cols; ++j)
    {
      // Print the entry or a zero
      const double* entry = locate(i, j);
      os << std::setw(w) << ((!entry)? 0.0 : *entry);
    }
    os << std::endl;
  }
  os << std::endl;

  os.flags(old_flags);
  os.precision(old_precision);
}


std::string
SparseMatrix::str(const bool scientific,
                  const unsigned int precision,
                  const unsigned int width) const
{
  std::stringstream ss;
  print_formatted(ss, scientific, precision, width);
  return ss.str();
}


//################################################## Methods


SparseMatrix
pdes::Math::operator*(const double factor, const SparseMatrix& A)
{ return SparseMatrix(A) *= factor; }


SparseMatrix
pdes::Math::operator*(const SparseMatrix& A, const double factor)
{ return SparseMatrix(A) *= factor; }


SparseMatrix
pdes::Math::operator/(const SparseMatrix& A, const double factor)
{ return SparseMatrix(A) /= factor; }


SparseMatrix
pdes::Math::operator+(const SparseMatrix& A, const SparseMatrix& B)
{
  SparseMatrix C(A);
  C.add(B);
  return C;
}


SparseMatrix
pdes::Math::operator-(const SparseMatrix& A, const SparseMatrix& B)
{
  SparseMatrix C(A);
  C.add(B, -1.0);
  return C;
}


void
pdes::Math::vmult(const SparseMatrix& A, const Vector& x, Vector& y)
{ A.vmult(x, y); }


Vector
pdes::Math::vmult(const SparseMatrix& A, const Vector& x)
{ return A.vmult(x); }


void
pdes::Math::Tvmult(const SparseMatrix& A, const Vector& x, Vector& y)
{ A.Tvmult(x, y); }


Vector
pdes::Math::Tvmult(const SparseMatrix& A, const Vector& x)
{ return A.Tvmult(x); }


std::ostream&
pdes::Math::operator<<(std::ostream& os, const SparseMatrix& A)
{ return os << A.str(); }

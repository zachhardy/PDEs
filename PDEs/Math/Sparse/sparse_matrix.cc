#include "sparse_matrix.h"
#include "matrix.h"
#include "vector.h"
#include "macros.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>


using namespace Math;


SparseMatrix::Iterator::Entry::
Entry(const size_t row, const size_t column, double& value) :
  row(row), column(column), value(value)
{}


SparseMatrix::Iterator::
Iterator(SparseMatrix* matrix, const size_t row, const unsigned int index) :
  matrix(matrix), row(row), index(index)
{}


SparseMatrix::Iterator::
Iterator(SparseMatrix* matrix) :
  matrix(matrix), row(-1), index(-1)
{}


SparseMatrix::Iterator&
SparseMatrix::Iterator::operator++()
{
  advance();
  return *this;
}


SparseMatrix::Iterator
SparseMatrix::Iterator::operator++(int)
{
  auto it = *this;
  advance();
  return it;
}


SparseMatrix::Iterator::Entry
SparseMatrix::Iterator::operator*() const
{
  return {row,
          matrix->colnums[row][index],
          matrix->values[row][index]};
}


bool
SparseMatrix::Iterator::operator==(const Iterator& other) const
{
  return (matrix == other.matrix &&
          row == other.row &&
          index == other.index);
}


bool
SparseMatrix::Iterator::operator!=(const Iterator& other) const
{
  return !(*this == other);
}


bool
SparseMatrix::Iterator::operator<(const Iterator& other) const
{
  assert(matrix == other.matrix);

  // Automatically false if at the end
  if (row == -1)
    return false;
  assert(row < matrix->n_rows());

  // Automatically true if other iterator is at the end
  if (other.row == -1)
    return true;
  assert(other.row < other.matrix->n_rows());

  // Otherwise, perform the comparison
  return (row < other.row ||
          (row == other.row && index < other.index));
}


void
SparseMatrix::Iterator::advance()
{
  ++index; // increment the index

  // Handle the end of a row
  if (index == matrix->row_length(row))
  {
    // If more rows exist, move to the next row
    if (row + 1 < matrix->n_rows())
      *this = matrix->begin(row + 1);

    // If no more rows exist, set to invald
    else
      *this = matrix->end();
  }
}


SparseMatrix::ConstIterator::ConstEntry::
ConstEntry(const size_t row, const size_t column, const double& value) :
    row(row), column(column), value(value)
{}


SparseMatrix::ConstIterator::
ConstIterator(const SparseMatrix* matrix,
              const size_t row,
              const unsigned int index) :
    matrix(matrix), row(row), index(index)
{}


SparseMatrix::ConstIterator::
ConstIterator(const SparseMatrix* matrix) :
    matrix(matrix), row(-1), index(-1)
{}


SparseMatrix::ConstIterator&
SparseMatrix::ConstIterator::operator++()
{
  advance();
  return *this;
}


SparseMatrix::ConstIterator
SparseMatrix::ConstIterator::operator++(int)
{
  auto it = *this;
  advance();
  return it;
}


SparseMatrix::ConstIterator::ConstEntry
SparseMatrix::ConstIterator::operator*() const
{
  return {row,
          matrix->colnums[row][index],
          matrix->values[row][index]};
}


bool
SparseMatrix::ConstIterator::
operator==(const ConstIterator& other) const
{
  return (matrix == other.matrix &&
          row == other.row &&
          index == other.index);
}


bool
SparseMatrix::ConstIterator::
operator!=(const ConstIterator& other) const
{
  return !(*this == other);
}


bool
SparseMatrix::ConstIterator::
operator<(const ConstIterator& other) const
{
  assert(matrix == other.matrix);

  // Automatically false if at the end
  if (row == -1)
    return false;
  assert(row < matrix->n_rows());

  // Automatically true if other iterator is at the end
  if (other.row == -1)
    return true;
  assert(other.row < other.matrix->n_rows());

  // Otherwise, perform the comparison
  return (row < other.row ||
          (row == other.row && index < other.index));
}


void
SparseMatrix::ConstIterator::advance()
{
  ++index; // increment the index

  // Handle the end of a row
  if (index == matrix->row_length(row))
  {
    // If more rows exist, move to the next row
    if (row + 1 < matrix->n_rows())
      *this = matrix->begin(row + 1);

      // If no more rows exist, set to invald
    else
      *this = matrix->end();
  }
}


SparseMatrix::RowIterator::
RowIterator(SparseMatrix* matrix, const size_t row) :
  matix(matrix), row(row)
{}


SparseMatrix::Iterator
SparseMatrix::RowIterator::begin()
{
  return matix->begin(row);
}


SparseMatrix::Iterator
SparseMatrix::RowIterator::end()
{
  return matix->end(row);
}


SparseMatrix::ConstRowIterator::
ConstRowIterator(const SparseMatrix* matrix, const size_t row) :
    matix(matrix), row(row)
{}


SparseMatrix::ConstIterator
SparseMatrix::ConstRowIterator::begin()
{
  return matix->begin(row);
}


SparseMatrix::ConstIterator
SparseMatrix::ConstRowIterator::end()
{
  return matix->end(row);
}


SparseMatrix::SparseMatrix()
{
  reinit(0, 0);
}


SparseMatrix::SparseMatrix(const size_t n_rows, const size_t n_cols)
{
  reinit(n_rows, n_cols);
}


void
SparseMatrix::reinit(const size_t n_rows, const size_t n_cols)
{
  has_entries = false;
  rows = n_rows;
  cols = n_cols;

  colnums.clear();
  colnums.resize(n_rows);

  values.clear();
  values.resize(n_rows);
}


SparseMatrix&
SparseMatrix::operator=(const double value)
{
  for (auto& row : values)
    for (auto& el : row)
      el = value;
  return *this;
}


void
SparseMatrix::copy_from(const Matrix& matrix)
{
  rows = matrix.n_rows();
  cols = matrix.n_cols();

  colnums.resize(n_rows());
  values.resize(n_rows());

  for (size_t row = 0; row < rows; ++row)
    for (size_t col = 0; col < cols; ++col)
      if (matrix(row, col) != 0.0)
        set(row, col, matrix(row, col));
  has_entries = true;
}


void
SparseMatrix::copy_from(const SparseMatrix& matrix)
{
  rows = matrix.n_rows();
  cols = matrix.n_cols();

  colnums = matrix.colnums;
  values = matrix.values;
  has_entries = matrix.has_entries;
}


size_t
SparseMatrix::n_rows() const
{
  return rows;
}


size_t
SparseMatrix::n_cols() const
{
  return cols;
}


size_t
SparseMatrix::n_nonzero_entries() const
{
  return std::accumulate(colnums.begin(), colnums.end(), 0,
                         [](size_t v, const std::vector<size_t>& row)
                         { return v + row.size(); });
}


unsigned int
SparseMatrix::row_length(const size_t row) const
{
  assert(row < rows);
  assert(colnums[row].size() == values[row].size());
  return colnums[row].size();
}


size_t
SparseMatrix::column(const size_t row, const unsigned int index) const
{
  assert(row < rows);
  assert(index < row_length(row));
  return colnums[row][index];
}


unsigned int
SparseMatrix::index(const size_t row, const size_t column) const
{
  assert(row < rows);
  assert(column < cols);

  // Binary search for the
  auto it = std::lower_bound(colnums[row].begin(),
                             colnums[row].end(), column);
  if (it != colnums[row].end() && *it == column)
    return it - colnums[row].begin();
  else
    return -1;
}


bool
SparseMatrix::empty() const
{
  return (rows == 0 && cols == 0 &&
          colnums.empty() && values.empty());
}


bool
SparseMatrix::exists(const size_t row, const size_t column) const
{
  return index(row, column) != -1;
}


bool
SparseMatrix::operator==(const SparseMatrix& other) const
{
  return (has_entries == other.has_entries &&
          rows == other.rows && cols == other.cols &&
          colnums == other.colnums && values == other.values);
}


bool
SparseMatrix::operator!=(const SparseMatrix& other) const
{ return !(*this == other); }


SparseMatrix::Iterator
SparseMatrix::begin()
{
  return (rows > 0)? begin(0) : end();
}


SparseMatrix::Iterator
SparseMatrix::end()
{
  return {this};
}


SparseMatrix::Iterator
SparseMatrix::begin(const size_t row)
{
  assert(row < rows);

  // If empty, return end
  if (!has_entries)
    return end();

  // Find the next non-empty row
  size_t r = row;
  while (r != row && colnums[r].size() == 0)
    ++r;

  // Return the appropriate row, or end
  if (r == rows) return {this};
  else return {this, r, 0};
}


SparseMatrix::Iterator
SparseMatrix::end(const size_t row)
{
  assert(row < rows);
  return (row + 1 == rows)? end() : begin(row + 1);
}


SparseMatrix::ConstIterator
SparseMatrix::begin() const
{
  return (rows > 0)? begin(0) : end();
}


SparseMatrix::ConstIterator
SparseMatrix::end() const
{
  return {this};
}


SparseMatrix::ConstIterator
SparseMatrix::begin(const size_t row) const
{
  assert(row < rows);

  // If empty, return end
  if (!has_entries)
    return end();

  // Find the next non-empty row
  size_t r = row;
  while (r != row && colnums[r].size() == 0)
    ++r;

  // Return the appropriate row, or end
  if (r == rows) return {this};
  else return {this, r, 0};
}


SparseMatrix::ConstIterator
SparseMatrix::end(const size_t row) const
{
  assert(row < rows);
  return (row + 1 == rows)? end() : begin(row + 1);
}


SparseMatrix::RowIterator
SparseMatrix::row_iterator(const size_t row)
{
  assert(row < rows);
  return {this, row};
}


SparseMatrix::ConstRowIterator
SparseMatrix::row_iterator(const size_t row) const
{
  assert(row < rows);
  return {this, row};
}


double&
SparseMatrix::operator()(const size_t i, const size_t j)
{
  const unsigned int idx = index(i, j);
  assert(idx != -1);
  return values[i][idx];
}


const double&
SparseMatrix::operator()(const size_t i, const size_t j) const
{
  const unsigned int idx = index(i, j);
  assert(idx != -1);
  return values[i][idx];
}


double
SparseMatrix::el(const size_t i, const size_t j)
{
  const unsigned int idx = index(i, j);
  return (idx != -1) ? values[i][idx] : 0.0;
}


double&
SparseMatrix::diag(const size_t i)
{
  const unsigned int idx = index(i, i);
  assert(idx != -1);
  return values[i][idx];
}


const double&
SparseMatrix::diag(const size_t i) const
{
  const unsigned int idx = index(i, i);
  assert(idx != -1);
  return values[i][idx];
}


double
SparseMatrix::diag_el(const size_t i) const
{
  const unsigned int idx = index(i, i);
  return (idx != -1) ? values[i][idx] : 0.0;
}


void
SparseMatrix::clear()
{
  has_entries = false;
  rows = cols = 0;
  colnums.clear();
  values.clear();
}


void
SparseMatrix::set(const size_t row, const size_t column, const double value)
{
  assert(row < rows);
  assert(column < cols);
  if (value == 0.0)
    return;

  has_entries = true;

  // If the row is empty or the column number is larger than all
  // current entries on the row, add to the back of the row.
  if (colnums[row].size() == 0 || colnums[row].back() < column)
  {
    colnums[row].push_back(column);
    values[row].push_back(value);
    return;
  }

  // Find the index to insert which maintains column ordering
  auto it = std::lower_bound(colnums[row].begin(), colnums[row].end(), column);
  size_t idx = it - colnums[row].begin();

  // Override the column if it already exists
  if (*it == column)
  {
    values[row][idx] = value;
    return;
  }

  // If it does not, insert it
  colnums[row].insert(it, column);
  values[row].insert(values[row].begin() + idx, value);
}


void
SparseMatrix::add(const size_t row, const size_t column, const double value)
{
  assert(row < rows);
  assert(column < cols);
  if (value == 0)
    return;

  has_entries = true;

  const unsigned int idx = index(row, column);
  if (idx == -1) set(row, column, value);
  else values[row][idx] += value;
}


void
SparseMatrix::swap_row(const size_t i, const size_t k)
{
  assert(i < rows);
  assert(k < rows);
  colnums[i].swap(colnums[k]);
  values[i].swap(values[k]);
}


void
SparseMatrix::swap(SparseMatrix& other)
{
  std::swap(has_entries, other.has_entries);
  std::swap(rows, other.rows);
  std::swap(cols, other.cols);

  colnums.swap(other.colnums);
  values.swap(other.values);
}


SparseMatrix&
SparseMatrix::scale(const double factor)
{
  for (auto& row_vals : values)
    for (auto& el : row_vals)
      el *= factor;
  return *this;
}


SparseMatrix&
SparseMatrix::operator-()
{
  return scale(-1.0);
}


SparseMatrix
SparseMatrix::operator-() const
{
  return SparseMatrix(*this).scale(-1.0);
}


SparseMatrix&
SparseMatrix::operator*=(const double factor)
{
  return scale(factor);
}


SparseMatrix
SparseMatrix::operator*(const double factor) const
{
  return SparseMatrix(*this).scale(factor);
}


SparseMatrix&
SparseMatrix::operator/=(const double factor)
{
  assert(factor != 0.0);
  return scale(1.0/factor);
}


SparseMatrix
SparseMatrix::operator/(const double factor) const
{
  assert(factor != 0.0);
  return SparseMatrix(*this).scale(1.0/factor);
}


SparseMatrix&
SparseMatrix::sadd(const double a, const double b, const SparseMatrix& B)
{
  assert(colnums == B.colnums);
  for (size_t row = 0; row < rows; ++row)
  {
    double* a_ij = &values[row][0];
    const double* b_ij = &B.values[row][0];
    const double* const eor = a_ij + row_length(row);

    for (; a_ij != eor; ++a_ij)
      *a_ij = a * *a_ij + b * *b_ij;
  }
  return *this;
}


SparseMatrix&
SparseMatrix::sadd(const double a,
                   const SparseMatrix& B)
{
  return sadd(a, 1.0, B);
}



SparseMatrix&
SparseMatrix::add(const double b, const SparseMatrix& B)
{
  return sadd(1.0, b, B);
}


SparseMatrix&
SparseMatrix::operator+=(const SparseMatrix& B)
{
  return add(1.0, B);
}


SparseMatrix
SparseMatrix::operator+(const SparseMatrix& B) const
{
  return SparseMatrix(*this).add(1.0, B);
}


SparseMatrix&
SparseMatrix::operator-=(const SparseMatrix& B)
{
  return add(-1.0, B);
}


SparseMatrix
SparseMatrix::operator-(const SparseMatrix& B) const
{
  return SparseMatrix(*this).add(-1.0, B);
}


void
SparseMatrix::vmult(const Vector& x, Vector& y, const bool adding) const
{
  assert(x.size() == cols);
  assert(y.size() == rows);

  double* dst_ptr = &y[0];

  for (size_t row = 0; row < rows; ++row)
  {
    const size_t* col_ptr = &colnums[row][0];
    const double* a_ij = &values[row][0];
    const double* const eor = a_ij + row_length(row);

    double val = adding? *dst_ptr : 0.0;
    while (a_ij != eor)
      val += *a_ij++*x[*col_ptr++];
    *dst_ptr++ = val;
  }
}


void
SparseMatrix::vmult_add(const Vector& x, Vector& y) const
{
  vmult(x, y, true);
}


Vector
SparseMatrix::vmult(const Vector& x) const
{
  Vector y(rows);
  vmult(x, y);
  return y;
}


Vector
SparseMatrix::operator*(const Vector& x) const
{
  return vmult(x);
}


void
SparseMatrix::Tvmult(const Vector& x, Vector& y, const bool adding) const
{
  assert(x.size() == cols);
  assert(y.size() == rows);

  double* dst_ptr = &y[0];
  const double* x_ptr = &x[0];

  for (size_t row = 0; row < rows; ++row, ++x_ptr)
  {
    const size_t* col_ptr = &colnums[row][0];
    const double* a_ij = &values[row][0];
    const double* const eor = a_ij + row_length(row);

    while (a_ij != eor)
      dst_ptr[*col_ptr++] += *a_ij++**x_ptr;
  }
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
{
  vmult(x, y, true);
}


std::string
SparseMatrix::str(const bool formatted,
                  const bool scientific,
                  const unsigned int precision,
                  const unsigned int width) const
{
  assert(!empty());
  std::stringstream ss;

  unsigned int w = width;
  if (scientific)
    ss.setf(std::ios::scientific, std::ios::floatfield);
  else
    ss.setf(std::ios::fixed, std::ios::floatfield);
  if (!width)
    w = scientific? precision + 9 : precision + 4;

  if (formatted)
  {
    for (size_t row = 0; row < rows; ++row)
    {
      for (size_t col = 0; col < cols; ++col)
      {
        const unsigned int idx = index(row, col);
        if (idx != -1)
          ss << std::setw(w) << values[row][idx] << " ";
        else
          ss << std::setw(w) << std::to_string(0.0) << " ";
      }
      ss << std::endl;
    }
  }
  else
  {
    ss << "(row, column):\tvalue" << std::endl
       << "-----------------------" << std::endl;
    for (size_t row = 0; row < rows; ++row)
      for (unsigned int index = 0; index < row_length(row); ++index)
        ss << "(" << row << ", " << colnums[row][index] << ")\t"
           << values[row][index] << std::endl;
    ss << std::endl;
  }
  return ss.str();
}



void
SparseMatrix::print(std::ostream& os,
                    const bool formatted,
                    const bool scientific,
                    const unsigned int precision,
                    const unsigned int width) const
{
  os << str(formatted, scientific, precision, width);
}


void
SparseMatrix::print_row(const size_t row,
                        std::ostream& os,
                        const bool scientific,
                        const unsigned int precision) const
{
  assert(!empty());

  // setup output stream format
  std::ios::fmtflags old_flags = os.flags();
  unsigned int old_precision = os.precision(precision);
  if (scientific)
    os.setf(std::ios::scientific, std::ios::floatfield);
  else
    os.setf(std::ios::fixed, std::ios::floatfield);

  os << "(row,col):\tvalue" << std::endl
     << "----------------" << std::endl;
  for (const auto& el : row_iterator(row))
    os << "(" << el.row
       << ", " << el.column
       << ")\t" << el.value << std::endl;
  os << std::endl;

  // reset output stream format
  os.precision(old_precision);
  os.flags(old_flags);
}



std::ostream&
Math::operator<<(std::ostream& os, const SparseMatrix& A)
{
  return os << A.str();
}

#include "sparse_matrix.h"
#include "matrix.h"
#include "vector.h"
#include "macros.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>


using namespace Math;
using namespace SparseMatrixIterators;


//######################################################################


Iterator::Accessor::
Accessor(SparseMatrix* matrix,
         const size_t row,
         const size_t index) :
  matrix(matrix),
  current_row(row),
  current_index(index)
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
}


Iterator::Accessor::
Accessor(SparseMatrix* matrix) :
  matrix(matrix),
  current_row(SparseMatrix::invalid_entry),
  current_index(SparseMatrix::invalid_entry)
{}


size_t
Iterator::Accessor::row() const
{
  assert(current_row < matrix->rows);
  return current_row;
}


size_t
Iterator::Accessor::column() const
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
  return matrix->colnums[current_row][current_index];
}


size_t
Iterator::Accessor::index() const
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
  return current_index;
}


double&
Iterator::Accessor::value() const
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
  return matrix->vals[current_row][current_index];
}


bool
Iterator::Accessor::operator==(const Accessor& other) const
{
  return (matrix == other.matrix &&
          current_row == other.current_row &&
          current_index == other.current_index);
}


bool
Iterator::Accessor::operator!=(const Accessor& other) const
{ return !(*this == other); }


bool
Iterator::Accessor::operator<(const Accessor& other) const
{
  assert(matrix == other.matrix);

  // If at end on matrix, nothing is less
  if (current_row == SparseMatrix::invalid_entry)
    return false;
  assert(current_row < matrix->rows);

  // If other is at end, everything is less
  if (other.current_row == SparseMatrix::invalid_entry)
    return true;
  assert(other.current_row < matrix->rows);

  // If here, do the real comparison
  return (current_row < other.current_row ||
          (current_row == other.current_row &&
           current_index < other.current_index));
}


void
Iterator::Accessor::advance()
{
  // Increment the position on the current row
  ++current_index;

  // Handle the end of a row
  if (current_index == matrix->row_length(current_row))
  {
    /* If the next row is not the last row, move the start of it,
     * otherwise, construct the invalid accessor for the end iterator. */
    if (current_row + 1 < matrix->rows)
      *this = *matrix->begin(current_row + 1);
    else
      *this = *matrix->end();
  }
}


Iterator::Iterator(SparseMatrix* matrix,
                   const size_t row,
                   const size_t index) :
  accessor(matrix, row, index)
{}


Iterator::Iterator(SparseMatrix* matrix) :
  accessor(matrix)
{}


Iterator&
Iterator::operator++()
{
  accessor.advance();
  return *this;
}


Iterator
Iterator::operator++(int)
{
  auto it = *this;
  accessor.advance();
  return it;
}


const Iterator::Accessor&
Iterator::operator*() const
{ return accessor; }


const Iterator::Accessor*
Iterator::operator->() const
{ return &accessor; }


bool
Iterator::operator==(const Iterator& other) const
{ return accessor == other.accessor; }


bool
Iterator::operator!=(const Iterator& other) const
{ return accessor != other.accessor; }


bool
Iterator::operator<(const Iterator& other) const
{ return accessor < other.accessor; }


//######################################################################


ConstIterator::Accessor::
Accessor(const SparseMatrix* matrix,
         const size_t row,
         const size_t index) :
  matrix(matrix),
  current_row(row),
  current_index(index)
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
}


ConstIterator::Accessor::
Accessor(const SparseMatrix* matrix) :
  matrix(matrix),
  current_row(SparseMatrix::invalid_entry),
  current_index(SparseMatrix::invalid_entry)
{}


size_t
ConstIterator::Accessor::row() const
{
  assert(current_row < matrix->rows);
  return current_row;
}


size_t
ConstIterator::Accessor::column() const
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
  return matrix->colnums[current_row][current_index];
}


size_t
ConstIterator::Accessor::index() const
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
  return current_index;
}


const double&
ConstIterator::Accessor::value() const
{
  assert(current_row < matrix->rows);
  assert(current_index < matrix->row_length(current_row));
  return matrix->vals[current_row][current_index];
}


bool
ConstIterator::Accessor::
operator==(const Accessor& other) const
{
  return (matrix == other.matrix &&
          current_row == other.current_row &&
          current_index == other.current_index);
}


bool
ConstIterator::Accessor::
operator!=(const Accessor& other) const
{ return !(*this == other); }


bool
ConstIterator::Accessor::
operator<(const Accessor& other) const
{
  assert(matrix == other.matrix);

  // If at end on matrix, nothing is less
  if (current_row == SparseMatrix::invalid_entry)
    return false;
  assert(current_row < matrix->rows);

  // If other is at end, everything is less
  if (other.current_row == SparseMatrix::invalid_entry)
    return true;
  assert(other.current_row < matrix->rows);

  // If here, do the real comparison
  return (current_row < other.current_row ||
          (current_row == other.current_row &&
           current_index < other.current_index));
}


void
ConstIterator::Accessor::advance()
{
  // Increment the position on the current row
  ++current_index;

  // Handle the end of a row
  if (current_index == matrix->row_length(current_row))
  {
    /* If the next row is not the last row, move the start of it,
     * otherwise, construct the invalid accessor for the end iterator. */
    if (current_row + 1 < matrix->rows)
      *this = *matrix->begin(current_row + 1);
    else
      *this = *matrix->end();
  }
}


ConstIterator::
ConstIterator(const SparseMatrix* matrix,
              const size_t row,
              const size_t index) :
  accessor(matrix, row, index)
{}


ConstIterator::
ConstIterator(const SparseMatrix* matrix) :
  accessor(matrix)
{}


ConstIterator&
ConstIterator::operator++()
{
  accessor.advance();
  return *this;
}


ConstIterator
ConstIterator::operator++(int)
{
  auto it = *this;
  accessor.advance();
  return it;
}


const ConstIterator::Accessor&
ConstIterator::operator*() const
{ return accessor; }


const ConstIterator::Accessor*
ConstIterator::operator->() const
{ return &accessor; }


bool
ConstIterator::operator==(const ConstIterator& other) const
{ return accessor == other.accessor; }


bool
ConstIterator::operator!=(const ConstIterator& other) const
{ return accessor != other.accessor; }


bool
ConstIterator::operator<(const ConstIterator& other) const
{ return accessor < other.accessor; }


//######################################################################


RowIterator::RowIterator(SparseMatrix* matrix,
                         const size_t row) :
  matrix(matrix), row(row)
{ assert(row < matrix->rows); }


Iterator
RowIterator::begin()
{ return matrix->begin(row); }


Iterator
RowIterator::end()
{ return matrix->end(row); }


//######################################################################


ConstRowIterator::
ConstRowIterator(const SparseMatrix* matrix,
                 const size_t row) :
  matrix(matrix), row(row)
{ assert(row < matrix->rows); }


ConstIterator
ConstRowIterator::begin() const
{ return matrix->begin(row); }


ConstIterator
ConstRowIterator::end() const
{ return matrix->end(row); }


//######################################################################


SparseMatrix::SparseMatrix() :
  has_entries(false),
  rows(0),
  cols(0),
  colnums(0),
  vals(0)
{}


SparseMatrix::SparseMatrix(const size_t n) :
  SparseMatrix()
{ reinit(n, n); }


SparseMatrix::SparseMatrix(const size_t n_rows,
                           const size_t n_cols) :
  SparseMatrix()
{ reinit(n_rows, n_cols); }


void
SparseMatrix::reinit(const size_t n_rows,
                     const size_t n_cols)
{
  has_entries = false;
  rows = n_rows;
  cols = n_cols;

  colnums.clear();
  colnums.resize(n_rows);

  vals.clear();
  vals.resize(n_rows);
}


SparseMatrix&
SparseMatrix::operator=(const value_type value)
{
  for (auto& row : vals)
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
  vals.resize(n_rows());

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
  vals = matrix.vals;
  has_entries = matrix.has_entries;
}


size_t
SparseMatrix::n_rows() const
{ return rows; }


size_t
SparseMatrix::n_cols() const
{ return cols; }


size_t
SparseMatrix::nnz() const
{
  return std::accumulate(colnums.begin(), colnums.end(), 0,
                         [](size_t v, const std::vector<size_t>& row)
                         { return v + row.size(); });
}


size_t
SparseMatrix::row_length(const size_t row) const
{
  assert(row < rows);
  assert(colnums[row].size() == vals[row].size());
  return colnums[row].size();
}


size_t
SparseMatrix::column(const size_t row, const size_t index) const
{
  assert(row < rows);
  assert(index < row_length(row));
  return colnums[row][index];
}


size_t
SparseMatrix::index(const size_t row, const size_t col) const
{
  assert(row < rows);
  assert(col < cols);

  // Binary search for the
  auto it = std::lower_bound(colnums[row].begin(),
                             colnums[row].end(), col);
  if (it != colnums[row].end() && *it == col)
    return it - colnums[row].begin();
  else
    return invalid_entry;
}


bool
SparseMatrix::empty() const
{
  return (rows == 0 && cols == 0 &&
          colnums.empty() && vals.empty());
}


bool
SparseMatrix::exists(const size_t row, const size_t col) const
{ return (index(row, col) != invalid_entry)? true : false; }


bool
SparseMatrix::operator==(const SparseMatrix& other) const
{
  return (has_entries == other.has_entries &&
          rows == other.rows && cols == other.cols &&
          colnums == other.colnums && vals == other.vals);
}


bool
SparseMatrix::operator!=(const SparseMatrix& other) const
{ return !(*this == other); }


Iterator
SparseMatrix::begin()
{ return (rows > 0)? begin(0) : end(); }


Iterator
SparseMatrix::end()
{ return {this}; }


Iterator
SparseMatrix::begin(const size_t row)
{
  assert(row < rows);
  if (!has_entries)
    return end();

  size_t r = row;
  while (r != row && colnums[r].size() == 0)
    ++r;

  if (r == rows) return {this};
  else return {this, r, 0};
}


Iterator
SparseMatrix::end(const size_t row)
{
  assert(row < rows);
  return (row + 1 == rows)? end() : begin(row + 1);
}


ConstIterator
SparseMatrix::begin() const
{ return (rows > 0)? begin(0) : end(); }


ConstIterator
SparseMatrix::end() const
{ return {this}; }


ConstIterator
SparseMatrix::begin(const size_t row) const
{
  assert(row < rows);
  if (!has_entries)
    return end();

  size_t r = row;
  while (r != row && colnums[r].size() == 0)
    ++r;

  if (r == rows) return {this};
  else return {this, r, 0};
}


ConstIterator
SparseMatrix::end(const size_t row) const
{
  assert(row < rows);
  return (row + 1 == rows)? end() : begin(row + 1);
}


RowIterator
SparseMatrix::row_iterator(const size_t row)
{
  assert(row < rows);
  return {this, row};
}


ConstRowIterator
SparseMatrix::row_iterator(const size_t row) const
{
  assert(row < rows);
  return {this, row};
}


double&
SparseMatrix::operator()(const size_t i, const size_t j)
{
  const size_t idx = index(i, j);
  assert(idx != invalid_entry);
  return vals[i][idx];
}


const double&
SparseMatrix::operator()(const size_t i, const size_t j) const
{
  const size_t idx = index(i, j);
  assert(idx != invalid_entry);
  return vals[i][idx];
}


double
SparseMatrix::el(const size_t i, const size_t j)
{
  const size_t idx = index(i, j);
  return (idx != invalid_entry)? vals[i][idx] : 0.0;
}


double&
SparseMatrix::diag(const size_t i)
{
  const size_t idx = index(i, i);
  assert(idx != invalid_entry);
  return vals[i][idx];
}


const double&
SparseMatrix::diag(const size_t i) const
{
  const size_t idx = index(i, i);
  assert(idx != invalid_entry);
  return vals[i][idx];
}


double
SparseMatrix::diag_el(const size_t i) const
{
  const size_t idx = index(i, i);
  return (idx != invalid_entry)? vals[i][idx] : 0.0;
}


void
SparseMatrix::clear()
{
  has_entries = false;
  rows = cols = 0;
  colnums.clear();
  vals.clear();
}


void
SparseMatrix::set(const size_t row,
                  const size_t col,
                  const value_type value)
{
  assert(row < rows);
  assert(col < cols);
  if (value == 0.0)
    return;

  has_entries = true;

  /* If the row is empty or the column number is larger than all
   * current entries on the row, add to the back of the row. */
  if (colnums[row].size() == 0 || colnums[row].back() < col)
  {
    colnums[row].push_back(col);
    vals[row].push_back(value);
    return;
  }

  // Find the index to insert which maintains column ordering
  auto it = std::lower_bound(colnums[row].begin(), colnums[row].end(), col);
  size_t idx = it - colnums[row].begin();

  // Override the column if it already exists
  if (*it == col)
  {
    vals[row][idx] = value;
    return;
  }

  // If it does not, insert it
  colnums[row].insert(it, col);
  vals[row].insert(vals[row].begin() + idx, value);
}


void
SparseMatrix::add(const size_t row,
                  const size_t col,
                  const value_type value)
{
  assert(row < rows);
  assert(col < cols);
  if (value == 0)
    return;

  has_entries = true;

  const size_t idx = index(row, col);
  if (idx == invalid_entry) set(row, col, value);
  else vals[row][idx] += value;
}


void
SparseMatrix::swap_row(const size_t i, const size_t k)
{
  assert(i < rows);
  assert(k < rows);
  colnums[i].swap(colnums[k]);
  vals[i].swap(vals[k]);
}


void
SparseMatrix::swap(SparseMatrix& other)
{
  std::swap(has_entries, other.has_entries);
  std::swap(rows, other.rows);
  std::swap(cols, other.cols);

  colnums.swap(other.colnums);
  vals.swap(other.vals);
}


SparseMatrix&
SparseMatrix::scale(const value_type factor)
{
  for (auto& row_vals : vals)
    for (auto& el : row_vals)
      el *= factor;
  return *this;
}


SparseMatrix&
SparseMatrix::add(const SparseMatrix& B,
                  const value_type factor)
{
  assert(colnums == B.colnums);
  for (size_t row = 0; row < rows; ++row)
  {
    value_type* a_ij = &vals[row][0];
    const value_type* b_ij = &B.vals[row][0];
    const value_type* const eor = a_ij + row_length(row);

    while (a_ij != eor)
      *a_ij++ += factor**b_ij++;
  }
  return *this;
}


SparseMatrix&
SparseMatrix::sadd(const value_type a,
                   const SparseMatrix& B)
{
  assert(colnums == B.colnums);
  for (size_t row = 0; row < rows; ++row)
  {
    value_type* a_ij = &vals[row][0];
    const value_type* b_ij = &B.vals[row][0];
    const value_type* const eor = a_ij + row_length(row);

    for (; a_ij != eor; ++a_ij)
      *a_ij = a**a_ij + *b_ij;
  }
  return *this;
}


SparseMatrix&
SparseMatrix::sadd(const value_type a,
                   const value_type b,
                   const SparseMatrix& B)
{
  assert(colnums == B.colnums);
  for (size_t row = 0; row < rows; ++row)
  {
    value_type* a_ij = &vals[row][0];
    const value_type* b_ij = &B.vals[row][0];
    const value_type* const eor = a_ij + row_length(row);

    for (; a_ij != eor; ++a_ij)
      *a_ij = a**a_ij + b**b_ij;
  }
  return *this;
}


void
SparseMatrix::vmult(const Vector& x,
                    Vector& y,
                    const bool adding) const
{
  assert(x.size() == cols);
  assert(y.size() == rows);

  value_type* dst_ptr = &y[0];

  for (size_t row = 0; row < rows; ++row)
  {
    const size_t* col_ptr = &colnums[row][0];
    const value_type* a_ij = &vals[row][0];
    const value_type* const eor = a_ij + row_length(row);

    value_type val = adding? *dst_ptr : 0.0;
    while (a_ij != eor)
      val += *a_ij++*x[*col_ptr++];
    *dst_ptr++ = val;
  }
}


void
SparseMatrix::Tvmult(const Vector& x,
                     Vector& y,
                     const bool adding) const
{
  assert(x.size() == cols);
  assert(y.size() == rows);

  value_type* dst_ptr = &y[0];
  const value_type* x_ptr = &x[0];

  for (size_t row = 0; row < rows; ++row, ++x_ptr)
  {
    const size_t* col_ptr = &colnums[row][0];
    const value_type* a_ij = &vals[row][0];
    const value_type* const eor = a_ij + row_length(row);

    while (a_ij != eor)
      dst_ptr[*col_ptr++] += *a_ij++**x_ptr;
  }
}


Vector
SparseMatrix::vmult(const Vector& x) const
{
  Vector y(rows);
  vmult(x, y);
  return y;
}


Vector
SparseMatrix::Tvmult(const Vector& x) const
{
  Vector y(cols);
  Tvmult(x, y);
  return y;
}


void
SparseMatrix::vmult_add(const Vector& x, Vector& y) const
{ vmult(x, y, true); }


void
SparseMatrix::Tvmult_add(const Vector& x, Vector& y) const
{ vmult(x, y, true); }


SparseMatrix&
SparseMatrix::operator-()
{ return scale(-1.0); }


SparseMatrix&
SparseMatrix::operator*=(const value_type factor)
{ return scale(factor); }


SparseMatrix&
SparseMatrix::operator/=(const value_type factor)
{
  assert(factor != 0.0);
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


void
SparseMatrix::print(std::ostream& os,
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

  os << "(row,col): value" << std::endl
     << "----------------" << std::endl;
  for (size_t row = 0; row < rows; ++row)
    for (size_t index = 0; index < row_length(row); ++index)
      os << "(" << row << "," << colnums[row][index] << ") "
         << vals[row][index] << std::endl;
  os << std::endl;

  // reset output stream format
  os.precision(old_precision);
  os.flags(old_flags);
}


void
SparseMatrix::print_formatted(std::ostream& os,
                              const bool scientific,
                              const unsigned int precision,
                              const unsigned int width) const
{
  assert(!empty());

  // setup output stream format
  unsigned int w = width;
  std::ios::fmtflags old_flags = os.flags();
  unsigned int old_precision = os.precision(precision);
  if (scientific)
    os.setf(std::ios::scientific, std::ios::floatfield);
  else
    os.setf(std::ios::fixed, std::ios::floatfield);
  if (!width)
    w = scientific? precision + 9 : precision + 4;

  for (size_t row = 0; row < rows; ++row)
  {
    for (size_t col = 0; col < cols; ++col)
    {
      const size_t idx = index(row, col);
      if (idx != invalid_entry)
        os << std::setw(w) << vals[row][idx] << " ";
      else
        os << std::setw(w) << std::to_string(0.0) << " ";
    }
    os << std::endl;
  }

  // reset output stream format
  os.precision(old_precision);
  os.flags(old_flags);
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


std::ostream&
Math::operator<<(std::ostream& os, const SparseMatrix& A)
{ return os << A.str(); }

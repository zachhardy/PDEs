#include "sparse_matrix.h"

using namespace pdes::Math;


//################################################## Entry

SparseMatrix::entry::
entry(const size_t& i, const size_t& j, double& val) :
    row(i), column(j), value(val)
{}


std::string
SparseMatrix::entry::
str() const
{
  std::stringstream ss;
  ss << "A(" << row << ", " << column << ") = " << value;
  return ss.str();
}

//################################################## Const Entry

SparseMatrix::const_entry::
const_entry(const size_t& i, const size_t& j, const double& val) :
    row(i), column(j), value(val)
{}


std::string
SparseMatrix::const_entry::
str() const
{
  std::stringstream ss;
  ss << "A(" << row << ", " << column << ") = " << value;
  return ss.str();
}


//################################################## Iterator

SparseMatrix::iterator::
iterator(SparseMatrix* sparse_matrix, const size_t row) :
    sparse_matrix_ptr(sparse_matrix),
    current_row(row),
    col_ptr(sparse_matrix_ptr->colnums[current_row].begin()),
    val_ptr(sparse_matrix_ptr->coeffs[current_row].begin())
{}


SparseMatrix::iterator::
iterator(SparseMatrix* sparse_matrix) :
    sparse_matrix_ptr(sparse_matrix),
    current_row(-1), col_ptr(), val_ptr()
{}


void
SparseMatrix::iterator::
advance()
{
  // Increment along the current row
  ++col_ptr; ++val_ptr;

  // If at the end of a row, handle it
  if (col_ptr == sparse_matrix_ptr->colnums[current_row].end())
  {
    // Increment the row to the next non-empty row.
    ++current_row;
    while (current_row < sparse_matrix_ptr->rows &&
           sparse_matrix_ptr->colnums.empty())
      ++current_row;

    /* Set the pointers to the next row, for valid rows, or to the invalid
     * iterator, if at the end of the matrix. */
    if (current_row < sparse_matrix_ptr->rows)
      *this = sparse_matrix_ptr->begin(current_row);
    else
      *this = sparse_matrix_ptr->end();
  }
}


SparseMatrix::iterator&
SparseMatrix::iterator::
operator++()
{ advance(); return *this; }


SparseMatrix::iterator
SparseMatrix::iterator::
operator++(int)
{ auto it = *this; advance(); return it; }


SparseMatrix::entry
SparseMatrix::iterator::
operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


bool
SparseMatrix::iterator::
operator==(const iterator& other) const
{
  return (sparse_matrix_ptr == other.sparse_matrix_ptr &&
          current_row == other.current_row &&
          col_ptr == other.col_ptr &&
          val_ptr == other.val_ptr);
}


bool
SparseMatrix::iterator::
operator!=(const iterator& other) const
{ return !(*this == other); }


//################################################## Constant Iterator


SparseMatrix::const_iterator::
const_iterator(const SparseMatrix* sparse_matrix, const size_t row) :
    sparse_matrix_ptr(sparse_matrix),
    current_row(row),
    col_ptr(sparse_matrix_ptr->colnums[current_row].begin()),
    val_ptr(sparse_matrix_ptr->coeffs[current_row].begin())
{}


SparseMatrix::const_iterator::
const_iterator(const SparseMatrix* sparse_matrix) :
    sparse_matrix_ptr(sparse_matrix),
    current_row(-1), col_ptr(), val_ptr()
{}


SparseMatrix::row_accessor::iterator::
iterator(SparseMatrix::row_accessor& accessor, const size_t elem) :
    current_row(accessor.row_num),
    col_ptr(accessor.colnums.begin() + elem),
    val_ptr(accessor.coeffs.begin() + elem)
{}


SparseMatrix::const_row_accessor::const_iterator::
const_iterator(const const_row_accessor& accessor, const size_t elem) :
    current_row(accessor.row_num),
    col_ptr(accessor.colnums.begin() + elem),
    val_ptr(accessor.coeffs.begin() + elem)
{}


void
SparseMatrix::const_iterator::
advance()
{
  // Increment along the current row
  ++col_ptr; ++val_ptr;

  // If at the end of a row, handle it
  if (col_ptr == sparse_matrix_ptr->colnums[current_row].end())
  {
    // Increment the row to the next non-empty row.
    ++current_row;
    while (current_row < sparse_matrix_ptr->rows &&
           sparse_matrix_ptr->colnums.empty())
      ++current_row;

    /* Set the pointers to the next row, for valid rows, or to the invalid
     * iterator, if at the end of the matrix. */
    if (current_row < sparse_matrix_ptr->rows)
      *this = sparse_matrix_ptr->begin(current_row);
    else
      *this = sparse_matrix_ptr->end();
  }
}


SparseMatrix::const_iterator&
SparseMatrix::const_iterator::
operator++()
{ advance(); return *this; }


SparseMatrix::const_iterator
SparseMatrix::const_iterator::
operator++(int)
{ auto it = *this; advance(); *this; return it; }


SparseMatrix::const_entry
SparseMatrix::const_iterator::
operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


bool
SparseMatrix::const_iterator::
operator==(const const_iterator& other) const
{
  return (sparse_matrix_ptr == other.sparse_matrix_ptr &&
          current_row == other.current_row &&
          col_ptr == other.col_ptr &&
          val_ptr == other.val_ptr);
}


bool
SparseMatrix::const_iterator::
operator!=(const const_iterator& other) const
{ return !(*this == other); }


//################################################## Row Accessor


SparseMatrix::row_accessor::
row_accessor(SparseMatrix& sparse_matrix, const size_t i) :
    row_num(i), colnums(sparse_matrix.colnums[row_num]),
    coeffs(sparse_matrix.coeffs[row_num])
{}


SparseMatrix::row_accessor::iterator&
SparseMatrix::row_accessor::iterator::
operator++()
{ ++col_ptr; ++val_ptr; return *this; }


SparseMatrix::row_accessor::iterator
SparseMatrix::row_accessor::iterator::
operator++(int)
{ auto it = *this; ++(*this); return it; }


SparseMatrix::entry
SparseMatrix::row_accessor::iterator::
operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


bool
SparseMatrix::row_accessor::iterator::
operator==(const iterator& other) const
{
  return (current_row == other.current_row &&
          col_ptr == other.col_ptr &&
          val_ptr == other.val_ptr);
}


bool
SparseMatrix::row_accessor::iterator::
operator!=(const iterator& other) const
{ return !(*this == other); }


SparseMatrix::row_accessor::iterator
SparseMatrix::row_accessor::
begin()
{ return {*this, 0}; }


SparseMatrix::row_accessor::iterator
SparseMatrix::row_accessor::
end()
{ return {*this, colnums.size()}; }


//################################################## Contant Row Accessor


SparseMatrix::const_row_accessor::
const_row_accessor(const SparseMatrix& sparse_matrix, const size_t i) :
    row_num(i), colnums(sparse_matrix.colnums[row_num]),
    coeffs(sparse_matrix.coeffs[row_num])
{}


SparseMatrix::const_row_accessor::const_iterator&
SparseMatrix::const_row_accessor::const_iterator::
operator++()
{ ++col_ptr; ++val_ptr; return *this; }


SparseMatrix::const_row_accessor::const_iterator
SparseMatrix::const_row_accessor::const_iterator::
operator++(int)
{ auto it = *this; ++(*this); return it; }


SparseMatrix::const_entry
SparseMatrix::const_row_accessor::const_iterator::
operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


bool
SparseMatrix::const_row_accessor::const_iterator::
operator==(const const_iterator& other) const
{
  return (current_row == other.current_row &&
          col_ptr == other.col_ptr &&
          val_ptr == other.val_ptr);
}


bool
SparseMatrix::const_row_accessor::const_iterator::
operator!=(const const_iterator& other) const
{ return !(*this == other); }


SparseMatrix::const_row_accessor::const_iterator
SparseMatrix::const_row_accessor::
begin() const
{ return {*this, 0}; }


SparseMatrix::const_row_accessor::const_iterator
SparseMatrix::const_row_accessor::
end() const
{ return {*this, colnums.size()}; }

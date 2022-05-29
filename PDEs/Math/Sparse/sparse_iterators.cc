#include "sparse_matrix.h"

using namespace pdes::Math;

//################################################## Constructors

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


//################################################## Advance

void
SparseMatrix::iterator::advance()
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


void
SparseMatrix::const_iterator::advance()
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


//################################################## Increment Operators

SparseMatrix::iterator&
SparseMatrix::iterator::operator++()
{ advance(); return *this; }


SparseMatrix::iterator
SparseMatrix::iterator::operator++(int)
{ auto it = *this; advance(); return it; }


SparseMatrix::const_iterator&
SparseMatrix::const_iterator::operator++()
{ advance(); return *this; }


SparseMatrix::const_iterator
SparseMatrix::const_iterator::operator++(int)
{ auto it = *this; advance(); *this; return it; }


SparseMatrix::row_accessor::iterator&
SparseMatrix::row_accessor::iterator::
operator++()
{ ++col_ptr; ++val_ptr; return *this; }


SparseMatrix::row_accessor::iterator
SparseMatrix::row_accessor::iterator::
operator++(int)
{ auto it = *this; ++(*this); return it; }


SparseMatrix::const_row_accessor::const_iterator&
SparseMatrix::const_row_accessor::const_iterator::
operator++()
{ ++col_ptr; ++val_ptr; return *this; }


SparseMatrix::const_row_accessor::const_iterator
SparseMatrix::const_row_accessor::const_iterator::
operator++(int)
{ auto it = *this; ++(*this); return it; }


//################################################## Dereference Operators


SparseMatrix::entry
SparseMatrix::iterator::operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


SparseMatrix::const_entry
SparseMatrix::const_iterator::operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


SparseMatrix::entry
SparseMatrix::row_accessor::iterator::operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


SparseMatrix::const_entry
SparseMatrix::const_row_accessor::const_iterator::operator*()
{ return {current_row, *col_ptr, *val_ptr}; }


//################################################## Comparison Operators


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


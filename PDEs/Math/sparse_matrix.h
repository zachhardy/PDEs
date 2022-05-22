#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "vector.h"
#include "matrix.h"
#include "macros.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cinttypes>


namespace pdes::Math
{


class SparseMatrix
{
public:
  using value_type = double;

private:
  size_t rows;  ///< The number of rows.
  size_t cols;  ///< The number of columns.

  /**
   * The row-wise nonzero column indices.
   */
  std::vector<std::vector<size_t>> colnums;

  /**
   * The row-wise nonzero data entries.
   */
  std::vector<std::vector<value_type>> coeffs;

public:
  /**
   * Default contructor.
   */
  SparseMatrix() :
      rows(0), cols(0), colnums(), coeffs()
  {}

  /**
   * Construct a sparse matrix with \p n_rows and \p n_cols with
   * \p default_row_length nonzero entries per row.
   */
  explicit
  SparseMatrix(const size_t n_rows,
               const size_t n_cols,
               const size_t default_row_length)
      : rows(n_rows), cols(n_cols),
        colnums(n_rows, std::vector<size_t>(default_row_length)),
        coeffs(n_rows, std::vector<value_type>(default_row_length))
  {}

  /**
   * Construct a square sparse matrix with dimension \p n with
   * \p default_row_length nonzero entries per row.
   */
  explicit
  SparseMatrix(const size_t n,
               const size_t default_row_length)
      : SparseMatrix(n, n, default_row_length)
  {}

  /**
   * Construct a sparse matrix from a sparsity pattern. In this context, a
   * sparsity pattern is defined by a list of lists. Each inner list stores
   * the nonzero column indices for a given row.
   */
  SparseMatrix(std::vector<std::vector<size_t>> sparsity_pattern)
      : rows(sparsity_pattern.size()), colnums(sparsity_pattern),
        coeffs(sparsity_pattern.size())
  {
    cols = 0;
    for (size_t i = 0; i < rows; ++i)
    {
      // Sort column indices and allocate data
      std::sort(colnums[i].begin(), colnums[i].end());
      coeffs[i].resize(colnums[i].size());

      // Set number of columns
      for (const auto& col : colnums[i])
        cols = (col > cols) ? col : cols;
    }
  }

  /**
   * Copy construction from a dense matrix.
   */
  SparseMatrix(const Matrix& other)
    : rows(other.n_rows()), cols(other.n_cols()),
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

  /**
   * Assignment with a dense matrix.
   */
  SparseMatrix&
  operator=(const Matrix& other)
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

  /**
   * Assignment with a scalar value. This routine sets all allocated data
   * to the specified value.
   */
  SparseMatrix&
  operator=(const value_type value)
  {
    Assert(!empty(), "Cannot set an empty matrix to a scalar.");
    for (auto elem : *this)
      elem.value = value;
  }

  /**
   * Test the equality of two sparse matrices.
   */
  bool
  operator==(const SparseMatrix& other) const
  {
    return (rows    == other.rows &&
            cols    == other.cols &&
            colnums == other.colnums &&
            coeffs == other.coeffs);
  }

  /**
   * Test the inequality of two sparse matrices.
   * \param other
   * \return
   */
  bool
  operator!=(const SparseMatrix& other) const
  { return !(*this == other); }

  /** \name Characteristices */
  // @{

  /**
   * Return the number of rows.
   */
  size_t
  n_rows() const
  { return rows; }

  /**
   * Return the number of columns.
   */
  size_t
  n_cols() const
  { return cols; }

  /**
   * Return the number of non-zero entries.
   */
  size_t
  nnz() const
  {
    size_t count = 0;
    for (size_t i = 0; i < rows; ++i)
      count += row_length(i);
    return count;
  }

  /**
   * Return the length of row \p i.
   */
  size_t
  row_length(const size_t i) const
  {
    Assert(i < rows, "Out of range error.");
    return colnums[i].size();
  }

  /**
   * Return whether the sparse matrix is empty.
   */
  bool
  empty() const
  {
    return (rows == 0 && cols == 0 &&
            colnums.empty() && coeffs.empty());
  }

  // @}
  /** \name Iterators */
  // @{

  /**
   * A struct defining a mutable entry in the sparse matrix. This acts as a
   * triplet containing a row, column, and value.
   */
  struct entry
  {
    const size_t& row, column;
    value_type& value;

    entry(const size_t& i,
          const size_t& j,
          value_type& val) :
      row(i), column(j), value(val) {}
  };


  /**
   * A custom mutable iterator over the sparse matrix entries. This iterator
   * essentially marches through each entry by storing iterators to the column
   * indices and values. Each increment advances each iterator and when the
   * end of a row is encountered, the pointers are reset to the start of the
   * subsequent row. When the end of the matrix is reached, an invalid iterator
   * with empty column and data iterators is returned.
   */
  class iterator
  {
  private:
    using ConstColumnIterator = std::vector<size_t>::const_iterator;
    using CoeffIterator = std::vector<value_type>::iterator;

  private:
    SparseMatrix*       sparse_matrix_ptr;
    size_t              current_row;
    ConstColumnIterator col_ptr;
    CoeffIterator       coeff_pointer;

    void
    advance()
    {
      // Increment along the current row
      ++col_ptr; ++coeff_pointer;

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

  public:
    iterator(SparseMatrix* sparse_matrix,
             const size_t row) :
        sparse_matrix_ptr(sparse_matrix),
        current_row(row),
        col_ptr(sparse_matrix_ptr->colnums[current_row].begin()),
        coeff_pointer(sparse_matrix_ptr->coeffs[current_row].begin()) {}

    iterator(SparseMatrix* sparse_matrix) :
        sparse_matrix_ptr(sparse_matrix),
        current_row(-1), col_ptr(), coeff_pointer() {}

    iterator&
    operator++()
    { advance(); return *this; }

    iterator
    operator++(int)
    { auto it = *this; advance(); return it; }

    entry
    operator*()
    { return {current_row, *col_ptr, *coeff_pointer}; }

    bool
    operator==(const iterator& other) const
    {
      return (sparse_matrix_ptr == other.sparse_matrix_ptr &&
              current_row == other.current_row &&
              col_ptr == other.col_ptr &&
              coeff_pointer == other.coeff_pointer);
    }

    bool
    operator!=(const iterator& other) const
    { return !(*this == other); }
  };


  /**
   * Return a mutable iterator to the first row of the sparse matrix.
   */
  iterator
  begin()
  { return (!empty()) ? begin(0) : end();  }

  /**
   * Return a mutable iterator to the end of the sparse matrix.
   */
  iterator
  end()
  { return {this}; }

  /**
   * Return a mutable iterator to the start of row \p i.
   */
  iterator
  begin(const size_t i)
  {
    Assert(i < rows, "Out of range error.");
    return {this, i};
  }

  /**
   * Return a mutable iterator to the end of row \p i.
   */
  iterator
  end(const size_t i)
  {
    Assert(i < rows, "Out of range error.");
    if (i + 1 == rows) return end();
    else return begin(i + 1);
  }

  /**
   * A struct defining a constant entry in the sparse matrix.
   * \see SparseMatrix::entry
   */
  struct const_entry
  {
    const size_t& row, column;
    const value_type& value;

    const_entry(const size_t& i,
                const size_t& j,
                const value_type& val) :
        row(i), column(j), value(val) {}
  };


  /**
   * A constant iterator over the elements of the sparse matrix.
   * \see SparseMatrix::iterator
   */
  class const_iterator
  {
  private:
    using ConstColumnIterator = std::vector<size_t>::const_iterator;
    using ConstCoeffIterator = std::vector<value_type>::const_iterator;

  private:
    const SparseMatrix*  sparse_matrix_ptr;
    size_t            current_row;
    ConstColumnIterator  col_ptr;
    ConstCoeffIterator   val_ptr;

    void
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

  public:
    const_iterator(const SparseMatrix* sparse_matrix,
                   const size_t row) :
        sparse_matrix_ptr(sparse_matrix),
        current_row(row),
        col_ptr(sparse_matrix_ptr->colnums[current_row].begin()),
        val_ptr(sparse_matrix_ptr->coeffs[current_row].begin()) {}

    const_iterator(const SparseMatrix* sparse_matrix) :
      sparse_matrix_ptr(sparse_matrix),
      current_row(-1), col_ptr(), val_ptr() {}

    const_iterator&
    operator++()
    { advance(); return *this; }

    const_iterator
    operator++(int)
    { auto it = *this; advance(); return it; }

    const_entry
    operator*()
    { return {current_row, *col_ptr, *val_ptr}; }

    bool
    operator==(const const_iterator& other) const
    {
      return (sparse_matrix_ptr == other.sparse_matrix_ptr &&
              current_row == other.current_row &&
              col_ptr == other.col_ptr &&
              val_ptr == other.val_ptr);
    }

    bool
    operator!=(const const_iterator& other) const
    { return !(*this == other); }
  };


  /**
   * Return a constant iterator to the start of the sparse matrix.
   */
  const_iterator
  begin() const
  { return (!empty()) ? begin(0) : end();  }

  /**
   * Return a constant iterator to the end of the sparse matrix.
   */
  const_iterator
  end() const
  { return {this}; }

  /**
   * Return a constant iterator to the start of row \p i.
   */
  const_iterator
  begin(const size_t i) const
  {
    Assert(i < rows, "Out of range error.");
    return {this, i};
  }

  /**
   * Return a constant iterator to the end of row \p i.
   */
  const_iterator
  end(const size_t i) const
  {
    Assert(i < rows, "Out of range error.");
    if (i + 1 == rows) return end();
    else return begin(i + 1);
  }


  /**
   * A struct defining a mutable row of the sparse matrix. This is used as
   * an interface to define range-based iterators over a particular row.
   */
  class row
  {
  private:
    SparseMatrix* sparse_matrix_ptr;
    const size_t row_num;

  public:
    row(SparseMatrix* sparse_matrix, const size_t i) :
      sparse_matrix_ptr(sparse_matrix), row_num(i) {}

    iterator
    begin()
    { return sparse_matrix_ptr->begin(row_num); }

    iterator
    end()
    { return sparse_matrix_ptr->end(row_num); }
  };


  /**
   * A convenience function which allows for mutable range-based iteration
   * over the specified row \p i.
   * \see SparseMatrix::const_row_iterator
   */
  row
  row_iterator(const size_t i)
  { return {this, i}; }

  /**
   * A struct defining a mutable row of the sparse matrix.
   * \see SparseMatrix::row
   */
  class const_row
  {
  private:
    const SparseMatrix* sparse_matrix_ptr;
    const size_t row_num;

  public:
    const_row(const SparseMatrix* sparse_matrix,
        const size_t i) :
        sparse_matrix_ptr(sparse_matrix), row_num(i) {}

    const_iterator
    begin()
    { return sparse_matrix_ptr->begin(row_num); }

    const_iterator
    end()
    { return sparse_matrix_ptr->end(row_num); }
  };


  /**
   * A convenience function which allows for constant range-based iteration
   * over the specified row \p i.
   * \see SparseMatrix::row_iterator
   */
  const_row
  const_row_iterator(const size_t i)
  { return {this, i}; }

  // @}
  /** \name Accessors */
  // @{

  /**
   * Return the column index located at relative position \p jr of row \p i.
   */
  const size_t&
  column(const size_t i, const size_t jr) const
  {
    Assert(i < rows, "Out of range error.");
    Assert(jr < row_length(i), "Relative index exceeds row length.");
    return colnums[i][jr];
  }


  /**
   * Read and write access to the element located at relative position \p jr
   * of row \p i.
   */
  value_type&
  value(const size_t i, const size_t jr)
  {
    Assert(i < rows, "Out of range error.");
    Assert(jr < row_length(i), "Relative index exceeds row length.");
    return coeffs[i][jr];
  }

  /**
   * Read access to the element located at relative position \p jr of row \p i.
   */
  const value_type&
  value(const size_t i, const size_t jr) const
  {
    Assert(i < rows, "Out of range error.");
    Assert(jr < row_length(i), "Relative index exceeds row length.");
    return coeffs[i][jr];
  }

  /**
   * Return a pointer to element <tt>(i, j)</tt>. If no element exists,
   * \p nullptr is returned.
   */
  value_type*
  locate(const size_t i, const size_t j)
  {
    Assert(i < rows && j < cols, "Out of range error.");
    for (size_t jr = 0; jr < row_length(i); ++jr)
      if (colnums[i][jr] == j)
        return &coeffs[i][jr];
    return nullptr;
  }

  /**
   * Return a constant pointer to element <tt>(i, j)</tt>. If no element
   * exists, \p nullptr is returned.
   */
  const value_type*
  locate(const size_t i, const size_t j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");
    for (uint64_t jr = 0; jr < row_length(i); ++jr)
      if (colnums[i][jr] == j)
        return &coeffs[i][jr];
    return nullptr;
  }


  /**
   * Read and write access to element <tt>(i, j)</tt>.
   * \throw If column \p j does not exist on row \p i.
   */
  value_type&
  operator()(const size_t i, const size_t j)
  {
    Assert(i < rows && j < cols, "Out of range error.");
    value_type* value_ptr = locate(i, j);
    Assert(value_ptr != nullptr, "Cannot access an uninitialized element.");
    return *value_ptr;
  }

  /**
   * Read access to element <tt>(i, j)</tt>.
   * \throw If column \p j does not exist on row \p i.
   */
  const value_type&
  operator()(const size_t i, const size_t j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");
    const value_type* value_ptr = locate(i, j);
    Assert(value_ptr != nullptr, "Cannot access an uninitialized element.");
    return *value_ptr;
  }

  /**
   * Return a pointer to the diagonal element of row \p i. If no diagonal
   * element exists, \p nullptr is returned.
   */
  value_type*
  diagonal(const size_t i)
  { return locate(i, i); }

  /**
   * Return a constant pointer to the diagonal element of row \p i. If no
   * diagonal element exists, \p nullptr is returned.
   */
  const value_type*
  diagonal(const size_t i) const
  { return locate(i, i); }

  /**
   * Return a vector of pointers to the diagonal elements of the sparse matrix.
   * If any particular diagonal is not present, its value is \p nullptr.
   */
   std::vector<value_type*>
   diagonal()
  {
     size_t min_dim = std::min(rows, cols);
     std::vector<value_type*> diag(min_dim);
     for (size_t i = 0; i < min_dim; ++i)
       diag[i] = diagonal(i);
     return diag;
  }

  /**
   * Return a vector of constant pointers to the diagonal elements of the
   * sparse matrix. If any particular diagonal is not present, its value is
   * \p nullptr.
   */
  std::vector<const value_type*>
  diagonal() const
  {
    size_t min_dim = std::min(rows, cols);
    std::vector<const value_type*> diag(min_dim);
    for (size_t i = 0; i < min_dim; ++i)
      diag[i] = diagonal(i);
    return diag;
  }

  // @}
  /** \name Modifiers */
  // @{

  /**
   * Return the sparse matrix to an uninitialized state.
   */
  void
  clear()
  {
    rows = cols = 0;
    colnums.clear();
    coeffs.clear();
  }

  /**
   * Reinitialize the sparse matrix with \p n_rows and \p n_cols with
   * \p default_row_length nonzero entries per row.
   */
  void
  reinit(const size_t n_rows,
         const size_t n_cols,
         const size_t default_row_length)
  {
    clear();
    rows = n_rows;
    cols = n_cols;
    colnums.resize(n_rows, std::vector<size_t>(default_row_length));
    coeffs.resize(n_rows, std::vector<value_type>(default_row_length));
  }

  /**
   * Reinitialize the sparse matrix with dimension \p n with
   * \p default_row_length nonzero entries per row.
   */
  void
  reinit(const size_t n, const size_t default_row_length)
  { reinit(n, n, default_row_length); }

  /**
   * Set element <tt>(i, j)</tt> to \p value. If the element is initialized,
   * override the value. If it is not, initialize it.
   */
  void
  set(const size_t i, const size_t j, const value_type value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

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

  /**
   * Add \p value to element <tt>(i, j)</tt>. If the element is not initialized,
   * the <code> set(const size_t i, const size_t j) </code> is called.
   */
  void
  add(const size_t i, const size_t j, const value_type value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // Locate the element
    value_type* value_ptr = locate(i, j);

    // Set if uninitialized, otherwise, add.
    if (locate(i, j) == nullptr) set(i, j, value);
    else *value_ptr += value;
  }

  /**
   * Swap the elements of two rows.
   */
  void
  swap_row(const size_t i, const size_t k)
  {
    Assert(i < rows && k < rows, "Out of range error.");
    colnums[i].swap(colnums[k]);
  }

  /**
   * Swap the elements of this sparse matrix with another.
   */
  void
  swap(SparseMatrix& other)
  {
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
    colnums.swap(other.colnums);
    coeffs.swap(other.coeffs);
  }

  // @}
  /** \name Scalar Operations */
  // @{

  /**
   * Negate the elements of the sparse matrix. This is computed via
   * \f$ \boldsymbol{A} = -\boldsymbol{A} = -a_{ij}, ~ \forall i, j \f$.
   */
  SparseMatrix&
  operator-()
  {
    for (entry elem : *this)
      elem.value = -elem.value;
    return *this;
  }

  /**
   * Multiply the elements of the sparse matrix by a scalar factor. This is
   * computed via
   * \f$ \boldsymbol{A} = \alpha \boldsymbol{A}
   *                    = \alpha a_{ij}, ~ \forall i, j
   * \f$.
   */
  SparseMatrix&
  operator*=(const value_type factor)
  {
    for (entry elem : *this)
      elem.value *= factor;
    return *this;
  }

  /**
   * Divide the elements of the sparse matrix by a scalar factor. This is
   * computed via
   * \f$ \boldsymbol{A} = \frac{1}{\alpha} \boldsymbol{A}
   *                    = \frac{\alpha}{a_{ij}}, ~ \forall i, j
   * \f$.
   */
  SparseMatrix&
  operator/=(const value_type factor)
  {
    Assert(factor != 0.0, "Zero division error.");

    for (entry elem : *this)
      elem.value /= factor;
    return *this;
  }

  // @}
  /** \name Linear Algebra */
  // @{

  /**
   * Add a sparse matrix multiplied by a scalar.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A} + \alpha \boldsymbol{B}
   *                    = \sum_{i=0}^{n} \sum_{j=0}^{n} a_{ij} + \alpha b_{ij},
   *                    ~ \forall i,j
   * \f$.
   */
  void
  add(const SparseMatrix& B,
      const value_type factor = 1.0)
  {
    Assert(colnums == B.colnums,
           "Adding two different sparse matrices with different "
           "sparsity patterns is not allowed.");

    std::vector<value_type>* row_ptr = &coeffs[0];
    const std::vector<value_type>* matrix_row_ptr = &B.coeffs[0];
    const std::vector<value_type>* const end_row_ptr = row_ptr + rows;

    while (row_ptr != end_row_ptr)
    {
      value_type* val_ptr = &(*row_ptr)[0];
      const value_type* matrix_ptr = &(*matrix_row_ptr)[0];
      const value_type* const end_ptr = val_ptr + row_ptr->size();

      while (val_ptr != end_ptr)
        *val_ptr++ += factor * *matrix_ptr++;

      row_ptr++; matrix_row_ptr++;
    }
  }

  /**
   * Compute a matrix-vector product. This is computed via
   * \f$ \vec{y} = \boldsymbol{A} \vec{x}
   *             = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i
   * \f$.
   */
  void
  vmult(const Vector& x, Vector& y,
        const bool adding = false) const
  {
    Assert(x.size() == cols, "Dimension mismatch error.");
    Assert(y.size() == rows, "Dimension mismatch error.");

    if (!adding) y = 0.0;
    for (const const_entry elem : *this)
      y[elem.row] += elem.value * x[elem.column];
  }

  /**
   * Return a matrix-vector product.
   * \see SparseMatrix::vmult
   */
  Vector
  vmult(const Vector& x) const
  {
    Vector y(rows);
    vmult(x, y);
    return y;
  }

  /**
   * Compute a matrix-vector product and add to the destination vector.
   * This is computed via
   * \f$ \vec{y} = \vec{y} + \boldsymbol{A} \vec{x}
   *             = y_i + \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i
   * \f$.
   */
  void
  vmult_add(const Vector& x, Vector& y)
  { vmult(x, y, true); }

  /**
   * Compute a transpose matrix-vector product. This is computed via
   * \f$ \vec{y} = \boldsymbol{A}^T \vec{x}
   *             = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall i
   * \f$.
   */
  void Tvmult(const Vector& x, Vector& y,
              const bool adding = false) const
  {
    Assert(x.size() == rows, "Dimension mismatch error.");
    Assert(y.size() == cols, "Dimension mismatch error.");

    if (!adding) y = 0.0;
    for (const const_entry elem : *this)
      y[elem.column] += elem.value * x[elem.row];
  }

  /**
   * Return a transpose matrix-vector product.
   * \see SparseMatrix::Tvmult
   */
  Vector
  Tvmult(const Vector& x)
  {
    Vector y(cols);
    Tvmult(x, y);
    return y;
  }

  /**
   * Compute a transpose matrix-vector product and to the destination vector.
   * This is computed via
   * \f[ \vec{y} = \vec{y} + \boldsymbol{A}^T \vec{x} \\
   *             = y_i + \sum_{i=1}^{n} a_{ji} x_i, ~ \forall i
   * \f]
   */
  void Tvmult_add(const Vector& x, Vector& y)
  { Tvmult(x, y, true); }

  /**
   * Compute a matrix-vector product.
   * \see SparseMatrix::vmult
   */
  Vector
  operator*(const Vector& x) const
  { return vmult(x); }

  // @}
  /** \name Printing Utilities */
  // @{

  std::string
  str(const bool scientific = false,
      const unsigned int precision = 3,
      const unsigned int width = 0) const
  {
    std::stringstream ss;
    print_formatted(ss, scientific, precision, width);
    return ss.str();
  }

  /**
   *  Print the sparse matrix as triplets of nonzero entries.
   *
   * \param os The output stream to print the matrix to.
   * \param scientific A flag for scientific notation.
   * \param precision The precision to display to sparse matrix elements.
   */
  void
  print(std::ostream& os = std::cout,
        const bool scientific = false,
        const unsigned int precision = 3) const
  {
    unsigned int w                   = precision + 7;
    std::ios::fmtflags old_flags     = os.flags();
    unsigned int       old_precision = os.precision(precision);

    if (scientific)
      os.setf(std::ios::scientific, std::ios::floatfield);
    else
      os.setf(std::ios::fixed, std::ios::floatfield);

    os.setf(std::ios::left, std::ios::adjustfield);
    os.width(w);

    os << "Row" << "Column" << "Value" << std::endl;
    os << "---" << "------" << "-----" << std::endl;

    for (const auto elem : *this)
      os << elem.row << elem.column << elem.value << std::endl;

    os.flags(old_flags);
    os.precision(old_precision);
  }

  /**
   * Print the matrix as if it were dense, filling in zero values.
   *
   * \param os The output stream to print the matrix to.
   * \param scientific A flag for scientific notation.
   * \param precision The precision to display to sparse matrix elements.
   * \param width The spacing between entries.
   */
  void
  print_formatted(std::ostream& os = std::cout,
                  const bool scientific = false,
                  const unsigned int precision = 3,
                  const unsigned int width = 0) const
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
    os.width(w);

    // Loop over rows and columns and print element-wise
    for (size_t i = 0; i < rows; ++i)
    {
      for (size_t j = 0; j < cols; ++j)
      {
        // Print the entry or a zero
        const value_type * entry = locate(i, j);
        os << ((!entry)? 0.0 : *entry);
      }
      os << std::endl;
    }
    os << std::endl;

    os.flags(old_flags);
    os.precision(old_precision);
  }

  // @}
};


/*-------------------- Inline Implementations --------------------*/


/**
 * Multiply a sparse matrix by a scalar factor.
 */
inline SparseMatrix
operator*(const double factor, const SparseMatrix& A)
{ return SparseMatrix(A) *= factor; }


/**
 * Multiply a sparse matrix by a scalar factor.
 */
inline SparseMatrix
operator*(const SparseMatrix& A, const double factor)
{ return SparseMatrix(A) *= factor; }


/**
 * Divide a sparse matrix by a scalar factor.
 */
inline SparseMatrix
operator/(const SparseMatrix& A, const double factor)
{ return SparseMatrix(A) /= factor; }


/**
 * Insert a sparse matrix into an output stream
 */
inline std::ostream&
operator<<(std::ostream& os, const SparseMatrix& A)
{ return os << A.str(); }

}
#endif //SPARSE_MATRIX_H

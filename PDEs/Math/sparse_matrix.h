#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "vector.h"
#include "matrix.h"
#include "exceptions.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cinttypes>


namespace math
{

template<typename number>
class SparseMatrixBase
{
public:
  virtual number&
  operator()(const uint64_t, const uint64_t) = 0;

  virtual const number&
  operator()(const uint64_t, const uint64_t) const = 0;
};


template<typename number>
class SparseMatrix : public SparseMatrixBase<number>
{
public:
  using value_type = number;
  using size_type  = uint64_t;

private:
  uint64_t rows;  ///< The number of rows.
  uint64_t cols;  ///< The number of columns.

  /// The row-wise nonzero column indices.
  std::vector<std::vector<size_type>> colnums;

  /// The row-wise nonzero data entries.
  std::vector<std::vector<value_type>> values;

public:

  class iterator
  {
  private:
    SparseMatrix* sparse_matrix_ptr;

    size_type current_row;
    size_type current_index;
    std::vector<size_type>::const_iterator col_ptr;
    typename std::vector<value_type>::iterator val_ptr;

    void
    advance()
    {
      // Increment the column and value pointer
      ++current_index; ++col_ptr; ++val_ptr;

      /* If the column pointer is at the end of the row, increment the row
       * and move the pointers to the start of the next row. */
      if (col_ptr == sparse_matrix_ptr->colnums[current_row].end())
      {
        ++current_row;

        /* If not the last row, move to the next row. If the last row, nothing
         * need be done since all attributes are equivalent to the stop
         * criteria. */
        if (current_row < sparse_matrix_ptr->n_rows())
        {
          current_index = 0;
          col_ptr = sparse_matrix_ptr->colnums[current_row].begin();
          val_ptr = sparse_matrix_ptr->values[current_row].begin();
        }
      }
    }

  public:
    iterator(SparseMatrix* sparse_matrix,
             const size_type row,
             const size_type index,
             std::vector<size_type>::const_iterator col,
             typename std::vector<value_type>::iterator val)
      : sparse_matrix_ptr(sparse_matrix),
        current_row(row), current_index(index),
        col_ptr(col), val_ptr(val)
    {}

    iterator&
    operator++()
    { advance(); return *this; }

    iterator
    operator++(int)
    { auto it = *this; advance(); return *this; }

    bool
    operator==(const iterator& other) const
    {
      return (sparse_matrix_ptr == other.sparse_matrix_ptr &&
              current_row       == other.current_row       &&
              current_index     == other.current_index     &&
              col_ptr           == other.col_ptr           &&
              val_ptr           == other.val_ptr);
    }

    bool
    operator!=(const iterator& other) const
    { return !(*this == other); }

    bool
    operator<(const iterator& other) const
    {
      if (current_row < other.current_row) return true;
      if (current_row == other.current_row)
        if (*col_ptr < *other.col_ptr) return true;
      return false;
    }

    bool
    operator<=(const iterator& other) const
    {
      if (*this == other) return true;
      return (*this < other);
    }

    bool
    operator>(const iterator& other) const
    { return !(*this <= other); }

    bool
    operator>=(const iterator& other) const
    { return !(*this < other); }

    iterator&
    operator*()
    { return *this; }

    iterator
    operator->()
    { return this; }

    const size_type&
    row() const
    { return current_row; }

    const size_type&
    index() const
    { return current_index; }

    const size_type&
    column() const
    { return *col_ptr; }

    value_type&
    value() const
    { return *val_ptr; }
  };


  class const_iterator
  {
  private:
    const SparseMatrix* sparse_matrix_ptr;

    size_type current_row;
    size_type current_index;
    std::vector<size_type>::const_iterator col_ptr;
    typename std::vector<value_type>::const_iterator val_ptr;

    void
    advance()
    {
      // Increment the column and value pointer
      ++current_index; ++col_ptr; ++val_ptr;

      /* If the column pointer is at the end of the row, increment the row
       * and move the pointers to the start of the next row. */
      if (col_ptr == sparse_matrix_ptr->colnums[current_row].end())
      {
        ++current_row;

        /* If not the last row, move to the next row. If the last row, nothing
         * need be done since all attributes are equivalent to the stop
         * criteria. */
        if (current_row < sparse_matrix_ptr->n_rows())
        {
          current_index = 0;
          col_ptr = sparse_matrix_ptr->colnums[current_row].begin();
          val_ptr = sparse_matrix_ptr->values[current_row].begin();
        }
      }
    }

  public:
    const_iterator(const SparseMatrix* sparse_matrix,
                   const size_type row,
                   const size_type index,
                   std::vector<size_type>::const_iterator col,
                   typename std::vector<value_type>::const_iterator val)
        : sparse_matrix_ptr(sparse_matrix),
          current_row(row), current_index(index),
          col_ptr(col), val_ptr(val)
    {}

    const_iterator&
    operator++()
    { advance(); return *this; }

    const_iterator
    operator++(int)
    { auto it = *this; advance(); return *this; }

    bool
    operator==(const const_iterator& other) const
    {
      return (sparse_matrix_ptr == other.sparse_matrix_ptr &&
              current_row       == other.current_row       &&
              current_index     == other.current_index     &&
              col_ptr           == other.col_ptr           &&
              val_ptr           == other.val_ptr);
    }

    bool
    operator!=(const const_iterator& other) const
    { return !(*this == other); }

    bool
    operator<(const const_iterator& other) const
    {
      if (current_row < other.current_row) return true;
      if (current_row == other.current_row)
        if (*col_ptr < *other.col_ptr) return true;
      return false;
    }

    bool
    operator<=(const const_iterator& other) const
    {
      if (*this == other) return true;
      return (*this < other);
    }

    bool
    operator>(const const_iterator& other) const
    { return !(*this <= other); }

    bool
    operator>=(const const_iterator& other) const
    { return !(*this < other); }

    const_iterator&
    operator*()
    { return *this; }

    const_iterator
    operator->()
    { return this; }

    const size_type&
    row() const
    { return current_row; }

    const size_type&
    index() const
    { return current_index; }

    const size_type&
    column() const
    { return *col_ptr; }

    const value_type&
    value() const
    { return *val_ptr; }
  };


public:
  /// Default constructor.
  SparseMatrix()
    : rows(0), cols(0), colnums(), values()
  {}

  /// Copy constructor.
  SparseMatrix(const SparseMatrix& other)
    : rows(other.rows), cols(other.cols),
      colnums(other.colnums), values(other.values)
  {}

  /// Move constructor.
  SparseMatrix(SparseMatrix&& other)
    : rows(other.rows), cols(other.cols),
      colnums(std::move(other.colnums)),
      values(std::move(other.values))
  {}

  /**
   * Construct a sparse matrix with \p n_rows and \p n_cols with
   * \p default_row_length nonzero entries per row.
   */
  explicit
  SparseMatrix(const size_type n_rows,
               const size_type n_cols,
               const size_type default_row_length)
      : rows(n_rows), cols(n_cols),
        colnums(n_rows, std::vector<size_type>(default_row_length)),
        values(n_rows, std::vector<value_type>(default_row_length))
  {}

  /**
   * Construct a square sparse matrix with dimension \p n with
   * \p default_row_length nonzero entries per row.
   */
  explicit
  SparseMatrix(const size_type n,
               const size_type default_row_length)
      : SparseMatrix(n, n, default_row_length)
  {}

  /**
   * Construct a sparse matrix from a sparsity pattern. In this context, a
   * sparsity pattern is defined by a list of lists. Each inner list stores
   * the nonzero column indices for a given row.
   */
  SparseMatrix(std::vector<std::vector<size_type>> sparsity_pattern)
      : rows(sparsity_pattern.size()), colnums(sparsity_pattern),
        values(sparsity_pattern.size())
  {
    cols = 0;
    for (size_type i = 0; i < rows; ++i)
    {
      // Sort column indices and allocate data
      std::sort(colnums[i].begin(), colnums[i].end());
      values[i].resize(colnums[i].size());

      // Set number of columns
      for (const auto& col : colnums[i])
        cols = (col > cols) ? col : cols;
    }
  }

  /// Construct from a dense matrix.
  SparseMatrix(const Matrix<number>& other)
    : rows(other.n_rows()), cols(other.n_cols()),
      colnums(other.n_rows()), values(other.n_rows())
  {
    for (size_type i = 0; i < rows; ++i)
      for (size_type j = 0; j < cols; ++j)
        if (other[i][j] != 0.0)
        {
          colnums[i].push_back(j);
          values[i].push_back(other[i][j]);
        }
  }

  /// Assignment operator.
  SparseMatrix&
  operator=(const SparseMatrix& other)
  {
    rows    = other.rows;
    cols    = other.cols;
    colnums = other.colnums;
    values  = other.values;
    return *this;
  }

  /// Assignment with a full matrix.
  SparseMatrix&
  operator=(const Matrix<value_type>& other)
  {
    colnums.clear();
    values.clear();

    rows = other.n_rows();
    cols = other.n_cols();
    colnums.resize(rows);
    values.resize(rows);
    for (size_type i = 0; i < rows; ++i)
      for (size_type j = 0; j < cols; ++j)
        if (other[i][j] != 0.0)
        {
          colnums[i].push_back(j);
          values[i].push_back(other[i][j]);
        }
    return *this;
  }

  /// Assignment with a scalar value.
  SparseMatrix&
  operator=(const value_type value)
  {
    Assert(!empty(), "Cannot set an empty matrix to a scalar.");
    for (iterator& elem : *this)
      elem.value() = value;
  }

  /// Equality comparison operator.
  bool
  operator==(const SparseMatrix& other) const
  {
    return (rows    == other.rows &&
            cols    == other.cols &&
            colnums == other.colnums &&
            values  == other.values);
  }

  /// Inequality comparison operator.
  bool
  operator!=(const SparseMatrix& other) const
  { return !(*this == other); }

  /** \name Information */
  // @{

  /// Return the number of rows.
  size_type
  n_rows() const
  { return rows; }

  /// Return the number of columns.
  size_type
  n_cols() const
  { return cols; }

  /// Return the number of nonzero elements.
  size_type
  nnz() const
  {
    size_type count = 0;
    for (size_type i = 0; i < rows; ++i)
      count += row_length(i);
    return count;
  }

  /// Return the number of nonzero entries in row \p i.
  size_type
  row_length(const size_type i) const
  {
    Assert(i < rows, "Out of range error.");
    return colnums[i].size();
  }

  /// Return the column index for nonzero entry \p jr of row \p i.
  size_type
  column(const size_type i, const size_type jr) const
  {
    Assert(i < rows, "Out of range error.");
    Assert(jr < row_length(i), "Relative index exceeds row length.");
    return colnums[i][jr];
  }

  /**
   * Return a pointer to the element at index <tt>(i, j)</tt>.
   * If no element exists null is returned.
   */
  value_type*
  locate(const size_type i, const size_type j)
  {
    Assert(i < rows && j < cols, "Out of range error.");
    for (size_type jr = 0; jr < row_length(i); ++jr)
      if (colnums[i][jr] == j)
        return &values[i][jr];
    return nullptr;
  }

  /**
   * Return a constant pointer to the element at index \p(i, \p j). If no
   * element exists, null is returned.
   */
  const value_type*
  locate(const size_type i, const size_type j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");
    for (uint64_t jr = 0; jr < row_length(i); ++jr)
      if (colnums[i][jr] == j)
        return &values[i][jr];
    return nullptr;
  }

  /// Return whether the SparseMatrix is empty.
  bool
  empty() const
  {
    return (rows == 0 && cols == 0 &&
            colnums.empty() && values.empty());
  }

  // @}
  /** \name Iterators */
  // @{

  /// Mutable iterator to the first row of the SparseMatrix.
  iterator
  begin()
  { return {this, 0, 0, colnums[0].begin(), values[0].begin()};  }

  /// Constant iterator to the first row of the SparseMatrix.
  const_iterator
  begin() const
  { return {this, 0, 0, colnums[0].begin(), values[0].begin()};  }

  /// Mutable iterator to the end of the SparseMatrix.
  iterator
  end()
  {
    return {this, rows, colnums[rows - 1].size(),
            colnums[rows - 1].end(), values[rows - 1].end()};
  }

  /// Constant iterator to the end of the SparseMatrix.
  const_iterator
  end() const
  {
    return {this, rows, colnums[rows - 1].size(),
            colnums[rows - 1].end(), values[rows - 1].end()};
  }

  /// Mutable iterator to the first entry of row \p i.
  iterator
  begin(const size_type i)
  {
    Assert(i < rows && !empty(), "Out of range error.");
    return {this, i, 0, colnums[i].begin(), values[i].begin()};
  }

  /// Constant iterator to the first entry of row \p i.
  const_iterator
  begin(const size_type i) const
  {
    Assert(i < rows, "Out of range error.");
    if (empty()) return end();
    return {this, i, 0, colnums[i].begin(), values[i].begin()};
  }

  /// Mutable iterator to the end of row \p i.
  iterator
  end(const size_type i)
  {
    Assert(i < rows, "Out of range error.");
    if (i + 1) return end();
    return {this, i, colnums[i].size(), colnums[i].end(), values[i].end()};
  }

  /// Constant iterator to the end of row \p i.
  const_iterator
  end(const size_type i) const
  {
    Assert(i < rows, "Out of range error.");
    if (i + 1 == rows) return end();
    else return {this, i, colnums[i].size(),
                 colnums[i].end(), values[i].end()};
  }

  // @}
  /** \name Accessors */
  // @{

  /// Read/write access to element <tt>(i, j)</tt>.
  value_type&
  operator()(const size_type i, const size_type j) override
  {
    Assert(i < rows && j < cols, "Out of range error.");
    value_type* value_ptr = locate(i, j);
    Assert(value_ptr != nullptr, "Cannot access an uninitialized element.");
    return *value_ptr;
  }

  /// Read access to element <tt>(i, j)</tt>.
  const value_type&
  operator()(const size_type i, const size_type j) const override
  {
    Assert(i < rows && j < cols, "Out of range error.");
    const value_type* value_ptr = locate(i, j);
    Assert(value_ptr != nullptr, "Cannot access an uninitialized element.");
    return *value_ptr;
  }

  Vector<value_type>
  diagonal() const
  {
    Assert(!empty(), "Empty matrix error.");

    Vector<value_type> x(std::min(rows, cols), 0.0);
    for (const_iterator& elem : *this)
      if (elem.row() == elem.column())
        x[elem.row()] = elem.value();
    return x;
  }

  // @}
  /** \name Modifiers */
  // @{

  /// Return the SparseMatrix to an uninitialized state.
  void
  clear()
  {
    rows = cols = 0;
    colnums.clear();
    values.clear();
  }

  /**
   * Reinit the sparse matrix with \p n_rows and \p n_cols with
   * \p default_row_length nonzero entries per row.
   */
  void
  reinit(const size_type n_rows,
         const size_type n_cols,
         const size_type default_row_length)
  {
    clear();
    rows = n_rows;
    cols = n_cols;
    colnums.resize(n_rows, std::vector<size_type>(default_row_length));
    values.resize(n_rows, std::vector<value_type>(default_row_length));
  }

  /**
   * Reinit the sparse matrix with dimension \p n with \p default_row_length
   * nonzero entries per row.
   */
  void
  reinit(const size_type n, const size_type default_row_length)
  { reinit(n, n, default_row_length); }

  /**
   * Set element <tt>(i, j)</tt> to \p value. If the element is initialized,
   * override the value. If it is not, initialize it.
   */
  void
  set(const size_type i, const size_type j, const value_type value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    /* If the row is empty or the column number is larger than all current
     * entries on the row, add to the back of the row. */
    if (colnums[i].size() == 0 or colnums[i].back() < j)
    {
      colnums[i].push_back(j);
      values[i].push_back(value);
      return;
    }

    // Find the index to insert the entry which maintains sorting.
    auto it = std::lower_bound(colnums[i].begin(), colnums[i].end(), j);
    uint64_t jr = it - colnums[i].begin();

    // If this points to an existing column, override the value
    if (*it == j)
    {
      values[i][jr] = value;
      return;
    }

    // Insert the entries into the data structures maintaining sorting
    colnums[i].insert(it, j);
    values[i].insert(values[i].begin() + jr, value);
  }

  /**
   * Add \p value to element <tt>(i, j)</tt>. If the element is not initialized,
   * the <code> set(const size_type i, const size_type j) </code> is called.
   */
  void
  add(const size_type i, const size_type j, const value_type value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // Locate the element
    number* value_ptr = locate(i, j);

    // Set if uninitialized, otherwise, add.
    if (locate(i, j) == nullptr) set(i, j, value);
    else *value_ptr += value;
  }

  /// Swap the elements of two rows.
  void
  swap_row(const size_type i, const size_type k)
  {
    Assert(i < rows && k < rows, "Out of range error.");
    colnums[i].swap(colnums[k]);
  }

  /// Swap the elements of this SparseMatrix with another.
  void
  swap(SparseMatrix& other)
  {
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
    colnums.swap(other.colnums);
    values.swap(other.values);
  }

  // @}
  /** \name Scalar Operations */
  // @{

  /// Element-wise negation in-place. */
  SparseMatrix&
  operator-()
  {
    for (iterator& elem : *this)
      elem.value() = -elem.value();
    return *this;
  }

  /// Element-wise multiplication by a scalar in-place.
  SparseMatrix&
  operator*=(const value_type value)
  {
    for (iterator& elem : *this)
      elem.value() *= value;
    return *this;
  }

  /// Element-wise multiplication by a scalar in-place.
  SparseMatrix&
  operator/=(const value_type value)
  {
    Assert(value != 0.0, "Zero division error.");

    for (iterator& elem : *this)
      elem.value() /= value;
    return *this;
  }

  // @}
  /** \name Linear Algebra */
  // @{

  /**
   * Add another SparseMatrix multiplied by a scalar to this in-place.
   *
   * This is computed via
   * \f[ \boldsymbol{C} = \boldsymbol{A} + \alpha \boldsymbol{B} \\
   *     c_{ik} = \sum_{i=0}^{n} \sum_{j=0}^{n} a_{ij} + \alpha b_{ij}.
   * \f]
   */
  void
  add(const SparseMatrix& B,
      const value_type factor = 1.0)
  {
    Assert(colnums == B.colnums,
           "Adding two different sparse matrices with different "
           "sparsity patterns is not allowed.");

    std::vector<value_type>* row_ptr = &values[0];
    const std::vector<value_type>* matrix_row_ptr = &B.values[0];
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
   * Compute a matrix-vector product.
   *
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  void
  vmult(const Vector<value_type>& x,
        Vector<value_type>& y,
        const bool adding = false) const
  {
    Assert(x.size() == cols, "Dimension mismatch error.");
    Assert(y.size() == rows, "Dimension mismatch error.");

    if (!adding) y = 0.0;
    for (const_iterator& elem : *this)
      y[elem.row()] += elem.value() * x[elem.column()];
  }

  /// See \ref vmult.
  Vector<value_type>
  vmult(const Vector<value_type>& x) const
  {
    Vector<value_type> y(rows);
    vmult(x, y);
    return y;
  }

  /**
   * Compute a matrix-vector product and add to the destination vector.
   *
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  void
  vmult_add(const Vector<value_type>& x,
            Vector<value_type>& y)
  { vmult(x, y, true); }

  /**
   * Compute a transpose matrix-vector product.
   *
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A}^T \vec{x} \\
   *     y_i = \sum_{i=1}^{n} a_{ji} x_i, \hspace{0.25cm} \forall i
   * \f]
   */
  void Tvmult(const Vector<value_type>& x,
              Vector<value_type>& y,
              const bool adding = false) const
  {
    Assert(x.size() == rows, "Dimension mismatch error.");
    Assert(y.size() == cols, "Dimension mismatch error.");

    if (!adding) y = 0.0;
    for (const_iterator& elem : *this)
      y[elem.column()] += elem.value() * x[elem.row()];
  }

  /// See \ref Tvmult.
  Vector<value_type>
  Tvmult(const Vector<value_type>& x)
  {
    Vector<value_type> y(cols);
    Tvmult(x, y);
    return y;
  }

  /**
   * Compute a transpose matrix-vector product and to the destination vector.
   *
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A}^T \vec{x} \\
   *     y_i = \sum_{i=1}^{n} a_{ji} x_i, \hspace{0.25cm} \forall i
   * \f]
   */
  void Tvmult_add(const Vector<value_type>& x,
                  Vector<value_type>& y)
  { Tvmult(x, y, true); }

  /// Compute a matrix-vector product. \see vmult
  Vector<value_type>
  operator*(const Vector<value_type>& x) const
  {
    Assert(cols == x.size(), "Dimension mismatch error.");

    Vector<value_type> y(rows);
    for (const_iterator& elem : *this)
      y[elem.row()] += elem.value() * x[elem.column()];
    return y;
  }

  // @}
  /** \name Printing Utilities */
  // @{

  /// Print the sparse matrix as triplets of nonzero entries.
  void
  print(std::ostream& os = std::cout) const
  {
    os.setf(std::ios::left, std::ios::adjustfield);

    os.fill(' '); os.precision(6);
    os << std::setw(8) << "Row"
       << std::setw(8) << "Column"
       << std::setw(10) << "Value" << std::endl;
    os << std::setw(8) << "---"
       << std::setw(8) << "------"
       << std::setw(10) << "-----" << std::endl;

    for (const_iterator& elem : *this)
        os << std::setw(8)  << elem.row()
           << std::setw(8)  << elem.column()
           << std::setw(10) << elem.value() << std::endl;
  }

  /// Print the sparse matrix as a normal matrix with zero fill ins.
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

    // Loop over rows and columns and print element-wise
    for (uint64_t i = 0; i < rows; ++i)
    {
      for (uint64_t j = 0; j < cols; ++j)
      {
        // Print the entry or a zero
        const number* entry = locate(i, j);
        os << std::setw(w) << ((!entry)? 0.0 : *entry);
      }
      os << std::endl;
    }
  }

  // @}
};


/*-------------------- Inline Implementations --------------------*/


/** Element-wise multiplication by a scalar value. */
template<typename number>
inline SparseMatrix<number>
operator*(const number value, const SparseMatrix<number>& A)
{
  return A * value;
}


}
#endif //SPARSE_MATRIX_H

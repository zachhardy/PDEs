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
class SparseMatrix
{
private:
  uint64_t rows;  ///< The number of rows.
  uint64_t cols;  ///< The number of columns.

  /** The row-wise nonzero column indices. */
  std::vector<std::vector<uint64_t>> colnums;

  /** The row-wise nonzero data entries. */
  std::vector<std::vector<number>> data;

public:
  /** Default constructor. */
  SparseMatrix()
      : rows(0), cols(0), colnums(), data()
  {}

  /** Construct a sparse matrix with \p n_rows and \p n_cols with
   *  \p default_row_length nonzero entries per row. */
  explicit SparseMatrix(const uint64_t n_rows,
                        const uint64_t n_cols,
                        const uint64_t default_row_length)
      : rows(n_rows), cols(n_cols),
        colnums(n_rows, std::vector<uint64_t>(default_row_length)),
        data(n_rows, std::vector<number>(default_row_length))
  {}

  /** Construct a square sparse matrix with dimension \p n with
   *  \p default_row_length nonzero entries per row. */
  explicit SparseMatrix(const uint64_t n,
                        const uint64_t default_row_length)
      : SparseMatrix(n, n, default_row_length)
  {}

  /** Construct a sparse matrix from a sparsity pattern. In this context, a
   *  sparsity pattern is defined by a list of lists. Each inner list stores
   *  the nonzero column indices for a given row. */
  SparseMatrix(std::vector<std::vector<uint64_t>> sparsity_pattern)
      : rows(sparsity_pattern.size()), colnums(sparsity_pattern),
        data(sparsity_pattern.size())
  {
    cols = 0;
    for (uint64_t i = 0; i < rows; ++i)
    {
      // Sort column indices and allocate data
      std::sort(colnums[i].begin(), colnums[i].end());
      data[i].resize(colnums[i].size());

      // Set number of columns
      for (const auto& col : colnums[i])
        cols = (col > cols) ? col : cols;
    }
  }

  /** Copy constructor. */
  SparseMatrix(const SparseMatrix& other)
      : rows(other.rows), cols(other.cols),
        colnums(other.colnums), data(other.data)
  {}

  /** Move constructor. */
  SparseMatrix(SparseMatrix&& other)
      : rows(other.rows), cols(other.cols),
        colnums(std::move(other.colnums)),
        data(std::move(other.data))
  {}

  /** Construct from a dense matrix. */
  SparseMatrix(const Matrix <number>& other)
      : rows(other.n_rows()), cols(other.n_cols()),
        colnums(other.n_rows()), data(other.n_rows())
  {
    for (uint64_t i = 0; i < rows; ++i)
      for (uint64_t j = 0; j < cols; ++j)
        if (other[i][j] != 0.0)
        {
          colnums[i].push_back(j);
          data[i].push_back(other[i][j]);
        }
  }

  /** Assignment operator. */
  SparseMatrix& operator=(const SparseMatrix& other)
  {
    rows = other.rows;
    cols = other.cols;
    colnums = other.colnums;
    data = other.data;
    return *this;
  }

  /** Equality operator. */
  bool operator==(const SparseMatrix& other) const
  {
    return (rows == other.rows && cols == other.cols &&
            colnums == other.colnums && data == other.data);
  }

  /** Inequality operator. */
  bool operator!=(const SparseMatrix& other) const
  {
    return !(*this == other);
  }

  //###########################################################################
  /** \name Information */
  // @{

  /** Return the number of rows. */
  uint64_t n_rows() const
  { return rows; }

  /** Return the number of columns. */
  uint64_t n_cols() const
  { return cols; }

  /** Return the number of nonzero elements. */
  uint64_t nnz() const
  {
    uint64_t count = 0;
    for (uint64_t i = 0; i < rows; ++i)
      count += row_length(i);
    return count;
  }

  /** Return the number of nonzero entries in row \p i. */
  uint64_t row_length(const uint64_t i) const
  {
    Assert(i < rows, "Out of range error.");
    return colnums[i].size();
  }

  /** Return the column index for nonzero entry \p jr of row \p i. */
  uint64_t column(const uint64_t i, const uint64_t jr) const
  {
    Assert(i < rows, "Out of range error.");
    Assert(jr < row_length(i), "Relative index exceeds row length.");
    return colnums[i][jr];
  }

  /** Return a pointer to the element at index \p(i, \p j). */
  number* locate(const uint64_t i, const uint64_t j)
  {
    Assert(i < rows && j < cols, "Out of range error.");
    for (uint64_t jr = 0; jr < row_length(i); ++jr)
      if (colnums[i][jr] == j)
        return &data[i][jr];
    return nullptr;
  }

  /** Return a constant pointer to the element at index \p(i, \p j). */
  const number* locate(const uint64_t i, const uint64_t j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");
    for (uint64_t jr = 0; jr < row_length(i); ++jr)
      if (colnums[i][jr] == j)
        return &data[i][jr];
    return nullptr;
  }

  /** Return whether the sparse matrix is empty. */
  bool empty() const
  {
    return (rows == 0 && cols == 0 &&
            colnums.empty() && data.empty());
  }

  // @}

  //###########################################################################
  /** \name Iterators */
  // @{

  /** Mutable iterator to the elements of the sparse matrix. */
  class iterator
  {
  private:
    SparseMatrix*   ref_sparse_matrix;
    uint64_t        current_row;
    uint64_t        current_index;

    void advance()
    {
      Assert(current_row < ref_sparse_matrix->n_rows(), "Out of range error.");

      // Increment the current index on the current row.
      ++current_index;

      /* If the current index is at the end of the row, move to the next row
       * or set the attributes to those of the invalid iterator. */
      if (current_index == ref_sparse_matrix->row_length(current_row))
      {
        if (current_row + 1 < ref_sparse_matrix->n_rows())
        {
          ++current_row;
          current_index = 0;
        }
        else
        {
          current_row = -1;
          current_index = -1;
        }
      }
    }

  public:
    iterator(SparseMatrix* sparse_matrix,
             const uint64_t row,
             const uint64_t index)
        : ref_sparse_matrix(sparse_matrix),
          current_row(row), current_index(index)
    {}

    iterator(SparseMatrix* sparse_matrix)
        : ref_sparse_matrix(sparse_matrix),
          current_row(-1), current_index(-1)
    {}

    iterator& operator++() { advance(); return *this; }
    iterator operator++(int) { auto it = *this; advance(); return it; }

    bool operator==(const iterator& other) const
    {
      return (ref_sparse_matrix == other.ref_sparse_matrix &&
              current_row == other.current_row &&
              current_index == other.current_index);
    }

    bool operator!=(const iterator& other) const
    { return !(*this == other); }

    iterator& operator*() { return *this; }

    uint64_t row() const { return current_row; }
    uint64_t index() const { return current_index; }
    uint64_t column() const
    { return ref_sparse_matrix->colnums[current_row][current_index]; }

    number& value()
    { return ref_sparse_matrix->data[current_row][current_index]; }
  };


  /** Mutable iterator to the start of row \p i. */
  iterator begin(const uint64_t i)
  {
    Assert(i < rows, "Out of range error.");

    // If empty, return the invalid end iterator.
    if (rows == 0)
      return {this};

    // Increment the row until a non-empty row is encountered.
    uint64_t row = i;
    while (row < rows && colnums[row].size() == 0)
      ++row;

    // Return the start of a row, or the invalid end iterator.
    if (row == rows) return {this};
    else return {this, row, 0};
  }

  /** Mutable iterator to the end of row \p i. */
  iterator end(const uint64_t i)
  {
    Assert(i < rows, "Out of range error.");
    if (i + 1 == rows) return {this};
    else return begin(i + 1);
  }

  /** Mutable iterator to the start of the sparse matrix. */
  iterator begin() { return begin(0); }

  /** Mutable iterator to the end of the sparse matrix. */
  iterator end() { return {this}; }


  /** Constant iterator over the elements of the sparse matrix. */
  class const_iterator
  {
  private:
    const SparseMatrix*   ref_sparse_matrix;
    uint64_t        current_row;
    uint64_t        current_index;

    void advance()
    {
      Assert(current_row < ref_sparse_matrix->n_rows(), "Out of range error.");

      // Increment the current index on the current row.
      ++current_index;

      /* If the current index is at the end of the row, move to the next row
       * or set the attributes to those of the invalid iterator. */
      if (current_index == ref_sparse_matrix->row_length(current_row))
      {
        if (current_row + 1 < ref_sparse_matrix->n_rows())
        {
          current_row += 1;
          current_index = 0;
        }
        else
        {
          current_row = -1;
          current_index = -1;
        }
      }
    }

  public:
    const_iterator(const SparseMatrix* sparse_matrix,
                   const uint64_t row,
                   const uint64_t index)
        : ref_sparse_matrix(sparse_matrix),
          current_row(row), current_index(index)
    {}

    const_iterator(const SparseMatrix* sparse_matrix)
        : ref_sparse_matrix(sparse_matrix),
          current_row(-1), current_index(-1)
    {}

    const_iterator& operator++() { advance(); return *this; }
    const_iterator operator++(int) { auto it = *this; advance(); return it; }

    bool operator==(const const_iterator& other) const
    {
      return (ref_sparse_matrix == other.ref_sparse_matrix &&
              current_row == other.current_row &&
              current_index == other.current_index);
    }

    bool operator!=(const const_iterator& other) const
    { return !(*this == other); }

    const_iterator& operator*() { return *this; }

    uint64_t row() const { return current_row; }
    uint64_t index() const { return current_index; }
    uint64_t column() const
    { return ref_sparse_matrix->colnums[current_row][current_index]; }

    number value() const
    { return ref_sparse_matrix->data[current_row][current_index]; }
  };


  /** Constant iterator to the start of row \p i. */
  const_iterator begin(const uint64_t i) const
  {
    Assert(i < rows, "Out of range error.");

    // If empty, return the invalid end iterator.
    if (rows == 0)
      return {this};

    // Increment the row until a non-empty row is encountered.
    uint64_t row = i;
    while (row < rows && colnums[row].size() == 0)
      ++row;

    // Return the start of a row, or the invalid end iterator.
    if (row == rows) return {this};
    else return {this, row, 0};
  }

  /** Constant iterator to the end of row \p i. */
  const_iterator end(const uint64_t i) const
  {
    Assert(i < rows, "Out of range error.");
    if (i + 1 == rows) return {this};
    else return begin(i + 1);
  }

  /** Constant iterator to the start of the sparse matrix. */
  const_iterator begin() const { return begin(0); }

  /** Constant iterator to the end of the sparse matrix. */
  const_iterator end() const { return {this}; }

  // @}

  //###########################################################################
  /** \name Accessors */
  // @{

  /** Read access to element \p(i, \p j). */
  number operator()(const uint64_t i, const uint64_t j) const
  {
    Assert(i < rows && j < cols, "Out of range error.");
    number* value_ptr = locate(i, j);
    Assert(value_ptr != nullptr, "Cannot access an uninitialized element.");
    return *value_ptr;
  }

  /** Read/write access to element \p(i, \p j). */
  number& operator()(const uint64_t i, const uint64_t j)
  {
    Assert(i < rows && j < cols, "Out of range error.");
    number* value_ptr = locate(i, j);
    Assert(value_ptr != nullptr, "Cannot access an uninitialized element.");
    return *value_ptr;
  }

  // @}

  //###########################################################################
  /** \name Modifiers */
  // @{

  /** Clear all data from the sparse matrix. */
  void clear()
  {
    rows = cols = 0;
    colnums.clear();
    data.clear();
  }

  /** Reinit the sparse matrix with \p n_rows and \p n_cols with
   *  \p default_row_length nonzero entries per row. */
  void reinit(const uint64_t n_rows, const uint64_t n_cols,
              const uint64_t default_row_length)
  {
    clear();
    rows = n_rows;
    cols = n_cols;
    colnums.resize(n_rows, std::vector<uint64_t>(default_row_length));
    data.resize(n_rows, std::vector<number>(default_row_length));
  }

  /** Reinit the sparse matrix with dimension \p n with \p default_row_length
   *  nonzero entries per row. */
  void reinit(const uint64_t n, const uint64_t default_row_length)
  {
    reinit(n, n, default_row_length);
  }

  /** Set the element at row \p i and column \p j to \p value. If the element
   *  is initialized, override the value. If it is not, initialize it. */
  void set(const uint64_t i, const uint64_t j, const number value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    /* If the row is empty or the column number is larger than all current
     * entries on the row, add to the back of the row. */
    if (colnums[i].size() == 0 or colnums[i].back() < j)
    {
      colnums[i].push_back(j);
      data[i].push_back(value);
      return;
    }

    // Find the index to insert the entry which maintains sorting.
    auto it = std::lower_bound(colnums[i].begin(), colnums[i].end(), j);
    uint64_t jr = it - colnums[i].begin();

    // If this points to an existing column, override the value
    if (*it == j)
    {
      data[i][jr] = value;
      return;
    }

    // Insert the entries into the data structures maintaining sorting
    colnums[i].insert(it, j);
    data[i].insert(data[i].begin() + jr, value);
  }

  /** Add \p value to the element at row \p i and column \p j. If the element
   *  is not initialized, the <code> set(const uint64_t i, const uint64_t j)
   *  </code> is called. */
  void add(const uint64_t i, const uint64_t j, const number value)
  {
    Assert(i < rows && j < cols, "Out of range error.");

    // Locate the element
    number* value_ptr = locate(i, j);

    // Set if uninitialized, otherwise, add.
    if (locate(i, j) == nullptr) set(i, j, value);
    else *value_ptr += value;
  }

  /** Swap the elements of two rows. */
  void swap_row(const uint64_t i0, const uint64_t i1)
  {
    Assert(i0 < rows && i1 < rows, "Out of range error.");
    colnums[i0].swap(colnums[i1]);
  }

  /** Swap the elements of this sparse matrix with another. */
  void swap(SparseMatrix& other)
  {
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
    colnums.swap(other.colnums);
    data.swap(other.data);
  }

  // @}

  //###########################################################################
  /** \name Scalar Operations */
  // @{

  /** Element-wise multiplication by a scalar. */
  SparseMatrix operator*(const number value) const
  {
    SparseMatrix A(*this);
    for (iterator& entry : A)
      entry.value() *= value;
    return A;
  }

  /** Element-wise multiplication by a scalar in-place. */
  SparseMatrix& operator*=(const number value)
  {
    for (iterator& entry : *this)
      entry.value() *= value;
    return *this;
  }

  /** Element-wise multiplication by a scalar. */
  SparseMatrix operator/(const number value) const
  {
    Assert(value != 0.0, "Zero division error.");
    SparseMatrix A(*this);
    for (iterator& entry : A)
      entry.value() /= value;
    return A;
  }

  /** Element-wise multiplication by a scalar in-place. */
  SparseMatrix& operator/=(const number value) const
  {
    Assert(value != 0.0, "Zero division error.");
    for (iterator& entry : *this)
      entry.value() /= value;
    return *this;
  }

  // @}

  //###########################################################################
  /** \name Linear Algebra */
  // @{

  /** Element-wise addition of two sparse matrices. */
  SparseMatrix operator+(const SparseMatrix& other) const
  {
    Assert(rows == other.rows && cols == other.cols,
           "Dimension mismatch error.");
    SparseMatrix A(*this);
    for (const_iterator& entry : other)
      A.add(entry.row(), entry.column(), entry.value());
    return A;
  }

  /** Element-wise addition of two sparse matrices. */
  SparseMatrix& operator+=(const SparseMatrix& other) const
  {
    Assert(rows == other.rows && cols == other.cols,
           "Dimension mismatch error.");
    for (const_iterator& entry : other)
      add(entry.row(), entry.column(), -entry.value());
    return *this;
  }

  /** Element-wise subtraction of two sparse matrices. */
  SparseMatrix operator-(const SparseMatrix& other) const
  {
    Assert(rows == other.rows && cols == other.cols,
           "Dimension mismatch error.");
    SparseMatrix A(*this);
    for (const_iterator& entry : other)
      A.add(entry.row(), entry.column(), entry.value());
    return A;
  }

  /** Element-wise subtraction of two sparse matrices. */
  SparseMatrix& operator-=(const SparseMatrix& other) const
  {
    Assert(rows == other.rows && cols == other.cols,
           "Dimension mismatch error.");
    for (const_iterator& entry : other)
      add(entry.row(), entry.column(), entry.value());
    return *this;
  }

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f[ \vec{y} = \boldsymbol{A} \vec{x} \\
   *     y_i = \sum_{j=1}^{n} a_{ij} x_j, \hspace{0.25cm} \forall i
   * \f]
   */
  Vector<number> operator*(const Vector<number>& x) const
  {
    Assert(cols == x.size(), "Dimension mismatch error.");
    Vector<number> y(rows);
    for (const_iterator& entry : *this)
      y[entry.row()] += entry.value() * x[entry.column()];
    return y;
  }

  /**
  * Compute a matrix-matrix product.
  * This is computed via
  * \f[ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}
  *     c_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj}, \hspace{0.25cm} \forall i, j
  * \f]
  */
  SparseMatrix operator*(const SparseMatrix& other) const
  {
    Assert(cols == other.n_rows(), "Dimension mismatch error.");

    SparseMatrix A(rows, other.n_cols());
    for (uint64_t i = 0; i < rows; ++i)
    {
      // Loop over relative column indices for row i of this
      for (const auto a_ij = this->begin(i); a_ij != this->end(i); ++a_ij)
      {
        uint64_t j = a_ij.column();
        number a_ij_val = a_ij.value();

        /* Loop over relative column indices of row j of other sparse matrix
         * and compute c_ik += a_ij b_jk. */
        for (const auto b_jk = other.begin(j); b_jk != other.end(j); ++b_jk)
          A.add(i, b_jk.column(), a_ij_val * b_jk.value());
      }
    }
  }

  /**
   * Return the transpose of a matrix.
   * This is computed via
   * \f[ \boldsymbol{B} = \boldsymbol{A}^T \\
   *     b_{ij} = a_{ji}, \hspace{0.25cm} \forall i, j
   * \f]
   */
  SparseMatrix transpose() const
  {
    SparseMatrix A(rows, cols, 0);
    for (const_iterator& entry : *this)
      A.set(entry.column(), entry.row(), entry.value());
    return A;
  }

  // @}

  //###########################################################################
  /** \name Printing */
  // @{

  /** Print the sparse matrix as triplets of nonzero entries. */
  void print(std::ostream& os) const
  {
    os.setf(std::ios::left, std::ios::adjustfield);

    os.fill(' '); os.precision(6);
    os << std::setw(8) << "Row"
       << std::setw(8) << "Column"
       << std::setw(10) << "Value" << std::endl;
    os << std::setw(8) << "---"
       << std::setw(8) << "------"
       << std::setw(10) << "-----" << std::endl;

    for (const auto entry : *this)
      os << std::setw(8)  << entry.row()
         << std::setw(8)  << entry.column()
         << std::setw(10) << entry.value() << std::endl;
  }

  /** Print the sparse matrix as a normal matrix with zero fill ins. */
  void print_formatted(std::ostream& os,
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

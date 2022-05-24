#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "vector.h"
#include "matrix.h"

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace pdes::Math
{


class SparseMatrix
{
public:
  using value_type = double;
  using SparsityPattern = std::vector<std::vector<size_t>>;

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
  SparseMatrix();

  /**
   * Construct a sparse matrix with \p n_rows and \p n_cols with
   * \p default_row_length nonzero entries per row.
   */
  explicit
  SparseMatrix(const size_t n_rows,
               const size_t n_cols,
               const size_t default_row_length);

  /**
   * Construct a square sparse matrix with dimension \p n with
   * \p default_row_length nonzero entries per row.
   */
  explicit
  SparseMatrix(const size_t n,
               const size_t default_row_length);

  /**
   * Construct a sparse matrix from a sparsity pattern. In this context, a
   * sparsity pattern is defined by a list of lists. Each inner list stores
   * the nonzero column indices for a given row.
   */
  SparseMatrix(SparsityPattern sparsity_pattern);


  /**
   * Copy construction from a dense matrix.
   */
  SparseMatrix(const Matrix& other);

  /**
   * Assignment with a dense matrix.
   */
  SparseMatrix&
  operator=(const Matrix& other);

  /**
   * Assignment with a scalar value. This routine sets all allocated data
   * to the specified value.
   */
  SparseMatrix&
  operator=(const value_type value);

  /**
   * Test the equality of two sparse matrices.
   */
  bool
  operator==(const SparseMatrix& other) const;

  /**
   * Test the inequality of two sparse matrices.
   * \param other
   * \return
   */
  bool
  operator!=(const SparseMatrix& other) const;

  /** \name Characteristics */
  // @{

  /**
   * Return the number of rows.
   */
  size_t
  n_rows() const;

  /**
   * Return the number of columns.
   */
  size_t
  n_cols() const;

  /**
   * Return the number of non-zero entries.
   */
  size_t
  nnz() const;

  /**
   * Return the length of row \p i.
   */
  size_t
  row_length(const size_t i) const;

  /**
   * Return whether the sparse matrix is empty.
   */
  bool
  empty() const;

  // @}
  /** \name Iterators */
  // @{

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
    CoeffIterator       coeff_ptr;

    void
    advance();

  public:
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
            value_type& val);
    };

  public:
    iterator(SparseMatrix* sparse_matrix, const size_t row);
    iterator(SparseMatrix* sparse_matrix);

    iterator&
    operator++();

    iterator
    operator++(int);

    entry
    operator*();

    bool
    operator==(const iterator& other) const;

    bool
    operator!=(const iterator& other) const;
  };


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
    row(SparseMatrix* sparse_matrix, const size_t i);

    iterator
    begin();

    iterator
    end();
  };


  /**
   * Return a mutable iterator to the first row of the sparse matrix.
   */
  iterator
  begin();

  /**
   * Return a mutable iterator to the end of the sparse matrix.
   */
  iterator
  end();

  /**
   * Return a mutable iterator to the start of row \p i.
   */
  iterator
  begin(const size_t i);

  /**
   * Return a mutable iterator to the end of row \p i.
   */
  iterator
  end(const size_t i);

  /**
   * A convenience function which allows for mutable range-based iteration
   * over the specified row \p i.
   * \see SparseMatrix::const_row_iterator
   */
  row
  row_iterator(const size_t i);


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
    size_t               current_row;
    ConstColumnIterator  col_ptr;
    ConstCoeffIterator   coeff_ptr;

    void
    advance()
    {
      // Increment along the current row
      ++col_ptr; ++coeff_ptr;

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
                  const value_type& val);
    };

  public:
    const_iterator(const SparseMatrix* sparse_matrix, const size_t row);
    const_iterator(const SparseMatrix* sparse_matrix);

    const_iterator&
    operator++();

    const_iterator
    operator++(int);

    const_entry
    operator*();

    bool
    operator==(const const_iterator& other) const;

    bool
    operator!=(const const_iterator& other) const;
  };


  /**
   * A struct defining a mutable row of the sparse matrix.
   * \see SparseMatrix::row
   */
  class const_row
  {
  private:
    const SparseMatrix* sparse_matrix_ptr;
    const size_t        row_num;

  public:
    const_row(const SparseMatrix* sparse_matrix, const size_t i);

    const_iterator
    begin();

    const_iterator
    end();
  };

  /**
   * Return a constant iterator to the start of the sparse matrix.
   */
  const_iterator
  begin() const;

  /**
   * Return a constant iterator to the end of the sparse matrix.
   */
  const_iterator
  end() const;

  /**
   * Return a constant iterator to the start of row \p i.
   */
  const_iterator
  begin(const size_t i) const;

  /**
   * Return a constant iterator to the end of row \p i.
   */
  const_iterator
  end(const size_t i) const;

  /**
   * A convenience function which allows for constant range-based iteration
   * over the specified row \p i.
   * \see SparseMatrix::row_iterator
   */
  const_row
  const_row_iterator(const size_t i) const;

  // @}
  /** \name Accessors */
  // @{

  /**
   * Return the column index located at relative position \p jr of row \p i.
   */
  const size_t&
  column(const size_t i, const size_t jr) const;


  /**
   * Read and write access to the element located at relative position \p jr
   * of row \p i.
   */
  value_type&
  value(const size_t i, const size_t jr);

  /**
   * Read access to the element located at relative position \p jr of row \p i.
   */
  const value_type&
  value(const size_t i, const size_t jr) const;

  /**
   * Return a pointer to element <tt>(i, j)</tt>. If no element exists,
   * \p nullptr is returned.
   */
  value_type*
  locate(const size_t i, const size_t j);

  /**
   * Return a constant pointer to element <tt>(i, j)</tt>. If no element
   * exists, \p nullptr is returned.
   */
  const value_type*
  locate(const size_t i, const size_t j) const;


  /**
   * Read and write access to element <tt>(i, j)</tt>.
   * \throw If column \p j does not exist on row \p i.
   */
  value_type&
  operator()(const size_t i, const size_t j);

  /**
   * Read access to element <tt>(i, j)</tt>.
   * \throw If column \p j does not exist on row \p i.
   */
  const value_type&
  operator()(const size_t i, const size_t j) const;

  /**
   * Return a pointer to the diagonal element of row \p i. If no diagonal
   * element exists, \p nullptr is returned.
   */
  value_type*
  diagonal(const size_t i);

  /**
   * Return a constant pointer to the diagonal element of row \p i. If no
   * diagonal element exists, \p nullptr is returned.
   */
  const value_type*
  diagonal(const size_t i) const;

  /**
   * Return a vector of pointers to the diagonal elements of the sparse matrix.
   * If any particular diagonal is not present, its value is \p nullptr.
   */
   std::vector<value_type*>
   diagonal();

  /**
   * Return a vector of constant pointers to the diagonal elements of the
   * sparse matrix. If any particular diagonal is not present, its value is
   * \p nullptr.
   */
  std::vector<const value_type*>
  diagonal() const;

  // @}
  /** \name Modifiers */
  // @{

  /**
   * Return the sparse matrix to an uninitialized state.
   */
  void
  clear();

  /**
   * Reinitialize the sparse matrix with \p n_rows and \p n_cols with
   * \p default_row_length nonzero entries per row.
   */
  void
  reinit(const size_t n_rows,
         const size_t n_cols,
         const size_t default_row_length);

  /**
   * Reinitialize the sparse matrix with dimension \p n with
   * \p default_row_length nonzero entries per row.
   */
  void
  reinit(const size_t n, const size_t default_row_length);

  /**
   * Set element <tt>(i, j)</tt> to \p value. If the element is initialized,
   * override the value. If it is not, initialize it.
   */
  void
  set(const size_t i, const size_t j, const value_type value);

  /**
   * Add \p value to element <tt>(i, j)</tt>. If the element is not initialized,
   * the <code> set(const size_t i, const size_t j) </code> is called.
   */
  void
  add(const size_t i, const size_t j, const value_type value);

  /**
   * Swap the elements of two rows.
   */
  void
  swap_row(const size_t i, const size_t k);

  /**
   * Swap the elements of this sparse matrix with another.
   */
  void
  swap(SparseMatrix& other);

  // @}
  /** \name Scalar Operations */
  // @{

  /**
   * Negate the elements of the sparse matrix. This is computed via
   * \f$ \boldsymbol{A} = -\boldsymbol{A} = -a_{ij}, ~ \forall i, j \f$.
   */
  SparseMatrix&
  operator-();

  /**
   * Return a sparse matrix with the negated elements.
   * \see SparseMatrix::operator-()
   */
  SparseMatrix
  operator-() const;

  /**
   * Multiply the elements of the sparse matrix by a scalar factor. This is
   * computed via
   * \f$ \boldsymbol{A} = \alpha \boldsymbol{A}
   *                    = \alpha a_{ij}, ~ \forall i, j
   * \f$.
   */
  SparseMatrix&
  operator*=(const value_type factor);

  /**
   * Divide the elements of the sparse matrix by a scalar factor. This is
   * computed via
   * \f$ \boldsymbol{A} = \frac{1}{\alpha} \boldsymbol{A}
   *                    = \frac{\alpha}{a_{ij}}, ~ \forall i, j
   * \f$.
   */
  SparseMatrix&
  operator/=(const value_type factor);

  // @}
  /** \name Linear Algebra */
  // @{

  /**
   * Add a sparse matrix multiplied by a scalar. This is only allowed for
   * sparse matrices with the same sparsity pattern and computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A} + \alpha \boldsymbol{B}
   *                    = \sum_{i=0}^{n} \sum_{j=0}^{n} a_{ij} + \alpha b_{ij},
   *                    ~ \forall i,j
   * \f$.
   */
  void
  add(const SparseMatrix& B,
      const value_type factor = 1.0);

  /**
   * Add another sparse matrix.
   * \see SparseMatrix::add(const SparseMatrix& const value_type)
   */
  SparseMatrix&
  operator+=(const SparseMatrix& B);

  /**
   * Subtract another matrix. his is only allowed for sparse matrices
   * with the same sparsity pattern and is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} - \boldsymbol{B}
   *                    = a_{ij} - b_{ij}, ~ \forall i, j
   * \f$.
   */
  SparseMatrix&
  operator-=(const SparseMatrix& B);

  /**
   * Compute a matrix-vector product. This is computed via
   * \f$ \vec{y} = \boldsymbol{A} \vec{x}
   *             = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i
   * \f$.
   */
  void
  vmult(const Vector& x, Vector& y,
        const bool adding = false) const;

  /**
   * Return a matrix-vector product.
   * \see SparseMatrix::vmult
   */
  Vector
  vmult(const Vector& x) const;

  /**
   * Compute a matrix-vector product and add to the destination vector.
   * This is computed via
   * \f$ \vec{y} = \vec{y} + \boldsymbol{A} \vec{x}
   *             = y_i + \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i
   * \f$.
   */
  void
  vmult_add(const Vector& x, Vector& y) const;

  /**
   * Compute a transpose matrix-vector product. This is computed via
   * \f$ \vec{y} = \boldsymbol{A}^T \vec{x}
   *             = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall i
   * \f$.
   */
  void Tvmult(const Vector& x, Vector& y,
              const bool adding = false) const;

  /**
   * Return a transpose matrix-vector product.
   * \see SparseMatrix::Tvmult
   */
  Vector
  Tvmult(const Vector& x) const;

  /**
   * Compute a transpose matrix-vector product and to the destination vector.
   * This is computed via
   * \f[ \vec{y} = \vec{y} + \boldsymbol{A}^T \vec{x} \\
   *             = y_i + \sum_{i=1}^{n} a_{ji} x_i, ~ \forall i
   * \f]
   */
  void Tvmult_add(const Vector& x, Vector& y) const;

  /**
   * Compute a matrix-vector product.
   * \see SparseMatrix::vmult
   */
  Vector
  operator*(const Vector& x) const;

  // @}
  /** \name Printing Utilities */
  // @{

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
        const unsigned int precision = 3) const;

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
                  const unsigned int width = 0) const;

  /**
   * Return the sparse matrix as a string.
   * \param scientific A flag for scientific notation.
   * \param precision The precision to display to sparse matrix elements.
   * \param width The spacing between entries.
   */
  std::string
  str(const bool scientific = false,
      const unsigned int precision = 3,
      const unsigned int width = 0) const;

  // @}
};

/*-------------------- Inline Implementations --------------------*/

/**
 * Multiply a sparse matrix by a scalar factor.
 * \see SparseMatrix::operator*=
 */
SparseMatrix
operator*(const double factor, const SparseMatrix& A);

/**
 * Multiply a sparse matrix by a scalar factor.
 * \see SparseMatrix::operator*=
 */
SparseMatrix
operator*(const SparseMatrix& A, const double factor);

/**
 * Divide a sparse matrix by a scalar factor.
 * \see SparseMatrix::operator/=
 */
SparseMatrix
operator/(const SparseMatrix& A, const double factor);

/**
 * Add two sparse matrices together.
 * \see SparseMatrix::operator+=
 *      SparseMatrix::add(const SparseMatrix&, const double)
 */
SparseMatrix
operator+(const SparseMatrix& A, const SparseMatrix& B);

/**
 * Subtract two sparse matrices.
 * \see SparseMatrix::operator-=
 *      SparseMatrix::add(const SparseMatrix&, const double)
 */
SparseMatrix
operator-(const SparseMatrix& A, const SparseMatrix& B);

/**
 * Compute a matrix-vector product.
 * \see SparseMatrix::vmult
 */
void
vmult(const SparseMatrix& A, const Vector& x, Vector& y);

/**
 * Return a matrix-vector product.
 * \see SparseMatrix::vmult
 */
Vector
vmult(const SparseMatrix& A, const Vector& x);

/**
 * Compute a transpose matrix-vector product.
 * \see SparseMatrix::vmult
 */
void
Tvmult(const SparseMatrix& A, const Vector& x, Vector& y);

/**
 * Return a transpose matrix-vector product.
 * \see SparseMatrix::vmult
 */
Vector
Tvmult(const SparseMatrix& A, const Vector& x);

/**
 * Insert a sparse matrix into an output stream
 */
std::ostream&
operator<<(std::ostream& os, const SparseMatrix& A);

}
#endif //SPARSE_MATRIX_H

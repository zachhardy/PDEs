#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "utilities.h"

#include <iostream>
#include <sstream>
#include <cstddef>
#include <vector>


namespace Math
{
  //########## Forward declarations
  class Vector;
  class Matrix;
  class SparseMatrix;


  namespace SparseMatrixIterators
  {
    class Iterator
    {
    private:
      class Accessor
      {
      public:

        Accessor(SparseMatrix* matrix,
                 const size_t row,
                 const size_t index);

        Accessor(SparseMatrix* matrix);

        size_t row() const;
        size_t column() const;
        size_t index() const;
        double& value() const;

        bool operator==(const Accessor& other) const;
        bool operator!=(const Accessor& other) const;
        bool operator<(const Accessor& other) const;

      protected:
        void advance();


        SparseMatrix* matrix;
        size_t current_row;
        size_t current_index;

        friend class Iterator;
      };

    public:
      Iterator(SparseMatrix* matrix,
               const size_t row,
               const size_t index);

      Iterator(SparseMatrix* matrix);

      Iterator& operator++();
      Iterator operator++(int);

      const Accessor& operator*() const;
      const Accessor* operator->() const;

      bool operator==(const Iterator& other) const;
      bool operator!=(const Iterator& other) const;
      bool operator<(const Iterator& other) const;

    private:
      Accessor accessor;
    };


    class ConstIterator
    {
    private:
      class Accessor
      {
      public:
        Accessor(const SparseMatrix* matrix,
                 const size_t row,
                 const size_t index);

        Accessor(const SparseMatrix* matrix);

        size_t row() const;
        size_t column() const;
        size_t index() const;
        const double& value() const;

        bool operator==(const Accessor& other) const;
        bool operator!=(const Accessor& other) const;
        bool operator<(const Accessor& other) const;

      protected:
        void advance();


        const SparseMatrix* matrix;
        size_t current_row;
        size_t current_index;

        friend class ConstIterator;
      };

    public:
      ConstIterator(const SparseMatrix* matrix,
                    const size_t row,
                    const size_t index);

      ConstIterator(const SparseMatrix* matrix);

      ConstIterator& operator++();
      ConstIterator operator++(int);

      const Accessor& operator*() const;
      const Accessor* operator->() const;

      bool operator==(const ConstIterator& other) const;
      bool operator!=(const ConstIterator& other) const;
      bool operator<(const ConstIterator& other) const;

    private:
      Accessor accessor;
    };


    class RowIterator
    {
    private:
      SparseMatrix* matrix;
      const size_t row;

    public:
      RowIterator(SparseMatrix* matrix,
                  const size_t row);

      Iterator begin();
      Iterator end();
    };


    class ConstRowIterator
    {
    private:
      const SparseMatrix* matrix;
      const size_t row;

    public:
      ConstRowIterator(const SparseMatrix* matrix,
                       const size_t row);

      ConstIterator begin() const;
      ConstIterator end() const;
    };

  }


  /** Implementation of a list of lists sparse matrix. */
  class SparseMatrix
  {
  public:
    using value_type = double;

    using iterator = SparseMatrixIterators::Iterator;
    using const_iterator = SparseMatrixIterators::ConstIterator;

    using RowIterator = SparseMatrixIterators::RowIterator;
    using ConstRowIterator = SparseMatrixIterators::ConstRowIterator;

    static const size_t invalid_entry = numbers::invalid_size_t;

  private:
    bool has_entries;

    size_t rows;
    size_t cols;

    std::vector<std::vector<size_t>> colnums;
    std::vector<std::vector<value_type>> vals;

  public:
    //################################################## Constructors

    /** \name Constructors and Initializers */
    // @{

    /** Construct an empty sparse matrix. */
    SparseMatrix();

    /** Construct a square sparse matrix with dimension \p n. */
    SparseMatrix(const size_t n);

    /** Construct a sparse matrix with \p n_rows and \p n_cols. */
    SparseMatrix(const size_t n_rows, const size_t n_cols);

    /** Assign all entries in the sparse matrix to a scalar value. */
    SparseMatrix& operator=(const value_type value);

    // @}

    //################################################## Information

    /** \name Information */
    // @{

    size_t n_rows() const;
    size_t n_cols() const;
    size_t nnz() const;

    size_t row_length(const size_t row) const;

    size_t column(const size_t row, const size_t index) const;

    /**
     * Return the relative index within the \p row of the specified \p col.
     * If the column does not exist, then the invalid entry value is returned.
     */
    size_t index(const size_t row, const size_t col) const;

    bool empty() const;
    bool exists(const size_t row, const size_t col) const;

    bool operator==(const SparseMatrix& other) const;
    bool operator!=(const SparseMatrix& other) const;

    // @}

    //################################################## Iterators

    /** \name Iterators */
    // @{

    iterator begin();
    iterator end();

    iterator begin(const size_t row);
    iterator end(const size_t row);

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator begin(const size_t row) const;
    const_iterator end(const size_t row) const;

    RowIterator row_iterator(const size_t row);
    ConstRowIterator row_iterator(const size_t row) const;


    // @}

    //################################################## Accessors

    /** \name Accessors */
    // @{

    /**
     * Return a reference to the value of the sparse matrix at row \p i and
     * column \p j. If the entry does not exist, an error is thrown. When
     * potential uninitialized elements may be encountered, the \ref el method should
     * be used.
     */
    value_type& operator()(const size_t i, const size_t j);

    /**
     * Return a const reference to the value of the sparse matrix at row \p i
     * and column \p j. If the entry does not exist, an error is thrown. When
     * potential uninitialized elements may be encountered, the \ref el method
     * should be used.
     */
    const value_type& operator()(const size_t i, const size_t j) const;

    /**
     * Return the value of the sparse matrix at row \p i and column \p j. If
     * the entry is uninitialized, zero is returned.
     */
    value_type el(const size_t i, const size_t j);

    /**
     * Return a reference to the <tt>i</tt>'th diagonal entry of the sparse
     * matrix. If the entry does not exist, an error is thrown. The \ref
     * diag_el method should be used if the diagonal may be zero.
     */
    value_type& diag(const size_t i);

    /**
     * Return a const reference to the <tt>i</tt>'th diagonal entry of the
     * sparse matrix. If the entry does not exist, an error is thrown. The \ref
     * diag_el method should be used if the diagonal may be zero.
     */
    const value_type& diag(const size_t i) const;

    /**
     * Return the value of the <tt>i</tt>'th diagonal entry of the sparse
     * matrix. If the entry does not exist, zero is returned.
     */
    value_type diag_el(const size_t i) const;

    // @}

    //================================================== Modifiers

    /** \name Modifiers */
    // @{

    /** Set the sparse matrix to an uninitialized state. */
    void clear();

    /** Reinitialize the sparse matrix with \p n_rows and \p n_cols. */
    void reinit(const size_t n_rows, const size_t n_cols);

    /** Copy the non-zero contents of a dense matrix. */
    void copy_from(const Matrix& matrix);

    /**
     * Set element <tt>(row, col)</tt> to \p value. If the element is
     * initialized, override the value. If it is not, initialize it.
     */
    void set(const size_t row,
             const size_t col,
             const value_type value);

    /**
     * Add \p value to element <tt>(row, col)</tt>. If the element is not
     * initialized, then set it, otherwise perform an add operation.
     */
    void add(const size_t row,
             const size_t col,
             const value_type value);


    void swap_row(const size_t i, const size_t k);
    void swap(SparseMatrix& other);

    // @}

    //================================================== Linear Algebra

    /** \name Linear Algebra */
    // @{

    /** Multiply the entries of the matrix by a scalar \p factor. */
    SparseMatrix& scale(const value_type factor);

    /**
     * Add another sparse matrix scaled by a scalar \p factor.
     *
     * \note To avoid expensive modifications to underlying structrue, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& add(const SparseMatrix& B,
                      const value_type factor = 1.0);

    /**
     * Scale this sparse matrix by the scalar value \p a and add another
     * sparse matrix.
     *
     * \note To avoid expensive modifications to underlying structrue, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& sadd(const value_type a, const SparseMatrix& B);

    /**
     * Scale this sparse matrix by the scalar value \p a and add another
     * sparse matrix scaled by the scalar value \p b.
     *
     * \note To avoid expensive modifications to underlying structrue, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& sadd(const value_type a,
                       const value_type b,
                       const SparseMatrix& B);

    /**
     * Compute a matrix-vector product. If the \p adding flag is set to \p true,
     * the matrix-vector product is added to the existing data within the
     * destination vector \f$ \vec{y} \f$.
     */
    void vmult(const Vector& x,
               Vector& y,
               const bool adding = false) const;

    /**
     * Compute a transpose matrix-vector product. If the \p adding flag is set
     * to \p true, the matrix-vector product is added to the existing data
     * within the destination vector \f$ \vec{y} \f$.
     */
    void Tvmult(const Vector& x,
                Vector& y,
                const bool adding = false) const;


    /** Return the result of a matrix-vector product. */
    Vector vmult(const Vector& x) const;

    /** Return the result of a transpose matrix-vector product. */
    Vector Tvmult(const Vector& x) const;

    /**
     * Add the result of a matrix-vector product to the destination vector
     * \f$ \vec{y} \f$ such that \f$ \vec{y} += \boldsymbol{A} \vec{x} \f$.
     */
    void vmult_add(const Vector& x, Vector& y) const;

    /**
     * Add the result of a transpose matrix-vector product to the destination
     * vector \f$ \vec{y} \f$ such that \f$ \vec{y} += \boldsymbol{A}^T \vec{x}
     * \f$.
     */
    void Tvmult_add(const Vector& x, Vector& y) const;

    /** Negate the entries of the sparse matrix. */
    SparseMatrix& operator-();

    /** Multiply the entries of the sparse matrix by a scalar \p factor. */
    SparseMatrix& operator*=(const value_type factor);

    /** Divide the entries of the sparse matrix by a scalar \p factor. */
    SparseMatrix& operator/=(const value_type factor);

    /**
     * Add another sparse matrix.
     *
     * \note To avoid expensive modifications to underlying structrue, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& operator+=(const SparseMatrix& B);

    /**
     * Subtract another sparse matrix.
     *
     * \note To avoid expensive modifications to underlying structrue, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& operator-=(const SparseMatrix& B);

    /** Compute a matrix-vector product. */
    Vector operator*(const Vector& x) const;

    // @}

    //================================================== Print Utilities

    /** \name Print Utilities */
    // @{

    /**
     *  Print the sparse matrix as triplets of nonzero entries.
     *
     * \param os The output stream to print the matrix to.
     * \param scientific A flag for scientific notation.
     * \param precision The precision to display to sparse matrix elements.
     */
    void print(std::ostream& os = std::cout,
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
    void print_formatted(std::ostream& os = std::cout,
                         const bool scientific = false,
                         const unsigned int precision = 3,
                         const unsigned int width = 0) const;

    /**
     * Return the sparse matrix as a string.
     *
     * \param scientific A flag for scientific notation.
     * \param precision The precision to display to sparse matrix elements.
     * \param width The spacing between entries.
     */
    std::string str(const bool scientific = false,
                    const unsigned int precision = 3,
                    const unsigned int width = 0) const;

    friend std::ostream&
    operator<<(std::ostream& os, const SparseMatrix& A);

    // @}

    friend class SparseMatrixIterators::Iterator;
    friend class SparseMatrixIterators::RowIterator;

    friend class SparseMatrixIterators::ConstIterator;
    friend class SparseMatrixIterators::ConstRowIterator;
  };


  /**
   * Insert a sparse matrix into an output stream
   *
   * \see SparseMatrix::str SparseMatrix::print
   */
  std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);

}
#endif //SPARSE_MATRIX_H

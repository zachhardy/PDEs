#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

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
    /** An iterator over all of the elements in a SparseMatrix. */
    class Iterator
    {
    private:
      /** A utility for accessing entries of the SparseMatrix. */
      class Accessor
      {
      protected:
        SparseMatrix* matrix;
        size_t current_row;
        size_t current_index;

      public:
        /** Construct for a particular entry. */
        Accessor(SparseMatrix* matrix,
                 const size_t row,
                 const size_t index);

        /** Construct an invalid accessor. */
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

        friend class Iterator;
      };

    private:
      Accessor accessor;

    public:
      /** Construct an iterator to a particular entry. */
      Iterator(SparseMatrix* matrix,
               const size_t row,
               const size_t index);

      /** Construct the end iterator. */
      Iterator(SparseMatrix* matrix);

      Iterator& operator++();
      Iterator operator++(int);

      const Accessor& operator*() const;
      const Accessor* operator->() const;

      bool operator==(const Iterator& other) const;
      bool operator!=(const Iterator& other) const;
      bool operator<(const Iterator& other) const;
    };

    /** A constant iterator over all the elements in a SparseMatrix. */
    class ConstIterator
    {
    private:
      /** A utility for accessing entries of the SparseMatrix. */
      class Accessor
      {
      protected:
        const SparseMatrix* matrix;
        size_t current_row;
        size_t current_index;

      public:
        /** Construct for a particular entry. */
        Accessor(const SparseMatrix* matrix,
                 const size_t row,
                 const size_t index);

        /** Construct an invalid accessor. */
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

        friend class ConstIterator;
      };

    private:
      Accessor accessor;

    public:
      /** Construct a constant iterator to a particular entry. */
      ConstIterator(const SparseMatrix* matrix,
                    const size_t row,
                    const size_t index);

      /** Construct the end constant iterator. */
      ConstIterator(const SparseMatrix* matrix);

      ConstIterator& operator++();
      ConstIterator operator++(int);

      const Accessor& operator*() const;
      const Accessor* operator->() const;

      bool operator==(const ConstIterator& other) const;
      bool operator!=(const ConstIterator& other) const;
      bool operator<(const ConstIterator& other) const;
    };


    /** An iterator over a single row of a SparseMatrix. */
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

    /** A constant iterator over a single row of a SparseMatrix. */
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
    using iterator = SparseMatrixIterators::Iterator;
    using const_iterator = SparseMatrixIterators::ConstIterator;

    using RowIterator = SparseMatrixIterators::RowIterator;
    using ConstRowIterator = SparseMatrixIterators::ConstRowIterator;

    static const size_t invalid_entry = -1;

  private:
    bool has_entries;

    size_t rows;
    size_t cols;

    std::vector<std::vector<size_t>> colnums;
    std::vector<std::vector<double>> vals;

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

    /** Element-wise assignment to a scalar. */
    SparseMatrix& operator=(const double value);

    // @}

    //################################################## Information

    /** \name Information */
    // @{

    size_t n_rows() const;
    size_t n_cols() const;

    /** Return the number of nonzero elements. */
    size_t nnz() const;

    size_t row_length(const size_t row) const;

    size_t column(const size_t row, const size_t index) const;

    /**
     * Return the nonzero index on \p row that corresponds to \p col. If the
     * column does not exist, an invalid value is returned.
     */
    size_t index(const size_t row, const size_t col) const;

    bool empty() const;

    /** Return whether the element (\p row, \p col) has been allocated. */
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
     * potential uninitialized elements may be encountered, the \ref el method
     * should be used.
     */
    double& operator()(const size_t i, const size_t j);

    /**
     * Return a const reference to the value of the sparse matrix at row \p i
     * and column \p j. If the entry does not exist, an error is thrown. When
     * potential uninitialized elements may be encountered, the \ref el method
     * should be used.
     */
    const double& operator()(const size_t i, const size_t j) const;

    /**
     * Return the value of the sparse matrix at row \p i and column \p j. If
     * the entry is uninitialized, zero is returned.
     */
    double el(const size_t i, const size_t j);

    /**
     * Return a reference to the <tt>i</tt>'th diagonal entry of the sparse
     * matrix. If the entry does not exist, an error is thrown. The \ref
     * diag_el method should be used if the diagonal may be zero.
     */
    double& diag(const size_t i);

    /**
     * Return a const reference to the <tt>i</tt>'th diagonal entry of the
     * sparse matrix. If the entry does not exist, an error is thrown. The \ref
     * diag_el method should be used if the diagonal may be zero.
     */
    const double& diag(const size_t i) const;

    /**
     * Return the value of the <tt>i</tt>'th diagonal entry of the sparse
     * matrix. If the entry does not exist, zero is returned.
     */
    double diag_el(const size_t i) const;

    // @}

    //================================================== Modifiers

    /** \name Modifiers */
    // @{

    void clear();

    /**
     * Reinitialize the sparse matrix with \p n_rows and \p n_cols. This
     * clears the existing data.
     */
    void reinit(const size_t n_rows, const size_t n_cols);

    /** Copy the non-zero contents of a dense matrix. */
    void copy_from(const Matrix& matrix);

    /** Copy from another sparse matrix. */
    void copy_from(const SparseMatrix& matrix);

    /**
     * Set element (\p row, \p col) to \p value. If the element is initialized,
     * override the value. If it is not, initialize it.
     */
    void set(const size_t row,
             const size_t col,
             const double value);

    /**
     * Add \p value to element (\p row, \p col). If the element is not
     * initialized, then set it, otherwise perform an add operation.
     */
    void add(const size_t row,
             const size_t col,
             const double value);

    void swap_row(const size_t i, const size_t k);
    void swap(SparseMatrix& other);

    // @}

    //================================================== Linear Algebra

    /** \name Linear Algebra */
    // @{

    /** Element-wise multiplication by a scalar in place. */
    SparseMatrix& scale(const double factor);

    /**
     * Element-wise addition with a SparseMatrix scaled by a scalar \p factor
     * in place.
     *
     * \note To avoid expensive modifications to underlying structure, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& add(const SparseMatrix& B,
                      const double factor = 1.0);

    /**
     * Element-wise multiplication by a scalar \p a followed by element-wise
     * addition with a SparseMatrix in place.
     *
     * \note To avoid expensive modifications to underlying structure, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& sadd(const double a, const SparseMatrix& B);

    /**
     * Element-wise multiplication by a scalar \p a followed by element-wise
     * addition with a SparseMatrix scaled by \p b in place.
     *
     * \note To avoid expensive modifications to underlying structure, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& sadd(const double a, const double b, const SparseMatrix& B);

    /**
     * Compute a matrix-vector product.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     * \param adding A flag for adding to or setting the destination Vector.
     */
    void vmult(const Vector& x, Vector& y,
               const bool adding = false) const;

    /**
     * Compute a transpose matrix-vector product.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     * \param adding A flag for adding to or setting the destination Vector.
     */
    void Tvmult(const Vector& x, Vector& y,
                const bool adding = false) const;


    /** Return a matrix-vector product. */
    Vector vmult(const Vector& x) const;

    /** Return a transpose matrix-vector product. */
    Vector Tvmult(const Vector& x) const;

    /**
     * Add a matrix-vector product to the destination Vector.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     */
    void vmult_add(const Vector& x, Vector& y) const;

    /**
     * Add a transpose matrix-vector product to the destination Vector.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     */
    void Tvmult_add(const Vector& x, Vector& y) const;

    /** Element-wise negation in place. */
    SparseMatrix& operator-();

    /** Element-wise multiplication by a scalar in place. */
    SparseMatrix& operator*=(const double factor);

    /** Element-wise division by a scalar in place. */
    SparseMatrix& operator/=(const double factor);

    /**
     * Element-wise addition with a SparseMatrix in place.
     *
     * \note To avoid expensive modifications to underlying structure, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& operator+=(const SparseMatrix& B);

    /**
     * Element-wise subtraction with a SparseMatrix in place.
     *
     * \note To avoid expensive modifications to underlying structure, the
     *       matrices must have the same sparsity pattern.
     */
    SparseMatrix& operator-=(const SparseMatrix& B);

    /** Return a matrix-vector product. */
    Vector operator*(const Vector& x) const;

    // @}

    //================================================== Print Utilities

    /** \name Print Utilities */
    // @{

    /**
     *  Print the sparse matrix as triplets of nonzero entries.
     *
     * \param scientific A flag for scientific notation.
     * \param precision The precision to display to sparse matrix elements.
     * \param os The output stream to print the matrix to.
     */
    void print(const bool scientific = true,
               const unsigned int precision = 3,
               std::ostream& os = std::cout) const;

    /**
     * Print a row of the sparse matrix.
     *
     * \param row The row to print.
     * \param scientific A flag for scientific notation.
     * \param precision The precision to display to sparse matrix elements.
     * \param os The output stream to print the matrix to.
     */
     void print_row(const size_t row,
                    const bool scientific = true,
                    const unsigned int precision = 3,
                    std::ostream& os = std::cout) const;

    /**
     * Print the matrix as if it were dense, filling in zero values.
     *
     * \param scientific A flag for scientific notation.
     * \param precision The precision to display to sparse matrix elements.
     * \param width The spacing between entries.
     * \param os The output stream to print the matrix to.
     */
    void print_formatted(const bool scientific = true,
                         const unsigned int precision = 3,
                         const unsigned int width = 0,
                         std::ostream& os = std::cout) const;

    /**
     * Return the sparse matrix as a string.
     *
     * \param scientific A flag for scientific notation.
     * \param precision The precision to display to sparse matrix elements.
     * \param width The spacing between entries.
     */
    std::string str(const bool scientific = true,
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


  std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);



}
#endif //SPARSE_MATRIX_H

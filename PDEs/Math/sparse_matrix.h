#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <iostream>
#include <sstream>
#include <cstddef>
#include <vector>


namespace PDEs
{
  namespace Math
  {
    //forward declarations
    class Vector;
    class Matrix;


    /**
     * Implementation of a list of lists sparse matrix.
     */
    class SparseMatrix
    {
    public:
      /**
       * Iterator implementation for the elements of the sparse matrix.
       */
      class Iterator
      {
      private:
        /**
         * A struct that describes the sparse matrix entry the iterator is
         * pointing to and provides a mutable reference to the underlying value.
         */
        struct Entry
        {
          const size_t row;
          const size_t column;
          double& value;

          /** Default constructor for a sparse matrix entry. */
          Entry(const size_t row,
                const size_t column,
                double& value);
        };

      private:
        SparseMatrix* matrix;
        size_t row;
        unsigned int index;

      public:
        /**
         * Construct an iterator that references a particular row and non-zero
         * index on that row of a sparse matrix.
         */
        Iterator(SparseMatrix* matrix,
                 const size_t row,
                 const unsigned int index);

        /**
         * Construct an iterator to designate the end of the sparse matrix
         * entries. This sets the row and index to an invalid element.
         */
        Iterator(SparseMatrix* matrix);

        /** Pre-increment operator. See \ref advance(). */
        Iterator& operator++();

        /** Post-increment operator. See \ref advance(). */
        Iterator operator++(int);

        /** Dereference operator. Return an Entry object. */
        Entry operator*() const;

        /** Return whether two iterators reference the same entry. */
        bool operator==(const Iterator& other) const;

        /** Return whether two iterators reference different entries. */
        bool operator!=(const Iterator& other) const;

        /**
         * Return whether this iterator references an element on a previous
         * index of the current row or a previous row.
         */
        bool operator<(const Iterator& other) const;

      private:
        /**
         * A utility function to increment the iterator.
         *
         * If the next entry is on the end of a row, the iterator is moved to
         * the next row. If the iterator is currently on the last entry of the
         * last row, the invalid iterator is returned.
         */
        void advance();
      };


      /**
       * An implementation of a constant iterator over the elements of the
       * sparse matrix.
       */
      class ConstIterator
      {
      private:
        //################################################## ConstEntry

        /**
         * A struct that describes the sparse matrix entry the constant
         * iterator is pointing to and provides a constant reference to the
         * underlying value.
         */
        struct ConstEntry
        {
          const size_t row;
          const size_t column;
          const double& value;

          /** Default constructor for a constant sparse matrix entry. */
          ConstEntry(const size_t row,
                     const size_t column,
                     const double& value);
        };

      private:
        const SparseMatrix* matrix;
        size_t row;
        unsigned int index;

      public:
        /**
         * Construct a constant iterator that references a particular row and
         * non-zero index on that row of a sparse matrix.
         */
        ConstIterator(const SparseMatrix* matrix,
                      const size_t row,
                      const unsigned int index);

        /**
         * Construct a constant iterator to designate the end of the sparse
         * matrix elements. This sets the row and index to an invalid element.
         */
        ConstIterator(const SparseMatrix* matrix);

        /** Pre-increment operator. See \ref advance(). */
        ConstIterator& operator++();

        /** Post-increment operator. See \ref advance(). */
        ConstIterator operator++(int);

        /** Dereference operator. Return a ConstEntry object. */
        ConstEntry operator*() const;

        /** Return whether two constant iterators reference the same entry. */
        bool operator==(const ConstIterator& other) const;

        /**
         * Return whether two constant iterators reference different entries.
         */
        bool operator!=(const ConstIterator& other) const;

        /**
         * Return whether this constant iterator references an entries on a
         * previous index of the current row or a previous row.
         */
        bool operator<(const ConstIterator& other) const;

      private:
        /**
         * A utility function to increment the constant iterator. If the next
         * entry is on the end of a row, the constant iterator is moved to the
         * next row. If the constant iterator is currently on the last entry
         * of the last row, the invalid constant iterator is returned.
         */
        void advance();
      };


      /**
       * An implementation of an iterator over a row of the sparse matrix.
       * This is just a specialization of the Iterator which implements
       * a \p begin and \p end routine.
       */
      class RowIterator
      {
      private:
        SparseMatrix* matix;
        const size_t row;

      public:
        /** Construct a row iterator pointing to the specified row. */
        RowIterator(SparseMatrix* matrix, const size_t row);

        /** Return an iterator to the first entry of the \p row. */
        Iterator begin();

        /** Return the entry to the end of the \p row. */
        Iterator end();
      };


      /**
       * An implementation of a constant iterator over a row of the sparse
       * matrix. This is just a specialization of the ConstIterator which
       * implements a \p begin and \p end routine.
       */
      class ConstRowIterator
      {
      private:
        const SparseMatrix* matix;
        const size_t row;

      public:
        /** Construct a constant row iterator pointing to the specified row. */
        ConstRowIterator(const SparseMatrix* matrix, const size_t row);

        /** Return a constant iterator to the first entry of the \p row. */
        ConstIterator begin();

        /** Return the constant iterator to the end of the \p row. */
        ConstIterator end();
      };

    private:
      bool has_entries;

      size_t rows;
      size_t cols;

      /**
       * Row-wise storage of the column indices which have non-zero entries on
       * the particular row.
       */
      std::vector<std::vector<size_t>> colnums;

      /**
       * Row-wise storage of the values of the non-zero entries of the
       * particular row.
       */
      std::vector<std::vector<double>> values;

    public:
      //################################################## Constructors

      /** Default constructor. Construct an empty sparse matrix. */
      SparseMatrix();

      /** Construct a sparse matrix with \p n_rows and \p n_cols. */
      SparseMatrix(const size_t n_rows, const size_t n_cols);

      /** Entry-wise assignment of the sparse matrix entries to a scalar. */
      SparseMatrix& operator=(const double value);

      /**
       * Reinitialize the sparse matrix with \p n_rows and \p n_cols.
       *
       * This clears the existing data and reallocates memory.
       */
      void reinit(const size_t n_rows, const size_t n_cols);

      /** Copy the non-zero contents of a dense matrix. */
      void copy_from(const Matrix& matrix);

      /** Copy from another sparse matrix. */
      void copy_from(const SparseMatrix& matrix);

      //################################################## Capacity

      /** Return the number of rows. */
      size_t n_rows() const;

      /** Return the number of columns. */
      size_t n_cols() const;

      /** Return the number of non-zero entries. */
      size_t
      n_nonzero_entries() const;

      /** Return the length (number of non-zero entries) of a \p row. */
      unsigned int row_length(const size_t row) const;

      /** Return whether the sparse matrix is empty or not. */
      bool empty() const;

      /**
       * Return whether an entry for a \p row and \p column has been allocated.
       */
      bool exists(const size_t row, const size_t column) const;

      //################################################## Data Access

      /**
       * Return a reference to the value of the sparse matrix at row \p i and
       * column \p j.
       *
       * If the entry does not exist, an error is thrown. When
       * potential uninitialized entries may be encountered, the \ref el method
       * should be used.
       */
      double& operator()(const size_t i, const size_t j);

      /**
       * Return a constant reference to the value of the sparse matrix at row
       * \p i and column \p j. If the entry does not exist, an error is thrown.
       * When potential uninitialized entries may be encountered, the \ref el
       * method should be used.
       */
      const double& operator()(const size_t i, const size_t j) const;

      /**
       * Return the value of the sparse matrix at row \p i and column \p j.
       *
       * If the entry is uninitialized, zero is returned.
       */
      double el(const size_t i, const size_t j);

      /**
       * Return a reference to the <tt>i</tt>'th diagonal entry of the sparse
       * matrix.
       *
       * If the entry does not exist, an error is thrown. The \ref diag_el
       * method should be used if the diagonal may be zero.
       */
      double& diag(const size_t i);

      /**
       * Return a constant reference to the <tt>i</tt>'th diagonal entry of
       * the sparse matrix.
       *
       * If the entry does not exist, an error is thrown. The \ref diag_el
       * method should be used if the diagonal may be zero.
       */
      const double& diag(const size_t i) const;

      /**
       * Return the value of the <tt>i</tt>'th diagonal entry of the sparse
       * matrix.
       *
       * If the entry does not exist, zero is returned.
       */
      double diag_el(const size_t i) const;

      /**
       * Return the column index for a particular non-zero \p index of a \p row.
       */
      size_t column(const size_t row, const unsigned int index) const;

      /**
       * Return the non-zero index on a \p row that corresponds to a \p column.
       *
       * If the \p column does not exist, an invalid value is returned.
       */
      unsigned int index(const size_t row, const size_t column) const;

      /** Return an iterator to the first entry of the first row. */
      Iterator begin();

      /** Return an iterator that designates the end of the sparse matrix. */
      Iterator end();

      /** Return an iterator to the first entry of \p row. */
      Iterator begin(const size_t row);

      /** Return an iterator to the end of \p row. */
      Iterator end(const size_t row);

      /** Return a constant iterator to the first entry of the first row. */
      ConstIterator begin() const;

      /**
       * Return a constant iterator that designates the end of the sparse
       * matrix.
       */
      ConstIterator end() const;

      /** Return a constant iterator to the first entry of \p row. */
      ConstIterator begin(const size_t row) const;

      /** Return a constant iterator to the end of \p row. */
      ConstIterator end(const size_t row) const;

      /** Return a row iterator for the specified \p row. */
      RowIterator row_iterator(const size_t row);

      /** Return a constant row iterator for \p row. */
      ConstRowIterator row_iterator(const size_t row) const;

      //################################################## Modifiers

      /** Delete the contents of the sparse matrix. */
      void clear();

      /**
       * Set the entry for the specified \p row and \p column to \p value.
       *
       * If the entry is initialized, override the value. If it is not,
       * initialize  it.
       */
      void set(const size_t row, const size_t column, const double value);

      /**
       * Add \p value to the entry at the specified \p row and \p column.
       *
       * If the entry is not initialized, then set it, otherwise add to it.
       */
      void add(const size_t row, const size_t column, const double value);

      /** Swap the contents of two rows of the sparse matrix. */
      void swap_row(const size_t i, const size_t k);

      /** Swap the contents of two sparse matrices. */
      void swap(SparseMatrix& other);

      //################################################## Scaling Operations

      /** Element-wise multiplication by a scalar. */
      SparseMatrix& scale(const double factor);

      /** Element-wise negation. */
      SparseMatrix& operator-();

      /** Return a sparse matrix with the negated entries. */
      SparseMatrix operator-() const;

      /** Element-wise multiplication by a scalar. */
      SparseMatrix& operator*=(const double factor);

      /** Return a sparse matrix with the entries multiplied by a scalar. */
      SparseMatrix operator*(const double factor) const;

      /** Element-wise division by a non-zero scalar. */
      SparseMatrix& operator/=(const double factor);

      /**
       * Return a sparse matrix with the entries divided by a non-zero scalar.
       */
      SparseMatrix operator/(const double factor) const;

      //################################################## Matrix Addition and
      //                                                   Subtraction

      /**
       * Element-wise multiplication by a scalar and addition by another
       * scaled sparse matrix, i.e. \f$ A = a A + b B \f$.
       *
       * \note To avoid expensive modifications to underlying structure, the
       *       matrices must have the same sparsity pattern.
       */
      SparseMatrix&
      sadd(const double a,
           const double b,
           const SparseMatrix& B);

      /**
       * Element-wise multiplication by a scalar and addition by another sparse
       * matrix, i.e. \f$ A = a A + B \f$.
       *
       * \note To avoid expensive modifications to underlying structure, the
       *       matrices must have the same sparsity pattern.
       */
      SparseMatrix& sadd(const double a, const SparseMatrix& B);

      /**
       * Element-wise addition by a scaled sparse matrix, i.e. \f$ A = A + b B
       * \f$.
       *
       * \note To avoid expensive modifications to underlying structure, the
       *       matrices must have the same sparsity pattern.
       */
      SparseMatrix& add(const double b, const SparseMatrix& B);

      /**
       * Element-wise addition by another sparse matrix.
       *
       * \note To avoid expensive modifications to underlying structure, the
       *       matrices must have the same sparsity pattern.
       */
      SparseMatrix& operator+=(const SparseMatrix& B);

      /** Return the sum of two sparse matrices. */
      SparseMatrix operator+(const SparseMatrix& B) const;

      /**
       * Element-wise subtraction by another sparse matrix.
       *
       * \note To avoid expensive modifications to underlying structure, the
       *       matrices must have the same sparsity pattern.
       */
      SparseMatrix& operator-=(const SparseMatrix& B);

      /** Return the difference between two sparse matrices. */
      SparseMatrix operator-(const SparseMatrix& B) const;

      //################################################## Matrix-Vector
      //                                                   Multiplication

      /**
       * Compute a matrix-vector product, i.e. \f$ y = A x \f$.
       *
       * The optional \p adding flag dictates whether to write or add to the
       * destination vector \p y.
       *
       * \note It is acceptable for the vectors \f$ x \f$ and \f$ y \f$ to be
       *    the same for square matrices.
       */
      void
      vmult(Vector& y,
            const Vector& x,
            const bool adding = false) const;

      /**
       * Add a matrix-vector product to the destination vector.
       *
       * \note It is acceptable for the vectors \f$ x \f$ and \f$ y \f$ to be
       *    the same for square matrices.
       */
      void
      vmult_add(Vector& y, const Vector& x) const;

      /** Return a matrix-vector product. */
      Vector
      operator*(const Vector& x) const;

      /**
       * Compute a transpose matrix-vector product via \f$ y = A^T x \f$
       *
       * The optional \p adding flag dictates whether to write or add to the
       * destination vector \p y.
       *
       * \note It is acceptable for the vectors \f$ x \f$ and \f$ y \f$ to be
       *    the same for square matrices.
       */
      void
      Tvmult(Vector& y,
             const Vector& x,
             const bool adding = false) const;

      /** Add a transpose matrix-vector product to the destination vector. */
      void
      Tvmult_add(Vector& y, const Vector& x) const;

      //################################################## Print Utilities

      /** Return the sparse matrix as a string with the specified formatting. */
      std::string
      str(const bool formatted = true,
          const bool scientific = true,
          const unsigned int precision = 3,
          const unsigned int width = 0) const;

      /**
       *  Print the sparse matrix as triplets of nonzero elements to an output
       *  stream with the specified formatting.
       */
      void
      print(std::ostream& os = std::cout,
            const bool formatted = false,
            const bool scientific = true,
            const unsigned int precision = 3,
            const unsigned int width = 0) const;

      /**
       * Print a row of the sparse matrix to an output stream with the
       * specified formatting.
       */
      void
      print_row(const size_t row,
                std::ostream& os = std::cout,
                const bool scientific = true,
                const unsigned int precision = 3) const;


      //################################################## Comparison

      /** Return whether all entries of two sparse matrices are equivalent. */
      bool operator==(const SparseMatrix& other) const;

      /** Return whether any entries of two sparse matrices are equivalent. */
      bool operator!=(const SparseMatrix& other) const;

      //################################################## Friends

      friend SparseMatrix
      operator*(const double factor, const SparseMatrix& A);

      friend std::ostream&
      operator<<(std::ostream& os, const SparseMatrix& A);
    };


    /** Multiply a sparse matrix by a scalar. */
    SparseMatrix operator*(const double factor, const SparseMatrix& A);


    /** Insert a sparse matrix into an output stream. */
    std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);
  }
}

#endif //SPARSE_MATRIX_H

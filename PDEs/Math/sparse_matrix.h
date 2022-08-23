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
    // forward declarations
    class Vector;

    class Matrix;


    /**
     * Implementation of a list of lists sparse matrix.
     */
    class SparseMatrix
    {
    public:
      /**
       * An implementation of an iterator over the entries of the sparse matrix.
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
          /**
           * The row index of the sparse matrix entry.
           */
          const size_t row;

          /**
           * The column index of the sparse matrix entry.
           */
          const size_t column;

          /**
           * Mutable reference to the sparse matrix entry.
           */
          double& value;

          /**
           * Construct the entry.
           */
          Entry(const size_t row,
                const size_t column,
                double& value);
        };

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
         * entries. This sets the row and index to invalid entries.
         */
        Iterator(SparseMatrix* matrix);

        /**
         * Pre-increment operator. This advances the iterator to the next entry,
         * then returns a reference to it. See \ref advance for the mechanics of
         * incrementing the iterator.
         */
        Iterator&
        operator++();

        /**
         * Post-increment operator. This stores a copy of the current iterator
         * state, advances the current iterator to the next entry, and then
         * returns the previously made copy. See \ref advance for the mechanics of
         * incrementing the iterator.
         */
        Iterator
        operator++(int);

        /**
         * Dereference operator. This returns an Entry object.
         */
        Entry
        operator*() const;

        /**
         * Return whether two iterators reference the same entry.
         */
        bool
        operator==(const Iterator& other) const;

        /**
         * Return whether two iterators reference different entries.
         */
        bool
        operator!=(const Iterator& other) const;

        /**
         * Return whether this iterator references an entry on a previous index
         * of the current row or a previous row.
         */
        bool
        operator<(const Iterator& other) const;


      private:
        /**
         * A utility function to increment the iterator. If the next entry is on
         * the end of a row, the iterator is moved to the next row. If the
         * iterator is currently on the last entry of the last row, the invalid
         * iterator is returned.
         */
        void
        advance();

        /**
         * A pointer to the sparse matrix associated with the iterator.
         */
        SparseMatrix* matrix;

        /**
         * The current row the iterator references.
         */
        size_t row;

        /**
         * The current non-zero index of the row the iterator references.
         */
        unsigned int index;
      };


      /**
       * An implementation of a constant iterator over the entries of the
       * sparse matrix.
       */
      class ConstIterator
      {
      private:
        /**
         * A struct that describes the sparse matrix entry the constant iterator
         * is pointing to and provides a constant reference to the underlying
         * value.
         */
        struct ConstEntry
        {
          /**
           * The row index of the sparse matrix entry.
           */
          const size_t row;

          /**
           * The column index of the sparse matrix entry.
           */
          const size_t column;

          /**
           * Constant reference to the sparse matrix entry.
           */
          const double& value;

          /**
           * Construct the entry.
           */
          ConstEntry(const size_t row,
                     const size_t column,
                     const double& value);
        };

      public:
        /**
         * Construct a constant iterator that references a particular row and
         * non-zero index on that row of a sparse matrix.
         */
        ConstIterator(const SparseMatrix* matrix,
                      const size_t row,
                      const unsigned int index);

        /**
         * Construct a constant iterator to designate the end of the sparse matrix
         * entries. This sets the row and index to invalid entries.
         */
        ConstIterator(const SparseMatrix* matrix);

        /**
         * Pre-increment operator. This advances the constant iterator to the next
         * entry, then returns a reference to it. See \ref advance for the
         * mechanics of incrementing the constant iterator.
         */
        ConstIterator&
        operator++();

        /**
         * Post-increment operator. This stores a copy of the current constant
         * iterator state, advances the current constant iterator to the next
         * entry, and then  returns the previously made copy. See \ref advance
         * for the mechanics of incrementing the constant iterator.
         */
        ConstIterator
        operator++(int);

        /**
         * Dereference operator. This returns a ConstEntry object.
         */
        ConstEntry
        operator*() const;

        /**
         * Return whether two constant iterators reference the same entry.
         */
        bool
        operator==(const ConstIterator& other) const;

        /**
         * Return whether two constant iterators reference different entries.
         */
        bool
        operator!=(const ConstIterator& other) const;

        /**
         * Return whether this constant iterator references an entry on a previous
         * index of the current row or a previous row.
         */
        bool
        operator<(const ConstIterator& other) const;


      private:
        /**
         * A utility function to increment the constant iterator. If the next
         * entry is on the end of a row, the constant iterator is moved to the
         * next row. If the constant iterator is currently on the last entry of
         * the last row, the invalid constant iterator is returned.
         */
        void
        advance();

        /**
         * A constant pointer to the sparse matrix associated with the iterator.
         */
        const SparseMatrix* matrix;

        /**
         * The current row the constant iterator references.
         */
        size_t row;

        /**
         * The current non-zero index of the row the constant iterator references.
         */
        unsigned int index;
      };


      /**
       * An implementation of an iterator over a row of the sparse matrix.
       * This is just a specialization of the Iterator which implements
       * a \p begin and \p end routine.
       */
      class RowIterator
      {
      public:
        /**
         * Construct a row iterator pointing to the specified row.
         */
        RowIterator(SparseMatrix* matrix, const size_t row);

        /**
         * Return an iterator to the first entry of the \p row.
         */
        Iterator
        begin();

        /**
         * Return the iterator to the end of the \p row.
         */
        Iterator
        end();

      private:
        /**
         * A pointer to the sparse matrix the row iterator references.
         */
        SparseMatrix* matix;

        /**
         * The row being iterated over.
         */
        const size_t row;
      };


      /**
       * An implementation of a constant iterator over a row of the sparse matrix.
       * This is just a specialization of the ConstIterator which implements
       * a \p begin and \p end routine.
       */
      class ConstRowIterator
      {
      public:
        /**
         * Construct a constant row iterator pointing to the specified row.
         */
        ConstRowIterator(const SparseMatrix* matrix, const size_t row);

        /**
         * Return a constant iterator to the first entry of the \p row.
         */
        ConstIterator
        begin();

        /**
         * Return the constant iterator to the end of the \p row.
         */
        ConstIterator
        end();

      private:
        /**
         * A constant pointer to the sparse matrix the row iterator references.
         */
        const SparseMatrix* matix;

        /**
         * The row being iterated over.
         */
        const size_t row;
      };

    public:
      /**
       * \name Constructors and initializers
       */
      /* @{ */

      /**
       * Default constructor. Construct an empty sparse matrix.
       */
      SparseMatrix();

      /**
       * Construct a sparse matrix with \p n_rows and \p n_cols.
       */
      SparseMatrix(const size_t n_rows, const size_t n_cols);

      /**
       * Element-wise assignment of the sparse matrix entries to a scalar.
       */
      SparseMatrix&
      operator=(const double value);

      /**
       * Reinitialize the sparse matrix with \p n_rows and \p n_cols. This
       * clears the existing data.
       */
      void
      reinit(const size_t n_rows, const size_t n_cols);

      /* @} */
      /**
       * \name Information about the sparse matrix.
       */
      /* @{ */

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
      n_nonzero_entries() const;

      /**
       * Return the length (number of non-zero entries) of a \p row.
       */
      unsigned int
      row_length(const size_t row) const;

      /**
       * Return the column index for a particular non-zero \p index of a \p row.
       */
      size_t
      column(const size_t row, const unsigned int index) const;

      /**
       * Return the non-zero index on a \p row that corresponds to a \p column.
       * If the \p column does not exist, an invalid value is returned.
       */
      unsigned int
      index(const size_t row, const size_t column) const;

      /**
       * Return whether the sparse matrix is empty or not.
       */
      bool
      empty() const;

      /**
       * Return whether an entry for a \p row and \p column has been allocated.
       */
      bool
      exists(const size_t row, const size_t column) const;

      /**
       * Return whether all entries of two sparse matrices are equivalent.
       */
      bool
      operator==(const SparseMatrix& other) const;

      /**
       * Return whether any entries of two sparse matrices are equivalent.
       */
      bool
      operator!=(const SparseMatrix& other) const;

      /* @} */
      /**
       * \name Accessors and iterators
       * */
      // @{

      /**
       * Return a reference to the value of the sparse matrix at row \p i and
       * column \p j. If the entry does not exist, an error is thrown. When
       * potential uninitialized elements may be encountered, the \ref el method
       * should be used.
       */
      double&
      operator()(const size_t i, const size_t j);

      /**
       * Return a constant reference to the value of the sparse matrix at row \p i
       * and column \p j. If the entry does not exist, an error is thrown. When
       * potential uninitialized elements may be encountered, the \ref el method
       * should be used.
       */
      const double&
      operator()(const size_t i, const size_t j) const;

      /**
       * Return the value of the sparse matrix at row \p i and column \p j. If
       * the entry is uninitialized, zero is returned.
       */
      double
      el(const size_t i, const size_t j);

      /**
       * Return a reference to the <tt>i</tt>'th diagonal entry of the sparse
       * matrix. If the entry does not exist, an error is thrown. The \ref
       * diag_el method should be used if the diagonal may be zero.
       */
      double&
      diag(const size_t i);

      /**
       * Return a constant reference to the <tt>i</tt>'th diagonal entry of the
       * sparse matrix. If the entry does not exist, an error is thrown. The \ref
       * diag_el method should be used if the diagonal may be zero.
       */
      const double&
      diag(const size_t i) const;

      /**
       * Return the value of the <tt>i</tt>'th diagonal entry of the sparse
       * matrix. If the entry does not exist, zero is returned.
       */
      double
      diag_el(const size_t i) const;


      /**
       * Return an iterator to the first entry of the first row.
       */
      Iterator
      begin();

      /**
       * Return an iterator that designates the end of the sparse matrix.
       */
      Iterator
      end();

      /**
       * Return an iterator to the first entry of \p row.
       */
      Iterator
      begin(const size_t row);

      /**
       * Return an iterator to the end of \p row.
       */
      Iterator
      end(const size_t row);

      /**
       * Return a constant iterator to the first entry of the first row.
       */
      ConstIterator
      begin() const;

      /**
       * Return a constant iterator that designates the end of the sparse matrix.
       */
      ConstIterator
      end() const;

      /**
       * Return a constant iterator to the first entry of \p row.
       */
      ConstIterator
      begin(const size_t row) const;

      /**
       * Return a constant iterator to the end of \p row.
       */
      ConstIterator
      end(const size_t row) const;

      /**
       * Return a row iterator for the specified \p row.
       */
      RowIterator
      row_iterator(const size_t row);

      /**
       * Return a constant row iterator for \p row.
       */
      ConstRowIterator
      row_iterator(const size_t row) const;

      /* @} */
      /**
       * \name Modifying the sparse matrix.
       */
      /* @{ */

      /**
       * Delete the contents of the sparse matrix.
       */
      void
      clear();

      /**
       * Copy the non-zero contents of a dense matrix.
       */
      void
      copy_from(const Matrix& matrix);

      /**
       * Copy from another sparse matrix.
       */
      void
      copy_from(const SparseMatrix& matrix);

      /**
       * Set entry for the specified \p row and \p column to \p value. If the
       * entry is initialized, override the value. If it is not, initialize it.
       */
      void
      set(const size_t row, const size_t column, const double value);

      /**
       * Add \p value to the entry at the specified \p row and \p column. If the
       * element is not initialized, then set it, otherwise perform an add
       * operation.
       */
      void
      add(const size_t row, const size_t column, const double value);

      /**
       * Swap the contents of two rows of the sparse matrix.
       */
      void
      swap_row(const size_t i, const size_t k);

      /**
       * Swap the contents of two sparse matrices.
       */
      void
      swap(SparseMatrix& other);

      /* @} */
      /**
       * \name Scaling operations
       */
      /* @{ */

      /**
       * Multiply the entries of sparse matrix by a scalar factor such that \f$
       * A = a A \f$.
       */
      SparseMatrix&
      scale(const double factor);

      /**
       * Negate the entries of the sparse matrix such that \f$ A = -A \f$.
       * This is equivalent to scaling by -1.0. See \ref scale.
       */
      SparseMatrix&
      operator-();

      /**
       * Return a sparse matrix containing the negated entries of this sparse
       * matrix. See \ref scale.
       */
      SparseMatrix
      operator-() const;

      /**
       * Multiply the entries of the sparse matrix by a scalar. See \ref scale.
       */
      SparseMatrix&
      operator*=(const double factor);

      /**
       * Return a sparse matrix containing the entries of this sparse matrix
       * multiplied by a scalar. See \ref scale.
       */
      SparseMatrix
      operator*(const double factor) const;

      /**
       * Divide the entries of the sparse matrix by a non-zero scalar.
       */
      SparseMatrix&
      operator/=(const double factor);

      /**
       * Return a sparse matrix containing the entries of this sparse matrix
       * divided by a non-zero scalar. See \ref scale.
       */
      SparseMatrix
      operator/(const double factor) const;

      /* @} */
      /**
       * \name Addition and subtraction operations
       */
      /* @{ */

      /**
       * Scale this sparse matrix by a scalar and add another scaled sparse matrix
       * to it such that \f$ A = a A + b B \f$.
       *
       * \note To avoid expensive modifications to underlying structure, the
       *       matrices must have the same sparsity pattern.
       */
      SparseMatrix&
      sadd(const double a, const double b, const SparseMatrix& B);

      /**
       * Scale this sparse matrix and add another. This is equivalent to calling
       * \ref sadd with <tt>b = 1.0</tt>. See \ref sadd.
       */
      SparseMatrix&
      sadd(const double a, const SparseMatrix& B);

      /**
       * Add a scaled sparse matrix to this one. This is equivalent to calling
       * \ref sadd with <tt>a = 1.0 </tt>. See \ref sadd.
       */
      SparseMatrix&
      add(const double b, const SparseMatrix& B);

      /**
       * Add a sparse matrix to this one. This is equivalent to calling \ref add
       * with <tt>b = 1.0</tt>. See \ref add.
       */
      SparseMatrix&
      operator+=(const SparseMatrix& B);

      /**
       * Return a sparse matrix containing the sum of this sparse matrix and
       * another. See \ref add.
       */
      SparseMatrix
      operator+(const SparseMatrix& B) const;

      /**
       * Subtract a sparse matrix from this one. This is equivalent to calling
       * \ref add with <tt>b = -1.0</tt>. See \ref add.
       */
      SparseMatrix&
      operator-=(const SparseMatrix& B);

      /**
       * Return a sparse matrix containing the difference between this sparse
       * matrix and another. See \ref add.
       */
      SparseMatrix
      operator-(const SparseMatrix& B) const;

      /* @} */
      /**
       * \name Matrix-vector products
       */
      /* @} */

      /**
       * Compute a matrix-vector product via \f$ y = A x = \sum_j a_{ij} x_j ~
       * \forall i \f$.
       *
       * \param[in] x The multiplying vector.
       * \param[out] y The destination vector.
       * \param adding A flag for adding to or setting the destination vector.
       */
      void
      vmult(const Vector& x, Vector& y, const bool adding = false) const;

      /**
       * Add a matrix-vector product to the destination vector. See \ref vmult.
       *
       * \param[in] x The multiplying vector.
       * \param[out] y The destination vector.
       */
      void
      vmult_add(const Vector& x, Vector& y) const;

      /**
       * Return a matrix-vector product. See \ref vmult.
       */
      Vector
      vmult(const Vector& x) const;

      /**
       * Return a matrix-vector product. See \ref vmult.
       */
      Vector
      operator*(const Vector& x) const;

      /**
       * Compute a transpose matrix-vector product via \f$ y = A^T x = \sum_j
       * a_{ji} x_i ~ \forall i
       *
       * \param[in] x The multiplying vector.
       * \param[out] y The destination vector.
       * \param adding A flag for adding to or setting the destination vector.
       */
      void
      Tvmult(const Vector& x, Vector& y, const bool adding = false) const;

      /**
       * Add a transpose matrix-vector product to the destination vector. See
       * \ref Tvmult.
       *
       * \param[in] x The multiplying vector.
       * \param[out] y The destination vector.
       */
      void
      Tvmult_add(const Vector& x, Vector& y) const;

      /**
       * Return a transpose matrix-vector product. See \ref Tvmult.
       */
      Vector
      Tvmult(const Vector& x) const;

      /* @} */
      /**
       * \name Print utilities
       */
      /* @{ */

      /**
       * Return the sparse matrix as a string.
       *
       * \param scientific A flag for scientific notation.
       * \param precision The precision to display to sparse matrix elements.
       * \param width The spacing between entries.
       */
      std::string
      str(const bool formatted = true,
          const bool scientific = true,
          const unsigned int precision = 3,
          const unsigned int width = 0) const;

      /**
       *  Print the sparse matrix as triplets of nonzero entries.
       *
       * \param os The output stream to print the matrix to.
       * \param formatted A flag to print in dense matrix or sparse format.
       * \param scientific A flag for scientific notation.
       * \param precision The precision to display to sparse matrix elements.
       */
      void
      print(std::ostream& os = std::cout,
            const bool formatted = false,
            const bool scientific = true,
            const unsigned int precision = 3,
            const unsigned int width = 0) const;

      /**
       * Print a row of the sparse matrix.
       *
       * \param row The row to print.
       * \param os The output stream to print the matrix to.
       * \param scientific A flag for scientific notation.
       * \param precision The precision to display to sparse matrix elements.
       */
      void
      print_row(const size_t row,
                std::ostream& os = std::cout,
                const bool scientific = true,
                const unsigned int precision = 3) const;

      /* @} */

    private:
      /**
       * A flag for whether entries have been allocated.
       */
      bool has_entries;

      size_t rows; ///< The number of rows.
      size_t cols; ///< The number of columns.

      /**
       * Row-wise storage of the column indices which have non-zero entries on
       * the particular row.
       */
      std::vector<std::vector<size_t>> colnums;

      /**
       * Row-wise storage of the values of the non-zero entries of the particular
       * row.
       */
      std::vector<std::vector<double>> values;
    };


    /**
     * Insert a sparse matrix into an output stream. This uses default parameters
     * from to \ref SparseMatrix::str routine.
     */
    std::ostream&
    operator<<(std::ostream& os, const SparseMatrix& A);
  }
}

#endif //SPARSE_MATRIX_H

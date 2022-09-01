#ifndef MATRIX_H
#define MATRIX_H

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


    /**
     * Implementation of a general linear algebra dense matrix.
     */
    class Matrix
    {
    public:
      using iterator = std::vector<Vector>::iterator;
      using const_iterator = std::vector<Vector>::const_iterator;

    protected:
      std::vector<Vector> values;

    public:
      //################################################## Constructors

      /** Default constructor. Create an empty matrix. */
      Matrix() = default;

      /** Copy constructor. Copy the internal data from another matrix. */
      Matrix(const Matrix& other) = default;

      /** Move constructor. Steal the internal data from another matrix. */
      Matrix(Matrix&& other) = default;

      /** Construct a matrix with \p n_rows and \p n_cols. */
      Matrix(const size_t n_rows,
             const size_t n_cols);

      /**
       * Construct a matrix with \p n_rows and \p n_cols whose entries are
       * set to \p value.
       */
      Matrix(const size_t n_rows,
             const size_t n_cols,
             const double value);

      /** Construct a matrix from initializer lists. */
      Matrix(const std::initializer_list<std::initializer_list<double>>& list);

      /**
       * Construct with the values pointed to by iterators in the range
       * <tt>[first, last)</tt>.
       *
       * Each iterator must point to another object with defined iterators
       * which all have identical distances from first to last.
       */
      template<typename InputIterator>
      Matrix(const InputIterator first, const InputIterator last);

      /**
       * Construct a matrix with \p n_rows and \p n_cols from <tt>n_rows *
       * n_cols</tt> contiguously stored entries.
       */
      Matrix(const size_t n_rows,
             const size_t n_cols,
             const double* value_ptr);

      /** Copy construction of a diagonal matrix from a Vector. */
      Matrix(const Vector& diagonal);

      /** Move construction of a diagonal matrix from a Vector. */
      Matrix(Vector&& diagonal);

      /** Construct a diagonal matrix from an initializer list. */
      Matrix(const std::initializer_list<double>& diagonal);

      /** Copy assignment from another matrix. */
      Matrix& operator=(const Matrix& other);

      /** Move assignment from another matrix. */
      Matrix& operator=(Matrix&& other);

      /** Entry-wise assignment to a scalar value. */
      Matrix& operator=(const double value);

      //################################################## Capacity

      /** Return the number of rows. */
      size_t
      n_rows() const;

      /** Return the number of columns. */
      size_t n_cols() const;

      /** Return the number of entries. */
      size_t size() const;

      /** Return the number of non-zero entries. */
      size_t n_nonzero_entries() const;

      /** Return whether the matrix is empty (no allocated entries) or not. */
      bool empty() const;

      //################################################## Data Access

      /** Read and write access for row \p i. */
      Vector& operator[](const size_t i);

      /** Read access for row \p i. */
      const Vector& operator[](const size_t i) const;

      /** Read and write access for row \p i. */
      Vector& operator()(const size_t i);

      /** Read access for row \p i. */
      const Vector& operator()(const size_t i) const;

      /** Read and write access for row \p i column \p j. */
      double& operator()(const size_t i, const size_t j);

      /** Read access for row \p i column \p j. */
      const double& operator()(const size_t i, const size_t j) const;

      /** Return a pointer the underlying rows. */
      Vector* data();

      /** Return a constant pointer to the underlying rows. */
      const Vector* data() const;

      /** Return a pointer to the underlying data on row \p i. */
      double* data(const size_t i);

      /** Return a constant pointer to the underlying data on row \p i. */
      const double* data(const size_t i) const;

      /** Return an iterator to the first row. */
      iterator begin();

      /** Return an iterator that designates the end of the rows. */
      iterator end();

      /** Return a constant iterator to the first row of the matrix.*/
      const_iterator begin() const;

      /** Return a constant iterator that designates the end of the rows. */
      const_iterator end() const;

      /** Return an iterator to the first entry of row \p i.*/
      std::vector<double>::iterator
      begin(const size_t i);

      /** Return an iterator that designates the end of row \p i. */
      std::vector<double>::iterator
      end(const size_t i);

      /** Return a constant iterator to the first entry of row \p i. */
      std::vector<double>::const_iterator
      begin(const size_t i) const;

      /** Return a constant iterator that designates the end of row \p i. */
      std::vector<double>::const_iterator
      end(const size_t i) const;

      //################################################## Modifiers

      /** Delete the contents of the matrix. */
      void clear();

      /**
       * Resize the matrix to \p n_rows and \p n_cols, setting new entries
       * to \p value.
       *
       * If either dimension is less than its current size, entries are
       * deleted from the back. If either is greater, new entries are
       * allocated at the back of that dimension and set to \p value.
       */
      void resize(const size_t n_rows, const size_t n_cols);

      /**
       * Resize the matrix to \p n_rows and \p n_cols, setting new entrie
       * to \p value.
       *
       * If either dimension is less than its current size, entries are
       * deleted from the back. If either is greater, new entries are
       * allocated at the back of that dimension and set to \p value.
       */
      void resize(const size_t n_rows,
                  const size_t n_cols,
                  const double value);

      /** Swap two rows of the matrix. */
      void swap_row(const size_t i, const size_t k);

      /** Swap two columns of the matrix. */
      void swap_column(const size_t j, const size_t k);

      /** Swap the contents of two matrices. */
      void swap(Matrix& other);

      /**
       * Set the diagonal of the matrix by copying data from a vector.
       *
       * If the matrix is empty, this creates a diagonal matrix. If the matrix
       * is not empty, the diagonal vector must be of the same size as the
       * minimum dimension of the matrix.
       */
      void set_diagonal(const Vector& diagonal);

      /**
       * Set the diagonal of the matrix by stealing data from a vector.
       *
       * If the matrix is empty, this creates a diagonal matrix. If the matrix
       * is not empty, the diagonal vector must be of the same size as the
       * minimum dimension of the matrix.
       */
       void set_diagonal(Vector&& diagonal);

       /**
        * Set the diagonal of the matrix from an initializer list.
        *
        * If the matrix is empty, this creates a diagonal matrix. If the matrix
        * is not empty, the diagonal vector must be of the same size as the
        * minimum dimension of the matrix.
        */
       void set_diagonal(const std::initializer_list<double>& diagonal);

      /**
       * Set the diagonal of the matrix to a single scalar value.
       *
       * If the matrix is empty, this initializes a matrix with a single entry
       * and sets it to \p value. If not empty, each diagonal entry is set to
       * \p value.
       */
      void set_diagonal(const double value);

      //################################################## Scaling Operations

      /** Entry-wise multiplication by a scalar. */
      Matrix& scale(const double factor);

      /** Entry-wise negation. */
      Matrix& operator-();

      /** Return a matrix with the negated entries. */
      Matrix
      operator-() const;

      /** Entry-wise multiplication by a scalar. */
      Matrix& operator*=(const double factor);

      /**
       * Return a matrix with the entries multiplied by a scalar.
       */
      Matrix operator*(const double factor) const;

      /** Entry-wise division by a non-zero scalar. */
      Matrix& operator/=(const double factor);

      /** Return a matrix with the entries divided by a non-zero scalar. */
      Matrix operator/(const double factor) const;

      //################################################## Matrix Addition and
      //                                                   Subtraction

      /**
       * Entry-wise multiplication by a scalar and addition by another
       * scale matrix, i.e. \f$ A = a A + b B \f$.
       */
      Matrix& sadd(const double a, const double b, const Matrix& B);

      /**
       * Entry-wise multiplication by a scalar and addition by another
       * matrix, i.e. \f$ A = a A + B \f$.
       */
      Matrix& sadd(const double a, const Matrix& B);

      /**
       * Entry-wise addition by a scaled matrix, i.e. \f$ A = A + b B \f$.
       */
      Matrix& add(const double b, const Matrix& B);

      /** Entry-wise addition by another matrix. */
      Matrix& operator+=(const Matrix& B);

      /** Return the sum of two matrices. */
      Matrix operator+(const Matrix& B) const;

      /** Entry-wise subtraction by another matrix. */
      Matrix& operator-=(const Matrix& B);

      /** Return the difference between two matrices. */
      Matrix operator-(const Matrix& B) const;

      /**
       * Entry-wise multiplication by a scalar and addition by a scaled
       * transpose matrix, i.e. \f$ A = a A + b B^T \f$.
       *
       * The transpose operation is applied within this routine by swapping
       * the row and column indices when querying the matrix \p B. The original
       * matrix \p B, and not its transpose, should be passed to this routine.
       */
      Matrix& sTadd(const double a, const double b, const Matrix& B);

      /**
       * Entry-wise multiplication by a scalar and addition by a transpose
       * matrix, i.e. \f$ A = a A + B \f$.
       *
       * The transpose operation is applied within this routine by swapping
       * the row and column indices when querying the matrix \p B. The original
       * matrix \p B, and not its transpose, should be passed to this routine.
       */
      Matrix& sTadd(const double a, const Matrix& B);

      /**
       * Entry-wise addition of a scaled transpose matrix, i.e. \f$ A = A +
       * b B^T \f$.
       *
       * The transpose operation is applied within this routine by swapping
       * the row and column indices when querying the matrix \p B. The original
       * matrix \p B, and not its transpose, should be passed to this routine.
       */
      Matrix& Tadd(const double b, const Matrix& B);

      //################################################## Matrix-Matrix
      //                                                   Multiplication

      /**
       * Compute a matrix-matrix product, i.e. \f$ C = A B \f$.
       *
       * The optional \p adding flag dictates whether to write or add to the
       * destination matrix \p C.
       */
      void
      mmult(Matrix& C,
            const Matrix& B,
            const bool adding = false) const;

      /** Return a matrix-matrix product. */
      Matrix operator*(const Matrix& B) const;

      /**
       * Compute a transpose matrix-matrix product via \f$ C = A^T B \f$.
       *
       * The optional \p adding flag dictates whether to write or add to the
       * destination matrix \p C.
       */
      void
      Tmmult(Matrix& C,
             const Matrix& B,
             const bool adding = false) const;

      /**
       * Compute a matrix-transpose matrix product, i.e. \f$ C = A B^T \f$.
       *
       * The optional \p adding flag dictates whether to write or add to the
       * destination matrix \p C.
       *
       * The transpose operation is applied within this routine by swapping
       * the row and column indices when querying the matrix \p B. The original
       * matrix \p B, and not its transpose, should be passed to this routine.
       */
      void
      mTmult(Matrix& C,
             const Matrix& B,
             const bool adding = false) const;

      /**
       * Compute a transpose matrix-transpose matrix product via \f$ C = A^T B^T
       * \f$.
       *
       * The optional \p adding flag dictates whether to write or add to the
       * destination matrix \p C.
       *
       * The transpose operation is applied within this routine by swapping
       * the row and column indices when querying the matrix \p B. The original
       * matrix \p B, and not its transpose, should be passed to this routine.
       */
      void
      TTmult(Matrix& C,
             const Matrix& B,
             const bool adding = false) const;

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

      /** Add a matrix-vector product to the destination vector \f$ y \f$. */
      void vmult_add(Vector& y, const Vector& x) const;

      /** Return a matrix-vector product. */
      Vector operator*(const Vector& x) const;

      /**
       * Compute a transpose matrix-vector product, i.e. \f$ y = A^T x \f$.
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

      /**
       * Add a transpose matrix-vector product to the destination vector \f$
       * y \f$.
       */
      void Tvmult_add(Vector& y, const Vector& x);

      //################################################## Print Utilities

      /** Return the matrix as a string with the specified formatting. */
      std::string
      str(const bool scientific = true,
          const unsigned int precision = 3,
          const unsigned int width = 0) const;

      /** Print the matrix to an output stream with the specified formatting. */
      void
      print(std::ostream& os = std::cout,
            const bool scientific = true,
            const unsigned int precision = 3,
            const unsigned int width = 0) const;

      //################################################## Comparison

      /** Return whether all entries of two matrices are equivalent. */
      bool operator==(const Matrix& other) const;

      /** Return whether any entries of two matrices are different. */
      bool operator!=(const Matrix& other) const;

      //################################################## Friends

      friend Matrix
      operator*(const double factor, const Matrix& A);

      friend std::ostream&
      operator<<(std::ostream& os, const Matrix& A);
    };


    /** Multiplication by a scalar. */
    Matrix operator*(const double factor, const Matrix& A);


    /** Insert the matrix into an output stream. */
    std::ostream& operator<<(std::ostream& os, const Matrix& A);
  }
}

#endif //MATRIX_H

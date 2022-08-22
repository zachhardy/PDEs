#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


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
    /**
     * Alias for an iterator over an STL vector of Vector objects.
     */
    using iterator = std::vector<Vector>::iterator;

    /**
     * Alias for a constant iterator over an STL vector of Vector objects.
     */
    using const_iterator = std::vector<Vector>::const_iterator;

  public:

    /**
     * \name Constructors and assignment
     */
    /* @{ */

    /**
     * Default constructor. Create an empty matrix.
     */
    Matrix() = default;

    /**
     * Copy constructor. Copy the internal data from another matrix.
     */
    Matrix(const Matrix& other) = default;

    /**
     * Move constructor. Steal the internal data from another matrix.
     */
    Matrix(Matrix&& other) = default;

    /**
     * Construct a matrix with \p n_rows and \p n_cols.
     */
    Matrix(const size_t n_rows,
           const size_t n_cols,
           const double value = 0.0);

    /**
     * Copy construction with nested initializer lists.
     */
    Matrix(const std::initializer_list<std::initializer_list<double>>& list);

    /**
     * Construct a matrix from iterators.
     */
    template<typename InputIterator>
    Matrix(const InputIterator first, const InputIterator last);

    /**
     * Construct a matrix with \p n_rows and \p n_cols from contiguously stored
     * elements.
     */
    Matrix(const size_t n_rows,
           const size_t n_cols,
           const double* value_ptr);

    /**
     * Copy construction of a diagonal matrix from a Vector.
     */
    Matrix(const Vector& diagonal);

    /**
     * Move construction of a diagonal matrix from a Vector.
     */
    Matrix(Vector&& diagonal);

    /**
     * Construct a diagonal matrix from an initializer list of doubles.
     */
    Matrix(const std::initializer_list<double>& diagonal);

    /**
     * Clear the Matrix and reinitialize it with \p n_rows and \p n_cols,
     * optionally set to \p value.
     */
    void
    reinit(const size_t n_rows,
           const size_t n_cols,
           const double value = 0.0);

    /**
     * Copy assignment from another matrix.
     */
    Matrix&
    operator=(const Matrix& other);

    /**
     * Move assignment from another matrix.
     */
    Matrix&
    operator=(Matrix&& other);

    /**
     * Element-wise assignment to a scalar value.
     */
    Matrix&
    operator=(const double value);

    /* @} */
    /**
     * \name Information about the matrix.
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
     * Return the number of elements.
     */
    size_t
    size() const;

    /**
     * Return the number of non-zero elements.
     */
    size_t
    n_nonzero_elements() const;

    /**
     * Return whether the matrix is empty (no allocated elements) or not.
     */
    bool
    empty() const;

    /**
     * Return the transpose of the matrix.
     */
    Matrix
    transpose() const;

    /**
     * Return whether all elements of two matrices are equivalent.
     */
    bool
    operator==(const Matrix& other) const;

    /**
     * Return whether any elements of two matrices are different.
     */
    bool
    operator!=(const Matrix& other) const;

    /* @} */
    /**
     * \name Accessors and iterators
     */
    /* @{ */

    /**
     * Read and write access for row \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    Vector&
    operator[](const size_t i);

    /**
     * Read access for row \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    const Vector&
    operator[](const size_t i) const;

    /**
     * Read and write access for row \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    Vector&
    operator()(const size_t i);

    /**
     * Read access for row \p i.
     * \note No bounds checking is performed. See \ref at.
     */
    const Vector&
    operator()(const size_t i) const;

    /**
     * Read and write access for row \p i with bounds checking.
     */
    Vector&
    at(const size_t i);

    /**
     * Read access for row \p i with bounds checking.
     */
    const Vector&
    at(const size_t i) const;

    /**
     * Read and write access for row \p i column \p j.
     * \note No bounds checking is performed. See \ref at.
     */
    double&
    operator()(const size_t i, const size_t j);

    /**
     * Read access for row \p i column \p j.
     * \note No bounds checking is performed. See \ref at.
     */
    const double&
    operator()(const size_t i, const size_t j) const;

    /**
     * Read and write access for row \p i column \p j with bounds checking.
     */
    double&
    at(const size_t i, const size_t j);

    /**
     * Read access for row \p i column \p j with bounds checking.
     */
    const double&
    at(const size_t i, const size_t j) const;

    /**
     * Read and write access for diagonal element \p i with bounds checking.
     */
    double&
    diag(const size_t i);

    /**
     * Read access for diagonal element \p i with bounds checking.
     */
    const double&
    diag(const size_t i) const;

    /**
     * Return a pointer the underlying rows.
     */
    Vector*
    data();

    /**
     * Return a constant pointer to the underlying rows.
     */
    const Vector*
    data() const;

    /**
     * Return a pointer to the underlying data on row \p i.
     */
    double*
    data(const size_t i);

    /**
     * Return a constant pointer to the underlying data on row \p i.
     */
    const double*
    data(const size_t i) const;

    /**
     * Return an iterator to the first row.
     */
    iterator
    begin();

    /**
     * Return an iterator that designates the end of the rows.
     */
    iterator
    end();

    /**
     * Return a constant iterator to the first row of the matrix.
     */
    const_iterator
    begin() const;

    /**
     * Return a constant iterator that designates the end of the rows.
     */
    const_iterator
    end() const;

    /**
     * Return an iterator to the first element of row \p i.
     */
    std::vector<double>::iterator
    begin(const size_t i);

    /**
     * Return an iterator that designates the end of row \p i.
     */
    std::vector<double>::iterator
    end(const size_t i);

    /**
     * Return a constant iterator to the first element of row \p i.
     */
    std::vector<double>::const_iterator
    begin(const size_t i) const;

    /**
     * Return a constant iterator that designates the end of row \p i.
     */
    std::vector<double>::const_iterator
    end(const size_t i) const;

    /* @} */
    /**
     * \name Modifying the matrix
     */
    /* @} */

    /**
     * Delete the contents of the matrix.
     */
    void
    clear();

    /**
     * Add a row to the back of the matrix. The input Vector have the size as
     * the number of columns in the matrix.
     */
    void
    push_back(const Vector& row);

    /**
     * Add a row to the back of the matrix. The input Vector have the size as
     * the number of columns in the matrix.
     * \note This operation clears the contents of the input Vector.
     */
    void
    push_back(Vector&& row);

    /**
     * Remove the last row from the matrix.
     */
    void
    pop_back();

    /**
     * Resize the matrix to have \p n_rows and \p n_cols. If either dimension
     * is less than the its current size, entries are deleted from the back.
     * If either is greater, new elements are allocated.
     */
    void
    resize(const size_t n_rows,
           const size_t n_cols,
           const double value = 0.0);

    /**
     * Swap two rows of the matrix.
     */
    void
    swap_row(const size_t i, const size_t k);

    /**
     * Swap two columns of the matrix.
     */
    void
    swap_column(const size_t j, const size_t k);

    /**
     * Swap the contents of two matrices.
     */
    void
    swap(Matrix& other);

    /**
     * Set the diagonal of the matrix by copying from a Vector.
     * If the matrix is empty, this creates a diagonal matrix with the
     * specified values. If the matrix is not empty, the diagonal entries must
     * be the same size as the minimum dimension of the matrix.
     */
    void
    set_diag(const Vector& diagonal);

    /**
     * Set the diagonal of the matrix by moving from a Vector.
     * See \ref set_diag.
     */
    void
    set_diag(Vector&& diagonal);

    /**
     * Set the diagonal of the matrix with an initializer list of doubles.
     * See \ref set_diag.
     */
    void
    set_diag(const std::initializer_list<double>& diagonal);

    /**
     * Set the diagonal of the matrix to a single scalar value. If the matrix
     * is empty, this initializes a matrix with one row and one column whose
     * entry is set to \p value. If not empty, each diagonal element is set to
     * \p value.
     */
    void
    set_diag(const double value);

    /* @} */
    /**
     * \name Scaling operations
     */
    /* @{ */

    /**
     * Multiply the matrix by a scalar such that \f$ \boldsymbol{A} = a
     * \boldsymbol{A} \f$.
     */
    Matrix&
    scale(const double factor);

    /**
     * Negate the elements of the matrix such that \f$ \boldsymbol{A} = -
     * \boldsymbol{A} \f$. This is equivalent to scaling by -1.0. See \ref
     * scale.
     */
    Matrix&
    operator-();

    /**
     * Return a matrix containing the negated elements of this matrix. See
     * \ref scale.
     */
    Matrix
    operator-() const;

    /**
     * Multiply the elements of the matrix by a scalar. See \ref scale.
     */
    Matrix&
    operator*=(const double factor);

    /**
     * Return a matrix containing the elements of this matrix multiplied by
     * a scalar. See \ref scale.
     */
    Matrix
    operator*(const double factor) const;

    /**
     * Divide the elements of the matrix by a non-zero scalar. See \ref scale.
     */
    Matrix&
    operator/=(const double factor);

    /**
     * Return a matrix containing the elements of this matrix divided by a
     * non-zero scalar. See \ref scale.
     */
    Matrix
    operator/(const double factor) const;

    /* @} */
    /**
     * \name Addition and subtraction operations
     */
    /* @{ */

    /**
     * Addition by a scaled matrix to this one such that  \f$ \boldsymbol{A} =
     * \boldsymbol{A} + b \boldsymbol{B} \f$. The dimensions of each matrix
     * must agree for this operation to be permissible.
     */
    Matrix&
    add(const double b, const Matrix& B);

    /**
     * Add another matrix to this one. This is equivalent to adding a matrix
     * scaled by 1.0. See \ref add.
     */
    Matrix&
    operator+=(const Matrix& B);

    /**
     * Return the sum of this matrix and another. See \ref add.
     */
    Matrix
    operator+(const Matrix& B) const;

    /**
     * Subtract another matrix from this one. This is equivalent to adding a
     * matrix scaled by -1.0. See \ref add.
     */
    Matrix&
    operator-=(const Matrix& B);

    /**
     * Return the difference between this matrix and another.
     */
    Matrix
    operator-(const Matrix& B) const;

    /**
     * Add the scaled transpose of another matrix to this one such that
     * \f$ \boldsymbol{A} = \boldsymbol{A} + b \boldsymbol{B}^T \f$. The
     * dimensions of this matrix and the transpose of the other must agree in
     * order for this operation to be permissible.
     */
    Matrix&
    Tadd(const double b, const Matrix& B);

    /**
     * Multiply this matrix by a scalar and add another scaled matrix to it
     * such that \f$ \boldsymbol{A} = a \boldsymbol{A} + b \boldsymbol{B} \f$.
     */
    Matrix&
    sadd(const double a, const double b, const Matrix& B);

    /**
     * Multiply this matrix by a scalar and add another to it such that \f$
     * \boldsymbol{A} = a \boldsymbol{A} + \boldymbol{B}. This is equivalent
     * to scaling the other matrix by 1.0 using the more general \ref sadd
     * routine. See \ref sadd.
     */
    Matrix&
    sadd(const double a, const Matrix& B);

    /**
     * Multiply this matrix by a scalar and add the scaled transpose of another
     * to it such that \f$ \boldsymbol{A} = a \boldymbol{A} + b \boldsymbol{B}^T
     * \f$.
     */
    Matrix&
    sTadd(const double a, const double b, const Matrix& B);

    /**
     * Multiply this matrix by a scalar and add the transpose of another to it
     * such that \boldsymbol{A} = a \boldsymbol{A} + \boldsymbol{B}^T \f$. This
     * is equivalent to scaling the other matrix by 1.0 using the more general
     * \ref sTadd routine. See \ref sTadd.
     */
    Matrix&
    sTadd(const double a, const Matrix& B);

    /* @} */
    /**
     * \name Matrix-matrix products
     */
    /* @{ */

    /**
     * Compute a matrix-matrix product via \f$ \boldsymbol{C} = \boldsymbol{A}
     * \boldsymbol{B} = \sum_{k=0}^{n} a_{ik} b_{kj}, ~ \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix.
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void
    mmult(const Matrix& B, Matrix& C, const bool adding = false) const;

    /**
     * Return a matrix-matrix product. See \ref mmult.
     */
    Matrix
    mmult(const Matrix& B) const;

    /**
     * Compute a transpose matrix-matrix product via \f$ \boldsymbol{C} =
     * \boldsymbol{A}^T \boldsymbol{B} = \sum_{k=0}^{n} a_{ki} b_{kj}, ~
     * \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix.
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void
    Tmmult(const Matrix& B, Matrix& C, const bool adding = false) const;

    /**
     * Return a transpose matrix-matrix product. See \ref Tmmult.
     */
    Matrix
    Tmmult(const Matrix& B) const;

    /**
     * Compute a matrix-transpose matrix product via \f$ \boldsymbol{C} =
     * \boldsymbol{A} \boldsymbol{B}^T = \sum_{k=1}^{n} a_{ik} b_{jk}, ~
     * \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix (not transposed).
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void
    mTmult(const Matrix& B, Matrix& C, const bool adding = false) const;

    /**
     * Return a matrix-transpose matrix product. See \ref mTmult.
     */
    Matrix
    mTmult(const Matrix& B) const;

    /**
     * Compute a transpose matrix-transpose matrix product via \f$
     * \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T = \sum_{k=1}^{n}
     * a_{ki} b_{jk}, ~ \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix (not transposed).
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void
    TTmult(const Matrix& B, Matrix& C, const bool adding = false) const;

    /**
     * Return a transpose matrix-transpose matrix product. See \ref TTmult.
     */
    Matrix
    TTmult(const Matrix& B) const;

    /* @} */
    /**
     * \name Matrix-vector products
     */
    /* @{ */

    /**
     * Compute a matrix-vector product via \f$ \vec{y} = \boldsymbol{A} \vec{x}
     * = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i \f$.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     * \param adding A flag for adding to or setting the destination Vector.
     */
    void
    vmult(const Vector& x, Vector& y, const bool adding = false) const;

    /**
     * Return a matrix-vector product. See \ref vmult.
     */
    Vector
    vmult(const Vector& x) const;

    /**
     * Add a matrix-vector product to the destination Vector. See \ref vmult.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     */
    void
    vmult_add(const Vector& x, Vector& y) const;

    /**
     * Return a matrix-vector product. See \ref vmult.
     */
    Vector
    operator*(const Vector& x) const;

    /**
     * Compute a transpose matrix-vector product via \f$ \vec{y} =
     * \boldsymbol{A}^T \vec{x} = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall  i \f$.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     * \param adding A flag for adding to or setting the destination Vector.
     */
    void
    Tvmult(const Vector& x, Vector& y, const bool adding = false) const;

    /**
     * Return a transpose matrix-vector product. See \ref Tvmult.
     */
    Vector
    Tvmult(const Vector& x) const;

    /**
     * Add a transpose matrix-vector product to the destination vector. See
     * \ref Tvmult.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     */
    void
    Tvmult_add(const Vector& x, Vector& y);

    /* @} */
    /**
     * \name Print Utilities
     */
    /* @{ */

    /**
     * Return the matrix as a string.
     *
     * \param scientific A flag for scientific notation.
     * \param precision The precision of the digits to display.
     * \param width The spacing between entries.
     */
    std::string
    str(const bool scientific = true,
        const unsigned int precision = 3,
        const unsigned int width = 0) const;

    /**
     * Return the matrix as a string. See \ref str.
     *
     * \param os The output stream to print the matrix to.
     * \param scientific A flag for scientific notation.
     * \param precision The precision of the digits to display.
     * \param width The spacing between entries.
     */
    void
    print(std::ostream& os = std::cout,
          const bool scientific = true,
          const unsigned int precision = 3,
          const unsigned int width = 0) const;

    /* @} */

  protected:
    /**
     * The underlying matrix data as an STL vector of Vector objects.
     */
    std::vector<Vector> values;
  };


  /**
   * Compute a matrix-matrix product. See \ref Matrix::mmult.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void
  mmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Return a matrix-matrix product. See \ref Matrix::mmult.
   */
  Matrix
  mmult(const Matrix& A, const Matrix& B);

  /**
   * Compute a transpose matrix-matrix product. See \ref Matrix::Tmmult.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void
  Tmmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Return a transpose matrix-matrix product. See \ref Matrix::Tmmult.
   */
  Matrix
  Tmmult(const Matrix& A, const Matrix& B);

  /**
   * Compute a matrix-transpose matrix product. See \ref Matrix::mTmult.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void
  mTmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Return a matrix-transpose matrix product. See \ref Matrix::mTmult.
   */
  Matrix
  mTmult(const Matrix& A, const Matrix& B);

  /**
   * Compute a transpose matrix-transpose matrix product.
   * See \ref Matrix::TTmult.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void
  TTmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Return a transpose matrix-transpose matrix product.
   * See \ref Matrix::TTmult.
   */
  Matrix
  TTmult(const Matrix& A, const Matrix& B);

  /**
   * Compute a matrix-vector product. See \ref Matrix::vmult.
   *
   * \param[in] A The Matrix.
   * \param[in] x The multiplying Vector.
   * \param[out] y The destination Vector.
   */
  void
  vmult(const Matrix& A, const Vector& x, Vector& y);

  /**
   * Return a matrix-vector product. See \ref Matrix::vmult.
   */
  Vector
  vmult(const Matrix& A, const Vector& x);

  /**
   * Compute a transpose matrix-vector product. See \ref Matrix::Tvmult.
   *
   * \param[in] A The Matrix.
   * \param[in] x The multiplying Vector.
   * \param[out] y The destination Vector.
   *
   */
  void
  Tvmult(const Matrix& A, const Vector& x, Vector& y);

  /**
   * Return a transpose matrix-vector product. See \ref Matrix::Tvmult.
   */
  Vector
  Tvmult(const Matrix& A, const Vector& x);

  /**
   * Insert the matrix into an output stream.
   */
  std::ostream&
  operator<<(std::ostream& os, const Matrix& A);
}

#endif //MATRIX_H

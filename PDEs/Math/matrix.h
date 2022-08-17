#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace Math
{
  //########## Forward declarations
  class Vector;


  /** Implementation of a general linear algebra dense matrix. */
  class Matrix
  {
  public:
    using iterator = std::vector<Vector>::iterator;
    using const_iterator = std::vector<Vector>::const_iterator;

    using STLMatrix = std::vector<std::vector<double>>;
    using InitializerMatrix = std::initializer_list<std::initializer_list<double>>;

  protected:
    std::vector<Vector> vals;

  public:

    //################################################## Constructors

    /** \name Constructors and Initialization */
    // @{

    Matrix() = default;

    /** Construct a matrix with \p n_rows and \p n_cols. */
    explicit Matrix(const size_t n_rows, const size_t n_cols);

    /** Construct a Matrix with \p n_rows and \p n_cols set to \p value. */
    explicit Matrix(const size_t n_rows,
                    const size_t n_cols,
                    const double value);

    Matrix(const STLMatrix& other);
    Matrix(STLMatrix&& other);

    Matrix(const InitializerMatrix list);

    Matrix& operator=(const STLMatrix& other);
    Matrix& operator=(STLMatrix&& other);

    /** Element-wise assignment to a scalar. */
    Matrix& operator=(const double value);

    // @}

    //################################################## Characteristics

    /** \name Characteristics */
    // @{

    size_t n_rows() const;
    size_t n_cols() const;
    size_t size() const;

    /** Return the number of non-zero elements in the matrix. */
    size_t nnz() const;

    bool empty() const;
    bool all_zero() const;

    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;

    // @}

    //################################################## Accessors

    /** \name Accessors */
    // @{

    Vector& operator[](const size_t i);
    const Vector& operator[](const size_t i) const;

    Vector& operator()(const size_t i);
    const Vector& operator()(const size_t i) const;

    Vector& at(const size_t i);
    const Vector& at(const size_t i) const;

    double& operator()(const size_t i, const size_t j);
    const double& operator()(const size_t i, const size_t j) const;

    double& at(const size_t i, const size_t j);
    const double& at(const size_t i, const size_t j) const;

    double& diag(const size_t i);
    const double& diag(const size_t i) const;

    Vector* data();
    const Vector* data() const;

    double* data(const size_t i);
    const double* data(const size_t i) const;

    // @}

    //################################################## Iterators

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    std::vector<double>::iterator begin(const size_t i);
    std::vector<double>::iterator end(const size_t i);

    std::vector<double>::const_iterator begin(const size_t i) const;
    std::vector<double>::const_iterator end(const size_t i) const;

    // @}

    //################################################## Modifiers

    /** \name Modifiers */
    // @{

    void clear();

    void pop_back();

    void push_back(const Vector& row);
    void push_back(Vector&& row);

    /**
     * Resize the Matrix to have \p n_rows and \p n_cols. If either dimension 
     * is less than the its current size, entries are deleted from the back.
     * If either is greater, new elements are unininitialized.
     */
    void resize(const size_t n_rows,
                const size_t n_cols);

    /**
     * Resize the Matrix to have \p n_rows and \p n_cols. If new elements are
     * created, they are set to \p value.
     */
    void resize(const size_t n_rows,
                const size_t n_cols,
                const double value);

    /**
     * Clear the Matrix and reinitialize it with \p n_rows and \p n_cols whose
     * elements are uninitialized.
     */
    void reinit(const size_t n_rows,
                const size_t n_cols);

    /** Clear the Matrix and reinitialize it with \p n_rows and \p n_cols whose
     * elements are set to \p value.
     */
    void reinit(const size_t n_rows,
                const size_t n_cols,
                const double value);

    void swap_row(const size_t i, const size_t k);
    void swap_column(const size_t j, const size_t k);
    void swap(Matrix& other);

    void set_diag(const Vector& diag);
    void set_diag(const double value);

    // @}

    //################################################## Linear Algebra

    /** \name Linear Algebra */
    // @{

    /** Element-wise multiplication by a scalar in place. */
    Matrix& scale(const double factor);

    /** Element-wise addition of a scaled Matrix in place. */
    Matrix& add(const Matrix& B, const double a = 1.0);

    /**
     * Element-wise addition of the transpose of a scaled Matrix in place.
     */
    Matrix& Tadd(const Matrix& B, const double a = 1.0);

    /**
     * Element-wise multiplication by a scalar followed by element-wise addition
     * of a Matrix in place.
     */
    Matrix& sadd(const double a, const Matrix& B);

    /**
     * Element-wise multiplication by a scalar \p a followed by element-wise
     * addition of a Matrix scaled by \p b in place.
     */
    Matrix& sadd(const double a, const double b, const Matrix& B);

    /**
     * Element-wise multiplication by a scalar \p a followed by element-wise
     * addition of the transpose of a Matrix in place.
     */
    Matrix& sTadd(const double a, const Matrix& B);

    /**
     * Element-wise multiplication by a scalar \p a followed by element-wise
     * addition of a Matrix scaled by \p b in place.
     */
    Matrix& sTadd(const double a, const double b, const Matrix& B);

    /**
     * Compute a matrix-matrix product via \f$ \boldsymbol{C} = \boldsymbol{A}
     * \boldsymbol{B} = \sum_{k=0}^{n} a_{ik} b_{kj}, ~ \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix.
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void mmult(const Matrix& B, Matrix& C,
               const bool adding = false) const;

    /**
     * Compute a transpose matrix-matrix product via \f$ \boldsymbol{C} =
     * \boldsymbol{A}^T \boldsymbol{B} = \sum_{k=0}^{n} a_{ki} b_{kj}, ~
     * \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix.
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void Tmmult(const Matrix& B, Matrix& C,
                const bool adding = false) const;

    /**
     * Compute a matrix-transpose matrix product via \f$ \boldsymbol{C} =
     * \boldsymbol{A} \boldsymbol{B}^T = \sum_{k=1}^{n} a_{ik} b_{jk}, ~
     * \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix (not transposed).
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void mTmult(const Matrix& B, Matrix& C,
                const bool adding = false) const;

    /**
     * Compute a transpose matrix-transpose matrix product via \f$
     * \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T = \sum_{k=1}^{n}
     * a_{ki} b_{jk}, ~ \forall i, j \f$.
     *
     * \param[in] B The multiplying Matrix (not transposed).
     * \param[out] C The destination Matrix.
     * \param adding A flag for adding to or setting the destination Matrix.
     */
    void TTmult(const Matrix& B, Matrix& C,
                const bool adding = false) const;

    /** Return a matrix-matrix product. */
    Matrix mmult(const Matrix& B) const;

    /** Return a transpose matrix-matrix product. */
    Matrix Tmmult(const Matrix& B) const;

    /** Return a matrix-transpose matrix product. */
    Matrix mTmult(const Matrix& B) const;

    /** Return a transpose matrix-transpose matrix product. */
    Matrix TTmult(const Matrix& B) const;

    /**
     * Compute a matrix-vector product via \f$ \vec{y} = \boldsymbol{A} \vec{x}
     * = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i \f$.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     * \param adding A flag for adding to or setting the destination Vector.
     */
    void vmult(const Vector& x, Vector& y,
               const bool adding = false) const;

    /**
     * Compute a transpose matrix-vector product via \f$ \vec{y} =
     * \boldsymbol{A}^T \vec{x} = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall  i \f$.
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
     * Add a transpose matrix-vector product to the destination vector.
     *
     * \param[in] x The multiplying Vector.
     * \param[out] y The destination Vector.
     */
    void Tvmult_add(const Vector& x, Vector& y);

    /** Element-wise negation in place. */
    Matrix& operator-();

    /** Return a Matrix with the negated elements. */
    Matrix operator-() const;

    /** Element-wise multiplication by a scalar in place. */
    Matrix& operator*=(const double factor);

    /** Element-wise division by a scalar in place. */
    Matrix& operator/=(const double factor);

    /** Element-wise addition with a Matrix in place. */
    Matrix& operator+=(const Matrix& B);

    /** Element-wise subtraction with a Matrix in place. */
    Matrix& operator-=(const Matrix& B);

    /** Return a matrix-vector product. */
    Vector operator*(const Vector& x) const;

    /** Return the transpose. */
    Matrix transpose() const;

    // @}

    //================================================== Print Utilities

    /** \name Print Utilities */
    // @{

    /**
     * Return the matrix as a string.
     *
     * \param os The output stream to print the matrix to.
     * \param scientific A flag for scientific notation.
     * \param precision The precision of the digits to display.
     * \param width The spacing between entries.
     */
    void print(std::ostream& os = std::cout,
               const bool scientific = true,
               const unsigned int precision = 3,
               const unsigned int width = 0) const;

    /**
     * Return the matrix as a string.
     *
     * \param scientific A flag for scientific notation.
     * \param precision The precision of the digits to display.
     * \param width The spacing between entries.
     */
    std::string str(const bool scientific = true,
                    const unsigned int precision = 3,
                    const unsigned int width = 0) const;

    // @}

  private:
    static bool valid_dimensions(const STLMatrix& A);
  };

  //################################################## Methods

  /** Element-wise addition. */
  Matrix operator+(const Matrix& A, const Matrix& B);

  /** Element-wise subtraction. */
  Matrix operator-(const Matrix& A, const Matrix& B);

  /** Return a matrix-matrix product. */
  Matrix operator*(const Matrix& A, const Matrix& B);

  /**
   * Compute a matrix-matrix product.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void mmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Compute a transpose matrix-matrix product.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void Tmmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Compute a matrix-transpose matrix product.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void mTmult(const Matrix& A, const Matrix& B, Matrix& C);

  /**
   * Compute a transpose matrix-transpose matrix product.
   *
   * \param[in] A The left Matrix.
   * \param[in] B The right multiplying Matrix.
   * \param[out] C The destination Matrix.
   */
  void TTmult(const Matrix& A, const Matrix& B, Matrix& C);

  /** Return a matrix-matrix product. */
  Matrix mmult(const Matrix& A, const Matrix& B);

  /** Return a transpose matrix-matrix product. */
  Matrix Tmmult(const Matrix& A, const Matrix& B);

  /** Return a matrix-transpose matrix product. */
  Matrix mTmult(const Matrix& A, const Matrix& B);

  /** Return a transpose matrix-transpose matrix product. */
  Matrix TTmult(const Matrix& A, const Matrix& B);

  /**
   * Compute a matrix-vector product.
   *
   * \param[in] A The Matrix.
   * \param[in] x The multiplying Vector.
   * \param[out] y The destination Vector.
   */
  void vmult(const Matrix& A, const Vector& x, Vector& y);

  /**
   * Compute a transpose matrix-vector product.
   *
   * \param[in] A The Matrix.
   * \param[in] x The multiplying Vector.
   * \param[out] y The destination Vector.
   *
   */
  void Tvmult(const Matrix& A, const Vector& x, Vector& y);

  /** Return a matrix-vector product. */
  Vector vmult(const Matrix& A, const Vector& x);

  /** Return a transpose matrix-vector product. */
  Vector Tvmult(const Matrix& A, const Vector& x);

  std::ostream& operator<<(std::ostream& os, const Matrix& A);

}

#endif //MATRIX_H

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
    using value_type = double;

    using iterator = std::vector<Vector>::iterator;
    using const_iterator = std::vector<Vector>::const_iterator;

    using STLMatrix = std::vector<std::vector<value_type>>;
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
                    const value_type value);

    /** Copy construction with nested STL vectors. */
    Matrix(const STLMatrix& other);

    /** Move construction from nested STL vectors. */
    Matrix(STLMatrix&& other);

    /** Construction from nested initializer lists. */
    Matrix(const InitializerMatrix list);

    /** Copy assignment from nested STL vectors. */
    Matrix& operator=(const STLMatrix& other);

    /** Move assignment from nested STL vectors. */
    Matrix& operator=(STLMatrix&& other);

    /** Set the elements to a scalar value. */
    Matrix& operator=(const value_type value);

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

    value_type& operator()(const size_t i, const size_t j);
    const value_type& operator()(const size_t i, const size_t j) const;

    value_type& at(const size_t i, const size_t j);
    const value_type& at(const size_t i, const size_t j) const;

    value_type& diag(const size_t i);
    const value_type& diag(const size_t i) const;

    Vector* data();
    const Vector* data() const;

    value_type* data(const size_t i);
    const value_type* data(const size_t i) const;

    // @}

    //################################################## Iterators

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    std::vector<value_type>::iterator begin(const size_t i);
    std::vector<value_type>::iterator end(const size_t i);

    std::vector<value_type>::const_iterator begin(const size_t i) const;
    std::vector<value_type>::const_iterator end(const size_t i) const;

    // @}

    //################################################## Modifiers

    /** \name Modifiers */
    // @{

    void clear();

    void pop_back();

    void push_back(const Vector& row);
    void push_back(Vector&& row);

    void resize(const size_t n_rows,
                const size_t n_cols);

    void resize(const size_t n_rows,
                const size_t n_cols,
                const value_type value);

    void reinit(const size_t n_rows,
                const size_t n_cols);

    void reinit(const size_t n_rows,
                const size_t n_cols,
                const value_type value);

    void swap_row(const size_t i, const size_t k);
    void swap_column(const size_t j, const size_t k);
    void swap(Matrix& other);

    void set_diag(const Vector& diag);
    void set_diag(const value_type value);

    // @}

    //################################################## Linear Algebra

    /** \name Linear Algebra */
    // @{

    /**
     * Scale the matrix by a scalar value. This is computed via \f$
     * \boldsymbol{A} = \alpha \boldsymbol{A} = \alpha a_{ij}, ~ \forall i, j \f$.
     */
    Matrix& scale(const value_type factor);

    /**
     * Add a multiple of another matrix. This is computed via \f$ \boldsymbol{A} =
     * \boldsymbol{A} + \alpha \boldsymbol{B} = a_{ij} \alpha b_{ij} \f$. The
     * default behavior is \f$ \alpha = 1.0 \f$.
     */
    Matrix& add(const Matrix& B, const value_type a = 1.0);

    /**
     * Add a multiple of the transpose of another matrix. This is computed via
     * \f$ \boldsymbol{A} = \boldsymbol{A} + \alpha \boldsymbol{B}^T = a_{ij} +
     * \alpha b{ji}, ~ \forall i, j \f$. The default behavior uses \f$ \alpha =
     * 1.0 \f$.
     */
    Matrix& Tadd(const Matrix& B, const value_type a = 1.0);

    /**
     * Scale this matrix and then add another. This is computed via \f$
     * \boldsymbol{A} = \alpha \boldsymbol{A} + \boldsymbol{B} = \alpha a_{ij} +
     * b_{ij}, ~ \forall i, j \f$.
     */
    Matrix& sadd(const value_type a, const Matrix& B);

    /**
     * Scale this matrix and then add a multiple of another matrix. This is
     * computed via \f$ \boldsymbol{A} = \alpha \boldsymbol{A} + \beta
     * \boldsymbol{B} = \alpha a_{ij} + \beta b_{ij}, ~ \forall i, j \f$.
     */
    Matrix& sadd(const value_type a, const value_type b, const Matrix& B);

    /**
     * Scale this matrix and then add the transpose of another. This is computed
     * via \f$ \boldsymbol{A} = \alpha \boldsymbol{A} + \boldsymbol{B}^T = \alpha
     * a_{ij} + b_{ji}, ~ \forall i, j \f$.
     */
    Matrix& sTadd(const value_type a, const Matrix& B);

    /**
     * Scale this matrix and then add a multiple of the transpose another. This
     * is computed via \f$ \boldsymbol{A} = \alpha \boldsymbol{A} + \beta
     * \boldsymbol{B}^T = \alpha a_{ij} + \beta \beta b_{ji}, ~ \forall i, j \f$.
     */
    Matrix& sTadd(const value_type a, const value_type b, const Matrix& B);

    /**
     * Compute a matrix-matrix product. This is computed via
     * \f$ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}  = \sum_{k=0}^{n}
     * a_{ik} b_{kj}, ~ \forall i, j \f$. If desired, setting the \p adding flag
     * to \p true will add the matrix-matrix product to the destination
     * matrix \f$ \boldsymbol{C} \f$.
     */
    void mmult(const Matrix& B, Matrix& C,
               const bool adding = false) const;

    /**
     * Compute a transpose matrix-matrix product. This is computed via
     * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B} = \sum_{k=0}^{n}
     * a_{ki} b_{kj}, ~ \forall i, j \f$. If desired, setting the \p adding flag
     * to \p true will add the matrix-matrix product to the destination
     * matrix \f$ \boldsymbol{C} \f$.
     */
    void Tmmult(const Matrix& B, Matrix& C,
                const bool adding = false) const;

    /**
     * Compute a matrix-transpose matrix product. This is computed via
     * \f$ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}^T = \sum_{k=1}^{n}
     * a_{ik} b_{jk}, ~ \forall i, j \f$. If desired, setting the \p adding flag
     * to \p true will add the matrix-matrix product to the destination
     * matrix \f$ \boldsymbol{C} \f$.
     */
    void mTmult(const Matrix& B, Matrix& C,
                const bool adding = false) const;

    /**
     * Compute a transpose matrix-transpose matrix product. This is
     * computed via \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T =
     * \sum_{k=1}^{n} a_{ki} b_{jk}, \forall i, j \f$. If desired, setting the
     * \p adding flag to \p true will add the matrix-matrix product to
     * the destination matrix \f$ \boldsymbol{C} \f$.
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
     * Compute a matrix-vector product. This is computed via \f$ \vec{y} =
     * \boldsymbol{A} \vec{x} = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i \f$.
     * If desired, setting the \p adding flag to \p true will add the
     * matrix-vector product to the destination vector \f$ \vec{y} \f$.
     */
    void vmult(const Vector& x, Vector& y,
               const bool adding = false) const;

    /**
     * Compute a transpose matrix-vector product. This is computed via \f$
     * \vec{y} = \boldsymbol{A}^T \vec{x} = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall
     * i \f$. If desired, setting the \p adding flag to \p true will add the
     * matrix-vector product to the destination vector \f$ \vec{y} \f$.
     */
    void Tvmult(const Vector& x, Vector& y,
                const bool adding = false) const;

    /** Return a matrix-vector product. */
    Vector vmult(const Vector& x) const;

    /** Return a transpose matrix-vector product. */
    Vector Tvmult(const Vector& x) const;

    /**
     * Add a matrix-vector product to the destination vector \f$ \vec{y} \f$.
     *
     * \see Matrix::vmult
     */
    void vmult_add(const Vector& x, Vector& y) const;

    /**
     * Add a transpose matrix-vector product to the destination vector
     * \f$ \vec{y} \f$
     *
     * \see Matrix::Tvmult
     */
    void Tvmult_add(const Vector& x, Vector& y);

    /**
     * Negate the elements of the matrix. This is computed via \f$ \boldsymbol{A}
     * = -\boldsymbol{A} = -a_{ij}, ~ \forall i, j \f$.
     */
    Matrix& operator-();

    /**
     * Return a matrix containing the negated elements of this matrix.
     *
     * \see Matrix::operator-()
     */
    Matrix operator-() const;

    /**
     * Multiply the elements of the matrix by a scalar factor. This is computed via
     * \f$ \boldsymbol{A} = \alpha \boldsymbol{A} = \alpha a_{ij} ~ \forall i, j
     * \f$.
     *
     * \see Matrix::scale
     */
    Matrix& operator*=(const value_type factor);

    /**
     * Divide the elements of the matrix by a scalar value. This is computed via
     * \f$ \boldsymbol{A} = \frac{\boldsymbol{A}}{\alpha} = \frac{a_{ij}}{\alpha},
     * ~ \forall i, j \f$.
     *
     * \see Matrix::scale
     */
    Matrix& operator/=(const value_type factor);

    /**
     * Add another matrix. This is computed via \f$ \boldsymbol{A} =
     * \boldsymbol{A} + \boldsymbol{B} = a_{ij} + b_{ij}, ~ \forall i, j \f$.
     *
     * \see Matrix::add
     */
    Matrix& operator+=(const Matrix& B);

    /**
     * Subtract another matrix. This is computed via \f$ \boldsymbol{A} =
     * \boldsymbol{A} - \boldsymbol{B} = a_{ij} - b_{ij}, ~ \forall i, j \f$.
     *
     * \see Matrix::add
     */
    Matrix& operator-=(const Matrix& B);

    /**
     * Compute a matrix-vector product.
     *
     * \see Matrix::vmult
     */
    Vector operator*(const Vector& x) const;

    /**
     * Return the transpose of the matrix.
     */
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

  /** Add two matrices together. */
  Matrix operator+(const Matrix& A, const Matrix& B);

  /** Subtract two matrices. */
  Matrix operator-(const Matrix& A, const Matrix& B);

  /** Return a matrix-matrix product. */
  Matrix operator*(const Matrix& A, const Matrix& B);

  /** Compute a matrix-matrix product. */
  void mmult(const Matrix& A, const Matrix& B, Matrix& C);

  /** Compute a transpose matrix-matrix product. */
  void Tmmult(const Matrix& A, const Matrix& B, Matrix& C);

  /** Compute a matrix-transpose matrix product. */
  void mTmult(const Matrix& A, const Matrix& B, Matrix& C);

  /** Compute a transpose matrix-transpose matrix product. */
  void TTmult(const Matrix& A, const Matrix& B, Matrix& C);

  /** Return a matrix-matrix product. */
  Matrix mmult(const Matrix& A, const Matrix& B);

  /** Return a transpose matrix-matrix product. */
  Matrix Tmmult(const Matrix& A, const Matrix& B);

  /** Return a matrix-transpose matrix product. */
  Matrix mTmult(const Matrix& A, const Matrix& B);

  /** Return a transpose matrix-transpose matrix product. */
  Matrix TTmult(const Matrix& A, const Matrix& B);

  /** Compute a matrix-vector product. */
  void vmult(const Matrix& A, const Vector& x, Vector& y);

  /** Compute a transpose matrix-vector product. */
  void Tvmult(const Matrix& A, const Vector& x, Vector& y);

  /** Return a matrix-vector product. */
  Vector vmult(const Matrix& A, const Vector& x);

  /** Return a transpose matrix-vector product. */
  Vector Tvmult(const Matrix& A, const Vector& x);

  /** Insert a matrix into an output stream. */
  std::ostream& operator<<(std::ostream& os, const Matrix& A);

}

#endif //MATRIX_H

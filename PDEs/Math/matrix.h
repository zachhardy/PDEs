#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <sstream>

#include <cstddef>
#include <vector>


namespace pdes::Math
{

//########## Forward declarations
class Vector;

/**
 * Implementation of a general linear algebra dense matrix.
 */
class Matrix
{
public:
  using value_type = double;
  using STLMatrix = std::vector<std::vector<value_type>>;
  using InitializerMatrix =
      std::initializer_list<std::initializer_list<double>>;

protected:
  std::vector<Vector> coeffs;

public:

  //================================================== Constructors

  /**
   * Default constructor. */
  Matrix() = default;

  /**
   * Construct a matrix with \p n_rows and \p n_cols. */
  explicit
  Matrix(const size_t n_rows, const size_t n_cols);

  /**
   * Construct a Matrix with \p n_rows and \p n_cols set to \p value. */
  explicit
  Matrix(const size_t n_rows,
         const size_t n_cols,
         const value_type value);

  /**
   * Copy construction with nested STL vectors.
   */
  Matrix(const STLMatrix& other);

  /**
   * Move construction from nested STL vectors.
   */
  Matrix(STLMatrix&& other);

  /**
   * Construction from nested initializer lists.
   */
  Matrix(const InitializerMatrix list);

  //================================================== Assignment

  /**
   * Copy assignment from nested STL vectors.
   */
  Matrix&
  operator=(const STLMatrix& other);

  /**
   * Move assignment from nested STL vectors.
   */
  Matrix&
  operator=(STLMatrix&& other);

  /**
   * Assignment to a scalar value.
   */
  Matrix&
  operator=(const value_type value);

  //================================================== Comparison

  /**
   * Test the equality of two matrices.
   */
  bool
  operator==(const Matrix& other) const;

  /**
   * Test the inequality of two matrics.
   */
  bool
  operator!=(const Matrix& other) const;

  //================================================== Characteristics

  /** \name Characteristics */
  // @{

  /**
   * Return the number of rows in the matrix.
   */
  size_t
  n_rows() const;

  /**
   * Return the number of columns in the matrix.
   */
  size_t 
  n_cols() const;

  /**
   * Return the number of elements in the matrix.
   */
  size_t
  size() const;

  /**
   * Return the number of non-zero elements in the matrix.
   */
  size_t
  nnz() const;

  /**
   * Return whether the matrix is empty.
   */
  bool
  empty() const;

  /**
   * Return whether the matrix is uniformly zero.
   */
  bool
  all_zero() const;

  // @}

  //================================================== Iterators

  /** \name Iterators */
  // @{

  /**
   * Mutable iterator to the start of the matrix.
   */
  std::vector<Vector>::iterator
  begin();

  /**
   * Mutable iterator to the end of the matrix.
   */
   std::vector<Vector>::iterator
   end();

  /**
   * Constant iterator to the start of the matrix.
   */
  std::vector<Vector>::const_iterator
  begin() const;

  /**
   * Constant iterator to the end of the matrix.
   */
  std::vector<Vector>::const_iterator
  end() const;

  /**
   * Mutable iterator to the start of row \p i of the matrix.
   */
  std::vector<value_type>::iterator
  begin(const size_t i);

  /**
   * Mutable iterator to the end of row \p i of the matrix.
   */
  std::vector<double>::iterator
  end(const size_t i);

  /**
   * Constant iterator to the start of row \p i of the matrix.
   */
  std::vector<double>::const_iterator
  begin(const size_t i) const;

  /**
   * Constant iterator to the end of row \p i of the matrix.
   */
  std::vector<double>::const_iterator
  end(const size_t i) const;

  // @}

  //================================================== Accessors

  /** \name Accessors */
  // @{

  /**
   * Read and write access for row \p i of the matrix.
   */
  Vector&
  operator[](const size_t i);

  /**
   * Read access for row \p i of the matrix.
   */
  const Vector&
  operator[](const size_t i) const;

  /**
   * Read and write access for row \p i of the matrix.
   */
  Vector&
  operator()(const size_t i);

  /**
   * Read access for row \p i of the matrix.
   */
  const Vector&
  operator()(const size_t i) const;

  /**
   * Read and write access for row \p i of the matrix with bounds checking.
   */
  Vector&
  at(const size_t i);

  /**
   * Read access for row \p i of the matrix with bounds checking.
   */
  const Vector&
  at(const size_t i) const;

  /**
   * Read and write access for element <tt>(i, j)</tt>.
   */
  value_type&
  operator()(const size_t i, const size_t j);

  /**
   * Read and write access for element <tt>(i, j)</tt>.
   */
  const value_type&
  operator()(const size_t i, const size_t j) const;

  /**
   * Read and write access for element <tt>(i, j)</tt> with bounds checking.
   */
  value_type&
  at(const size_t i, const size_t j);

  /**
   * Read and write access for element <tt>(i, j)</tt> with bounds checking.
   */
  const value_type&
  at(const size_t i, const size_t j) const;

  /**
   * Read and write access to the <tt>i</tt>'th diagonal element.
   */
   value_type&
   diagonal(const size_t i);

  /**
   * Read access to the <tt>i</tt>'th diagonal element.
   */
  const value_type&
  diagonal(const size_t i) const;

  /**
   * Return the diagonal of the matrix.
   */
  Vector
  diagonal() const;

  /**
   * Return a mutable pointer to the underlying matrix data.
   */
  Vector*
  data();

  /**
   * Return a constant pointer to the underlying matrix data.
   */
  const Vector*
  data() const;

  /**
   * Return a pointer to the underlying data for row \p i of the matrix.
   */
  value_type*
  data(const size_t i);

  /**
   * Return a pointer to the underlying data for row \p i of the matrix.
   */
  const value_type*
  data(const size_t i) const;

  // @}

  //================================================== Modifiers

  /** \name Modifiers */
  // @{

  /**
   * Return the matrix to its uninitialized state.
   */
  void
  clear();

  /**
   * Remove the last row from the matrix.
   */
  void
  pop_back();

  /**
   * Add a row to the back of the matrix.
   */
  void
  push_back(const Vector& row);

  /**
   * Move a row to the back of the matrix.
   */
  void
  push_back(Vector&& row);

  /**
   * Resize the matrix to \p n_rows and \p n_cols.
   */
  void
  resize(const size_t n_rows, const size_t n_cols);

  /**
   * Resize to \p n_rows and \p n_cols, setting new elements to \p value.
   */
  void
  resize(const size_t n_rows,
         const size_t n_cols,
         const value_type value);

  /**
   * Alias to \ref resize.
   */
  void
  reinit(const size_t n_rows, const size_t n_cols);

  /**
   * Alias to \ref resize.
   */
  void
  reinit(const size_t n_rows,
         const size_t n_cols,
         const value_type value);

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
   * Swap the contents of this matrix with another.
   */
  void
  swap(Matrix& other);

  /**
   * Set the diagonal of the matrix with a vector.
   */
  void
  set_diagonal(const Vector& diag);

  /**
   * Set the diagonal of the matrix with a fixed scalar value.
   */
  void
  set_diagonal(const value_type value);

  // @}

  //================================================== Scalar Operations

  /** \name Scalar Operations */
  // @{



  /**
   * Negate the elements of the matrix. This is compute via
   * \f$ \boldsymbol{A} = -\boldsymbol{A} = -a_{ij}, ~ \forall i, j \f$.
   */
  Matrix&
  operator-();

  /**
   * Return a matrix with the negated elements.
   * \see Matrix::operator-()
   */
  Matrix
  operator-() const;

  /**
   * Multiply the elements of the matrix by a scalar value.
   * This is computed via
   * \f$ \boldsymbol{A} = \alpha \boldsymbol{A}
   *                    = \alpha a_{ij} ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator*=(const value_type value);

  /**
   * Divide the elements of the matrix by a scalar value.
   * This is computed via
   * \f$ \boldsymbol{A} = \frac{\boldsymbol{A}}{\alpha}
   *                    = \frac{a_{ij}}{\alpha}, ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator/=(const value_type value);

  //@}

  //================================================== Linear Algebra

  /** \name Linear Algebra */
  // @{

  /**
   * Addition of a scaled matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} + \alpha \boldsymbol{B}
   *                    = a_{ij} + \alpha b_{ij}, ~ \forall i, j
   * \f$.
   */
  void add(const Matrix& B, const value_type factor = 1.0);


  /**
   * Addition of a scaled transpose Matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} + \alpha \boldsymbol{B}^T
   *                    = a_{ij} + \alpha b{ji}, ~ \forall i, j
   * \f$.
   */
  void Tadd(const Matrix& B, const value_type factor = 1.0);

  /**
   * Add another matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} + \boldsymbol{B}
   *                    = a_{ij} + b_{ij}, ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator+=(const Matrix& B);

  /**
   * Subtract another matrix. This is computed via
   * \f$ \boldsymbol{A} = \boldsymbol{A} - \boldsymbol{B}
   *                    = a_{ij} - b_{ij}, ~ \forall i, j
   * \f$.
   */
  Matrix&
  operator-=(const Matrix& B);

  /**
   * Compute a matrix-matrix multiplication. This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}
   *                    = \sum_{k=0}^{n} a_{ik} b_{kj}, ~ \forall i, j
   * \f$.
   */
  void
  mmult(const Matrix& B, Matrix& C,
        const bool adding = false) const;

  /**
   * Return a matrix-matrix product.
   * \see Matrix::mmult
   */
  Matrix
  mmult(const Matrix& B) const;

  /**
   * Compute a transpose matrix-matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B} \\
   *                    = \sum_{k=0}^{n} a_{ki} b_{kj}, ~ \forall i, j
   * \f$.
   */
  void
  Tmmult(const Matrix& B, Matrix& C,
         const bool adding = false) const;

  /**
   * Return a transpose matrix-matrix product.
   * \see Matrix::Tmmult
   */
  Matrix
  Tmmult(const Matrix& B) const;

  /**
   * Compute a matrix-transpose matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A} \boldsymbol{B}^T \\
   *                    = \sum_{k=1}^{n} a_{ik} b_{jk}, ~ \forall i, j
   * \f$.
   */
  void
  mTmult(const Matrix& B, Matrix& C,
         const bool adding = false) const;

  /**
   * Return a matrix-transpose matrix multiplication.
   * \see Matrix::mTmult
   */
  Matrix
  mTmult(const Matrix& B) const;

  /**
   * Compute a transpose matrix-transpose matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T \\
   *                    = \sum_{k=1}^{n} a_{ki} b_{jk}, \forall i, j
   * \f$.
   */
  void
  TTmult(const Matrix& B, Matrix& C,
         const bool adding = false) const;

  /**
   * Return a transpose matrix-transpose matrix multiplication.
   * This is computed via
   * \f$ \boldsymbol{C} = \boldsymbol{A}^T \boldsymbol{B}^T
   *                    = \sum_{k=0}^{n} a_{ki} b_{jk}, ~ \forall i, j
   * \f$.
   */
  Matrix
  TTmult(const Matrix& B) const;

  /**
   * Compute a matrix-vector product.
   * This is computed via
   * \f$ \vec{y} = \boldsymbol{A} \vec{x} \\
   *             = \sum_{j=1}^{n} a_{ij} x_j, ~ \forall i
   * \f$.
   */
  void
  vmult(const Vector& x, Vector& y,
        const bool adding = false) const;

  /**
   * Return a matrix-vector product.
   * \see Matrix::vmult
   */
  Vector
  vmult(const Vector& x) const;

  /**
   * Add a matrix-vector product to the destination vector.
   * \see Matrix::vmult
   */
  void
  vmult_add(const Vector& x, Vector& y) const;

  /**
   * Compute a transpose matrix-vector product.
   * This is computed via
   * \f$ \vec{y} = \boldsymbol{A}^T \vec{x}
   *             = \sum_{i=1}^{n} a_{ji} x_i, ~ \forall i
   * \f$.
   */
  void
  Tvmult(const Vector& x, Vector& y,
         const bool adding = false) const;

  /**
   * Return a transpose matrix-vector product.
   * \see Matrix::Tvmult
   */
  Vector
  Tvmult(const Vector& x) const;

  /**
   * Add a transpose matrix-vector product to the destination vector.
   * \see Matrix::Tvmult
   */
  void
  Tvmult_add(const Vector& x, Vector& y);

  /**
   * Compute a matrix-vector product.
   * \see Matrix::vmult
   */
  Vector
  operator*(const Vector& x) const;

  /**
   * Return the transpose of the matrix.
   */
  Matrix
  transpose() const;

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
  void
  print(std::ostream& os = std::cout,
        const bool scientific = false,
        const unsigned int precision = 3,
        const unsigned int width = 0) const;

  /**
   * Return the matrix as a string.
   *
   * \param scientific A flag for scientific notation.
   * \param precision The precision of the digits to display.
   * \param width The spacing between entries.
   */
  std::string
  str(const bool scientific = false,
      const unsigned int precision = 3,
      const unsigned int width = 0) const;

  // @}

private:
  static bool
  valid_dimensions(const STLMatrix& A);

};

//================================================== Methods

/**
 * Add two matrices together.
 * \see Matrix::add Matrix::operator+=
 */
Matrix
operator+(const Matrix& A, const Matrix& B);

/**
 * Subtract two matrices.
 * \see Matrix::operator-=
 */
Matrix
operator-(const Matrix& A, const Matrix& B);

/**
 * Return a matrix-matrix multiplication.
 * \see Matrix::mmult
 */
Matrix
operator*(const Matrix& A, const Matrix& B);

/**
 * Compute a matrix-matrix multiplication.
 * \see Matrix::mmult
 */
void
mmult(const Matrix& A, const Matrix& B, Matrix& C);

/**
 * Return a matrix-matrix multiplication.
 * \see Matrix::mmult
 */
Matrix
mmult(const Matrix& A, const Matrix& B);

/**
 * Compute a transpose matrix-matrix multiplication.
 * \see Matrix::Tmmult
 */
void
Tmmult(const Matrix& A, const Matrix& B, Matrix& C);

/**
 * Return a transpose matrix-matrix multiplication.
 * \see Matrix::Tmmult
 */
Matrix
Tmmult(const Matrix& A, const Matrix& B);

/**
 * Compute a matrix-transpose matrix multiplication.
 * \see Matrix::mTmult
 */
void
mTmult(const Matrix& A, const Matrix& B, Matrix& C);

/**
 * Return a matrix-transpose matrix multiplication.
 * \see Matrix::mTmult
 */
Matrix
mTmult(const Matrix& A, const Matrix& B);

/**
 * Compute a transpose matrix-transpose matrix multiplication.
 * \see Matrix::TTmmult
 */
void
TTmult(const Matrix& A, const Matrix& B, Matrix& C);

/**
 * Return a transpose matrix-transpose matrix multiplication.
 * \see Matrix::TTmmult
 */
Matrix
TTmult(const Matrix& A, const Matrix& B);

/**
 * Compute a matrix-vector product.
 * \see Matrix::vmult
 */
void
vmult(const Matrix& A, const Vector& x, Vector& y);

/**
 * Return a matrix-vector product.
 * \see Matrix::vmult
 */
Vector
vmult(const Matrix& A, const Vector& x);

/**
 * Compute a transpose matrix-vector product.
 * \see Matrix::vmult
 */
void
Tvmult(const Matrix& A, const Vector& x, Vector& y);

/**
 * Return a transpose matrix-vector product.
 * \see Matrix::vmult
 */
Vector
Tvmult(const Matrix& A, const Vector& x);

/**
 * Insert a matrix into an output stream.
 */
std::ostream&
operator<<(std::ostream& os, const Matrix& A);

}

#endif //MATRIX_H

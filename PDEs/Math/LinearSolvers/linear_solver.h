#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "matrix.h"
#include "sparse_matrix.h"
#include "vector.h"

namespace math
{

enum class LinearSolverType
{
  LU = 0,
  CHOLESKY = 1,
};

//######################################################################
/**
 * An abstract base class for solving the linear system
 * \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
template<typename value_type>
class LinearSolver
{
protected:
  bool initialized = false;
  Matrix<value_type>& A;

public:
  LinearSolver(Matrix<value_type>& matrix) : A(matrix)
  {
    if (A.n_rows() != A.n_cols())
    {
      std::stringstream err;
      err << "LinearSystem::" << __FUNCTION__ << ": "
          << "Invalid inputs. The matrix must be square and of the same "
          << "dimension as the right-hand side vector.";
      throw std::runtime_error(err.str());
    }
  }

  void set_matrix(Matrix<value_type>& matrix)
  {
    A = matrix;
    initialized = false;
  }

public:
  /** Abstract setup method. */
  virtual void setup() = 0;

  /** Abstract solve method. */
  virtual Vector<value_type> solve(const Vector<value_type>& b) = 0;
};


//######################################################################
/**
 * An abstract base class for solving the linear system
 * \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
template<typename value_type>
class SparseLinearSolver
{
protected:
  bool initialized = false;
  SparseMatrix<value_type>& A;

public:
  SparseLinearSolver(SparseMatrix<value_type>& matrix) : A(matrix)
  {
    if (A.n_rows() != A.n_cols())
    {
      std::stringstream err;
      err << "LinearSystem::" << __FUNCTION__ << ": "
          << "Invalid inputs. The matrix must be square and of the same "
          << "dimension as the right-hand side vector.";
      throw std::runtime_error(err.str());
    }
  }

  void set_matrix(SparseMatrix<value_type>& matrix)
  {
    A = matrix;
    initialized = false;
  }

public:
  /** Abstract setup method. */
  virtual void setup() = 0;

  /** Abstract solve method. */
  virtual Vector<value_type> solve(const Vector<value_type>& b) = 0;
};




}
#endif //LINEAR_SOLVER_H

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
template<typename number>
class LinearSolver
{
protected:
  bool initialized = false;
  Matrix<number>& A;

public:
  LinearSolver(Matrix<number>& matrix) : A(matrix)
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

  void set_matrix(Matrix<number>& matrix)
  {
    A = matrix;
    initialized = false;
  }

public:
  /** Abstract setup method. */
  virtual void setup() = 0;

  /** Abstract solve method. */
  virtual Vector<number> solve(const Vector<number>& b) = 0;
};


//######################################################################
/**
 * An abstract base class for solving the linear system
 * \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
 */
template<typename number>
class SparseLinearSolver
{
protected:
  bool initialized = false;
  SparseMatrix<number>& A;

public:
  SparseLinearSolver(SparseMatrix<number>& matrix) : A(matrix)
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

  void set_matrix(SparseMatrix<number>& matrix)
  {
    A = matrix;
    initialized = false;
  }

public:
  /** Abstract setup method. */
  virtual void setup() = 0;

  /** Abstract solve method. */
  virtual Vector<number> solve(const Vector<number>& b) = 0;
};




}
#endif //LINEAR_SOLVER_H

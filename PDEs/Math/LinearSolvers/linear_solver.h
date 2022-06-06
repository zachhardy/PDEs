#ifndef LINEAR_SOLVER_BASE_H
#define LINEAR_SOLVER_BASE_H

#include <cstddef>
#include <string>

//########## Forward declarations
namespace Math
{
  class Vector;
  class Matrix;
  class SparseMatrix;
}

namespace Math::LinearSolver
{

  struct Options
  {
    size_t verbosity = 0;

    double tolerance = 1.0e-6;
    size_t max_iterations = 500;
  };


  /**
   * Base class from which all linear solvers must derive.
   */
  template<class MatrixType>
  class LinearSolverBase
  {
  public:
    /**
     * Abstract method for solving a linear system.
     */
    virtual void
    solve(Vector& x, const Vector& b) const = 0;

    /**
     * Return the solution to \f$ \boldsymbol{A} \vec{x} = \vec{b} \f$.
     */
    Vector
    solve(const Vector& b) const;

    /**
     * Attach a matrix to the solver.
     */
    virtual void
    set_matrix(const MatrixType& matrix);
  };


  //###########################################################################


  /**
   * Base class for direct solvers.
   */
  template<class MatrixType>
  class DirectSolverBase : public LinearSolverBase<MatrixType>
  {
  protected:
    MatrixType A;
    bool factorized = false;

  public:
    /**
     * Default constructor.
     */
    DirectSolverBase();

    /**
     * Abstract routine for factorizing the matrix.
     */
    virtual void
    factorize() = 0;

    /**
     * Attach a matrix to the solver.
     */
    virtual void
    set_matrix(const MatrixType& matrix) override;
  };


  //###########################################################################


  /**
   * Base class for iterative solvers.
   */
  class IterativeSolverBase : public LinearSolverBase<SparseMatrix>
  {
  protected:
    const std::string solver_name;
    size_t verbosity = 0;

    const SparseMatrix* A;

    double tolerance;
    size_t max_iterations;

  public:
    /**
     * Default constructor using Options struct to set parameters.
     *
     * \note The Options struct contains all parameters necessary for all
     *       implemented iterative solvers. Each derived class should set
     *       any and all appropriate parameters from this.
     */
    IterativeSolverBase(const Options& opts = Options(),
                        const std::string name = "Undefined");

    /**
     * Attach a matrix to the solver.
     */
    void set_matrix(const SparseMatrix& matrix) override;

  protected:
    /**
     * Check whether the solver has converged.
     */
    virtual bool
    check(const size_t iteration, const double value) const;

    /**
     * Throw an error when convergence criteria is not met.
     */
    void
    throw_convergence_error(const size_t iteration,
                            const double value) const;
  };

}
#endif //LINEAR_SOLVER_BASE_H

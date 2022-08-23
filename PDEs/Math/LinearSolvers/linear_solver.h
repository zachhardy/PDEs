#ifndef LINEAR_SOLVER_BASE_H
#define LINEAR_SOLVER_BASE_H

#include <cstddef>
#include <string>


namespace PDEs
{
  namespace Math
  {
    // forward declarations
    class Vector;
    class Matrix;
    class SparseMatrix;


    namespace LinearSolvers
    {

      /**
       * A struct containing solver options for linear solvers. This is used to
       * minimize the number of arguments in iterative solver constructors.
       */
      struct Options
      {
        /**
         * The convergence tolerance.
         */
        double tolerance = 1.0e-6;

        /**
         * The maximum number of iterations allowed.
         */
        unsigned int max_iterations = 500;

        /**
         * The level of screen output for the linear solver.
         */
        unsigned int verbosity = 0;

        /**
         * Default constructor.
         */
        Options(const double tolerance = 1.0e-6,
                const unsigned int max_iterations = 500,
                const unsigned int verbosity = 0);
      };


      /**
       * A base class from which all linear solvers must derive. This is templated
       * on the MatrixType in order to accommodate both dense matrices and sparse
       * matrices.
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
         * Return the solution to \f$ A x = b \f$.
         */
        Vector
        solve(const Vector& b) const;

        /**
         * Attach a matrix to the solver.
         */
        virtual void
        set_matrix(const MatrixType& matrix);
      };


      /**
       * Base class for direct solvers.
       */
      template<class MatrixType>
      class DirectSolverBase : public LinearSolverBase<MatrixType>
      {
      public:
        using LinearSolverBase<MatrixType>::solve;

        /**
         * Default constructor.
         */
        DirectSolverBase();

        /**
         * An abstract routine for applying a particular factorization to the
         * attached matrix.
         */
        virtual void
        factorize() = 0;

        /**
         * Attach a matrix to the direct solver.
         */
        virtual void
        set_matrix(const MatrixType& matrix) override;

      protected:
        /**
         * The matrix that describes the linear system being solved. This matrix
         * will be factorized, and therefore modified.
         */
        MatrixType A;

        /**
         * A flag specifying whether the matrix is in its factorized form or still
         * requires factorization before \ref solve can be called.
         */
        bool factorized = false;
      };


      /**
       * Base class for iterative solvers. This is defaulted to sparse matrices.
       */
      class IterativeSolverBase : public LinearSolverBase<SparseMatrix>
      {
      public:
        using LinearSolverBase<SparseMatrix>::solve;

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
         * Attach the sparse matrix to the iterative linear solver.
         */
        void
        set_matrix(const SparseMatrix& matrix) override;


      protected:
        /**
         * Check whether the solver has converged. This is an abstract class that
         * should be overridden with the appropriate check for convergence for the
         * particular iterative solver.
         */
        virtual bool
        check(const unsigned int iteration, const double value) const;

        /**
         * Throw an error when convergence criteria is not met.
         */
        void
        throw_convergence_error(const unsigned int iteration,
                                const double value) const;

      protected:
        /**
         * A constant pointer to the sparse matrix of the linear system.
         */
        const SparseMatrix* A;

        /**
         * The iteration convergence tolerance.
         */
        double tolerance;

        /**
         * The maximum number of iterations allowed.
         */
        unsigned int max_iterations;

        /**
         * The level of screen output from the linear solver.
         */
        unsigned int verbosity = 0;

        /**
         * A string that is used to describe the iterative solver when outputting
         * results. This is set by the derived classes.
         */
        const std::string solver_name;
      };

    }
  }
}
#endif //LINEAR_SOLVER_BASE_H

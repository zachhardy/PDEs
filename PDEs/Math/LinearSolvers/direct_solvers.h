#ifndef DIRECT_SOLVERS_H
#define DIRECT_SOLVERS_H

#include "linear_solver.h"

#include "matrix.h"
#include "Math/sparse_matrix.h"

#include <vector>
#include <cstddef>


namespace PDEs
{
  namespace Math
  {
    namespace LinearSolvers
    {

      /**
       * Implementation of an LU decomposition solver. See more at
       * <a href="https://en.wikipedia.org/wiki/LU_decomposition">
       * https://en.wikipedia.org/wiki/LU_decomposition</a>.
       */
      class LU : public DirectSolverBase<Matrix>
      {
      public:
        /**
         * Default constructor. Construct an direct LU solver, optionally with
         * row pivoting.
         */
        LU(const bool pivot = true);

        /**
         * Attach a dense matrix to the solver.
         */
        void
        set_matrix(const Matrix& matrix) override;

        /**
         * Perform an LU factorization of the matrix \f$ A \f$ in
         * place.
         *
         * An LU factorization defines the relationship
         * \f$ A = L U \f$ where  \f$ L \f$ is a lower triangular matrix and
         * \f$ U \f$ is an upper triangular matrix. The factorization is
         * performed in place.
         *
         * The algorithm used to do perform this factorization is an extension
         * of the formation of a row-echelon form matrix in that the upper
         * triangular matrix is identical to the row-echelon form. The lower
         * triangular matrix then contains the row operations used to form upper
         * triangular system.
         */
        void
        factorize() override;

        /**
         * Solve an LU factored linear system.
         *
         * The LU factored linear system is define by \f$ A x = L U x =b \f$.
         * The solution \f$ x \f$ is obtained in a two step process. First,
         * define \f$ y = U x \f$ and plug this in to obtain \f$ L y = b \f$.
         * The vector \f$ y \f$ can be obtained using forward-substitution after
         * reordering the right-hand side vector \f$ b \f$ according to the
         * pivot mapping vector. Next, the solution \f$ x \f$ is computed using
         * the previous definition  \f$ U x = y \f$ where \f$ y \f$ is now the
         * source term. This system can be solved using back-substitution.
         */
        void
        solve(Vector& x, const Vector& b) const override;

      private:
        /**
         * A flag designating whether or not to perform row pivots.
         */
        bool pivot_flag;

        /**
         * The pivot mapping vector. The index corresponds to the initial row
         * number and the value to the pivoted row number. This is used to map the
         * right-hand side vector to the correct row when solving.
         */
        std::vector<size_t> row_pivots;
      };


      /**
       * Implementation of a Cholesky decomposition solver. See more at
       * <a href="https://en.wikipedia.org/wiki/Cholesky_decomposition">
       * https://en.wikipedia.org/wiki/Cholesky_decomposition</a>.
       */
      class Cholesky : public DirectSolverBase<Matrix>
      {
      public:
        /**
         * Default constructor.
         */
        Cholesky();

        /**
         * Perform a Cholesky factorization on the matrix \f$ A \f$.
         *
         * Cholesky factorization is meant for symmetric positive definite
         * matrices into a lower triangular matrix and its transpose. This method
         * is more efficient than the LU decomposition when applicable.
         *
         * \note Checks are not performed to ensure symmetric positive
         *       definiteness.  The user is responsible for ensuring the matrix
         *       fits this criteria.
         */
        void
        factorize() override;

        /**
         * Solve the Cholesky factored linear system.
         *
         * The Cholesky solve is a specialization of the LU solve in that
         * \f$ U = L^T \f$.
         *
         * \see LU::solve
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };


      /**
       * Implementation of a sparse LU solver. See more at
       * <a href="https://en.wikipedia.org/wiki/LU_decomposition">
       * https://en.wikipedia.org/wiki/LU_decomposition</a>.
       */
      class SparseLU : public DirectSolverBase<SparseMatrix>
      {
      public:
        /**
         * Default constructor. Construct a sparse LU solver, optionally with
         * row pivoting.
         */
        SparseLU(const bool pivot = true);

        /**
         * Attach a matrix to the solver.
         */
        void
        set_matrix(const SparseMatrix& matrix) override;

        /**
         * Perform an LU factorization on the matrix \f$ A \f$ in
         * place. See \ref LU::factorize
         */
        void
        factorize() override;

        /** Solve the LU factored linear system. See \ref LU::solve */
        void
        solve(Vector& x, const Vector& b) const override;

      private:
        /**
         * A flag designating whether or not to perform row pivots.
         */
        bool pivot_flag = true;

        /**
         * The pivot mapping vector. The index corresponds to the initial row
         * number and the value to the pivoted row number. This is used to map the
         * right-hand side vector to the correct row when solving.
         */
        std::vector<size_t> row_pivots;
      };


      /**
       * Implementation of a sparse Cholesky solver. See more at
       * <a href="https://en.wikipedia.org/wiki/Cholesky_decomposition">
       * https://en.wikipedia.org/wiki/Cholesky_decomposition</a>.
       */
      class SparseCholesky : public DirectSolverBase<SparseMatrix>
      {
      public:
        SparseCholesky();

        /**
         * Perform a Cholesky factorization on the matrix \f$ A \f$.
         * See \ref Cholesky::solve
         */
        void
        factorize() override;

        /**
         * Solve the Cholesky factored linear system. See \ref Cholesky::solve
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };

    }
  }
}
#endif //DIRECT_SOLVERS_H

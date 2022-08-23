#ifndef LU_H
#define LU_H

#include "../linear_solver.h"

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
       * Implementation of an LU decomposition solver.
       *
       * An LU solver is a direct solver for general matrices. This method
       * first performs a factorization of the matrix \f$ A \f$ into an upper
       * and lower triangular matrix such that \f$ A = L U \f$. This is often
       * performed in place. The factorization algorithm follows that of
       * the Gaussian elimination algorithm. In fact, upper triangular matrix of
       * an LU factorization is the row-echelon form matrix obtained from
       * Gaussian elimination. The lower triangular matrix contains the row
       * operations used to form upper triangular system.
       *
       * The LU factored linear system is defined by \f$ A x = L U x = b \f$.
       * The solution \f$ x \f$ is obtained in a two step process. First,
       * define \f$ y = U x \f$ and plug this in to obtain \f$ L y = b \f$.
       * The vector \f$ y \f$ can be obtained using forward-substitution after
       * reordering the right-hand side vector \f$ b \f$ according to the
       * #pivot_mapping vector. Next, the solution \f$ x \f$ is computed using
       * the previous definition  \f$ U x = y \f$ where \f$ y \f$ is now the
       * source term. This system can be solved using back-substitution.
       *
       * Optionally, row-pivoting can be used. This functionality improves the
       * numerical stability of the algorithm. As the factorizations sweep
       * along the columns, division by the diagonal is performed at each
       * subsequent step. When pivoting is on, rows will be swapped when there
       * is a sub-diagonal entry of the current column that is larger than the
       * current diagonal. By doing this, division by the largest available
       * number is always performed.
       *
       * See more at https://en.wikipedia.org/wiki/LU_decomposition.
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
         * Perform an LU factorization of the matrix \f$ A \f$ in place.
         */
        void
        factorize() override;

        /**
         * Solve an LU factored linear system.
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
         * number and the value to the pivoted row number. This is used to map
         * the right-hand side vector to the correct row when solving.
         */
        std::vector<size_t> row_pivots;
      };


      /**
       * Implementation of a sparse LU solver. For descriptions of the LU
       * decomposition solver see \ref LU.
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

    }
  }
}


#endif //LU_H

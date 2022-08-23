#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "Math/LinearSolvers/linear_solver.h"

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
       * Implementation of a Cholesky decomposition solver.
       *
       * Cholesky decompositions are only valid for symmetric positive definite
       * matrices. Similar to the LU decomposition, this method first performs
       * a factorization of \f$ A \f$ such that \f$ A = L L^T \f$. This special
       * case of an LU decomposition should be used when applicable as that it
       * is more efficient than its more general LU decomposition counterpart.
       *
       * The Cholesky factored system can be solved in the same way as the LU
       * decomposition noting that \f$ U = L^T \f$. See \ref LU for more
       * details.
       *
       * See more at https://en.wikipedia.org/wiki/Cholesky_decomposition.
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
         * \note Checks are not performed to ensure symmetric positive
         *       definiteness.  The user is responsible for ensuring the matrix
         *       fits this criteria.
         */
        void
        factorize() override;

        /**
         * Solve the Cholesky factored linear system. See \ref LU::solve
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };


      /**
       * Implementation of a sparse Cholesky solver. For descriptions of the
       * Cholesky decomposition solver see \erf Cholesky.
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
#endif //CHOLESKY_H

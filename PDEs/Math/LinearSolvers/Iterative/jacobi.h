#ifndef JACOBI_H
#define JACOBI_H

#include "LinearSolvers/linear_solver.h"


namespace PDEs
{
  namespace Math
  {
    namespace LinearSolvers
    {

      /**
       * Implementation of the Jacobi iterative method.
       *
       * This is a splitting method which decomposes the matrix into diagonal
       * and off-diagonal components such that \f$ A = D + L + U \f$. The
       * resulting linear system is then \f$ A x = (D + L + U) x = b \f$.
       * The fixed-point problem is then formed by lagging the off-diagonal
       * terms:
       * \f[ x^{\ell + 1} = D^{-1} \left( b - (L + U) x^\ell \right). \f]
       *
       * See more at https://en.wikipedia.org/wiki/Jacobi_method.
       */
      class Jacobi : public IterativeSolverBase
      {
      public:
        /**
         * Default constructor.
         */
        Jacobi(const Options& opts = Options());

        /**
         * Iteratively solve the system using the Jacobi method.
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };

    }
  }
}

#endif //JACOBI_H

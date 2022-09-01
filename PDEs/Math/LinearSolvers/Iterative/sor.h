#ifndef SOR_H
#define SOR_H

#include "LinearSolvers/linear_solver.h"


namespace PDEs
{
  namespace Math
  {
    namespace LinearSolvers
    {
      /**
       * Implementation of the successive over-relaxation (SOR) iterative
       * method.
       *
       * Similar to the Jacobi method, SOR decomposes the matrix such that \f$
       * A = D + L + U \f$. The resulting linear system is multiplied by the
       * relaxation parameter \f$ \omega \f$ and rearranged to obtain
       * \f[
       *    (D + \omega L) x = \left( \omega b + \left[ (1 - \omega) D - U
       *    \right] x \right).
       * \f]
       * Using forward-substitution, the fixed point problem becomes
       * \f[
       *    x^{\ell+1} = (1 - \omega) x^\ell + \omega D^{-1} \left( b - L
       *    x^{\ell+1} - U x^{\ell} \right).
       * \f]
       *
       * See more at https://en.wikipedia.org/wiki/Successive_over-relaxation.
       */
      class SOR : public IterativeSolverBase
      {
      protected:
        const double omega; ///< The relaxation parameter

      public:
        /** Default constructor. */
        SOR(const double omega = 1.5,
            const Options& opts = Options(),
            const std::string solver_name = "SOR");

        /** Solve the system using the SOR iterative method. */
        virtual void solve(Vector& x, const Vector& b) const override;
      };


      /**
       * Implementation of the Gauss-Seidel iterative method. Gauss-Seidel is
       * a specialization of the SOR iterative method where \f$ \omega = 1.0
       * \f$.
       *
       * See more at https://en.wikipedia.org/wiki/Gauss-Seidel_method.
       */
      class GaussSeidel : public SOR
      {
      public:
        /** Default constructor. */
        GaussSeidel(const Options& opts = Options());
      };


      /**
       * Implementation of the symmetric successive over-relaxation (SSOR)
       * iterative method. This is a specialization of the SOR method for
       * symmetric matrices which performs a forward and backward SOR sweep
       * per iteration.
       */
      class SSOR : public SOR
      {
      public:
        /** Default constructor. */
        SSOR(const double omega = 1.5, const Options& opts = Options());

        /** Solve the system using the SSOR method. */
        void
        solve(Vector& x, const Vector& b) const override;
      };
    }
  }
}
#endif //SOR_H

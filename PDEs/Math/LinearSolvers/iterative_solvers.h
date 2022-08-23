#ifndef ITERATIVE_SOLVERS_H
#define ITERATIVE_SOLVERS_H

#include "LinearSolvers/linear_solver.h"


namespace PDEs
{
  namespace Math
  {
    namespace LinearSolvers
    {
      /**
       * Implementation of the Jacobi iterative method. See more at
       * <a href="https://en.wikipedia.org/wiki/Jacobi_method">
       * https://en.wikipedia.org/wiki/Jacobi_method</a>.
       */
      class Jacobi : public IterativeSolverBase
      {
      public:
        /**
         * Default constructor.
         */
        Jacobi(const Options& opts = Options());

        /**
         * Solve the system using the Jacobi iterative method. The method splits
         * the matrix into diagonal, strictly lower triangular, and strictly upper
         * triangular matrices such that \f$ A = D + L + U \f$. The off-diagonal
         * terms are lagged, leading to the update function \f$ x^{\ell+1} =
         * D^{-1} \left( b - (L + U) x^\ell \right) \f$.
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };


      /**
       * Implementation of the successive over-relaxation (SOR) iterative method.
       * See more at
       * <a href="https://en.wikipedia.org/wiki/Successive_over-relaxation">
       * https://en.wikipedia.org/wiki/Successive_over-relaxation</a>.
       */
      class SOR : public IterativeSolverBase
      {
      public:
        /**
         * Default constructor.
         *
         * \param omega The relaxation parameter.
         */
        SOR(const double omega = 1.5,
            const Options& opts = Options(),
            const std::string solver_name = "SOR");

        /**
         * Solve the system using the SOR iterative method. This method splits
         * the matrix into diagonal, strictly lower triangular, and strictly
         * upper triangular matrices such that \f$ A = D + L + U \f$. The system
         * \f$ A x = (D + L + U) x = b \f$ is then multiplied by \f$ \omega \f$
         * to obtain \f$ \left( \omega (D + L + U) - D + D \right) x = \omega b
         * \f$. Lagging the upper triangular terms yields the update function \f$
         * x^{\ell + 1} = (1 - \omega) x^\ell + \omega D^{-1} \left( b - L
         * x^{\ell + 1} - U x^\ell \right) \f$. This can be solved using
         * forward-substitution.
         */
        virtual void
        solve(Vector& x, const Vector& b) const override;

      protected:
        const double omega;
      };


      /**
       * Implementation of the Gauss-Seidel iterative method. Gauss-Seidel is
       * a specialization of the SOR iterative method where \f$ \omega = 1.0 \f$.
       * See <a href="https://en.wikipedia.org/wiki/Gauss-Seidel_method">
       * https://en.wikipedia.org/wiki/Gauss-Seidel_method</a>.
       */
      class GaussSeidel : public SOR
      {
      public:
        /**
         * Default constructor.
         */
        GaussSeidel(const Options& opts = Options());
      };


      /**
       * Implementation of the symmetric successive over-relaxation (SSOR)
       * iterative method. This is a specialization of the SOR method for
       * symmetric matrices.
       */
      class SSOR : public SOR
      {
      public:
        /**
         * Default constructor.
         * \param omega  The relaxation parameter.
         */
        SSOR(const double omega = 1.5, const Options& opts = Options());

        /**
         * Solve the system using the SSOR method. This performs a forward and
         * backward sweep using the SOR method.
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };


      /**
       * Implementation of the conjugate gradient (CG) method. This method is
       * only valid for symmetric positive definite matrices. See more at
       * <a href="https://en.wikipedia.org/wiki/Conjugate_gradient_method">
         * https://en.wikipedia.org/wiki/Conjugate_gradient_method</a>.
       */
      class CG : public IterativeSolverBase
      {
      public:
        /**
         * Default constructor.
         */
        CG(const Options& opts = Options());

        /**
         * Solve the system using the CG method. This method is similar to
         * gradient descent, however, an additional constraint is imposed that
         * requires all subsequent search directions to be orthogonal to those
         * prior.
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };
    }
  }
}
#endif //ITERATIVE_SOLVERS_H

#ifndef CG_H
#define CG_H

#include "LinearSolvers/linear_solver.h"


namespace PDEs
{
  namespace Math
  {
    namespace LinearSolvers
    {

      /**
       * Implementation of the conjugate gradient (CG) method.
       *
       * This method is similar to gradient descent, however, an addition
       * constraint is imposed which requires that each subsequent search
       * direction be orthogonal to those prior.
       *
       * This method is only applicable for symmetric positive definite
       * matrices.
       *
       * See more at https://en.wikipedia.org/wiki/Conjugate_gradient_method.
       */
      class CG : public IterativeSolverBase
      {
      public:
        /**
         * Default constructor.
         */
        CG(const Options& opts = Options());

        /**
         * Solve the system using the CG method.
         */
        void
        solve(Vector& x, const Vector& b) const override;
      };
    }
  }
}
#endif //CG_H

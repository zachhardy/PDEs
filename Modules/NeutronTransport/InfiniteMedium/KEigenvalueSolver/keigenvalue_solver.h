#ifndef KEIGENVALUE_SOLVER_H
#define KEIGENVALUE_SOLVER_H

#include "../SteadyStateSolver/steadystate_solver.h"


namespace NeutronTransport
{
  namespace InfiniteMedium
  {

    /**
     * Implementation of a \f$ k \f$-eigenvalue solver.
     */
    class KEigenvalueSolver : public SteadyStateSolver
    {
    public:
      /**
       * The convergence tolerance for the outer iterations.
       */
      double outer_tolerance = 1.0e-8;

      /**
       * The maximum number of outer iterations allowed.
       */
      unsigned int max_outer_iterations = 1000;


      using SteadyStateSolver::initialize;

      /**
       * Execute the multi-group diffusion \f$ k \f$-eigenvalue solver.
       */
      virtual void
      execute() override;

    protected:
      /**
       * Compute the total neutron production rate.
       */
      double
      compute_production();

      /**
       * The current estimate of the k-eigenvalue.
       */
      double k_eff = 1.0;
    };

  }
}


#endif //KEIGENVALUE_SOLVER_H

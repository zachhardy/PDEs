#ifndef KEIGENVALUE_SOLVER_H
#define KEIGENVALUE_SOLVER_H

#include "InfiniteMedium/SteadyStateSolver/steadystate_solver.h"


namespace InfiniteMedium
{
  /**
   * Implementation of a \f$ k \f$-eigenvalue solver.
   */
  class KEigenvalueSolver : public SteadyStateSolver
  {
  public:
    double outer_tolerance = 1.0e-8;
    unsigned int max_outer_iterations = 1000;

  protected:
    double k_eff = 1.0;

  public:
    using SteadyStateSolver::initialize;
    virtual void execute() override;

  protected:
    /** Compute the total neutron production rate. */
    double compute_production();
  };
}


#endif //KEIGENVALUE_SOLVER_H

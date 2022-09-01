#ifndef KEIGENVALUE_SOLVER_H
#define KEIGENVALUE_SOLVER_H

#include "../SteadyStateSolver/steadystate_solver.h"


namespace NeutronDiffusion
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
    /** The current estimate of the \f$ k \f$-eigenvalue. */
    double k_eff = 1.0;

  public:
    virtual void execute() override;

    /**
     * Write the result of the simulation to an output file with the
     * specified prefix. This creates a file named <tt><file_prefix>.data</tt>
     * in the \p directory and file named \f$ k_eff.txt \f$ with the
     * converged eigenvalue.
     */
    virtual void
    write(const std::string directory,
          const std::string file_prefix) const override;

  protected:

    /** Implementation of the power method. */
    void power_method();

    /** Compute the total neutron production rate. */
    double compute_production();
  };

}


#endif //KEIGENVALUE_SOLVER_H

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
    /**
     * The current estimate of the k-eigenvalue.
     */
    double k_eff = 1.0;

    /**
     * The convergence tolerance for the outer iterations.
     */
    double outer_tolerance = 1.0e-8;

    /**
     * The maximum number of outer iterations allowed.
     */
    unsigned int max_outer_iterations = 1000;

    /**
     * Execute the multi-group diffusion \f$ k \f$-eigenvalue solver.
     */
    virtual void
    execute() override;

    /**
     * Write the result of the simulation to an output file. This writes the
     * same file as that in the SteadyStateSolver, but also a file named
     * \p k_eff.txt which contains the converged eigenvalue.
     *
     * \param output_directory The directory where the output should be placed.
     * \param file_prefix The name of the file without a suffix. By default,
     *      the suffix \p .data will be added to this input.
     */
    virtual void
    write(const std::string& output_directory,
          const std::string& file_prefix) const override;

  protected:
    /**
     * Implementation of the power method algorithm.
     */
    void
    power_method();

    /**
     * Compute the total neutron production rate.
     */
    double
    compute_production();
  };

}


#endif //KEIGENVALUE_SOLVER_H

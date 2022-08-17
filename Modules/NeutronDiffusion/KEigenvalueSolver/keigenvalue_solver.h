#ifndef KEIGENVALUE_SOLVER_H
#define KEIGENVALUE_SOLVER_H

#include "../SteadyStateSolver/steadystate_solver.h"


namespace NeutronDiffusion
{

  /** Implementation of a k-eigenvalue solver. */
  class KEigenvalueSolver : public SteadyStateSolver
  {
  public:
    double k_eff = 1.0;

    unsigned int max_outer_iterations = 1000;
    double outer_tolerance = 1.0e-8;

    /*-------------------- Interface Routines ----------*/
  public:
    virtual void execute() override;

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

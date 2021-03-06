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

    double tolerance = 1.0e-8;
    size_t max_iterations = 1000;

  public:
    /** Execute the k-eigenvalue solver. */
    virtual void execute() override;

  protected:
    /** Implementation of the power method algorithm. */
    void power_method();

    /** Compute the total neutron production rate. */
    double compute_production();

  public:
    virtual void write(const std::string& output_directory,
                       const std::string& file_prefix) const override;
  };

}


#endif //KEIGENVALUE_SOLVER_H

#include "steadystate_solver.h"

/// Run the steady state multigroup diffusion simulation.
void neutron_diffusion::SteadyStateSolver::execute()
{
  assemble_matrix();
  assemble_rhs();

  linear_solver->setup();
  phi = linear_solver->solve(system_rhs);
}

//######################################################################

/** \brief Assemble the matrix with the routine associated with the
 *  discretization type. */
void neutron_diffusion::SteadyStateSolver::assemble_matrix()
{
  switch (discretization->type)
  {
    case DiscretizationMethod::FINITE_VOLUME: assemble_fv_matrix();
    default: return;
  }
}

/** \brief Assemble the right-hand side vector with the routine associated with
 *  the discretization type. */
void neutron_diffusion::SteadyStateSolver::assemble_rhs()
{
  switch (discretization->type)
  {
    case DiscretizationMethod::FINITE_VOLUME: assemble_fv_rhs();
    default: return;
  }
}

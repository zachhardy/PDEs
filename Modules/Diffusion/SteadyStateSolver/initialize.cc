#include "steadystate_solver.h"


/**
 * \brief Initialize the solver.
 *
 * Check the inputs of the solver to ensure that they exist and are
 * appropriately specified. Perform basic operations to parse the inputs and
 * set up to solver.
 */
void diffusion::SteadyStateSolver::initialize()
{
  std::cout << "Initializing the solver." << std::endl;

  check_inputs();

  initialize_materials();

  initialize_boundaries();

  phi.resize(discretization->n_dofs(n_groups));
  if (use_precursors)
    precursors.resize(discretization->n_dofs(max_precursors_per_material));

  system_rhs.resize(discretization->n_dofs(n_groups));
  system_matrix.resize(discretization->n_dofs(n_groups),
                       discretization->n_dofs(n_groups));
}


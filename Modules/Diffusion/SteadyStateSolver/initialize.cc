#include "steadystate_solver.h"


/**
 * Check the inputs of the solver to ensure that they exist and are
 * appropriately spcified. Perform basic operations to parse the inputs and
 * set up to solver.
 */
void diffusion::SteadyStateSolver::initialize()
{
  std::cout << "Initializing the solver." << std::endl;

  // Preliminary input checks
  check_inputs();

  // Initialize the materials
  initialize_materials();

  // Initialize the boundaries
  initialize_boundaries();

  // Initialize data storage
  phi.resize(discretization->n_dofs(n_groups));
  if (use_precursors)
    precursors.resize(discretization->n_dofs(max_precursors_per_material));

  system_rhs.resize(discretization->n_dofs(n_groups));
  system_matrix.resize(discretization->n_dofs(n_groups),
                       discretization->n_dofs(n_groups));
}


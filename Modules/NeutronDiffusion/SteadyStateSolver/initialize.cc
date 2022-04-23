#include "steadystate_solver.h"

/// Initialize the multigroup diffusion solver.
void neutron_diffusion::SteadyStateSolver::initialize()
{
  std::cout << "Initializing the solver...\n";

  input_checks();
  initialize_materials();
  initialize_boundaries();

}


/// Validate the general setup of the simulation.
void neutron_diffusion::SteadyStateSolver::input_checks()
{
  // Check that the mesh
  if (mesh == nullptr)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "No mesh attached to the solver.";
    throw std::runtime_error(err.str());
  }
  if (mesh->dim > 1)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "Only 1D problems have been implemented.";
    throw std::runtime_error(err.str());
  }

  // Check the discretization
  if (discretization == nullptr)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "No discretization attached to the solver.";
    throw std::runtime_error(err.str());
  }
  if (discretization->type != DiscretizationMethod::FINITE_VOLUME)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "Only the finite volume method has been implemented.";
    throw std::runtime_error(err.str());
  }
}



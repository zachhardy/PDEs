#include "steadystate_solver.h"

/// Check the solver to ensure valid inputs.
void diffusion::SteadyStateSolver::check_inputs()
{
  // Check the mesh
  if (mesh == nullptr)
  {
    std::stringstream err;
    err << solver_string << "::" << __FUNCTION__ << ": "
        << "No mesh available to the solver.";
    throw std::runtime_error(err.str());
  }
  if (mesh->dim > 1)
  {
    std::stringstream err;
    err << solver_string << "::" << __FUNCTION__ << ": "
        << "Only 1D mesh capabilities are implemented.";
    throw std::runtime_error(err.str());
  }

  // Check the discretization
  if (discretization == nullptr)
  {
    std::stringstream err;
    err << solver_string << "::" << __FUNCTION__ << ": "
        << "No discretization available to the solver.";
    throw std::runtime_error(err.str());
  }
  if (discretization->type != SDMethod::FINITE_VOLUME)
  {
    std::stringstream err;
    err << solver_string << "::" << __FUNCTION__ << ": "
        << "Only finite volume discretizations are available.";
    throw std::runtime_error(err.str());
  }
}


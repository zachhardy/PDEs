#include "steadystate_solver.h"


/// Create a boundary condition for each boundary and each group.
void neutron_diffusion::SteadyStateSolver::initialize_boundaries()
{
  std::cout << "\nInitializing simulation boundaries...\n";

  if (mesh->dim == 1 and boundary_info.size() != 2)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "1D problems must have 2 boundary conditions";
    throw std::runtime_error(err.str());
  }

  for (const auto& boundary : boundary_info)
  {
    std::vector<BndryPtr> mg_bcs;
    for (size_t g = 0; g < n_groups; ++g)
    {
      BndryPtr bc;
      switch (boundary.first)
      {
        case BoundaryType::ZERO_FLUX:
        { bc = std::make_shared<DirichletBoundary>(); break; }
        case BoundaryType::REFLECTIVE:
        { bc = std::make_shared<NeumannBoundary>(); break; }
        case BoundaryType::VACUUM:
        { bc = std::make_shared<RobinBoundary>(); break; }
        case BoundaryType::DIRICHLET:
        {
          auto& bval = boundary_values[boundary.second][g][0];
          bc = std::make_shared<DirichletBoundary>(bval);
          break;
        }
        case BoundaryType::NEUMANN:
        {
          auto& bval = boundary_values[boundary.second][g][0];
          bc = std::make_shared<NeumannBoundary>(bval);
          break;
        }
        case BoundaryType::MARSHAK:
        {
          auto& bval = boundary_values[boundary.second][g][0];
          bc = std::make_shared<RobinBoundary>(bval);
          break;
        }
        case BoundaryType::ROBIN:
        {
          auto& bval = boundary_values[boundary.second][g];
          bc = std::make_shared<RobinBoundary>(bval[0], bval[1], bval[2]);
          break;
        }
      }
      mg_bcs.emplace_back(bc);
    }
    boundaries.emplace_back(mg_bcs);
  }

  std::cout << "Boundaries initialized: " << boundaries.size() << "\n";
  std::cout << "Done initializing boundaries.\n";
}

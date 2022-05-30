#include "steadystate_solver.h"
#include "macros.h"

void
NeutronDiffusion::SteadyStateSolver::
initialize_boundaries()
{
  std::cout << "Initializing simulation boundaries.\n";

  //============================== Check number of boundaries
  if (mesh->dim == 1)
    Assert(boundary_info.size() == 2,
           "1D problems must have 2 boundary conditions.");

  //============================== Check specified boundary values
  for (const auto& bndry_vals : boundary_values)
  {
    Assert(bndry_vals.size() == n_groups,
           "Specified boundary values must have as many entries as groups "
           "in the simulation.");
  }

  //================================================== Loop over boundaries
  for (const auto& boundary : boundary_info)
  {
    std::vector<BndryPtr> mg_bcs;
    for (size_t g = 0; g < n_groups; ++g)
    {
      BndryPtr bc;
      switch (boundary.first)
      {
        case BoundaryType::ZERO_FLUX:
        {
          bc = std::make_shared<DirichletBoundary>();
          break;
        }

        case BoundaryType::REFLECTIVE:
        {
          bc = std::make_shared<NeumannBoundary>();
          break;
        }

        case BoundaryType::VACUUM:
        {
          bc = std::make_shared<RobinBoundary>();
          break;
        }

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
          const auto& bval = boundary_values[boundary.second][g];
          Assert(bval.size() == 3,
                 "Fully specified Robin boundaries must have 3 "
                 "values for each group.");

          bc = std::make_shared<RobinBoundary>(bval[0], bval[1], bval[2]);
          break;
        }
      }
      mg_bcs.emplace_back(bc);
    }
    boundaries.emplace_back(mg_bcs);
  }

  std::cout << "Boundaries initialized: " << boundaries.size() << "\n";
}

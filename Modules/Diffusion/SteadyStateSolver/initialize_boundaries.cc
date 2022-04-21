#include "steadystate_solver.h"

/**
 * \brief Initialize the boundary conditions.
 *
 * Boundary conditions are stored first boundary-wise, then group-wise. This
 * allows for the complete specification of unique boundary conditions across
 * all groups. This routine essentially parses the specialized input structures
 * to create appropriate objects.
 */
void diffusion::SteadyStateSolver::initialize_boundaries()
{
  // Loop over boundaries
  int boundary_id = 0;
  for (const auto& bndry : boundary_info)
  {
    std::vector<BndryPtr> mg_boundary;
    const auto& btype = bndry.first;
    const auto& bvals = bndry.second;

    if (btype == BoundaryType::ZERO_FLUX)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<DirichletBoundary>();
        mg_boundary.emplace_back(bndry_g);
      }
    }
    else if (btype == BoundaryType::DIRICHLET)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<DirichletBoundary>(bvals[g][0]);
        mg_boundary.emplace_back(bndry_g);
      }
    }
    else if (btype == BoundaryType::REFLECTIVE)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<NeumannBoundary>();
        mg_boundary.emplace_back(bndry_g);
      }
    }
    else if (btype == BoundaryType::NEUMANN)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<NeumannBoundary>(bvals[g][0]);
        mg_boundary.emplace_back(bndry_g);
      }
    }
    else if (btype == BoundaryType::VACUUM)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<RobinBoundary>();
        mg_boundary.emplace_back(bndry_g);
      }
    }
    else if (btype == BoundaryType::MARSHAK)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<RobinBoundary>(bvals[g][0]);
        mg_boundary.emplace_back(bndry_g);
      }
    }
    else if (btype == BoundaryType::ROBIN)
    {
      for (int g = 0; g < n_groups; ++g)
      {
        auto bndry_g = std::make_shared<RobinBoundary>(bvals[g][0],
                                                       bvals[g][1],
                                                       bvals[g][2]);
        mg_boundary.emplace_back(bndry_g);
      }
    }
    boundaries.emplace_back(mg_boundary);
  }
}
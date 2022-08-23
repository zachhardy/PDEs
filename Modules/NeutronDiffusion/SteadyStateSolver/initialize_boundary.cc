#include "steadystate_solver.h"

#include <cassert>


using namespace NeutronDiffusion;


void
NeutronDiffusion::SteadyStateSolver::
initialize_boundaries()
{
  std::cout << "Initializing boundary conditions.\n";

  //============================================================
  // Check the number of boundaries
  //============================================================

  if (mesh->dimension == 1)
    assert(boundary_info.size() == 2);
  else if (mesh->dimension == 2)
    assert(boundary_info.size() == 4);

  //============================================================
  // Check the boundary values
  //============================================================

  for (const auto& bndry_vals : boundary_values)
    assert(bndry_vals.size() == n_groups);

  //============================================================
  // Create the appropriate boundary conditions based on inputs
  //============================================================

  for (const auto& boundary : boundary_info)
  {
    std::vector<BndryPtr> mg_bcs;

    for (unsigned int g = 0; g < n_groups; ++g)
    {
      BndryPtr bc;
      const auto btype = boundary.first;

      if (btype == BoundaryType::ZERO_FLUX)
        bc = std::make_shared<DirichletBoundary>();
      else if (btype == BoundaryType::REFLECTIVE)
        bc = std::make_shared<NeumannBoundary>();
      else if (btype == BoundaryType::VACUUM)
        bc = std::make_shared<RobinBoundary>();
      else
      {
        const auto& bvals = boundary_values[boundary.second][g];
        if (btype == BoundaryType::DIRICHLET)
          bc = std::make_shared<DirichletBoundary>(bvals[0]);
        else if (btype == BoundaryType::NEUMANN)
          bc = std::make_shared<NeumannBoundary>(bvals[0]);
        else if (btype == BoundaryType::MARSHAK)
          bc = std::make_shared<RobinBoundary>(bvals[0]);
        else
        {
          assert(bvals.size() == 3);
          bc = std::make_shared<RobinBoundary>(bvals[0], bvals[1], bvals[2]);
        }
      }
      mg_bcs.emplace_back(bc);
    }
    boundaries.emplace_back(mg_bcs);
  }
  std::cout << "Boundaries initialized: " << boundaries.size() << "\n";
}

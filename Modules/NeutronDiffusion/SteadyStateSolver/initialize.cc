#include "steadystate_solver.h"

#include "Discretization/FiniteVolume/fv.h"
#include "LinearSolvers/direct_solvers.h"

#include <algorithm>
#include <set>
#include <cassert>


using namespace NeutronDiffusion;


void
SteadyStateSolver::
initialize()
{
  std::cout
    << "\n*********************************************\n"
    <<   "Initializing the multi-group diffusion solver"
    << "\n*********************************************\n"
    << "Number of groups:  " << groups.size() << "\n\n";

  //============================================================
  // Check the mesh
  //============================================================

  assert(mesh != nullptr);
  assert(mesh->dimension < 3);

  //============================================================
  // Initialize the discretization
  //============================================================

  if (discretization_method == SDMethod::FINITE_VOLUME)
  {
    std::cout << "Initializing spatial discretization.\n";
    discretization = std::make_shared<FiniteVolume>(mesh);
  }
  else
    throw std::runtime_error("Invalid spatial discretization method.");

  //============================================================
  // Check the groups, initialize the material properties
  //============================================================

  // Check that there are groups, sort them, check for duplicataes
  assert(!groups.empty());
  std::sort(groups.begin(), groups.end());
  std::set<unsigned int> grps(groups.begin(), groups.end());
  assert(groups.size() == grps.size());

  // Initialize the materials
  initialize_materials();

  //============================================================
  // Initialize the boundary conditions
  //============================================================

  initialize_boundaries();

  //============================================================
  // Initialize data storage
  //============================================================

  size_t n_phi_dofs = discretization->n_dofs(n_groups);
  size_t n_precursor_dofs = discretization->n_dofs(max_precursors);

  phi.resize(n_phi_dofs, 0.0);
  precursors.resize(n_precursor_dofs, 0.0);

  A.reinit(n_phi_dofs, n_phi_dofs);
  b.resize(n_phi_dofs, 0.0);

  std::cout
    << "\nSimulation Summary:\n"
    <<   "-------------------\n"
    << "# of cells             :  " << mesh->cells.size() << std::endl
    << "# of nodes             :  " << discretization->n_nodes() << std::endl
    << "# of phi unknowns      :  " << n_phi_dofs << std::endl
    << "# of precursor unknowns:  " << n_precursor_dofs << std::endl;
}

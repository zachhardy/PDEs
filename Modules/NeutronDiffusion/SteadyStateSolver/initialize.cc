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
  std::cout << "Initializing the diffusion solver...\n";

  //============================================================
  // Check the mesh
  //============================================================

  assert(mesh != nullptr);
  assert(mesh->dimension < 3);

  //============================================================
  // Initialize the discretization
  //============================================================

  if (discretization_method == SDMethod::FINITE_VOLUME)
    discretization = std::make_shared<FiniteVolume>(mesh);
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

  size_t n_nodes = discretization->n_nodes();

  phi.resize(n_groups* n_nodes, 0.0);
  precursors.resize(max_precursors * n_nodes, 0.0);

  A.reinit(n_groups * n_nodes, n_groups * n_nodes);
  b.resize(n_groups * n_nodes, 0.0);
}

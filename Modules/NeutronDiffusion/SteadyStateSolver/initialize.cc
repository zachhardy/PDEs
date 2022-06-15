#include "steadystate_solver.h"
#include "../groupset.h"

#include "Discretization/FiniteVolume/fv.h"

#include "LinearSolvers/Direct/sparse_lu.h"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include <algorithm>
#include <set>
#include <cassert>


using namespace Math;
using namespace Math::LinearSolver;

using namespace NeutronDiffusion;


void
SteadyStateSolver::initialize()
{
  std::cout << "Initializing diffusion solver.\n";

  // If the full system is being solved, only use one groupset.
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
  {
    std::cout << "Solution technique set to full system.\n";

    groupsets.clear();
    groupsets.emplace_back(0);
    for (const size_t group : groups)
      groupsets.front().groups.emplace_back(group);

    linear_solver = std::make_shared<PETScSolver>(KSPGMRES, PCLU);
  }

  input_checks();

  //================================================== Initialize objects

  if (discretization_method == DiscretizationMethod::FINITE_VOLUME)
    discretization = std::make_shared<FiniteVolume>(mesh);
  else
    throw std::runtime_error("Invalid spatial discretization method.");

  initialize_materials();
  initialize_boundaries();

  //================================================== Initialize data storage
  size_t n_nodes = discretization->n_nodes();

  phi.resize(n_groups*n_nodes, 0.0);
  phi_ell.resize(phi.size(), 0.0);
  precursors.resize(max_precursors*n_nodes, 0.0);

  //================================================== Initialize groups

  // Initialize the groupsets
  for (auto& groupset : groupsets)
  {
    const size_t n_gsg = groupset.groups.size();

    groupset.A.reinit(n_gsg*n_nodes, n_gsg*n_nodes);
    groupset.b.resize(n_gsg*n_nodes, 0.0);
  }//for groupset
}

//######################################################################


void
SteadyStateSolver::input_checks()
{
  //================================================== Check the groups
  // Ensure groups and groupsets were added
  assert(!groups.empty());
  assert(!groupsets.empty());

  // Sort the groups
  std::sort(groups.begin(), groups.end());
  for (size_t g = 1; g < n_groups; ++g)
    assert(groups[g] - groups[g - 1] == 1);

  // Check that the groupsets contain all groups, no duplicates
  std::set<size_t> groupset_groups;
  for (auto& groupset : groupsets)
  {
    assert(!groupset.groups.empty());

    // Sort the groupset
    std::sort(groupset.groups.begin(), groupset.groups.end());
    for (size_t g = 1; g < groupset.groups.size(); ++g)
      assert(groupset.groups[g] - groupset.groups[g - 1] == 1);

    // Collect all groupset groups
    for (const auto& group : groupset.groups)
    {
      assert(groupset_groups.find(group) == groupset_groups.end());
      groupset_groups.insert(group);
    }
  }
  std::set<size_t> groups_set(groups.begin(), groups.end());
  assert(groupset_groups == groups_set);

  //================================================== Check the mesh
  assert(mesh != nullptr);
  assert(mesh->dim == 1);
}

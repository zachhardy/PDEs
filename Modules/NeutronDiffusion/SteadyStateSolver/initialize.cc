#include "steadystate_solver.h"
#include "NeutronDiffusion/Groupset/groupset.h"

#include "Discretization/FiniteVolume/fv.h"

#include "macros.h"

#include <set>


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
      groupsets[0].groups.emplace_back(group);

    // Set solver type to LU
    linear_solver_type = LinearSolverType::LU;
  }

  input_checks();

  //================================================== Initialize objects

  if (discretization_method == SDMethod::FINITE_VOLUME)
    discretization = std::make_shared<FiniteVolume>(mesh);
  else
    throw std::runtime_error("Invalid spatial discretization method.");

  initialize_materials();
  initialize_boundaries();

  //================================================== Initialize data storage
  size_t n_nodes = discretization->n_nodes();

  phi.resize(n_groups * n_nodes, 0.0);
  phi_ell.resize(phi.size(), 0.0);
  precursors.resize(max_precursors * n_nodes, 0.0);

  //================================================== Initialize groupsets
  for (auto& groupset : groupsets)
  {
    const size_t n_gsg = groupset.groups.size();
    groupset.matrix.reinit(n_gsg * n_nodes, n_gsg * n_nodes);
    groupset.rhs.resize(n_gsg * n_nodes, 0.0);
  }//for groupset
}

//######################################################################


void
SteadyStateSolver::input_checks()
{
  //================================================== Check the groups
  // Ensure groups and groupsets were added
  Assert(!groups.empty(), "Groups must be added to the solver.");
  Assert(!groupsets.empty(), "Groupsets must be added to the solver.");

  // Check that the groupsets contain all groups, no duplicates
  std::set<size_t> groupset_groups;
  for (const auto& groupset : groupsets)
  {
    Assert(!groupset.groups.empty(),
           "Groups must be added to each groupset.");

    // Collect all groupset groups
    for (const auto& group : groupset.groups)
    {
      Assert(groupset_groups.find(group) == groupset_groups.end(),
             "Duplicate group found.");

      groupset_groups.insert(group);
    }
  }

  std::set<size_t> groups_set(groups.begin(), groups.end());
  Assert(groupset_groups == groups_set,
         "The groupsets must contain all specified groups.");

  //================================================== Check the mesh
  Assert(mesh != nullptr, "No mesh found.");
  Assert(mesh->dim == 1, "Only 1D problems are implemented.");
}

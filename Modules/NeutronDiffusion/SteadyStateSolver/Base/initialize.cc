#include "steadystate_solver.h"

#include "LinearSolvers/Direct/lu.h"
#include "LinearSolvers/Direct/cholesky.h"

#include <set>

/**
 * \brief Initializer the solver.
 *
 * This routine ensures that the specified setup is valid. For example, the
 * mesh, groups, groupsets, materials, and boundaries are all checked and
 * relevant properties are initialized.
 */
void neutron_diffusion::SteadyStateSolver::initialize()
{
  std::cout << "Initializing solver...\n";

  // If the full system is being solved, only use one groupset.
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
  {
    std::cout << "Solution technique set to full system.\n";

    groupsets.clear();
    groupsets.emplace_back(0);
    for (const uint64_t group : groups)
      groupsets[0].groups.emplace_back(group);
  }

  //================================================== Initialize objects
  input_checks();
  initialize_discretization();
  initialize_materials();
  initialize_boundaries();

  //================================================== Initialize data storage
  uint64_t n_nodes = discretization->n_nodes();

  phi.resize(n_groups * n_nodes, 0.0);
  phi_ell.resize(phi.size(), 0.0);
  precursors.resize(max_precursors_per_material * n_nodes, 0.0);

  //================================================== Initialize groupsets
  for (auto& groupset : groupsets)
  {
    const uint64_t n_gsg = groupset.groups.size();
    groupset.matrix.resize(n_gsg * n_nodes, n_gsg * n_nodes);
    groupset.rhs.resize(n_gsg * n_nodes, 0.0);
  }//for groupset

  std::cout << "Done initializing solver.\n";
}

//######################################################################

/// Validate the general setup of the simulation.
void neutron_diffusion::SteadyStateSolver::input_checks()
{
  //================================================== Check the groups
  // Ensure groups and groupsets were added
  if (groups.empty() or groupsets.empty())
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "Groups and groupsets must be added to the solver "
        << "before initialization.";
    throw std::runtime_error(err.str());
  }

  // Check that the groupsets contain all groups, no duplicates
  std::set<uint64_t> groupset_groups;
  for (const auto& groupset : groupsets)
  {
    // Groupsets must have groups
    if (groupset.groups.empty())
    {
      std::stringstream err;
      err << solver_string << __FUNCTION__ << ": "
          << "No groups added to groupset " << groupset.id << ".";
      throw std::runtime_error(err.str());
    }

    // Collect all groupset groups
    for (const auto& group : groupset.groups)
    {
      // No duplicate groups are allowed
      if (groupset_groups.find(group) != groupset_groups.end())
      {
        std::stringstream err;
        err << solver_string << __FUNCTION__ << ": "
            << "Duplicate group found in groupset " << groupset.id << ".";
        throw std::runtime_error(err.str());
      }
      groupset_groups.insert(group);
    }
  }

  std::set<uint64_t> groups_set(groups.begin(), groups.end());
  if (groupset_groups != groups_set)
  {
    std::stringstream err;
    err << solver_string << __FUNCTION__ << ": "
        << "The groupsets must contain all specified groups.";
    throw std::runtime_error(err.str());
  }

  //================================================== Check the mesh
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
}



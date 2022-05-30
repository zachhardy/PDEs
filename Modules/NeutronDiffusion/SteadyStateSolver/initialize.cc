#include "steadystate_solver.h"
#include "NeutronDiffusion/Groupset/groupset.h"

#include "LinearSolvers/Direct/sparse_lu.h"
#include "LinearSolvers/Direct/sparse_cholesky.h"
#include "LinearSolvers/Iterative/jacobi.h"
#include "LinearSolvers/Iterative/gauss_seidel.h"
#include "LinearSolvers/Iterative/sor.h"
#include "LinearSolvers/Iterative/ssor.h"
#include "LinearSolvers/Iterative/cg.h"

#include "macros.h"

#include <set>

using namespace pdes::Math;
using namespace pdes::Math::LinearSolver;

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
  }

  //================================================== Initialize objects
  input_checks();
  initialize_discretization();
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

//######################################################################

std::shared_ptr<LinearSolver::LinearSolverBase>
SteadyStateSolver::initialize_linear_solver(Groupset& groupset)
{
  switch (linear_solver_type)
  {
    case LinearSolverType::LU:
      return std::make_shared<SparseLU>(groupset.matrix);

    case LinearSolverType::CHOLESKY:
      return std::make_shared<SparseCholesky>(groupset.matrix);

    case LinearSolverType::JACOBI:
      return std::make_shared<Jacobi>(groupset.matrix, linear_solver_opts);

    case LinearSolverType::GAUSS_SEIDEL:
      return std::make_shared<GaussSeidel>(groupset.matrix, linear_solver_opts);

    case LinearSolverType::SOR:
      return std::make_shared<SOR>(groupset.matrix, linear_solver_opts);

    case LinearSolverType::SSOR:
      return std::make_shared<SSOR>(groupset.matrix, linear_solver_opts);

    case LinearSolverType::CG:
      return std::make_shared<CG>(groupset.matrix, linear_solver_opts);

    default:
      throw std::runtime_error("Invalid linear solver type encountered.");
  }
}
#include "steadystate_solver_fv.h"

#include <cmath>

void
NeutronDiffusion::SteadyStateSolver_FV::
scoped_transfer(const Groupset& groupset, const Math::Vector& x,
                Math::Vector& dst)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();
  const auto n_gsg = groupset.groups.size();

  for (const auto& cell : mesh->cells)
  {
    const size_t gs_uk_map = cell.id * n_gsg;
    const size_t full_uk_map = cell.id * n_groups;
    for (size_t g = gs_i; g <= gs_f; ++g)
      dst[full_uk_map + g] = x[gs_uk_map + g];
  }
}

//######################################################################

void
NeutronDiffusion::SteadyStateSolver_FV::
scoped_copy(const Groupset& groupset, const Math::Vector& x,
            Math::Vector& dst)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  for (const auto& cell : mesh->cells)
  {
    const size_t uk_map = cell.id * n_groups;
    for (size_t g = gs_i; g <= gs_f; ++g)
      dst[uk_map + g] = x[uk_map + g];
  }
}

//######################################################################

double
NeutronDiffusion::SteadyStateSolver_FV::
compute_change(const Groupset& groupset)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  double norm = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const size_t uk_map = cell.id * n_groups;
    for (size_t g = gs_i; g <= gs_f; ++g)
    {
      const size_t dof = uk_map + g;
      double delta = std::fabs(phi[dof] - phi_ell[dof]);
      norm += delta * delta;
    }
  }
  return std::sqrt(norm);
}

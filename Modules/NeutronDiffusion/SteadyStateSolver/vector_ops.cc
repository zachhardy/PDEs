#include "steadystate_solver.h"

#include <cmath>


void
NeutronDiffusion::SteadyStateSolver::
scoped_transfer(const Groupset& groupset,
                const Vector& x, Vector& dst)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();
  const auto n_gsg = groupset.groups.size();
  const auto npc = discretization->nodes_per_cell();

  for (const auto& cell : mesh->cells)
    for (size_t i = 0; i < npc; ++i)
    {
      const size_t gs_uk_map = cell.id*npc*n_gsg + i*n_gsg;
      const size_t mg_uk_map = cell.id*npc*n_groups + i*n_groups;

      for (size_t g = gs_i; g <= gs_f; ++g)
        dst[mg_uk_map + g] = x[gs_uk_map + g];
    }
}

//###########################################################################

void
NeutronDiffusion::SteadyStateSolver::
scoped_copy(const Groupset& groupset,
            const Vector& x, Vector& dst)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();
  const auto npc = discretization->nodes_per_cell();

  for (const auto& cell : mesh->cells)
    for (size_t i = 0; i < npc; ++i)
    {
      const size_t uk_map = cell.id*npc*n_groups + i*n_groups;

      for (size_t g = gs_i; g <= gs_f; ++g)
        dst[uk_map + g] = x[uk_map + g];
    }
}

//###########################################################################

double
NeutronDiffusion::SteadyStateSolver::
compute_change(const Groupset& groupset)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();
  const auto npc = discretization->nodes_per_cell();

  double diff_norm = 0.0, ref_norm = 0.0;
  for (const auto& cell : mesh->cells)
    for (size_t i = 0; i < npc; ++i)
    {
      const size_t uk_map = cell.id*npc*n_groups + i*n_groups;
      for (size_t g = gs_i; g <= gs_f; ++g)
      {
        const size_t dof = uk_map + g;
        double delta = std::fabs(phi[dof] - phi_ell[dof]);
        diff_norm += delta*delta;
        ref_norm += std::fabs(phi[dof]*phi[dof]);
      }
    }
  return std::sqrt(diff_norm) / std::sqrt(ref_norm);
}

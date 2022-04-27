#include "steadystate_solver_fv.h"

void neutron_diffusion::SteadyStateSolver_FV::
scoped_copy(const Groupset& groupset, const math::Vector& x,
            math::Vector& destination)
{
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();
  const auto n_gsg = groupset.groups.size();

  for (const auto& cell : mesh->cells)
  {
    const size_t i = cell->id * n_gsg;
    const size_t uk_map = cell->id * n_groups;

    for (size_t g = gs_i; g <= gs_f; ++g)
      destination[uk_map + g] = x[uk_map + g];
  }
}

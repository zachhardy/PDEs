#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::compute_fission_rate()
{
  //==================== Loop over cells
  fission_rate = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const size_t uk_map = cell.id * n_groups;
    const double* nu_sigf = &xs->nu_sigma_f[0];

    for (const auto& g : groups)
      fission_rate[cell.id] += nu_sigf[g] * phi[uk_map + g];
  }
}


double
TransientSolver::compute_reactor_power()
{
  //==================== Loop over cells
  double p = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (xs->is_fissile)
      p += fission_rate[cell.id] * cell.volume;
   }
  return energy_per_fission * p;
}




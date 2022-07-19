#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::compute_fission_rate()
{
  fission_rate = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const size_t uk_map = cell.id * n_groups;
    const double* sig_f = &xs->sigma_f[0];

    for (const auto& g : groups)
      fission_rate[cell.id] += sig_f[g] * phi[uk_map + g];
  }
}


void
TransientSolver::compute_power()
{
  double p = 0.0, p_max = 0.0, volume = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (xs->is_fissile)
    {
      p += energy_per_fission * fission_rate[cell.id] * cell.volume;
      volume += cell.volume;
      p_max = std::max(p_max, fission_rate[cell.id]);
    }
  }
  power = p;
  average_power_density = power / volume;
  peak_power_density = energy_per_fission * p_max;
}





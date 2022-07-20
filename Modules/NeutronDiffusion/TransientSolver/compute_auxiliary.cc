#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::update_fission_rate()
{
  fission_rate = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const size_t uk_map = cell.id * n_groups;
    const double* sig_f = &xs->sigma_f[0];

    for (unsigned int gr = 0; gr < n_groups; ++gr)
      fission_rate[cell.id] += sig_f[groups[gr]] * phi[uk_map + gr];
  }
}


void
TransientSolver::update_temperature()
{
  const double alpha = conversion_factor;
  const double eff_dt = effective_time_step();
  average_fuel_temperature = 0.0; double V = 0.0;
  for (const auto& cell : mesh->cells)
    temperature[cell.id] =
        temperature_old[cell.id] +
        eff_dt * alpha * fission_rate[cell.id];
}


void
TransientSolver::compute_bulk_properties()
{
  const double Ef = energy_per_fission;
  const double alpha = conversion_factor;

  double V = 0.0;
  power = 0.0; average_fuel_temperature = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const double& T = temperature[cell.id];
    const double& Sf = fission_rate[cell.id];
    const double& volume = cell.volume;

    V += volume;
    power += Ef * Sf * volume;
    average_fuel_temperature += T * volume;

    peak_power_density = std::max(peak_power_density, Ef * Sf);
    peak_fuel_temperature = std::max(peak_fuel_temperature, T);
  }
  average_power_density = power / V;
  average_fuel_temperature /= V;
}

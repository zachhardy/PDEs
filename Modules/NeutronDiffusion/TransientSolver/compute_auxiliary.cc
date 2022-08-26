#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::
update_precursors()
{
  // Get effective time step
  const auto eff_dt = effective_time_step();

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const auto uk_map_g = n_groups * cell.id;
    const auto uk_map_j = max_precursors * cell.id;

    // Compute delayed fission rate
    const auto* nud_sigf = xs->nu_delayed_sigma_f.data();

    double f = 0.0;
    for (unsigned int gr = 0; gr < n_groups; ++gr)
      f += nud_sigf[groups[gr]] * phi[uk_map_g + gr];

    // Update the precursors
    const auto* lambda = xs->precursor_lambda.data();
    const auto* gamma = xs->precursor_yield.data();

    for (unsigned int j = 0; j < xs->n_precursors; ++j)
    {
      const auto coeff = 1.0 + eff_dt*lambda[j];
      const auto c_old = precursors_old[uk_map_j + j];
      precursors[uk_map_j + j] = (c_old + eff_dt*gamma[j] * f) / coeff;
    }//for precursor
  }//for cell
}


void
TransientSolver::
update_fission_rate()
{
  fission_rate = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const auto uk_map = n_groups * cell.id;
    const auto* sig_f = xs->sigma_f.data();

    for (unsigned int gr = 0; gr < n_groups; ++gr)
      fission_rate[cell.id] += sig_f[groups[gr]] * phi[uk_map + gr];
  }
}


void
TransientSolver::
update_temperature()
{
  const auto alpha = conversion_factor;
  const auto eff_dt = effective_time_step();
  average_fuel_temperature = 0.0; double V = 0.0;
  for (const auto& cell : mesh->cells)
    temperature[cell.id] =
        temperature_old[cell.id] +
        eff_dt * alpha * fission_rate[cell.id];
}


void
TransientSolver::
compute_bulk_properties()
{
  const auto Ef = energy_per_fission;

  double V = 0.0;
  power = 0.0; average_fuel_temperature = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const auto& T = temperature[cell.id];
    const auto& Sf = fission_rate[cell.id];
    const auto& volume = cell.volume;

    V += volume;
    power += Ef * Sf * volume;
    average_fuel_temperature += T * volume;

    peak_power_density = std::max(peak_power_density, Ef * Sf);
    peak_fuel_temperature = std::max(peak_fuel_temperature, T);
  }
  average_power_density = power / V;
  average_fuel_temperature /= V;
}

#include "transient_solver.h"

#include <filesystem>
#include <cassert>

using namespace NeutronDiffusion;


void
TransientSolver::initialize()
{
  KEigenvalueSolver::initialize();
  KEigenvalueSolver::execute();

  //==================== Check temporal parameters
  if (output_frequency < 0.0) output_frequency = dt;
  if (dt > output_frequency) dt = output_frequency;

  //==================== Clear output directory
  typedef std::filesystem::directory_iterator DirectoryIterator;

  if (write_outputs)
    for (const auto& entry : DirectoryIterator(output_directory))
      std::filesystem::remove_all(entry.path());

  //==================== Initialize auxiliary vector
  fission_rate.resize(mesh->cells.size(), 0.0);
  temperature.resize(mesh->cells.size(), 0.0);

  compute_initial_values();
}


void
TransientSolver::compute_initial_values()
{
  //==================== Normalize fission cross sections
  if (normalize_fission_xs)
    for (auto& xs : material_xs)
      for (size_t g = 0; g < n_groups; ++g)
      {
        xs->sigma_f[g] /= k_eff;
        xs->nu_sigma_f[g] /= k_eff;
        xs->nu_prompt_sigma_f[g] /= k_eff;
        xs->nu_delayed_sigma_f[g] /= k_eff;
      }

  //==================== Compute fuel volume
  reactor_volume = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (xs->is_fissile)
      reactor_volume += cell.volume;
  }

  //========== Normalize the scalar flux
  compute_fission_rate();
  phi.scale(power / compute_reactor_power());

  //========== Recompute auxiliary quantities
  compute_fission_rate();
  power = compute_reactor_power();
  average_power_density = power / reactor_volume;
  if (use_precursors)
    compute_precursors();

  //=========== Set old solution vectors
  step_solutions();
}

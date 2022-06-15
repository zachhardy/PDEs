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
  // Evaluate initial condition functions, if provided
  if (not initial_conditions.empty())
    evaluate_initial_conditions();

  // Otherwise, use the k-eigenvalue solver result
  else
  {
    // Normalize fission cross sections
    if (normalize_fission_xs)
      for (auto& xs : material_xs)
        for (size_t g = 0; g < n_groups; ++g)
        {
          xs->sigma_f[g] /= k_eff;
          xs->nu_sigma_f[g] /= k_eff;
          xs->nu_prompt_sigma_f[g] /= k_eff;
          xs->nu_delayed_sigma_f[g] /= k_eff;
        }
  }

  // Normalize the scalar flux
  if (normalization_method != NormalizationMethod::NONE)
  {
    double initial_power = power;
    compute_fission_rate();
    compute_power();

    if (normalization_method == NormalizationMethod::TOTAL_POWER)
      phi.scale(initial_power / power);
    else if (normalization_method == NormalizationMethod::AVERAGE_POWER)
      phi.scale(initial_power / average_power_density);
  }

  // Compute initial auxiliary quantities
  compute_fission_rate();
  compute_power();
  if (use_precursors)
  {
    if (initial_conditions.empty())
      compute_precursors();
    else
      precursors = 0.0;
  }

  // Set old solution vectors
  step_solutions();
}


void
TransientSolver::evaluate_initial_conditions()
{
  const size_t npc = discretization->nodes_per_cell();
  for (const auto& cell : mesh->cells)
  {
    const auto& nodes = discretization->nodes(cell);
    assert(nodes.size() == npc);

    for (size_t i = 0; i < npc; ++i)
    {
      const auto& node = nodes[i];
      const size_t uk_map = cell.id*npc*n_groups + i*n_groups;

      for (const auto& ic : initial_conditions)
      {
        const size_t g = ic.first;
        assert(g >= groups.front() && g <= groups.back());

        const auto f = ic.second;
        phi[uk_map + g] = f(node);
      }
    }//for node
  }//for cell
}

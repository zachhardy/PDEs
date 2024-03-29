#include "transient_solver.h"

#include <fstream>
#include <iomanip>
#include <filesystem>
#include <cassert>

using namespace NeutronDiffusion;


void
TransientSolver::initialize()
{

  KEigenvalueSolver::initialize();
  KEigenvalueSolver::execute();

  // Check temporal parameters
  if (output_frequency < 0.0) output_frequency = dt;
  if (dt > output_frequency) dt = output_frequency;

  // Clear output directory
  if (write_outputs)
  {
    std::cout
      << "Setting up output directories at " << output_directory << std::endl;

    if (not std::filesystem::is_directory(output_directory))
      std::filesystem::create_directory(output_directory);

    using DirectoryIterator = std::filesystem::directory_iterator;
    for (const auto& entry: DirectoryIterator(output_directory))
        std::filesystem::remove_all(entry.path());

    // Write the geometry file
    discretization->write(output_directory);

    // Initialize the summary file
    std::string summary_filepath = output_directory + "/summary.txt";
    std::ofstream file(summary_filepath,
                       std::ofstream::out | std::ofstream::trunc);

    file.setf(std::ios::left);
    file << "# "
         << std::setw(18) << "Time "
         << std::setw(18) << "Power "
         << std::setw(18) << "Peak Pwr Density "
         << std::setw(18) << "Avg Power Density "
         << std::setw(18) << "Peak Fuel Temp "
         << std::setw(18) << "Avg Fuel Temp" << std::endl;
    file.close();

  }

  // Check for non-static xs
  for (const auto& xs : material_xs)
    if (xs->sigma_a_function)
    {
      has_dynamic_xs = true;
      break;
    }

  // Initialize auxiliary vectors
  fission_rate.resize(mesh->cells.size(), 0.0);
  temperature.resize(mesh->cells.size(), initial_temperature);

  compute_initial_values();
}


void
TransientSolver::compute_initial_values()
{
  std::cout << "Computing initial conditions.\n";

  // Evaluate initial condition functions, if provided
  if (not initial_conditions.empty())
  {
    std::cout << "Evaluating specified initial condition functions.\n";

    for (const auto& cell : mesh->cells)
    {
      const auto& nodes = discretization->nodes(cell);
      for (unsigned int i = 0; i < nodes.size(); ++i)
      {
        const auto& node = nodes[i];
        const auto uk_map = n_groups*nodes.size()*cell.id + n_groups*i;

        for (const auto& ic : initial_conditions)
        {
          const auto g = ic.first;
          assert(g < n_groups);

          const auto f = ic.second;
          phi[uk_map + g] = f(node);
        }
      }//for node
    }//for cell
  }

  // Otherwise, use the k-eigenvalue solver result
  else
  {
    std::cout
      << "Using the k-eigenvalue problem result for the initial condition.\n";

    // Normalize fission cross sections
    if (normalize_fission_xs)
    {
      std::cout
        << "Normalizing the fission cross-sections to the k-eigenvalue "
        << "(" << k_eff << ").\n";

      for (auto &xs: material_xs)
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          xs->sigma_f[g] /= k_eff;
          xs->nu_sigma_f[g] /= k_eff;
          xs->nu_prompt_sigma_f[g] /= k_eff;
          xs->nu_delayed_sigma_f[g] /= k_eff;
        }
    }
  }

  // Normalize the scalar flux
  if (normalization_method != NormMethod::NONE)
  {
    std::cout
      << "Normalizing the initial condition to the specified "
      << (normalization_method == NormMethod::TOTAL_POWER? "total" : "average")
      << " power "
      << (normalization_method == NormMethod::TOTAL_POWER? "" : "density ")
      << "(" << initial_power << ")\n";

    update_fission_rate();
    compute_bulk_properties();

    if (normalization_method == NormMethod::TOTAL_POWER)
      phi.scale(initial_power / power);
    else if (normalization_method == NormMethod::AVERAGE_POWER)
      phi.scale(initial_power / average_power_density);
  }

  // Compute initial auxiliary quantities
  update_fission_rate();
  compute_bulk_properties();
  if (use_precursors)
  {
    if (initial_conditions.empty())
      compute_precursors();
    else
      precursors = 0.0;
  }

  // Set old solution vectors
  step_solutions();
  compute_bulk_properties();
}

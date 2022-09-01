#include "transient_solver.h"

#include <cmath>
#include <iomanip>
#include <numeric>


using namespace NeutronDiffusion;


void
TransientSolver::execute_time_step(bool reconstruct_matrices)
{
  // Update cross sections
  if (has_dynamic_xs)
  {
    const auto eff_dt = effective_time_step();
    for (const auto& cell : mesh->cells)
      cellwise_xs[cell.id].update({time + eff_dt,
                                   temperature[cell.id]});
    reconstruct_matrices = true;
  }

  if (reconstruct_matrices)
    rebuild_matrix();

  // Solve for the scalar flux
  if (algorithm == Algorithm::DIRECT)
  {
    b = 0.0;
    set_transient_source(APPLY_MATERIAL_SOURCE | APPLY_BOUNDARY_SOURCE);
    linear_solver->solve(phi, b);
  }
  else
    iterative_time_step_solve(APPLY_MATERIAL_SOURCE | APPLY_BOUNDARY_SOURCE |
                              APPLY_SCATTER_SOURCE | APPLY_FISSION_SOURCE);

    // Update the temperature and precursors
    update_fission_rate();
    update_temperature();
    if (use_precursors)
      update_precursors();

  if (time_stepping_method == TimeSteppingMethod::CRANK_NICHOLSON)
  {
    phi.sadd(2.0, -1.0, phi_old);
    temperature.sadd(2.0, -1.0, temperature_old);
    if (use_precursors)
      precursors.sadd(2.0, -1.0, precursors_old);
    update_fission_rate();
  }
}


void
TransientSolver::iterative_time_step_solve(SourceFlags source_flags)
{
  unsigned int nit;
  double change;
  bool converged;

  // Start iterations
  phi_ell = phi_old;
  for (nit = 0; nit < max_inner_iterations; ++nit)
  {
    // Compute the RHS and solve
    b = 0.0;
    set_transient_source(source_flags);
    linear_solver->solve(phi, b);

    // Convergence check, finalize iteration
    change = l1_norm(phi - phi_ell);
    converged = change < inner_tolerance;
    phi_ell = phi;

    // Print iteration information
    if (verbosity > 1)
      std::cout
        << std::left << "inner::"
        << "Iteration  " << std::setw(3) << nit << "  "
        << "Change  " << std::setw(8) << change
        << (converged? "  CONVERGED" : "  ")
        << std::endl;

    if (converged) break;
  }//for nit
}


void
TransientSolver::refine_time_step()
{
  double dP = std::fabs(power - power_old)/std::fabs(power_old);
  while (dP > refine_threshold)
  {
    dt /= 2.0;
    dt = (dt > dt_min)? dt : dt_min;

    execute_time_step(true);
    compute_bulk_properties();

    dP = std::fabs(power - power_old)/std::fabs(power_old);
    if (dt == dt_min) break;
  }
}


void
TransientSolver::coarsen_time_step()
{
  double dP = std::fabs(power - power_old)/std::fabs(power_old);
  if (dP < coarsen_threshold)
  {
    dt *= 2.0;
    if (dt > output_frequency)
      dt = output_frequency;
    rebuild_matrix();
  }
}


void
TransientSolver::step_solutions()
{
  power_old = power;
  phi_old = phi;
  temperature_old = temperature;
  if (use_precursors)
    precursors_old = precursors;
}


double
TransientSolver::effective_time_step()
{
  switch (time_stepping_method)
  {
    case TimeSteppingMethod::BACKWARD_EULER:
      return dt;
    case TimeSteppingMethod::CRANK_NICHOLSON:
      return 0.5 * dt;
    default:
      throw std::runtime_error("Invalid time stepping method.");
  }
}

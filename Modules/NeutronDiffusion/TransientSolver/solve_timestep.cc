#include "transient_solver.h"

#include <cmath>
#include <iomanip>
#include <numeric>


using namespace NeutronDiffusion;


void
TransientSolver::solve_time_step()
{
  phi_ell = phi_old;

  // Update cross sections
  if (has_dynamic_xs)
  {
    const auto eff_dt = effective_time_step();
    for (const auto& cell : mesh->cells)
      cellwise_xs[cell.id].update(time + eff_dt);
    assemble_matrices();
  }

  // Solve for the scalar flux
  if (solution_technique == SolutionTechnique::FULL_SYSTEM)
    solve_full_system_time_step(APPLY_MATERIAL_SOURCE);
  else
    for (auto& groupset: groupsets)
      solve_groupset_time_step(
        groupset,
        APPLY_MATERIAL_SOURCE |
        APPLY_WGS_SCATTER_SOURCE | APPLY_WGS_FISSION_SOURCE |
        APPLY_AGS_SCATTER_SOURCE | APPLY_AGS_FISSION_SOURCE);

  // Update the precursors
  if (use_precursors)
    update_precursors();

  if (time_stepping_method == TimeSteppingMethod::CRANK_NICHOLSON)
  {
    phi.sadd(2.0, -1.0, phi_old);
    if (use_precursors)
      precursors.sadd(2.0, -1.0, precursor_old);
    compute_fission_rate();
  }
}

//######################################################################

void
TransientSolver::solve_groupset_time_step(Groupset& groupset,
                                          SourceFlags source_flags)
{
  Vector& b = groupset.b;

  size_t nit;
  double change;
  bool converged = false;

  // Start iterations
  for (nit = 0; nit < groupset.max_iterations; ++nit)
  {
    // Compute the RHS and solve
    b = 0.0;
    set_transient_source(groupset, source_flags);
    auto x = linear_solver->solve(b);

    // Convergence check, finalize iteration
    scoped_transfer(groupset, x, phi);
    change = compute_change(groupset);
    scoped_copy(groupset, phi, phi_ell);

    if (change < groupset.tolerance)
      converged = true;

    // Print iteration information
    if (verbosity > 1)
    {
      std::stringstream iter_info;
      iter_info
        << std::left << "SourceIteration::"
        << "Step  " << std::setw(3) << nit << "   "
        << "Value  " << change;
      if (converged) iter_info << "  CONVERGED";
      std::cout << iter_info.str() << std::endl;
    }

    if (converged) break;
  }//for nit
  compute_fission_rate();
}

//######################################################################

void
TransientSolver::solve_full_system_time_step(SourceFlags source_flags)
{
  groupsets.front().b = 0.0;
  set_transient_source(groupsets.front(), source_flags);
  phi = linear_solver->solve(groupsets.front().b);
  compute_fission_rate();
}

//######################################################################

void
TransientSolver::refine_time_step()
{
  double dP = std::fabs(power - power_old)/std::fabs(power_old);
  while (dP > refine_threshold)
  {
    dt /= 2.0;
    assemble_matrices();

    solve_time_step();
    compute_power();

    dP = std::fabs(power - power_old)/std::fabs(power_old);
  }
}

//######################################################################

void
TransientSolver::coarsen_time_step()
{
  double dP = std::fabs(power - power_old)/std::fabs(power_old);
  if (dP < coarsen_threshold)
  {
    dt *= 2.0;
    if (dt > output_frequency)
      dt = output_frequency;
    assemble_matrices();
  }
}

//######################################################################

void
TransientSolver::step_solutions()
{
  power_old = power;
  phi_old = phi_ell = phi;
  temperature_old = temperature;
  if (use_precursors)
    precursor_old = precursors;
}

//######################################################################

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

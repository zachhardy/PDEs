#include "transient_solver.h"

#include <iomanip>
#include <cmath>

using namespace NeutronDiffusion;


void
TransientSolver::execute()
{
  const double eps = 1.0e-10;

  // Output settings
  size_t output = 0;
  double next_output = output_frequency;
  if (write_outputs)
    write(output++);

  // Initialize matrices
  assemble_matrices();

  const double dt_initial = dt;
  bool reconstruct_matrices;

  // Time stepping loop
  time = t_start;
  size_t step = 0;
  while (time < t_end - eps)
  {
    /* This flag is used to tell the solve_time_step routine whether or not
     * to reconstruct the matrices or not. This gets set to true when the
     * time step changes or the cross-sections are modified. */
    reconstruct_matrices = false;

    //==================================================
    // Modify time steps to coincide with output times and
    // the end of the simulation.
    //==================================================

    if (write_outputs && time + dt > next_output + eps)
    {
      dt = next_output - time;
      reconstruct_matrices = true;
    }

    if (time + dt > t_end + eps)
    {
      dt = t_end - time;
      reconstruct_matrices = true;
    }

    //==================================================
    // Solve the time step
    //==================================================

    solve_time_step(reconstruct_matrices);
    compute_bulk_properties();

    if (adaptivity)
      refine_time_step();

    //==================================================
    // Postprocess the time step
    //==================================================

    time += dt;
    ++step;

    // Output solutions
    if (std::fabs(time - next_output) < eps)
    {
      write(output++);
      next_output += output_frequency;
      if (next_output > t_end || std::fabs(next_output - t_end) < eps)
        next_output = t_end;
    }

    // Check to see if time steps should be coarsened.
    if (adaptivity)
      coarsen_time_step();

    // If no adaptivity, reset the time step to the original when
    // modified for outputting purposes
    else if (dt != dt_initial)
    {
      dt = dt_initial;
      assemble_matrices();
    }

    // Move the solutions to the next time step
    step_solutions();

    std::cout
      << "\n***** Time Step " << step << " *****\n"
      << "Simulation Time:        " << time << " s\n"
      << "Time Step Size :        " << dt << " s\n"
      << "Reactor Power  :        " << power << " W\n"
      << "Average Power Density:  " << average_power_density << " W/cm^3\n"
      << "Average Temperature:    " << average_fuel_temperature << " K\n";
  }

  // Reset dt to see the initial time step size
  dt = dt_initial;
}

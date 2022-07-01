#include "transient_solver.h"

#include <iomanip>
#include <cmath>

using namespace NeutronDiffusion;


void
TransientSolver::execute()
{
  // Output settings
  size_t output = 0;
  double next_output = output_frequency;
  if (write_outputs)
    write_snapshot(output++);

  // Initialize matrices
  assemble_matrices();

  double dt_prev = dt;
  const double dt_initial = dt;

  // Time stepping loop
  time = t_start;
  size_t step = 0;
  while (time < t_end - 1.0e-12)
  {
    //========================================
    // Modify time steps to coincide with output
    // times and the end of the simulation.
    //========================================

    if (write_outputs && time + dt > next_output)
    {
      dt_prev = dt;
      dt = next_output - time;
    }

    if (time + dt > t_end)
      dt = t_end - time;

    //========================================
    // Solve the time step
    //========================================

    solve_time_step();
    compute_power();

    //========================================
    // Postprocess the time step
    //========================================

    time += dt;
    ++step;
    step_solutions();

    // Output solutions
    if (std::fabs(time - next_output) < 1.0e-12)
    {
      write_snapshot(output++);
      next_output += output_frequency;
      if (next_output > t_end ||
          std::fabs(next_output - t_end) < 1.0e-12)
        next_output = t_end;
    }

    std::cout
      << "\n***** Time Step " << step << " *****\n"
      << "Simulation Time:   " << time << " s\n"
      << "Time Step Size :   " << dt << " s\n"
      << "Reactor Power  :   " << power << " W\n";

    // Reset dt if changed for an output
    dt = dt_prev;
  }
  // Reset dt to see the initial time step size
  dt = dt_initial;
}
#include "transient_solver.h"

#include <iomanip>

using namespace NeutronDiffusion;


void
TransientSolver::execute()
{
  // Initialize matrices
  assemble_matrices();

  double next_output = output_frequency;
  const double dt_initial = dt;

  // Time stepping loop
  time = t_start;
  size_t step = 0, output = 0;
  while (time < t_end - 1.0e-12)
  {
    // Force coincidence with output times
    if (write_outputs && time + dt > next_output)
      dt = next_output - time;

    // Force coincidence with t_end
    if (time + dt > t_end)
      dt = t_end - time;

    solve_time_step();
    compute_power();

    time += dt;
    ++step;

    step_solutions();

    if (verbosity > 0)
      std::cout
        << "\n***** Time Step " << step << " *****\n"
        << "Simulation Time:   " << time << " s\n"
        << "Time Step Size :   " << dt << " s\n"
        << "Reactor Power  :   " << power << " W\n";
  }
}
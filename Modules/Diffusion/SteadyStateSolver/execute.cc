#include "steadystate_solver.h"

#include "LinearSolvers/gauss_elimination.h"

#include <iomanip>

void diffusion::SteadyStateSolver::execute()
{
  std::cout << "Executing the solver." << std::endl;

  assemble_matrix();
  assemble_rhs_vector();

  GaussElimination linear_solver(system_matrix, system_rhs, false, true);
  linear_solver.setup();
  phi = linear_solver.solve();
}
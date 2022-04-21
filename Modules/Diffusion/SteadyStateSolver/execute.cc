#include "steadystate_solver.h"
#include "LinearSolvers/lu.h"
#include <iomanip>

void diffusion::SteadyStateSolver::execute()
{
  std::cout << "Executing the solver." << std::endl;

  assemble_matrix();
  assemble_rhs_vector();

  math::LU solver(system_matrix);
  phi = solver.solve(system_rhs);

  std::cout << "Solution:\t" << phi.to_string();
}
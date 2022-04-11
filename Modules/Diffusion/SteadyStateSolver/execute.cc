#include "steadystate_solver.h"

#include <iomanip>

void diffusion::SteadyStateSolver::execute()
{
  std::cout << "Executing the solver." << std::endl;

  assemble_matrix();
  assemble_rhs_vector();

  std::cout << "Matrix:\n" << system_matrix.to_string();
  std::cout << "Vector:\n" << system_rhs.to_string();
}
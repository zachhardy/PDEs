#include "steadystate_solver.h"

#include <iomanip>

void diffusion::SteadyStateSolver::execute()
{
  std::cout << "Executing the solver." << std::endl;

  assemble_matrix();
  assemble_rhs_vector();

}
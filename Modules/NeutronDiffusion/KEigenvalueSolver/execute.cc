#include "keigenvalue_solver.h"


using namespace NeutronDiffusion;


void
KEigenvalueSolver::execute()
{
  std::cout
      << "\n*******************************************************\n"
      <<   "Executing the multi-group diffusion k-eigenvalue solver"
      << "\n*******************************************************\n";

  if (algorithm == Algorithm::DIRECT)
    assemble_matrix(ASSEMBLE_SCATTER);
  else
    assemble_matrix(NO_ASSEMBLER_FLAGS);
  linear_solver->set_matrix(A);

  power_method();

  if (use_precursors)
  {
    compute_precursors();
    precursors /= k_eff;
  }
}
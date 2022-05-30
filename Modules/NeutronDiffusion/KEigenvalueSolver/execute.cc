#include "keigenvalue_solver.h"


using namespace NeutronDiffusion;


void
KEigenvalueSolver::execute()
{
  for (auto& groupset : groupsets)
  {
    if (solution_technique == SolutionTechnique::GROUPSET_WISE)
      assemble_matrix(groupset);
    else
      assemble_matrix(groupset, ASSEMBLE_SCATTER);
  }

  power_method();

  if (use_precursors)
  {
    compute_precursors();
    precursors /= k_eff;
  }


}
#include "keigenvalue_solver.h"


using namespace NeutronDiffusion;


void
KEigenvalueSolver::execute()
{
  solution_technique = SolutionTechnique::GROUPSET_WISE;
  for (auto& groupset : groupsets)
    assemble_matrix(groupset);

  power_method();

  if (use_precursors)
  {
    compute_precursors();
    precursors /= k_eff;
  }


}
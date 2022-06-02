#include "steadystate_solver.h"


using namespace NeutronDiffusion;


void
SteadyStateSolver::
assemble_matrix(Groupset& groupset, AssemblerFlags assembler_flags)
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    {
      fv_assemble_matrix(groupset, assembler_flags);
      break;
    }
    default:throw std::runtime_error("Invalid discretization method.");
  }
}


void
SteadyStateSolver::
set_source(Groupset& groupset, SourceFlags source_flags)
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    {
      fv_set_source(groupset, source_flags);
      break;
    }
    default:throw std::runtime_error("Invalid discretization method.");
  }
}


void SteadyStateSolver::
compute_precursors()
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    {
      fv_compute_precursors();
      break;
    }
    default:throw std::runtime_error("Invalid discretization method.");
  }
}
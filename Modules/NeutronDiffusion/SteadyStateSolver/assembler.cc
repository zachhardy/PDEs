#include "steadystate_solver.h"


using namespace NeutronDiffusion;


void
SteadyStateSolver::assemble_matrix(Groupset& groupset)
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    { fv_assemble_matrix(groupset); break; }
    default:
      throw std::runtime_error("Unimplemented discretization method.");
  }
}


void
SteadyStateSolver::
set_source(Groupset& groupset, SourceFlags source_flags)
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    { fv_set_source(groupset, source_flags); break; }
    default:
      throw std::runtime_error("Unimplemented discretization method.");
  }
}


void SteadyStateSolver::
compute_precursors()
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
    { fv_compute_precursors(); break; }
    default:
      throw std::runtime_error("Unimplemented discretization method.");
  }
}
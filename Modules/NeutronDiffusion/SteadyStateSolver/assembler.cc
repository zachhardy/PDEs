#include "steadystate_solver.h"


using namespace NeutronDiffusion;


void
SteadyStateSolver::assemble_matrix(Groupset& groupset)
{ fv_assemble_matrix(groupset); }


void
SteadyStateSolver::
set_source(Groupset& groupset, SourceFlags source_flags)
{ fv_set_source(groupset, source_flags); }


void SteadyStateSolver::
compute_precursors()
{ fv_compute_precursors(); }
#include "steadystate_solver_fv.h"

#include "macros.h"


void
NeutronDiffusion::SteadyStateSolver_FV::
compute_precursors()
{
  Assert(use_precursors,
         "Cannot call compute_precursors when use_precursors "
         "flag is set to false.");

  precursors = 0.0;

  //======================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.id;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    if (!xs->is_fissile)
      continue;

    const double* lambda = xs->precursor_lambda.data();
    const double* gamma = xs->precursor_yield.data();
    const double* nud_sigf = xs->nu_delayed_sigma_f.data();
    const double* x = &phi[cell.id * n_groups];
    const size_t uk_map = cell.id * max_precursors;

    for (size_t j = 0; j < xs->n_precursors; ++j)
      for (size_t g = 0; g < n_groups; ++g)
        precursors[uk_map + j] += gamma[j]/lambda[j]*nud_sigf[g] * x[g];
  }
}
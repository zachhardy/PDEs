#include "steadystate_solver.h"

#include <cassert>


using namespace NeutronDiffusion;


void
SteadyStateSolver::compute_precursors()
{
  assert(use_precursors);

  precursors = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const double* lambda = xs->precursor_lambda.data();
    const double* gamma = xs->precursor_yield.data();
    const double* nud_sigf = xs->nu_delayed_sigma_f.data();

    const size_t uk_map = cell.id*n_groups;
    const size_t prec_uk_map = cell.id*max_precursors;

    // Loop over precursors
    for (size_t j = 0; j < xs->n_precursors; ++j)
    {
      double value = 0.0;
      const double coeff = gamma[j]/lambda[j];
      for (size_t gr = 0; gr < n_groups; ++gr)
      {
        const size_t g = groups[gr];
        value += coeff*nud_sigf[g] * phi[uk_map + g];
      }
      precursors[prec_uk_map + j] = value;
    }
  }
}
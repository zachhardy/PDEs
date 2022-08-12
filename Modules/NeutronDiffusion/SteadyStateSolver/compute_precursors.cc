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

    const auto* lambda = xs->precursor_lambda.data();
    const auto* gamma = xs->precursor_yield.data();
    const auto* nud_sigf = xs->nu_delayed_sigma_f.data();

    const auto uk_map_g = n_groups * cell.id;
    const auto uk_map_j = max_precursors * cell.id;

    // Loop over precursors
    for (unsigned int j = 0; j < xs->n_precursors; ++j)
    {
      double value = 0.0;
      const auto coeff = gamma[j]/lambda[j];
      for (unsigned int gr = 0; gr < n_groups; ++gr)
        value += coeff * nud_sigf[groups[gr]] * phi[uk_map_g + gr];
      precursors[uk_map_j + j] = value;
    }
  }
}
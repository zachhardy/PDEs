#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::update_precursors()
{
  // Get effective time step
  const double eff_dt = effective_time_step();

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    if (not xs->is_fissile)
      continue;

    const size_t uk_map = cell.id * n_groups;
    const size_t prec_uk_map = cell.id * max_precursors;

    // Compute delayed fission rate
    const double* nud_sigf = xs->nu_delayed_sigma_f.data();

    double f = 0.0;
    for (const auto& g : groups)
      f += nud_sigf[g] * phi[uk_map + g];

    // Update the precursors
    const double* lambda = xs->precursor_lambda.data();
    const double* gamma = xs->precursor_yield.data();

    for (size_t j = 0; j < xs->n_precursors; ++j)
    {
      const double coeff = 1.0 + eff_dt*lambda[j];
      const double c_old = precursor_old[prec_uk_map + j];

      precursors[prec_uk_map + j] =
        (c_old + eff_dt*gamma[j] * f) / coeff;
    }//for precursor
  }//for cell
}

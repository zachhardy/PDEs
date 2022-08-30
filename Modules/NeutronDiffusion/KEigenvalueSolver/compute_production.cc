#include "keigenvalue_solver.h"


using namespace NeutronDiffusion;


double
KEigenvalueSolver::compute_production()
{
  double production = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const auto volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    if (xs->is_fissile)
    {
      const auto uk_map = n_groups * cell.id;
      const auto* nu_sigf = xs->nu_sigma_f.data();

      double cell_production = 0.0;
      for (unsigned int g = 0; g < n_groups; ++g)
        cell_production += nu_sigf[g] * phi[uk_map + g];
      production += cell_production*volume;
    }
  }
  return production;
}

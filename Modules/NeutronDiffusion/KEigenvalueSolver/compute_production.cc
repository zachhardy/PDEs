#include "keigenvalue_solver.h"


using namespace NeutronDiffusion;


double
KEigenvalueSolver::compute_production()
{
  switch (discretization_method)
  {
    case DiscretizationMethod::FINITE_VOLUME:
      return fv_compute_production();
    default:
      throw std::runtime_error("Invalid discretization method.");
  }
}


double
KEigenvalueSolver::fv_compute_production()
{
  double production = 0.0;
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    if (xs->is_fissile)
    {
      const size_t uk_map = cell.id*n_groups;
      const double* x = &phi[uk_map];
      const double* nu_sigf = &xs->nu_sigma_f[0];

      double cell_production = 0.0;
      for (const auto& g : groups)
        cell_production += nu_sigf[g]*x[g];
      production += cell_production*volume;
    }
  }
  return production;
}
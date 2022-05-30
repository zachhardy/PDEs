#include "keigenvalue_solver.h"

#include <cmath>
#include <iomanip>


using namespace NeutronDiffusion;


void
KEigenvalueSolver::power_method()
{
  std::cout << "********** Solving the k-eigenvalue problem "
            << "using the Power Method.\n";

  phi = phi_ell = 1.0;

  double f_prev = 1.0;
  double k_eff_prev = k_eff;
  double f_new, k_eff_change, rho;

  size_t nit;
  bool converged = false;

  for (nit = 0; nit < max_iterations; ++nit)
  {
    //==================== Solve each groupset
    for (auto& groupset : groupsets)
    {
      // Precompute the fission source
      groupset.rhs = 0.0;
      set_source(groupset,
                 APPLY_WGS_FISSION_SOURCE |
                 APPLY_AGS_FISSION_SOURCE);
      groupset.rhs /= k_eff;

      // Converge the scattering source
      solve_groupset(groupset,
                     APPLY_WGS_SCATTER_SOURCE |
                     APPLY_AGS_SCATTER_SOURCE);
    }

    //==================== Recompute the k-eigenvalue
    f_new = compute_production();
    k_eff *= f_new / f_prev;
    rho = (k_eff - 1.0)/k_eff;

    k_eff_change = std::fabs(k_eff - k_eff_prev) / k_eff;
    f_prev = f_new;
    k_eff_prev = k_eff;

    if (k_eff_change <= tolerance)
      converged = true;

    // Print iteration information
    std::stringstream iter_info;
    iter_info << "Iteration: " << std::setw(3) << nit << "   "
              << "k_eff: " << std::setw(6) << k_eff << "   "
              << "k_eff Change: " << std::setw(6) << k_eff_change << "   "
              << "Reactivity: " << std::setw(6) << rho * 1.0e5;
    if (converged) iter_info << "   CONVERGED\n";
    std::cout << iter_info.str() << "\n";

    if (converged) break;
  }

  std::cout << "\n"
            << "          Final k-eigenvalue     :     "
            <<  std::setprecision(6) << k_eff << "\n"
            << "          Final Change           :     "
            << std::setprecision(6) << k_eff_change << "\n";
}
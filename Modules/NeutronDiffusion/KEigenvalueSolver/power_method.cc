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

  double production_ell = 1.0;
  double k_eff_ell = k_eff;
  double production, rho;
  double k_eff_change, phi_change;

  size_t nit;
  bool converged = false;

  for (nit = 0; nit < max_iterations; ++nit)
  {
    //==================== Precompute the fission source
    for (auto& groupset : groupsets)
    {
      groupset.b = 0.0;
      set_source(groupset, APPLY_WGS_FISSION_SOURCE |
                           APPLY_AGS_FISSION_SOURCE);
      groupset.b /= k_eff;
    }

    //==================== Solve the system
    if (solution_technique == SolutionTechnique::GROUPSET_WISE)
      for (auto& groupset : groupsets)
      {
        // Converge the scattering source
        solve_groupset(groupset,
                       APPLY_WGS_SCATTER_SOURCE |
                       APPLY_AGS_SCATTER_SOURCE);
      }
    else
      solve_full_system(NO_SOURCE_FLAGS);


    //==================== Recompute the k-eigenvalue
    production = compute_production();
    k_eff *= production / production_ell;
    rho = (k_eff - 1.0)/k_eff;

    k_eff_change = std::fabs(k_eff - k_eff_ell) / k_eff;
    phi_change = l2_norm(phi - phi_ell);

    production_ell = production;
    k_eff_ell = k_eff;
    phi_ell = phi;

    if (k_eff_change <= tolerance &&
        phi_change <= tolerance)
      converged = true;

    // Print iteration information
    if (verbosity > 0)
    {
      std::stringstream iter_info;
      iter_info
        << std::left << "k-Eigenvalue::"
        << "Step  " << std::setw(4) << nit
        << "k_eff  " << std::setw(8) << k_eff
        << "k_eff Change  " << std::setw(14) << k_eff_change
        << "Reactivity  " << rho * 1.0e5;
      if (converged) iter_info << "   CONVERGED";
      std::cout << iter_info.str() << std::endl;
    }

    if (converged) break;
  }

  std::stringstream summary;
  if (converged)
    summary << "\n***** k-Eigenvalue Solver Converged! *****\n";
  else
    summary << "\n!!*!! WARNING: k-Eigenvalue Solver NOT Converged !!*!!\n";

  summary << "Final k-Eigenvalue:         "
          << std::left << std::setw(6) << k_eff << "\n"
          << "Final k-Eigenvalue Change:  "
          << std::left << std::setw(6) << k_eff_change << "\n"
          << "Final Phi Change:           "
          << std::left << std::setw(6) << phi_change << "\n"
          << "# of Iterations:            " << nit;
  std::cout << summary.str() << std::endl;
}
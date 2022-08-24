#include "keigenvalue_solver.h"

#include <cmath>
#include <iomanip>


using namespace NeutronDiffusion;


void
KEigenvalueSolver::power_method()
{
  std::cout << "\n********** Solving the k-eigenvalue problem "
            << "using the Power Method.\n\n";

  phi = 1.0;
  auto x = phi;

  auto production = compute_production();
  auto production_ell = production;
  auto k_eff_ell = k_eff;

  double rho;
  double k_eff_change, phi_change;

  unsigned int nit;
  bool converged = false;

  for (nit = 0; nit < max_outer_iterations; ++nit)
  {
    //========================================
    // Precompute the fission source
    //========================================

    b = 0.0;
    set_source(APPLY_FISSION_SOURCE | APPLY_BOUNDARY_SOURCE);
    b /= k_eff;

    //========================================
    // Solve the system
    //========================================
    unsigned int inner_nit;
    if (algorithm == Algorithm::DIRECT)
      linear_solver->solve(phi, b);
    else
      inner_nit = iterative_solve(APPLY_SCATTER_SOURCE);

    //========================================
    // Recompute the k-eigenvalue
    //========================================

    production = compute_production();
    k_eff *= production/production_ell;
    rho = (k_eff - 1.0)/k_eff;

    k_eff_change = std::fabs(k_eff - k_eff_ell)/k_eff;
    phi_change = l1_norm(phi - x) / l1_norm(phi);

    production_ell = production;
    k_eff_ell = k_eff;
    x = phi;

    converged = (k_eff_change < outer_tolerance &&
                 phi_change < outer_tolerance);

    // Print iteration information
    if (verbosity > 0)
    {
      std::stringstream iter_info;
      iter_info
        << std::left << "outer::"
        << "Iteration  " << std::setw(4) << nit
        << "k_eff  " << std::setw(10) << k_eff
        << "k_eff Change  " << std::setw(14) << k_eff_change
        << "Phi Change  " << std::setw(14) << phi_change
        << "Reactivity  " << rho*1.0e5;
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
          << std::left << std::setw(6) << k_eff << std::endl
          << "Final k-Eigenvalue Change:  "
          << std::left << std::setw(6) << k_eff_change << std::endl
          << "Final Phi Change:           "
          << std::left << std::setw(6) << phi_change << std::endl
          << "# of Iterations:            " << nit << std::endl;
  std::cout << summary.str() << std::endl;
}
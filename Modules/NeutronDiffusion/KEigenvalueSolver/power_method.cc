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
  phi_ell = phi;
  auto phi_tmp = phi_ell;

  auto production = compute_production();
  auto production_ell = production;
  auto k_eff_ell = k_eff;

  unsigned int nit;
  bool converged = false;
  double k_eff_change, phi_change;
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
    std::pair<unsigned int, double> inner_result;
    if (algorithm == Algorithm::DIRECT)
      linear_solver->solve(phi, b);
    else
      inner_result = iterative_solve(APPLY_SCATTER_SOURCE);

    //========================================
    // Recompute the k-eigenvalue
    //========================================

    production = compute_production();
    k_eff *= production/production_ell;

    k_eff_change = std::fabs(k_eff - k_eff_ell)/k_eff;
    phi_change = l1_norm(phi - phi_tmp) / l1_norm(phi);

    production_ell = production;
    k_eff_ell = k_eff;
    phi_tmp = phi;

    converged = (k_eff_change < outer_tolerance &&
                 phi_change < outer_tolerance);

    // Print iteration information
    if (verbosity > 0)
      std::cout
        << std::left << "outer::"
        << "Iteration  " << std::setw(4) << nit
        << "k_eff  "
        << std::setprecision(6) << std::setw(10) << k_eff
        << "k_eff Change  "
        << std::setprecision(6) << std::setw(14) << k_eff_change
        << "Phi Change  "
        << std::setprecision(6) << std::setw(14) << phi_change
        << (converged? "CONVERGED" : "") << std::endl;

    if (converged) break;
  }

  std::cout
    << (converged?
        "\n***** k-Eigenvalue Solver Converged! *****\n" :
        "\n!!*!! WARNING: k-Eigenvalue Solver NOT Converged !!*!!\n")
    << "Final k-Eigenvalue:         "
    << std::left << std::setw(6) << k_eff << std::endl
    << "Iterations:                 "
    << std::left << std::setw(3) << nit << std::endl
    << "Final k-Eigenvalue Change:  "
    << std::left << std::setw(6) << k_eff_change << std::endl
    << "Final Phi Change:           "
    << std::left << std::setw(6) << phi_change << std::endl
    << std::endl;
}
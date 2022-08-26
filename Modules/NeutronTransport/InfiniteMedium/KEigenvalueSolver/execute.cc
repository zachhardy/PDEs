#include "keigenvalue_solver.h"

#include <cmath>
#include <iomanip>
#include <cassert>


using namespace NeutronTransport;
using namespace InfiniteMedium;


void KEigenvalueSolver::execute()
{
  assert(xs->is_fissile);

  phi = 1.0;
  phi_ell = phi;
  auto phi_tmp = phi_ell;

  auto production = compute_production();
  auto production_ell = production;
  auto k_eff_ell = k_eff;

  unsigned int nit;
  double k_eff_change, phi_change;
  for (nit = 0; nit < max_outer_iterations; ++nit)
  {
    q_moments = 0.0;
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int gp = 0; gp < n_groups; ++gp)
        q_moments[g] += xs->chi[g]*xs->nu_sigma_f[gp]*phi[gp];
    q_moments /= k_eff;

    auto result = source_iterations();
    if (result.first == max_inner_iterations)
      std::cout
        << "!!*!! WARNING: Inner iterations did not converge... "
        << "Final difference was " << result.second << ". !!*!!\n";

    production = compute_production();
    k_eff *= production/production_ell;
    double rho = (k_eff - 1.0)/k_eff;

    k_eff_change = std::fabs(k_eff - k_eff_ell)/k_eff;
    phi_change = l1_norm(phi - phi_tmp)/l1_norm(phi);

    production_ell = production;
    k_eff_ell = k_eff;
    phi_tmp = phi = phi_ell;

    bool converged = (k_eff_change < outer_tolerance &&
                      phi_change < outer_tolerance);

    if (verbosity > 0)
      std::cout
        << std::left << "outer::"
        << "Iteration  " << std::setw(4) << nit
        << "k_eff  " << std::setw(12) << k_eff
        << "k_eff Change  " << std::setw(12) << k_eff_change
        << "Phi Change  " << std::setw(12) << phi_change
        << "Reactivity  " << rho*1.0e5
        << (converged? "  CONVERGED" : "  ")
        << std::endl;

    if (converged) break;
  }//for nit

  std::cout
    << (nit < max_outer_iterations?
        "\n***** k-Eigenvalue Solver Converged! *****\n" :
        "\n!!*!! WARNING: k-Eigenvalue Solver NOT Converged !!*!!\n")
    << "Final k-Eigenvalue:         "
    << std::left << std::setw(6) << k_eff << std::endl
    << "# of Iterations:            " << nit << std::endl
    << "Final k-Eigenvalue Change:  "
    << std::left << std::setw(6) << k_eff_change << std::endl
    << "Final Phi Change:           "
    << std::left << std::setw(6) << phi_change << std::endl
    << std::endl;
}


double
KEigenvalueSolver::compute_production()
{
  double production = 0.0;
  for (unsigned int g = 0; g < n_groups; ++g)
    production += xs->nu_sigma_f[g] * phi[g];
  return production;
}

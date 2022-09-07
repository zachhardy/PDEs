#include "steadystate_solver.h"

#include <cassert>
#include <iomanip>


using namespace PDEs;
using namespace Math;
using namespace InfiniteMedium;


void
SteadyStateSolver::execute()
{
  std::cout
      << "\n**********************************************************\n"
      << "Executing the multi-group infinite medium transport solver"
      << "\n**********************************************************\n";

  auto result = source_iterations(APPLY_MATERIAL_SOURCE |
                                  APPLY_SCATTER_SOURCE |
                                  APPLY_FISSION_SOURCE);

  bool converged = result.first != max_inner_iterations;

  if (verbosity > 0)
    std::cout
        << (converged ?
            "\n***** Steady State Solver CONVERGED *****\n" :
            "\n!!*!! WARNING: Steady State Solver NOT CONVERGED *****\n")
        << "# of Iterations:  "
        << std::left << std::setw(6) << result.first << std::endl
        << "Final Phi Change  "
        << std::left << std::setw(6) << result.second << std::endl
        << std::endl;
}


std::pair<unsigned int, double>
SteadyStateSolver::source_iterations(SourceFlags source_flags)
{
  phi_ell = phi;
  auto q_moments_init = q_moments;

  unsigned int nit;
  double change;
  for (nit = 0; nit < max_inner_iterations; ++nit)
  {
    q_moments = q_moments_init;
    set_source(source_flags);
    sweep();
    if (use_dsa) dsa();
    update_neutron_density();

    change = l1_norm(phi - phi_ell);
    double balance = check_balance();
    phi_ell = phi;

    if (verbosity > 1)
      std::cout
          << std::left << "inner::"
          << "Iteration  " << std::setw(3) << nit << "   "
          << "Value  " << std::setw(10) << change << "  "
          << "Balance  " << std::setw(10) << balance << "  "
          << (change < inner_tolerance ? "CONVERGED" : "")
          << std::endl;

    if (change < inner_tolerance)
      break;
  }
  return {nit, change};
}


void
SteadyStateSolver::sweep()
{
  phi = 0.0;
  for (unsigned int n = 0; n < n_angles; ++n)
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      // Solve for the angular flux
      double x = 0.0;
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        x += moment_to_discrete[ell][n] * q_moments[ell * n_groups + g];
      x /= density * xs->sigma_t[g];

      // Store the angular flux
      psi[n * n_groups + g] = x;

      // Accumulate flux moment
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        phi[ell * n_groups + g] += discrete_to_moment[ell][n] * x;
    }//for g
}


void
SteadyStateSolver::set_source(SourceFlags source_flags)
{
  const bool apply_mat_src = (source_flags & APPLY_MATERIAL_SOURCE);
  const bool apply_scatter_src = (source_flags & APPLY_SCATTER_SOURCE);
  const bool apply_fission_src = (source_flags & APPLY_FISSION_SOURCE);

  // Loop over moments
  for (unsigned int ell = 0; ell < n_moments; ++ell)
  {
    // The first DoF for this moment
    const size_t uk_map = ell * n_groups;

    // Loop over groups
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double rhs = 0.0; // increment into this to avoid many accessor calls

      // Inhomogeneous isotropic source
      if (apply_mat_src && ell == 0 && src)
        rhs += src->values[g];

      // Scattering term
      if (apply_scatter_src)
      {
        const auto* sig_ell = xs->transfer_matrices[ell][g].data();
        for (unsigned int gp = 0; gp < n_groups; ++gp)
          rhs += density * *sig_ell++ * phi[uk_map + gp];
      }

      // Fission term
      if (apply_fission_src)
      {
        const auto chi = xs->chi[g];
        const auto* nusig_f = xs->nu_sigma_f.data();
        for (unsigned int gp = 0; gp < n_groups; ++gp)
          rhs += chi * density * *nusig_f++ * phi[uk_map + gp];
      }
      q_moments[uk_map + g] += rhs;
    }//for group
  }//for moment
}


void
SteadyStateSolver::dsa()
{
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    double dphi = phi[g] - phi_ell[g];
    phi[g] += xs->sigma_s[g] / xs->sigma_a[g] * dphi;
  }
}


void
SteadyStateSolver::update_neutron_density()
{
  density = 0.0;
  for (unsigned int g = 0; g < n_groups; ++g)
    density += phi[g] * xs->inv_velocity[g];
}


double
SteadyStateSolver::check_balance()
{
  double balance = 0.0;
  for (unsigned int g = 0; g < n_groups; ++g)
    balance += src->values[g] - xs->sigma_a[g] * phi[g];
  return balance;
}


const Vector&
SteadyStateSolver::get_angular_flux() const
{
  return psi;
}


const Vector&
SteadyStateSolver::get_flux_moments() const
{
  return phi;
}

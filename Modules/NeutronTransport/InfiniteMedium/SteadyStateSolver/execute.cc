#include "steadystate_solver.h"

#include <cassert>
#include <iomanip>


using namespace PDEs;
using namespace Math;
using namespace NeutronTransport;
using namespace InfiniteMedium;


void
SteadyStateSolver::execute()
{
  std::cout
    << "\n**********************************************************\n"
    << "Executing the multi-group infinite medium transport solver"
    << "\n**********************************************************\n";

  auto result = source_iterations();
  bool converged = result.second < inner_tolerance;

  if (verbosity > 0)
    std::cout
      << (converged?
          "\n***** Steady State Solver CONVERGED *****\n" :
          "\n!!*!! WARNING: Steady State Solver NOT CONVERGED *****\n")
      << "# of Iterations:  "
      << std::left << std::setw(6) << result.first << std::endl
      << "Final Phi Change  "
      << std::left << std::setw(6) << result.second << std::endl
      << std::endl;
}


std::pair<unsigned int, double>
SteadyStateSolver::source_iterations()
{
  phi_ell = phi;
  auto q_moments_init = q_moments;

  unsigned int nit;
  double change;
  for (nit = 0; nit < max_inner_iterations; ++nit)
  {
    q_moments = q_moments_init;
    set_source();
    sweep();
    if (use_dsa)
      dsa();

    change = l1_norm(phi - phi_ell);
    double balance = check_balance();
    phi_ell = phi;

    if (verbosity > 1)
      std::cout
        << std::left << "inner::"
        << "Iteration  " << std::setw(3) << nit << "   "
        << "Value  " << std::setw(10) << change << "  "
        << "Balance  " << std::setw(10) << balance << "  "
        << (change < inner_tolerance? "CONVERGED" : "")
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
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      // Solve for the angular flux
      double x = 0.0;
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        x += moment_to_discrete[ell][n]*q_moments[ell*n_groups + g];
      x /= xs->sigma_t[g];

      // Store the angular flux
      psi[n*n_groups + g] = x;

      // Accumulate flux moment
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        phi[ell*n_groups + g] += discrete_to_moment[ell][n]*x;
    }//for g
  }//for n
}


void
SteadyStateSolver::set_source()
{
  // Loop over moments
  const auto* S_ptr = xs->transfer_matrices.data();
  for (unsigned int ell = 0; ell < n_moments; ++ell, ++S_ptr)
  {
    // The first DoF for this moment
    const size_t uk_map = ell*n_groups;

    // Loop over groups
    const auto* Sg_ptr = S_ptr->data();
    for (unsigned int g = 0; g < n_groups; ++g, ++Sg_ptr)
    {
      double rhs = 0.0; // increment into this to avoid many accessor calls

      // Inhomogeneous isotropic source
      if (ell == 0 && src)
        rhs += src->values[g];

      // Scattering term
      const double* sig_ell = Sg_ptr->data();
      for (unsigned int gp = 0; gp < n_groups; ++gp, ++sig_ell)
        rhs += *sig_ell * phi[uk_map + gp];

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
    phi[g] += xs->sigma_s[g]/xs->sigma_a[g]*dphi;
  }
}


double
SteadyStateSolver::check_balance()
{
  double balance = 0.0;
  for (unsigned int g = 0; g < n_groups; ++g)
    balance += src->values[g] - xs->sigma_a[g]*phi[g];
  return balance;
}

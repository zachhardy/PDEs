#include "infinite_medium.h"

#include "Math/Polynomials/polynomial.h"

#include <numeric>
#include <cassert>
#include <iomanip>


using namespace PDEs;
using namespace Math;
using namespace Polynomials;
using namespace NeutronTransport;


void
InfiniteMedium::initialize()
{
  assert(quadrature->size() % 2  == 0);
  assert(src == nullptr || src->values.size() == xs->n_groups);

  std::cout
      << "\n*********************************************************\n"
      <<   "Initializing infinite medium multi-group transport solver"
      << "\n*********************************************************\n";

  n_groups = xs->n_groups;
  n_moments = xs->scattering_order + 1;
  n_angles = quadrature->size();

  psi.resize(n_groups * n_angles);
  phi.resize(n_groups * n_moments);
  q_moments.resize(n_groups * n_moments);

  compute_moment_to_discrete_operator();
  compute_discrete_to_moment_operator();

  std::cout
      << "\n------------------------------\n"
      <<   "--- Simulation Information ---"
      << "\n------------------------------\n"
      << "Groups:                 " << n_groups << std::endl
      << "Moments:                " << n_moments << std::endl
      << "Angles:                 " << n_angles << std::endl
      << "Angular flux unknowns:  " << psi.size() << std::endl
      << "Flux moment unknowns:   " << phi.size() << std::endl;
}


void
InfiniteMedium::execute()
{
  std::cout
      << "\n**********************************************************\n"
      <<   "Executing the multi-group infinite medium transport solver"
      << "\n**********************************************************\n";

  auto phi_ell = phi;
  for (unsigned int nit = 0; nit < max_iterations; ++nit)
  {
    set_source();
    solve();

    double change = l1_norm(phi - phi_ell);
    phi_ell = phi;

    std::stringstream iter_info;
    iter_info
        << std::left << "inner::"
        << "Iteration  " << std::setw(3) << nit << "   "
        << "Value  " << change;
    if (change < tolerance) iter_info << "  CONVERGED";
    std::cout << iter_info.str() << std::endl;

    if (change < tolerance)
      break;
  }
}


void
InfiniteMedium::solve()
{
  phi = 0.0;
  for (unsigned int n = 0; n < n_angles; ++n)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      // Solve for the angular flux
      double x = 0.0;
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        x += moment_to_discrete[ell][n] * q_moments[ell*n_groups + g];
      x /= xs->sigma_t[g];

      // Store the angular flux
      psi[n*n_groups + g] = x;

      // Accumulate flux moment
      for (unsigned int ell = 0; ell < n_moments; ++ell)
        phi[ell*n_groups + g] += discrete_to_moment[ell][n] * x;
    }//for g
  }//for n
}


void
InfiniteMedium::set_source()
{
  // Reset the source moments vector to zero
  q_moments = 0.0;

  // Loop over moments
  const auto* S_ptr = xs->transfer_matrices.data();
  for (unsigned int ell = 0; ell < n_moments; ++ell, ++S_ptr)
  {
    // The first DoF for this moment
    const size_t uk_map = ell * n_groups;

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
InfiniteMedium::compute_moment_to_discrete_operator()
{
  moment_to_discrete.reinit(n_moments, n_angles);

  const auto& wgts = quadrature->get_weights();
  const auto& qpoints = quadrature->get_quadrature_points();
  const double norm = std::accumulate(wgts.begin(), wgts.end(), 0.0);

  for (unsigned int ell = 0; ell < n_moments; ++ell)
  {
    double f = (2.0*ell + 1)/norm;
    for (unsigned int n = 0; n < n_angles; ++n)
      moment_to_discrete[ell][n] = f * legendre(ell, qpoints[n].z());
  }

//  std::cout << "\nMoment to Discrete Operator:\n";
//  for (unsigned int n = 0; n < n_angles; ++n)
//  {
//    std::cout << std::left << std::setw(5) << n;
//    for (unsigned int ell = 0; ell < n_moments; ++ell)
//      std::cout
//          << std::setw(15) << std::left << std::fixed
//          << std::setprecision(10) << moment_to_discrete[ell][n]
//          << (ell + 1 < n_moments? " " : "\n");
//  }
}


void
InfiniteMedium::compute_discrete_to_moment_operator()
{
  discrete_to_moment.reinit(n_moments, n_angles);

  const auto& wgts = quadrature->get_weights();
  const auto& qpoints = quadrature->get_quadrature_points();

  for (unsigned int ell = 0; ell < n_moments; ++ell)
    for (unsigned int n = 0; n < n_angles; ++n)
      discrete_to_moment[ell][n] = wgts[n] * legendre(ell, qpoints[n].z());

//  std::cout << "\nDiscrete to Moment Operator:\n";
//  for (unsigned int n = 0; n < n_angles; ++n)
//  {
//    std::cout << std::left << std::setw(5) << n;
//    for (unsigned int ell = 0; ell < n_moments; ++ell)
//      std::cout
//          << std::setw(15) << std::left << std::fixed
//          << std::setprecision(10) << discrete_to_moment[ell][n]
//          << (ell + 1 < n_moments ? " " : "\n");
//  }
}


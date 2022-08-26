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

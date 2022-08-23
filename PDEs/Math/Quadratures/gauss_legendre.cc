#include "gauss_legendre.h"
#include "Polynomials/polynomial.h"

#include <cmath>
#include <limits>
#include <algorithm>
#include <cassert>


using namespace PDEs;
using namespace Math;
using namespace Polynomials;


GaussLegendreQuadrature::
GaussLegendreQuadrature(const unsigned int n, const bool verbose) :
    Quadrature(n)
{
  // Compute the quadrature points
  auto qpoints = this->legendre_roots(n);
  for (const auto& qpoint : qpoints)
    quadrature_points.emplace_back(0.0, 0.0, qpoint);

  // Compute the weights
  for (unsigned int q = 0; q < qpoints.size(); ++q)
  {
    const double& qp = qpoints[q];
    double dPn = dlegendre(n, qpoints[q]);
    weights.push_back(2.0 / (1.0 - qp*qp) / (dPn*dPn));
  }

  std::cout << "Quadrature Points: [";
  for (unsigned int i = 0; i < n; ++i)
    std::cout << qpoints[i] << (i + 1 < n ? "  " : "]\n");

  std::cout << "Quadrature Weights:[";
  for (unsigned int i = 0; i < n; ++i)
    std::cout << weights[i] << (i + 1 < n ? "  " : "]\n");
}


std::vector<double>
GaussLegendreQuadrature::legendre_roots(const unsigned int n)
{
  assert(n > 0);

  // Define the tolerance
  const double tol = 4.0 * std::numeric_limits<double>::epsilon();

  // Approximate the zeros by detecting sign changes in the Legendre
  // polynomial. This is done by using a partitioning of the domain into
  // search intervals.
  unsigned int n_bins = (n < 128) ? 1000 : 10000;
  const double delta = 2.0 / n_bins;  //bin width

  size_t count = 0;
  std::vector<double> roots(n, 0.0);
  for (size_t i = 0; i < n_bins; ++i)
  {
    const double x_i = -1.0 + i * delta;
    const double x_ip = std::min(x_i + delta, 1.0);
    if (legendre(n, x_i) * legendre(n, x_ip) < 0)
      roots[count++] = 0.5 * (x_i + x_ip);
  }

  // For each initial guess, apply Newton's method
  for (unsigned int k = 0; k < n; ++k)
  {
    for (unsigned int nit = 0; nit < 1000; ++nit)
    {
      double x0 = roots[k];
      roots[k] = x0 - legendre(n, x0) / dlegendre(n, x0);
      if (std::fabs(roots[k] - x0) < tol)
        break;
    }
  }
  std::stable_sort(roots.begin(), roots.end());
  return roots;
}

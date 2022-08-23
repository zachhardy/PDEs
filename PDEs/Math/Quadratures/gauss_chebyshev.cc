#include "gauss_chebyshev.h"
#include "Polynomials/polynomial.h"

#include <cmath>

#include <iomanip>
#include <cassert>


using namespace PDEs;
using namespace Math;
using namespace Polynomials;


GaussChebyshevQuadrature::
GaussChebyshevQuadrature(const unsigned int n, const bool verbose) :
    Quadrature(n)
{
  if (verbose)
    std::cout << "\nInitializing Gauss-Chebyshev quadratures...\n"
              << "Number of Points:\t" << n << std::endl
              << "---------------------------" << std::endl
              << std::setw(3) << "n"
              << std::setw(12) << "Point"
              << std::setw(12) << "Weight" << std::endl
              << "---------------------------" << std::endl;

  // Compute the quadrature points
  auto qpoints = this->chebyshev_roots(n);
  for (const auto& qpoint: qpoints)
    quadrature_points.emplace_back(0.0, 0.0, qpoint);

  // Compute the weights
  weights.assign(n, M_PI / n);

  if (verbose)
    for (unsigned int i = 0; i < n; ++i)
      std::cout << std::setw(3) << i
                << std::setw(12) << qpoints[i]
                << std::setw(12) << weights[i] << std::endl;
}


std::vector<double>
GaussChebyshevQuadrature::chebyshev_roots(const unsigned int n)
{
  std::vector<double> roots(n, 0.0);
  for (unsigned int i = 0; i < n; ++i)
    roots[i] = -std::cos((2.0 * i + 1.0) / (2.0 * n) * M_PI);
  return roots;
}

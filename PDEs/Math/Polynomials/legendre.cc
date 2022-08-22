#include "legendre.h"

#include <cmath>
#include <limits>
#include <algorithm>
#include <cassert>


double
Math::legendre(const unsigned int n, const double x)
{
  if (n == 0) return 1.0;
  if (n == 1) return x;

  double P_m = 1.0, P = x, P_p;
  for (unsigned int p = 2; p <= n; ++p)
  {
    P_p = (2.0*p - 1.0)/p * x * P - (p - 1.0)/p * P_m;
    P_m = P;
    P = P_p;
  }
  return P_p;
}


double
Math::dlegendre(const unsigned int n, const double x)
{
  if (n == 0) return 0.0;
  if (n == 1) return 1.0;

  // Invalid at x = -1 and x = 1, add or subtract epsilon
  double x_s = x;
  const double e_mach = std::numeric_limits<double>::epsilon();
  if (x == -1.0) x_s += 4.0 * e_mach;
  if (x == 1.0)  x_s -= 4.0 * e_mach;

  return n/(x_s*x_s - 1.0) * (x_s*legendre(n, x_s) - legendre(n - 1, x_s));
}


double
Math::d2legendre(const unsigned int n, const double x)
{
  if (n == 0) return 0.0;
  if (n == 1) return 0.0;

  // Invalid at x = -1 and x = 1, add or subtract epsilon
  double x_s = x;
  if (x == -1.0) x_s += 4.0*std::numeric_limits<double>::epsilon();
  if (x == 1.0)  x_s -= 4.0*std::numeric_limits<double>::epsilon();

  // Compute an epsilon for the numerical derivative
  const double e_mach = std::numeric_limits<double>::epsilon();
  const double eps = (1.0 + x) * std::sqrt(e_mach);

  // Define evaluation points
  const double x_m = std::min(x - eps, -1.0);
  const double x_p = std::max(x + eps, 1.0);

  // Compute the numerical derivative
  const double v_p = dlegendre(n, x_m);
  const double v_m = dlegendre(n, x_p);
  return (v_p - v_m) / (2.0 * eps);
}


std::vector<double>
Math::legendre_roots(const unsigned int n)
{
  assert(n > 0);

  // Define the tolerance
  const double tol = 4.0 * std::numeric_limits<double>::epsilon();

  // Approximate the zeros by detecting sign changes in the Legendre
  // polynomial. This is done by using a partitioning of the domain into
  // search intervals.
  unsigned int n_bins = (n < 128)? 1000 : 10000;
  const double delta = 2.0/n_bins;  //bin width

  size_t count = 0;
  std::vector<double> roots(n, 0.0);
  for (size_t i = 0; i < n_bins; ++i)
  {
    const double x_i = -1.0 + i*delta;
    const double x_ip = std::min(x_i + delta, 1.0);
    if (legendre(n, x_i) * legendre(n, x_ip) < 0)
      roots[count++] = 0.5*(x_i + x_ip);
  }

  // For each initial guess, apply Newton's method
  for (unsigned int k = 0; k < n; ++k)
  {
    for (unsigned int nit = 0; nit < 1000; ++nit)
    {
      double x0 = roots[k];
      roots[k] = x0 - legendre(n, x0)/dlegendre(n, x0);
      if (std::fabs(roots[k] - x0) < tol)
        break;
    }
  }
  std::stable_sort(roots.begin(), roots.end());
  return roots;
}

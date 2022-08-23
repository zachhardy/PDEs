#include "polynomial.h"

#include <cmath>
#include <limits>
#include <algorithm>
#include <cassert>


using namespace PDEs;
using namespace Math;
using namespace Polynomials;


double
Polynomials::legendre(const unsigned int n, const double x)
{
  if (n == 0) return 1.0;
  if (n == 1) return x;

  double P_m = 1.0, P = x, P_p;
  for (unsigned int p = 2; p <= n; ++p)
  {
    P_p = (2.0 * p - 1.0) / p * x * P - (p - 1.0) / p * P_m;
    P_m = P;
    P = P_p;
  }
  return P_p;
}


double
Polynomials::dlegendre(const unsigned int n, const double x)
{
  if (n == 0) return 0.0;
  if (n == 1) return 1.0;

  // Invalid at x = -1 and x = 1, add or subtract epsilon
  double x_s = x;
  const double e_mach = std::numeric_limits<double>::epsilon();
  if (x == -1.0) x_s += 4.0 * e_mach;
  if (x == 1.0) x_s -= 4.0 * e_mach;

  return n / (x_s * x_s - 1.0) * (x_s * legendre(n, x_s) - legendre(n - 1, x_s));
}


double
Polynomials::d2legendre(const unsigned int n, const double x)
{
  if (n == 0) return 0.0;
  if (n == 1) return 0.0;

  // Invalid at x = -1 and x = 1, add or subtract epsilon
  double x_s = x;
  if (x == -1.0) x_s += 4.0 * std::numeric_limits<double>::epsilon();
  if (x == 1.0) x_s -= 4.0 * std::numeric_limits<double>::epsilon();

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


double
Polynomials::chebyshev(const unsigned int n, const double x)
{
  if (n == 0) return 1.0;
  if (n == 1) return x;

  double T_m = 1.0, T = x, T_p;
  for (unsigned int p = 2; p <= n; ++p)
  {
    T_p = 2.0 * x * T - T_m;
    T_m = T;
    T = T_p;
  }
  return T_p;
}

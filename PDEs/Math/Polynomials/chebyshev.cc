#include "chebyshev.h"

#include <cmath>
#include <limits>
#include <algorithm>
#include <cassert>


double
Math::chebyshev(const unsigned int n, const double x)
{
  if (n == 0) return 1.0;
  if (n == 1) return x;

  double T_m = 1.0, T = x, T_p;
  for (unsigned int p = 2; p <= n; ++p)
  {
    T_p = 2.0*x * T - T_m;
    T_m = T;
    T = T_p;
  }
  return T_p;
}


std::vector<double>
Math::chebyshev_roots(const unsigned int n)
{
  assert(n > 0);
  std::vector<double> roots(n);
  for (unsigned int k = 0; k < n; ++k)
    roots[k] = -std::cos(M_PI * (k + 0.5)/n);
  return roots;
}

#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <vector>


namespace Math
{
  /**
   * Evaluate the order \p n Chebyshev polynomial.
   *
   * \param n The Chebyshev polynomial order.
   * \param x The evaluation point. This must be in the range [-1, 1].
   */
  double chebyshev(const unsigned int n, const double x);

  std::vector<double> chebyshev_roots(const unsigned int n);
}

#endif //CHEBYSHEV_H

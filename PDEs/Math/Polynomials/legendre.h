#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <vector>


namespace Math
{
  /**
   * Evaluate the order \p n Legendre polynomial.
   *
   * \param n The Legendre polynomial order.
   * \param x The evaluation point. This must be in the range [-1, 1].
   */
  double legendre(const unsigned int n, const double x);

  /**
   * Evaluate the derivative of the order \p n Legendre polynomial.
   *
   * \param n The Legendre polynomial order.
   * \param x The evaluation point. This must be in the range [-1, 1].
   */
  double dlegendre(const unsigned int n, const double x);


  /**
   * Evaluate the second derivative of the order \p n Legendre polynomial.
   *
   * \param n The Legendre polynomial order.
   * \param x The evaluation point. This must be in the range [-1, 1].
   */
  double d2legendre(const unsigned int n, const double x);


  std::vector<double> legendre_roots(const unsigned int n);
}

#endif //LEGENDRE_H

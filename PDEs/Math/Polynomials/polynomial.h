#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <vector>


namespace PDEs
{
  namespace Math
  {
    namespace Polynomials
    {
      /**
       * Evaluate the order \p n Legendre polynomial.
       *
       * \param n The Legendre polynomial order.
       * \param x The evaluation point. This must be in the range [-1, 1].
       */
      double
      legendre(const unsigned int n, const double x);

      /**
       * Evaluate the derivative of the order \p n Legendre polynomial.
       *
       * \param n The Legendre polynomial order.
       * \param x The evaluation point. This must be in the range [-1, 1].
       */
      double
      dlegendre(const unsigned int n, const double x);


      /**
       * Evaluate the second derivative of the order \p n Legendre polynomial.
       *
       * \param n The Legendre polynomial order.
       * \param x The evaluation point. This must be in the range [-1, 1].
       */
      double
      d2legendre(const unsigned int n, const double x);

      /**
       * Return the roots of the Legendre polynomial of degree \p n. This is
       * used to define the quadrature points for Gauss-Legendre quadratures.
       *
       * \param n The Legendre polynomial order.
       */
      std::vector<double>
      legendre_roots(const unsigned int n);

      /**
       * Evaluate the order \p n Chebyshev polynomial.
       *
       * \param n The Chebyshev polynomial order.
       * \param x The evaluation point. This must be in the range [-1, 1].
       */
      double
      chebyshev(const unsigned int n, const double x);

      /**
       * Return the roots of the Chebyshev polynomial of degree \p n. This is
       * used to define the quadrature points for Gauss-Chebyshev quadratures.
       *
       * \param n The Chebyshev polynomial order.
       */
      std::vector<double>
      chebyshev_roots(const unsigned int n);
    }
  }
}

#endif //CHEBYSHEV_H

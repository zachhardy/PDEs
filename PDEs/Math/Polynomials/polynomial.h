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
       * Evaluate the order \p n Legendre polynomial at \p x.
       *
       * \note \p x should be in the range [-1, 1].
       */
      double legendre(const unsigned int n, const double x);

      /**
       * Evaluate the derivative of the order \p n Legendre polynomial at \p x.
       *
       * \note \p x should be in the range [-1, 1].
       */
      double dlegendre(const unsigned int n, const double x);


      /**
       * Evaluate the second derivative of the order \p n Legendre polynomial
       * at \p x.
       *
       * \note \p x should be in the range [-1, 1].
       */
      double d2legendre(const unsigned int n, const double x);

      /**
       * Evaluate the order \p n Chebyshev polynomial at \p x.
       *
       * \note \p x should be in the range [-1, 1].
       */
      double chebyshev(const unsigned int n, const double x);
    }
  }
}

#endif //CHEBYSHEV_H

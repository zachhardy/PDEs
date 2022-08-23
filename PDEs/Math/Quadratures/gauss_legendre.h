#ifndef PDES_GAUSS_LEGENDRE_H
#define PDES_GAUSS_LEGENDRE_H

#include "quadrature.h"


namespace PDEs
{
  namespace Math
  {
    /**
     * Implementation of Gauss-Legendre 1D quadratures. See more at
     * <a href="https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature">
     * https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature</a>.
     */
    class GaussLegendreQuadrature : public Quadrature
    {
    public:
      /**
       * Construct an \p n point Gauss-Legendre quadrature set. This routine
       * initializes the quadrature points using a root-finding algorithm for
       * the degree \p n Legendre polynomial, then initializes the associated
       * quadrature weights.
       */
      GaussLegendreQuadrature(const unsigned int n,
                              const bool verbose = false);

      /**
       * Return the roots of the degree \p n Legendre polynomial.
       */
      static std::vector<double>
      legendre_roots(const unsigned int n);

    };
  }// Math
}// PDEs


#endif //PDES_GAUSS_LEGENDRE_H

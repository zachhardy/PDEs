#ifndef PDES_GAUSS_LEGENDRE_H
#define PDES_GAUSS_LEGENDRE_H

#include "quadrature.h"


namespace PDEs
{
  namespace Math
  {
    /**
     * Implementation of Gauss-Legendre 1D quadratures.
     *
     * Gauss-Legendre quadratures points are defined by the zeros of the
     * Legendre polynomials. At construction for an \p n point quadrature
     * set, the degree \p n Legendre polynomial zeros are computed using a
     * root-finding algorithm along with the associate quadrature weights.
     *
     * See more at
     * https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature.
     */
    class GaussLegendreQuadrature : public Quadrature
    {
    public:
      /** Construct an \p n point Gauss-Legendre quadrature set. */
      GaussLegendreQuadrature(const unsigned int n,
                              const bool verbose = false);

      /** Return the roots of the degree \p n Legendre polynomial. */
      static std::vector<double> legendre_roots(const unsigned int n);

    };
  }// Math
}// PDEs


#endif //PDES_GAUSS_LEGENDRE_H

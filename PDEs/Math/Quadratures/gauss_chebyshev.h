#ifndef GAUSS_CHEBYSHEV_H
#define GAUSS_CHEBYSHEV_H

#include "quadrature.h"


namespace PDEs
{
  namespace Math
  {

    /**
     * Implementation of Gauss-Chebyshev 1D quadratures.
     *
     * Gauss-Chebyshev quadratures points are defined by the zeros of the
     * Chebyshev polynomials. The zeros of these polynomials have an analytic
     * form and the quadrature weights uniform. These are populated at
     * construction.
     *
     * See more at
     * https://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature.
     */
    class GaussChebyshevQuadrature : public Quadrature
    {
    public:
      /**
       * Construct an \p n point Gauss-Chebyshev quadrature set. This routine
       * initializes the quadrature points with the analytic roots of the
       * Chebyshev polynomials, then initializes the associated quadrature
       * weights.
       */
      GaussChebyshevQuadrature(const unsigned int n,
                               const bool verbose = false);

      /**
       * Return the roots of the degree \p n Chebyshev polynomial.
       */
      static std::vector<double>
      chebyshev_roots(const unsigned int n);
    };

  } // Math
} // PDEs

#endif //GAUSS_CHEBYSHEV_H

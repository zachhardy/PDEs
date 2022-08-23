#ifndef PDES_QUADRATURE_H
#define PDES_QUADRATURE_H

#include "Grid/cartesian_vector.h"

#include <vector>


namespace PDEs
{
  using namespace Grid;


  namespace Math
  {
    /**
     * A base class for quadratures.
     */
    class Quadrature
    {
    public:
      /**
       * Construct a quadrature set with \p n quadrature points.
       */
      explicit Quadrature(const unsigned int n);

    protected:
      /**
       * The number of quadrature points.
       */
      const unsigned int n_quadrature_points;

      /**
       * The quadrature points in the quadrature sets. Note that for 1D
       * quadratures the points reside in the z-coordinate.
       */
      std::vector<Point> quadrature_points;

      /**
       * The associated quadrature weights.
       */
      std::vector<double> weights;
    };

  } // PDEs
} // Math

#endif //PDES_QUADRATURE_H

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "Grid/cartesian_vector.h"

#include <vector>


namespace PDEs
{
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

      /**
       * Return the size of the quadrature set.
       */
      unsigned int
      size() const;

      /**
       * Return a the <tt>i</tt>'th quadrature point.
       */
      const Grid::Point&
      quadrature_point(const unsigned int i) const;

      /**
       * Return a the <tt>i</tt>'th quadrature weight.
       */
      const double&
      weight(const unsigned int i) const;

      /**
       * Return the set of quadrature points.
       */
      const std::vector<Grid::Point>&
      get_quadrature_points() const;

      /**
       * Return the set of quadrature weights.
       */
      const std::vector<double>&
      get_weights() const;

    protected:
      const unsigned int n_quadrature_points;

      /**
       * The quadrature points in the quadrature sets. Note that for 1D
       * quadratures the points reside in the z-coordinate.
       */
      std::vector<Grid::Point> quadrature_points;

      /**
       * The associated quadrature weights.
       */
      std::vector<double> weights;
    };
  } // PDEs
} // Math

#endif //QUADRATURE_H

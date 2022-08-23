#include "quadrature.h"


using namespace PDEs;
using namespace Math;


Quadrature::Quadrature(const unsigned int n) :
    n_quadrature_points(n)
{}


unsigned int
Quadrature::size() const
{
  return n_quadrature_points;
}


const Grid::Point&
Quadrature::quadrature_point(const unsigned int i) const
{
  return quadrature_points[i];
}


const double&
Quadrature::weight(const unsigned int i) const
{
  return weights[i];
}


const std::vector<Grid::Point>&
Quadrature::get_quadrature_points() const
{
  return quadrature_points;
}


const std::vector<double>&
Quadrature::get_weights() const
{
  return weights;
}




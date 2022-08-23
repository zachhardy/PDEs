#include "quadrature.h"


using namespace PDEs;
using namespace Math;


Quadrature::Quadrature(const unsigned int n) :
  n_quadrature_points(n), quadrature_points(n), weights(n)
{}

#include "boundaries.h"

using namespace NeutronDiffusion;


Boundary::Boundary(BoundaryType type)
  : bndry_type(type)
{}


BoundaryType
Boundary::type() const
{ return bndry_type; }

//######################################################################

DirichletBoundary::DirichletBoundary()
  : Boundary(BoundaryType::DIRICHLET)
{}


DirichletBoundary::DirichletBoundary(const double value)
  : Boundary(BoundaryType::DIRICHLET), value(value)
{}

//######################################################################

NeumannBoundary::NeumannBoundary()
  : Boundary(BoundaryType::NEUMANN)
{}


NeumannBoundary::NeumannBoundary(const double value)
  : Boundary(BoundaryType::NEUMANN), value(value)
{}

//######################################################################

RobinBoundary::RobinBoundary()
  : Boundary(BoundaryType::ROBIN)
{}


RobinBoundary::RobinBoundary(const double j_inc)
  : Boundary(BoundaryType::ROBIN), f(j_inc)
{}


RobinBoundary::RobinBoundary(const double a, const double b, const double f)
  : Boundary(BoundaryType::ROBIN), a(a), b(b), f(f)
{}




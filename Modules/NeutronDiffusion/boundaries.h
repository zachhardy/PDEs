#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <sstream>
#include <vector>


namespace NeutronDiffusion
{

/// Types of diffusion boundary conditions.
enum class BoundaryType
{
  DIRICHLET   = 1,  ///< Dirichlet boundary.
  ZERO_FLUX   = 2,  ///< Homogeneous Dirichlet boundary.
  NEUMANN     = 3,  ///< Neumann boundary.
  REFLECTIVE  = 4,  ///< Homogeneous Neumann boundary.
  ROBIN       = 5,  ///< Robin boundary.
  VACUUM      = 6,  ///< Homogeneous Robin boundry with `a, b = 0.25, 0.5`.
  MARSHAK     = 7   ///< Robin boundary with `a, b = 0.25, 0.5`.
};

//######################################################################

/// Abstract base class for diffusion boundaries.
class Boundary
{
public:
  const BoundaryType type;

public:
  explicit Boundary(BoundaryType boundary_type) : type(boundary_type) {}
};

//######################################################################

/// Dirichlet boundary given by \f$ u_b = f^d \f$.
class DirichletBoundary : public Boundary
{
public:
  double value = 0.0;

public:
  /// Construct a zero flux boundary.
  DirichletBoundary() : Boundary(BoundaryType::DIRICHLET) {}

  explicit DirichletBoundary(double boundary_value)
    : Boundary(BoundaryType::DIRICHLET), value(boundary_value)
  {}
};

//######################################################################

/// Neumann boundary given by \f$ \partial_{\hat{n}_b} u = f^n \f$.
class NeumannBoundary : public Boundary
{
public:
  double value = 0.0;

public:
  /// Construct a reflective boundary.
  NeumannBoundary() : Boundary(BoundaryType::NEUMANN) {}

  explicit NeumannBoundary(double boundary_value)
    : Boundary(BoundaryType::NEUMANN), value(boundary_value)
  {}
};

//######################################################################

/// Robin boundary given by \f$ a u_b + b \partial_{\hat{n}_b} u = f^r \f$.
class RobinBoundary : public Boundary
{
public:
    double a = 0.25;
    double b = 0.5;
    double f = 0.0;

public:
  /// Contruct a vacuum boundary.
  RobinBoundary() : Boundary(BoundaryType::ROBIN) {}

  /// Construct a Marshak boundary from an incident partial current.
  explicit RobinBoundary(double j_inc)
    : Boundary(BoundaryType::ROBIN), f(j_inc)
  {}

  explicit RobinBoundary(double a_value, double b_value, double f_value)
      : Boundary(BoundaryType::ROBIN), a(a_value), b(b_value), f(f_value)
  {}
};

}
#endif //BOUNDARIES_H

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <sstream>
#include <vector>


namespace NeutronDiffusion
{

/** Types of diffusion boundary conditions. */
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

/** Abstract base class for diffusion boundaries. */
class Boundary
{
private:
  const BoundaryType bndry_type;

public:
  /** Default constructor. */
  explicit
  Boundary(BoundaryType type) : bndry_type(type) {}

  /** Get the boundary type. */
  BoundaryType
  type() const;
};

//######################################################################

/** Dirichlet boundary given by \f$ u_b = f^d \f$. */
class DirichletBoundary : public Boundary
{
public:
  double value = 0.0;

public:
  /** Construct a zero flux boundary. */
  DirichletBoundary();

  /** Construct a general Dirichlet boundary. */
  explicit
   DirichletBoundary(const double value);
};

//######################################################################

/** Neumann boundary given by \f$ \partial_{\hat{n}_b} u = f^n \f$. */
class NeumannBoundary : public Boundary
{
public:
  double value = 0.0;

public:
  /** Construct a reflective boundary. */
  NeumannBoundary();

  /** Construct a general Neumann boundary. */
  explicit
  NeumannBoundary(const double value)
    : Boundary(BoundaryType::NEUMANN), value(value)
  {}
};

//######################################################################

/** Robin boundary given by \f$ a u_b + b \partial_{\hat{n}_b} u = f^r \f$. */
class RobinBoundary : public Boundary
{
public:
    double a = 0.25;
    double b = 0.5;
    double f = 0.0;

public:
  /** Contruct a vacuum boundary. */
  RobinBoundary();

  /** Construct a Marshak boundary from an incident partial current. */
  explicit
  RobinBoundary(const double j_inc);

  /** Construct a general Robin boundary. */
  explicit
  RobinBoundary(const double a,  const double b, const double f);
};

}
#endif //BOUNDARIES_H

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


  /**
   * Abstract base class for diffusion boundaries.
   */
  class Boundary
  {
  private:
    const BoundaryType bndry_type;

  public:
    explicit Boundary(BoundaryType type);

    BoundaryType type() const;
  };

  //######################################################################

  /**
   * Dirichlet boundary given by \f$ u_b = f^d \f$.
   */
  class DirichletBoundary : public Boundary
  {
  public:
    double value;

  public:
    /** Construct a Dirichlet boundary condition. */
    explicit DirichletBoundary(const double value = 0.0);
  };

  //######################################################################

  /**
   * Neumann boundary given by \f$ \partial_{\hat{n}_b} u = f^n \f$.
   */
  class NeumannBoundary : public Boundary
  {
  public:
    double value;

  public:
    /** Construct a Neumann boundary condition. */
    explicit NeumannBoundary(const double value = 0.0);
  };

  //######################################################################

  /**
   * Robin boundary given by \f$ a u_b + b \nabla u \cdot \hat{n}_b  = f^r \f$.
   */
  class RobinBoundary : public Boundary
  {
  public:
    double a = 0.25; ///< The coefficient for the unknown
    double b = 0.5; ///< The coefficient for the gradient term.
    double f = 0.0; ///< The boundary value.

  public:
    /** Construct a vacuum boundary. */
    RobinBoundary();

    /** Construct a Marshak boundary from an incident partial current. */
    explicit RobinBoundary(const double j_inc);

    /** Construct a general Robin boundary. */
    RobinBoundary(const double a,  const double b, const double f);
  };

}
#endif //BOUNDARIES_H

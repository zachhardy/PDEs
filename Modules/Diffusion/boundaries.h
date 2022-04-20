#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "diffusion.h"

#include <sstream>
#include <vector>


namespace diffusion
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

}


//######################################################################
/// Abstract base class for diffusion boundaries.
class diffusion::Boundary
{
public:
  /// An identifier for the boundary type.
  const BoundaryType type;

public:
  /// Default constructor.
  explicit Boundary(BoundaryType boundary_type) : type(boundary_type) {}
};


//######################################################################
/**
 * \brief Derived Boundary class for Dirichlet boundaries.
 *
 * These boundaries are given by \f$ u_b = f^d \f$.
 */
class diffusion::DirichletBoundary : public Boundary
{
public:
  double value = 0.0; ///< The boundary value.

public:
  /// Construct a homogeneous Dirichlet boundary.
  DirichletBoundary() : Boundary(BoundaryType::DIRICHLET) {}

   /// Construct a Dirichlet boundary.
  explicit DirichletBoundary(double boundary_value)
    : Boundary(BoundaryType::DIRICHLET), value(boundary_value)
  {}
};


//######################################################################
/**
 * \brief Derived Boundary class for Neumann boundaries.
 *
 * These boundaries are given by \f$ \partial_{\hat{n}_b} u = f^n \f$.
 */
class diffusion::NeumannBoundary : public Boundary
{
public:
  double value = 0.0; ///< The boundary derivative.

public:
  /// Construct a reflective boundary.
  NeumannBoundary() : Boundary(BoundaryType::NEUMANN) {}

  /// Construct a Neumann boundary.
  explicit NeumannBoundary(double boundary_value)
    : Boundary(BoundaryType::NEUMANN), value(boundary_value)
  {}
};


/**
 * \brief Derived Boundary class for Robin boundaries.
 *
 * These boundaries are given by
 * \f$ a u_b + b \partial_{\hat{n}_b} u = f^r \f$.
 */
class diffusion::RobinBoundary : public Boundary
{
public:
    double a = 0.25; ///< The scaling factor for the boundary value.
    double b = 0.5;  ///< The scaling factor for the boundary derivative.
    double f = 0.0;  ///< The boundary source term.

public:
  /// Contruct a vacuum boundary.
  RobinBoundary() : Boundary(BoundaryType::ROBIN) {}

  /**
   * \brief Construct a Marshak boundary.
   * \param j_inc The incident partial current.
   */
  explicit RobinBoundary(double j_inc)
    : Boundary(BoundaryType::ROBIN), f(j_inc)
  {}

  /**
   * \brief Construct a fully specified Robin boundary
   * \param a_value The scaling factor for the boundary value.
   * \param b_value The scaling factor for the boundary derivative.
   * \param f_value The boundary source term.
   */
  explicit RobinBoundary(double a_value, double b_value, double f_value)
      : Boundary(BoundaryType::ROBIN), a(a_value), b(b_value), f(f_value)
  {}
};

#endif //BOUNDARIES_H

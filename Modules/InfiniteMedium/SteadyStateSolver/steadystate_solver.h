#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "Math/vector.h"
#include "Math/matrix.h"

#include "Math/Quadratures/gauss_legendre.h"

#include "Physics/material.h"
#include "CrossSections/cross_sections.h"


namespace InfiniteMedium
{
    using namespace PDEs;
    using namespace Math;
    using namespace Physics;

  /**
   * Bitwise source flags for right-hand side vector construction.
   */
  enum SourceFlags : int
  {
    NO_SOURCE_FLAGS = 0,
    APPLY_MATERIAL_SOURCE = (1 << 0),
    APPLY_SCATTER_SOURCE = (1 << 1),
    APPLY_FISSION_SOURCE = (1 << 2),
  };


  inline SourceFlags
  operator|(const SourceFlags f1, const SourceFlags f2)
  {
    return static_cast<SourceFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
  }

  //######################################################################

  /**
   * An infinite medium multi-group discrete ordinates solver.
   */
  class SteadyStateSolver
  {
  public:

    bool use_dsa = false;

    double inner_tolerance = 1.0e-6;
    unsigned int max_inner_iterations = 100;

    unsigned int verbosity = 0;

    std::shared_ptr<CrossSections> xs;
    std::shared_ptr<GaussLegendreQuadrature> quadrature;
    std::shared_ptr<IsotropicMultiGroupSource> src;

  protected:

    unsigned int n_groups;
    unsigned int n_moments;
    unsigned int n_angles;

    /**
     * The multi-group angular flux vector.
     *
     * This vector stores groups contiguously. This means that for each angle,
     * all groups for that angle are contiguous.
     */
    Vector psi;

    /**
     * The multi-group flux moments vector.
     *
     * This vector stores groups contiguously. This means that for each angle,
     * all groups for that angle are contiguous.
     */
    Vector phi;
    Vector phi_ell; ///< The multi-group flux moments last iteration.
    Vector q_moments; ///< The multi-group source moments.

    Matrix discrete_to_moment; ///< The mapping from angle to moments.
    Matrix moment_to_discrete; ///< The mapping from moment to angles.

  public:
    virtual void initialize();
    virtual void execute();

    const Vector& get_angular_flux() const;
    const Vector& get_flux_moments() const;

    void
    write_angular_flux(const std::string file_prefix,
                       const std::string output_directory = ".") const;

    void
    write_flux_moments(const std::string file_prefix,
                       const std::string output_directory = ".") const;


  protected:
    /** Apply source iterations on the sources within \p source_flags. */
    std::pair<unsigned int, double>
    source_iterations(SourceFlags source_flags);

    /**
     * Set the right-hand side source moments vector. This is an additive
     * routine which will only add the specified sources to the source
     * moments vector. Source options include the inhomogeneous source,
     * scattering source, and fission source.
     */
    void
    set_source(SourceFlags source_flags = NO_SOURCE_FLAGS);

    /**
     * Sweep over the energy angles and energy groups to solve for the
     * next angular flux and flux moment updates.
     */
    void sweep();

    /** Apply an infinite medium version of DSA for acceleration.*/
    void dsa();

    void compute_moment_to_discrete_operator();
    void compute_discrete_to_moment_operator();

    double check_balance();
  };
}

#endif //STEADYSTATE_SOLVER_H

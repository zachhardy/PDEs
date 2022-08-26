#ifndef STEADYSTATE_SOLVER_H
#define STEADYSTATE_SOLVER_H

#include "Math/vector.h"
#include "Math/matrix.h"

#include "Math/Quadratures/gauss_legendre.h"

#include "Physics/material.h"
#include "CrossSections/cross_sections.h"


namespace NeutronTransport
{
  namespace InfiniteMedium
  {
    using namespace PDEs;
    using namespace Math;
    using namespace Physics;

    /**
     * An infinite medium multi-group discrete ordinates solver.
     */
    class SteadyStateSolver
    {
    public:

      /**
       * A flag for using diffusion synthetic acceleration.
       */
      bool use_dsa = false;

      /**
       * Inner iteration convergence inner_tolerance.
       */
      double inner_tolerance = 1.0e-6;

      /**
       * Maximum number of inner iterations to attempt.
       */
      unsigned int max_inner_iterations = 100;

      /**
       * A flag for the level of screen output.
       */
      unsigned int verbosity = 0;

      /**
       * The material cross-sections.
       */
      std::shared_ptr<CrossSections> xs;

      /**
       * The angular quadrature set.
       */
      std::shared_ptr<GaussLegendreQuadrature> quadrature;

      /**
       * An optional isotropic multi-group source for the problem.
       */
      std::shared_ptr<IsotropicMultiGroupSource> src;

      /**
       * Initialize the infinite medium transport solver.
       */
      virtual void
      initialize();

      /**
       * Execute the infinite medium transport solver.
       */
      virtual void
      execute();

    protected:
      /**
       * Implementation of the source iterations routine.
       */
      std::pair<unsigned int, double>
      source_iterations();

      /**
       * Compute the source moments.
       */
      void
      set_source();

      /**
       * Solve the system for the next angular flux and flux moment update.
       */
      void
      sweep();

      /**
       * Implementation of an infinite medium DSA algorithm.
       */
      void
      dsa();

      /**
       * Compute the moment-to-discrete operator.
       */
      void
      compute_moment_to_discrete_operator();

      /**
       * Compute the discrete-to-moment operator.
       */
      void
      compute_discrete_to_moment_operator();


      /**
       * The number of energy groups. This is obtained from the cross-sections.
       */
      unsigned int n_groups;

      /**
       * The number of flux moments. This is obtained from the cross-sections and
       * is defined as <tt>scattering_order + 1</tt>.
       */
      unsigned int n_moments;


      /**
       * The number of discrete angles. This is derived from the quadrature set.
       */
      unsigned int n_angles;

      /**
       * The angular flux vector. This vector is group-contiguous, or \f$ \psi =
       * [\psi_{n_1, g_1}, \ldots, \psi_{n_1, g_G}, \ldots, \psi_{n_N, g_1},
       * \ldots, \psi_{n_N, g_G}] \f$.
       */
      Vector psi;

      /**
       * The flux moments vector. This vector is group-contiguous, or \f$ \phi =
       * [\phi_{m_1, g_1}, \ldots, \phi_{m_1, g_G}, \ldots, \phi_{n_N, g_1,
       * \ldots, \phi_{n_N, g_G}] \f$.
       */
      Vector phi;

      /**
       * The flux moments from last iteration.
       */
      Vector phi_ell;

      /**
       * The source moments vector. This is ordered the same as #phi.
       */
      Vector q_moments;

      /**
       * The moment-to-discrete operator.
       */
      Matrix discrete_to_moment;

      /**
       * The discrete-to-moment operator.
       */
      Matrix moment_to_discrete;

      /**
       * Compute the particle balance.
       */
      double
      check_balance();
    };
  }
}

#endif //STEADYSTATE_SOLVER_H

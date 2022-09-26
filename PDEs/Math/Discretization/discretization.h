#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "mesh.h"

#include <cinttypes>


namespace PDEs
{
  namespace Math
  {

    /**
     * Available spatial discretization methods.
     */
    enum class SpatialDiscretizationMethod
    {
      FINITE_VOLUME = 1,  ///< Finite volume (FV)
      PIECEWISE_LINEAR_CONTINUOUS = 2,  ///< Linear continuous (PWLC)
      PIECEWISE_LINEAR_DISCONTINIOUS = 3,  ///< Linear dicontinuous (PWLD)
      LAGRANGE_CONTINUOUS = 4,  ///< Continuous finite elements (CFEM)
      LAGRANGE_DISCONTINUOUS = 5   ///< Discontinuous finite elements (DFEM)
    };

    //######################################################################

    /**
     * Abstracted base class for spatial discretizations.
     *
     * A discretization is built upon a mesh. Derived classes are meant to
     * contain all members and routines that are necessary to define the discrete
     * representation of a solution and the operations necessary for the
     * construction a linear system to solve.
     */
    class Discretization
    {
    public:
      const std::shared_ptr<Grid::Mesh> mesh;
      const SpatialDiscretizationMethod method;

    public:
      Discretization(const std::shared_ptr<Grid::Mesh> reference_mesh,
                     const SpatialDiscretizationMethod discretization_method) :
          mesh(reference_mesh), method(discretization_method)
      {}

      /**Return the number of nodes in the discretization. */
      virtual size_t n_nodes() const = 0;

      /** Return the number of nodes per cell in the discretization. */
      virtual unsigned int nodes_per_cell() const = 0;

      /**
       * Return the number of degrees of freedom (DoFs) in the discretization.
       * The number of DoFs is defined as the number of nodes multiplied by the
       * number of solution components.
       */
      virtual size_t n_dofs(const unsigned int n_components) const = 0;

      /**
       * Return the number of degrees of freedom (DoFs) per cell in the
       * discretization. The number of DoFs per cell is defined as the number of
       * nodes per cell multiplied by the number of solution components.
       */
      virtual unsigned int
      dofs_per_cell(const unsigned int n_components) const = 0;

      /** Return the coordinates of the nodes on the specified \p cell. */
      virtual std::vector<Grid::Node>
      nodes(const Grid::Cell& cell) const = 0;

      /**
       * Define the sparsity pattern. This routine defines the column indices of
       * non-zero entries per row for a problem with the specified number of
       * components. If the \p is_coupled flag is set to \p true, it is assumed
       * that all components are coupled to one another, otherwise, it is
       * assumed that the system is uncoupled in all components. The resulting
       * sparsity pattern is written into \p pattern.
       */
      virtual void
      make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                            const unsigned int n_components = 1,
                            const bool is_coupled = false) const = 0;

      /** Write the discretization to a file. */
      void write(const std::string directory,
                 const std::string file_prefix = "geom") const;
    };

  }
}
#endif //SPATIAL_DISCRETIZATION_H

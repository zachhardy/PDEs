#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "mesh.h"

#include <cinttypes>


namespace Math
{

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
   * A discretization is built upon a Mesh object. Derived classes
   * are meant to contain all members and routines that are necessary to define
   * the discrete representation of a solution and the operations necessary for
   * aiding in the construction a linear system to solve.
   */
  class Discretization
  {
  public:
    const std::shared_ptr<Grid::Mesh> mesh;
    const SpatialDiscretizationMethod type;

    explicit
    Discretization(const std::shared_ptr<Grid::Mesh> reference_mesh,
                   const SpatialDiscretizationMethod discretization_type)
      : mesh(reference_mesh), type(discretization_type)
    {}

    virtual size_t n_nodes() const = 0;
    virtual unsigned int nodes_per_cell() const = 0;

    virtual size_t n_dofs(const unsigned int n_components) const = 0;

    virtual unsigned int
    dofs_per_cell(const unsigned int n_components) const = 0;

    /** Get the location of the nodes on a cell. */
    virtual std::vector<Grid::CartesianVector> nodes(const Grid::Cell& cell) const = 0;

    /**
     * Define the sparsity pattern. This routine defines the column indices of
     * non-zero entries per row for a problem with the specified number of
     * components. If the \p is_coupled flag is set to \p true, it is assumed
     * that all components are coupled to one another, otherwise, it is assumed
     * that the system is uncoupled in all components.
     *
     * \param[out] pattern The column indices per row to allocate for.
     * \param n_components The number of components in the solution.
     * \param is_coupled A flag for allocating storage for coupling between
     *   solution components.
     */
    virtual void
    make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                          const unsigned int n_components = 1,
                          const bool is_coupled = false) const = 0;
  };

}
#endif //SPATIAL_DISCRETIZATION_H

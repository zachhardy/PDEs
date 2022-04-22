#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "Grid/grid_structs.h"
#include "Grid/mesh.h"

#include <memory>

namespace discretization
{

enum class DiscretizationMethod
{
  FINITE_VOLUME = 1,  ///< Finite volume (FV)
  PIECEWISE_LINEAR_CONTINUOUS = 2,  ///< Linear continuous (PWLC)
  PIECEWISE_LINEAR_DISCONTINIOUS = 3,  ///< Linear dicontinuous (PWLD)
  LAGRANGE_CONTINUOUS = 4,  ///< Continuous finite elements (CFEM)
  LAGRANGE_DISCONTINUOUS = 5   ///< Discontinuous finite elements (DFEM)
};

//######################################################################

/**
 * \brief Abstracted base class for spatial discretizations.
 *
 * A SpatialDiscretization is built upon a Mesh object. Derived classes
 * are meant to contain all members and routines that are necessary to define
 * the discrete representation of a solution and the operations necessary for
 * aiding in the construction a linear system to solve.
 */
class SpatialDiscretization
{
public:
  const std::shared_ptr<grid::Mesh> mesh;
  const DiscretizationMethod type;

public:
  explicit
  SpatialDiscretization(const std::shared_ptr<grid::Mesh> reference_mesh,
                        const DiscretizationMethod discretization_type)
      : mesh(reference_mesh), type(discretization_type)
  {}

public:
  /**
   * \brief Get the number of nodes in the discretization.
   *
   * This method should be overridden in derived classses.
   */
  virtual size_t n_nodes() const
  { return 0; }

  /**
   * \brief Get the number of DoFs in the spatial discretization.
   *
   * This method should be overridden in derived classes.
   */
  virtual size_t n_dofs(const size_t n_components) const
  { return 0; }

  /**
   * \brief Get the number of nodes per cell.
   *
   * This method should be overridden in derived classes.
   */
  virtual size_t nodes_per_cell() const
  { return 0; }

  /**
   * \brief Get the number of DoFs per cell.
   *
   * This method should be overriden in derived classes.
   */
  virtual size_t dofs_per_cell(const size_t n_components) const
  { return 0; }

  /**
   * \brief Get the location of the nodes across the spatial domain.
   *
   * This method should be overridden in derived classes.
   */
  virtual std::vector<grid::Point> nodes() const = 0;

};

}
#endif //SPATIAL_DISCRETIZATION_H

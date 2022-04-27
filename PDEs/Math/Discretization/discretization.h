#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "Grid/grid_structs.h"
#include "Grid/mesh.h"

#include <memory>

namespace math
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
 * A Discretization is built upon a Mesh object. Derived classes
 * are meant to contain all members and routines that are necessary to define
 * the discrete representation of a solution and the operations necessary for
 * aiding in the construction a linear system to solve.
 */
class Discretization
{
public:
  const std::shared_ptr<grid::Mesh> mesh;
  const DiscretizationMethod type;

public:
  explicit
  Discretization(const std::shared_ptr<grid::Mesh> reference_mesh,
                 const DiscretizationMethod discretization_type)
      : mesh(reference_mesh), type(discretization_type)
  {}

public:
  /// Get the number of nodes in the discretization.
  virtual size_t n_nodes() const { return 0; }

  /// Get the number of DoFs in the spatial discretization.
  virtual size_t n_dofs(const size_t n_components) const { return 0; }

  /// Get the number of nodes per cell.
  virtual size_t nodes_per_cell() const { return 0; }

  /// Get the number of DoFs per cell.
  virtual size_t dofs_per_cell(const size_t n_components) const { return 0; }

   /// Get the location of the nodes on a cell.
  virtual std::vector<grid::Point> nodes(const grid::Cell& cell) const = 0;
};

}
#endif //SPATIAL_DISCRETIZATION_H

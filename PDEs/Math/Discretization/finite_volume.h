#ifndef FINITE_VOLUME_H
#define FINITE_VOLUME_H

#include "Discretization/discretization.h"
#include "Grid/cell.h"

namespace math
{

/**
 * \brief A class for finite volume (FV) discretizations.
 *
 * Finite volume discretizations use cell-averaged unknowns, colloquially
 * located at the cell centers. This method is generally derived by assuming
 * this form of the solution, plugging it into the governing law, and
 * integrating over each Cell of the Mesh. This method primarily relies on
 * geometric information contained within the Cell and Face objects and does
 * not require many additional routines.
 */
class FiniteVolume : public Discretization
{
public:   /*---------- Constructors, Destructors, and Assignment ----------*/

  /// Default constructor.
  explicit FiniteVolume(std::shared_ptr<grid::Mesh> reference_mesh)
      : Discretization(reference_mesh, DiscretizationMethod::FINITE_VOLUME)
  {}

public: /*---------- Routines ----------*/

  /**
   * \brief Get the number of nodes in the FV discretization.
   *
   * For FV methods, this is the number of Cell objects on the Mesh.
   */
  size_t n_nodes() const override
  { return mesh->cells.size(); }

  /**
   * \brief Get the number of DoFs in the FV discretization.
   *
   * For FV methods, this the number of Cell objects on the Mesh multiplied by
   * the number of solution components \p n_components, such that:
   * \f[ n_{dofs} = n_{comp} n_{nodes} = n_{comp} n_{cells}. \f]
   *
   * \param n_components The number of solution components.
   */
  size_t n_dofs(const size_t n_components) const override
  { return n_components * n_nodes(); }

  /**
   * \brief Get the number of nodes per cell in the FV discretization.
   *
   * For FV methods, this is 1.
   */
  size_t nodes_per_cell() const override
  { return 1; }

  /**
   * \brief Get the number of DoFs per cell in the FV discretization.
   *
   * For FV methods, this is the number of solution componentd \p n_components.
   *
   * \param n_components The number of solution components.
   */
  size_t dofs_per_cell(const size_t n_components) const override
  { return n_components * nodes_per_cell(); }

  /**
   * \brief Get the nodes in the FV discretization.
   *
   * For FV methods, this is the Centroid of each Cell object.
   */
  std::vector<grid::Point> nodes(const grid::Cell& cell) const override
  { return std::vector<grid::Point>(1, cell.centroid); }

};

}
#endif //FINITE_VOLUME_H

#ifndef FV_H
#define FV_H

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

  /** Return the number of nodes per cell.
   *  For FV discretizations, there is 1 node per cell. */
  size_t nodes_per_cell() const override
  { return 1; }

  /** Return the number of DoFs per cell.
   *  This returns the number of nodes per cell multiplied by the number of
   *  solution components. For FV discretizations, this returns the number of
   *  solution components. See \ref nodes_per_cell. */
  size_t dofs_per_cell(const size_t n_components) const override
  { return n_components * nodes_per_cell(); }

  /** Return the number of nodes in the discretization.
   * For FV discretizations, this is equivalent to the number of cells. */
  size_t n_nodes() const override
  { return mesh->cells.size(); }

  /** Return the number of DoFs in the discretization.
   *  This returns the number of nodes multiplied by the number of solution
   *  components. For FV discretizations, this returns the number of cells times
   *  the number of solution components. See \ref n_nodes.*/
  size_t n_dofs(const size_t n_components) const override
  { return n_components * n_nodes(); }

  /** Return the location of the nodes on the specified Cell. For FV
   *  discretizations, this returns the Cell centroid. */
  std::vector<grid::Point> nodes(const grid::Cell& cell) const override
  { return std::vector<grid::Point>(1, cell.centroid); }

  void
  make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                        const size_t n_components = 1,
                        const bool is_coupled = false) const override;
};

}
#endif //FV_H

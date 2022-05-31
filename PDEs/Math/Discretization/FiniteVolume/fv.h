#ifndef FV_H
#define FV_H

#include "Discretization/discretization.h"
#include "cell.h"

#include <cinttypes>


namespace Math
{

  /**
   * A class for finite volume (FV) discretizations.
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
  public:

    /**
     * Default constructor.
     */
    explicit
    FiniteVolume(std::shared_ptr<Grid::Mesh> reference_mesh);

  public:

    /**
     * Return the number of nodes per cell. For FV discretizations, there is 1
     * node per cell.
     */
    size_t
    nodes_per_cell() const override;

    /**
     * Return the number of DoFs per cell. This returns the number of nodes per
     * cell multiplied by the number of solution components. For FV
     * discretizations, this returns the number of solution components.
     * \see FiniteVolume::nodes_per_cell */
    size_t
    dofs_per_cell(const size_t n_components) const override;

    /**
     * Return the number of nodes in the discretization. For FV discretizations,
     * this is equivalent to the number of cells. */
    size_t
    n_nodes() const override;

    /**
     * Return the number of DoFs in the discretization. This returns the number
     * of nodes multiplied by the number of solution components. For FV
     * discretizations, this returns the number of cells times the number of
     * solution components.
     * \see FiniteVolume::n_nodes
     */
    size_t
    n_dofs(const size_t n_components) const override;

    /**
     * Return the location of the nodes on the specified Cell. For FV
     * discretizations, this returns the Cell centroid.
     */
    std::vector<Grid::Point>
    nodes(const Grid::Cell& cell) const override;

    /**
     * Define the sparsity pattern. This routine defines the column indices of
     * non-zero entries per row for a problem with the specified number of
     * components. If the \p is_coupled flag is insert to \p true, it is assumed
     * that all components are coupled to one another, otherwise, it is assumed
     * that the system is uncoupled in all components. For FV discretizations,
     * each cell's node is couple to only the neighboring cell nodes.
     *
     * \param pattern The column indices per row to allocate for.
     * \param n_components The number of components in the solution.
     * \param is_coupled A flag for allocating storage for coupling between
     *   solution components.
     */
    void
    make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                          const size_t n_components = 1,
                          const bool is_coupled = false) const override;
  };

}
#endif //FV_H

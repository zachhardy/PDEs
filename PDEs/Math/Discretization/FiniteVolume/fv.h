#ifndef FV_H
#define FV_H

#include "Discretization/discretization.h"
#include "cell.h"

#include <cstddef>


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
    explicit FiniteVolume(std::shared_ptr<Grid::Mesh> reference_mesh);

    size_t n_nodes() const override;
    unsigned int nodes_per_cell() const override;

    size_t n_dofs(const unsigned int n_components) const override;
    unsigned int dofs_per_cell(const unsigned int n_components) const override;

    /** For FV, the only Node is the Centroid of the Cell. */
    std::vector<Grid::Point> nodes(const Grid::Cell& cell) const override;

    /**
     * Define the sparsity pattern. This routine defines the column indices of
     * non-zero entries per row for a problem with the specified number of
     * components. If the \p is_coupled flag is set to \p true, it is assumed
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
                          const unsigned int n_components = 1,
                          const bool is_coupled = false) const override;
  };

}
#endif //FV_H

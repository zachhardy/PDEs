#ifndef FV_H
#define FV_H

#include "Discretization/discretization.h"
#include "cell.h"

#include <cstddef>


namespace PDEs
{
  namespace Math
  {

    /**
     * Implementation of a finite volume (FV) discretization.
     *
     * A finite volume discretization uses cell-averaged unknowns, located at the
     * cell centers. This type is generally derived by assuming this form of the
     * solution, plugging it into the governing law, and  integrating over each
     * cell of the mesh. This method primarily relies on geometric information
     * contained within the cell and face objects and does not require many
     * additional routines.
     */
    class FiniteVolume : public Discretization
    {
    public:
      explicit FiniteVolume(std::shared_ptr<Grid::Mesh> reference_mesh);

      size_t
      n_nodes() const override;

      unsigned int
      nodes_per_cell() const override;

      size_t
      n_dofs(const unsigned int n_components) const override;

      unsigned int
      dofs_per_cell(const unsigned int n_components) const override;

      std::vector<Grid::Node>
      nodes(const Grid::Cell& cell) const override;

      void
      make_sparsity_pattern(std::vector<std::vector<size_t>> pattern,
                            const unsigned int n_components = 1,
                            const bool is_coupled = false) const override;
    };

  }
}
#endif //FV_H
